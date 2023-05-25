% CPPAnalysis
% run GetCPP first

if ~exist('filter_tf','file'); eeglab; close all; end
%%

clear; clc; close all;

%%%%%% set options
useCSD = 1;
excludeBadPps = 1; % remove pps with <640 good trials?
excludeTooFew = 1; % remove pps with <10 per cond*cert
excludeByRT = 1; % remove trials outside [100 1500] ms
doFilt = 1; % whether to plot with loPass filtered data
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty

if useCSD
    cppName = 'GetCPPCSD.mat';
    cppTopoName = 'GetCPPTopoCSD.mat';
    mapLims = [-15 15];
else
    cppName = 'GetCPP.mat';
    cppTopoName = 'GetCPPTopo.mat';
    mapLims = [-3 3];
end

% load stuff
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));
load(fullfile(outFolder, cppName),'cpp','cppStim','eeg','stimWin','preRespWindow','chans5','cppFilt','cppStimFilt');
load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');
load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels','behTab');
load(fullfile('./Saves/',cppTopoName),'cppMeanTopo','wins','cppRTTopo','cppTopoAmplWindow');


%% remove stuff - flagged, outliers, bad pps, too few trials per bin

toRemove = sq(any(isnan(cpp),2)); % flagged trials or cpp outliers

% ppsToExclude = [2 16]; % can manually remove these?
% ppsToExclude(16)=1; % noisy low numbers of low conf3
if excludeBadPps
    toRemove(ppsToExclude,:) = 1;     
end

if excludeByRT
    rtLims = [100 1540];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

if excludeCoMFromCert
    %%%% set any certainty with CoM to NaN - will be ignored everywhere
    behData.certainty(behData.CoM==1) = NaN; % will have to update in behTab
    toRemove(behData.CoM==1) = 1;
end

if excludeTooFew
    
    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);

    toRemove(toRemove2) = 1;
% 
%     % remove those people entirely, not just for that cond
%     toRemove(any(CountUnique(groupMeans(behData.certainty,2,behData.cond,'dim'),3)<10,[2 3]),:) = 1;
  
end


% actually remove from data now
cpp(repmat(permute(toRemove,[1,3,2]),1,size(cpp,2),1)) = NaN;
cppFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppFilt,2),1)) = NaN;
cppStim(repmat(permute(toRemove,[1,3,2]),1,size(cppStim,2),1)) = NaN;
cppStimFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppStimFilt,2),1)) = NaN;
cppMeanTopo(repmat(toRemove,1,1,128,size(cppMeanTopo,4))) = NaN;
cppRTTopo(repmat(toRemove,1,1,128)) = NaN;
cppTopoAmplWindow(repmat(toRemove,1,1,128,2)) = NaN;

% also need to exclude from behData, as will use that to exclude low cell counts
behNames = fieldnames(behData);
for i = 1:length(behNames)
    behData.(behNames{i})(toRemove) = NaN;
end

%% pre-response baseline?

% % re-baseline to pre-response
% baselineInds = isBetween(eeg.respTimes, [-100 0]);
% cppFilt = cppFilt - nanmean(cppFilt(:,baselineInds,:,:),2);
% cpp = cpp - nanmean(cpp(:,baselineInds,:,:),2);

%% define windows and labels

amplWindows = [-140 -50; 700 1000]; % ms
slopeWindows = [-310 -100; 700 1000];

cmap = crameri('vik');

factors = {'acc','certainty','CoM','confInR1','conf3'};
nF = 3;%length(factors);
labels.diffNames = {'Error - Correct' , 'Low - High', 'CoM - NoCoM','Low - High','Low - High' };
diffFlips = [-1 -1 1 -1 -1]; % invert some to match order

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
        [flipud(crameri('roma',6)); .2 .2 .2];
        [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];};

%% filter out 30Hz SSVEP? just for plotting, not analysis
if ~doFilt

    cppFilt = cpp;
    cppStimFilt = cppStim;

end

%% get means + slopes

cppVars.ampl = sq(nanmean(cpp(:,isBetween(eeg.respTimes,amplWindows(1,:)),:),2));
cppVars.amplPost = sq(nanmean(cpp(:,isBetween(eeg.respTimes,amplWindows(2,:)),:),2));

cppVars.slope = FitCPPSlope(cpp, slopeWindows(1,:), eeg.respTimes); % get slopes
cppVars.slopePost = FitCPPSlope(cpp, slopeWindows(2,:), eeg.respTimes); % get slopes

cppVarNames = fieldnames(cppVars);
nCppVars = length(cppVarNames);

varsByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),cppVars,'UniformOutput',0); % split by evidence condition
nF=2;

%% make lme table

if excludeCoMFromCert % exclude from this table
    behTab.certainty(behTab.CoM>0) = NaN;
end

% remove trials from loaded behdata table
behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));

% ampl + slope into table
cppTab = struct2table(structfun(@(x) col(x), cppVars, 'UniformOutput',0));

regTab = horzcat(behTab, cppTab); % combine
regNames = regTab.Properties.VariableNames;

isLogistic = cellRegexpi(regNames, 'Logistic')>0; % which are for logistic regr, don't zscore those
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

nF = 4 - 2*excludeCoMFromCert;%length(factors); % reset

%% statistics in the paper:

% switch to 'slope' to get the stats for the supplementary materials
cppName = 'ampl';

% 1st/2nd order acc + RT
% make table with 1st/2nd order in one column, repeat behData, combine pre/post choice CPP 
regTab2 = [ [ones(size(behTab,1),1); ones(size(behTab,1),1)*2], repmat(table2array(behTab),2,1), reshape(table2array(cppTab),[],2) ];
isLogistic2 = [true ~isLogistic([1:end-4 end-3 end-1])];
regTab2(:,isLogistic2) = nanzscore(regTab2(:,isLogistic2));
regTab2 = array2table(regTab2, 'VariableNames', [{'time'}, regNames(1:end-4) cppVarNames([1 3])']);
% get acc+RT from each resp into one
regTab2.accLogisticBoth = [regTab.accLogistic; regTab.confAccLogistic];
regTab2.RTLogBoth = nanzscore([behTab.RTLog; behTab.confRTLog]);


% this will apply logistic/normal regression to a cell array, depending on
% whether 'Logistic' is in the formula
fitglmeCell = @(x) fitglme(regTab, x, glmeArgs{ ~isempty(regexp(x,'Logistic'))+ 1}{:});
fitglmeCell2 = @(x) fitglme(regTab2, x, glmeArgs{ ~isempty(regexp(x,'Logistic'))+ 1}{:}); % uses regTab2
% for eachcond sep, i=0 | 1
fitglmeCellSep = @(x, i) fitglme(regTab(regTab.cond>0 == i,:), x, glmeArgs{ ~isempty(regexp(x,'Logistic'))+ 1}{:}); 

% first-order/pre choice CPP
%%%%%% NOTE: if excludeCoMFromCert=1, these are in the main text of the
%%%%%% paper, if it is 0, they are in the supplementary
paperFormulaePre = {
     '%s ~ 1 + acc*cond + (1 | pp)'; % cpp ampl
     '%s ~ 1 + certainty*cond + (1 | pp)';
     '%s ~ 1 + acc*certainty*cond+ (1 | pp)'; % 3-way 
     };
paperFormulaePre = cellfun(@(x) sprintf(x, cppName), paperFormulaePre, 'UniformOutput',0);
paperPreFits = cellfun(fitglmeCell, paperFormulaePre,'UniformOutput',0);

% run these separate LME in each ev cond separately
% excludeCoMFromCert
paperFormulaePreSep = {'%s ~ 1 + certainty + (1 | pp)';};
paperFormulaePreSep = cellfun(@(x) sprintf(x, cppName), paperFormulaePreSep, 'UniformOutput',0);

paperPreSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaePreSep,2,1), {0;1}, 'UniformOutput',0);


if ~excludeCoMFromCert % only do second-order if CoM trials are included, 

    % does acc & RT change between 1st/2nd resp? and how does post-dec ev
    % affect it?
    % don't excludeCoMFromCert
    paperFormulaeBeh2 = {'accLogisticBoth ~ 1 + cond*time + (1 | pp)';
        'RTLogBoth ~ 1 + acc*cond*time + (1 | pp)'};
    paperBehFits2 = cellfun(fitglmeCell2, paperFormulaeBeh2,'UniformOutput',0);
    
    
    % don't excludeCoMFromCert for this
    % how do Ev & 1st-Acc affect certainty,CoM and final RT?
    paperFormulaeBeh = {'confInR1 ~ 1 + cond*acc + (1 | pp)';
         'CoMLogistic ~ 1 + cond*acc + (1 | pp)';
         'RTLog ~ 1 + cond*acc + (1 | pp)';};
    paperBehFits = cellfun(fitglmeCell, paperFormulaeBeh,'UniformOutput',0);

    % post-choice/2nd order CPP
    % don't excludeCoMFromCert
    paperFormulaePost = {'%sPost ~ 1 + cond*confInR1 + (1 | pp)';
        '%sPost ~ 1 + cond*conf3+ (1 | pp)';
        '%sPost ~ 1 + cond*certainty + (1 | pp)';
        };
    paperFormulaePost = cellfun(@(x) sprintf(x, cppName), paperFormulaePost, 'UniformOutput',0);
    paperPostFits = cellfun(fitglmeCell, paperFormulaePost,'UniformOutput',0);
    
    
    
    % don't excludeCoMFRomCert
    paperFormulaePostSep = {'%sPost ~ 1 + confInR1 + (1 | pp)'; 
                            '%sPost ~ 1 + conf3 + (1 | pp)';
                            '%sPost ~ 1 + certainty + (1 | pp)';};
    paperFormulaePostSep = cellfun(@(x) sprintf(x, cppName), paperFormulaePostSep, 'UniformOutput',0);
    paperPostSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaePostSep,1,2), {0 1; 0 1; 0 1}, 'UniformOutput',0);
    
    
end

% post-hoc tests: make certainty a categorical, re-run and use coefTest to
% compare levels of Certainty. Continued only
paperFormulaePosthoc = { '%s ~ 1 + certCat + (1 | pp)';  % in Interrupted
    '%sPost ~ 1 + certCat+ (1 | pp)'}; % in Continued
paperFormulaePosthoc = cellfun(@(x) sprintf(x, cppName), paperFormulaePosthoc, 'UniformOutput',0);

regTab.certCat = categorical(regTab.certainty);
certContrasts = [0 1 0; % low vs high
        0 0 1; % med vs high
        0 1 -1];  % low vs med

% just do one of these
for j = 2-excludeCoMFromCert
    f1 = fitglme(regTab(regTab.cond>=0 == j-1,:), paperFormulaePosthoc{j}, 'DummyVarCoding','effects');
    for i = 1:size(certContrasts,1)
        [postHocStats(i,1,j), postHocStats(i,2,j), postHocStats(i,3,j), postHocStats(i,4,j)] = coefTest(f1, certContrasts(i,:));
    end
end

% u = unique(regTab.confInR1);u(7:end)=[];
% regTab.confInR1_cat = [regTab.confInR1 == u'] * [1;2;3;4;5;6];
% regTab.confInR1_cat(isnan(regTab.confInR1)) = NaN;
% regTab.confInR1_cat = categorical(regTab.confInR1_cat);
% f1 = fitglme(regTab(regTab.cond>=0,:), 'amplPost ~ 1 + confInR1_cat + (1 | pp)', 'DummyVarCoding','effects');


%% save now

saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
saveName = sprintf('CPPAnalysis_%s_%s.mat', saveOpts{1,useCSD+1}, saveOpts{2, excludeCoMFromCert+1});

save(fullfile(outFolder, saveName));

disp('saved');


%% start plots - don't run these, instead call PlotPaperFigures.m
% %% plot per pp
% 
% figure();
% [r,c] = GetSubPlotShape(fileInfo.nPP);
% for iPP = 1:fileInfo.nPP
%     
%     subplot(r,c,iPP);
%     errorBarPlot(permute(cppFilt(iPP,:,:),[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes);
%     
%     xline(0,':k');
%     yline(0, ':k');
%     
%     title(sprintf('%s: %s', fileInfo.ppID{iPP}));%, cppChans{iPP}));
%     
% end
% % makeSubplotScalesEqual(r,c,1:fileInfo.nPP);
% 
% 
% % %% all on one - includes tobeexcluded pps
% % 
% % figure();
% % c = get(gca,'ColorOrder');
% % for i = 1:2
% %     hold on; 
% %     h{i} = errorBarPlot(permute(cpp(ppsToExclude==(i-1),:,:),[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes,'color',c(i,:));
% % end
% % 
% % xline(0,':k');
% % yline(0, ':k');
% % ylabel('\muV'); xlabel('time from resp 1 (ms)');
% % box off;
% % 
% % title('grand mean');
% % legend([h{1}{1,2}, h{2}{1,2}], {'Include','Exclude'},'Location','Best');
% 
% %% plot all included pps together
% 
% figure();
% h = errorBarPlot(permute(cppFilt,[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes);
% 
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% 
% title('grand mean');
% 
% 
% 
% %% average over pps
% 
% figure();
% 
% errorBarPlot(nanmean(cppFilt,3),'area',1,'xaxisvalues',eeg.respTimes);
% 
% xline(0,':k');
% yline(0, ':k');
% xlines(amplWindows, '-k');
% xlines(slopeWindows, '-r');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% 
% title('grand mean');
% 
%% average over pps - split by cond

figure();
c = get(gca,'ColorOrder');
set(gca,'ColorOrder', [c(1:2,:); 1 0 0], 'nextplot','replacechildren');

resplocked = groupMeans(cppFilt,3, repmat(permute(behData.cond,[1,3,2]),1,length(eeg.respTimes)));
resplocked(:,:,3) = diff(resplocked,[],3);
h = errorBarPlot(resplocked,'area',1,'xaxisvalues',eeg.respTimes);
h{3,1}.LineStyle = '--';

% h{3,2}.Visible = 'off'; % just line
hold on;
fill(col(repmat(amplWindows(2,:),2,1)), [ylim fliplr(ylim)], 'k','FaceAlpha',.2,'LineStyle','none')

xline(0,':k');
yline(0, ':k');
% xlines(amplWindows, '-k');
% xlines(slopeWindows, '-r');
ylabel('CPP \muV/m^2'); xlabel('time from initial response (ms)');
box off;
legend([h{:,1}], [labels.cond, 'difference wave'],'Location','Best');
xlim([-500 1000]);
title('grand mean');


% %% plot CPP by RT bin
% 
% % 2 or 3 speed bins? uncomment one below to choose
% binNames = {'fast','slow'};
% % binNames = {'fast','medium','slow'};
% p = linspace(0,100, length(binNames)+1);
% p(1) = [];
% prc = prctile(behData.RT, p,2);
% binInds = sum(behData.RT > permute(prc,[1,3,2]),3);
% 
% cppByRT = groupMeans(cppFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(cpp,2),1),'dim'); %[npp t bins tr]
% 
% figure();
% subplot(1,2,1);
% cppStimByRT = groupMeans(cppStimFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(cppStim,2),1),'dim'); %[npp t bins tr]
% h = errorBarPlot(nanmean(cppStimByRT,4),'area',1,'xaxisvalues',eeg.stimTimes);
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from stim (ms)');
% box off;
% title('evidence locked');
% 
% subplot(1,2,2)
% h = errorBarPlot(nanmean(cppByRT,4),'area',1,'xaxisvalues',eeg.respTimes);
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% title('response locked');
% legend([h{:,1}], binNames,'Location','Best');
% 
% %% plot CPP by RT bin, and acc
% 
% % 2 or 3 speed bins? uncomment one below to choose
% accByRT = repmat(permute(groupMeans(behData.acc, 2, binInds,'dim'),[1 4 2 3]),1,size(cpp,2),1,1); %[npp t bins acc tr]
% cppStimByRTAcc = groupMeans(cppStimByRT, 4, accByRT(:,1:size(cppStim,2),:,:)); %[npp t bins tr]
% cppByRTAcc = groupMeans(cppByRT,4,accByRT); %[pp t rt acc]
% 
% figure();
% for i = 1:2
%     subplot(2,2,i*2-1);
%     
%     h = errorBarPlot(cppStimByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.stimTimes);
%     xline(0,':k');
%     yline(0, ':k');
%     ylabel('\muV'); xlabel('time from stim (ms)');
%     box off;
%     title(['evidence locked, ' labels.acc{i}]);
%     
%     subplot(2,2,i*2)
%     h = errorBarPlot(cppByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.respTimes);
%     xline(0,':k');
%     yline(0, ':k');
%     ylabel('\muV'); xlabel('time from resp 1 (ms)');
%     box off;
%     title(['response locked, ' labels.acc{i}]);
%     legend([h{:,1}], binNames,'Location','Best');
% end
% %% indiv plot RT splits
% 
% figure();
% for iPP = 1:fileInfo.nPP
%     subplot(6,5,iPP);
%     h = errorBarPlot(permute(cppByRT(iPP,:,:,:),[4,2,3,1]),'area',1,'xaxisvalues',eeg.respTimes);
%     xline(0,':k');
%     yline(0, ':k');
%     ylabel('\muV'); xlabel('time from resp 1 (ms)');
%     box off;
%     title(iPP);
% end    
% legend([h{:,1}], binNames,'Location','Best');
% 
%% topoplot at each window
chanInds = 0;

nWins = length(wins)-1;
[r,c] = GetSubPlotShape(nWins);
figure();
for i = 1:nWins
    subplot(r,c,i); 
    topoplot(sq(nanmean(nanmean(cppMeanTopo(:,:,:,i),2),1)),...
        eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

    title(wins(i+1));
end
colorbar

%% topoplot at each window - cont-inter

% get cont - int
cppMeanTopoByCond = groupMeans(cppMeanTopo,2,repmat(behData.cond,1,1,128,20)); %[pp cond chan time]
cppMeanTopoByCond(:,3,:,:) = diff(cppMeanTopoByCond(:,1:2,:,:),[],2); % diff

nWins = length(wins)-1;
[r,c] = GetSubPlotShape(nWins);
for j = 1:3
    figure();
    for i = 1:nWins
        subplot(r,c,i); 
        topoplot(sq(nanmean(cppMeanTopoByCond(:,j,:,i),1)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});
        title(wins(i+1));
    end
    colorbar
    if j < 3; SuperTitle(labels.cond(j)); else SuperTitle('Continued-Interrupted');end
end

% 
% %% show just a few times, [int cont diff]
% 
% winInds = find(any(wins(2:end) == (100:200:1000)'));
% nWins = length(winInds);
% figure();
% for i = 1:nWins
%     for j = 1:3
%         subplot(3,nWins,(j-1)*nWins+i); 
%         topoplot(sq(nanmean(cppMeanTopoByCond(:,j,:,winInds(i)),1)),...
%             eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
%             'emarker',{'.','k',10,1});
%         if j==1; title(wins(winInds(i))); end
%         if i==1
%             if j<3; ylabel(labels.cond(j)); else ylabel('Cont-Int');end
%         end
%     end
% end
% colorbar('Location','North','Position',[0.3631 0.3333 0.4387 0.0535]);
% % SuperTitle('Cont-Interr');
% 
%% plot topo - just low-conf continued

chanInds = 0;

cppTopoContLow = cppMeanTopo;
cppTopoContLow(repmat(behData.certainty>1 | behData.cond==1, 1,1,eeg.nChans,nWins)) = NaN; % remove all but continued-low conf

nWins = length(wins)-1;
[r,c] = GetSubPlotShape(nWins);
figure();
for i = 1:nWins
    subplot(r,c,i); 
    topoplot(sq(nanmean(nanmean(cppTopoContLow(:,:,:,i),2),1)),...
        eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

    title(wins(i+1));
end
colorbar
SuperTitle('Continued, Low Confidence');

% 
% 
% %% Fig 3.4
% figure();
% 
% subplot(2,2,1); % topoplot of activity at RT
% 
% 
% topoplot(sq(nanmean(nanmean(cppRTTopo,2),1)),...
%     eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',mapLims,...
%     'emarker',{'.','k',10,1});
% title('at RT')
% 
% subplot(2,2,2); % stimlocked trace by cond
% stimlocked = groupMeans(cppStimFilt,3,repmat(permute(behData.cond2,[1,3,2]),1,size(cppStim,2)),'dim'); %[pp t cond tr]
% eeg.stimTimes = eeg.epochTimes(isBetween(eeg.epochTimes, stimWin));
% 
% win = [-200 1200]; % for plotting
% h = plot(eeg.stimTimes, sq(nanmean(nanmean(stimlocked(:,isBetween(eeg.stimTimes,stimWin),:,:),4),1)),'LineWidth',2);
% xline(0);
% legend(h, labels.cond2,'Location','Best');
% title('evidence onset');
% xlim(win);
% 
% subplot(2,2,3);
% data = groupMeans(cppFilt,3,repmat(permute(behData.cond2,[1,3,2]),1,size(cpp,2)),'dim'); %[pp t cond tr]
% 
% win = [-500 1000];
% h = plot(eeg.respTimes(isBetween(eeg.respTimes,win)), sq(nanmean(nanmean(data(:,isBetween(eeg.respTimes,win),:,:),4),1)),'LineWidth',2);
% xline(0);
% legend(h, labels.cond2,'Location','Best');
% title('first-order response');
% 
% % do a topoplot - of 2nd cpp
% % axes('Position',[.7 .7 .2 .2])
% subplot(2,2,4);
% box off
% win = [500 800];
% cppMeanTopo500 = nanmean(cppMeanTopo(:,:,:,17:19),4); % average over win
% topoplot(sq(nanmean(nanmean(cppMeanTopo500,2),1)),...
%     eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',mapLims,...
%     'emarker',{'.','k',10,1});
% title([num2str(win) ' ms']);

%% make fig3.5

% split [stimlocked cpp] by [conf; acc; CoM], do a topoplot of 600-800ms
% post decision

cppMeanTopo600 = nanmean(cppMeanTopo(:,:,:,18:19),4); % average over 600-800ms
stimWin = [-200 1200];
respWin = [-500 1000];
ylims = [-1.5 4; -1.4 5];

nF=2;
figure();
for i = 1:nF
    subplot(nF,3,i*3-2);
    
    % split stimlocked by cond
    stimlocked = groupMeans(cppStimFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(cppStim,2))); %[pp t cond]
    stimlocked(:,:,end+1) = diff(stimlocked(:,:,[1 end]),[],3) * diffFlips(i); % get diff
    
    % plot
    h = plot(eeg.stimTimes(isBetween(eeg.stimTimes,stimWin)), sq(nanmean(stimlocked(:,isBetween(eeg.stimTimes,stimWin),:,:),1)), 'LineWidth',2);
    for j = 1:length(h); h(j).Color = cols{i}(j,:);end
    
    xline(0);yline(0);
    legend(h, [labels.(factors{i}) labels.diffNames{i}],'Location','Best');
    if i==1; title('evidence onset');end
    xlim(stimWin); %ylim(ylims(1,:));
    ylabel(factors{i});
    
    
    
    subplot(nF,3,i*3-1);
        % split cpp by cond
    resplocked = groupMeans(cppFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(cpp,2))); %[pp t cond]
    resplocked(:,:,end+1) = diff(resplocked(:,:,[1 end]),[],3) * diffFlips(i); % get diff

    
    h = plot(eeg.respTimes(isBetween(eeg.respTimes,respWin)), sq(nanmean(resplocked(:,isBetween(eeg.respTimes,respWin),:,:),1)), 'LineWidth',2);
    for j = 1:length(h); h(j).Color = cols{i}(j,:);end
    xline(0);yline(0);
%     xlines(amplWindows,'-k');xlines(slopeWindows, '-r');
    if i==1; title('first-order response'); end
    xlim(respWin); %ylim(ylims(2,:));
    
    
    % add a topoplot
    subplot(nF,3,i*3);

    % split by cond
    respTopo = groupMeans(cppMeanTopo600,2,repmat(behData.(factors{i}),1,1,128));
    respTopo = diff(respTopo(:,[1 end],:),[],2) * diffFlips(i);
    topoplot(sq(nanmean(respTopo,1)),...
        eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
        'emarker',{'.','k',10,1});
    title('600:800ms post dec');
    
end


%% plot resplocked by factors, separately for cont/int



% split cppFilt by Cond
cppFiltCond = groupMeans(cppFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(cpp,2)),'dim'); %[pp t cond tr]
figure();
for i = 1:nF
%     figure();
   
    factorByCond = repmat(permute(groupMeans(behData.(factors{i}),2,behData.cond,'dim'),[1 4 2 3]),1,length(eeg.respTimes),1,1); % to match cppFiltCond

    % split cpp by cond
    resplocked = groupMeans(cppFiltCond,4,factorByCond); %[pp t cond factor]
    resplocked(:,:,:,end+1) = diff(resplocked(:,:,:,[1 end]),[],4) * diffFlips(i); % get diff

    
    for j = 1:2
        subplot(nF,2,(i-1)*2+j);
%         subplot(1,2,j);
        set(gca,'ColorOrder', cols{i}, 'nextplot','replacechildren');
        h = plot(eeg.respTimes(isBetween(eeg.respTimes,respWin)), sq(nanmean(resplocked(:,isBetween(eeg.respTimes,respWin),j,1:end),1)), 'LineWidth',2);
%         h = errorBarPlot(sq(resplocked(:,isBetween(eeg.respTimes,respWin),j,1:end-1)), 'area',1,'xaxisvalues',eeg.respTimes(isBetween(eeg.respTimes,respWin)));
%         for j = 1:length(h); h(j).Color = cols{i}(j,:);end
        xline(0);yline(0);
    %     xlines(amplWindows,'-k');xlines(slopeWindows, '-r');
        if i%==1; 
            title(labels.cond(j)); 
        end

        if j==1
            ylabel('CPP \muV'); 
            legend(h, [labels.(factors{i}) labels.diffNames{i}],'Location','NorthEast');
%             legend([h{:,1}], [labels.(factors{i}) labels.diffNames{i}],'Location','NorthEast');
        end

        xlim(respWin); %ylim(ylims(2,:));
    end    
    makeSubplotScalesEqual(nF,2)
         
end

%% can also perm plot stuff

% cppFiltCond = groupMeans(cppFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(cpp,2)),'dim'); %[pp t cond tr]
% 
% fig = figure();
% % 
% % for i = 1:5
% % %     subplotInds = [5 2 i*2-1; 5 2 i*2]; % column per fig type
% %     figure(fig);
% %     subplotInds = [3,2,i]; % column per fig type
% %     myPermPlots(cppFiltCond, factors{i}, eeg.respTimes, behData, diffFlips(i), cols{i}, respWin, labels, labels.diffNames(i),subplotInds);
% % end
% 
% clusterArgs = {'cluster',1,'clustermethod','mass','nperms',10000};
% n = {'acc','certainty'};%,'conf3','acc','CoM','confInR1'};
% [r,c] = GetSubPlotShape(length(n));
% for i = 1:length(n)
%     j = find(strcmp(factors, n(i)));
% %     subplotInds = [5 2 i*2-1; 5 2 i*2]; % column per fig type
%     figure(fig);
%     subplotInds = [r,c,i]; % column per fig type
%     myPermPlots(cppFiltCond, factors{j}, eeg.respTimes, behData, diffFlips(j), cols{j}, respWin, labels, labels.diffNames(j),subplotInds,[],0,clusterArgs);
% end

%% plot acc*fac for conditions separately
clusterArgs = {'cluster',1,'clustermethod','mass','nperms',5000};

respWin = [-500 1000];
tInds = isBetween(eeg.respTimes, respWin);

% split by acc too
cppFiltCond = groupMeans(cppFilt(:,tInds,:),3,repmat(permute(behData.cond,[1,3,2]),1,sum(tInds)),'dim'); %[pp t cond tr]
accByCond = repmat(permute(groupMeans(behData.acc,2,behData.cond,'dim'),[1 4 2 3]),1,size(cppFiltCond,2),1,1); % to match cppFiltCond
cppFiltCondAcc = groupMeans(cppFiltCond,4,accByCond,'dim'); %[pp time cond acc tr] 


lineStyles = {':','-'};
for iF = 4
    figure();
    factor = factors{iF};
    facByCond  = repmat(permute(groupMeans(behData.(factor),2,behData.cond,'dim'),[1 4 2 3]),1,size(cppFiltCond,2),1,1); % to match cppFiltCond
    facByCondAcc = groupMeans(facByCond,4,accByCond,'dim');
    cppFiltCondAccFac = groupMeans(cppFiltCondAcc, 5, facByCondAcc); %[pp time cond acc fac]

    % permute here to swap cond/acc
%     cppFiltCondAccFac = permute(cppFiltCondAccFac,[1,2,4,3,5]);

    % get diff waves - 
%     cppFiltCondAccFac(:,:,:,:,end+1) = diff(cppFiltCondAccFac(:,:,:,:,[1 end]),[],5) .* diffFlips(iF);
    for j = 1 % acc
%         subplot(1,2,j);
%         set(gca,'ColorOrder', cols{iF}, 'nextplot','replacechildren');

        for k = 1:2 
            subplot(1,2,k);
%              subplot(2,2,(j-1)*2+k);
            set(gca,'ColorOrder', cols{iF}, 'nextplot','replacechildren');
                    hold on;

%             h = plot(eeg.respTimes(tInds), sq(nanmean(cppFiltCondAccFac(:,:,j,k,:),1)), 'LineWidth',2,'LineStyle',lineStyles{2});
            h = errorBarPlot(sq(cppFiltCondAccFac(:,:,j,k,:)),'area',1,'xaxisvalues',eeg.respTimes(tInds),'plotargs',{'LineStyle',lineStyles{2}});
            for ii=1:size(h,1)
%                 h{ii,2}.Visible='off'; % uncomment to remove err bars
                h{ii,2}.FaceAlpha = .1;
%                 h{ii,1}.LineWidth = 2;
            end
        
            title([labels.cond{j}, ',' labels.acc{k}]);
            xline(0);yline(0);
            ylabel(factor); 
    
%             % perm testing - get corr-Err, Com-none, interaction
%             clear diffWaves;
%             diffWaves(:,:,1) = sq(diff(nanmean(cppFiltCondAccFac(:,:,j,:,1:2),5),[],4)); % av over fac, effect of acc
%             diffWaves(:,:,2) = sq(diff(nanmean(cppFiltCondAccFac(:,:,j,:,1:2),4),[],5)) .* diffFlips(iF); % av over acc, effect of fac
%             diffWaves(:,:,3) = sq(diff(diff(cppFiltCondAccFac(:,:,j,:,1:2),[],5),[],4)) .* diffFlips(iF); % diff over fac, then acc
%             p = [];
%             for k = 1:3
%                 [~,p(k,:)] = permutationOLS(diffWaves(:,:   ,k),[],[],[],clusterArgs{:});
%                 pbar(p(k,:), 'xVals', eeg.respTimes(tInds), 'yVal',  min(ylim) + k/5, 'plotargs',{'LineWidth',5,'Color',circshift([1 0 0],k-1)});
%             end
        
        end
    end
%     clear legStyles;
%     for i = 1:length(h)-1
%         legStyles{1,i} = {'-','Color',cols{iF}(i,:)};
%     end
%     legStyles = [legStyles {{lineStyles{1},'Color', [.5 .5 .5]}, {lineStyles{2},'Color', [.5 .5 .5]}, {'-k'}} ];
%     legNames = [labels.(factor), labels.acc labels.diffNames{iF}];
%     emptyLegend(length(legNames), legStyles, {'LineWidth',2}, legNames, {'Location','Best'});   
    legend([h{:,1}], labels.(factors{iF}),'Location','Best')
end

% makeSubplotScalesEqual(1,2);

%% split by conf*acc*com*cert?

iF = find(strcmp(factors, 'confInR1')); % confInR1, will split into cert*CoM later

facByCond  = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.cond,'dim'),[1 4 2 3]),1,size(cppFiltCond,2),1,1); % to match cppFiltCond
facByCondAcc = groupMeans(facByCond,4,accByCond,'dim');
cppFiltCondAccFac = groupMeans(cppFiltCondAcc, 5, facByCondAcc); %[pp time cond acc fac]
% need to change order of conf6
cppFiltCondAccFac = cppFiltCondAccFac(:,:,:,:,[4 5 6 3 2 1]);
cppFiltCondAccCertCoM = reshape(cppFiltCondAccFac,30,[],2,2,3,2); %[pp time cond acc cert com


iF = 2; % reset to cert
figure();
for i = 2
    for j = 1:2
        for k = 1:2
            subplot(2,2,(j-1)*2+k);
            set(gca,'ColorOrder',cols{iF}(1:end-1,:),'nextplot','replacechildren');
            hold on;

            h = errorBarPlot(sq(cppFiltCondAccCertCoM(:,:,i,j,:,k)), 'area',1,'xaxisvalues',eeg.respTimes(tInds),'plotargs',{'LineStyle',lineStyles{k}});
%             h = plot(eeg.respTimes(tInds), sq(nanmean(cppFiltCondAccCertCoM(:,:,i,j,:,k),1)), 'LineWidth',2,'LineStyle',lineStyles{k});
        
        
            title([labels.cond{i}, ', ' labels.acc{j}]);
            xline(0);yline(0);
            ylabel(factors{iF}); 
        end
%         % perm testing - get corr-Err, Com-none, interaction
%         clear diffWaves;
%         diffWaves(:,:,1) = sq(diff(nanmean(cppFiltCondAccCertCoM(:,:,j,:,1:2),5),[],4)); % av over fac, effect of acc
%         diffWaves(:,:,2) = sq(diff(nanmean(cppFiltCondAccCertCoM(:,:,j,:,1:2),4),[],5)) .* diffFlips(iF); % av over acc, effect of fac
%         diffWaves(:,:,3) = sq(diff(diff(cppFiltCondAccCertCoM(:,:,j,:,1:2),[],5),[],4)) .* diffFlips(iF); % diff over fac, then acc
%         p = [];
%         for k = 1:3
%             [~,p(k,:)] = permutationOLS(diffWaves(:,:   ,k),[],[],[],clusterArgs{:});
%             pbar(p(k,:), 'xVals', eeg.respTimes(tInds), 'yVal',  min(ylim) + k/5, 'plotargs',{'LineWidth',5,'Color',circshift([1 0 0],k-1)});
%         end
        
    end
end

clear legStyles;
for i = 1:length(h)
    legStyles{1,i} = {'-','Color',cols{iF}(i,:)};
end
legStyles = [legStyles {{lineStyles{1},'Color', [.5 .5 .5]}, {lineStyles{2},'Color', [.5 .5 .5]}, {'-k'}} ];
legNames = [labels.(factors{iF}), labels.CoM labels.diffNames{iF}];
emptyLegend(length(legNames)-1, legStyles, {'LineWidth',2}, legNames, {'Location','Best'});   

makeSubplotScalesEqual(2,2);

% will draw the means within window below


%% correlate ampl/slopes
% 
% figure;
% scatterArgs = {'pearson',0,'plot_ci',1,'text',2};
% for i = 1:2
%     subplot(2,2,i);
% 
%     scatterRegress(col(cppVars.(cppVarNames{i})), col(cppVars.(cppVarNames{i+2})), scatterArgs{:});
%     xlabel(cppVarNames{i});
%     ylabel(cppVarNames{i+2});
% 
%     subplot(2,2,i + 2);
%     conditionalPlot(cppVars.(cppVarNames{i})', cppVars.(cppVarNames{i+2})');
%     xlabel(cppVarNames{i});
%     ylabel(cppVarNames{i+2});
% 
% end
% 
% % and correlate pre/post
% 
% figure;
% scatterArgs = {'pearson',0,'plot_ci',1,'text',2};
% for i = 1:2
%     subplot(2,2,i);
% 
%     scatterRegress(col(cppVars.(cppVarNames{i*2-1})), col(cppVars.(cppVarNames{i*2})), scatterArgs{:});
%     xlabel(cppVarNames{i*2-1});
%     ylabel(cppVarNames{i*2});
% 
%     subplot(2,2,i + 2);
%     conditionalPlot(cppVars.(cppVarNames{i*2-1})', cppVars.(cppVarNames{i*2})');
%     xlabel(cppVarNames{i*2-1});
%     ylabel(cppVarNames{i*2});
% 
% end
%% split by cond and then cert/acc/com



figure();
for iF = 1:nF
    % split by cond %[pp cond 3 tr]
    factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');

    for iV = 1:nCppVars

        % split by factor
        data = groupMeans(varsByCond.(cppVarNames{iV}), 3, factor); %[pp cond 3 factor]
        
        subplot(4,nF,(iV-1)*nF+iF)
        h = errorBarPlot(data, 'type','bar');
        set(gca,'XTickLabel',labels.cond);
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(i,:);
        end
        if iV==1
            legend(h, labels.(factors{iF}),'Location','Best'); 
            title(factors{iF});
        end
        if iF==1; ylabel(cppVarNames{iV}); end
        
        
    end
    
end

%% draw pre/post separately

for iP = 2:-1:1
    figure();
    for iF = 1:nF
        % split by cond %[pp cond 3 tr]
        factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');
    
        for iV = 1:2
            subplot(2,nF,(iV-1)*nF+iF)

            if iF==1 && iP==2 % no cond for acc
                data = groupMeans(cppVars.(cppVarNames{iV*2-(iP-1)}), 2, behData.(factors{iF})); %[pp factor]
            else % include cond
                % split by factor
                data = groupMeans(varsByCond.(cppVarNames{iV*2-(iP-1)}), 3, factor); %[pp cond 3 factor]
            end
            
            h = errorBarPlot(data, 'type','bar');
            if iF==1 && iP==2
                h.FaceColor = 'flat';
                h.CData = cols{iF}(1:end-1,:);
                set(gca,'XTickLabel',labels.(factors{iF}));
            else
                set(gca,'XTickLabel',labels.cond);
                for i = 1:length(h)
                    h(i).FaceColor = cols{iF}(i,:);
                end
            end


            
            if iV==1
                if iF==2 && iP==2
                    legend(h, labels.(factors{iF}),'Location','Best'); 
                end
                title(factors{iF});
            end
            if iF==1; ylabel(cppVarNames{iV*2-(iP-1)}); end
            
            
        end
    end    
    makeSubplotScalesEqual(2,2,1:2);
    makeSubplotScalesEqual(2,2,3:4);

end

%% split acc/err


behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
varsByCondAcc = structfun(@(x) groupMeans(x, 3, behDataByCond.acc,'dim'),varsByCond,'UniformOutput',0); %[pp cond acc tr]

% figure();
for iF = 2
    figure();
    factor = groupMeans(behDataByCond.(factors{iF}),3,behDataByCond.acc, 'dim');

    for iV = 1:2%nCppVars    
        

        % split by factor
        data = groupMeans(varsByCondAcc.(cppVarNames{iV*2}), 4, factor); %[pp cond acc factor]
%         data = permute(data,[1 3 2 4]); % need to permute as errbpl reshapes it, which messes order
%         data = reshape(data,30,4,[]);

        for iA = 1:2
            subplot(2,2,(iV-1)*2+iA);
            set(gca,'ColorOrder',cols{iF}(1:end-1,:),'nextplot','replacechildren');
            h = errorBarPlot(sq(data(:,:,iA,:)), 'type','bar');
            set(gca,'XTickLabel', labels.cond);
%             set(gca,'XTickLabel',col(strcat(repmat(labels.cond',1,2), repmat(labels.acc,2,1))'));
    %         [h(1:length(h)/2).FaceAlpha] = deal(.5); % make errors translucent
    %         for i = 1:length(h)
    %             h(i).FaceColor = cols{iF}(i,:);
    %         end
            if iV==1
                if iA==1
                    legend(h, [labels.(factors{iF})],'Location','Best'); 
                end
                title(labels.acc{iA})
            end
            ylabel(cppVarNames{iV*2}); 
        
        end 
        makeSubplotScalesEqual(2,2,iV*2-1:iV*2);
    end
    
end


%% split by certainty,com,acc,cond
% massively unequal trials numbers per cell (min mean is <1 int,corr,high,Change, max mean is
% 183)

behDataByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), behDataByCond, 'UniformOutput',0);
cppVarsByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), varsByCond,'UniformOutput',0);

amplPostByCondAccCertCoM = groupMeans(cppVarsByCondAcc.amplPost, 4, behDataByCondAcc.confInR1); % split by conf6
amplPostByCondAccCertCoM = amplPostByCondAccCertCoM(:,:,:,[4 5 6 3 2 1]); % invert these to make go low-high low-high
amplPostByCondAccCertCoM = reshape(amplPostByCondAccCertCoM,30,2,2,3,2); % split by CoM

% amplPostByCondAccCertCoM = sq(nanmean(cppFiltCondAccCertCoM(:,isBetween(eeg.respTimes(tInds), [800 950]),:,:,:,:),2));

figure();

for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        set(gca,'ColorOrder',cols{2}(1:end-1,:),'nextplot','replacechildren');

        h = errorBarPlot(permute(amplPostByCondAccCertCoM(:,i,j,:,:),[1 5 4 2 3]),'type','bar');

        title([labels.cond{i}, ', ' labels.acc{j}]);
        xlabel(factors{4}); 
        xticks(1:2); xticklabels(labels.CoM);
        ylabel('post-dec CPP ampl');

    end
end

legend(h, labels.certainty, 'Location','Best');
makeSubplotScalesEqual(2,2);

%% now for slopes

slopePostByCondAccCertCoM = groupMeans(cppVarsByCondAcc.slopePost, 4, behDataByCondAcc.confInR1); % split by conf6
slopePostByCondAccCertCoM = slopePostByCondAccCertCoM(:,:,:,[4 5 6 3 2 1]); % invert these to make go low-high low-high
slopePostByCondAccCertCoM = reshape(slopePostByCondAccCertCoM,30,2,2,3,2); % split by CoM

% amplPostByCondAccCertCoM = sq(nanmean(cppFiltCondAccCertCoM(:,isBetween(eeg.respTimes(tInds), [800 950]),:,:,:,:),2));

figure();

for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        set(gca,'ColorOrder',cols{2}(1:end-1,:),'nextplot','replacechildren');

        h = errorBarPlot(permute(slopePostByCondAccCertCoM(:,i,j,:,:),[1 5 4 2 3]),'type','bar');

        title([labels.cond{i}, ', ' labels.acc{j}]);
        xlabel(factors{4}); 
        xticks(1:2); xticklabels(labels.CoM);
        ylabel('post-dec CPP slope');

    end
end

legend(h, labels.certainty, 'Location','Best');
makeSubplotScalesEqual(2,2);

%% now both together

figure();

for i = 1:2
    if i==1
        data = amplPostByCondAccCertCoM;
    else
        data = slopePostByCondAccCertCoM;
    end
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        set(gca,'ColorOrder',cols{2}(1:end-1,:),'nextplot','replacechildren');

        h = errorBarPlot(reshape(permute(data(:,:,:,:,j),[1 2 3 4]), fileInfo.nPP, 4, 3),'type','bar');

        title(labels.CoM{j});
        set(gca,'XTickLabels',labels.acc);
        xlabel({'interrupt          continue';'accuracy'});
        ylabel(['post-dec CPP ' cppVarNames{i*2}]);

    end
    makeSubplotScalesEqual(2,2,i*2-1:i*2);
end

legend(h, labels.certainty, 'Location','Best');


%% line plot

figure();

for iCond = 1:2
    subplot(1,2,iCond)
    set(gca,'ColorOrder',cols{3}(1:2,:),'nextplot','replacechildren');
    h = errorBarPlot(reshape(permute(amplPostByCondAccCertCoM(:,iCond,:,:,:),[1,4,3,5,2]),30,[],2),'type','line','plotargs',{'LineStyle','none','LineWidth',2});
    JoinUpErrorPoints(h,reshape(1:length(h(1).XData),3,[])');
    xticks(1:6)
    xticklabels(repmat(labels.certainty,1,2))
    xtickangle(0);
    xlim([.5 6.5]);
    xlabel({'error              correct', 'certainty'});
    ylabel('post-dec CPP ampl');

    hold on
    for i =1:3
        h1(i,:) = plot([i i+3], sq(nanmean(amplPostByCondAccCertCoM(:,iCond,:,i,:),1)), 'o', 'Color', cols{2}(i,:),'LineWidth',2);
    end
    title(labels.cond(iCond));
end

legend([h h1(:,1)'], [labels.CoM labels.certainty],'Location','Best');
makeSubplotScalesEqual(1,2);


%% 2-way cond*factor

[anovas2, anovas2Tabs] = deal(struct());
for iF = 1:nF
    factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');

    for iV = 1:nCppVars
        
        data = groupMeans(varsByCond.(cppVarNames{iV}), 3, factor); %[pp cond 3 factor]
                
        anovas2.(factors{iF}).(cppVarNames{iV}) = rmanova(data, {'pp','cond',factors{iF}},'categorical',2:3, 'DummyVarCoding','effects');
        
    end
    anovas2Tabs.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovas2.(factors{iF}),'UniformOutput',0),'rowNames',anovas2.(factors{iF}).ampl.Term(2:end));
end

%% cont/inter separately

[anovasSep, anovasSepTabs] = deal(struct());
for iF = 1:nF
    factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');

    for iV = 1:nCppVars
        
        data = groupMeans(varsByCond.(cppVarNames{iV}), 3, factor); %[pp cond 3 factor]

        for j = 1:2        
            anovasSep.(labels.cond{j}).(factors{iF}).(cppVarNames{iV}) = rmanova(sq(data(:,j,:)), {'pp',factors{iF}},'categorical',[], 'DummyVarCoding','effects');
        end
        
    end
    for j = 1:2
        anovasSepTabs.(labels.cond{j}).(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasSep.(labels.cond{j}).(factors{iF}),'UniformOutput',0),'rowNames',anovasSep.(labels.cond{j}).(factors{iF}).ampl.Term(2:end));
    end
end

%% 3-way: pre/post * cond * [cert | acc | CoM]

% anovas3 = struct();
% for iF = 1:nF
%     factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');
% 
%     for iV = 1:2
%         
%         dataPre = groupMeans(varsByCond.(varNames{iV*2-1}), 3, factor); %[pp cond 3 factor]
%         dataPost = groupMeans(varsByCond.(varNames{iV*2}), 3, factor); %[pp cond 3 factor]
%         
%         data = cat(4, dataPre, dataPost);
%         
%         anovas3.(factors{iF}).(varNames{iV*2-1}) = rmanova(data, {'pp','cond',factors{iF},'time'},'categorical',2:4, 'DummyVarCoding','effects');
%         
%     end
%     anovas3Tabs.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovas3.(factors{iF}),'UniformOutput',0),'rowNames',anovas3.(factors{iF}).ampl.Term(2:end));
% 
% end



%% does single-trial lme have better sensitivity?



% %% look at init-dec things split by factors
% labels.confAcc = labels.acc;
% labels.confResp = {'def left','probs left','maybe left','maybe right','probs right','def right'};
% labels.scoretemp = labels.confResp;
% 
% behVarNames2 = {'acc','RT','confAcc','confRT','confInR1','certainty','conf3','CoM'};
% nBehVars2 = length(behVarNames2);
% 
% alpha = .01;
% 
% form = '%s ~ 1 + %s + (1 | pp)';
% fit = {};
% 
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = [1 3]
%     figure();
% 
%     for iB = 1:length(behVarNames2)
%         subplot(r,c,iB);
%         if any(strcmp(factors, behVarNames2{iB}))
%             set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%         end
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             % split by factor
%             data = groupMeans(cppVars.(cppVarNames{iV}), 2, behData.(behVarNames2{iB})); %[pp cond 3 factor]
%         
%         
%             h = errorBarPlot(data, 'type','bar');
%             ylabel(cppVarNames{iV});
%             xlabel(behVarNames2{iB})
% 
%         
%         else % cond plot
%             [~,~,~,~,h] = conditionalPlot(behData.(behVarNames2{iB})', cppVars.(cppVarNames{iV})');
%             xlabel(behVarNames2{iB});
%             ylabel(cppVarNames{iV});
% 
%         end
% 
%         SuperTitle(cppVarNames{iV});
% 
%         % do a regression and put star
%         fit{iB} = fitglme(regTab, sprintf(form, cppVarNames{iV}, behVarNames2{iB}));
%         if fit{iB}.Coefficients.pValue(strcmp(fit{iB}.CoefficientNames, behVarNames2{iB})) < alpha
%              text(0.5, 0.9, '*','Units','Normalized','FontSize',16);
%         end
%         
%         
%     end
%     cppBehRegs.beh.(cppVarNames{iV}) = StatsTableFromLMECells(fit,behVarNames2);
%     
% end
% 
% %% look at init-dec things split by factors + acc
% 
% 
% behDataByAcc = rmfield(behData, behNames(~ismember(behNames, [behVarNames2 {'cond'}]))); 
% behDataByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), behDataByAcc, 'UniformOutput',0);
% 
% cppVarsByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), cppVars,'UniformOutput',0);
% 
% form = '%s ~ 1 + acc*%s + (1 | pp)';
% fit = {};
% 
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = [1 3]
%     figure();
% 
%     for iB = 1:nBehVars2
%         subplot(r,c,iB);
%         if any(strcmp(factors, behVarNames2{iB}))
%             set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%         end
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             % split by AccFac
% 
%             data = groupMeans(cppVarsByAcc.(cppVarNames{iV}), 3, behDataByAcc.(behVarNames2{iB})); %[pp acc fac tr]
%         
%             
%             h = errorBarPlot(data, 'type','bar');
%             xticklabels(labels.acc);
%             ylabel(cppVarNames{iV});
%             if iB==1 || iB==nBehVars2; legend(h, labels.(behVarNames2{iB}),'Location','Best'); end
% 
%         
%         else % cond plot
%             
%             [~,~,~,~,h] = conditionalPlot(permute(behDataByAcc.(behVarNames2{iB}),[3,1,2]), permute(cppVarsByAcc.(cppVarNames{iV}),[3 1 2]));
%             xlabel(behVarNames2{iB});
%             ylabel(cppVarNames{iV});
%             legend(h(2:2:end), labels.acc,'Location','Best');
% 
%         end
% 
%         title(behVarNames2{iB});
% 
%         % do a regression and put star
%         if ~strcmp(behVarNames2{iB},'acc')
%             fit{iB} = fitglme(regTab, sprintf(form, cppVarNames{iV}, behVarNames2{iB}));
%             for j = 1:3
%                 if fit{iB}.Coefficients.pValue(j+1) < alpha
%                     t = text(0.5, (1-j/10), '*','Units','Normalized','FontSize',16,'Color',[1 2 3]==j);
%                 end
%             end
%         end
% 
%     end
%     cppBehRegs.behAcc.(cppVarNames{iV}) = StatsTableFromLMECells(fit,behVarNames2);
% end
% 
% %% look at post-dec things split by factors + cond
% 
% behDataByCond = rmfield(behData, behNames(~ismember(behNames, [behVarNames2 {'cond'}]))); 
% behDataByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), behDataByCond, 'UniformOutput',0);
% 
% cppVarsByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), cppVars,'UniformOutput',0);
% 
% form = '%s ~ 1 + cond*%s + (1 | pp)';
% fit = {};
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = [2 4]
%     figure();
% 
%     for iB = 1:nBehVars2
%         subplot(r,c,iB);
%         if any(strcmp(factors, behVarNames2{iB}))
%             set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%         end
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             % split by AccFac
% 
%             data = groupMeans(cppVarsByCond.(cppVarNames{iV}), 3, behDataByCond.(behVarNames2{iB})); %[pp acc fac tr]
%         
%             
%             h = errorBarPlot(data, 'type','bar');
%             xticklabels(labels.cond);
%             ylabel(cppVarNames{iV});
%             if iB==1 || iB==nBehVars2; legend(h, labels.(behVarNames2{iB}),'Location','Best'); end
% 
%         
%         else % cond plot
%             
%             [~,~,~,~,h] = conditionalPlot(permute(behDataByCond.(behVarNames2{iB}),[3,1,2]), permute(cppVarsByCond.(cppVarNames{iV}),[3 1 2]));
%             xlabel(behVarNames2{iB});
%             ylabel(cppVarNames{iV});
%             legend(h(2:2:end), labels.cond,'Location','Best');
% 
%         end
% 
%         title(behVarNames2{iB});
% 
%         % do a regression and put star
%         fit{iB} = fitglme(regTab, sprintf(form, cppVarNames{iV}, behVarNames2{iB}));
%         for j = 1:3
%             if fit{iB}.Coefficients.pValue(j+1) < alpha
%                 t = text(0.5, (1-j/10), '*','Units','Normalized','FontSize',16,'Color',[1 2 3]==j);
%             end
%         end
%     end
%     cppBehRegs.behCond.(cppVarNames{iV}) = StatsTableFromLMECells(fit,behVarNames2);
% 
% end
% 
% %% look at post-dec things split by factors + cond + acc
% 
% warning('OFF','NANCAT:emptyCells')
% behDataByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), behDataByCond, 'UniformOutput',0);
% 
% cppVarsByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), cppVarsByCond,'UniformOutput',0);
% warning('ON', 'NANCAT:emptyCells');
% 
% textCols = {'r','b','g','m','c','y','k'};
% 
% form = '%s ~ 1 + acc*cond*%s + (1 | pp)';
% fit = {};
%         
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = [2 4]
%     figure();
%     colOrder = get(gca,'ColorOrder');
%     colOrder = colOrder(1:6,:);
%     
%     for iB = 1:nBehVars2
%         subplot(r,c,iB);
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             if any(strcmp(factors, behVarNames2{iB}))
%                 set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%             else
%                 set(gca,'ColorOrder',colOrder(1:length(labels.(behVarNames2{iB})),:),'nextplot','replacechildren');
%             end
%             % split by AccFac
% 
%             data = groupMeans(cppVarsByCondAcc.(cppVarNames{iV}), 4, behDataByCondAcc.(behVarNames2{iB})); %[pp cond acc fac]
%             data = cat(2, sq(data(:,1,:,:)), sq(data(:,2,:,:)));
% 
%             h = errorBarPlot(data, 'type','bar');
% %             [h(1:length(h)/2).FaceAlpha] = deal(.5);
% 
% %             xticklabels(labels.cond);
%             xticklabels(col(strcat(repmat(labels.cond',1,2), repmat(labels.acc,2,1))'));
%             ylabel(cppVarNames{iV});
%             legend(h, labels.(behVarNames2{iB}),'Location','Best');
% %             if iB==1 || iB==nBehVars2; legend([h(length(h)/2+1:end) h(1)], [labels.(behVarNames2{iB}) 'error...'],'Location','Best'); end
% 
%         
%         else % cond plot
%             
%             [~,~,~,~,h] = conditionalPlot(reshape(permute(behDataByCondAcc.(behVarNames2{iB}),[4,1,2,3]),[],30,4), reshape(permute(cppVarsByCondAcc.(cppVarNames{iV}),[4,1,2,3]),[],30,4));
%             xlabel(behVarNames2{iB});
%             ylabel(cppVarNames{iV});
%             if iB==2; legend(h(2:2:end), col(strcat(repmat(labels.cond,2,1), repmat(labels.acc',1,2))),'Location','Best'); end
% 
%         end
% 
%         title(behVarNames2{iB});
% 
%         % do a regression and put star
%         if ~strcmp(behVarNames2{iB},'acc')
%             fit{iB} = fitglme(regTab, sprintf(form, cppVarNames{iV}, behVarNames2{iB}));
%             for j = 1:length(fit{iB}.CoefficientNames)-1
%                 if fit{iB}.Coefficients.pValue(j+1) < alpha
%                     t = text(0.5, (1-j/10), '*','Units','Normalized','FontSize',16,'Color',textCols{j});
%                 end
%             end
% 
%         end
%     end
%     cppBehRegs.behCondAcc.(cppVarNames{iV}) = StatsTableFromLMECells(fit,behVarNames2);
% 
% end
% 
% 
%% do regressions across each time step - all factors?

% do some times
times = -500:30:1000;
% times = eeg.respTimes(isBetween(eeg.respTimes, [-500 1000]));
formula = 'amplWin ~ 1 + cond*confInR1+ (1 | pp)';
% run one reg to get number of terms
regTab.amplWin = nanzscore(col(cpp(:,1,:))); % use mean
fit = fitglme(regTab(regTab.cond>=-10,:), formula);
stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);
for i = 1:length(times)-1
    % get resplockcumul at t
    % add to regTab
    % regress

%     regTab.amplWin = nanzscore(col(cpp(:,find(eeg.respTimes<=times(i),1,'last'),:))); % use single point

    regTab.amplWin = nanzscore(col(nanmean(cppFilt(:,isBetween(eeg.respTimes, times([i i+1])),:),2))); % use mean


    if sum(~isnan(regTab.amplWin)) > 100
        fit = fitglme(regTab(regTab.cond>=-10,:), formula);
        stats(:,i,1) = fit.Coefficients.Estimate(2:end);
        stats(:,i,2) = fit.Coefficients.pValue(2:end);
    end
end

% alpha = .05 / (length(times)-1); % bonferroni for times - seems too strict, removes effects that survive permutation testing below
% aboveAlpha = stats(:,:,2) > alpha;
% stats(repmat(aboveAlpha,1,1,2)) = NaN;

figure();
subplot(2,1,1)
imagep(stats(:,:,2),fit.CoefficientNames(2:end),times(2:end));
xlabel('time')
ylabel('p-value')

subplot(2,1,2); 
imagesc(stats(:,:,1), [-1 1].*max(abs(stats(:,:,1)),[],'all'));
set(gca, 'ColorMap', crameri('vik'));
set(gca,'YTick',1:20,'YTickLabels',fit.CoefficientNames(2:end), 'XTick', 1:length(times)-1, 'XTickLabels', times(2:end));
xlabel('time')
ylabel('beta')
colorbar;

figure();
isSig = find(any(stats(:,:,2) < .05, 2));
h = plot(times(2:end), stats(isSig,:,1)');
yline(0,':k'); xline(0, ':k');
hold on;
yVal = [min(ylim) abs(diff(ylim))]; 
for j = 1:length(isSig)
    pbar(stats(isSig(j),:,2), 'xVals',times(2:end), 'yVal',yVal(1),'plotargs',{'Color',h(j).Color,'LineWidth',3});
    yVal(1) = yVal(1) - yVal(2)/20; 
end
ylabel('beta');
xlabel('time');
legend(h, fit.CoefficientNames(isSig+1),'Location','Best');

% %% how does that compare to perm testing each one?
% 
% clusterArgs = {'cluster',1,'clustermethod','mass','nperms',5000};
% % behData.hiCert = behData.certainty;
% % behData.hiCert(behData.hiCert==2) = NaN;
% 
% factor = 'CoM';%'CoM'; % must be 2 levels only
%     
% p = [];
% [~,p] = permutationOLS(diff(groupMeans(cpp,3,repmat(permute(behData.cond,[1,3,2]),1,1025,1)),[],3),[],[],[],clusterArgs{:});
% [~,p(2,:)] = permutationOLS(diff(groupMeans(cpp,3,repmat(permute(behData.acc,[1,3,2]),1,1025,1)),[],3),[],[],[],clusterArgs{:});
% [~,p(3,:)] = permutationOLS(diff(groupMeans(cpp,3,repmat(permute(behData.(factor),[1,3,2]),1,1025,1)),[],3),[],[],[],clusterArgs{:});
% 
% % inter
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.cond + 2*behData.acc,[1,3,2]),1,1025,1)),30,[],2,2);
% y = diff(diff(y,[],4),[],3);
% [~,p(4,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% 
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.cond + 2*behData.(factor),[1,3,2]),1,1025,1)),30,[],2,2);
% y = diff(diff(y,[],4),[],3);
% [~,p(5,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.acc + 2*behData.(factor),[1,3,2]),1,1025,1)),30,[],2,2);
% y = diff(diff(y,[],4),[],3);
% [~,p(6,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% 
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.cond + 2*behData.acc + 4*behData.(factor),[1,3,2]),1,1025,1)),30,[],2,2,2);
% y = diff(diff(diff(y,[],5),[],4),[],3);
% [~,p(7,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% 
% %
% terms = {'cond','acc',factor,'cond:acc',['cond:' factor],['acc:' factor], ['cond:acc:' factor]};
% figure();
% imagep(p, terms);
% xticks(1:100:1000);
% xticklabels(round(eeg.respTimes(xticks)));


%% 1-sample t-tests on post-choice slopes

% Given how subtle the build-up is now for the high confidence continued 
% evidence trace it's probably worth doing a one-sample t-test on its slope 
% to show that it's reliably positive

varsByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),cppVars,'UniformOutput',0);
behDataByCond = structfun(@(x) groupMeans(x,2,behData.cond,'dim'), behData, 'UniformOutput',0);

varsByCondCert = structfun(@(x) groupMeans(x,3,behDataByCond.certainty,'dim'), varsByCond, 'UniformOutput',0);

p = zeros(2,3);
for iC = 1:2
    for iCert = 1:3
        [~,p(iC,iCert)] = ttest(nanmean(varsByCondCert.slopePost(:,iC,iCert,:),4));
    end
end

%% plot means in diff windows

% behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]

names.confInR1 = 'confidence-in-initial-choice';
names.certainty = 'confidence-in-final-choice';
names.conf3 = 'confidence-in-initial-choice: binned';
names.CoM = 'change-of-mind';
names.acc = 'initial accuracy';

labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};
labels.conf3 = {'certain/probably CoM', 'maybe CoM/no-CoM', 'probably/certain no-CoM'};
labels.certainty = {'maybe CoM/no-CoM', 'probably CoM/no-CoM', 'certain CoM/no-CoM'};

win = [200 350]; %ms after resp
nW = size(win,1);
figure();


iFs = [4];
for j = 1:length(iFs)
    iF=iFs(j);
    factor = factors{iF};
    
    for i = 1:nW
        newAmpl = sq(nanmean(cpp(:,isBetween(eeg.respTimes,win(i,:)),:),2));
        amplByCond = groupMeans(newAmpl,2,behData.cond,'dim');
        
        amplByCondFac = groupMeans(amplByCond,3,behDataByCond.(factor));
        
        subplot(length(iFs),nW,(i-1)*length(iFs)+j);
    
        set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
        h = errorBarPlot(amplByCondFac,'type','bar');
        xticks(1:2); xticklabels(labels.cond);
        ylabel(sprintf('mean Pe amplitude%d:%dms', win(i,1), win(i,2)));
        title(names.(factors{iF}));
    
        regTab.f = nanzscore(col(newAmpl));
        for iC = 1:2
            f = fitglme(regTab((regTab.cond>0) == (iC-1),:), sprintf('f ~ 1 + %s + (1 | pp)', factor));
            if f.Coefficients.pValue(2)<.05
                hold on;
                text(iC, max(ylim)/2, '*','FontSize',20);
            end
        end
        
    
        if i==1; legend(h, labels.(factor),'Location','Best'); end
    
    end

end

%% 
regTab.preEv=regTab.f;
% now do preresp-baseline, store as regTab.f
% regTab.f = nanzscore(col(newAmpl));
r2=[regTab;regTab];
% newAmpl = sq(nanmean(dataResp.cppFiltCond(:, isBetween(eeg.respTimes, amplWindows(2,:)),iC,:),2)); %[pp tr]
% newAmpl2 = sq(nanmean(cppFiltCond(:, isBetween(eeg.respTimes, amplWindows(2,:)),iC,:),2)); %[pp tr]
% r2.f = nanzscore([col(newAmpl); col(newAmpl2)]);
r2.Pe = nanzscore([regTab.preEv; regTab.f]);
r2.bs = nanzscore([zeros(numel(newAmpl),1); ones(numel(newAmpl),1)]);
fitglme(r2, sprintf('Pe ~ 1 + bs*%s + (1 | bs)', factors{iF}))

