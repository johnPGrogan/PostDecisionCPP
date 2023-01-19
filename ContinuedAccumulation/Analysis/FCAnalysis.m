% FCAnalysis
% run GetFC first

if ~exist('filter_tf','file'); eeglab; close all; end

%%

clear; clc; close all;

%%%%%% set options
useCSD = 1;
excludeBadPps = 1; % remove pps with <640 good trials?
excludeTooFew = 1; % remove pps with <20 per conf3
excludeByRT = 1; % remove trials outside [100 1500] ms
doFilt = 1; % whether to plot with loPass filtered data
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty

if useCSD
    fcName = 'GetFCCSD.mat';
    fcTopoName = 'GetCPPTopoCSD.mat'; % topo uses same data as cpp
    mapLims = [-15 15];
else
    fcName = 'GetFC.mat';
    fcTopoName = 'GetCPPTopo.mat';
    mapLims = [-3 3];
end

% load stuff
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));
load(fullfile(outFolder, fcName),'fc','fcStim','eeg','stimWin','chans5','fcFilt','fcStimFilt');
load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');
load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels','behTab');
load(fullfile('./Saves/',fcTopoName),'cppMeanTopo','wins','cppRTTopo','cppTopoAmplWindow');


%% remove stuff - flagged, outliers, bad pps, too few trials per bin

toRemove = sq(any(isnan(fc),2)); % flagged trials or cpp outliers

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;     
end

if excludeCoMFromCert
    %%%% set any certainty with CoM to NaN - will be ignored everywhere
    behData.certainty(behData.CoM==1) = NaN;
end

if excludeByRT
    rtLims = [100 1540];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

if excludeTooFew

    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.cond);

    toRemove(toRemove2) = 1;

end

% actually remove from data now
fc(repmat(permute(toRemove,[1,3,2]),1,size(fc,2),1)) = NaN;
fcFilt(repmat(permute(toRemove,[1,3,2]),1,size(fcFilt,2),1)) = NaN;
fcStim(repmat(permute(toRemove,[1,3,2]),1,size(fcStim,2),1)) = NaN;
fcStimFilt(repmat(permute(toRemove,[1,3,2]),1,size(fcStimFilt,2),1)) = NaN;
cppMeanTopo(repmat(toRemove,1,1,128,size(cppMeanTopo,4))) = NaN;
cppRTTopo(repmat(toRemove,1,1,128)) = NaN;
cppTopoAmplWindow(repmat(toRemove,1,1,128,2)) = NaN;


%% define windows and labels

amplWindows = [-140 -50; 700 1000]; % ms
slopeWindows = [-310 -100; 400 700];

cmap = crameri('vik');

factors = {'acc','certainty','CoM','confInR1','conf3'};
nF = length(factors);
labels.diffNames = {'Error - Correct' ,'Low - High',  'CoM - NoCoM','Low - High','Low - High'};
diffFlips = [-1 -1 1 -1 -1]; % invert some to match order

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2]; };

stimWin = [-200 1200];
respWin = [-500 1000];

%% filter out 30Hz SSVEP?
if ~doFilt

    fcFilt = fc;
    fcStimFilt = fcStim;

end

%% get means + slopes

fcVars.ampl = sq(nanmean(fc(:,isBetween(eeg.respTimes,amplWindows(1,:)),:),2));
fcVars.amplPost = sq(nanmean(fc(:,isBetween(eeg.respTimes,amplWindows(2,:)),:),2));

fcVars.slope = FitCPPSlope(fc, slopeWindows(1,:), eeg.respTimes); % get slopes
fcVars.slopePost = FitCPPSlope(fc, slopeWindows(2,:), eeg.respTimes); % get slopes

fcVarNames = fieldnames(fcVars);
nFcVars = 2;%length(fcVarNames);

varsByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),fcVars,'UniformOutput',0);


%% does single-trial lme have better sensitivity?

if excludeCoMFromCert
    behTab.certainty(behTab.CoM>0) = NaN;
end

% remove trials
behTab(col(toRemove),:) = array2table(NaN(sum(toRemove,'all'),width(behTab)));

% vars
fcTab = struct2table(structfun(@(x) col(x), fcVars, 'UniformOutput',0));
% matVars = nanzscore(matVars);

regTab = horzcat(behTab, fcTab);
regNames = regTab.Properties.VariableNames;

isLogistic = cellRegexpi(regNames, 'Logistic')>0;
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));
% 
% glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

nF=3;


%% paper statistics (in supplementary)

% this will apply regression to a cell array
fitglmeCell = @(x) fitglme(regTab, x);
% for eachcond sep, i=0 | 1
fitglmeCellSep = @(x, i) fitglme(regTab(regTab.cond>0 == i,:), x); 



if excludeCoMFromCert
    % first-order/pre choice CPP. excludeComFromCert == 1
    paperFormulaePre = {
         'ampl ~ 1 + acc*cond + (1 | pp)'; % cpp ampl
         'ampl ~ 1 + certainty*cond + (1 | pp)';
         'ampl ~ 1 + acc*certainty*cond+ (1 | pp)'; % 3-way 
         };
    paperPreFits = cellfun(fitglmeCell, paperFormulaePre,'UniformOutput',0);
    
    % run these separate LME in each ev cond separately
    % excludeCoMFromCert
    paperFormulaePreSep = {'ampl ~ 1 + certainty + (1 | pp)';};
    paperPreSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaePreSep,2,1), {0;1}, 'UniformOutput',0);


else
    
    % post-choice/2nd order CPP
    % don't excludeCoMFromCert
    paperFormulaePost = {'amplPost ~ 1 + cond*acc + (1 | pp)';
        'amplPost ~ 1 + cond*confInR1 + (1 | pp)';
        'amplPost ~ 1 + cond*conf3 + (1 | pp)';
        'amplPost ~ 1 + cond*certainty+ (1 | pp)';
        };
    paperPostFits = cellfun(fitglmeCell, paperFormulaePost,'UniformOutput',0);
    
    
    
    % don't excludeCoMFRomCert
    paperFormulaePostSep = {'amplPost ~ 1 + acc + (1 | pp)';
        'amplPost ~ 1 + confInR1 + (1 | pp)';
        'amplPost ~ 1 + certainty + (1 | pp)'};
    paperPostSepFits = cellfun(fitglmeCellSep, repmat(paperFormulaePostSep,1,2), {0 1; 0 1; 0 1}, 'UniformOutput',0);
    
    
end

%% save now

saveOpts = {'Volt','CSD'; '', 'ExclCoMFromCert'};
saveName = sprintf('FCAnalysis_%s_%s.mat', saveOpts{1,useCSD+1}, saveOpts{2, excludeCoMFromCert+1});

save(fullfile(outFolder, saveName));

%%%% don't run the figures below, instead run PlotPaperFigures.m

%% plot per pp

figure();
[r,c] = GetSubPlotShape(fileInfo.nPP);
for iPP = 1:fileInfo.nPP
    
    subplot(r,c,iPP);
    errorBarPlot(permute(fcFilt(iPP,:,:),[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes);
    xlines(slopeWindows);
    xline(0,':k');
    yline(0, ':k');
    
    title(sprintf('%s: %s', fileInfo.ppID{iPP}));%, fcChans{iPP}));
    
end
% makeSubplotScalesEqual(r,c,1:fileInfo.nPP);


% %% all on one - includes tobeexcluded pps
% 
% figure();
% c = get(gca,'ColorOrder');
% for i = 1:2
%     hold on; 
%     h{i} = errorBarPlot(permute(fc(ppsToExclude==(i-1),:,:),[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes,'color',c(i,:));
% end
% 
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% 
% title('grand mean');
% legend([h{1}{1,2}, h{2}{1,2}], {'Include','Exclude'},'Location','Best');

%% plot all included pps together
% 
figure();
% h = errorBarPlot(permute(fcFilt,[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes);
plot(eeg.respTimes, sq(nanmean(fcFilt,3)),'LineWidth',2);
xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;

title('grand mean');



%% average over pps

figure();

errorBarPlot(nanmean(fcFilt,3),'area',1,'xaxisvalues',eeg.respTimes);

xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;

title('grand mean');
% 
%% average over pps - split by cond

figure();

h = errorBarPlot(groupMeans(fcFilt,3, repmat(permute(behData.cond,[1,3,2]),1,length(eeg.respTimes))),'area',1,'xaxisvalues',eeg.respTimes);
xline(0,':k');
yline(0, ':k');
% xlines(amplWindows, '-k');
% xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;
legend([h{:,1}], labels.cond,'Location','Best');

title('grand mean');


%% plot FC by RT bin
% 
% % 2 or 3 speed bins? uncomment one below to choose
% binNames = {'fast','slow'};
% % binNames = {'fast','medium','slow'};
% p = linspace(0,100, length(binNames)+1);
% p(1) = [];
% prc = prctile(behData.RT, p,2);
% binInds = sum(behData.RT > permute(prc,[1,3,2]),3);
% 
% fcByRT = groupMeans(fcFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(fc,2),1),'dim'); %[npp t bins tr]
% 
% figure();
% subplot(1,2,1);
% fcStimByRT = groupMeans(fcStimFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(fcStim,2),1),'dim'); %[npp t bins tr]
% h = errorBarPlot(nanmean(fcStimByRT,4),'area',1,'xaxisvalues',eeg.stimTimes);
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from stim (ms)');
% box off;
% title('evidence locked');
% 
% subplot(1,2,2)
% h = errorBarPlot(nanmean(fcByRT,4),'area',1,'xaxisvalues',eeg.respTimes);
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% title('response locked');
% legend([h{:,1}], binNames,'Location','Best');
% 
% %% plot fc by RT bin, and acc
% 
% % 2 or 3 speed bins? uncomment one below to choose
% accByRT = repmat(permute(groupMeans(behData.acc, 2, binInds,'dim'),[1 4 2 3]),1,size(fc,2),1,1); %[npp t bins acc tr]
% fcStimByRTAcc = groupMeans(fcStimByRT, 4, accByRT(:,1:size(fcStim,2),:,:)); %[npp t bins tr]
% fcByRTAcc = groupMeans(fcByRT,4,accByRT); %[pp t rt acc]
% 
% figure();
% for i = 1:2
%     subplot(2,2,i*2-1);
%     
%     h = errorBarPlot(fcStimByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.stimTimes);
%     xline(0,':k');
%     yline(0, ':k');
%     ylabel('\muV'); xlabel('time from stim (ms)');
%     box off;
%     title(['evidence locked, ' labels.acc{i}]);
%     
%     subplot(2,2,i*2)
%     h = errorBarPlot(fcByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.respTimes);
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
%     h = errorBarPlot(permute(fcByRT(iPP,:,:,:),[4,2,3,1]),'area',1,'xaxisvalues',eeg.respTimes);
%     xline(0,':k');
%     yline(0, ':k');
%     ylabel('\muV'); xlabel('time from resp 1 (ms)');
%     box off;
%     title(iPP);
% end    
% legend([h{:,1}], binNames,'Location','Best');
% 

%% make fig3.5

% split [stimlocked fc] by [conf; acc; CoM], do a topoplot of 600-800ms
% post decision

cppMeanTopo600 = nanmean(cppMeanTopo(:,:,:,18:19),4); % average over 600-800ms
ylims = [-1.5 4; -1.4 5];

nF=3;
figure();
for i = 1:nF
    subplot(nF,3,i*3-2);
    
    % split stimlocked by cond
    stimlocked = groupMeans(fcStimFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(fcStim,2))); %[pp t cond]
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
        % split fc by cond
    resplocked = groupMeans(fcFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(fc,2))); %[pp t cond]
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



% split fcFilt by Cond
fcFiltCond = groupMeans(fcFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(fc,2)),'dim'); %[pp t cond tr]
figure();
for i = 1:nF
%     figure();
   
    factorByCond = repmat(permute(groupMeans(behData.(factors{i}),2,behData.cond,'dim'),[1 4 2 3]),1,length(eeg.respTimes),1,1); % to match fcFiltCond

    % split fc by cond
    resplocked = groupMeans(fcFiltCond,4,factorByCond); %[pp t cond factor]
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
            ylabel('FCP \muV'); 
            legend(h, [labels.(factors{i}) labels.diffNames{i}],'Location','NorthEast');
%             legend([h{:,1}], [labels.(factors{i}) labels.diffNames{i}],'Location','NorthEast');
        end

        xlim(respWin); %ylim(ylims(2,:));
    end    
         
end
makeSubplotScalesEqual(nF,2)

%% can also perm plot stuff

fcFiltCond = groupMeans(fcFilt,3,repmat(permute(behData.cond,[1,3,2]),1,size(fc,2)),'dim'); %[pp t cond tr]

fig = figure();
% 
% for i = 1:5
% %     subplotInds = [5 2 i*2-1; 5 2 i*2]; % column per fig type
%     figure(fig);
%     subplotInds = [3,2,i]; % column per fig type
%     myPermPlots(fcFiltCond, factors{i}, eeg.respTimes, behData, diffFlips(i), cols{i}, respWin, labels, labels.diffNames(i),subplotInds);
% end

clusterArgs = {'cluster',1,'clustermethod','mass','nperms',10000};
n = {'acc','certainty'};%,'conf3','acc','CoM','confInR1'};
[r,c] = GetSubPlotShape(length(n));
for i = 1:length(n)
    j = find(strcmp(factors, n(i)));
%     subplotInds = [5 2 i*2-1; 5 2 i*2]; % column per fig type
    figure(fig);
    subplotInds = [r,c,i]; % column per fig type
    myPermPlots(fcFiltCond, factors{j}, eeg.respTimes, behData, diffFlips(j), cols{j}, respWin, labels, labels.diffNames(j),subplotInds,[],0,clusterArgs);
end

%% plot acc*fac for conditions separately
clusterArgs = {'cluster',1,'clustermethod','mass','nperms',5000};

respWin = [-500 1000];
tInds = isBetween(eeg.respTimes, respWin);

% split by acc too
fcFiltCond = groupMeans(fcFilt(:,tInds,:),3,repmat(permute(behData.cond,[1,3,2]),1,sum(tInds)),'dim'); %[pp t cond tr]
accByCond = repmat(permute(groupMeans(behData.acc,2,behData.cond,'dim'),[1 4 2 3]),1,size(fcFiltCond,2),1,1); % to match fcFiltCond
fcFiltCondAcc = groupMeans(fcFiltCond,4,accByCond,'dim'); %[pp time cond acc tr] 


lineStyles = {':','-'};
for iF = 2
    figure();
    factor = factors{iF};
    facByCond  = repmat(permute(groupMeans(behData.(factor),2,behData.cond,'dim'),[1 4 2 3]),1,size(fcFiltCond,2),1,1); % to match fcFiltCond
    facByCondAcc = groupMeans(facByCond,4,accByCond,'dim');
    fcFiltCondAccFac = groupMeans(fcFiltCondAcc, 5, facByCondAcc); %[pp time cond acc fac]

    % permute here to swap cond/acc
%     fcFiltCondAccFac = permute(fcFiltCondAccFac,[1,2,4,3,5]);

    % get diff waves - 
%     fcFiltCondAccFac(:,:,:,:,end+1) = diff(fcFiltCondAccFac(:,:,:,:,[1 end]),[],5) .* diffFlips(iF);
    for j = 1:2 % acc
%         subplot(1,2,j);
%         set(gca,'ColorOrder', cols{iF}, 'nextplot','replacechildren');

        for k = 1:2    
             subplot(2,2,(j-1)*2+k);
            set(gca,'ColorOrder', cols{iF}, 'nextplot','replacechildren');
                    hold on;

%             h = plot(eeg.respTimes(tInds), sq(nanmean(fcFiltCondAccFac(:,:,j,k,:),1)), 'LineWidth',2,'LineStyle',lineStyles{2});
            h = errorBarPlot(sq(fcFiltCondAccFac(:,:,j,k,:)),'area',1,'xaxisvalues',eeg.respTimes(tInds),'plotargs',{'LineStyle',lineStyles{2}});
        
        
            title([labels.cond{j}, ',' labels.acc{k}]);
            xline(0);yline(0);
            ylabel(factor); 
    
%             % perm testing - get corr-Err, Com-none, interaction
%             clear diffWaves;
%             diffWaves(:,:,1) = sq(diff(nanmean(fcFiltCondAccFac(:,:,j,:,1:2),5),[],4)); % av over fac, effect of acc
%             diffWaves(:,:,2) = sq(diff(nanmean(fcFiltCondAccFac(:,:,j,:,1:2),4),[],5)) .* diffFlips(iF); % av over acc, effect of fac
%             diffWaves(:,:,3) = sq(diff(diff(fcFiltCondAccFac(:,:,j,:,1:2),[],5),[],4)) .* diffFlips(iF); % diff over fac, then acc
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

%     legend(h, labels.(factors{iF}),'Location','Best')
    legend([h{:,1}], labels.(factors{iF}),'Location','Best')
end

% makeSubplotScalesEqual(1,2);

%% split by conf*acc*com*cert?

iF = find(strcmp(factors, 'confInR1')); % confInR1, will split into cert*CoM later

facByCond  = repmat(permute(groupMeans(behData.(factors{iF}),2,behData.cond,'dim'),[1 4 2 3]),1,size(fcFiltCond,2),1,1); % to match fcFiltCond
facByCondAcc = groupMeans(facByCond,4,accByCond,'dim');
fcFiltCondAccFac = groupMeans(fcFiltCondAcc, 5, facByCondAcc); %[pp time cond acc fac]
% need to change order of conf6
fcFiltCondAccFac = fcFiltCondAccFac(:,:,:,:,[4 5 6 3 2 1]);
fcFiltCondAccCertCoM = reshape(fcFiltCondAccFac,30,[],2,2,3,2); %[pp time cond acc cert com


iF = 2; % reset to cert
figure();
for i = 2
    for j = 1:2
        for k = 1:2
            subplot(2,2,(j-1)*2+k);
            set(gca,'ColorOrder',cols{iF}(1:end-1,:),'nextplot','replacechildren');
            hold on;

            h = errorBarPlot(sq(fcFiltCondAccCertCoM(:,:,i,j,:,k)), 'area',1,'xaxisvalues',eeg.respTimes(tInds),'plotargs',{'LineStyle',lineStyles{k}});
%             h = plot(eeg.respTimes(tInds), sq(nanmean(fcFiltCondAccCertCoM(:,:,i,j,:,k),1)), 'LineWidth',2,'LineStyle',lineStyles{k});
        
        
            title([labels.cond{i}, ', ' labels.acc{j}]);
            xline(0);yline(0);
            ylabel(factors{iF}); 
        end
%         % perm testing - get corr-Err, Com-none, interaction
%         clear diffWaves;
%         diffWaves(:,:,1) = sq(diff(nanmean(fcFiltCondAccCertCoM(:,:,j,:,1:2),5),[],4)); % av over fac, effect of acc
%         diffWaves(:,:,2) = sq(diff(nanmean(fcFiltCondAccCertCoM(:,:,j,:,1:2),4),[],5)) .* diffFlips(iF); % av over acc, effect of fac
%         diffWaves(:,:,3) = sq(diff(diff(fcFiltCondAccCertCoM(:,:,j,:,1:2),[],5),[],4)) .* diffFlips(iF); % diff over fac, then acc
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
%     scatterRegress(col(fcVars.(fcVarNames{i})), col(fcVars.(fcVarNames{i+2})), scatterArgs{:});
%     xlabel(fcVarNames{i});
%     ylabel(fcVarNames{i+2});
% 
%     subplot(2,2,i + 2);
%     conditionalPlot(fcVars.(fcVarNames{i})', fcVars.(fcVarNames{i+2})');
%     xlabel(fcVarNames{i});
%     ylabel(fcVarNames{i+2});
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
%     scatterRegress(col(fcVars.(fcVarNames{i*2-1})), col(fcVars.(fcVarNames{i*2})), scatterArgs{:});
%     xlabel(fcVarNames{i*2-1});
%     ylabel(fcVarNames{i*2});
% 
%     subplot(2,2,i + 2);
%     conditionalPlot(fcVars.(fcVarNames{i*2-1})', fcVars.(fcVarNames{i*2})');
%     xlabel(fcVarNames{i*2-1});
%     ylabel(fcVarNames{i*2});
% 
% end
%% split by cond and then cert/acc/com



figure();
for iF = 1:nF
    % split by cond %[pp cond 3 tr]
    factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');

    for iV = 1:nFcVars

        % split by factor
        data = groupMeans(varsByCond.(fcVarNames{iV}), 3, factor); %[pp cond 3 factor]
        
        subplot(nFcVars,nF,(iV-1)*nF+iF)
        h = errorBarPlot(data, 'type','bar');
        set(gca,'XTickLabel',labels.cond);
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(i,:);
        end
        if iV==1
            legend(h, labels.(factors{iF}),'Location','Best'); 
            title(factors{iF});
        end
        if iF==1; ylabel(fcVarNames{iV}); end
        
        
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
                data = groupMeans(fcVars.(fcVarNames{iV*2-(iP-1)}), 2, behData.(factors{iF})); %[pp factor]
            else % include cond
                % split by factor
                data = groupMeans(varsByCond.(fcVarNames{iV*2-(iP-1)}), 3, factor); %[pp cond 3 factor]
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
            if iF==1; ylabel(fcVarNames{iV*2-(iP-1)}); end
            
            
        end
    end    
    makeSubplotScalesEqual(2,nF,1:nF);
    makeSubplotScalesEqual(2,nF,(nF+1):2*nF);

end

%% split acc/err


behDataByCond = structfun(@(x) groupMeans(x, 2, behData.cond,'dim'),behData,'UniformOutput',0); %[pp cond tr]
varsByCondAcc = structfun(@(x) groupMeans(x, 3, behDataByCond.confAcc,'dim'),varsByCond,'UniformOutput',0); %[pp cond acc tr]

% figure();
for iF = 2
    figure();
    factor = groupMeans(behDataByCond.(factors{iF}),3,behDataByCond.confAcc, 'dim');

    for iV = 1:2%nFcVars    
        

        % split by factor
        data = groupMeans(varsByCondAcc.(fcVarNames{iV*2}), 4, factor); %[pp cond acc factor]
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
            ylabel(fcVarNames{iV*2}); 
        
        end 
        makeSubplotScalesEqual(2,2,iV*2-1:iV*2);
    end
    
end


%% split by certainty,com,acc,cond
% massively unequal trials numbers per cell (min mean is <1 int,corr,high,Change, max mean is
% 183)

behDataByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), behDataByCond, 'UniformOutput',0);
fcVarsByCondAcc = structfun(@(x) groupMeans(x,3,behDataByCond.acc,'dim'), varsByCond,'UniformOutput',0);

amplPostByCondAccCertCoM = groupMeans(fcVarsByCondAcc.amplPost, 4, behDataByCondAcc.confInR1); % split by conf6
amplPostByCondAccCertCoM = amplPostByCondAccCertCoM(:,:,:,[4 5 6 3 2 1]); % invert these to make go low-high low-high
amplPostByCondAccCertCoM = reshape(amplPostByCondAccCertCoM,30,2,2,3,2); % split by CoM

% amplPostByCondAccCertCoM = sq(nanmean(fcFiltCondAccCertCoM(:,isBetween(eeg.respTimes(tInds), [800 950]),:,:,:,:),2));

figure();

for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        set(gca,'ColorOrder',cols{2}(1:end-1,:),'nextplot','replacechildren');

        h = errorBarPlot(permute(amplPostByCondAccCertCoM(:,i,j,:,:),[1 5 4 2 3]),'type','bar');

        title([labels.cond{i}, ', ' labels.acc{j}]);
        xlabel(factors{4}); 
        xticks(1:2); xticklabels(labels.CoM);
        ylabel('post-dec FCP ampl');

    end
end

legend(h, labels.certainty, 'Location','Best');
makeSubplotScalesEqual(2,2);

%% now for slopes

slopePostByCondAccCertCoM = groupMeans(fcVarsByCondAcc.slopePost, 4, behDataByCondAcc.confInR1); % split by conf6
slopePostByCondAccCertCoM = slopePostByCondAccCertCoM(:,:,:,[4 5 6 3 2 1]); % invert these to make go low-high low-high
slopePostByCondAccCertCoM = reshape(slopePostByCondAccCertCoM,30,2,2,3,2); % split by CoM

% amplPostByCondAccCertCoM = sq(nanmean(fcFiltCondAccCertCoM(:,isBetween(eeg.respTimes(tInds), [800 950]),:,:,:,:),2));

figure();

for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        set(gca,'ColorOrder',cols{2}(1:end-1,:),'nextplot','replacechildren');

        h = errorBarPlot(permute(slopePostByCondAccCertCoM(:,i,j,:,:),[1 5 4 2 3]),'type','bar');

        title([labels.cond{i}, ', ' labels.acc{j}]);
        xlabel(factors{4}); 
        xticks(1:2); xticklabels(labels.CoM);
        ylabel('post-dec FCP slope');

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
        ylabel(['post-dec FCP ' fcVarNames{i*2}]);

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
    ylabel('post-dec FCP ampl');

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

    for iV = 1:nFcVars
        
        data = groupMeans(varsByCond.(fcVarNames{iV}), 3, factor); %[pp cond 3 factor]
                
        anovas2.(factors{iF}).(fcVarNames{iV}) = rmanova(data, {'pp','cond',factors{iF}},'categorical',2:3, 'DummyVarCoding','effects');
        
    end
    anovas2Tabs.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovas2.(factors{iF}),'UniformOutput',0),'rowNames',anovas2.(factors{iF}).ampl.Term(2:end));
end

%% cont/inter separately

[anovasSep, anovasSepTabs] = deal(struct());
for iF = 1:nF
    factor = groupMeans(behData.(factors{iF}),2,behData.cond, 'dim');

    for iV = 1:nFcVars
        
        data = groupMeans(varsByCond.(fcVarNames{iV}), 3, factor); %[pp cond 3 factor]

        for j = 1:2        
            anovasSep.(labels.cond{j}).(factors{iF}).(fcVarNames{iV}) = rmanova(sq(data(:,j,:)), {'pp',factors{iF}},'categorical',[], 'DummyVarCoding','effects');
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




%% do regressions across each time step - all factors?

% do some times
times = -500:100:1000;
% times = eeg.respTimes(isBetween(eeg.respTimes, [-500 1000]));
% f4 = '%s ~ 1 + cond*acc*certainty*CoM + (1 | pp)';
f4 = '%s ~ 1 + cond*acc*certainty + (1 | pp)';
fit = fitglme(regTab, sprintf(f4, 'amplPost'));

stats = NaN(length(fit.CoefficientNames)-1,length(times),2);
for i = 1:length(times)-1
    % get resplockcumul at t
    % add to regTab
    % regress

%     regTab.amplWin = nanzscore(col(fc(:,find(eeg.respTimes<=times(i),1,'last'),:))); % use single point

%     regTab.amplWin = nanzscore(col(nanmean(fc(:,isBetween(eeg.respTimes, times([i i+1])),:),2))); % use mean

    regTab.amplWin = nanzscore(col(FitCPPSlope(fc, times([i i+1]), eeg.respTimes))); % use slope, still call it mean for simplicity


    if sum(~isnan(regTab.amplWin)) > 100
        fit = fitglme(regTab, sprintf(f4,'amplWin'));
        stats(:,i,1) = fit.Coefficients.Estimate(2:end);
        stats(:,i,2) = fit.Coefficients.pValue(2:end);
    end
end

% alpha = .05 / (length(times)-1); % bonferroni for times - seems too strict, removes effects that survive permutation testing below
% aboveAlpha = stats(:,:,2) > alpha;
% stats(repmat(aboveAlpha,1,1,2)) = NaN;

figure();
imagep(stats(:,:,2),fit.CoefficientNames(2:end),times);
xlabel('time')
ylabel('beta')

