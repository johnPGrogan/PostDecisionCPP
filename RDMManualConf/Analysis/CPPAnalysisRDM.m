% CPPAnalysisRDM
% run GetCPP first

if ~exist('filter_tf','file'); eeglab; close all; end
%%

clear; clc; close all;

%%%%%% set options
useCSD = 1;
excludeBadPps = 1; % remove pps with <640 good trials?
excludeTooFew = 1; % remove pps with <20 per conf3
excludeByRT = 1; % remove trials outside [100 1500] ms
doFilt = 1; % whether to plot with loPass filtered data

if useCSD
    cppName = 'GetCPPCSD.mat';
    cppTopoName = 'GetCPPTopoCSD.mat';
    cppTopoNameRC = 'GetCPPRespCueTopoCSD.mat';
    mapLims = [-20 20];
else
    cppName = 'GetCPP.mat';
    cppTopoName = 'GetCPPTopo.mat';
    cppTopoNameRC = 'GetCPPRespCueTopo.mat';
    mapLims = [-3 3];
end

% load stuff
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));
load(fullfile(outFolder, cppName),'cpp','cppStim','eeg','stimWin','preRespWindow','cppChans','cppFilt','cppStimFilt','cppRespCue','cppRespCueFilt','chans5');
load(fullfile(outFolder, 'FlagArtefacts.mat'),'isFlagged','ppsToExclude');
load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels');
load(fullfile(outFolder, cppTopoName),'cppMeanTopo','wins','cppRTTopo','cppTopoAmplWindow');
load(fullfile(outFolder, cppTopoNameRC),'respCueMeanTopo');

behData = rmfield(behData, {'dataMat','dataMatNames'}); % remove these

%% remove stuff - flagged, outliers, bad pps, too few trials per bin

toRemove = sq(any(isnan(cpp),2)); % flagged trials already out

if excludeBadPps
    toRemove(ppsToExclude,:) = 1;     
end

if excludeByRT
    rtLims = [100 1500];
    toRemove(~isBetween(behData.RT, rtLims)) = 1;
end

if excludeTooFew
    % only remove those trials, not the entire person
    behData.certainty(toRemove) = NaN;
    toRemove2 = GetCellsWithLowN(10, behData.certainty, behData.stimDur);

    toRemove(toRemove2) = 1;

%     % remove those people entirely, not just for that cond
%     toRemove(any(CountUnique(groupMeans(behData.certainty,2,behData.stimDur,'dim'),3)<10,[2 3]),:) = 1;

end

% actually remove from data now
cpp(repmat(permute(toRemove,[1,3,2]),1,size(cpp,2),1)) = NaN;
cppFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppFilt,2),1)) = NaN;
cppStim(repmat(permute(toRemove,[1,3,2]),1,size(cppStim,2),1)) = NaN;
cppStimFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppStimFilt,2),1)) = NaN;
cppRespCue(repmat(permute(toRemove,[1,3,2]),1,size(cppRespCue,2),1)) = NaN;
cppRespCueFilt(repmat(permute(toRemove,[1,3,2]),1,size(cppRespCueFilt,2),1)) = NaN;
cppMeanTopo(repmat(toRemove,1,1,128,size(cppMeanTopo,4))) = NaN;
cppRTTopo(repmat(toRemove,1,1,128)) = NaN;
respCueMeanTopo(repmat(toRemove,1,1,128,size(respCueMeanTopo,4))) = NaN;
cppTopoAmplWindow(repmat(toRemove,1,1,128)) = NaN;

% remove from behData too
behNames = fieldnames(behData);
for i = 1:length(behNames)
    behData.(behNames{i})(toRemove) = NaN;
end


%% define windows and labels

amplWindows = [-150 -50;]; % ms, resplocked
slopeWindows = [-500 -200;];

cmap = crameri('vik');

factors = {'acc','certainty','confInCorr','confResp','stimDur'};
nF = length(factors);
labels.diffNames = {'Error - Correct' ,'Low - High', 'Low - High' ,'Low - High','Slow - Fast',};
diffFlips = [-1 -1 -1 -1 -1]; % invert some to match order

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2]
        };

nF=2;

%% filter out 30Hz SSVEP? - only for plotting
if ~doFilt

    cppFilt = cpp;
    cppStimFilt = cppStim;

end

%% get means + slopes

cppVars.ampl = sq(nanmean(cpp(:,isBetween(eeg.respTimes,amplWindows(1,:)),:),2));

cppVars.slope = FitCPPSlope(cpp, slopeWindows(1,:), eeg.respTimes); % get slopes

cppVarNames = fieldnames(cppVars);
nCppVars = length(cppVarNames);

%% does single-trial lme have better sensitivity?
behVarNames = {'pp','stimDur','corrLR','respLR','acc','RT','confResp','certainty','confInCorr'};
nBehVars = length(behVarNames);
behNames = fieldnames(behData); % names of all behData

behTabNames = behVarNames;
behTab = struct2table(structfun(@col, rmfield(behData, behNames(~ismember(behNames, behTabNames))), 'UniformOutput',0));

nF = length(factors); % do all

% vars
cppTab = struct2table(structfun(@(x) col(x), cppVars, 'UniformOutput',0));
% matVars = nanzscore(matVars);

regTab = horzcat(behTab, cppTab);
regNames = regTab.Properties.VariableNames;

% re-order names 
[~,b] = ismember([behTabNames, cppVarNames'],regNames);
regTab = regTab(:, b);

regTab.accLogistic = regTab.acc; % only if DV
regTab.RTLog = col(log(behData.RT));

regNames = regTab.Properties.VariableNames;

isLogistic = cellRegexpi(regNames, 'Logistic')>0;
regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV


%% statistics in the paper:

paperFormulae = {'accLogistic ~ 1 + stimDur + (1 | pp)';
     'certainty ~ 1 + stimDur + (1 | pp)';
     'RTLog ~ 1 + stimDur + (1 | pp)';
     'ampl ~ 1 + acc*stimDur + (1 | pp)'; % cpp ampl
     'ampl ~ 1 + certainty*stimDur + (1 | pp)';
     'ampl ~ 1 + acc*certainty*stimDur + (1 | pp)'; % 3-way 
     };
fitglmeCell = @(x) fitglme(regTab, x, glmeArgs{ ~isempty(regexp(x,'Logistic'))+ 1}{:});
paperFits = cellfun(fitglmeCell, paperFormulae, 'UniformOutput',0);

% post-hoc tests: make certainty a categorical, re-run and use coefTest to
% compare levels of Certainty
regTab.certCat = categorical(regTab.certainty);
f1 = fitglme(regTab, 'ampl ~ 1 + certCat*stimDur + (1 | pp)','DummyVarCoding','effects');
certContrasts = [0 0 1 0 0 0; % low vs high
    0 0 0 1 0 0; % med vs high
    0 0 1 -1 0 0];  % low vs med
for i = 1:size(certContrasts,1)
    [postHocStats(i,1), postHocStats(i,2), postHocStats(i,3), postHocStats(i,4)] = coefTest(f1, certContrasts(i,:));
end


%% save now

saveOpts = {'Volt','CSD'};
saveName = sprintf('CPPAnalysis_%s.mat', saveOpts{1,useCSD+1});

save(fullfile(outFolder, saveName));

% use ../../ContinuedAccumulation/Analysis/PlotPaperFigures.m to plot

%%
% 
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
% %% all on one - includes tobeexcluded pps
% 
% figure();
% c = get(gca,'ColorOrder');
% for i = 1:2
%     hold on; 
%     h{i} = errorBarPlot(permute(cpp(ppsToExclude==(i-1),:,:),[3,2,1]),'area',1,'xaxisvalues',eeg.respTimes,'color',c(i,:));
% end
% 
% xline(0,':k');
% yline(0, ':k');
% ylabel('\muV'); xlabel('time from resp 1 (ms)');
% box off;
% 
% title('grand mean');
% legend([h{1}{1,2}, h{2}{1,2}], {'Include','Exclude'},'Location','Best');
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

%% average over pps

figure();

errorBarPlot(nanmean(cppFilt,3),'area',1,'xaxisvalues',eeg.respTimes);

xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;

title('grand mean');

%% average over pps - split by stimDur

figure();
subplot(2,2,1);
h = errorBarPlot(groupMeans(cppStimFilt,3, repmat(permute(behData.stimDur,[1,3,2]),1,length(eeg.stimTimes))),'area',1,'xaxisvalues',eeg.stimTimes);
xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from evidence onset(ms)');
box off;
legend([h{:,1}], labels.stimDur,'Location','Best');
title('evidence locked');

subplot(2,2,2);
h = errorBarPlot(groupMeans(cppRespCueFilt,3, repmat(permute(behData.stimDur,[1,3,2]),1,length(eeg.respCueTimes))),'area',1,'xaxisvalues',eeg.respCueTimes);
xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from respcue (ms)');
box off;
legend([h{:,1}], labels.stimDur,'Location','Best');
title('response cue locked');
xlim([-1000 100])

subplot(2,2,3);
h = errorBarPlot(groupMeans(cppFilt,3, repmat(permute(behData.stimDur,[1,3,2]),1,length(eeg.respTimes))),'area',1,'xaxisvalues',eeg.respTimes);
xline(0,':k');
yline(0, ':k');
xlines(amplWindows, '-k');
xlines(slopeWindows, '-r');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;
legend([h{:,1}], labels.stimDur,'Location','Best');
title('response locked');
xlim([-1500 500]);

%% plot CPP by RT bin

% 2 or 3 speed bins? uncomment one below to choose
binNames = {'fast','slow'};
% binNames = {'fast','medium','slow'};
p = linspace(0,100, length(binNames)+1);
p(1) = [];
prc = prctile(behData.RT, p,2);
binInds = sum(behData.RT > permute(prc,[1,3,2]),3);

cppByRT = groupMeans(cppFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(cpp,2),1),'dim'); %[npp t bins tr]
cppStimByRT = groupMeans(cppStimFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(cppStim,2),1),'dim'); %[npp t bins tr]
cppRespCueByRT = groupMeans(cppRespCueFilt, 3, repmat(permute(binInds,[1,3,2]),1,size(cppRespCue,2),1),'dim'); %[npp t bins tr]

figure();
subplot(2,2,1);
h = errorBarPlot(nanmean(cppStimByRT,4),'area',1,'xaxisvalues',eeg.stimTimes);
xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from stim (ms)');
box off;
title('evidence locked');

subplot(2,2,2);
h = errorBarPlot(nanmean(cppRespCueByRT,4),'area',1,'xaxisvalues',eeg.respCueTimes);
xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from resp cue (ms)');
box off;
title('response cue locked');


subplot(2,2,3)
h = errorBarPlot(nanmean(cppByRT,4),'area',1,'xaxisvalues',eeg.respTimes);
xline(0,':k');
yline(0, ':k');
ylabel('\muV'); xlabel('time from resp 1 (ms)');
box off;
title('response locked');
legend([h{:,1}], binNames,'Location','Best');

%% plot median splits per stim dur? just for resp first

figure();
durByRT = groupMeans(behData.stimDur,2,binInds,'dim');
cppByRTDur = groupMeans(cppByRT, 4, repmat(permute(durByRT,[1,4,2,3]),1,size(cpp,2),1,1)); % pp time rt dur

for i = 1:3
    subplot(2,2,i);
    h = errorBarPlot(cppByRTDur(:,:,:,i),'area',1,'xaxisvalues',eeg.respTimes);
    xline(0,':k');
    yline(0, ':k');
    ylabel('\muV'); xlabel('time from resp (ms)');
    box off;
    title(labels.stimDur{i});
end
legend([h{:,1}], binNames,'Location','Best');

%% plot CPP by RT bin, and acc

% 2 or 3 speed bins? uncomment one below to choose
accByRT = repmat(permute(groupMeans(behData.acc, 2, binInds,'dim'),[1 4 2 3]),1,size(cpp,2),1,1); %[npp t bins acc tr]
cppStimByRTAcc = groupMeans(cppStimByRT, 4, accByRT(:,1:size(cppStim,2),:,:)); %[npp t bins tr]
cppRespCueByRTAcc = groupMeans(cppRespCueByRT, 4, accByRT(:,1:size(cppRespCue,2),:,:)); %[npp t bins tr]
cppByRTAcc = groupMeans(cppByRT,4,accByRT); %[pp t rt acc]

figure();
for i = 1:2
    subplot(2,3,i*3-2);
    h = errorBarPlot(cppStimByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.stimTimes);
    xline(0,':k');
    yline(0, ':k');
    ylabel('\muV'); xlabel('time from stim (ms)');
    box off;
    title(['evidence locked, ' labels.acc{i}]);
    
    subplot(2,3,i*3-1);
    h = errorBarPlot(cppRespCueByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.respCueTimes);
    xline(0,':k');
    yline(0, ':k');
    ylabel('\muV'); xlabel('time from resp cue (ms)');
    box off;
    title(['resp cue locked, ' labels.acc{i}]);
        
    subplot(2,3,i*3)
    h = errorBarPlot(cppByRTAcc(:,:,:,i),'area',1,'xaxisvalues',eeg.respTimes);
    xline(0,':k');
    yline(0, ':k');
    ylabel('\muV'); xlabel('time from resp 1 (ms)');
    box off;
    title(['response locked, ' labels.acc{i}]);
    legend([h{:,1}], binNames,'Location','Best');
end
%% indiv plot RT splits

figure();
for iPP = 1:fileInfo.nPP
    subplot(6,5,iPP);
    h = errorBarPlot(permute(cppByRT(iPP,:,:,:),[4,2,3,1]),'area',1,'xaxisvalues',eeg.respTimes);
    xline(0,':k');
    yline(0, ':k');
    ylabel('\muV'); xlabel('time from resp 1 (ms)');
    box off;
    title(iPP);
end    
legend([h{:,1}], binNames,'Location','Best');

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
c = colorbar('Position',[0.3024 0.7024 0.0315 0.2179]);


%% Fig 3.4
figure();

subplot(2,2,1); % topoplot of activity at RT
topoplot(sq(nanmean(nanmean(cppRTTopo,2),1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',mapLims,...
    'emarker',{'.','k',10,1});
title('at RT')


subplot(2,2,2); % stimlocked trace by stimDur
stimlocked = groupMeans(cppStimFilt,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(cppStim,2)),'dim'); %[pp t stimDur tr]
eeg.stimTimes = eeg.epochTimes(isBetween(eeg.epochTimes, stimWin));
h = plot(eeg.stimTimes, sq(nanmean(nanmean(stimlocked(:,isBetween(eeg.stimTimes,stimWin),:,:),4),1)),'LineWidth',2);
xline(0);
legend(h, labels.stimDur,'Location','Best');
title('evidence onset');
xlim(stimWin);

subplot(2,2,3); % respcue
respCueLocked = groupMeans(cppRespCueFilt,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(cppRespCue,2)),'dim'); %[pp t stimDur tr]
win = [-1000 200]; % for plotting
h = plot(eeg.respCueTimes(isBetween(eeg.respCueTimes,win)), sq(nanmean(nanmean(respCueLocked(:,isBetween(eeg.respCueTimes,win),:,:),4),1)),'LineWidth',2);
xline(0);
legend(h, labels.stimDur,'Location','Best');
title('resp cue onset');
xlim(win);


subplot(2,2,4);
data = groupMeans(cppFilt,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(cpp,2)),'dim'); %[pp t stimDur tr]
win = [-500 1000];
h = plot(eeg.respTimes(isBetween(eeg.respTimes,win)), sq(nanmean(nanmean(data(:,isBetween(eeg.respTimes,win),:,:),4),1)),'LineWidth',2);
xline(0);
legend(h, labels.stimDur,'Location','Best');
title('first-order response');



%% make fig3.5

% split [stimlocked cpp] by [conf; acc;], do a topoplot
% post decision

cppMeanTopo_100 = nanmean(cppMeanTopo(:,:,:,10),4); % average over -100:0ms

nF=2; % just these two

stimWin = [-200 1500];
respWin = [-1500 1500];
ylims = [-1.5 4; -1.4 5];
    
figure();
for i = 1:nF
    subplot(nF,4,i*4-3);
    
    % split stimlocked by stimDur
    stimlocked = groupMeans(cppStimFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(cppStim,2))); %[pp t stimDur]
    stimlocked(:,:,end+1) = diff(stimlocked(:,:,[1 end]),[],3) * diffFlips(i); % get diff
    
    % plot
    h = plot(eeg.stimTimes(isBetween(eeg.stimTimes,stimWin)), sq(nanmean(stimlocked(:,isBetween(eeg.stimTimes,stimWin),:,:),1)), 'LineWidth',2);
    for j = 1:length(h); h(j).Color = cols{i}(j,:);end
    
    xline(0);yline(0);
    legend(h, [labels.(factors{i}) labels.diffNames{i}],'Location','Best');
    if i==1; title('evidence onset');end
    xlim(stimWin); %ylim(ylims(1,:));
    ylabel(factors{i});
    
    
    subplot(nF,4,i*4-2);
    % split respCueLocked by stimDur
    respCueLocked = groupMeans(cppRespCueFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(cppRespCue,2))); %[pp t stimDur]
    respCueLocked(:,:,end+1) = diff(respCueLocked(:,:,[1 end]),[],3) * diffFlips(i); % get diff
    % plot
    win = [-1000 200];
    h = plot(eeg.respCueTimes(isBetween(eeg.respCueTimes,win)), sq(nanmean(respCueLocked(:,isBetween(eeg.respCueTimes,win),:,:),1)), 'LineWidth',2);
    for j = 1:length(h); h(j).Color = cols{i}(j,:);end
    
    xline(0);yline(0);
    if i==1; title('resp cue onset');end
    xlim(win); %ylim(ylims(1,:));
    
    
    subplot(nF,4,i*4-1);
        % split cpp by stimDur
    resplocked = groupMeans(cppFilt,3,repmat(permute(behData.(factors{i}),[1,3,2]),1,size(cpp,2))); %[pp t stimDur]
    resplocked(:,:,end+1) = diff(resplocked(:,:,[1 end]),[],3) * diffFlips(i); % get diff
    h = plot(eeg.respTimes(isBetween(eeg.respTimes,respWin)), sq(nanmean(resplocked(:,isBetween(eeg.respTimes,respWin),:,:),1)), 'LineWidth',2);
    for j = 1:length(h); h(j).Color = cols{i}(j,:);end
    xline(0);yline(0);
%     xlines(amplWindows,'-k');xlines(slopeWindows, '-r');
    if i==1; title('first-order response'); end
    xlim(respWin); %ylim(ylims(2,:));
    
    
    % add a topoplot
    subplot(nF,4,i*4);

    % split by stimDur
    respTopo = groupMeans(cppMeanTopo_100,2,repmat(behData.(factors{i}),1,1,128));
    respTopo = diff(respTopo(:,[1 end],:),[],2) * diffFlips(i);
    topoplot(sq(nanmean(respTopo,1)),...
        eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
        'emarker',{'.','k',10,1});
    if i==1; title('-100:0ms to response'); end
    
end
%% fig 2.6
% columns: respcuelocked + resplocked + topo
% rows: durations
% lines = certainty

i = 2; % certainty

behDataByDur = structfun(@(x) groupMeans(x, 2, behData.stimDur,'dim'),behData,'UniformOutput',0); %[pp dur tr]


% split by stimDur, then by factor
cppRespCueFiltDur = groupMeans(cppRespCueFilt,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(cppRespCue,2)),'dim'); %[pp t dur tr]
cppFiltDur = groupMeans(cppFilt,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(cpp,2)),'dim'); %[pp t dur tr]

 % split by dur, fac
respCueLocked = groupMeans(cppRespCueFiltDur,4,repmat(permute(behDataByDur.(factors{i}),[1,4,2,3]),1, size(cppRespCue,2),1,1)); %[pp t dur factor]
respCueLocked(:,:,:,end+1) = diff(respCueLocked(:,:,:,[1 end]),[],4) * diffFlips(i); % get diff
    
% split by dur, fac
respLocked = groupMeans(cppFiltDur,4,repmat(permute(behDataByDur.(factors{i}),[1,4,2,3]),1, size(cpp,2),1,1)); %[pp t dur factor]
respLocked(:,:,:,end+1) = diff(respLocked(:,:,:,[1 end]),[],4) * diffFlips(i); % get diff
    
respTopo = groupMeans(cppMeanTopo_100,2,repmat(behData.stimDur,1,1,128));
respCueTopo = groupMeans(respCueMeanTopo(:,:,:,10),2,repmat(behData.stimDur,1,1,128));

rcWin = [-800 0];
rWin = [-1200 0];
rcInds = isBetween(eeg.respCueTimes, rcWin);
rInds = isBetween(eeg.respTimes, rWin);
%

doErrBars = 1;
figure();
for iD = 1:3
    % respcue
    subplot(3,4,iD*4-3);
    set(gca,'ColorOrder',cols{i},'nextplot','replacechildren');
    h = errorBarPlot(sq(respCueLocked(:,rcInds,iD,1:end-1)),'area',1,'xaxisvalues',eeg.respCueTimes(rcInds));
    if ~doErrBars
        for k = 1:size(h,1); h{k,2}.Visible='off';end
    end
    xline(0);
    if iD==3; xlabel('time to response cue (ms)'); end
    ylabel(labels.stimDur{iD}); 
    if iD==1
        legend([h{:,1}], [labels.(factors{i}) labels.diffNames{i}],'Location','Best'); 
        title('resp cue locked');
    end
    
    
    % topoplot - no conditions
    subplot(3,4,iD*4-2);
    topoplot(sq(nanmean(respCueTopo(:,iD,:),1)),...
        eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
        'emarker',{'.','k',10,1});
    if iD==1; title('-100:0ms to resp cue'); end
    
    % resp
    subplot(3,4,iD*4-1);
    set(gca,'ColorOrder',cols{i},'nextplot','replacechildren');
    h = errorBarPlot(sq(respLocked(:,rInds,iD,1:end-1)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
    if ~doErrBars
        for k = 1:size(h,1); h{k,2}.Visible='off';end
    end
    xline(0);
    if iD==3; xlabel('time to response (ms)'); end
    ylabel(labels.stimDur{iD}); 
    if iD==1
%         legend([h{:,1}], [labels.(factors{i}) labels.diffNames{i}],'Location','Best'); 
        title('resp cue');
    end
    
    
    % topoplot - no conditions
    subplot(3,4,iD*4);
    topoplot(sq(nanmean(respTopo(:,iD,:),1)),...
        eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
        'emarker',{'.','k',10,1});
    if iD==1; title('-100:0ms to response'); end
end


%% split by duration + certainty (accuracy seems to not do much)
cppFiltDurCert = groupMeans(cppFiltDur,4,repmat(permute(behDataByDur.certainty,[1,4,2,3]),1,size(cpp,2),1,1)); % [pp t dur cert]

figure();
for iD = 1:3
    subplot(2,3,iD);
    set(gca,'ColorOrder',cols{i},'nextplot','replacechildren');
    h = errorBarPlot(sq(cppFiltDurCert(:,rInds,iD,:)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
    if ~doErrBars
        for k = 1:size(h,1); h{k,2}.Visible='off';end
    end
    xline(0);
    xlabel('time to response (ms)');
    ylabel(labels.stimDur{iD});
    if iD==1
        legend([h{:,1}], [labels.(factors{i})],'Location','Best');
    end
end
makeSubplotScalesEqual(2,3,1:3)

% now plot with dur as legend
% figure();
for iC = 1:3
    subplot(2,3,iC+3);
%     set(gca,'ColorOrder',cols{i},'nextplot','replacechildren');
    h = errorBarPlot(sq(cppFiltDurCert(:,rInds,:,iC)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
    if ~doErrBars
        for k = 1:size(h,1); h{k,2}.Visible='off';end
    end
    xline(0);
    xlabel('time to response (ms)');
    ylabel(labels.certainty{iC});
    if iC==1
        legend([h{:,1}], [labels.stimDur],'Location','Best');
    end
end
makeSubplotScalesEqual(2,3)

%% plot traces for dur*acc*cert

% just resplocked for now

% split resplocked by accuracy too
behDataByDurAcc = structfun(@(x) groupMeans(x, 3, behDataByDur.acc,'dim'),behDataByDur,'UniformOutput',0); %[pp dur tr]

cppFiltDurAcc = groupMeans(cppFiltDur,4,repmat(permute(behDataByDur.acc,[1,4,2,3]),1,size(cpp,2),1,1),'dim'); % [pp t dur acc]
cppFiltDurAccCert = groupMeans(cppFiltDurAcc,5,repmat(permute(behDataByDurAcc.certainty,[1,5,2,3,4]),1,size(cpp,2))); % [pp t dur acc cert]

figure();
for iD = 1:3
    for iAcc = 1:2
        subplot(3,2,(iD-1)*2+iAcc);
        set(gca,'ColorOrder',cols{i},'nextplot','replacechildren');
        h = errorBarPlot(sq(cppFiltDurAccCert(:,rInds,iD,iAcc,:)),'area',1,'xaxisvalues',eeg.respTimes(rInds));
        if ~doErrBars
            for k = 1:size(h,1); h{k,2}.Visible='off';end
        end
        xline(0);
        if iD==3; xlabel('time to response (ms)'); end
        if iAcc==1; ylabel(labels.stimDur{iD}); end
        if iD==1; title(labels.acc(iAcc)); end
        if iD==1 && iAcc==1
            legend([h{:,1}], [labels.(factors{i})],'Location','Best'); 
        end
    end
end
    

%% plot means by factors (don't split by duration

figure();
for iF = 1:2

    for iV = 1:nCppVars
        subplot(2,nF,(iV-1)*(nF)+iF)
        
        % split by factor
        data = groupMeans(cppVars.(cppVarNames{iV}), 2, behData.(factors{iF})); %[pp factor]
        h = errorBarPlot(data, 'type','bar');
        
        h.FaceColor = 'flat';
        h.CData = cols{iF}(1:end-1,:);
        
        set(gca,'XTickLabel',labels.(factors{iF}));
        xlabel(factors{iF});
        
        if iF==1; ylabel(cppVarNames{iV}); end
        
        
    end
    
end


%% split by dur and then factors

varsByDur = structfun(@(x) groupMeans(x, 2, behData.stimDur,'dim'),cppVars,'UniformOutput',0);

figure();
for iF = 1:nF
    % split by stimDur %[pp stimDur 3 tr]
    factor = groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim');

    for iV = 1:nCppVars

        % split by factor
        data = groupMeans(varsByDur.(cppVarNames{iV}), 3, factor); %[pp dur 3 factor]
        
        subplot(2,nF,(iV-1)*(nF)+iF)
        h = errorBarPlot(data, 'type','bar');
        set(gca,'XTickLabel',labels.stimDur);
        xlabel('duration');
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

%% split acc/err


behDataByDur = structfun(@(x) groupMeans(x, 2, behData.stimDur,'dim'),behData,'UniformOutput',0); %[pp dur tr]
varsByDurAcc = structfun(@(x) groupMeans(x, 3, behDataByDur.acc,'dim'),varsByDur,'UniformOutput',0); %[pp dur acc tr]

% change order on x-axis
% figure();
for iF = 2%1:nF
    figure();
    factor = groupMeans(behDataByDur.(factors{iF}),3,behDataByDur.acc, 'dim');

    for iV = 1:nCppVars

        % split by factor
        data = groupMeans(varsByDurAcc.(cppVarNames{iV}), 4, factor); %[pp stimDur acc factor]
%         data = permute(data,[1 2 3 4]); % need to permute as errbpl reshapes it, which messes order; [pp stimDur factor acc
%         data = reshape(data,fileInfo.nPP,6,[]);
        
%         subplot(2,1,iV);
        for iA=1:2
            subplot(2,2,(iV-1)*2+iA);
            set(gca,'ColorOrder',cols{iF}(1:end-1,:),'nextplot','replacechildren');
            h = errorBarPlot(sq(data(:,:,iA,:)), 'type','bar');
            set(gca,'XTickLabels', labels.stimDur, 'XTickLabelRotation',0);
            xlabel('duration');
            ylabel(cppVarNames{iV}); 
            
            title(labels.acc(iA));
            if iV==1 && iA==1
                legend(h, [labels.(factors{iF})],'Location','Best'); 
            end
            
        end
        
        makeSubplotScalesEqual(2,2,iV*2-1:iV*2);
    end
    
end

%% justby acc

figure();

varsByAcc = structfun(@(x) groupMeans(x, 2, behData.acc,'dim'),cppVars,'UniformOutput',0);

figure();
for iF = 2
    factor = groupMeans(behData.(factors{iF}),2,behData.acc, 'dim');

    for iV = 1:nCppVars

        % split by factor
        data = groupMeans(varsByAcc.(cppVarNames{iV}), 3, factor); %[pp dur 3 factor]
        
        subplot(2,1,iV);
%         subplot(2,nF,(iV-1)*(nF)+iF)
        h = errorBarPlot(data, 'type','bar');
        set(gca,'XTickLabel',labels.acc);
        xlabel('duration');
        for i = 1:length(h)
            h(i).FaceColor = cols{iF}(i,:);
        end
        if iV==1
            legend(h, labels.(factors{iF}),'Location','Best'); 
            title(factors{iF});
        end
        ylabel(cppVarNames{iV}); 
        
        
    end
    
end
%% 2-way dur*factor

[anovas2, anovas2Tabs] = deal(struct());
for iF = 1:nF-1
    factor = groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim');

    for iV = 1:nCppVars
        
        data = groupMeans(varsByDur.(cppVarNames{iV}), 3, factor); %[pp stimDur 3 factor]
                
        anovas2.(factors{iF}).(cppVarNames{iV}) = rmanova(data, {'pp','stimDur',factors{iF}},'categorical',2:3, 'DummyVarCoding','effects');
        
    end
    anovas2Tabs.(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovas2.(factors{iF}),'UniformOutput',0),'rowNames',anovas2.(factors{iF}).ampl.Term(2:end));
end

%% durations separately

[anovasSep, anovasSepTabs] = deal(struct());
for iF = 1:nF-1
    factor = groupMeans(behData.(factors{iF}),2,behData.stimDur, 'dim');

    for iV = 1:nCppVars
        
        data = groupMeans(varsByDur.(cppVarNames{iV}), 3, factor); %[pp stimDur 3 factor]

        for j = 1:3   
            anovasSep.(labels.durNames{j}).(factors{iF}).(cppVarNames{iV}) = rmanova(sq(data(:,j,:)), {'pp',factors{iF}},'categorical',[], 'DummyVarCoding','effects');
        end
        
    end
    for j = 1:3
        anovasSepTabs.(labels.durNames{j}).(factors{iF}) = struct2table(structfun(@(x) x.pValue(2:end), anovasSep.(labels.durNames{j}).(factors{iF}),'UniformOutput',0),'rowNames',anovasSep.(labels.durNames{j}).(factors{iF}).ampl.Term(2:end));
    end
end



%% look at init-dec things split by factors
% 
% behVarNames2 = {'acc','RT','confResp','confInCorr','certainty','corrLR','respLR','stimDur'};
% nBehVars2 = length(behVarNames2);
% 
% alpha = .01;
% 
% form = '%s ~ 1 + %s + (1 | pp)';
% fit = {};
% 
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = 1:2
%     figure();
% 
%     for iB = 1:length(behVarNames2)
%         subplot(r,c,iB);
%         if any(strcmp(factors, behVarNames2{iB}))
%             set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%         end
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             % split by factor
%             data = groupMeans(cppVars.(cppVarNames{iV}), 2, behData.(behVarNames2{iB})); %[pp stimDur 3 factor]
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
% behDataByAcc = rmfield(behData, behNames(~ismember(behNames, [behVarNames2 {'stimDur'}]))); 
% behDataByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), behDataByAcc, 'UniformOutput',0);
% 
% cppVarsByAcc = structfun(@(x) groupMeans(x,2,behData.acc,'dim'), cppVars,'UniformOutput',0);
% 
% form = '%s ~ 1 + acc*%s + (1 | pp)';
% fit = {};
% 
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = 1:2
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
% %% look at things split by factors + stimDur
% 
% behDataByDur = rmfield(behData, behNames(~ismember(behNames, [behVarNames2 {'stimDur'}]))); 
% behDataByDur = structfun(@(x) groupMeans(x,2,behData.stimDur,'dim'), behDataByDur, 'UniformOutput',0);
% 
% cppVarsByDur = structfun(@(x) groupMeans(x,2,behData.stimDur,'dim'), cppVars,'UniformOutput',0);
% 
% form = '%s ~ 1 + stimDur*%s + (1 | pp)';
% fit = {};
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = 1:2
%     figure();
% 
%     for iB = 1:nBehVars2-1
%         subplot(r,c,iB);
%         if any(strcmp(factors, behVarNames2{iB}))
%             set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%         end
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             % split by AccFac
% 
%             data = groupMeans(cppVarsByDur.(cppVarNames{iV}), 3, behDataByDur.(behVarNames2{iB})); %[pp acc fac tr]
%         
%             
%             h = errorBarPlot(data, 'type','bar');
%             xticklabels(labels.stimDur1);
%             ylabel(cppVarNames{iV});
%             if iB==1 || iB==nBehVars2; legend(h, labels.(behVarNames2{iB}),'Location','Best'); end
% 
%         
%         else % cond plot
%             
%             [~,~,~,~,h] = conditionalPlot(permute(behDataByDur.(behVarNames2{iB}),[3,1,2]), permute(cppVarsByDur.(cppVarNames{iV}),[3 1 2]));
%             xlabel(behVarNames2{iB});
%             ylabel(cppVarNames{iV});
%             legend(h(2:2:end), labels.stimDur1,'Location','Best');
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
%     cppBehRegs.behDur.(cppVarNames{iV}) = StatsTableFromLMECells(fit,behVarNames2(1:end-1));
% 
% end
% 
% %% look at post-dec things split by factors + stimDur + acc
% 
% warning('OFF','NANCAT:emptyCells')
% behDataByDurAcc = structfun(@(x) groupMeans(x,3,behDataByDur.acc,'dim'), behDataByDur, 'UniformOutput',0);
% 
% cppVarsByDurAcc = structfun(@(x) groupMeans(x,3,behDataByDur.acc,'dim'), cppVarsByDur,'UniformOutput',0);
% warning('ON', 'NANCAT:emptyCells');
% 
% textCols = {'r','b','g','m','c','y','k'};
% 
% form = '%s ~ 1 + acc*stimDur*%s + (1 | pp)';
% fit = {};
%         
% [r,c] = GetSubPlotShape(nBehVars2);
% for iV = 1:2
%     figure();
%     colOrder = get(gca,'ColorOrder');
%     colOrder = colOrder(1:6,:);
%     
%     for iB = 1:nBehVars2-1
%         subplot(r,c,iB);
%         if length(CountUnique(behData.(behVarNames2{iB}))) <= 6
%             if any(strcmp(factors, behVarNames2{iB}))
%                 set(gca,'ColorOrder',cols{strcmp(factors, behVarNames2{iB})}(1:end-1,:),'nextplot','replacechildren');
%             else
%                 set(gca,'ColorOrder',colOrder(1:length(labels.(behVarNames2{iB})),:),'nextplot','replacechildren');
%             end
%             % split by AccFac
% 
%             data = groupMeans(cppVarsByDurAcc.(cppVarNames{iV}), 4, behDataByDurAcc.(behVarNames2{iB})); %[pp stimDur acc fac]
%             data = cat(2, sq(data(:,1,:,:)), sq(data(:,2,:,:)),sq(data(:,3,:,:)));
% 
%             h = errorBarPlot(data, 'type','bar');
% %             [h(1:length(h)/2).FaceAlpha] = deal(.5);
% 
% %             xticklabels(labels.stimDur);
%             xticklabels(col(strcat(repmat(labels.stimDur1',1,2), repmat(labels.acc,3,1))'));
%             ylabel(cppVarNames{iV});
%             legend(h, labels.(behVarNames2{iB}),'Location','Best');
% %             if iB==1 || iB==nBehVars2; legend([h(length(h)/2+1:end) h(1)], [labels.(behVarNames2{iB}) 'error...'],'Location','Best'); end
% 
%         
%         else % cond plot
%             
%             [~,~,~,~,h] = conditionalPlot(reshape(permute(behDataByDurAcc.(behVarNames2{iB}),[4,1,2,3]),[],fileInfo.nPP,6), reshape(permute(cppVarsByDurAcc.(cppVarNames{iV}),[4,1,2,3]),[],fileInfo.nPP,6));
%             xlabel(behVarNames2{iB});
%             ylabel(cppVarNames{iV});
%             if iB==2; legend(h(2:2:end), col(strcat(repmat(labels.acc,3,1), repmat(labels.stimDur1',1,2))),'Location','Best'); end
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
%     cppBehRegs.behDurAcc.(cppVarNames{iV}) = StatsTableFromLMECells(fit,behVarNames2(1:end-1));
% 
% end


%% do regressions across each time step - all factors?

% do some times
times = -1000:100:1000;
% times = eeg.respTimes(isBetween(eeg.respTimes, [-500 1000]));
stats = NaN(7,length(times),2);
for i = 1:length(times)-1
    % get resplockcumul at t
    % add to regTab
    % regress

%     regTab.amplWin = nanzscore(col(cpp(:,find(eeg.respTimes<=times(i),1,'last'),:))); % use single point

    regTab.amplWin = nanzscore(col(nanmean(cpp(:,isBetween(eeg.respTimes, times([i i+1])),:),2))); % use mean
%     regTab.slopeWin = nanzscore(col(FitCPPSlope(cpp, times([i i+1]), eeg.respTimes)));


    if sum(~isnan(regTab.amplWin)) > 100
        fit = fitglme(regTab, 'amplWin ~ 1 + stimDur*acc*certainty + (1 | pp)');
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

%% how does that compare to perm testing each one?
% 
% clusterArgs = {'cluster',1,'clustermethod','mass','nperms',5000};
% behData.hiCert = behData.certainty;
% behData.hiCert(behData.hiCert==2) = NaN;
% behData.isLong = behData.stimDur;
% behData.isLong(behData.stimDur==500) = NaN;
% 
% factor = 'hiCert';%'CoM'; % must be 2 levels only
% nRespTimes = length(eeg.respTimes);
% 
% p = [];
% [~,p] = permutationOLS(diff(groupMeans(cpp,3,repmat(permute(behData.isLong,[1,3,2]),1,nRespTimes,1)),[],3),[],[],[],clusterArgs{:});
% [~,p(2,:)] = permutationOLS(diff(groupMeans(cpp,3,repmat(permute(behData.isLong,[1,3,2]),1,nRespTimes,1)),[],3),[],[],[],clusterArgs{:});
% [~,p(3,:)] = permutationOLS(diff(groupMeans(cpp,3,repmat(permute(behData.(factor),[1,3,2]),1,nRespTimes,1)),[],3),[],[],[],clusterArgs{:});
% 
% % inter
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.isLong + 2*behData.acc,[1,3,2]),1,nRespTimes,1)),fileInfo.nPP,[],2,2);
% y = diff(diff(y,[],4),[],3);
% [~,p(4,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% 
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.isLong + 2*behData.(factor),[1,3,2]),1,nRespTimes,1)),fileInfo.nPP,[],2,2);
% y = diff(diff(y,[],4),[],3);
% [~,p(5,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.acc + 2*behData.(factor),[1,3,2]),1,nRespTimes,1)),fileInfo.nPP,[],2,2);
% y = diff(diff(y,[],4),[],3);
% [~,p(6,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% 
% y = reshape(groupMeans(cpp,3,repmat(permute(behData.isLong + 2*behData.acc + 4*behData.(factor),[1,3,2]),1,nRespTimes,1)),fileInfo.nPP,[],2,2,2);
% y = diff(diff(diff(y,[],5),[],4),[],3);
% [~,p(7,:)] = permutationOLS(y,[],[],[],clusterArgs{:});
% 
% 
% %
% terms = {'stimDur','acc',factor,'stimDur:acc',['stimDur:' factor],['acc:' factor], ['stimDur:acc:' factor]};
% figure();
% imagep(p, terms);
% xticks(1:100:1000);
% xticklabels(round(eeg.respTimes(xticks)));



%% stimlocked plots by durs

% load('./Saves/CPPAnalysis_CSD.mat','cppStimFilt','eeg','behData','labels');

% split by duration and cert

eeg.nStims = length(eeg.stimTimes);
cppByDur = groupMeans(cppStimFilt, 3, repmat(permute(behData.stimDur,[1,3,2]),1,eeg.nStims),'dim'); %[pp t dur tr]
behDataByDur.certainty = groupMeans(behData.certainty,2,behData.stimDur,'dim'); %[pp dur tr]
cppByDurCert = groupMeans(cppByDur,4,repmat(permute(behDataByDur.certainty,[1,4,2,3]),1,eeg.nStims)); %[pp t dur cert


%% plot 3 time windows, each with approp durations

figure();
colours.certainty = [0    0.4470    0.8; .0    0.7510         0;  0.8500    0.3250    0.0980; .2 .2 .2];
lines.alphas.certainty = .2;
showErrBars = 0;

tWins = [-200 350; -200 500; -200 750];

for i = 1:3 % duration
    subplot(2,2,i);
   
    set(gca,'ColorOrder',colours.certainty,'nextplot','replacechildren');
    tInds = isBetween(eeg.stimTimes, tWins(i,:));
    hold on;
    h = errorBarPlot(sq(cppByDurCert(:,tInds,i,:)),'area',1,'xaxisvalues',eeg.stimTimes(tInds),'plotargs',{'LineWidth',2});
    if 0%showErrBars==0
        for j=1:size(h,1) % still show diff wave
            h{j,2}.Visible='off'; 
%             h{j,1}.LineWidth = 2;
        end % remove error bars
    else % make more transparent
        for j=1:size(h,1)
            h{j,2}.FaceAlpha = lines.alphas.certainty;
%             h{j,1}.LineWidth = 2;
        end
    end

    % regress?
    if 1%doStats % regress factor within this cond, at every 100ms bin
        u=[350 500 750];% stimdur values

        if ~exist('regTab','var')
            regTab = table;
            regTab.pp = nanzscore(col(behData.pp)); % participant for RE
            regTab.stimDur = col(behData.stimDur);
        end
        regTab.certainty = nanzscore(col(behData.certainty)); % add this factor
        times = [0:100:u(i) u(i)];
        formula = sprintf('amplWin ~ 1 + %s + (1 | pp)', 'certainty'); % formula to run, in this cond
        % run one reg to get number of terms
        regTab.amplWin = nanzscore(col(cppStimFilt(:,1,:))); % use mean
        fit = fitglme(regTab, formula);
        stats = NaN(length(fit.CoefficientNames)-1,length(times)-1,2);
        ind = strcmp(fit.CoefficientNames,'certainty');
        for iT = 1:length(times)-1
            % take mean within window
            regTab.amplWin = nanzscore(col(nanmean(cppStimFilt(:,isBetween(eeg.stimTimes, times([iT iT+1])),:),2))); % use mean
            if sum(~isnan(regTab.amplWin)) > 100 % if at least 100 data
                fit = fitglme(regTab(regTab.stimDur==u(i),:), formula);
                stats(:,iT,2) = fit.Coefficients.pValue(2:end); % p-value
            end
        end

        hold on;
        yVal = min(ylim); 

        pbar(col([stats(find(ind)-1,:,2);stats(find(ind)-1,:,2)])', 'xVals',col([times(1:end-1);times(2:end)]), 'yVal',yVal,'plotargs',{'Color','k','LineWidth',3});

    end
    
    xline(0,'k');
    yline(0,'k');
    xlabel('time from stimulus onset (ms)');
    ylabel('CPP \muV');
    xlim(max(tWins));
    title(labels.stimDur{i});
    if i==1
        legend(flip([h{:,1}]),flip(labels.certainty),'Location','Best');

%         hold on; emptyLegend(6, { {'Color', colours.certainty(1,:),'LineWidth',2},{'Color', colours.certainty(2,:),'LineWidth',2},{'Color', colours.certainty(3,:),'LineWidth',2},...
%             {'Color','k','LineWidth',3}, {'Color','k','LineWidth',2}, {'Color','k','LineWidth',1}}, {}, [labels.certainty';labels.stimDur'], {'Location','Best'});
    end

end