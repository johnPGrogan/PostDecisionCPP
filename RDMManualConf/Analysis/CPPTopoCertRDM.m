% CPPTopoCertRDM
% plot certainty effects topography i.e. certain minus maybe
% first at RT
% then maybe at different times in lead up
% and potentially during stim-locked too, or maybe at respCueOnset time?

clear; clc; close all;

%%%%%% set options
opts.useCSD = 1;
opts.excludeBadPps = 1; % remove pps with <640 good trials?
opts.excludeTooFew = 1; % remove pps with <20 per conf3
opts.excludeByRT = 1; % remove trials outside [100 1500] ms
opts.doFilt = 1; % whether to plot with loPass filtered data

saveOpts = {'Volt','CSD'};
saveName = sprintf('CPPAnalysis_%s.mat', saveOpts{1,opts.useCSD+1});

% load stuff
outFolder = './Saves';
opts1 = load(fullfile(outFolder, saveName), 'useCSD','excludeBadPps','excludeTooFew','excludeByRT','doFilt');
if ~equals(opts, opts1)
    error('loaded options do not match set options above');
end

load(fullfile(outFolder, saveName), 'eeg','cppFilt','cppStimFilt','cppRespCueFilt','chans5',...
    'behData','labels','cppMeanTopo','cppRTTopo','cppTopoAmplWindow',...
    'respCueMeanTopo','fileInfo');

if useCSD
    cppName = 'GetCPPCSD.mat';
    cppTopoName = 'GetCPPTopoCSD.mat';
    cppTopoNameRC = 'GetCPPRespCueTopoCSD.mat';
    stimTopoName = 'GetCPPStimTopoCSD.mat';
    mapLims = [-20 20];
else
    cppName = 'GetCPP.mat';
    cppTopoName = 'GetCPPTopo.mat';
    cppTopoNameRC = 'GetCPPRespCueTopo.mat';
    mapLims = [-3 3];
end

% load these ones separately
% load(fullfile(outFolder, 'ExtractEpochs.mat'));
% load(fullfile(outFolder, cppName),'eeg','cppFilt','cppStimFilt','cppRespCueFilt','chans5');
% load(fullfile(outFolder, 'FlagArtefacts.mat'),'isFlagged','ppsToExclude');
% load(fullfile(outFolder, 'BehDataLoad.mat'),'behData','labels');
% load(fullfile(outFolder, cppTopoName),'cppMeanTopo','wins','cppRTTopo','cppTopoAmplWindow');
load(fullfile(outFolder, cppTopoNameRC),'wins');
respCueWins = wins(2:end);
load(fullfile(outFolder, stimTopoName),'stimEndTopo','stimMeanTopo','wins');
stimEndTopo(:,:,129:end)=[]; % remove eye chans
stimWins = wins(2:end);

% behData = rmfield(behData, {'dataMat','dataMatNames'}); % remove these

%% remove stuff - flagged, outliers, bad pps, too few trials per bin

toRemove = sq(any(isnan(cppFilt),2)); % flagged trials already out

% remove from these that aren't in loaded file
stimEndTopo(repmat(toRemove,1,1,size(stimEndTopo,3))) = NaN;
stimMeanTopo(repmat(toRemove,1,1,128,size(stimEndTopo,4))) = NaN;



%% define windows and labels

cmap = crameri('vik');

factors = {'acc','certainty','confInCorr','confResp','stimDur'};
nF = length(factors);
labels.diffNames = {'Error - Correct' ,'Low - High', 'Low - High' ,'Low - High','Slow - Fast',};
diffFlips = [-1 -1 -1 -1 -1]; % invert some to match order

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2]
        };

nF=2;

legArgs = {'Location','Best','AutoUpdate','off','box','off'};
%% topoplot each certainty level

factors = {'certainty'};
iF = 1;

% first at RT

topoByFac = groupMeans(cppRTTopo, 2, repmat(behData.(factors{iF}),1,1,128)); %[pp cert chan
% respTopo = groupMeans(cppMeanTopo_100,2,repmat(behData.(factors{i}),1,1,128));

n = size(topoByFac,2);

[r,c] = GetSubPlotShape(n);

% chanInds = ismember(eeg.chanNames, chans5); % show CPP electrodes
[~, chanInds] = ismember(chans5, eeg.chanNames);

figure();
for i = 1:n
    subplot(r,c,i); 
    topoplot(sq(nanmean(topoByFac(:,i,:),1)),...
        eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

    title(labels.(factors{iF}){i});
end
c = colorbar('Position',[0.3024 0.7024 0.0315 0.2179]);

%% now plot difference

figure();
topoplot(sq(diff(nanmean(topoByFac(:,[1 end],:),1),[],2)),...
    eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

title(sprintf('%s - %s', labels.(factors{iF}){1}, labels.(factors{iF}){end}));

c = colorbar('Position',[0.3024 0.7024 0.0315 0.2179]);

%% plot difference across times?

nWins = length(wins)-1;
[r,c] = GetSubPlotShape(nWins);
figure();
for i = 1:nWins
    subplot(r,c,i); 
    topoDiff = groupMeans(cppMeanTopo(:,:,:,i), 2, repmat(behData.(factors{iF}),1,1,128)); %[pp factor chan]
    
    topoplot(sq(diff(nanmean(topoDiff(:,[1 end],:),1),[],2)),...
        eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

    title(wins(i+1));
end
c = colorbar('Position',[0.3024 0.7024 0.0315 0.2179]);

%% do respCue locked

% respCueWins = -700:100:1000;

nRCWins = length(respCueWins)-1;
[r,c] = GetSubPlotShape(nRCWins);
figure();
for i = 1:nRCWins
    subplot(r,c,i); 
    topoDiff = groupMeans(respCueMeanTopo(:,:,:,i), 2, repmat(behData.(factors{iF}),1,1,128)); %[pp factor chan]
    
    topoplot(sq(diff(nanmean(topoDiff(:,[1 end],:),1),[],2)),...
        eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

    title(respCueWins(i+1));
end
c = colorbar('Position',[0.3024 0.7024 0.0315 0.2179]);

%% do stimlocked

% stimWins = -1000:100:500;

nStimWins = length(stimWins)-1;
[r,c] = GetSubPlotShape(nStimWins);
figure();
for i = 1:nStimWins
    subplot(r,c,i); 
    topoDiff = groupMeans(stimMeanTopo(:,:,:,i), 2, repmat(behData.(factors{iF}),1,1,128)); %[pp factor chan]
    
    topoplot(sq(diff(nanmean(topoDiff(:,[1 end],:),1),[],2)),...
        eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),'mapLimits',mapLims,...
        'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});

    title(stimWins(i+1));
end
c = colorbar('Position',[0.3024 0.7024 0.0315 0.2179]);

%% make fig3.5

% split [stimlocked cpp] by [conf; acc;], do a topoplot
% post decision

factors = {'acc','certainty'};
labels.diffNames = {'correct - error', 'certain - maybe'};
cppMeanTopo_100 = nanmean(cppMeanTopo(:,:,:,10),4); % average over -100:0ms
nPP = size(cppStimFilt,1);
nF=1; % just certainty
iF = 2;

doPerms = 1;
clusterArgs = {'cluster',1,'clustermethod','mean','two_tailed',true};

stimWin = [-200 750];
respCueWin = [-500 1000];
respWin = [-1500 200];
% ylims = [-1.5 4; -1.4 5];
ylims = [-10 60];
    
figure();
subplot(2,3,1); % stim-locked

% split stimlocked by stimDur
stimlocked = groupMeans(cppStimFilt,3,repmat(permute(behData.(factors{iF}),[1,3,2]),1,size(cppStimFilt,2))); %[pp t stimDur]
stimlocked(:,:,end+1) = diff(stimlocked(:,:,[1 end]),[],3); % get diff

% plot
%     h = plot(eeg.stimTimes(isBetween(eeg.stimTimes,stimWin)), sq(nanmean(stimlocked(:,isBetween(eeg.stimTimes,stimWin),:,:),1)), 'LineWidth',2);
%     for j = 1:length(h); h(j).Color = cols{i}(j,:);end
set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
h = errorBarPlot(stimlocked(:,isBetween(eeg.stimTimes,stimWin),:,:), 'area',1,'xaxisvalues', eeg.stimTimes(isBetween(eeg.stimTimes, stimWin)));
xline(0);yline(0);
legend([h{:,1}], [labels.(factors{iF}) labels.diffNames{iF}],legArgs{:});
title('evidence onset');
xlim(stimWin); ylim(ylims);
ylabel('CPP \muV/cm^2');
if doPerms % permOLS 
    % difference wave? 
    d = stimlocked(:,isBetween(eeg.stimTimes, stimWin),end); %[pp t 1]
    [~,p] = permutationOLS(d, [], [], [], clusterArgs{:});
    hold on;
    pbar(p, 'xVals', eeg.stimTimes(isBetween(eeg.stimTimes, stimWin)), 'yVal', min(ylims));
end

% add a topoplot
subplot(2,3,4);
% split by stimDur
respTopo = groupMeans(stimEndTopo,2,repmat(behData.(factors{iF}),1,1,128));
respTopo = diff(respTopo(:,[1 end],:),[],2);
topoplot(sq(nanmean(respTopo,1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});
title('last 100ms of evidence (across stim durations)'); 


%%%%

subplot(2,3,2); % respCuelocked
% split respCueLocked by stimDur
respCueLocked = groupMeans(cppRespCueFilt,3,repmat(permute(behData.(factors{iF}),[1,3,2]),1,size(cppRespCueFilt,2))); %[pp t stimDur]
respCueLocked(:,:,end+1) = diff(respCueLocked(:,:,[1 end]),[],3); % get diff
% plot
% h = plot(eeg.respCueTimes(isBetween(eeg.respCueTimes,win)), sq(nanmean(respCueLocked(:,isBetween(eeg.respCueTimes,win),:,:),1)), 'LineWidth',2);
% for j = 1:length(h); h(j).Color = cols{i}(j,:);end
set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
h = errorBarPlot(respCueLocked(:,isBetween(eeg.respCueTimes,respCueWin),:,:), 'area',1,'xaxisvalues', eeg.respCueTimes(isBetween(eeg.respCueTimes, respCueWin)));
xline(0);yline(0);
title('response-cue onset');
xlim(respCueWin); ylim(ylims);
if doPerms % permOLS 
    % difference wave? 
    d = respCueLocked(:,isBetween(eeg.respCueTimes, respCueWin),end); %[pp t 1]
    [~,p] = permutationOLS(d, [], [], [], clusterArgs{:});
    hold on;
    pbar(p, 'xVals', eeg.respCueTimes(isBetween(eeg.respCueTimes, respCueWin)), 'yVal', min(ylims));
end


% add a topoplot
subplot(2,3,5); % at 100ms before resp cue appears
respTopo = groupMeans(respCueMeanTopo(:,:,:,9),2,repmat(behData.(factors{iF}),1,1,128));
respTopo = diff(respTopo(:,[1 end],:),[],2);
topoplot(sq(nanmean(respTopo,1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});
title('100:200ms after response cue'); 

subplot(2,3,3); % split cpp by stimDur
resplocked = groupMeans(cppFilt,3,repmat(permute(behData.(factors{iF}),[1,3,2]),1,size(cppFilt,2))); %[pp t stimDur]
resplocked(:,:,end+1) = diff(resplocked(:,:,[1 end]),[],3); % get diff
% h = plot(eeg.respTimes(isBetween(eeg.respTimes,respWin)), sq(nanmean(resplocked(:,isBetween(eeg.respTimes,respWin),:,:),1)), 'LineWidth',2);
% for j = 1:length(h); h(j).Color = cols{i}(j,:);end
set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
h = errorBarPlot(resplocked(:,isBetween(eeg.respTimes,respWin),:,:), 'area',1,'xaxisvalues', eeg.respTimes(isBetween(eeg.respTimes, respWin)));
xline(0);yline(0);
%     xlines(amplWindows,'-k');xlines(slopeWindows, '-r');
title('response');
xlim(respWin); ylim(ylims);
if doPerms % permOLS 
    % difference wave?
    d = resplocked(:,isBetween(eeg.respTimes, respWin),end); %[pp t 1]
    [~,p] = permutationOLS(d, [], [], [], clusterArgs{:});
    hold on;
    pbar(p, 'xVals', eeg.respTimes(isBetween(eeg.respTimes, respWin)), 'yVal', min(ylims));
end

% add a topoplot
subplot(2,3,6); % at 100ms before RT
respTopo = groupMeans(cppMeanTopo(:,:,:,10),2,repmat(behData.(factors{iF}),1,1,128));
respTopo = diff(respTopo(:,[1 end],:),[],2);
topoplot(sq(nanmean(respTopo,1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits', mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {chanInds, '.','k',10,1});
title('-100:0ms to response');


%% last 100ms of stim-locked is the same as 100ms pre-respcue

% instead plot stim-locked for each stimDur


% split stimlocked by stimDur
stimByDur = groupMeans(cppStimFilt,3,repmat(permute(behData.stimDur,[1,3,2]),1,size(cppStimFilt,2)),'dim'); %[pp t stimDur tr]
certByDur = groupMeans(behData.certainty, 2, behData.stimDur,'dim'); %[pp dur tr]
stimByDurCert = groupMeans(stimByDur, 4, repmat(permute(certByDur,[1 4 2 3]),1,size(cppStimFilt,2))); %[pp t dur cert]
stimByDurCert(:,:,:,4) = diff(stimByDurCert(:,:,:,[1 3]),[],4);

figure();
for i = 1:3
    subplot(2,3,i); % stim-locked
    
    set(gca,'ColorOrder',cols{iF},'nextplot','replacechildren');
    h = errorBarPlot(sq(stimByDurCert(:,isBetween(eeg.stimTimes,stimWin),i,:)), 'area',1,'xaxisvalues', eeg.stimTimes(isBetween(eeg.stimTimes, stimWin)));
    xline(0);yline(0);
    xline(str2num(labels.stimDur1{i}),'--k');
%     xlim([0 750]);
    legend([h{:,1}], [labels.(factors{iF}) labels.diffNames{iF}],legArgs{:});
    title('time from evidence onset (ms)');
    xlim(stimWin); %ylim(ylims(1,:));
    ylabel(factors{iF});
    title(labels.stimDur{i});
end