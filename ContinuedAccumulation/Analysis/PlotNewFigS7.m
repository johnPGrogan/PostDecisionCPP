function PlotNewFigS7()
% plot stim/respCue/resp locked expt 1 traces,
% and the certain-maybe difference topoplots
gf = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',14);
%%%%%% set options
opts.useCSD = 1;
opts.excludeBadPps = 1; % remove pps with <640 good trials?
opts.excludeTooFew = 1; % remove pps with <20 per conf3
opts.excludeByRT = 1; % remove trials outside [100 1500] ms
opts.doFilt = 1; % whether to plot with loPass filtered data

saveOpts = {'Volt','CSD'};
saveName = sprintf('CPPAnalysis_%s.mat', saveOpts{1,opts.useCSD+1});

% load stuff
% outFolder = './Saves';
outFolder = '../../RDMManualConf/Analysis/Saves';
opts1 = load(fullfile(outFolder, saveName), 'useCSD','excludeBadPps','excludeTooFew','excludeByRT','doFilt');
if ~equals(opts, opts1)
    error('loaded options do not match set options above');
end

load(fullfile(outFolder, saveName), 'eeg','cppFilt','cppStimFilt','cppRespCueFilt','chans5',...
    'behData','labels','cppMeanTopo','respCueMeanTopo');

if opts.useCSD

    stimTopoName = 'GetCPPStimTopoCSD.mat';
    mapLims = [-20 20];
else

    mapLims = [-3 3];
end




load(fullfile(outFolder, stimTopoName),'stimEndTopo');
stimEndTopo(:,:,129:end)=[]; % remove eye chans


[~, chanInds] = ismember(chans5, eeg.chanNames);

%% remove stuff - flagged, outliers, bad pps, too few trials per bin

toRemove = sq(any(isnan(cppFilt),2)); % flagged trials already out

% remove from these that aren't in loaded file
stimEndTopo(repmat(toRemove,1,1,size(stimEndTopo,3))) = NaN;



%% define windows and labels

cmap = crameri('vik');

factors = {'acc','certainty'};
labels.diffNames = {'correct - error', 'certain - maybe'};
iF = 2;

cols = {[0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
        [  0    0.4470    0.8; .0    0.7510   0;  0.8500    0.3250    0.0980; .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [crameri('roma',6); .2 .2 .2];
        [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; .2 .2 .2]
        };


legArgs = {'Location','Best','AutoUpdate','off','box','off'};


data = workspace2struct();

figure();
doPlots(data);

% reset
set(0,'DefaultAxesFontSize',gf); % reset

end

function doPlots(data)

struct2workspace(data); % unload

%% make fig3.5

% split [stimlocked cpp] by [conf; acc;], do a topoplot
% post decision



doPerms = 1;
clusterArgs = {'cluster',1,'clustermethod','mean','two_tailed',true};

stimWin = [-200 750];
respCueWin = [-500 1000];
respWin = [-1500 200];
% ylims = [-1.5 4; -1.4 5];
ylims = [-10 60];
    
% figure();
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


end