% PhaseTopoAnalysis
% load phase topographies, plot means, over time, look at indiv
% trying to decide whether to just use Oz, or pick diff per person

if ~exist('topoplot','file'); eeglab nogui; end

clc; clear; close all;

%% plot 30Hz signal phase on trials

useCSD = 1;

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat');
load('./Saves/BehDataLoad.mat','behData');
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');

if useCSD
    loadName = fullfile(outFolder, 'GetPhaseTaggingTopo_CSD.mat');
    saveName = fullfile(outFolder, 'PhaseTopoAnalysis_CSD.mat');
else
    loadName = fullfile(outFolder, 'GetPhaseTaggingTopo.mat');
    saveName = fullfile(outFolder, 'PhaseTopoAnalysis.mat');
end


load(loadName);
projWaveRespMeanWindows = real(projWaveRespMeanWindows);

cmap = crameri('vik');
mapLims = [-400 400];

time = 400; % end of 100ms time window to use
iT = find(respWins(2:end)==time);

%% plot meanat 300ms

visChanNames = {'A23','A24','A25','A26','A27','A13','A14','A22','A15','A28'};
[~, visChanInds] = ismember(visChanNames, eeg.chanNames);

figure();
topoplot(sq(nanmean(nanmean(projWaveRespMeanWindows(~ppsToExclude,iT,:,:),4),1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {visChanInds, '.','k',10,1});
c = colorbar('Limits', [0 mapLims(2)]);

%% pick 5

chans5 = {'A22','A23','A24','A25'}; % midline only
% chans5 =  {'A23','A24','A25','A26','A27'}; % mids + right
% chans5 = {'A23','A24','A25','A26','A27','A13','A14','A22','A15','A28'}; % all

[~, chans5Inds] = ismember(chans5, eeg.chanNames);

%%
figure();
for iPP = 1:fileInfo.nPP
    if 1%~ppsToExclude(iPP)
        subplot(5,6,iPP);
        topoplot(sq(nanmean(projWaveRespMeanWindows(iPP,iT,:,:),4)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap, 'mapLimits', mapLims,...
            'emarker',{'.','k',10,1}, 'emarker2', {chans5Inds, '.','k',10,1});
        title(iPP);
    end
end
        
%% pick maximal in that per person

means = sq(nanmean(projWaveRespMeanWindows(:,iT,chans5Inds,:),4))'; % get mean voltages per channel per pp   

% pick maximal amplitude surrounding response execution    
[m,j] = max(means);

visChans = chans5(j); % get name
visChanInds = chans5Inds(j); % index

%% save

save(saveName, 'visChans','visChanInds');
    


