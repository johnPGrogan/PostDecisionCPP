% GetSSVEPChans
% load CPP topo, pick best chan per pp
% load that channel and do FFT for 30Hz (normalised to 28+32)
% 

clc; clear; close all;

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat');
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');

useCSD = 1;
doBaseline = 1;

if useCSD
    loadName = 'GetSSVEPTopoCSD.mat';
    saveName1 = 'GetSSVEPChans_CSD.mat';
else
    loadName = 'GetSSVEPTopo.mat';
    saveName1 = 'GetSSVEPChans.mat';
end

load(fullfile(outFolder, loadName), 'ssvepTopo', 'stimWindows', 'ssvepTopoResp','respWindows')

window = [300 500]; % for mean
winInds = isBetween(stimWindows, window);

cmap = crameri('vik');
mapLims = [-.05 .05];


%% topoplot

chans= {'A23','A24','A22','A15','A14','A28','A27','A26','A25','A13','A16','A29',...
    'A29','A30','B6','B7','A21'};
[~, chansInds] = ismember(chans, eeg.chanNames);

figure();
topoplot(sq(nanmean(nanmean(nanmean(ssvepTopo(~ppsToExclude,:,winInds,:),3),4),1)),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'mapLimits',mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2',  {chansInds, '.','k',10,1});
colorbar

chans5 = {'A29','A30','B6','B7','A21','A22'};
[~, chans5Inds] = ismember(chans5, eeg.chanNames);

%% each pp

figure();
for iPP = 1:fileInfo.nPP
    if 1%~ppsToExclude(iPP)
        subplot(5,6,iPP);
        topoplot(sq(nanmean(nanmean(ssvepTopo(iPP,:,winInds,:),3),4)),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,...
            'emarker',{'.','k',10,1}, 'emarker2', {chans5Inds, '.','k',10,1});
        title(iPP);
    end
end


%% pick best per pp 

vals = nanmean(nanmean(ssvepTopo(:,chans5Inds,winInds,:),3),4)'; % get mean voltages per channel per pp   

% pick maximal amplitude surrounding response execution    
[m,j] = max(vals);

ssvepChans = chans5(j); % get name
ssvepChanInds = chans5Inds(j); % index
    
%% load each chan

ssvep.stim = NaN(fileInfo.nPP, length(stimWindows), fileInfo.maxTr);
ssvep.resp = NaN(fileInfo.nPP, length(respWindows), fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    
    ssvep.stim(iPP,:,:) = sq(ssvepTopo(iPP,ssvepChanInds(iPP),:,:));

    ssvep.resp(iPP,:,:) = sq(ssvepTopoResp(iPP,ssvepChanInds(iPP),:,:));

end


%% get stuff into structs to match GetSSVEP.mat


windows.stim = stimWindows;
windows.resp = respWindows;
nTimes.stim = length(stimWindows);
nTimes.resp = length(respWindows);
lockNames ={'stim','resp'};

%% remove flagged

for i = 1:2
    ssvep.(lockNames{i})(repmat(permute(isFlagged,[2,3,1]),1,nTimes.(lockNames{i}),1)) = NaN;
end

%% baseline?
baselineInds = find(windows.stim <= -250,1,'last');%isBetween(windows.stim, [-150 0]);
baseline = nanmean(ssvep.stim(:, baselineInds,:),2);
    
if doBaseline
    
    for i = 1:2
        ssvep.(lockNames{i}) = log(ssvep.(lockNames{i}) ./ baseline);
    end


end

%% save

save(fullfile('./Saves/', saveName1), 'ssvep','windows','nTimes','lockNames',...
    'baseline','doBaseline','baselineInds');
