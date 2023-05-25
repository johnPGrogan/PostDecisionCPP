function GetCPPStimTopoRDM
% get the mean CPP across topography - within each 100ms period
% relative to stimulus onset

load('./Saves/ExtractEpochs.mat');
eeg.stimWindowS = (-512:256) + 540; % to add on to evOnset time
eeg.stimTimes = ((-512:256) ./ 512) .* 1000; % times
eeg.nStimSamples = length(eeg.stimTimes);

useCSD = 1;
if useCSD
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\CSD\';
    intpSuffix = '_intp_csd.mat';
    saveName = 'GetCPPStimTopoCSD.mat';
else
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\Interp';
    intpSuffix = '_intp.mat';
    saveName = 'GetCPPStimTopo.mat';
end

load('./Saves/FlagArtefacts.mat','isFlagged');
load('./Saves/BehDataLoad.mat','behData');

%% 


wins = -500:100:500;
winInds = sum(eeg.stimTimes > wins'); % 1:21
% winInds(1) = 1; % include -1000 in this

tInds = isBetween(eeg.stimTimes, wins([1 end]));
winInds(find(tInds,1)) = 1; % -500 in first bin

% stimSlopeTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans);
stimMeanTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans, length(wins)-1);
grandMeanTopo = NaN(fileInfo.nPP, length(eeg.stimTimes), eeg.nChansTot);
stimEndTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChansTot, 1);
stimEndTimes = [250 350; 400 500; 650 750];

for iPP = 1:fileInfo.nPP
    
    disp(iPP);
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpSuffix]),'erp');
    
    % discard times outside wins

    stimMeanTopo(iPP,:,1:128,:) = permute(nanmean(groupMeans(data.erp(1:128,tInds,:),2,winInds(tInds),'dim'),1),[1,2,4,3]);
    
%     % get -50:50ms data too - will fit slopes to it
%     stimSlopeTopo(iPP,:,1:128) = FitCPPSlope(data.erp(1:128,:,:),[-50 50], eeg.stimTimes)';

    
    % get mean topo across trials [pp times chans]
    % remove flagged
    data.erp(:,:,isFlagged(:,iPP)) = NaN;
    grandMeanTopo(iPP,:,:) = nanmean(data.erp(:,isBetween(eeg.epochTimes, minMax(eeg.stimTimes)),:),3)'; % get eye chans too

    % also get last 100ms of evidence, differs by stimDur
    for iS = 1:3
        thisDur = behData.stimDur(iPP,:)==stimEndTimes(iS,2);
        stimEndTopo(iPP,thisDur,:) = sq(nanmean(data.erp(:,isBetween(eeg.epochTimes,stimEndTimes(iS,:)),thisDur),2))';
    end

    
end


save(fullfile('./Saves/', saveName), 'grandMeanTopo','stimMeanTopo','wins','stimEndTopo');


end