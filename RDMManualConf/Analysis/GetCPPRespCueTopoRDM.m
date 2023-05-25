function GetCPPRespCueTopoRDM
% get the mean CPP across topography - within each 100ms period
% relative to response cue onset (evidence offset)

load('./Saves/ExtractEpochs.mat');
eeg.respCueWindowS = (-512:512) + 540; % to add on to evOnset time
eeg.respCueTimes = ((-512:512) ./ 512) .* 1000; % times
eeg.nRespCueSamples = length(eeg.respCueTimes);

useCSD = 1;
if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\RespCueLockCSD\';
    saveName = 'GetCPPRespCueTopoCSD.mat';
else
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\RespLock';
    saveName = 'GetCPPRespCueTopo.mat';
end

load('./Saves/FlagArtefacts.mat','isFlagged');

%% 


wins = -700:100:1000;
winInds = sum(eeg.respCueTimes > wins'); % 1:21
% winInds(1) = 1; % include -1000 in this

tInds = isBetween(eeg.respCueTimes, wins([1 end]));
winInds(find(tInds,1)) = 1; % -500 in first bin

respCueSlopeTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans);
respCueMeanTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans, length(wins)-1);
grandMeanTopo = NaN(fileInfo.nPP, length(eeg.respCueTimes), eeg.nChansTot);

for iPP = 1:fileInfo.nPP
    
    disp(iPP);
    data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_respCue']),'respCueErp');
    
    % discard times outside wins

    respCueMeanTopo(iPP,:,1:128,:) = permute(nanmean(groupMeans(data.respCueErp(1:128,tInds,:),2,winInds(tInds),'dim'),1),[1,2,4,3]);
    
    % get -50:50ms data too - will fit slopes to it
    respCueSlopeTopo(iPP,:,1:128) = FitCPPSlope(data.respCueErp(1:128,:,:),[-50 50], eeg.respCueTimes)';

    
    % get mean topo across trials [pp times chans]
    % remove flagged
    data.respCueErp(:,:,isFlagged(:,iPP)) = NaN;
    grandMeanTopo(iPP,:,:) = nanmean(data.respCueErp,3)'; % get eye chans too
    
end


save(fullfile('./Saves/', saveName), 'respCueSlopeTopo', 'grandMeanTopo','respCueMeanTopo','wins');


end