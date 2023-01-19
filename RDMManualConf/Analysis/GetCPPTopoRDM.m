function GetCPPTopoRDM
% get the mean CPP across topography - within each 100ms period

load('./Saves/ExtractEpochs.mat');
eeg.respWindowS = (-768:768) + 540; % to add on to evOnset time
eeg.respTimes = ((-768:768) ./ 512) .* 1000; % times
eeg.nRespSamples = length(eeg.respTimes);

useCSD = 1;
if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\ResplockCSD\';
    saveName = 'GetCPPTopoCSD.mat';
else
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\Resplock';
    saveName = 'GetCPPTopo.mat';
end

load('./Saves/FlagArtefacts.mat','isFlagged');

%% 


wins = -1000:100:1000;
winInds = sum(eeg.respTimes > wins'); % 1:21
winInds(1) = 1; % include -1000 in this

tInds = isBetween(eeg.respTimes, wins([1 end]));

cppRTTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans); % at RT sample
cppMeanTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans, length(wins)); % in 100ms bins
grandMeanTopo = NaN(fileInfo.nPP, length(eeg.respTimes), eeg.nChansTot); % av over all trials
cppSlopeTopo = cppRTTopo; % slope around RT
cppTopoAmplWindow = cppRTTopo;% in amplWindow

amplWindow = [-150 -50]; %ms before response
amplWindowInds = isBetween(eeg.respTimes, amplWindow);

for iPP = 1:fileInfo.nPP
    
    disp(iPP);
    data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'respErp');

    % within mean ampl window
    cppTopoAmplWindow(iPP,:,1:128) =  sq(nanmean(data.respErp(1:128,amplWindowInds,:),2))'; % RT
    
    % mean in each bin
    cppMeanTopo(iPP,:,1:128,:) = permute(nanmean(groupMeans(data.respErp(1:128,tInds,:),2,winInds(tInds),'dim'),1),[1,2,4,3]);
    
    cppRTTopo(iPP,:,1:128) =  sq(nanmean(data.respErp(1:128,find(eeg.respTimes>0,1),:),2))'; % RT

    % get -50:50ms data too - will fit slopes to it
    cppSlopeTopo(iPP,:,1:128) = FitCPPSlope(data.respErp(1:128,:,:),[-50 50], eeg.respTimes)';
    
    % get mean topo across trials [pp times chans]
    % remove flagged
    data.respErp(:,:,isFlagged(:,iPP)) = NaN;
    grandMeanTopo(iPP,:,:) = nanmean(data.respErp,3)'; % get eye chans too
%     
end

cppMeanTopo(:,:,:,1) = []; % remove <1000ms bin

save(fullfile('./Saves/', saveName), 'cppTopoAmplWindow','amplWindow', 'cppSlopeTopo', 'grandMeanTopo','cppMeanTopo','wins','cppRTTopo');


end