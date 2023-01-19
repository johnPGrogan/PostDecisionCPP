function GetCPPTopo
% get the mean CPP across topography - within each 100ms period

load('./Saves/ExtractEpochs.mat');
eeg.respWindowS = 513 + (-512:512); % 513 is time 0 (evOnset)
eeg.respTimes = ((-512:512) ./ 512) .* 1000; % times

useCSD = 1;
if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\ResplockCSD_2\';
    saveName = 'GetCPPTopoCSD.mat';
else
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Resplock';
    saveName = 'GetCPPTopo.mat';
end

load('./Saves/FlagArtefacts3.mat','isFlagged');

%% 


wins = -1000:100:1000;
winInds = sum(eeg.respTimes > wins'); % 1:21
winInds(1) = 1; % include -1000 in this

cppRTTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans); % at time of RT
cppMeanTopo = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans, length(wins)-1); % in each 100ms window
grandMeanTopo = NaN(fileInfo.nPP, length(eeg.respTimes), eeg.nChansTot); % mean across all trials
cppSlopeTopo = cppRTTopo; % CPP slope around RT

cppTopoAmplWindow = repmat(cppRTTopo,1,1,1,2);% in both amplitude windows
amplWindows = [-140 -50; 500 800]; % ms
for i = 1:2 % get indices within these windows
    amplWindowInds(i,:) = isBetween(eeg.respTimes, amplWindows(i,:));
end


for iPP = 1:fileInfo.nPP
    
    disp(iPP);
    data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'respErp');

    % within mean ampl windows
    for i = 1:2
        cppTopoAmplWindow(iPP,:,1:128,i) =  sq(nanmean(data.respErp(1:128,amplWindowInds(i,:),:),2))'; % RT
    end
    
    % mean within each 100ms window
    cppMeanTopo(iPP,:,1:128,:) = permute(nanmean(groupMeans(data.respErp(1:128,:,:),2,winInds,'dim'),1),[1,2,4,3]);
    
    cppRTTopo(iPP,:,1:128) =  sq(nanmean(data.respErp(1:128,513,:),2))'; % at RT sample

    % get -50:50ms data too - fit slopes to it
    cppSlopeTopo(iPP,:,1:128) = FitCPPSlope(data.respErp(1:128,:,:),[-50 50], eeg.respTimes)';

    
    % get mean topo across trials [pp times chans]
    % remove flagged
    data.respErp(:,:,isFlagged(:,iPP)) = NaN;
    grandMeanTopo(iPP,:,:) = nanmean(data.respErp,3)'; % get eye chans too
    
end

save(fullfile('./Saves/', saveName), 'cppTopoAmplWindow','amplWindows','cppSlopeTopo', 'grandMeanTopo','cppMeanTopo','wins','cppRTTopo');


end