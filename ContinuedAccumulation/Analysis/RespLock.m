% RespLock
% re-epoch to the initial response -1000:1000ms


clear; close all;

useCSD = 1; % CSD or voltage?
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\ResplockCSD_2\';
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD\'; % my channels
    intpName = '_intp_csd.mat';
else
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Resplock\';
    intpName = '_intp.mat';
end

fileInfo.rawFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Raw\';

load('./Saves/FlagArtefacts3.mat','isFlagged');

%%

[respInd, respSamp] = deal(NaN(fileInfo.maxTr, fileInfo.nPP));

eeg.respWindowS = 513 + [-512:512]; % 513 is time 0 (evOnset)
eeg.respTimes = ((-512:512) ./ 512) .* 1000; % times

for iPP = 1:fileInfo.nPP
% try
    disp([fileInfo.ppID{iPP} '...'])
    
    
    % get trigger times
    a = load(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']), 'trigs','sTimes');
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp','RS');
    
    respSamp(:,iPP) = data.RS; % sample of first response
    
    respErp = NaN(eeg.nChansTot, 1025, fileInfo.maxTr);
    % get from -1000:1000 around response
    for i = 1:length(a.trigs)
        if ~isnan(respSamp(i,iPP))
            windowInds = respSamp(i,iPP) + eeg.respWindowS;
            
            respErp(:,:,i) = data.erp(:, windowInds, i);
        end
    end
    
    %% get per time step - across all trials
    
    respErp(:,:,isFlagged(:,iPP)) = NaN; % remove flagged
    preRespMeans = sq(nanmean(respErp,3));
    
    %% save the resp-locked
    
    save(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'-v7.3',...
        'preRespMeans','respErp','eeg');
% catch ME
%     disp(ME);
%     keyboard;
% end
end