% RespLockRDM
% re-epoch to RT, remove flagged trials
% takes around 30 mins for whole script


clear; close all;

useCSD = 1; % CSD or voltage?
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\ResplockCSD\';
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\CSD\'; % my channels
    intpName = '_intp_csd.mat';
else
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\Resplock\';
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\Interp\'; % not projectsHDD
    intpName = '_intp.mat';
end

fileInfo.rawFolder = 'D:\TCD\Projects\RDMManualConf\Data\Raw\';

load('./Saves/FlagArtefacts.mat','isFlagged','art');

%%

[respInd, respSamp] = deal(NaN(fileInfo.maxTr, fileInfo.nPP));

% 540 is 1050ms, i.e. sample of stimOnset
eeg.respWindowS = (-768:768) + 540; % to add on to evOnset time
eeg.respTimes = ((-768:768) ./ 512) .* 1000; % times
eeg.nRespSamples = length(eeg.respTimes);

for iPP = 1:fileInfo.nPP
% try
    disp([fileInfo.ppID{iPP} '...'])
    
    
    % get trigger times
    a = load(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']), 'trigs','sTimes','stimDurSample');
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp','RS');
    
    respSamp(:,iPP) = data.RS + a.stimDurSample; % evOnset is different in each trial
    
    respErp = NaN(eeg.nChansTot, eeg.nRespSamples, fileInfo.maxTr);
    % get from -1000:1000 around response
    for i = 1:length(a.trigs)
        if ~isnan(respSamp(i,iPP)) && isBetween(data.RS(i),  art.rtLimsSamp)
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