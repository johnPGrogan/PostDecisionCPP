% RespCueLockRDM
% re-epoch to the time of response cue onset (evOffset), remove flagged
% trials
% takes around 30 mins for whole script


clear; close all;

useCSD = 1; % CSD or voltage?
outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

if useCSD
    fileInfo.respCueFolder = 'D:\TCD\Projects\RDMManualConf\Data\RespCueLockCSD\';
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\CSD\'; % my channels
    intpName = '_intp_csd.mat';
else
    fileInfo.respCueFolder = 'D:\TCD\Projects\RDMManualConf\Data\RespCueLock\';
    intpName = '_intp.mat';
end

fileInfo.rawFolder = 'D:\TCD\Projects\RDMManualConf\Data\Raw\';

load('./Saves/FlagArtefacts.mat','isFlagged','art');

%%

[respCueSamp] = deal(NaN(fileInfo.maxTr, fileInfo.nPP));

% 540 is 1050ms, i.e. sample of stimOnset
eeg.respCueWindowS = (-512:512) + 540; % to add on to evOnset time
eeg.respCueTimes = ((-512:512) ./ 512) .* 1000; % times
eeg.nRespCueSamples = length(eeg.respCueTimes);

for iPP = 1:fileInfo.nPP
% try
    disp([fileInfo.ppID{iPP} '...'])
    
    
    % get trigger times
    a = load(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']), 'trigs','sTimes','stimDurSample');
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp');
    
    respCueSamp(:,iPP) = a.stimDurSample; % evOnset is different in each trial
    
    respCueErp = NaN(eeg.nChansTot, eeg.nRespCueSamples, fileInfo.maxTr);
    % get from -1000:1000 around response cue
    for i = 1:length(a.trigs)
        if ~isnan(respCueSamp(i,iPP))
            windowInds = respCueSamp(i,iPP) + eeg.respCueWindowS;
            
            respCueErp(:,:,i) = data.erp(:, windowInds, i);
        end
    end
    
    %% get per time step - across all trials
    
    respCueErp(:,:,isFlagged(:,iPP)) = NaN; % remove flagged
    preRespCueMeans = sq(nanmean(respCueErp,3));
    
    %% save the respCue-locked
    
    save(fullfile(fileInfo.respCueFolder, [fileInfo.ppID{iPP} '_respCue']),'-v7.3',...
        'preRespCueMeans','respCueErp','eeg');
% catch ME
%     disp(ME);
%     keyboard;
% end
end