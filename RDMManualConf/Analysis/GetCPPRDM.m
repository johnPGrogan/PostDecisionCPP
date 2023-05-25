% GetCPPRDM
% Topoplot the grand mean to choose 5 channels
% Then pick the maximal of those per person
% Outlier exclusion then

clear; close all;

loadUp = 1; % load up PreRespMeanVoltage file?
useCSD = 1;

outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\ResplockCSD\';
    fileInfo.respCueFolder = 'D:\TCD\Projects\RDMManualConf\Data\RespCueLockCSD\';
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\CSD\';
    loadName = 'PreRespMeanVoltageCSD.mat';
    saveName = 'GetCPPCSD.mat';
    intpSuffix = '_intp_csd.mat';
    mapLims = [-20 20];
else
    fileInfo.respFolder = 'D:\TCD\Projects\RDMManualConf\Data\Resplock\';
    fileInfo.respCueFolder = 'D:\TCD\Projects\RDMManualConf\Data\RespCueLock\';
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\Interp\'; % not projectsHDD
    loadName = 'PreRespMeanVoltage.mat';
    saveName = 'GetCPP.mat';
    intpSuffix = '_intp.mat';
    mapLims = [-3 3];

end

load(fullfile(outFolder, 'FlagArtefacts.mat'),'isFlagged','ppsToExclude');


eeg.respWindowS = -768:768; % to add on to evOnset time
eeg.respTimes = ((-768:768) ./ 512) .* 1000; % times
eeg.nRespSamples = length(eeg.respTimes);

% resp cue locked times
eeg.respCueWindowS = (-512:512) + 540; % to add on to evOnset time
eeg.respCueTimes = ((-512:512) ./ 512) .* 1000; % times
eeg.nRespCueSamples = length(eeg.respCueTimes);

%%

preRespWindow = [-170 -50]; % time window to use 
preRespInds = isBetween(eeg.respTimes, preRespWindow);

preRespMeans = NaN(eeg.nChans, fileInfo.nPP);

%% get grand mean to choose channels

if loadUp && exist(fullfile(fileInfo.respFolder, loadName), 'file')
    r = load(fullfile(fileInfo.respFolder, loadName),'preRespMeans','fileInfo');
    
    ppInds = 1:fileInfo.nPP;
    preRespMeans = r.preRespMeans(:,ppInds);
     
else

    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])

        % get resplocked data
        data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'preRespMeans');

        % mean in window
        preRespMeans(:,iPP) = nanmean(data.preRespMeans(1:eeg.nChans, preRespInds),2);

    end

    save(fullfile(fileInfo.respFolder, loadName), 'preRespMeans', 'preRespWindow', 'fileInfo');
end

%% topoplot that


% plot average topography from -150:-50ms before response
cppChanNames = {'A4','A19','A20','A21','A22','A32','A31','A30','A29','A16','A17','A18','A5'};
[~, cppChanInds] = ismember(cppChanNames, eeg.chanNames);


% pick 5 electrode sites of the positive signal centroparietally
if useCSD
    chans5 = {'A4','A19','A20','A21','A22'}; % CSD - from max
else
    chans5 = {'A4','A19','A20','A5','A32'}; % voltage - from max
end

[~, chans5Inds] = ismember(chans5, eeg.chanNames);

cmap = 'jet'; % crameri('vik');

figure();
topoplot(nanmean(preRespMeans(:,~ppsToExclude),2),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,...
    'emarker',{'.','k',10,1}, 'emarker2', {chans5Inds, '.','k',10,1});
colorbar



keyboard; %%%%%% edit these picked channels

%% topoplot per person

figure();
for iPP = 1:fileInfo.nPP
    if ~ppsToExclude(iPP)
        subplot(5,6,iPP);
        topoplot(preRespMeans(:,iPP),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',[-1 1]*max(abs(preRespMeans(:,iPP))),...
            'emarker',{'.','k',10,1}, 'emarker2', {chans5Inds, '.','k',10,1});
        title(iPP);
    end
end
        
%% pick maximal in that per person

cpps = preRespMeans(chans5Inds,:); % get mean voltages per channel per pp   

% pick maximal amplitude surrounding response execution    
[m,j] = max(cpps);

cppChans = chans5(j); % get name
cppChanInds = chans5Inds(j); % index
    

%% load that data again and store cpp chans, remove outliers
stimWin = [-200 2250];
stimWinInds = isBetween(eeg.epochTimes, stimWin);
eeg.stimTimes = eeg.epochTimes(stimWinInds);


cpp = NaN(fileInfo.nPP, length(eeg.respTimes), fileInfo.maxTr);
cppStim = NaN(fileInfo.nPP, sum(stimWinInds), fileInfo.maxTr);
cppRespCue = NaN(fileInfo.nPP, length(eeg.respCueTimes), fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'respErp');
    
    % store that channel as cpp
    cpp(iPP, :, :) = permute(data.respErp(cppChanInds(iPP), :, :),[4,2,3,1]);
    
    cpp(iPP,:,isFlagged(:,iPP)) = NaN; % remove flagged trials before mean calc
    
    
    % also get stimlocked- outlier remove by resplocked
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpSuffix]),'erp');

    % get CPP channel
    cppStim(iPP, :, :) = permute(data.erp(cppChanInds(iPP), stimWinInds, :),[4,2,3,1]);

    % remove flagged
    cppStim(iPP,:,isFlagged(:,iPP)) = NaN;
    
    
    % and get respCueLocked
    % also get stimlocked- outlier remove by resplocked
    data = load(fullfile(fileInfo.respCueFolder, [fileInfo.ppID{iPP} '_respCue']),'respCueErp');

    % get CPP channel
    cppRespCue(iPP, :, :) = permute(data.respCueErp(cppChanInds(iPP), :, :),[4,2,3,1]);

    % remove flagged
    cppRespCue(iPP,:,isFlagged(:,iPP)) = NaN;
    
    
end


%% save a filtered copy - just for plotting

loPass = 10; %Hz
% filtfilt in matlab2021b doesn't accept NaN, so set to zero then replace
isAllNaN = repmat(all(isnan(cpp),2),1,length(eeg.respTimes),1); 
cpp(isAllNaN) = 0;
% 10Hz lo-pass
cppFilt = eegfilt(reshape(permute(cpp,[1,3,2]),fileInfo.nPP*fileInfo.maxTr,length(eeg.respTimes)),eeg.fs,0, loPass)'; % low-pass
cppFilt = permute(reshape(cppFilt,length(eeg.respTimes),fileInfo.nPP,fileInfo.maxTr),[2,1,3]);
% set zeros back to NAN
cpp(isAllNaN) = NaN;
cppFilt(isAllNaN) = NaN;

isAllNaNStim = repmat(all(isnan(cppStim),2),1,length(eeg.stimTimes),1);
cppStim(isAllNaNStim) = 0;
cppStimFilt = eegfilt(reshape(permute(cppStim,[1,3,2]),fileInfo.nPP*fileInfo.maxTr,length(eeg.stimTimes)), eeg.fs, 0, loPass)'; % low-pass
cppStimFilt = permute(reshape(cppStimFilt,length(eeg.stimTimes),fileInfo.nPP,fileInfo.maxTr),[2,1,3]);
cppStim(isAllNaNStim) = NaN;
cppStimFilt(isAllNaNStim) = NaN;


isAllNaNStim = repmat(all(isnan(cppRespCue),2),1,length(eeg.respCueTimes),1);
cppRespCue(isAllNaNStim) = 0;
cppRespCueFilt = eegfilt(reshape(permute(cppRespCue,[1,3,2]),fileInfo.nPP*fileInfo.maxTr,length(eeg.respCueTimes)), eeg.fs, 0, loPass)'; % low-pass
cppRespCueFilt = permute(reshape(cppRespCueFilt,length(eeg.respCueTimes),fileInfo.nPP,fileInfo.maxTr),[2,1,3]);
cppRespCue(isAllNaNStim) = NaN;
cppRespCueFilt(isAllNaNStim) = NaN;


%% save

save(fullfile(outFolder, saveName), ...
    'cpp','cppChans','cppChanInds','chans5','preRespWindow','preRespInds','fileInfo','eeg',...
    'cppStim','stimWin','cppFilt','cppStimFilt',...
    'cppRespCue','cppRespCueFilt');