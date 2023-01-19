% GetCPP
% Topoplot the grand mean to choose 5 channels
% Then pick the maximal of those per person
% Outlier exclusion then

if ~exist('filter_tf','file'); eeglab; close all; end
clear; close all;

loadUp = 1; % load up PreRespMeanVoltage file?
useCSD = 1;

outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\ResplockCSD_2\';
    fileInfo.confRespFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSDConfResp\';
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD\';
    loadName = 'PreRespMeanVoltageCSD.mat';
    saveName = 'GetCPPCSD.mat';
    intpSuffix = '_intp_csd.mat';
    mapLims = [-15 15];
else
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Resplock\';
    loadName = 'PreRespMeanVoltage.mat';
    saveName = 'GetCPP.mat';
    intpSuffix = '_intp.mat';
    mapLims = [-3 3];

end

load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');
% cf = load(fullfile(outFolder, 'FlagArtefactsConfResp.mat'),'isFlagged','ppsToExclude');


eeg.respWindowS = 513 + (-512:512); % 513 is time 0 (evOnset)
eeg.respTimes = ((-512:512) ./ 512) .* 1000; % times


%%

preRespWindow = [-170 -50]; % ms before initial RT
preRespInds = isBetween(eeg.respTimes, preRespWindow);

confRespInds = isBetween(eeg.confRespTimes, preRespWindow);

preRespMeans = NaN(eeg.nChans, fileInfo.nPP);

%% get grand mean to choose channels

if loadUp && exist(fullfile(outFolder, loadName), 'file')
    r = load(fullfile(outFolder, loadName),'preRespMeans','fileInfo');
    
    ppInds = 1:30;
    preRespMeans = r.preRespMeans(:,ppInds);
     
else

    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])

        % get resplocked data
        data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'preRespMeans');
        
        % get mean over time
        preRespMeans(:,iPP) = nanmean(data.preRespMeans(1:eeg.nChans, preRespInds),2);
    
    end

    save(fullfile(outFolder, loadName), 'preRespMeans', 'preRespWindow', 'fileInfo');
end

%% topoplot that


% plot average topography from -150:-50ms before response
cppChanNames = {'A1','A2','A3','A4','A5','A6','A7','A18','A19','A20','A31','A32',...
    'B1','B2','B3','B4','B20','C1','D1','D14','D15','D16'};
[~, cppChanInds] = ismember(cppChanNames, eeg.chanNames);

cmap = 'jet'; % crameri('vik');

figure();
topoplot(nanmean(preRespMeans(:,~ppsToExclude),2),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {cppChanInds, '.','k',10,1});
colorbar


% pick 5 electrode sites of the positive signal centroparietally
if useCSD
    chans5 = {'A2','A3','A4','A19','D16'}; % CSD - from max
else
    chans5 = {'A2','A3','A4','D15','D16','B2'}; % Wouter's ones
end


[~, chans5Inds] = ismember(chans5, eeg.chanNames);

keyboard; %%%%%% edit these picked channels

%% topoplot per person

figure();
for iPP = 1:fileInfo.nPP
    if 1%~ppsToExclude(iPP)
        subplot(5,6,iPP);
        topoplot(preRespMeans(:,iPP),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,...
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
stimWin = [-500 1200];
stimWinInds = isBetween(eeg.epochTimes, stimWin);
eeg.stimTimes = eeg.epochTimes(stimWinInds);

cpp = NaN(fileInfo.nPP, length(eeg.respTimes), fileInfo.maxTr);
cppStim = NaN(fileInfo.nPP, sum(stimWinInds), fileInfo.maxTr);

% eeg.confRespTimes = ((-512:52) ./ 512) .* 1000; % times
confWin = minMax(eeg.confRespTimes,2);
% cppConfResp = NaN(fileInfo.nPP, length(eeg.confRespTimes), fileInfo.maxTr);

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
    

%     % and confResp locked
%     data = load(fullfile(fileInfo.confRespFolder, [fileInfo.ppID{iPP} '_intp_csd']),'erp');
%     
%     % get CPP channel
%     cppConfResp(iPP, :, :) = permute(data.erp(cppChanInds(iPP), :, :),[4,2,3,1]);
% 
%     % remove flagged
%     cppConfResp(iPP,:,cf.isFlagged(:,iPP)) = NaN;
%     
%     % redo means, as different flaggin used
%     cppMeans= sq(nanmean(cppConfResp(iPP, confRespInds,:),2)); % mean over window
%     
    
end



%% save a filtered copy - for plotting only

loPass = 10; %Hz

cppFilt = filterCPP(cpp, loPass, eeg, fileInfo);
cppStimFilt = filterCPP(cppStim, loPass, eeg, fileInfo);
cppConfRespFilt = filterCPP(cppConfResp, loPass, eeg, fileInfo);

%% save

save(fullfile(outFolder, saveName), ...
    'cpp','cppChans','cppChanInds','chans5','preRespWindow','preRespInds','fileInfo','eeg',...
    'cppConfResp','cppConfRespFilt','confWin',...
    'cppStim','stimWin','cppFilt','cppStimFilt');
    

%%%%%%%%%%%%%%%%%%%
%%

function cppFilt = filterCPP(cpp, loPass, eeg, fileInfo)
% replace nan with zero to allow matlab2021b to use filtfilt, 
    isAllNaN = repmat(all(isnan(cpp),2),1,size(cpp,2),1);
    cpp(isAllNaN) = 0;
    
    % filter
    cppFilt = eegfilt(reshape(permute(cpp,[1,3,2]),fileInfo.nPP*fileInfo.maxTr,size(cpp,2)),eeg.fs,0, loPass)'; % low-pass
    cppFilt = permute(reshape(cppFilt,size(cpp,2),fileInfo.nPP,fileInfo.maxTr),[2,1,3]);

    % set zeros back to NAN
    cppFilt(isAllNaN) = NaN;

end