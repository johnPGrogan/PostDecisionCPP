% GetFC
% Pick one frontal electrode per person from grand topography
% Topoplot the grand mean to choose 5 channels
% Then pick the maximal of those per person
% Outlier exclusion then

if ~exist('filter_tf','file'); eeglab; close all; end

clear; close all; clc

loadUp = 1; % load up PreRespMeanVoltage file?
useCSD = 1;
cmap =  crameri('vik');

outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

if useCSD
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\ResplockCSD_2\';
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD\';
    loadName = 'PostRespMeanVoltageCSD.mat';
    saveName = 'GetFCCSD.mat';
    intpSuffix = '_intp_csd.mat';
    mapLims = [-15 15];
else
    fileInfo.respFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Resplock\';
    loadName = 'PostRespMeanVoltage.mat';
    saveName = 'GetFC.mat';
    intpSuffix = '_intp.mat';
    mapLims = [-3 3];

end

load(fullfile(outFolder, 'FlagArtefacts3.mat'),'isFlagged','ppsToExclude');


eeg.respWindowS = 513 + (-512:512); % 513 is time 0 (evOnset)
eeg.respTimes = ((-512:512) ./ 512) .* 1000; % times


%%

postRespWindow = [850 1000]; % time window
postRespInds = isBetween(eeg.respTimes, postRespWindow);

postRespMeans = NaN(eeg.nChans, fileInfo.nPP);

%% get grand mean to choose channels

if loadUp && exist(fullfile(outFolder, loadName), 'file')
    r = load(fullfile(outFolder, loadName),'postRespMeans','fileInfo');

    ppInds = 1:30;
    postRespMeans = r.postRespMeans(:,ppInds);
  
else

    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])

        % get resplocked data
        data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'preRespMeans');

        % get mean
        postRespMeans(:,iPP) = nanmean(data.preRespMeans(1:eeg.nChans, postRespInds),2);

    end

    save(fullfile(outFolder, loadName), 'postRespMeans', 'postRespWindow', 'fileInfo');
end

%% topoplot that


% plot average topography from -150:-50ms before response
fcChanNames = {'C23','C22','C21','C25','C24','C11','C12',...
    'D2','D1','C1','C2','A1','B1','D15'};
[~, fcChanInds] = ismember(fcChanNames, eeg.chanNames);

figure();
topoplot(nanmean(postRespMeans(:,~ppsToExclude),2),...
    eeg.chanlocs, 'electrodes','off','colormap',cmap,'maplimits',mapLims,...
    'emarker',{'.','k',10,1}, 'emarker2', {fcChanInds, '.','k',10,1});
colorbar

% pick 5 electrode sites of the positive signal centroparietally
if useCSD
    chans5 = {'C23','C22','C1','D1','A1'};
else
    chans5 = {'A2','A3','A4','D15','D16','B2'}; % Wouter's ones
end

[~, chans5Inds] = ismember(chans5, eeg.chanNames);


%% topoplot per person

figure();
for iPP = 1:fileInfo.nPP
    if 1%~ppsToExclude(iPP)
        subplot(5,6,iPP);
        topoplot(postRespMeans(:,iPP),...
            eeg.chanlocs, 'electrodes','off','colormap',cmap,...
            'emarker',{'.','k',10,1}, 'emarker2', {chans5Inds, '.','k',10,1});
        title(iPP);
    end
end
        
%% pick maximal in that per person

keyboard; %%%%%% edit these picked channels

fcs = postRespMeans(chans5Inds,:); % get mean voltages per channel per pp   

% pick maximal amplitude surrounding response execution    
[m,j] = max(fcs);

fcChans = chans5(j); % get name
fcChanInds = chans5Inds(j); % index
    

%% load that data again and store fc chans, remove outliers
stimWin = [-500 1200];
stimWinInds = isBetween(eeg.epochTimes, stimWin);
eeg.stimTimes = eeg.epochTimes(stimWinInds);

fc = NaN(fileInfo.nPP, length(eeg.respTimes), fileInfo.maxTr);
fcStim = NaN(fileInfo.nPP, sum(stimWinInds), fileInfo.maxTr);

eeg.confRespTimes = ((-512:52) ./ 512) .* 1000; % times
confWin = minMax(eeg.confRespTimes,2);
% fcConfResp = NaN(fileInfo.nPP, length(eeg.confRespTimes), fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_resp']),'respErp');
    
    % store that channel as fc
    fc(iPP, :, :) = permute(data.respErp(fcChanInds(iPP), :, :),[4,2,3,1]);
    
    fc(iPP,:,isFlagged(:,iPP)) = NaN; % remove flagged trials before mean calc
    
    % also get stimlocked- outlier remove by resplocked
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpSuffix]),'erp');

    % get fc channel
    fcStim(iPP, :, :) = permute(data.erp(fcChanInds(iPP), stimWinInds, :),[4,2,3,1]);

    % remove flagged
    fcStim(iPP,:,isFlagged(:,iPP)) = NaN;

%     % and confResp locked
%     data = load(fullfile(fileInfo.respFolder, [fileInfo.ppID{iPP} '_confResp']),'respErp');
%     
%     % get fc channel
%     fcConfResp(iPP, :, :) = permute(data.respErp(fcChanInds(iPP), :, :),[4,2,3,1]);
% 
%     % remove flagged
%     fcConfResp(iPP,:,isFlagged(:,iPP)) = NaN;

end


%% save a filtered copy 

loPass = 10; %Hz

fcFilt = filterCPP(fc, loPass, eeg, fileInfo);
fcStimFilt = filterCPP(fcStim, loPass, eeg, fileInfo);
% fcConfRespFilt = filterCPP(fcConfResp, loPass, eeg, fileInfo);

%% save

save(fullfile(outFolder, saveName), ...
    'fc','fcChans','fcChanInds','chans5','postRespWindow','postRespInds','fileInfo','eeg',...
    'fcStim','stimWin','fcFilt','fcStimFilt');%,...
%     'fcConfResp','fcConfRespFilt','confWin');

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