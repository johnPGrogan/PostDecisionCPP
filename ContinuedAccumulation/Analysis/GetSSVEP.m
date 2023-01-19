% GetSSVEP
% load up each individual file
% extract ssvep (at Oz elec)
% save that

clc; clear all; close all;

load('./Saves/ExtractEpochs.mat')
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');
load('./Saves/BehDataLoad.mat','behData'); % for l/r

useCSD = 1; % CSD or voltage
doBaseline = 1; % SSVEP already present in baseline period

if useCSD
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT_CSD';
    load('./Saves/DoFFTCSD.mat'); % get FFT info
    saveName = 'GetSSVEPCSD.mat';

else
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT';
    load('./Saves/DoFFT.mat'); % get FFT info
    saveName = 'GetSSVEP.mat';

end

windows.stim = stimWindows;
windows.resp = respWindows;
nTimes.stim = length(stimWindows);
nTimes.resp = length(respWindows);

lockNames ={'stim','resp'};
%% load up all - get the grand average for good trials

for i = 1:2
    ssvep.(lockNames{i}) = NaN(fileInfo.nPP, nTimes.(lockNames{i}), fileInfo.maxTr);
end
for iPP = 1:fileInfo.nPP

    disp([fileInfo.ppID{iPP} '...'])

    % get resplocked data
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'ssvep');

    ssvep.stim(iPP,:,:) = data.ssvep;


    % get resplocked data
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'ssvep');

    ssvep.resp(iPP,:,:) = data.ssvep;
end

%% remove flagged

for i = 1:2
    ssvep.(lockNames{i})(repmat(permute(isFlagged,[2,3,1]),1,nTimes.(lockNames{i}),1)) = NaN;
end

%% baseline?
baselineInds = find(windows.stim <= -250,1,'last');%isBetween(windows.stim, [-150 0]);
baseline = nanmean(ssvep.stim(:, baselineInds,:),2);
    
if doBaseline
    
    for i = 1:2
        ssvep.(lockNames{i}) = log(ssvep.(lockNames{i}) ./ baseline);
    end


end
%% save'
save(fullfile('./Saves/', saveName), 'ssvep','windows','nTimes','lockNames',...
    'baseline','doBaseline','baselineInds');

