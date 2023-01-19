% GetMuBetaTopo
% use picked channels from GetMuBeta.m
% load each indiv data file
% get means within each 100ms bin

%% set up
clc; clear all; close all;

load('./Saves/ExtractEpochs.mat')
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');
load('./Saves/BehDataLoad.mat','behData'); % for l/r

loadUp = 1; % load up PreRespMeanVoltage file?
useCSD = 1; % CSD or voltage

if useCSD
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT_CSD';
    load('./Saves/DoFFTCSD.mat'); % get FFT info
    loadName = 'PreRespMeanBetaCSD.mat';
    saveName = 'GetMuBetaTopoCSD.mat';
    mapLims = [-1.3 1.3];
else
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT';
    load('./Saves/DoFFT.mat'); % get FFT info
    loadName = 'PreRespMeanBeta.mat';
    saveName = 'GetMuBetaTopo.mat';
    mapLims = [-.1 .1];
end

%% windows

% respWindows
wins = -1000:100:1000;
winInds = sum(respWindows > wins'); % 1:21
winInds(1) = 1; % include -1000 in this
nWins = length(wins)-1;


%% 

betaTopoWins = NaN(fileInfo.nPP, fileInfo.maxTr, eeg.nChans, nWins); %[pp tr chans win]
% data.STFT is [hemi, t, tr] per pp

for iPP = 1:fileInfo.nPP
    disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'STFT');
    
    % split into bins, av within them. make into [pp tr chans wins]
    betaTopoWins(iPP,:,:,:) = permute(nanmean(groupMeans(data.STFT(1:eeg.nChans,:,:),2,winInds,'dim'),1),[1,2,4,3]);

    betaTopoWins(iPP,isFlagged(:,iPP),:,:) = NaN; % remove flagged trials
    
end



%% save

save(fullfile('./Saves',saveName),...
    'betaTopoWins','wins')


