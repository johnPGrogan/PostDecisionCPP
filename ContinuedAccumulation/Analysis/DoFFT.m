% DoFFT - fast fourier transforms 

useCSD = 1;


clc; close all; clearvars -except useCSD;
load('./Saves/ExtractEpochs.mat');
if useCSD
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT_CSD\';
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD_ref_noBS\';
    saveName1 = 'DoFFTCSD.mat';
    intpName = '_intp_csd.mat';

else
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT';
    saveName1 = 'DoFFT.mat';
    intpName = '_intp.mat';

end
load('./Saves/BehDataLoad.mat','behData'); % for RS

%% SET FFT params

fs=eeg.fs;
% STFT parameters
nFreqs = round( (fs/30) * 15); % nFreqs is length of window, pick SSVEP freq multiples
freqs = ((0:nFreqs-1)*fs) / nFreqs; % frequency scale, given window length (remember resolution = 1/window-duration)
betaFreqInds = find( (freqs>=8 & freqs<=30) ); % beta power freqs

stimWindows = -250:20:1450; % centres of each window
windowSize = 256; % samples; 500ms; 4 cycles of 8Hz

ssvepFreqInds = find(isBetween(freqs, [28, 32])); % 30Hz SSVEP + adjacent frequencies
alphaFreqInds = find(isBetween(freqs, [8 14])); % alpha band
thetaFreqInds = find(isBetween(freqs, [5 7])); % theta band

respWindowS = 513 + (-640:640); % 513 is time 0 (evOnset). samples
respTimes = ((-640:640) ./ 512) .* 1000; % times
respWindows = -1000:20:1000; % ms

save(fullfile('./Saves/',saveName1), 'fs','nFreqs','freqs','betaFreqInds','stimWindows',...
    'windowSize','ssvepFreqInds','respWindowS','respTimes','respWindows')

for iPP = 1:fileInfo.nPP

if 1%isempty(whos('-file',fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim.mat']),'alpha')) % only do if not already done alpha
    disp(iPP);
    tic;
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp'); % interpolated, ev-locked epochs
    
    % compute short time fourier transform (STFT) for each single trial:
    
    [STFT, ssvep, alpha, theta] = myFFT(data.erp, nFreqs, stimWindows, betaFreqInds, ssvepFreqInds, eeg.epochTimes, alphaFreqInds, thetaFreqInds);

    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'-v7.3',...
        'STFT','ssvep','alpha','theta');
    
    
    %% also do resp locked
    % quicker to recalc this than load it up

    respErp = NaN(eeg.nChansTot, length(respTimes), fileInfo.maxTr);
    % get from -1250:1250 around response
    for i = 1:fileInfo.maxTr
        if ~isnan(behData.RS(iPP,i))
            windowInds = behData.RS(iPP,i) + respWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                respErp(:,:,i) = data.erp(:, windowInds, i); % whole epoch is used
            else
                % if the RT is < 250ms then will have to nanpad?
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                respErp(:,:,i) = [NaN(eeg.nChansTot,length(respTimes)-length(windowInds))  data.erp(:,windowInds,i)];
            end
        end
    end
    % run SFTF on this
    [STFT, ssvep, alpha, theta] = myFFT(respErp, nFreqs, respWindows, betaFreqInds, ssvepFreqInds, respTimes, alphaFreqInds, thetaFreqInds);
    
    % save
    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'-v7.3',...
        'STFT','ssvep','alpha','theta');
    
    t = toc
end    
end


function [STFT, ssvep, alpha,theta] = myFFT(erp, nFreqs, windowCentres, betaFreqInds, ssvepFreqInds, epochTimes, alphaFreqInds,thetaFreqInds)
% do FFT on each time window, get mean mu-beta, and ssvep (normalised), and
% theta and alpha if requested

[nChans,~,nTr] = size(erp);
STFT = NaN(nChans, length(windowCentres), nTr); % preallocate
ssvep = NaN(1,length(windowCentres), nTr);
if exist('alphaFreqInds', 'var') && ~isempty(alphaFreqInds)
    alpha = STFT;
else
    alpha = [];
end
if exist('thetaFreqInds', 'var') && ~isempty(thetaFreqInds)
    theta = STFT;
else
    theta = [];
end

for iT = 1:length(windowCentres) % for each STFT timepoint (window centre) 
    
    [~,samp] = min(abs(epochTimes-windowCentres(iT))); % find centrepoint of window
    spec = abs( fft( erp(:,samp-round(nFreqs/2)+(1:nFreqs),:) ,[],2) ) ./ (nFreqs/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
    
    % store mean beta
    STFT(:,iT,:) = nanmean(spec(:,betaFreqInds,:),2);
    
    % also get the SSVEP - 30Hz at Oz, divided by 28 + 32hz
    ssvep(1,iT,:) = spec(23, ssvepFreqInds(2), :) ./ nansum(spec(23,ssvepFreqInds([1 3]),:),2);

    if exist('alphaFreqInds', 'var') && ~isempty(alphaFreqInds)
        alpha(:,iT,:) = nanmean(spec(:,alphaFreqInds,:),2);
    end

    if exist('thetaFreqInds', 'var') && ~isempty(thetaFreqInds)
        theta(:,iT,:) = nanmean(spec(:,thetaFreqInds,:),2);
    end
    
end
end