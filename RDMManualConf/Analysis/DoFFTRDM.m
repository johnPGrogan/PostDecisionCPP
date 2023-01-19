% DoFFTRDM
% run per person, get beta, alpha, theta bands out
% save per peron
for useCSD = 1 %1:2 will do both


clc; close all; clearvars -except useCSD;
load('./Saves/ExtractEpochs.mat');
if useCSD
    fileInfo.fftFolder = 'D:\TCD\Projects\RDMManualConf\Data\FFT_CSD\';
    fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\CSD\';
    saveName1 = 'DoFFTCSD.mat';
    intpName = '_intp_csd.mat';

else
    fileInfo.fftFolder = 'D:\TCD\Projects\RDMManualConf\Data\FFT';
    saveName1 = 'DoFFT.mat';
    intpName = '_intp.mat';

end
load('./Saves/BehDataLoad.mat','behData'); % for RS

%% SET

fs=eeg.fs;
% STFT parameters
nFreqs = round( (fs/30) * 15); % length of window = num freqs,
freqs = ((0:nFreqs-1)*fs) / nFreqs; % frequency scale, given window length (remember resolution = 1/window-duration)
windowSize = 256; % samples; 500ms; 4 cycles of 8Hz

freqBands = [15 30; 8 14; 5 7]; % [beta, alpha, theta]

% stim windows, must be windowSize from edges of epoch
stimWindows = -400:20:1500; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time

% respcuelocked - min respcue is 350ms, so 
respCueWindows = -800:20:500;
respCueWindowS = (-666:512) + 540; % to add on to evOnset time
respCueTimes = ((-666:512) ./ 512) .* 1000; % times

% resplocked
respWindows = -1000:20:1000; % so take window 500ms either side
respWindowS = (-768:768) + 540; % to add on to evOnset time
respTimes = ((-768:768) ./ 512) .* 1000; % times


save(fullfile('./Saves/',saveName1), 'fs','nFreqs','freqs','freqBands','windowSize',...
    'stimWindows','respTimes','respWindows','respCueWindows','respCueTimes')

for iPP = 1:fileInfo.nPP
if 1%~exist( fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim.mat']), 'file')
    disp(iPP);
    tic;
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp'); % interpolated, ev-locked epochs
    
    % compute short time fourier transform (STFT) for each single trial:
    STFT = myFFT(data.erp, stimWindows, eeg.epochTimes, freqs, freqBands);
    
    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'-v7.3',...
        'STFT','freqBands');
        
    %% also do resp locked
    % quicker to recalc this than load it up

    raw = load(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw.mat']),'stimDurSample'); % interpolated, ev-locked epochs

    respCueErp = NaN(eeg.nChansTot, length(respCueTimes), fileInfo.maxTr);
    % get from -1250:1250 around response cue
    for i = 1:fileInfo.maxTr
        if ~isnan(raw.stimDurSample(i))
            windowInds = raw.stimDurSample(i) + respCueWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                respCueErp(:,:,i) = data.erp(:, windowInds, i); % whole epoch is used
            else
                % if the RT is < 250ms then will have to nanpad?
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                respCueErp(:,:,i) = [NaN(eeg.nChansTot,length(respCueTimes)-length(windowInds))  data.erp(:,windowInds,i)];
            end
        end
    end
    % run SFTF on this
    STFT = myFFT(respCueErp, respCueWindows, respCueTimes, freqs, freqBands);
    % save
    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_respCue']),'-v7.3',...
        'STFT','freqBands');
    
    
    %% also do resp locked
    % quicker to recalc this than load it up

    respErp = NaN(eeg.nChansTot, length(respTimes), fileInfo.maxTr);
    % get from -1250:1250 around response
    for i = 1:fileInfo.maxTr
        if ~isnan(behData.RS(iPP,i))
            windowInds = behData.RS(iPP,i) + raw.stimDurSample(i) + respWindowS;
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
    STFT = myFFT(respErp, respWindows, respTimes, freqs, freqBands);
    
    % save
    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'-v7.3',...
        'STFT','freqBands');
    
    t = toc
end    
end

end

function STFT = myFFT(erp, windowCentres, epochTimes, freqs, freqBands)
% STFT = myFFT(erp, windowCentres, epochTimes, freqs, freqBands)
% do FFT on each time window, get abs power, take mean within frequency
% bands
% Inputs:
%   erp = [nChans, nTimes, nTrials]
%   windowCentres = time points of centres of windows
%   epochTimes = time points of erp [1 nTimes]
%   freqs = vector of frequencies desired, will be used to pick window
%       widths (1 sample per frequency)
%   freqBands = [nBands 2]: [min max] freqs of each band, 1 row per band
% 
% Outputs:
%   STFT = [nChans, nWindows, nTrials, nBands] mean abs power in each band

nFreqs = length(freqs);
nBands = size(freqBands,1);

bandInds = false(nFreqs, nBands);
for iB = 1:nBands
    bandInds(:,iB) = isBetween(freqs, freqBands(iB,:));
end

[nChans,~,nTr] = size(erp);
STFT = NaN(nChans, length(windowCentres), nTr, nBands); % preallocate

for iT = 1:length(windowCentres) % for each STFT timepoint (window centre) we'll compute the FFT for all trials at once
    
    [~,samp] = min(abs(epochTimes-windowCentres(iT))); % find the sample point in the ERP epoch corresponding to the centre of the current FFT window
    spec = abs( fft( erp(:,samp-round(nFreqs/2)+(1:nFreqs),:) ,[],2) ) ./ (nFreqs/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
    
    for iB = 1:nBands
        % Save the result for each trial, just like the matrix 'erp'
        STFT(:,iT,:,iB) = nanmean(spec(:,bandInds(:,iB),:),2);
    end
    
end

end