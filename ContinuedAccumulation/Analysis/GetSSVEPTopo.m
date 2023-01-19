% GetSSVEPTopo
% do FFT on all channels (stim + resp locked), and get 30Hz signal
% (normalised to adjacents)

clc; clear; close all;

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat');
load('./Saves/FlagArtefacts3.mat','isFlagged');

useCSD = 1;
doBaseline = 1;

if useCSD
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT_CSD\';
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD_ref_noBS\';
    saveName1 = 'GetSSVEPTopoCSD.mat';
    intpName = '_intp_csd.mat';

else
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT';
    saveName1 = 'GetSSVEPTopo.mat';
    intpName = '_intp.mat';

end
load('./Saves/BehDataLoad.mat','behData'); % for RS

%% SET

ssvepHz = 30; % 30 or 15hz signal
fs=eeg.fs;
% STFT parameters
nFreqs = round( (fs/ssvepHz) * 15); % Window of how many sample points? If there is an SSVEP involved, whether or not you are interested in analyzing it, it is good to have all power related to the SSVEP isolated in a single frequency bin. This happens when you choose a window length that is an integer number of SSVEP cycles. In this example, there is a 17Hz SSVEP and I've made the widnow length = 5 cycles
freqs = ((0:nFreqs-1)*fs) / nFreqs; % frequency scale, given window length (remember resolution = 1/window-duration)
stimWindows = -300:100:1500; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time

windowSize = length(freqs); % samples; 500ms; 4 cycles of 8Hz

ssvepFreqInds = find(any(freqs == (ssvepHz .* [14/15 1 16/15]'))); % will get 14:16 or 28:2:32

respWindowS = 513 + (-640:640); % 513 is time 0 (evOnset)
respTimes = ((-640:640) ./ 512) .* 1000; % times
respWindows = -1000:100:1000;

ssvepTopo = zeros(30,128,length(stimWindows),1280);
ssvepTopoResp = zeros(30,128,length(respWindows),1280);

%%
for iPP = 1:fileInfo.nPP
    disp(iPP);
    tic;
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp'); % interpolated, ev-locked epochs
    
    % compute short time fourier transform (STFT) for each single trial:
    
    ssvepTopo(iPP,:,:,:) = myFFT(data.erp, nFreqs, stimWindows, ssvepFreqInds, eeg.epochTimes);
    
        
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
    ssvepTopoResp(iPP,:,:,:) = myFFT(respErp, nFreqs, respWindows, ssvepFreqInds, respTimes);
    
    t = toc % ~3 mins per pp
end

%% remove flagged

ssvepTopo(repmat(permute(isFlagged,[2,3,4,1]),1,128,size(ssvepTopo,3),1)) = NaN;
ssvepTopoResp(repmat(permute(isFlagged,[2,3,4,1]),1,128,size(ssvepTopoResp,3),1)) = NaN;

%% baseline?
baselineInds = find(stimWindows <= -250,1,'last');%isBetween(windows.stim, [-150 0]);
baseline = nanmean(ssvepTopo(:,:, baselineInds,:),3);
    
if doBaseline
    
    ssvepTopo = log(ssvepTopo ./ baseline);
    ssvepTopoResp = log(ssvepTopoResp ./ baseline);


end


%% save
save(fullfile(outFolder, saveName1), '-v7.3','ssvepTopo','ssvepTopoResp','ssvepFreqInds','eeg','baselineInds','baseline','doBaseline','stimWindows','respWindows','freqs');


%%%%%%%%%%%%%%%%%%%%%
%%

function ssvepTopo = myFFT(erp, nFreqs, windowCentres, ssvepFreqInds, epochTimes)
% do FFT on each time window, get mean mu-beta, and ssvep (normalised)

[nChans,~,nTr] = size(erp);
ssvepTopo = NaN(128, length(windowCentres), nTr);

for iT = 1:length(windowCentres) % for each STFT timepoint (window centre) we'll compute the FFT for all trials at once
    
    [~,samp] = min(abs(epochTimes-windowCentres(iT))); % find the sample point in the ERP epoch corresponding to the centre of the current FFT window
    spec = abs( fft( erp(:,samp-round(nFreqs/2)+(1:nFreqs),:) ,[],2) ) ./ (nFreqs/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
    
    
    % also get the SSVEP - 30Hz at Oz, divided by 28 + 32hz
    ssvepTopo(:,iT,:) = spec(1:128, ssvepFreqInds(2), :) ./ nansum(spec(1:128,ssvepFreqInds([1 3]),:),2);
    
end

end