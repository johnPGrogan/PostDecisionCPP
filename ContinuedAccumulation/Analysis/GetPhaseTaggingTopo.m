% GetPhaseTaggingTopo
% calculate the phase tagging metric on each channel per person
% pick best channel per person
% save the topographies

clc; clear; close all;

%% plot 30Hz signal phase on trials

useCSD = 1;

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat');
load('./Saves/BehDataLoad.mat','behData');
load('./Saves/NoiseData.mat','firstTilt');
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');

if useCSD
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD\'; % where to load from
    fileInfo.phaseTopoFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\PhaseTopoCSD\'; % where to load from
    intpName = '_intp_csd.mat';
    saveName = fullfile(outFolder, 'GetPhaseTaggingTopo_CSD.mat');
else
    fileInfo.phaseTopoFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\PhaseTopo\'; % where to load from
    intpName = '_intp.mat';
    saveName = fullfile(outFolder, 'GetPhaseTaggingTopo.mat');
end


%% 

% get whether the first tilt was the bright or dim one
brightFirst = double(firstTilt == behData.corrLR);
brightFirst(isnan(behData.corrLR)) = NaN;


refTimes = [100 400]; % [stim resp]

ozInd = find(strcmp(eeg.chanNames, 'A23')); % which is Oz channel

% 171samples = 334ms, which is 5 cycles of 15Hz
fftlen = 171; % 273;  % in sample points. 164 samples is 320 ms, almost exactly 6 cycles of the flicker.
% Make the complex sinusoid. It's rather like a wavelet but is not tapered like normal wavelets are:

fv = 15; %Hz
wvlt = exp(1i*2*pi*fv*[1:fftlen]/eeg.fs); % make wavelet


%% trim epoch

tLims = [-200 1500];
tInds = isBetween(eeg.epochTimes, tLims);
stimTimes = eeg.epochTimes(tInds);

iT = find(eeg.epochTimes>refTimes(1),1); % pick a time when there is a contrast difference, and the traces show a strong effect

stimWins = -200:100:1500;
stimWinInds = sum(stimTimes > stimWins'); % 1:21
stimWinInds(1) = 1; % include -1000 in this

% pre-resp locking (from stim ssvep)
preRespWins = -1500:100:100;
preRespLims = [-1500 100];
preRespWindowS = 513 + (-768:52);
preRespTimes = ((-768:52) ./ 512) .* 1000; % times
preRespWinInds = sum(preRespTimes > preRespWins'); % 1:21
preRespWinInds(1) = 1; % include -1000 in this


% for resplocking
respWins = -500:100:1100;
respLims = [-500 1250];
respWindowS = 513 + (-256:640); % 513 is time 0 (evOnset)
respTimes = ((-256:640) ./ 512) .* 1000; % times
respWinInds = sum(respTimes > respWins'); % 1:21
respWinInds(1) = 1; % include -1000 in this

iTResp = find(respTimes>refTimes(2),1); % pick a time when there is a contrast difference, and the traces show a strong effect


%% convolve signal with data, per person

% now convolve this with all erps at all channels:
projWaveStimMean = NaN(fileInfo.nPP, sum(tInds), eeg.nChans);
projWaveStimMeanWindows = NaN(fileInfo.nPP, max(stimWinInds), eeg.nChans, fileInfo.maxTr);
projWavePreRespMean = NaN(fileInfo.nPP, length(preRespTimes), eeg.nChans);
projWavePreRespMeanWindows = NaN(fileInfo.nPP, max(preRespWinInds), eeg.nChans, fileInfo.maxTr);
projWaveRespMean = NaN(fileInfo.nPP, length(respTimes), eeg.nChans);
projWaveRespMeanWindows = NaN(fileInfo.nPP, max(respWinInds), eeg.nChans, fileInfo.maxTr);


for iPP = 1:fileInfo.nPP
    disp(iPP);
    tic;
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp'); % interpolated, ev-locked epochs
    
    ssvepStim = NaN(eeg.nChans, eeg.nSamples, fileInfo.maxTr);
    for iCh = 1:eeg.nChans
        for iTr = 1:fileInfo.maxTr
            ssvepStim(iCh,:,iTr) = conv(data.erp(iCh,:,iTr),wvlt,'same'); % convolve
        end
    end

    % rotate + project
    % first split dim + bright trials
    ssvepStimByBright = permute(groupMeans(ssvepStim, 3, brightFirst(iPP,:), 'dim'),[4,1,3,2]); % [pp time bright tr]  
    refPhaseStim = angle(nanmean(diff(ssvepStimByBright(:,iT,:,:),[],3),4)); % bright-dim, mean over trials and pps 
    
    % make a reference wave, extending back and forwards from this time point
    refWave1 = exp(1i*((eeg.epochTimes-eeg.epochTimes(iT))/1000*fv*2*pi+refPhaseStim));
    
    projWaveStim = ssvepStim ./ refWave1; % project into the reference phase angle
    % real positive means towards that angle, negative means away from
    
    % dimFirst is opposite this angle, so flip that to align them
    projWaveStim(:,:,brightFirst(iPP,:)==0) = -projWaveStim(:,:,brightFirst(iPP,:)==0);

    % remove flagged trials
    projWaveStim(:,:,isFlagged(:,iPP)==1) = NaN;

    
    %%%%% also get resp-locked, from this stim   
    for iTr = 1:fileInfo.maxTr
        if ~isnan(behData.RS(iPP,iTr))
            windowInds = behData.RS(iPP,iTr) + preRespWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                for iCh = 1:eeg.nChans
                    projWavePreResp(iCh,:,iTr) = projWaveStim(iCh, windowInds, iTr); % whole epoch is used
                end
            else
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                for iCh = 1:eeg.nChans
                    projWavePreResp(iCh,:,iTr) = [NaN(1,length(preRespTimes)-length(windowInds))  projWaveStim(iCh,windowInds,iTr)];
                end
            end
        end
    end    
    
    % can now trim stim-locked
    projWaveStim = projWaveStim(:,tInds,:);

    % store mean across trials
    projWaveStimMean(iPP,:,:) = nanmean(real(projWaveStim),3)';

    % also take 20 time samples
    projWaveStimMeanWindows(iPP,:,:,:) = permute(nanmean(groupMeans(real(projWaveStim),2,stimWinInds,'dim'),1),[1,3,4,2]);

    % get means across trials and within windows
    projWavePreRespMean(iPP,:,:) = nanmean(real(projWavePreResp),3)';
    projWavePreRespMeanWindows(iPP,:,:,:) = permute(nanmean(groupMeans(real(projWavePreResp),2,preRespWinInds,'dim'),1),[1,3,4,2]);



    %%%%%% now resplock on ssvep, and calc for post
    ssvepResp = NaN(fileInfo.nPP, length(respTimes), fileInfo.maxTr);    
    for iTr = 1:fileInfo.maxTr
        if ~isnan(behData.RS(iPP,iTr))
            windowInds = behData.RS(iPP,iTr) + respWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                 for iCh = 1:eeg.nChans
                    ssvepResp(iCh,:,iTr) = ssvepStim(iCh, windowInds, iTr); % whole epoch is used
                 end
            else
                % if the RT is < 250ms then will have to nanpad?
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                 for iCh = 1:eeg.nChans
                    ssvepResp(iCh,:,iTr) = [NaN(1,length(respTimes)-length(windowInds))  ssvepStim(iCh,windowInds,iTr)];
                end
            end
        end
    end
    
    ssvepRespByBright = permute(groupMeans(ssvepResp, 3, brightFirst(iPP,:), 'dim'),[4,1,3,2]); % [pp time bright tr]  
    refPhaseResp = angle(nanmean(diff(ssvepRespByBright(:,iTResp,:,:),[],3),4)); % bright-dim, mean over trials and pps 
    
    % make a reference wave, extending back and forwards from this time point
    refWave2 = exp(1i*((respTimes-respTimes(iT))/1000*fv*2*pi+refPhaseResp));
    
    projWaveResp = ssvepResp ./ refWave2; % project into the reference phase angle
    
    % dimFirst is opposite this angle, so flip that to align them
    projWaveResp(:,:,brightFirst(iPP,:)==0) = -projWaveResp(:,:,brightFirst(iPP,:)==0);

    % store mean across trials
    projWaveRespMean(iPP,:,:) = nanmean(real(projWaveResp),3)';

    % also take 20 time samples
    projWaveRespMeanWindows(iPP,:,:,:) = permute(nanmean(groupMeans(real(projWaveResp),2,respWinInds,'dim'),1),[1,3,4,2]);

end

%% save

save(saveName, 'projWaveStimMean','projWaveStimMeanWindows','stimWins',...
    'projWaveRespMean','projWaveRespMeanWindows','respWins',...
    'projWavePreRespMean','projWavePreRespMeanWindows','preRespWins',...
    'respTimes','stimTimes','preRespTimes');
