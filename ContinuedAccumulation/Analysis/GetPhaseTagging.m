% GetPhaseTagging
% Load up just Oz channel from everyone (or chans picked from topog)
% do convolution
% project into the dominant phase + flip
% save the traces, to be analysed separately
% takes about 40mins if loadUp=0, otherwise ~1mins

clc; clear; close all;

%% settings

useCSD = 1;
usePickedChans = 1; % use maximal chans in 200-300ms window?
loadUp = 1; % 1 = load up previously loaded channels, 0 = re-get them

outFolder = './Saves';
load('./Saves/ExtractEpochs.mat');
load('./Saves/BehDataLoad.mat','behData');
load('./Saves/NoiseData.mat','firstTilt');
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');

if useCSD
    fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD\'; % where to load from
    intpName = '_intp_csd.mat';
    loadName = fullfile(outFolder, 'OzData_CSD.mat');
    saveName = fullfile(outFolder, 'GetPhaseTagging_CSD2.mat');
    topoName = fullfile(outFolder, 'PhaseTopoAnalysis_CSD.mat');
else
    intpName = '_intp.mat';
    loadName = fullfile('OzData.mat');
    saveName = fullfile(outFolder, 'GetPhaseTagging.mat');
    topoName = fullfile(outFolder, 'PhaseTopoAnalysis.mat');
end


%% 

% get whether the first tilt was the bright or dim one
brightFirst = double(firstTilt == behData.corrLR);
brightFirst(isnan(behData.corrLR)) = NaN;


%% times to use as reference phase

refTimes = [100 400]; % [stim resp]

% set times for epoching
stimLims = [-200 1500];
stimTimes = eeg.epochTimes(isBetween(eeg.epochTimes, stimLims));

% for resplocking
respLims = [-500 1250];
respWindowS = 513 + (-256:640); % 513 is time 0 (evOnset)
respTimes = ((-256:640) ./ 512) .* 1000; % times

% pre-resp locking (i.e. initial stimulus, but 0 is RT)
preRespLims = [-1500 100];
preRespWindowS = 513 + (-768:52);
preRespTimes = ((-768:52) ./ 512) .* 1000; % times


%% pick channels - Oz or picked ones from topographies

if usePickedChans % pick from topographies
    load(topoName, 'visChanInds');
    ozInds = visChanInds;
else
    ozInds = repmat(find(strcmp(eeg.chanNames, 'A23')), 1, fileInfo.nPP); % which is Oz channel
end

%% stim locked - load up Oz from each person

if exist(loadName, 'file') && loadUp
    d = load(loadName, 'ozData','ozInds','usePickedChans');
    if d.usePickedChans == usePickedChans && all(d.ozInds == ozInds)
        ozData = d.ozData;
    else
        error('loaded data does not match ozInds');
    end


else
    ozData = NaN(fileInfo.nPP, eeg.nSamples, fileInfo.maxTr);
    for iPP = 1:fileInfo.nPP
        disp(iPP);
        data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp'); % interpolated, ev-locked epochs

        ozData(iPP,:,:) = data.erp(ozInds(iPP),:,:);

    end
    save(loadName, 'ozData', 'ozInds','usePickedChans');

end

%% remove bad trials?

ozData(repmat(permute(isFlagged,[2,3,1]),1,eeg.nSamples)==1) = NaN;
ozData(ppsToExclude==1,:,:) = NaN;


%% simon's way - convolution on 15hz signal

% 171samples = 334ms, which is 5 cycles of 15Hz
fftlen = 171; % 273;  % in sample points. 164 samples is 320 ms, almost exactly 6 cycles of the flicker.
% Make the complex sinusoid. It's rather like a wavelet but is not tapered like normal wavelets are:

fv = 15; %Hz
wvlt = exp(1i*2*pi*fv*[1:fftlen]/eeg.fs);

% take entire epoch for now, so can resplock from this later
% now convolve this with all erps at all channels:
ssvepStim = NaN(fileInfo.nPP, eeg.nSamples, fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    for iTr = 1:fileInfo.maxTr
        ssvepStim(iPP,:,iTr) = conv(ozData(iPP,:,iTr),wvlt,'same'); % convolve
    end
end


%% project the phase into the phase when it is driven by contrast

% pick an average phase per person, at at time when it should be driven by
% contrast difference

% first split dim + bright trials
ssvepStimByBright = groupMeans(ssvepStim, 3, repmat(permute(brightFirst,[1,3,2]),1,size(ssvepStim,2)), 'dim'); % [pp time bright tr]

iT = find(eeg.epochTimes>refTimes(1),1); % pick a time when there is a contrast difference, and the traces show a strong effect
refPhaseStim = angle(nanmean(diff(ssvepStimByBright(:,iT,:,:),[],3),4)); % bright-dim, mean over trials and pps 

% make a reference wave, extending back and forwards from this time point
refWave1 = exp(1i*((eeg.epochTimes-eeg.epochTimes(iT))/1000*fv*2*pi+refPhaseStim));

projWaveStim = ssvepStim ./ refWave1; % project into the reference phase angle
% real positive means towards that angle, negative means away from

% dimFirst is opposite this angle, so flip that to align them
projWaveStim(repmat(permute(brightFirst,[1,3,2]),1,size(projWaveStim,2))==0) = -projWaveStim(repmat(permute(brightFirst,[1,3,2]),1,size(projWaveStim,2))==0);

%% re-align this to have zero be the response time
% so can look at pre-response SSVEP to initial stimulus

projWavePreResp =NaN(fileInfo.nPP, length(preRespTimes), fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    disp(iPP);
    for i = 1:fileInfo.maxTr
        if ~isnan(behData.RS(iPP,i))
            windowInds = behData.RS(iPP,i) + preRespWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                projWavePreResp(iPP,:,i) = projWaveStim(iPP, windowInds, i); % whole epoch is used
            else
                % if the RT is < 250ms then will have to nanpad?
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                projWavePreResp(iPP,:,i) = [NaN(1,length(preRespTimes)-length(windowInds))  projWaveStim(iPP,windowInds,i)];
            end
        end
    end
end

%% trim the stim times now 

projWaveStim = projWaveStim(:, isBetween(eeg.epochTimes, stimLims), :);

%% do resplocked also

% re-resplock, quicker than loading
respErp = NaN(fileInfo.nPP, length(respTimes), fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    disp(iPP);
    for i = 1:fileInfo.maxTr
        if ~isnan(behData.RS(iPP,i))
            windowInds = behData.RS(iPP,i) + respWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                respErp(iPP,:,i) = ozData(iPP, windowInds, i); % whole epoch is used
            else
                % if the RT is < 250ms then will have to nanpad?
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                respErp(iPP,:,i) = [NaN(1,length(respTimes)-length(windowInds))  ozData(iPP,windowInds,i)];
            end
        end
    end
end


%% do convolution 

% now convolve this with all erps at all channels:
ssvepResp = NaN(fileInfo.nPP, length(respTimes), fileInfo.maxTr);

for iPP = 1:fileInfo.nPP
    for iTr = 1:fileInfo.maxTr
        ssvepResp(iPP,:,iTr) = conv(respErp(iPP,:,iTr),wvlt,'same');
    end
end


%% get a reference phase - after rotating dimFirst by 180 to align

% split dim + bright first
ssvepRespByBright = groupMeans(ssvepResp,3,repmat(permute(brightFirst,[1,3,2]),1,size(ssvepResp,2)),'dim'); % [pp time bright tr]

iT = find(respTimes>refTimes(2),1);
refPhaseResp = angle(nanmean(diff(ssvepRespByBright(:,iT,:,:),[],3),4)); % bright-dim, mean over trials and pps 

% make reference wave
refWaveResp = exp(1i*((respTimes-respTimes(iT))/1000*fv*2*pi+refPhaseResp));

projWaveResp = ssvepResp ./ refWaveResp; % project

% invert dimFirst
projWaveResp(repmat(permute(brightFirst,[1,3,2]),1,size(projWaveResp,2))==0) = -projWaveResp(repmat(permute(brightFirst,[1,3,2]),1,size(projWaveResp,2))==0) ;

%% store into a structure to save

phase.stim = real(projWaveStim); % just keep real for now
phase.resp = real(projWaveResp);
phase.preResp = real(projWavePreResp);
phaseTimes.stim = stimTimes;
phaseTimes.resp = respTimes;
phaseTimes.preResp = preRespTimes;
nTimes = structfun(@length,  phaseTimes,'UniformOutput',0);
phaseNames = {'stim','preResp','resp'};
nPhase = length(phaseNames);



%% save

save(saveName, 'phase','phaseTimes','nTimes','phaseNames', 'nPhase',...
    'brightFirst','usePickedChans','ozInds');