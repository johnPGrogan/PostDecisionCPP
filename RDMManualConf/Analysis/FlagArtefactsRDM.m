% FlagArtefactsRDM
% load up all front + veog from all people, save. load that to speed up
% artefact stuff

clc; clear all; close all; 
outFolder = './Saves';
loadName = fullfile(outFolder, 'artStuff.mat');
load(fullfile(outFolder, 'ExtractEpochs.mat'));
load(fullfile(outFolder, 'BehDataLoad.mat'),'behData')

fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\Interp\'; % my channels

% frontal channels - some people had very noisy VEOG (probably bad contact)
% so I am using the frontal channel above the eyes to check them
frontInds = find(ismember(eeg.chanNames, {'C29','C17','C16'}));

% some people had different channel orders, so these are the VEOG channel
% indices per person
veogChans = [3 4; 1 5; 3 4; 3 4; 3 4;
    3 4; 3 4; 3 4; 3 4; 3 4;
    1 4; 3 4; 3 4; 3 4; 3 4;
    1 4; 3 4; 3 4; 3 4; 3 4;
    3 4; 3 4; 1 4; 3 4; 3 4;
    3 4; 3 4;
    ];

art.rtLimsSamp = [51 768]; % in samples
art.rtLimsTime = art.rtLimsSamp / 512 * 1000; % in ms relative to evidence onset

eeg.veogChans = [129 130]; % 131 is upper and 132 is lower - BUT CHECK THAT THESE ARE ACTUALLY YOUR VEOG CHANNELS - DIFFERENT PHD STUDENTS USE DIFFERENT CONVENTIONS!
art.blinkTh = 200; % threshold for detecting a blink in VEOG upper-lower difference signal (peak2peak moving)
art.blinkWindow = 200; % window (ms) for blink detection
art.blinkStep = 100; % stepsize (ms)
art.blinkTh2 = 100; % for eeglab func

art.saccTh = 50; % thresh for step-change (saccade detection)
art.saccWindow = 100; % ms
art.saccStep = 30; %ms

art.moveTh = [-200 200]; % [min max] change for detecting movement or other artefacts
art.artifTh = 100; % this is for SCALP channels - it's high because it is in CSD units, which take bigger values
art.chans = 1:128; % Artifact rejection channels - which electrodes should we reject based on?
art.windowTime = [-200 2250]; % in msec, the limits of the time window relative to the event in which you will check for an artifact
% will need to change this to run up until final cue/resp?
art.windowInds= isBetween(eeg.epochTimes, art.windowTime); % indices of the timepoints for the artifact check window whose limits (min and max) are defined above

art.goodThresh = fileInfo.maxTr/2; % need at least 50% good trials?

isArt = NaN(eeg.nChans, fileInfo.maxTr, fileInfo.nPP);
isMovement = NaN(fileInfo.maxTr, fileInfo.nPP);
veog = NaN(1, sum(art.windowInds), fileInfo.maxTr,fileInfo.nPP);
frontChans = NaN(3, sum(art.windowInds), fileInfo.maxTr,fileInfo.nPP);


%% load or get data

if exist(loadName,'file')
    load(loadName)
else
    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])
        load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} '_intp']),'erp','stimDur','RS')


        erpArt = NaN(eeg.nChansTot, sum(art.windowInds), fileInfo.maxTr);
        for iTr = 1:fileInfo.maxTr
            if ~isnan(RS(iTr))

                % take from -200:response
                respLims = [art.windowTime(1) stimDur(iTr)+ (RS(iTr)/eeg.fs*1000)]; % can add time on here, <200ms should not exclude too many more trials

                respInds = isBetween(eeg.epochTimes, respLims);
                erpArt(:,1:sum(respInds),iTr) = erp(:, respInds, iTr); % get resp epoch
            end
        end

        %% artefact reject?

        isArt(:,:,iPP) = sq(max(abs(erpArt(1:eeg.nChans,:,:)),[],2) > art.artifTh);

        veog(1,1:size(erpArt,2),:,iPP) = diff(erpArt(128+veogChans(iPP,:),:,:),[],1);

        frontChans(:,1:size(erpArt,2),:,iPP) = erpArt(frontInds,:,:);

%         isMovement(:,iPP) = sq(any(myArtextval(erpArt(1:eeg.nChans,:,:), art.moveTh))); %[chan 1 tr]. pop_artextval - min/max in window

    end

    save(loadName,'veog','frontChans','isMovement','isArt');
end


%% now flag
isBlink = sq(max(veog,[],2)-min(veog,[],2)) > art.blinkTh;

%% try other approaches
isBlink3 = sq(max(abs(veog),[],2) > art.blinkTh); % this will give a logical 1 if this threshold is exceeded on a certain trial, 0 otherwise
% isBlink2 = sq(myArtmwppth(permute(veog,[4,2,3,1]), art.blinkTh2, art.blinkWindow, art.blinkStep, eeg.fs))'; %[tr 1]
% isSaccade = sq(myArtstep(permute(veog,[4,2,3,1]), art.saccTh, art.saccWindow, art.saccStep, eeg.fs))';

isMissing = isnan(behData.RS)';
isBadRT = ~isBetween(behData.RS, art.rtLimsSamp)';

isFrontBlink = sq((max(frontChans(2,:,:,:),[],2) - min(frontChans(2,:,:,:),[],2)) > 100); % just C17 (AFz)
isFrontBlink2 = sq(myArtmwppth(permute(frontChans(2,:,:,:),[4,2,3,1]), 50, art.blinkWindow, art.blinkStep, eeg.fs))';
% isFrontSaccade = sq(myArtstep(permute(frontChans(2,:,:,:),[4,2,3,1]), art.saccTh, art.saccWindow, art.saccStep, eeg.fs))';

%% draw figs
for iPP = 1%:fileInfo.nPP

    %%
    figure(1);clf;
    bar(nansum(isArt(:,:,iPP),2)); xlabel('chan');ylabel('num arts');
    ylim([0 fileInfo.maxTr]);
    
    %% plot some stuff
    
    figure(2);clf;
    nTr=20;
    isOut = isBlink; %isBlink | sq(any(isArt,1)); %sq(any(isArt,1));
    t = {'accept','reject'};
    for i = 1:2
        subplot(2,1,i);
        plot(sq(veog(1,:,find(isOut(:,iPP) == i-1, nTr),iPP)));
        ylim([-500 500]);
        ylabel('\muV'); xlabel('time');
        %         xline(-200); xline(2500); yline(0);
        title(t{i});
    end
    
    drawnow;
    
    %% look at frontal channels only

    % plot some trials, then split by good and bad
    nTr = 20;
    figure(3);clf
    for i = 1:3
        subplot(3,3,i*3-2);
        plot(sq(frontChans(i,:,find(~isnan(frontChans(i,200,:,iPP)),nTr),iPP)));
        ylabel(eeg.chanNames(frontInds(i)));
        title('all')
    
        subplot(3,3,i*3-1);
        plot(sq(frontChans(i,:,find(~isFrontBlink2(:,iPP),nTr),iPP)),'b');
        title('include')
    
        subplot(3,3,i*3);
        plot(sq(frontChans(i,:,find(isFrontBlink2(:,iPP),nTr),iPP)),'r');
        title('reject')
    end
    
    %% plot the ones that differ between veog/frontals
   
    figure(4); clf
    d(:,1) = isBlink(:,iPP) & ~isFrontBlink(:,iPP);
    d(:,2) = ~isBlink(:,iPP) & isFrontBlink(:,iPP);
    
    subplot(2,2,1);
    plot(sq(veog(1,:,find(d(:,1),nTr),iPP)));
    ylabel('VEOG')
    title('veog blink, eeg missed')
    
    subplot(2,2,2);
    plot(sq(veog(1,:,find(d(:,2),nTr),iPP)));
    ylabel('VEOG')
    title('EEG blink, VEOG missed')
    
    subplot(2,2,3);
    plot(sq(frontChans(2,:,find(d(:,1),nTr),iPP)));
    ylabel('EEG')
    %     title('veog blink, eeg missed')
    
    subplot(2,2,4);
    plot(sq(frontChans(2,:,find(d(:,2),nTr),iPP)));
    ylabel('EEG')
    %     title('EEG blink, VEOG missed')
    
    %% error bar
    
    figure(5); clf
    subplot(2,2,1)
    h = errorBarPlot(permute(groupMeans(veog(1,:,:,iPP),3,isBlink(:,iPP),'dim'),[2,1,3]),'area',1);
    yline(0);
    legend([h{:,1}], {'accept', 'reject'}, 'Location','Best');
    ylabel('veog');
    title('veog blink')
    
    subplot(2,2,2)
    h = errorBarPlot(permute(groupMeans(veog(1,:,:,iPP),3,isFrontBlink(:,iPP),'dim'),[2,1,3]),'area',1);
    yline(0);
    ylabel('veog');
    title('EEG blink')
    
    subplot(2,2,3);
    h = errorBarPlot(permute(groupMeans(frontChans(2,:,:,iPP),3,isBlink(:,iPP),'dim'),[2,1,3]),'area',1);
    yline(0);
    ylabel('eeg')
    
    subplot(2,2,4);
    h = errorBarPlot(permute(groupMeans(frontChans(2,:,:,iPP),3,isFrontBlink(:,iPP),'dim'),[2,1,3]),'area',1);
    yline(0);
    ylabel('eeg')
    
%     makeSubplotScalesEqual(2,2);
    
    %%
end

%% pick different settings for diff people

useFrontBlink2 = [1 2 13 21 24]; % use frontal channels for these people as their VEOG was much worse
isFrontBlink2(:,13) = sq(myArtmwppth(frontChans(2,:,:,13), 80, art.blinkWindow, art.blinkStep, eeg.fs))'; % use bigger thresh here

blinksToUse = isBlink;
blinksToUse(:, useFrontBlink2) = isFrontBlink2(:,useFrontBlink2);

% should I be removing badRT here?
isFlagged = blinksToUse | sq(any(isArt==1,1)) ; % 1 = flagged for removal [tr pp]
isGood = isFlagged==0;

nGoodTrials = nansum(isGood);
ppsToExclude = nGoodTrials < art.goodThresh ; %4th dim


%% are they missing at random?

load(fullfile(outFolder, 'BehDataLoad.mat'));

% regression with means rejected by condition (+ acc?)

factors = {'corrLR','stimDur','respLR','acc','confResp'};
for i = 1:length(factors)
    means = groupMeans(isGood', 2, behData.(factors{i})); %[pp cond (interrupt/cont)]
    anovas.(factors{i}) = rmanova(means, {'pp','cond'},'categorical',2);
    disp(factors{i});
    disp(anovas.(factors{i}));
end



%% save

save(fullfile(outFolder, 'FlagArtefacts.mat'), ...
    'isFlagged','isBlink','isArt','isMovement','isBlink2',...
    'ppsToExclude','isBadRT', 'art','isSaccade','isBlink3','isMissing','nGoodTrials','isGood','blinksToUse','useFrontBlink2','isFrontBlink2');