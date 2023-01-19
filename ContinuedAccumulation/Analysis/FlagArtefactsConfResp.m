% function FlagArtefactsConfResp()
% flag on the conf-resp locked epochs

outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));

fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\InterpConfResp\'; % my channels

art.rtLimsSamp = [51 768]; % in samples
art.rtLimsTime = art.rtLimsSamp / 512 * 1000; % in ms relative to evidence onset

eeg.veogChans = [129 130]; % 131 is upper and 132 is lower - BUT CHECK THAT THESE ARE ACTUALLY YOUR VEOG CHANNELS - DIFFERENT PHD STUDENTS USE DIFFERENT CONVENTIONS!
art.blinkTh = 250; % threshold for detecting a blink in VEOG upper-lower difference signal (peak2peak moving)
art.blinkWindow = 200; % window (ms) for blink detection
art.blinkStep = 100; % stepsize (ms)
art.blinkTh2 = 100; % for eeglab func

art.saccTh = 50; % thresh for step-change (saccade detection)
art.saccWindow = 100; % ms
art.saccStep = 30; %ms

art.moveTh = [-200 200]; % [min max] change for detecting movement or other artefacts
art.artifTh = 100; % this is for SCALP channels - it's high because it is in CSD units, which take bigger values
art.chans = 1:128; % Artifact rejection channels - which electrodes should we reject based on?
art.windowTime = [-1000 0]; % in msec, the limits of the time window relative to the event in which you will check for an artifact
% will need to change this to run up until final cue/resp?
art.windowInds= isBetween(eeg.confRespTimes, art.windowTime); % indices of the timepoints for the artifact check window whose limits (min and max) are defined above

art.goodThresh = 640; % need at least 50% good trials?

[isBlink, isBlink2, nArtPerTr, isSaccade, isMovement, isBadRT,isBlink3, isMissing] = deal(NaN(fileInfo.maxTr, fileInfo.nPP));
isArt = NaN(eeg.nChans, fileInfo.maxTr, fileInfo.nPP);


for iPP = 1:fileInfo.nPP
    
    disp([fileInfo.ppID{iPP} '...'])
    
    load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} '_intp']),'erp','confCueS','RS')
    
%     
%     % average reference - already done
%     erp(1:eeg.nChans,:,:) = erp(1:eeg.nChans,:,:) - nanmean(erp(1:eeg.nChans,:,:));
    
%     erpArt = NaN(eeg.nChansTot, sum(art.windowInds), fileInfo.maxTr);
%     %% art window is -400: (RT+60samples or 1500ms if RT < 1000ms)
%     for iTr = 1:fileInfo.maxTr
%         if ~isnan(RS(iTr))
% %             if eeg.epochTimes(RS(iTr)) > 0 % if RT > 1sec
% %                 respLims = [art.windowTime(1) (RS(iTr)+60)/eeg.fs*1000]; % go until RT + 60 samples
% %             else 
% %                 respLims = art.windowTime; % else [-400 1500]
% %             end
% 
% %             % take from -400:confCueTime
% %             respLims = [art.windowTime(1) confCueS(iTr)/eeg.fs*1000];
% %             
% %     %         respSamples = round(respLims(1)/1000*512) : round(respLims(2)/1000*512); % convert to samples
% %             respInds = isBetween(eeg.epochTimes, respLims);
%             
%             respInds = art.windowInds; % -1000:0 to conf resp
%             erpArt(:,1:sum(respInds),iTr) = erp(:, respInds, iTr); % get resp epoch
%         end
%     end

    erpArt = erp(:,art.windowInds,:);
    %% artefact reject?

    isArt(:,:,iPP) = sq(max(abs(erpArt(1:eeg.nChans,:,:)),[],2) > art.artifTh);

    veog = erpArt(129,:,:) - erpArt(130,:,:);
    isBlink(:,iPP) = (max(veog,[],2)-min(veog,[],2)) > art.blinkTh;

    %% try other approaches
    isBlink3(:,iPP) = max(abs(veog),[],2) > art.blinkTh; % this will give a logical 1 if this threshold is exceeded on a certain trial, 0 otherwise
    isBlink2(:,iPP) = sq(myArtmwppth(veog, art.blinkTh2, art.blinkWindow, art.blinkStep, eeg.fs)); %[tr 1]
    isSaccade(:,iPP) = sq(myArtstep(veog, art.saccTh, art.saccWindow, art.saccStep, eeg.fs));
    isMovement(:,iPP) = sq(any(myArtextval(erpArt(1:eeg.nChansTot,:,:), art.moveTh))); %[chan 1 tr]. pop_artextval - min/max in window
                
    isMissing(:,iPP) = isnan(RS);
    isBadRT(:,iPP) = ~isBetween(RS, art.rtLimsSamp);
    
    figure(1);clf;
    bar(nansum(isArt(:,:,iPP),2)); xlabel('chan');ylabel('num arts');
    drawnow;
%     %% plot some stuff
%     
%     clf;
%     
%     isOut = isBlink; %isBlink | sq(any(isArt,1)); %sq(any(isArt,1));
%     t = {'accept','reject'};
%     for i = 1:2
%         subplot(2,1,i);
%         plot(eeg.epochTimes, sq(veog(1,:,isOut(:,iPP) == i-1)));
%         ylim([-500 500]);
%         ylabel('\muV'); xlabel('time');
%         xline(-200); xline(2500); yline(0);
%         title(t{i});
%     end
    
end

%%
isFlagged = isBlink==1 | sq(any(isArt==1,1)) ; % 1 = flagged for removal [tr pp]
isGood = isFlagged==0;

nGoodTrials = nansum(isGood);
ppsToExclude = nGoodTrials < art.goodThresh ; %4th dim

%% plot num good trials
figure();
hist(nGoodTrials);
xline(art.goodThresh);
xlabel('# good trials');
ylabel('# pps');

%% plot veog

figure();
h = errorBarPlot(permute(groupMeans(veog,3,f.isGood(:,iPP),'dim'),[2,1,3]),'area',1);
legend([h{:,1}], {'reject','include'},'Location','Best');
ylabel('VEOG');
xlabel('sample from baseline start');

%% are they missing at random?

load(fullfile(outFolder, 'BehDataLoad.mat'));

% regression with means rejected by condition (+ acc?)

factors = {'corrLR','cond','respLR1','acc','respLR2','certainty','CoM','conf3'};
for i = 1:length(factors)
    means = groupMeans(isGood', 2, behData.(factors{i})); %[pp cond (interrupt/cont)]
    anovas.(factors{i}) = rmanova(means, {'pp','cond'},'categorical',2);
    disp(factors{i});
    disp(anovas.(factors{i}));
end



%% save

save(fullfile(outFolder, 'FlagArtefactsConfResp.mat'), ...
    'isFlagged','isBlink','isArt','isMovement','isBlink2',...
    'ppsToExclude','isBadRT', 'art','isSaccade','isBlink3','isMissing','nGoodTrials','isGood');