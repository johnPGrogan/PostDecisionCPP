function ExtractEpochsRDM()
% Adapted from SimonPipeline/ExtractERPs.m, by John Grogan, 2021.
% For Wouter's task 2 - RDM manual + confidence
clc;
if ~exist('pop_biosig','file')
    eeglab; close all;
end


%% set up folders

% adjust these paths to where you saved your data
fileInfo.dataFolder = {'D:\TCD\Projects\WouterDataBackup\EEG data\'}; % 1-16
fileInfo.rawFolder = 'D:\TCD\Projects\RDMManualConf\Data\Raw';
fileInfo.interpFolder = 'D:\TCD\Projects\RDMManualConf\Data\Interp';
fileInfo.outFolder = './Saves'; % this is where output files go - make sure this folder exists

%% find files

f = []; % combine across Day 1 + 2
for i = 1:length(fileInfo.dataFolder)
    f = cat(1, f, dir(fileInfo.dataFolder{i}));
end
fNames = {f.name}';

toKeep = cellRegexpi(fNames, '\w\w') > 0; % folders with data

fNames(~toKeep) = []; % remove
f(~toKeep) = [];

fileInfo.ppNames = fNames;
fileInfo.ppID = fileInfo.ppNames;
fileInfo.ppNum = 1:length(fNames); % give them a number too

fileInfo.nPP = length(fileInfo.ppNames);

fileInfo.nBlocks = 8;
fileInfo.nTr = 72;
fileInfo.maxTr = fileInfo.nBlocks * fileInfo.nTr;


%% go through each folder and get the bdf files

f2 = cell(fileInfo.nPP, fileInfo.nBlocks);
for iPP = 1:fileInfo.nPP
    f1 = dir(fullfile(f(iPP).folder, f(iPP).name)); 
    toKeep = cellRegexpi({f1.name}, '.bdf')>0;
    f2(iPP,:) = {f1(toKeep).name};
end

blockNums = cellfun(@(x) str2double(x(6:7)), f2, 'UniformOutput',0); % block numbers
fileInfo.blockNums = cell2mat(blockNums); % remove from cells


fileInfo.folders = f; % store
fileInfo.files = f2;


%% set triggers

triggers.relevant = [5 16 21 22 23 31 32 33 36 101 102 103 104 105 106];

triggers.stimOnset = 5;
triggers.evOnset = 16; % evidence on
triggers.evOffset = [101 102 103 104 105 106]; % at offset
triggers.resp1 = [23 22 21 31 32 33]; % l/r

triggers.stimLR = [1 1 1 2 2 2]; % to match triggers.evOffset
triggers.stimDur = [350 500 750 350 500 750]; % for triggers.evOffset
triggers.respLR = [1 1 1 2 2 2]; % to match triggers.resp1
triggers.confResp = [1 2 3 4 5 6]; % to match resp1


%% EEG params

eeg.reref = 0; % use average reference
eeg.fs = 512; % sample rate
eeg.nChans = 128;  % number of EEG channels
eeg.nExts = 5;     % number of external channels (EOG etc)
eeg.nChansTot = eeg.nChans + eeg.nExts;

% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
chanInfo = load('chanlocsBioSemi128.mat'); % note this was made by calling>> readlocs('cap128.loc') , and, since the locations for the 8 externals (129:136) are all (1,0,0), getting rid of those by calling chanlocs = chanlocs(1:128)
eeg.chanlocs = chanInfo.chanlocs;
eeg.chanNames = {eeg.chanlocs.labels}';

% define the epoch to extract:
eeg.epochLimsMS = [-1050 3750];    % resp can then be [-1500:1500], and trim stim to [-500 1500]
eeg.epochSamples = round(eeg.epochLimsMS(1)/1000*eeg.fs):round(eeg.epochLimsMS(2)/1000*eeg.fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
eeg.epochTimes = eeg.epochSamples*1000/eeg.fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds

%define the baseline window:
eeg.blWin = [-200 0]; % to motion onset

eeg.nSamples = length(eeg.epochSamples); % num timepoints in epoch

%% init filters

filters.doDetrend = 1;     % detrend the data (1) or not (0)?
filters.doLPF = 1;    % 1 = low-pass filter the data, 0=don't.
filters.doHPF = 0;    % high-pass filter the data? Usually we don't do this unless the drift in the data is bad and not linear (so detrending alone won't get rid of it)

filters.loCutOff = 40;       % Low Pass Filter cutoff in Hz
filters.hiCutOff = 0;    % Cutoff frequency for the HIGHPASS filter (which is a 'low cut off' of the overall banpass we're effectively doing). Typically 0.05 Hz and very rarely above 0.1 Hz because it can result in distortions of the ERPs - see Steve Luck stuff online

%% save this info

save(fullfile(fileInfo.outFolder, 'ExtractEpochs.mat'),'fileInfo','filters','eeg','triggers');

%% extract and save each pp - looping through all available blocks

for iPP = 1:fileInfo.nPP
    extraction(iPP, fileInfo, eeg, filters, triggers);
end

end

function epochs = extraction(iPP, fileInfo, eeg, filters, triggers)
% epochs = extraction(iPP, fileInfo, eeg, filters, triggers)
% Loops through all blocks for 1 person, loads EEG, trims to relevant
% triggers, applies filters, extracts epochs, baselines, and gets some
% behavioural response data from epochs too. Saves outputs


%% extract

% Now loop through the subjects and extract the single trial ERPs and the other important informaiton about every trial:

    % preallocate
    [corrLR, stimDur, stimDurSample, respLR, RS, blockNum, confResp] = deal(NaN(1, fileInfo.maxTr));
    erp = NaN(eeg.nChansTot, eeg.nSamples, fileInfo.maxTr);
    blAmp = erp(:,1,:);
    [epochTrigs, epochSTimes] = deal(num2cell(NaN(1, fileInfo.maxTr)));
    dataMat = NaN(fileInfo.maxTr, 12);
    nChansOrig = NaN(1, fileInfo.nBlocks);
      
    fileNames = {};   
    for iB = 1:fileInfo.nBlocks
        
        try % see if the file for this block exists?
            fileNames{iB} = fullfile(fileInfo.folders(iPP).folder, fileInfo.folders(iPP).name, fileInfo.files{iPP, iB});
            fprintf('\nBlock %d: %s\n', iB, fileNames{iB});
            EEG = pop_biosig(fileNames{iB}); % read in EEG - this is an EEGLAB function that outputs a structure 'EEG' with lots of fileds, the most important of which is EEG.data - the actual EEG data!
        catch ME
            fprintf('\nfile not found: %s\n', fileNames{iB});
            disp(ME);
            continue; % skip to next block
        end

        if EEG.srate ~= eeg.fs 
            EEG = pop_resample(EEG,eeg.fs);
        end
        

        %%
        % some people have diff n chans
        nChansOrig(1,iB) = size(EEG.data,1);
        
        % remove empty chans
        EEG.data = EEG.data(1:eeg.nChansTot,:,:); % keep veog (+ 1 other external)?
        EEG.chanlocs = EEG.chanlocs(1:eeg.nChansTot);
        
        % Fish out the event triggers and times
        nEvents = length(EEG.event); % total number of events
        
        if nEvents == 0
            warning('no trigger events found in this block (duration is %d), skipping', length(EEG.times)/eeg.fs);
            continue; 
        end % no triggers in this block?
        
        trigs = {EEG.event.type}; % this is sometimes 'condition X' or just 'X' (where X is a string number)
        
        if iscellstr(trigs)
            % get triggers as numbers only
            trigsNum = NaN(1,nEvents);
            for i = 1:nEvents
                trigsNum(i) = str2num(trigs{i}(regexp(trigs{i}, '\d')));
            end
            trigs = trigsNum; % replace
        elseif all(isnumeric([trigs{:}]))
            trigs = [trigs{:}]; % 
        else
            keyboard;
        end

        sTimes = round([EEG.event.latency]);
        % so trigs and stimes are the exact same length; trigs has the trigger codes and stimes the SAMPLE POINTS in the
        % continuous EEG where those triggers happened.

        % one person had wrong trigger codes - they had 768 added for some
        % reason?
        % check and fix
        if ~any(ismember(trigs, triggers.relevant))
            keyboard;
            fprintf('\n no relevant trigger codes found, will try adjusting by 768 to fix known error...')
            if all(ismember(trigs-768, [triggers.relevant 0]))
                trigs = trigs - 768; % minus
                trigs(trigs==0) = 1; % the 'start' should be 1 not zero
                if all(ismember(triggers.relevant, trigs))
                    fprintf(' fixed!\n');
                else
                    fprintf(' NOT FIXED - manually adjust\n');
                    beep;
                    keyboard;
                end
            else
                fprintf(' NOT FIXED - manually adjust\n');
                beep;
                keyboard;
            end
        end
        
        if ~any(trigs(1) == [4 5])
            trigs(1)=1; % usually because recording was late starting
%             keyboard;
        end
        
        % detrend?
        if filters.doDetrend
            EEG.data = detrend(EEG.data')';
        end
        
        if filters.doHPF
            [B,A] = butter(3,filters.hiCutOff/(eeg.fs/2),'high'); % butter inputs are order, cutoff freq (in normalised scale where 1=half sample rate, and optional third argument telling it to be e.g. 'high'-pass)
            EEG.data = filtfilt(B,A,double(EEG.data)')';
        end
        
        % LP Filter
        if filters.doLPF
            EEG.data = eegfilt(EEG.data,eeg.fs,filters.hiCutOff,filters.loCutOff);
        end 
        
    
        %% load beh data
        matName = [fileNames{iB}(1:end-4) '.mat'];
        mat = load(matName);

        dataMatNames = {'pp','block','coherence','stimDur','direction','respLR','respT','perf','acc','acc2','confScale','scalevalue'};
%         mat.Confidence.scoretemp(mat.Confidence.scoretemp==0) = NaN;
        
        nTrThisBlock = length(mat.trial.condition);

        % check some stuff 
        if nTrThisBlock < fileInfo.nTr
            if iPP~=9 && iPP~=11
                keyboard; % iPP=9, first block is only 48 trials but last is 92 to make it up, I will correct this at the end.
                %iPP=11, first block is only 48 trials. others are 72
            end
        end
        if ~all(mat.trial.condition == mat.trial.trialCond)
            keyboard; 
        end


        % if final trial was a miss there will be no confscores - append NaN
        if length(mat.trial.Confidence.scalevalue) < nTrThisBlock
            mat.trial.Confidence.scalevalue(length(mat.trial.Confidence.scalevalue):nTrThisBlock) = NaN;
            mat.trial.Confidence.scaledconf(length(mat.trial.Confidence.scaledconf):nTrThisBlock) = NaN;
            mat.trial.Confidence.Keypressed{length(mat.trial.Confidence.Keypressed):nTrThisBlock} = [];
        end

        directions = [180 180 180 0 0 0]; % set these, so can get values per trial
        durations = [350 500 750 350 500 750];

        mat.respLR = (mat.trial.Confidence.scalevalue > 3) + 1; %1=left(s,d,f), 2=right(j,k,l)
        
        % check (pp1 has a leading zero here giving 73 trials in
        % mat.triggers.respLR, so ignore that)
        if length(mat.triggers.RespLR)<length(mat.respLR) || ~all(mat.respLR == mat.triggers.RespLR(end-nTrThisBlock+1:end))
%             keyboard; % sometimes the mat.triggers.RespLR has a zero, can ignore as will use mat.respLR
        end

        mat.confScale = mat.trial.Confidence.scaledconf .* (-1 + 2*mat.trial.Acc); % [-3 -2 -1 1 2 3] (neg are incorrect)

        % make into matrix
        origMat = [repmat(iPP,1,nTrThisBlock);
            repmat(iB,1,nTrThisBlock);
            mat.trial.coherence;
            durations(mat.trial.condition);
            directions(mat.trial.condition);
            mat.respLR; 
            mat.RTs*1000; 
            mat.triggers.perf;
            mat.trial.Acc;
            mat.Acc; % not sure what this is, ranges from 0-1. may be cumul acc or similar?
            mat.confScale; 
            mat.trial.Confidence.scalevalue;]';      

        
        %% find epochs
                
        % find indices of triggers corresponding to the events we want to
        % time-lock to (evidence onset)
%         targs = find(any(trigs == triggers.evOnset'));
        trialStarts = find(trigs==triggers.stimOnset);
        
        % fix known issues due to missing trials in EEG files etc. Will
        % match trials by RT later to validate them
        if all(length(trialStarts) ~= [nTrThisBlock  sum(~isnan(mat.RTs))])
            % missing triggers, usually at start, so add in trialStarts
            if iPP==3 && iB==1 % seems there's a missing trial in this block, plus a missing response trigger in trial 16 or 21?
                trigs = [5 trigs];
                sTimes = [1 sTimes];
            
            elseif iPP==3 && iB==1  % last trial is missing, so can ignore
                % do nothing

            elseif iPP==6 && iB==8 % last trial missing so can ignore
                % do nothing

            elseif iPP==18 && iB==4 % very first trigger missed, but rest of trial is there
                trigs = [5 trigs];
                sTimes = [sTimes(1)-512 sTimes]; % this is the mean number of samples between 5&16, so guess that.
                % shouldn't matter much as only differs by 1 frame (~14ms)
            
            elseif iPP==22 && iB==2 % 61 trials here, 72 in mat. 11 trials missing from beginning
                trigs = [repmat(5,1,11), trigs];
                sTimes = [1:11, sTimes];

            else
                keyboard;
            end
            trialStarts = find(trigs==triggers.stimOnset); % redo
        end
        
        % Now loop through the single trials and grow the trial-descriptor vectors and 'erp' matrix
        for iE = 1:length(trialStarts)
            
            iTr = (iB-1)*fileInfo.nTr + iE; % trial index in entire experiment

            % store this data now, so conditions are available even if no
            % response is found
            dataMat(iTr,:) = origMat(iE,:);

            % first, check whether the epoch relative to the current event is even contained within the bounds of the EEG data (maybe the recording was stopped mid-trial?) - if not, skip the trial ('continue')
            if sTimes(trialStarts(iE)) < abs(eeg.epochSamples(1)) || (sTimes(trialStarts(iE)+1) + eeg.epochSamples(end) + 512) > size(EEG.data,2)
                if iE==length(trialStarts) && (sTimes(trialStarts(iE)) + eeg.epochSamples(end)) <= size(EEG.data,2)
%                     keyboard; % this could potentially be included, only 1 trial found so far
                end
                continue; 
            end
            
            thisTr = trialStarts(iE); % find start of trial
            if ~ismember(trigs(thisTr+1), triggers.evOnset)
                keyboard;
            else
                thisTarg = thisTr+1; % index of this evOnset trigger
            end
            thisTrig = trigs(thisTarg); % trigger number
            
            % after the evOnset comes the evOffset. resp should only be
            % after this
            thisOffset = thisTarg + 1; % ev offset 
            if thisOffset <= length(trigs) && ismember(trigs(thisOffset),  triggers.evOffset)
                stimDurSample(iTr) = sTimes(thisOffset) - sTimes(thisTarg); % actual duration (samples)
                corrLR(iTr) = triggers.stimLR(trigs(thisOffset) == triggers.evOffset); % 1/2 l/r
                stimDur(iTr) = triggers.stimDur(trigs(thisOffset) == triggers.evOffset); % expected duration (ms)
            else
                keyboard;
            end

            nextRespInd = thisOffset + 1; % should be immediately after it. some RTs will be after 1500s so will be skipped
            if nextRespInd<=length(trigs) && ismember(trigs(nextRespInd), triggers.resp1) && ~isnan(mat.RTs(iE)) 
                RS(iTr) = sTimes(nextRespInd)-sTimes(thisOffset);
                respLR(iTr) = triggers.respLR(trigs(nextRespInd) == triggers.resp1); % 1/2 l/r
                confResp(iTr) = triggers.confResp(trigs(nextRespInd) == triggers.resp1); % 1-6
            else
                continue; % finish here
            end    
            
            % check RT roughly matches respT
            if abs( (RS(iTr)/eeg.fs*1000) - mat.RTs(iE)*1000 ) > 20 % ms
                keyboard;
            end
            
            if iE<length(trialStarts) % check in this epoch
                nextTr = trialStarts(iE+1);
            else
                nextTr = length(trigs); % if last, check until end
            end

            % are there any extra triggers left in this trial (may be
            % outside epoch below)
            if nextTr > (nextRespInd+1)
%                 keyboard;
                extraTrigs{iTr} = [trigs(nextRespInd:nextTr);  sTimes(nextRespInd:nextTr)];
            end
            
            % get trial info
            blockNum(iTr) = fileInfo.blockNums(iPP,iB);
            

            %% Now extract the epoch
            ep = EEG.data(:,sTimes(thisTarg) + eeg.epochSamples);
            %Baseline correction
            blAmp(:,1,iTr) = nanmean(ep(:,isBetween(eeg.epochTimes, eeg.blWin)),2); % store this
            ep = ep - blAmp(:,1,iTr);
            erp(:,:,iTr) = ep; % now add the current epoch onto the growing 'erp' matrix by concatenation along the 3rd dimension

            
            %% also store al the triggers in an epoch
            
            % get all sTimes within the epoch
            trigInds = isBetween(sTimes, sTimes(thisTarg) + minMax(eeg.epochSamples,2));
            
            epochTrigs{iTr} = trigs(trigInds);
            epochSTimes{iTr} = sTimes(trigInds);
            
        end
    end
    
    %% check eeg vs datamat?
    
    
    if iPP==9 % remove extra blank trials in iPP=9 from dataMat, so that checks will line up
        dataMat(all(isnan(dataMat),2),:) = [];
    end
        

    eegStuff = [stimDur' corrLR' respLR' confResp'];
    matStuff = [dataMat(:,4), (2-dataMat(:,5)./180), dataMat(:,6), dataMat(:,12)];
    rtStuff = RS'./eeg.fs*1000 - dataMat(:,7);
    if ~all(eegStuff ~= matStuff | isnan(eegStuff) | isnan(matStuff) | abs(rtStuff) <= 20,'all')
        keyboard;
    end
    
    
    %% It's worth taking a look at this subject's VEOG, just for the last block, to get a sense whether blinks will be well
%     % detected by a certain threshold:

%     isLastBlock = blockNum==max(blockNum);
%     clf;
%     plot(reshape(erp(131,:,isLastBlock)-erp(133,:,isLastBlock), 1,[]));
%     yline(250,'-k');
%     title(['VEOG ' fileInfo.ppID{iPP}])
%     
    %% Now that we have all trials from all blocks, save everything for this subject:
    
    epochs.erp = erp;
    epochs.baseline = blAmp;
    epochs.blockNum = blockNum;
    epochs.corrLR = corrLR;
    epochs.stimDur = stimDur;
    epochs.stimDurSample = stimDurSample;
    epochs.respLR = respLR;
    epochs.RS = RS; % maybe convert into ms from samples?
    epochs.acc = respLR == corrLR;
    epochs.confResp = confResp;
    epochs.trigs = epochTrigs;
    epochs.sTimes = epochSTimes;
    epochs.dataMat = dataMat;
    epochs.dataMatNames = dataMatNames;
    epochs.nChansOrig = nChansOrig;
    
    %% iPP==9. first block was only 48 trials, but last was 92 to make up for it
%     so I'm going to remove those 24 NaN trials from block 1, to make the
%     number of trials match up. blockNum is in dataMat(:,2)

    if iPP == 9
        toRemove = isBetween(1:size(epochs.erp,3), [49 72]); % only these 24 were NaN
       
        epochs2 = rmfield(epochs, {'baseline','erp','dataMat','dataMatNames','nChansOrig'});
        % trial is 2nd dim
        
        epochs2 = structfun(@(x) x(:,~toRemove), epochs2,'UniformOutput',0);
        
        % these have trials in other dims
        epochs2.baseline = epochs.baseline(:,:,~toRemove);
        epochs2.erp = epochs.erp(:,:,~toRemove);
%         epochs2.dataMat = epochs.dataMat(~toRemove,:); % this is already done above now
        
        epochs2.dataMatNames = epochs.dataMatNames; % doesnt have trials
        epochs2.nChansOrig = epochs.nChansOrig;
        
        epochs = epochs2; % replace
    end
    
    
    %%
    fprintf('\n Saving raw file...\n');
    %     %% append only dataMat
%     AppendIntoFiles(fileInfo.ppID{iPP}, epochs, 'D:\TCD\Projects\RDMManualConf\Data\')

    save(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']),...
        '-v7.3', '-struct','epochs');

end