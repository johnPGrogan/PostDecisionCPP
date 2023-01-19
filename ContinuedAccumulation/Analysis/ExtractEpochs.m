function ExtractEpochs()
% Adapted from SimonPipeline/ExtractERPs.m, by John Grogan, 2021.

if ~exist('pop_biosig','file')
    eeglab; close all;
end


%% set up folders

% change these folders to where you have saved the data
fileInfo.dataFolder = {'D:\TCD\Projects\WouterDataBackup\Continued accumulation\Full dataset'}; % 1-16
fileInfo.rawFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Raw'; % where saved to
fileInfo.interpFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\Interp';
fileInfo.outFolder = './Saves'; % this is where the new files will be saved in later scripts

%% find files

f = []; % combine across Day 1 + 2
for i = 1:length(fileInfo.dataFolder)
    f = cat(1, f, dir(fileInfo.dataFolder{i}));
end
fNames = {f.name}';

toKeep = cellRegexpi(fNames, 'PP\d\d_EXP\d\d.bdf') > 0; % bdfs for experimental blocks

fNames(~toKeep) = []; % remove
f(~toKeep) = [];

ppNames = cellfun(@(x) x(1:4), fNames, 'UniformOutput', 0); % get ppNames
blockNums = cellfun(@(x) str2double(x(9:10)), fNames, 'UniformOutput',0); % block numbers
blockNums = cat(1, blockNums{:}); % remove from cells

days = arrayfun(@(x) str2double(x.folder(end)), f); % get day folder 

% append day2 blockNums - if laoding from Day1 + 2 folders
blockNums(days==2) = blockNums(days==2) + 8;

fileInfo.files = f; % store
fileInfo.fNames = fNames;
fileInfo.ppNames = ppNames;
fileInfo.blockNums = blockNums;

fileInfo.ppID = unique(fileInfo.ppNames); % 1 per person

fileInfo.nPP = length(fileInfo.ppID);

fileInfo.nBlocks = 16;
fileInfo.nTr = 80;
fileInfo.maxTr = fileInfo.nBlocks * fileInfo.nTr;

%% set triggers

triggers.relevant = [12 4 5 8 9 13 16 17 20 21 24 25 28 29];

triggers.stimOnset = 12;
triggers.evOnset = [4 5 8 9]; % evidence on
triggers.resp1 = [16 17]; % l/r
triggers.confCueOnset = 13;
triggers.resp2 = [20 21 24 25 28 29]; % s,d,f,j,k,l

triggers.isLeft = [4 5]; % [4 5 8 9]
triggers.isCont = [4 8]; % [4 5 8 9]

%% EEG params

eeg.reref = 0; % use average reference
eeg.fs = 512; % sample rate
eeg.nChans = 128;  % number of EEG channels
eeg.nExts = 2;     % number of external channels (EOG etc)
eeg.nChansTot = eeg.nChans + eeg.nExts;

% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
chanInfo = load('chanlocsBioSemi128.mat'); % note this was made by calling>> readlocs('cap128.loc') , and, since the locations for the 8 externals (129:136) are all (1,0,0), getting rid of those by calling chanlocs = chanlocs(1:128)
eeg.chanlocs = chanInfo.chanlocs;
eeg.chanNames = {eeg.chanlocs.labels}';

% define the epoch to extract:
eeg.epochLimsMS = [-1000 3000];    % note the use of integer number of cycles of SSVEP (17Hz)
eeg.epochSamples = round(eeg.epochLimsMS(1)/1000*eeg.fs):round(eeg.epochLimsMS(2)/1000*eeg.fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
eeg.epochTimes = eeg.epochSamples*1000/eeg.fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds

%define the baseline window:
eeg.blWin = [-400 0]; % the unusual number -117.6 ms comes from taking exactly two cycles of the SSVEP flicker rate (17Hz). If there were no flicker, you might choose something rounder like [-100 0]

eeg.nSamples = length(eeg.epochSamples); % num timepoints in epoch

% uncomment this if locking to confidence response
% % conf resp locking
% eeg.confRespLimsMS = [-1500 500]; % around conf resp
% eeg.confRespSamples = round(eeg.confRespLimsMS(1)/1000*eeg.fs):round(eeg.confRespLimsMS(2)/1000*eeg.fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
% eeg.confRespTimes = eeg.confRespSamples*1000/eeg.fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds

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
    
    [corrLR, cond, respLR, RS, blockNum, confCueS, confRS, confResp] = deal(NaN(1, fileInfo.maxTr));
    erp = NaN(eeg.nChansTot, eeg.nSamples, fileInfo.maxTr);
    blAmp = erp(:,1,:);
    [epochTrigs, epochSTimes] = deal(num2cell(NaN(1, fileInfo.maxTr)));
    dataMat = NaN(fileInfo.maxTr, 15);
    
    
    % now list the eeg files to be used for each subject:
    thisPP = find(cellRegexpi(fileInfo.fNames, fileInfo.ppID(iPP)) == 1);
    thisBlocks = fileInfo.blockNums(thisPP);
    
    fileNames = {};   
    for iB = 1:fileInfo.nBlocks
        
        try % see if the file for this block exists?
            fileNames{iB} = fullfile(fileInfo.files(thisPP(thisBlocks(iB))).folder, fileInfo.files(thisPP(thisBlocks(iB))).name);
            fprintf('\nBlock %d: %s\n', fileInfo.blockNums(thisPP(thisBlocks(iB))), fileInfo.files(thisPP(thisBlocks(iB))).name);
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
        % remove empty chans
        EEG.data = EEG.data(1:eeg.nChansTot,:,:); % keep veog (+ 1 other external)?
        EEG.chanlocs = EEG.chanlocs(1:eeg.nChansTot);
        
        % Fish out the event triggers and times
        nEvents = length(EEG.event); % total number of events
        
        if nEvents == 0
            warning('no trigger events found in this block (duration is %d), skipping', length(EEG.times)/eeg.fs);
            continue; 
        end % no triggers in this block?
        
        trigs = [EEG.event.edftype];
        if length(trigs) < nEvents
            emptyInds = find(cellfun(@isempty, {EEG.event.edftype}));
            for i = 1:length(emptyInds)
                EEG.event(emptyInds(i)).edftype = str2num(EEG.event(emptyInds(i)).type);
            end

            trigs = [EEG.event.edftype];
            if length(trigs) < nEvents

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
            end
        end

        sTimes = round([EEG.event.latency]);

        % so trigs and stimes are the exact same length; trigs has the trigger codes and stimes the SAMPLE POINTS in the
        % continuous EEG where those triggers happened.

        % some people had wrong trigger codes - they had 768 added for some reason?
        % check and fix
        if ~any(ismember(trigs, triggers.relevant))
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
        
        if ~any(trigs(1) == [1 12])
            trigs(1)=1; % it's usually because the recording started late, so set this to 'start of 
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

        dataMatNames = {'pp','block','day','deltaC','cue','respLR','respT','respLRconf','respTconf','perf2','perf','conf','scoretemp','score','scalevalue'};
        mat.Confidence.scoretemp(mat.Confidence.scoretemp==0) = NaN;
        
        % if final trial was a miss there will be no confscores - append NaN
        if length(mat.Confidence.scoretemp)<80
            mat.Confidence.scoretemp(length(mat.Confidence.scoretemp):80) = NaN;
            mat.Confidence.Keypressed{length(mat.Confidence.scoretemp):80} = [];
        end
        origMat = [repmat(iPP,1,length(mat.respLR));....
            repmat(iB,1,length(mat.respLR));...
            repmat(ceil(iB/8),1,length(mat.respLR));...
            mat.par.dc; mat.par.cues';
            mat.respLR; mat.respT*1000;...% original choice
            mat.respLRconf; mat.respTconf*1000;... % conf choice
            mat.perf2; mat.perf;...
            mat.Confidence.confidence; mat.Confidence.scoretemp;...
            mat.Confidence.score; mat.Confidence.scalevalue]';      

        
        %% find epochs
                
        % find indices of triggers corresponding to the events we want to
        % time-lock to (evidence onset)
        trialStarts = find(trigs==triggers.stimOnset);
        
        % there are some missing or misaligned triggers (i.e. mismatched to
        % trials), I have manually fixed these known issues here, and am
        % matching trigger RT to behavioural RT below to check trials are
        % correctly aligned
        if all(length(trialStarts) ~= [length(mat.respT) sum(~isnan(mat.respT))])

            if iPP==9 && iB==11 % this person had no triggers recorded for the first trial
                trigs = [12 trigs]; %pretend stim onset
                sTimes = [1 sTimes]; % will exclude
            elseif iPP==13 && iB==3 % first 17 trials were missing triggers
                trigs = [repmat(12,1,17) trigs];
                sTimes = [1:17, sTimes];
            elseif iPP==16 && iB==2 % only final 12 trials have triggers captured
                trigs = [repmat(12,1,68) trigs];
                sTimes = [1:68, sTimes];
            elseif iPP==17 && iB==3 % first trial only has '8'
                trigs = [12 trigs];
                sTimes = [1 sTimes];
            elseif iPP==23 && iB==11 % this trialStart was received wrong for some reason
                trigs(trigs==44) = 12; 
            elseif iPP==25 && iB==13 % missing trials at start
                trigs = [repmat(12,1,5), trigs];
                sTimes = [1:5, sTimes];
            elseif iPP==29 && iB==7 
                % missing last 8 trials, so just continue
            else
                keyboard;
            end
            trialStarts = find(trigs==triggers.stimOnset);
        end
        
        % Now loop through the single trials and grow the trial-descriptor vectors and 'erp' matrix
        for iE = 1:length(trialStarts)
            
            iTr = (iB-1)*fileInfo.nTr + iE; % trial index in entire experiment

            %store beh data now, so will still be there even if trial is
            %invalid later
            dataMat(iTr,:) = origMat(iE,:); 
            
            % first, check whether the epoch relative to the current event is even contained within the bounds of the EEG data (maybe the recording was stopped mid-trial?) - if not, skip the trial ('continue')
            if sTimes(trialStarts(iE)) < abs(eeg.epochSamples(1)) 
                continue; 
            end
            if (sTimes(trialStarts(iE)) + eeg.epochSamples(end) > size(EEG.data,2))
                keyboard;
            end
            
            thisTr = trialStarts(iE); % start of this trial
            if ~ismember(trigs(thisTr+1), triggers.evOnset)
                % is probs a premature response - e.g. [12 17 12]
                if (thisTr+2)<=length(trigs) && ismember(trigs(thisTr+1), triggers.resp1) && ismember(trigs(thisTr+2), triggers.stimOnset)
                    continue;
                elseif (thisTr+1)<=length(trigs) && ~ismember(trigs(thisTr+1), triggers.relevant)
%                     keyboard; 
                    thisTarg = thisTr+2;
                    % if there's a random trigger between '12' and [4 5 8 9], can just skip it: thisTarg = thisTr+2;
                    % e.g. iPP=16,iB=16,iE=77, [12 32 5 17] 
                else % something else
                    % iPP==11 && iB==3 && iE=80, final trial in block was 
                    % [12 16] = [stimOnset, response], so can just continue
%                     fprintf('\nno target detected after stimulus: %d %d %d %d %d\n', trigs(thisTr+(-2:2)))
                    keyboard;
                    continue;
                end
            else
                thisTarg = thisTr+1; % index of this evOnset trigger
            end
            thisTrig = trigs(thisTarg); % trigger number
            
                                   
            nextRespInd = thisTarg + 1; % should be immediately after it           
            if nextRespInd<length(trigs) && ismember(trigs(nextRespInd), triggers.resp1) && ~isnan(mat.respT(iE))
                RS(iTr) = sTimes(nextRespInd)-sTimes(thisTarg);
                respLR(iTr) = trigs(nextRespInd)-15; % subtract 15 as a quick way to translate response trigger code to 1's and 2's for left/right
            else
                continue; % finish here
            end    
            
            % check RT roughly matches respT
            if abs( (RS(iTr)/eeg.fs*1000) - mat.respT(iE)*1000 ) > 10 % ms
                keyboard;
            end
            
            if iE<length(trialStarts) % check in this epoch
                nextTr = trialStarts(iE+1);
            else
                nextTr = length(trigs); % if last, check until end
            end

            % also get confResp + confRT
            % was there a confCue?
            confCueInd = find(ismember(trigs, triggers.confCueOnset) & sTimes>sTimes(thisTarg) & sTimes<=sTimes(nextTr),1);
            if ~isempty(confCueInd)
                confCueS(iTr) = sTimes(confCueInd) - sTimes(thisTarg);
                confRespInd = find(ismember(trigs, triggers.resp2) & sTimes>sTimes(thisTarg) & sTimes<=sTimes(nextTr),1);
                if ~isempty(confRespInd)
                    confRS(iTr) = sTimes(confRespInd) - sTimes(confCueInd);
                    confResp(iTr) = find(trigs(confRespInd)==triggers.resp2);
                    
                    if abs( (confRS(iTr)/eeg.fs*1000) - mat.respTconf(iE)*1000 ) > 10 % ms
                        keyboard;
                    end
                end
            end
            
                        
            % get trial info
            blockNum(iTr) = fileInfo.blockNums(thisPP(iB));
            corrLR(iTr) = 2 - any(thisTrig == triggers.isLeft); %1=left, 2=right
            cond(iTr) = 1 + any(thisTrig == triggers.isCont); %1=interrupt, 2=continued

            %% Now extract the epoch
            ep = EEG.data(:,sTimes(thisTarg) + eeg.epochSamples);
            %Baseline correction
            blAmp(:,1,iTr) = nanmean(ep(:,isBetween(eeg.epochTimes, eeg.blWin)),2); % store this
            ep = ep - blAmp(:,1,iTr);
            erp(:,:,iTr) = ep; % now add the current epoch onto the growing 'erp' matrix by concatenation along the 3rd dimension

            %%% uncomment this if doing conf-resp locking
%             %% only keep the bits around the confidence response times
%             % [-1500 500];
%             
%             % check there was a confresp
%             if ~isempty(confCueInd) && ~isempty(confRespInd)
% 
%                 if (sTimes(confRespInd) + eeg.confRespSamples(end)) <= size(EEG.data,2)
%                     ep = EEG.data(:, sTimes(confRespInd) + eeg.confRespSamples);
%                 else
%                     inds = eeg.confRespSamples <= (size(EEG.data,2) - sTimes(confRespInd));
%                     ep = EEG.data(:, sTimes(confRespInd) + eeg.confRespSamples(inds));
%                 end
%                 ep = ep - blAmp(:,1,iTr); % baseline using the ev-locked baseline
% 
%             else
%                 ep = NaN(eeg.nChansTot, length(eeg.confRespSamples));
%             end
% 
%             erp(:,1:size(ep,2),iTr) = ep; % now add the current epoch onto the growing 'erp' matrix by concatenation along the 3rd dimension
            
            %% also store al the triggers in an epoch
            
            %%% uncomment this if doing confResp locking
%             if ~isempty(confCueInd) && ~isempty(confRespInd) 
%                 % get all sTimes within the epoch
%                 trigInds = isBetween(sTimes, sTimes(confRespInd) + minMax(eeg.confRespSamples,2));

            trigInds = isBetween(sTimes, sTimes(thisTarg) + minMax(eeg.epochSamples,2));
            epochTrigs{iTr} = trigs(trigInds);
            epochSTimes{iTr} = sTimes(trigInds);
%             end
            

        end
    end
    
    
    
    %% It's worth taking a look at this subject's VEOG, just for the last block, to get a sense whether blinks will be well
%     % detected by a certain threshold:

    isLastBlock = blockNum==max(blockNum);
    clf;
    plot(reshape(erp(129,:,isLastBlock)-erp(130,:,isLastBlock), 1,[]));
    yline(250,'-k');
    title(['VEOG ' fileInfo.ppID{iPP}])
    
    %% Now that we have all trials from all blocks, save everything for this subject:
    
    epochs.erp = erp;
    epochs.baseline = blAmp;
    epochs.blockNum = blockNum;
    epochs.isLeft = corrLR;
    epochs.cond = cond;
    epochs.respLR = respLR;
    epochs.RS = RS; % maybe convert into ms from samples?
    epochs.corrLR = corrLR;
    epochs.acc = respLR == corrLR;
    epochs.confCueS = confCueS;
    epochs.confRS = confRS;
    epochs.confResp = confResp;
    epochs.trigs = epochTrigs;
    epochs.sTimes = epochSTimes;
    epochs.dataMat = dataMat;
    epochs.dataMatNames = dataMatNames;

%     %% append only dataMat
%     AppendIntoFiles(fileInfo.ppID{iPP}, epochs)
     
    %%
    fprintf('\n Saving raw file...\n');
    save(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']),...
        '-v7.3', '-struct','epochs');

end