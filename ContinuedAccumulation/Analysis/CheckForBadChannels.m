% CheckForBadChannels
% This script opens matfiles with single-trial ERPs that have been
% extracted - usually without any re-referencing or interpolation ('raw')
% but possibly with low-pass filter and detrending done - and examines the
% data for channels that are abnormally high in variance
clear; close all;

outFolder = './Saves';
load(fullfile(outFolder, 'ExtractEpochs.mat'));


loadOld = 0;
if loadOld
    % load this - will have ones already done
    c = load('ChansToInterp2.mat', 'chansToInterp');
    chansToInterp = c.chansToInterp; 
    
    % don't do the ones already done?
    ppsToDo = find(~ismember(fileInfo.ppID, fieldnames(c.chansToInterp)))';
%     ppsToDo = 1:fileInfo.nPP;
else
    chansToInterp = struct();
    ppsToDo = 1:fileInfo.nPP;
end



%% 

% We're going to plot the standard deviation of each channel and spot ones that stick out. It's worthwhile to particularly
% scrutinize electrodes that might be particularly important in the final analysis. Here, we want to use electrodes around
% electrode CPz and around left and right motor areas (see 'cap_128_layout_medium.jpg'):
importantChans = [2 3 4 19 84 85 86 52 53 54 110 114 115 104 103 94 93 81 80 72 71 59]; % so we can pay particular attention to those
reRefChan = 23; % pick a channel to re-reference the data to and get a second picture of variance across channels - 
% important just because sometimes the reference used to load the EEG might itself have been bad during the recording...
% This is Oz

% Cz=A1, Pz=A19, Oz=A23, FPz=C17, Fz=C21, C3=D19, T7=D23, C4=B22, T8=B26
% [1, 19, 23, 54, 58, 81, 85, 115, 119]
%%

% Manually go through each subject in turn and use the code snippets that are commented out near the end to investigate bad channels
for iPP = ppsToDo % change the number manually, before the %, and hit F5
    close all;
    data = load(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw'])); % contains epochs structure's fields
    
    data.erp = data.erp(1:eeg.nChans,:,:); % only EEG channels

    %%
    blockNums = unique(data.blockNum(~isnan(data.blockNum)));
    
    % concatenate all epochs and get Standard deviation (SD) per channel
%     [SD, SD2] = deal(zeros(eeg.nChansTot, length(unique(blockNum))));
%     for iBl = blockNums
%         iTrs = find(blockNum==iBl); % find the indices of the trials from the current block, which are not outliers
%         conc = reshape(erp(:,:,iTrs),[eeg.nChansTot,length(iTrs)*eeg.nSamples]); % concatenate all trials
%         conc2 = conc - repmat(conc(reRefChan,:),[eeg.nChansTot,1]); % Also look at SD when referenced to somewhere at the back - electrode Oz for example
%         for iCh = 1:eeg.nChansTot % include externals here as well, because it might be a good idea to check whether those are noisy too
%             SD(iCh,iBl) = std(conc(iCh,:)); % measure S.D. of each channel 
%             SD2(iCh,iBl) = std(conc2(iCh,:));
%         end
%     end

    %% remove very bad trials
    
    % % %     % can remove entire bad trials and recalc SD?
    badTrials = find(nanmean(nanstd(data.erp(1:eeg.nChans,:,:),[],2) > 100, 1)>.9);
    if any(badTrials)
        disp(badTrials);
        data.erp(:,:,badTrials) = NaN;
    end

    
    conc = reshape(data.erp, eeg.nChans, eeg.nSamples*fileInfo.nTr, fileInfo.nBlocks); % [chan time*tr block]
    
    SD = sq(nanstd(conc, [], 2));
    SDr = sq(nanstd(conc - conc(reRefChan,:,:),[],2));
%     SDr = sq(nanstd(conc - nanmean(conc,1),[],2));
    
    
    %% also get artefact counts
    
    isArt = sq(max(abs(data.erp(1:eeg.nChans, isBetween(eeg.epochTimes, [-200 2500]),:)),[],2) - 100 > 0);  % is uV over 100 in window?

    chanArts = sum(isArt,2); % arts per channel
    chanArtsPerBlock = sq(sum(reshape(isArt,eeg.nChans,80,16),2));
    % will interpolate if > 40? 30?
    
    %% are there any channels that stick out in terms of variance?
    
    
    figure();
    plot(0:100, prctile(col(SD), 0:100), '-x');
    xlabel('%'); ylabel('SD');
    ylim([0 200]);
    % find the inflection point, can use that to set sdLims
    sdLims = [1 100]  
    
%     keyboard; % adjust sdLims
    %%
    figure; 
    subplot(2,1,1); hold on; plot(SD); 
    ylabel('SD'); xlabel('channel'); 
    title(['pp ' num2str(iPP) ' ' fileInfo.ppID{iPP}])
    plot([1; 1]*importantChans, [0;1] * max(SD(importantChans,:),[],2)', '-k')
    
    legend(num2str(blockNums'),'Location','Best');
    yline(sdLims(2),':k');
    subplot(2,1,2); hold on; plot(SDr); ylabel('SD'); xlabel('channel');
    
    figure();
    imagesc(SD, sdLims .* [1 1]); colorbar
    ylabel('channels'); xlabel('block');
    % horizontal stripes are bad channels, vertical are bad blocks
    
    SDct = squeeze(nanstd(data.erp,[],2));
    figure; imagesc(SDct,sdLims.*[1 1]); ylabel('channel'); xlabel('trial'); colorbar
    hold on; plot(80:80:1280,eeg.nChans, 'kx')
%     horizontal stripes are bad channels, vertical are bad trials/times
%     
%     figure(); 
%     imagesc(SDct(importantChans,:), sdLims.*[1 1.2]); 
%     hold on; plot(80:80:1280,max(ylim), 'kx')
%     yticks(1:length(importantChans)); yticklabels(importantChans); 
%     xticks(1:fileInfo.nTr:fileInfo.maxTr); xticklabels(1:fileInfo.nBlocks);
%     ylabel('important channels'); xlabel('block'); 
%     now taking a look through the marked 'important' channels, I see 52 sort of sticks out in block (check [sortedSD,I] = sort(SD(52,:)))...5! So finally: 
%     
%     plot artefacts too
    figure();
    imagesc(chanArtsPerBlock, [0 80]); colorbar
    title('artefacts');
    ylabel('channel'); xlabel('block');
    
%     [bl,ch] = find(chanArtsPerBlock' > 40); sortrows([bl,ch])


    % other useful plotting code:
    % plot a selected bunch of electrodes, concatenated across trials:
%     elec2plot = [2 14 17 23 36 55 56]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))

    %% examine the plots
    
    ch2interp = cell(1,length(blockNums)); % makes an empty cell array for filling in the bad channels for this subject
   
    % look at the figures - the first figure (subplot) can be used to see
    % which blocks stick out - the second one (imagesc) also shows this
    
    % the next one (imagesc) shows SD on each trial, look for horizontal or
    % vertical stripes to decide if it is a bad channel or trial
    % you can click on a bright patch and use (get, gca, 'CurrentPoint') to
    % find the [x y z] coords - there will be some slight rounding errors
    % depending on how large the image is onscreen
    
    % the next one shows the same but just for the 'importantChans'
    
    % use 
    % [bl,ch] = find(SD(1:128,:)'>sdLims(2)); sortrows([bl,ch])
    % to find the blocks that are bad for a certain channel
    
    % or
    % [b, i] = sort(SD(13,:))
    
    % add them to the 'bad' variable below = {block [badChansInThatBlock]}
    % This will fill in chansToInterp, which is saved for the
    % InterpBadChannels.m file
    
    


%     [bl, ch] = find(~isBetween(SD', sdLims) | chanArtsPerBlock' > 30);
%     sortrows([bl, ch])
%     
%     nBadChans = length(unique(ch))% > 10%?
%     [a,b] = CountUnique(bl); [b,a]

    ch2interp = [] % enter numbers here - will be done for all blocks
    
    keyboard; % edit the thresholds above, or add others into [bl ch] to be included below
%% add chans to interpolate per block

%     for iBl = 1:fileInfo.nBlocks
%         nBadChansPerBlock = length(unique(ch(bl == iBl)));
%         if nBadChansPerBlock > eeg.nChans/10
%             keyboard;
%         end
%         ch2interp{iBl} = ch(bl == iBl);
%     end
    
    %% store
    chansToInterp.(fileInfo.ppID{iPP}) = ch2interp; % store each pp in a field
    chansToInterp = orderfields(chansToInterp);
    
    save('ChansToInterp2.mat', 'chansToInterp', 'importantChans');

    %%%% also store in the loaded file
    save(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']),...
        '-append','ch2interp','sdLims');
end


