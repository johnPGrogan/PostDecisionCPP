% GetMuBeta
% load FFTs
% remove flagged, grand topo pre-resp, pick elecs, outlier removal
% go back and apply the same to the stimlocked too

clc; clear all; close all;

load('./Saves/ExtractEpochs.mat')
load('./Saves/FlagArtefacts3.mat','isFlagged','ppsToExclude');
load('./Saves/BehDataLoad.mat','behData'); % for l/r

loadUp = 1; % load up PreRespMeanVoltage file?
useCSD = 1; % CSD or voltage

if useCSD
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT_CSD';
    load('./Saves/DoFFTCSD.mat'); % get FFT info
    loadName = 'PreRespMeanBetaCSD.mat';
    saveName = 'GetMuBetaCSD.mat';
    mapLims = [-1.3 1.3];

else
    fileInfo.fftFolder = 'D:\TCD\Projects\ContinuedAccumulation\Data\FFT';
    load('./Saves/DoFFT.mat'); % get FFT info
    loadName = 'PreRespMeanBeta.mat';
    saveName = 'GetMuBeta.mat';
    mapLims = [-.1 .1];

end

preRespWindow = [-150 -50]; % window before init RT
preRespInds = isBetween(respWindows, preRespWindow);


%% load up all - get the grand average for good trials

if loadUp && exist(fullfile(fileInfo.fftFolder, loadName), 'file')
    r = load(fullfile(fileInfo.fftFolder, loadName),'betaLat');
    betaLat = r.betaLat;

else

    betaLat = NaN(eeg.nChansTot, fileInfo.nPP);
    betaTopo = NaN(eeg.nChansTot, fileInfo.nPP);
    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])

        % get resplocked data
        data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'STFT');
     
        % remove flagged trials
        data.STFT(:,:,isFlagged(:,iPP)) = NaN;
        
        % lateralisation = left resp - right resp
        betas = nanmean(data.STFT(:,preRespInds,:),2); % average over time
        betaLat(:,iPP) = nanmean(betas(:,1,behData.respLR1(iPP,:)==1),3) - nanmean(betas(:,1,behData.respLR1(iPP,:)==2),3);
        
        % also store this for lateralisation data for topoplot at RT
        betaTopo(:,iPP) = nanmean(data.STFT(:,respWindows==0,behData.respLR1(iPP,:)==1),3) - nanmean(data.STFT(:,respWindows==0,behData.respLR1(iPP,:)==2),3);

    end

    save(fullfile(fileInfo.fftFolder, loadName), 'betaLat', 'preRespWindow','betaTopo');
end


%% topoplot that




figure();
topoplot(nanmean(betaLat(:,:),2),...
    eeg.chanlocs, 'electrodes','on','colormap',crameri('vik'),'mapLimits',mapLims,...
    'emarker',{'.','k',10,1});
colorbar


% % pick two clusters (left + right hemisphere)
chanClusters = {'D19','D28','D18'; 'B22','B21','B18'}'; % 
[~, chansInds] = ismember(col(chanClusters), eeg.chanNames);

keyboard; %%%%%% edit these picked channels

%% topoplot per person

figure();
for iPP = 1:fileInfo.nPP
    if 1%~ppsToExclude(iPP)
        subplot(5,6,iPP);
        topoplot(betaLat(:,iPP),...
            eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),...
            'emarker',{'.','k',10,1}, 'emarker2', {chansInds, '.','k',10,1});
        title(iPP);
    end
end

%% pick maximal in that per person

betas = betaLat(chansInds,:); % get mean voltages per channel per pp   
betas = reshape(betas, 3,2,[]); % left/ right

% pick maximal amplitude surrounding response execution    

flips = [1 -1]; % invert right side
% do per side
for i = 1:2
    [m,j] = max(betas(:,i,:) * flips(i));
    betasPicked(:,i)  = m *flips(i);
    betaChans(:,i) = chanClusters(j,i); % get name
    betaChanInds(:,i) = chansInds(j + (i-1)*3); % index
end
    

%% load that data again and store cpp chans, remove outliers

betas = NaN(fileInfo.nPP, length(respWindows), fileInfo.maxTr, 2); % 4th dim is left/right
betaStims = NaN(fileInfo.nPP, length(stimWindows), fileInfo.maxTr, 2); % 4th dim is left/right

for iPP = 1:fileInfo.nPP
    disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'STFT');
    
    % store that channel as cpp
    betas(iPP, :, :,:) = permute(data.STFT(betaChanInds(iPP,:), :, :),[4,2,3,1]);
    
    betas(iPP,:,isFlagged(:,iPP),:) = NaN; % remove flagged trials before mean calc
    
    % remove outliers
    betaMeans= -diff(nanmean(betas(iPP, preRespInds,:,:),2),[],4); % mean over window, lateralise
        
    
    %%%%%%% also do the stim-locked
%     disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'STFT');
    
    % store that channel as cpp
    betaStims(iPP, :, :,:) = permute(data.STFT(betaChanInds(iPP,:), :, :),[4,2,3,1]);
        
    betaStims(iPP,:,isFlagged(:,iPP),:) = NaN; % remove flagged trials before mean calc
       
end

% betaStims(ppsToExclude,:,:,:) = NaN; % make sure they are gone

%% make both into [ipsi contra ipsi-contra]

% flip those where respLR1==2
for iPP = 1:fileInfo.nPP
    betas(iPP,:,behData.respLR1(iPP,:)==2,:) = flip(betas(iPP,:,behData.respLR1(iPP,:)==2,:),4);
    betaStims(iPP,:,behData.respLR1(iPP,:)==2,:) = flip(betaStims(iPP,:,behData.respLR1(iPP,:)==2,:),4);
end

% contra - ipsi
betas(:,:,:,3) = diff(betas,[],4);
betaStims(:,:,:,3) = diff(betaStims,[],4);

%% save

save(fullfile('./Saves',saveName), ...
    'betas','betaChans','betaChanInds','chanClusters','preRespWindow','preRespInds','betaStims','betasPicked');