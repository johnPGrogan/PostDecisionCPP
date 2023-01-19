% BehDataLoadRDM
% Load up each eeg data file, extract behaviour, save as one file
% behData.certainty is main confidence metric: 1=low/maybe,
% 2=medium/probably, 3=high/certain
% save

clear; close all;

%% 

outFolder = './Saves';

% loop through each file again, and load the raw eeg file
load(fullfile(outFolder, 'ExtractEpochs.mat'));

for iPP = 1:fileInfo.nPP
    
    eegFile = fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw.mat']);
    
    if exist(eegFile, 'file')
        fprintf('\n loading %s\n', fileInfo.ppID{iPP});
        d = load(eegFile, 'corrLR','stimDur','respLR','RS','acc','confResp','dataMat','dataMatNames');
        
        behData.stimDur(iPP,:) = d.stimDur; % duration of stimulus
        behData.respLR(iPP,:) = d.respLR; % direction of coherent motion
        behData.RS(iPP,:) = d.RS; % RT in samples
        behData.corrLR(iPP,:) = d.corrLR; % correct response
        behData.acc(iPP,:) = d.acc; % accuracy
        behData.confResp(iPP,:) = d.confResp; % confidence response
        behData.dataMat(iPP,:,:) = permute(d.dataMat,[3,1,2]); % store matrix
        behData.pp(iPP,:) = repmat(iPP,1,fileInfo.maxTr); % ppID

    else
        fprintf('\n Not Found: %s\n', fileInfo.ppID{iPP});
        
    end
    
end
    
%% convert

behData.confResp(behData.confResp==0) = NaN; % missing is NaN

% certainty/uncertainty [maybe probably certain]
behData.certainty = abs(round(behData.confResp -3.5)); %1-6 -> [3 2 1 1 2 3]

behData.dataMatNames = d.dataMatNames;

% get conf in correct response
behData.confInCorr = behData.confResp;
behData.confInCorr(behData.corrLR == 1) = 7 - behData.confInCorr(behData.corrLR==1);

%% also get these

% get RT in ms
behData.RT = behData.RS/eeg.fs * 1000;

%% reset NaN in acc

behData.acc = double(behData.acc);
behData.acc(isnan(behData.RT)) = NaN;

%% labels for things

labels.durNames = {'short','medium','long'};
labels.stimDur1 = {'350','500','750'};
labels.stimDur = {'350ms','500ms','750ms'};
labels.corrLR = {'left', 'right'};
labels.respLR = {'left', 'right'};
labels.acc = {'error', 'correct'};
labels.confInCorr = {'defWrong','probsWrong','maybeWrong','maybeRight','probsRight','defRight'};
labels.certainty = {'maybe','probably','certain'};
labels.confResp = {'defLeft','probsLeft','maybeLeft','maybeRight','probsRight','defRight'};

%% save

clear d; % don't save this
% save(fullfile(fileInfo.rawFolder, 'BehDataLoad.mat'));
save(fullfile(fileInfo.outFolder, 'BehDataLoad.mat'));
        