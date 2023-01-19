% BehDataLoad
% Load up each eeg data file, extract behaviour, save as one file
% save
% behData.certainty is the main confidence metric used: 1=low(maybe),
% 2=medium (probably), 3=high (certain)
% behData.cond is post-decison evidence condition: 1=interrupted,
% 2=continued
% 

clear; close all;

%% 

outFolder = './Saves';

% loop through each file again, and load the raw eeg file
load(fullfile(outFolder, 'ExtractEpochs.mat'));

for iPP = 1:fileInfo.nPP
    
    eegFile = fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw.mat']);
    
    if exist(eegFile, 'file')
        fprintf('\n loading %s\n', fileInfo.ppID{iPP});
        d = load(eegFile, 'isLeft','cond','respLR','RS','corrLR','acc','confCueS','confRS','confResp','dataMat','dataMatNames');
        
        behData.isLeft(iPP,:) = d.isLeft; %which stimulus is correct 
        behData.cond(iPP,:) = d.cond; % post-dec evidence
        behData.respLR1(iPP,:) = d.respLR; % initial response button
        behData.RS(iPP,:) = d.RS; % initial response RT (in samples, from triggers)
        behData.corrLR(iPP,:) = d.corrLR; % initial response correct answer
        behData.acc(iPP,:) = d.acc; % initial response accuracy
        behData.confCueS(iPP,:) = d.confCueS; % confidence cue onset time (samples)
        behData.confRS(iPP,:) = d.confRS; % confidence response RT (samples)
        behData.confResp(iPP,:) = d.confResp; % confidence button response
        behData.pp(iPP,:) = repmat(iPP,1,fileInfo.maxTr); % ppID

        dataMat(iPP,:,:) = permute(d.dataMat,[3,1,2]); % store all other data in matrix
    else
        fprintf('\n Not Found: %s\n', fileInfo.ppID{iPP});
        
    end
    
end
    
dataMatNames = d.dataMatNames; 

%% convert

behData.confResp(behData.confResp==0) = NaN;
% make conf into left/right
behData.respLR2 = 1 + double(behData.confResp>3); % [l r] = [1 2]
behData.respLR2(isnan(behData.confResp)) = NaN;

% confidence response accuracy
behData.confAcc = double(behData.respLR2 == behData.corrLR);
behData.confAcc(isnan(behData.confResp)) = NaN;

% certainty/uncertainty [maybe probably certain]
behData.certainty = abs(round(behData.confResp -3.5)); %1-6 -> [3 2 1 1 2 3]

% did they change their mind?
behData.CoM = double(behData.respLR1 == (3-behData.respLR2)); % is opposite resp
behData.CoM(isnan(behData.confResp)) = NaN;

% get conf in respLR1: 6 = high, 1 =low (i.e. certain wrong)
behData.confInR1 = behData.confResp;
behData.confInR1(behData.respLR1 == 1) = 7 - behData.confInR1(behData.respLR1 == 1); % swap left ones from 1:6 to 6:1

% convert from [1 2 3 4 5 6] to [1 1 2 2 3 3]. not really using this
behData.conf3 =  ceil(behData.confInR1 / 2); 


% get confidence that it was left
behData.scoretemp = behData.confResp;
behData.scoretemp(behData.isLeft==1) = 7 - behData.scoretemp(behData.isLeft==1);

% get confidence in correct option
behData.confInCorr = behData.confResp;
behData.confInCorr(behData.corrLR == 1) = 7 - behData.confInCorr(behData.corrLR==1);

%% also get these

% get RT in time from samples
behData.RT = behData.RS/512 * 1000;
behData.confRT = behData.confRS / 512 * 1000;
behData.score = dataMat(:,:,14);

% invert cond
behData.cond2 = 3 - behData.cond; %[1 = cont, 2 = inter]

%% reset NaN in acc

behData.acc = double(behData.acc);
behData.acc(isnan(behData.RT)) = NaN;

%% labels for things

labels.cond = {'Extinguished','Continued'};
labels.cond2 = flip(labels.cond);
labels.corrLR = {'left', 'right'};
labels.respLR = {'left', 'right'};
labels.acc = {'error', 'correct'};
labels.CoM = {'noChange','change'};
labels.confInR1 = {'certain CoM', 'probably CoM', 'maybe CoM', 'maybe no-CoM', 'probably no-CoM', 'certain no-CoM'};
labels.conf3 = {'certain/probably CoM', 'maybe CoM/no-CoM', 'probably/certain no-CoM'};
labels.certainty = {'maybe CoM/no-CoM', 'probably CoM/no-CoM', 'certain CoM/no-CoM'};
labels.confInCorr = labels.confInR1;

%% make table

behVarNames = {'acc','RT','confAcc','confRT','confInR1','certainty','CoM','confInCorr','conf3'};
nBehVars = length(behVarNames);
behDataNames = fieldnames(behData); % names of all behData

behNames = ['pp','cond',behVarNames]; % not DVs
behTab = struct2table(structfun(@col, rmfield(behData, behDataNames(~ismember(behDataNames, behNames))), 'UniformOutput',0));
behNames =  behTab.Properties.VariableNames; % update

% make list of names with 'logistic' appended
logIVs = {'acc','confAcc','CoM'};
ivNames = behVarNames;
logInds = find(ismember(behVarNames, logIVs));
ivNames(logInds) = strcat(ivNames(logInds), 'Logistic');

% make copy before nanzscoring
for j = logInds
    behTab.(ivNames{j}) = behTab.(behNames{strcmp(behNames, behVarNames{j})});
end

% also get log RT 
behTab.RTLog = log(behTab.RT);
behTab.confRTLog = log(behTab.confRT);

behNames = behTab.Properties.VariableNames; % update

isLogistic = ismember(behNames, ivNames(logInds)); % which are [0 1] for logistic regressions

% don't z-score yet, will have to exclude later
% behTab(:,~isLogistic) = varfun(@nanzscore, behTab(:,~isLogistic));
glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

% may need to rezscore after excluding people/trials


%% save

clear d; % don't save this
save(fullfile(fileInfo.rawFolder, 'BehDataLoad.mat'));
save(fullfile(fileInfo.outFolder, 'BehDataLoad.mat'));
        