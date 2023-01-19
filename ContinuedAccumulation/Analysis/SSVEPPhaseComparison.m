% SSVEPPhaseComparison
% Loads up 15Hz dSSVEP and 30Hz SSVEP, checks settings match up, runs an
% LME comparing them

%% compare 30Hz and 15Hz signals with freq*acc*com LME
% after doing <10 exclusions etc

clc; clear; close all;

%% 
showErrBars = 1;
faceAlpha = .2;

colours.certainty = [0    0.4470    0.8; 0.8500    0.3250    0.0980; .0    0.7510         0; .2 .2 .2];
colours.cond = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; .2 .2 .2];
colours.CoM = [0.4660    0.6740    0.1880; 0.9290    0.6940    0.1250; .2 .2 .2]; % CoM
colours.confInR1 =  [crameri('roma',6); .2 .2 .2];
colours.topo = crameri('vik');


opts.useCSD = 1;
opts.excludeByRT = 1; % remove trials outside [100 1500] ms
opts.excludeBadPps = 1;
opts.excludeTooFew = 1;
excludeCoMFromCert = 0; % remove CoM trials from behData.certainty

% opts.outFolder = './Saves';
outFolder = '../../../../../OneDrive - TCDUD.onmicrosoft.com/ProjectsOD/ContinuedAccumulation/Analysis/Saves/';

%% load phase

opts.saveOpts = {'Volt','CSD'; '', 'ExclRT'};
opts.saveName = sprintf('PhaseAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.excludeByRT+1});

optsNames = fieldnames(opts);
data15 = load(fullfile(outFolder, opts.saveName), optsNames{:}, ...
    'meanPhase','phaseNames','nPhase', 'regTab', ...
    'behData','factors','nF', 'labels','toRemove');


% check things match

dataNames = fieldnames(data15);
data15b = rmfield(data15, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data15b)
    warning('loaded data and options do not match');
    keyboard;
end

% move into workspace
% struct2workspace(data15)

%% now load 30hz


opts.useCPPSSVEP = 0; % load up SSVEP in CPP chans

opts.saveOpts = {'Volt','CSD'; '', 'CPPChans'};
opts.saveName = sprintf('SSVEPAnalysis_%s_%s.mat', opts.saveOpts{1,opts.useCSD+1}, opts.saveOpts{2, opts.useCPPSSVEP+1});

optsNames = fieldnames(opts);
data30 = load(fullfile(outFolder, opts.saveName), optsNames{:}, ...
    'ssvepMeans','lockNames','behData','toRemove','regTab');


% check things match

dataNames = fieldnames(data30);
data30b = rmfield(data30, dataNames(~ismember(dataNames, optsNames)));

if ~equals(opts, data30)
    warning('loaded data and options do not match');
    keyboard;
end


%% 15hz has interrupted set to NaN, so do this for 30Hz too

data30.regTab(data30.regTab.cond<=0,:) = table(NaN);

% check tables match up

if ~all( (data15.regTab.accLogistic == data30.regTab.accLogistic) & ...
        (data15.regTab.CoMLogistic == data30.regTab.CoMLogistic) & ...
        (data15.regTab.confAccLogistic == data30.regTab.confAccLogistic) | isnan(data30.regTab.accLogistic))
    warning('regTabs do not match');
end

%% combine and analyse

regTab = data15.regTab;
regTab.resp30 = data30.regTab.resp;

% need to duplicate it 
regTab2 = vertcat(regTab, regTab);
regTab2.ssvep = nanzscore([regTab.resp; regTab.resp30]); % combine
regTab2.freq = [ones(height(regTab),1); ones(height(regTab),1)*2]; % freq
regTab2.freq(isnan(regTab2.ssvep)) = NaN;
regTab2.freq = nanzscore(regTab2.freq);


%% fit

formula = 'ssvep ~ 1 + freq*acc*confInR1 + (1 | pp)';

fit = fitglme(regTab2, formula)
% 3-way interaction is significant, meaning acc*CoM is diff in each signal

%% sep 

formulaSep = repmat({'ssvep ~ 1 + acc*confInR1 + (1 | pp)'},2,1);

fitglmeCellSep = @(x, i) fitglme(regTab2(regTab2.freq>0 == i,:), x); 

fitSep = cellfun(fitglmeCellSep, formulaSep, {0;1},'UniformOutput',0);
% these give the same results to the separate LME in each signal done in
% the separate scripts I loaded the data from

