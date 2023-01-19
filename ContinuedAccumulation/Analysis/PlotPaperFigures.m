function PlotPaperFigures()
% statistics are in the files that these scripts load up, saved as
% 'paperFits' or something similar
% The data files loaded from the following functions all have the name of
% the script used to create them in their name e.g. 'CPPAnalysis_CSD.mat'
% was created by 'CPPAnalysis.m' and so on
% 
% use print('-vector','-dsvg','FILENAME') to save topoplots as svg
% 
% make sure Crameri, matlib, eeglab, are on the path
addpath ../../Funcs/



%% main Figures

%% First-order CPP and behavioural data for task 1 & 2
PlotNewFig2; 
% requires:
% RDMManualConf/Analysis/Saves/CPPAnalysis_CSD.mat
% RDMManualConf/Analysis/Saves/BehDataLoad.mat
% RDMManualConf/Analysis/Saves/FlagArtefacts.mat
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat
% ContinuedAccumulation/Analysis/Saves/BehDataLoad.mat
% ContinuedAccumulation/Analysis/Saves/FlagArtefacts3.mat
% 


%% Second-order (post-decision) CPP for Task 2
PlotNewFig3; % also plots supplementary figure of initial accuracy effect
% Requires:
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat


%% Plot CPP peak latency + final-RTs
PlotNewFig6; % this is figure 5 now
% Requires:
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat

%% Post-decision phase-tagging (differential SSVEP) and beta (motor prep)
PlotNewFig5; % this is figure 6 + 7 now
% Requires:
% ContinuedAccumulation/Analysis/Saves/PhaseAnalysis_CSD_ExclRT.mat
% ContinuedAccumulation/Analysis/Saves/MuBetaAnalysis_CSD__.mat
% ContinuedAccumulation/Analysis/Saves/GetPhaseTaggingTopo_CSD.mat
% ContinuedAccumulation/Analysis/Saves/ExtractEpochs.mat
% ContinuedAccumulation/Analysis/Saves/BehDataLoad.mat
% ContinuedAccumulation/Analysis/Saves/FlagArtefacts3.mat
% ContinuedAccumulation/Analysis/Saves/DoFFTCSD.mat

%% Plot optional stopping illustration
OptStopExample2; % Figure 8


%% below are all the supplementary figures
%%%%%%%%%%

%% post-choice CPP between Continued & Extinguished - suppl fig 1
% and Initial Accuracy split also - suppl fig 5

PlotNewFig3(1); 
% Requires:
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat


%% Task 2 first-order CPP with CoM included 
PlotNewFig2(1); % suppl fig 2
% Requires:
% 
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_ExclCoMFromCert.mat

%% pre-resp vs pre-evidence baseline-corrected post-choice CPP
% suppl figure 3
PlotNewFig4; % plots post-choice CPP for Extinguished evidence trials, comparing two baseline-corrections
% Requires:
% ContinuedAccumulation/Analysis/Saves/CPPAnalysis_CSD_.mat

%% 30Hz SSVEP - Supple fig 4
PlotSSVEP30Fig;
% Requires 
% ContinuedAccumulation/Analysis/Saves/SSVEPAnalysis_CSD_.mat
% ContinuedAccumulation/Analysis/Saves/GetSSVEPTopo_CSD.mat
% ContinuedAccumulation/Analysis/Saves/ExtractEpochs.mat
% ContinuedAccumulation/Analysis/Saves/BehDataLoad.mat
% ContinuedAccumulation/Analysis/Saves/FlagArtefacts3.mat



end