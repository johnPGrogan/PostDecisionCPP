% RDMBasicPipelineSteps
% These are the scripts I ran for the analysis for the paper.
% The plotting is done with
% ../../ContinuedAccumulation/Analysis/PlotPaperFigs.m

%% Extraction

ExtractEpochsRDM % load up data, filter, epoch, save

%% find bad channels to mark for interpolation

CheckForBadChannelsRDM

%% interpolate

InterpBadChannels('./Saves/ExtrachEpochs.mat','D:\TCD\Projects\RDMManualConf\Data\Interp');

%% Flag artefacts

FlagArtefactsRDM

%% CSD

CSDTransform('D:\TCD\Projects\RDMManualConf\Data\Interp', 'D:\TCD\Projects\RDMManualConf\Data\CSD', 0);

%% resplock

RespLockRDM

% also lock to response cue onset
RespCueLockRDM

%% behavioural data

BehDataLoadRDM

%% CPP

GetCPPRDM
GetCPPTopoRDM
GetCPPRespCueTopoRDM
CPPAnalysisRDM

%% FFT

DoFFTRDM

%% beta

GetMuBetaRDM
MuBetaAnalysisRDM




