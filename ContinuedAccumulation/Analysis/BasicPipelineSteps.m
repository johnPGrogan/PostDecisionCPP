%% BasicPipeline


eeglab; close all;


%% find bdf files for pps, load, extract, filter, epoch, baseline
ExtractEpochs
    

%% bad channel checking
CheckForBadChannels
    

%% interpolate bad channels
InterpBadChannels('./Saves/ExtractEpochs.mat','D:\TCD\Projects\ContinuedAccumulation\Data\Interp')
    

%% artefact flagging 
FlagArtefacts.m
    
    
%% CSD transform

CSDTransform('D:\TCD\Projects\ContinuedAccumulation\Data\Interp', 'D:\TCD\Projects\ContinuedAccumulation\Data\CSD', 0);

%% response-lock

RespLock.m

%% get behavioural data

BehDataLoad.m

%% Get and analyse

GetCPP.m
GetCPPTopo.m
CPPAnalysis.m

%% run FFT

doFFT.m

%% mu-beta

GetMuBeta.m
GetMuBetaTopo.m
MuBetaAnalysis.m

%% do stim noise analysis (loadNoise must be run before SSVEP stuff)

LoadNoise;
NoiseAnalysis;

%% 30 Hz SSVEP

% GetSSVEPTopo.m % if you want to plot topographies of SSVEP
GetSSVEP.m
ssvepAnalysis.m

%% phase

GetPhaseTaggingTopo.m
GetPhaseTagging.m
PhaseAnalysis.m

SSVEPPhaseComparison.m % compare 30Hz and 15Hz dSSVEP

%% Fronto-central electrodes

GetFC.m
FCAnalysis.m


%% plot figures

PlotPaperFigures.m


