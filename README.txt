This repository contains the analysis code for the paper CITATION.
The task code, and anonymous pre-processed and analysed data are available at https://osf.io/4dqkz/

Written by John P Grogan, some code adapted from Simon Kelly's code.

Tasks
	RDMManualConf/ is Task 1 (random dot motion with simultaneous manual response of direction and confidence).
	ContinuedAccumulation/ is Task 2 (two-stage contrast discrimination with initial identity response, and delayed confidence response, with evidence interrupted or continued in the delay period).

Running
	Ensure the Funcs/ folder is on your path, and download the data from OSF. If running the RDMManualConf analysis, make sure the ContinuedAccumulation/Analysis folder is also on your path.
	You can run ContinuedAccumulation/Analysis/PlotPaperFigs.m to make the figures from the paper, and the files those scripts load up contain the statistics.
	There are also all the scripts needed pre-process and analyse the data for those figures, which can be run from ContinuedAccumulation/Analysis/BasicPipelineSteps.m and RDMManualConf/Analysis/BasicPipelineStepsRDM.m.

Requirements
	This analysis ran in Matlab R2021b, with with following packages:
		EEGLab v2021.1
		ERPLab v8.20
		BIOsig 3.7.9
		CDSToolbox v1.1
		Crameri (for the colormaps)
		Matlib (https://github.com/sgmanohar/matlib)
		barwitherr.m

	and Matlab toolboxes:
		Signal processing toolbox
		
	(My functions required (in Funcs/ folder):
		cellRegexpi.m
		col.m
		CountUnique.m
		emptyLegend.m
		FitCPPSlope.m
		GetCellsWithLowN.m
		isBetween.m
		minMax.m
		myArtmwppth.m
		pbar.m
		tail.m)


Data files required (see Data/ component in OSF project https://osf.io/4dqkz/)

	(these are the files actually loaded for the figures)
		BehDataLoad.mat
		CPPAnalysys_CSD_.mat
		CPPAnalysis_CSD_ExcluCoMFromCert.mat
		CPPPeakLatency.m
		DoFFTCSD.mat
		ExtractEpochs.mat
		FlagArtefacts3.mat
		GetMuBetaTopoCSD.mat		
		GetPhaseTaggingTopo_CSD.mat
		GetSSVEPTopoCSD.mat
		MuBetaAnalysis_CSD_.mat
		MuBetaAnalysis_CSD_ExclCoMFromCert.mat
		PhaseAnalysis_CSD_ExclRT.mat
		PhaseTopoAnalysis_CSD.mat
		SSVEPAnalysis_CSD_.mat
		
		(and these in RDMManualConf/Analysis/Saves/)
		BehDataLoad.mat
		CPPAnalysis_CSD.mat
		ExtractEpochs.mat
		FlagArtefacts.mat
		
	(these files are made using the pre-processed data, )
		ChansToInterp2.mat
		FlagArtefacts3.mat
		ExtractEpochs.mat
		BehDataLoad.mat
		DoFFTCSD.mat
		GetCPPCSD.mat
		GetCPPRespCueTopoCSD
		GetCPPTopoCSD.mat
		GetMuBetaCSD.mat
		GetMuBetaTopoCSD.mat
		GetPhaseTagging_CSD2.mat
		GetSSVEPCSD.mat
		PhaseTopoAnalysis_CSD.mat
		PreRespMeanBetaCSD.mat
		
	(plus these in RDMManualConf/Analysis/Saves)
		BehDataLoad.mat
		ExtractEpochs.mat
		DoFFTCSD.mat
		FlagArtefacts.mat
		GetCPPCSD.mat
		GetCPPTopoCSD.mat
		GetCPPRespCueTopoCSD.mat
		GetMuBetaCSD.mat
		PreRespMeanBetaCSD.mat
		
	(in Funcs/)
		chanlocsBioSemi128.mat
		CSD_coords.mat
		
