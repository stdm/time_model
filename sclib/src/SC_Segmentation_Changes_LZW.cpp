/**************************************************************************/
/*    Responsibility:																											*/
/*      - This class implements the (speaker) change detector introduced  */
/*        in 'Real-Time Unsupervised Speaker Change Detection', L.Lu,     */
/*        H.J.Zhang, 2002 (IEEE)                                          */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#include "SC_Segmentation_Changes_LZW.h"
#include "SC_ModelHandler.h"
#include "SC_MixtureModel_GMM.h"
#include "SC_GroundTruth_SCiVo.h"
#include "SC_Model_qGMM.h"
#include "SC_FeatureHandler.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_Changes_LZW::SC_Segmentation_Changes_LZW(SC_TweakableParameters* pTweak, int mode) : SC_Segmentation_Changes(pTweak, mode) {
  this->pMatrixFunc = new SC_MatrixFunctions();
  this->pDist = new SC_DistanceMeasures(this->pTweak, this->pMatrixFunc);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_Changes_LZW::~SC_Segmentation_Changes_LZW() {
  MFree_0D(this->pDist);
  MFree_0D(this->pMatrixFunc);
}

//====================================================================================================================
// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
// implementation
// The features needed for training are meant to reside in an array of feature-containers as returned by 
// SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
//====================================================================================================================
int SC_Segmentation_Changes_LZW::trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	int idxMFCC = sclib::bitPosition(sclib::featureMFCC), idxLSP = sclib::bitPosition(sclib::featureLSP), idxPitch = sclib::bitPosition(sclib::featurePitch);
	
	if (pFeatures[idxMFCC] == NULL || pFeatures[idxLSP] == NULL || pFeatures[idxPitch] == NULL ||
		  pFeatures[idxMFCC]->Hdr.frameSize != pFeatures[idxLSP]->Hdr.frameSize ||
			pFeatures[idxMFCC]->Hdr.frameStep != pFeatures[idxLSP]->Hdr.frameStep ||
			pFeatures[idxMFCC]->Hdr.sampleRate != pFeatures[idxLSP]->Hdr.sampleRate ||
			pFeatures[idxPitch]->Hdr.frameSize != pFeatures[idxLSP]->Hdr.frameSize ||
			pFeatures[idxPitch]->Hdr.frameStep != pFeatures[idxLSP]->Hdr.frameStep ||
			pFeatures[idxPitch]->Hdr.sampleRate != pFeatures[idxLSP]->Hdr.sampleRate) {
    return SVLIB_Fail; //pitch, mfccs and lsps must be extracted on the same frame-basis, otherwise return error
  } else {
		if (this->mode == sclib::modeSpeakerChange) {
			return trainLuZhangWindow(pGT, segmentStart, segmentEnd, pFeatures, sclib::atSpeech, sclib::atPause|sclib::atSilence, sclib::atSpeakerBoundary);
		} else {
			//return trainLuZhangWindow(pGT, segmentStart, segmentEnd, pFeatures, sclib::atNoise, sclib::noType, sclib::atNoiseBoundary);
			return SVLIB_Fail; //we typicallay have no groundtruth for acoustic changes and if we had, no methods to evaluate them just like speaker-boundary data
		}
  }
}

//====================================================================================================================
//	estimates the needed prior probabilities for the bayesian fusion engine from the labeled training data and prints
//  it to a file
//====================================================================================================================
int SC_Segmentation_Changes_LZW::trainLuZhangWindow(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, unsigned long int type, unsigned long int typesNot, unsigned long int boundaryType) {
	unsigned int f, featureIdx[3];
  long int windowCount, changeCount = 0;
  unsigned long int extra, y, windowLength, windowStep, maxLen, minLen, maxDim, minDim, *time2Change;
	double **oldCov[3], **actCov[3];
  double **distances; //array to hold: distance to predecessor: [x][0] MFCC-based, [x][1] LSP-based, [x][2] Pitch-based
	double **changeDist, **noChangeDist;
  unsigned long int **segments; //array to hold: [x][0] segStart; [x][1] segEnd (in correspondence with the same indices of distances[])
  SV_Data *pStrippedFeatures[3]; //copied-together features, [0] MFCCs, [1] LSPs, [2] Pitch
	double *changeMean, *noChangeMean, *changeVariance, *noChangeVariance, time2ChangeMean, time2ChangeVar;
	SV_Data *pPriors = new SV_Data(6, 3); // P_change          | time2ChangeMean  | time2ChangeVar
		                                    // P_noChange        | 0.0              | 0.0
	                                      // changeMean_MFCC   | changeMean_LSP   | changeMean_Pitch
	                                      // changeVar_MFCC    | changeVar_LSP    | changeVar_Pitch
	                                      // noChangeMean_MFCC | noChangeMean_LSP | noChangeMean_Pitch
	                                      // noChangeVar_MFCC  | noChangeVar_LSP  | noChangeVar_Pitch
	bool includesChangePoint;
	SC_FeatureHandler *pSaver = new SC_FeatureHandler(this->pTweak);
	SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak, true);
	int bestOrder;
	SV_Data *pTime2ChangeData = NULL;
	SC_Model *pTime2ChangeModel;
	
	//some initializations
	for (f = 0; f < 3; f++) {
		oldCov[f] = NULL;
		actCov[f] = NULL;
	}
	featureIdx[0] = sclib::bitPosition(sclib::featureMFCC);
	featureIdx[1] = sclib::bitPosition(sclib::featureLSP);
	featureIdx[2] = sclib::bitPosition(sclib::featurePitch);

  //collect all desired frames in this segment, with their frameNr in the zeroth col.
	//because all 3 feature-sets must have equal frame sizes and steps, the ones from mfccs are used here and later on without loss of generality
  segmentEnd = sclib::min(segmentEnd, pGT->getConverter()->audioFrame2sample(pGT->getLastAudioFrameNrInSegment(segmentStart, segmentEnd, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep), pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd));
  pStrippedFeatures[0] = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureMFCC)], segmentStart, segmentStart, segmentEnd, type, typesNot, 1, 1, true);
	pStrippedFeatures[1] = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureLSP)], segmentStart, segmentStart, segmentEnd, type, typesNot, 1, 1, true);
	pStrippedFeatures[2] = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featurePitch)], segmentStart, segmentStart, segmentEnd, type, typesNot, 1, 1, true);
  if (pStrippedFeatures[0] == NULL) { //no frames of desired type here
		MFree_0D(pPriors);
		MFree_0D(pSaver);
		MFree_0D(pModelHandler);
    return SVLIB_Fail;
  }

  //determine window parameters
 	windowLength = pGT->getConverter()->ms2audioFrame(this->pTweak->segmentationChangesLz.detectorWindowLength, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd); //measured in audioFrames
	windowStep = pGT->getConverter()->ms2audioFrame(this->pTweak->segmentationChangesLz.detectorWindowStep, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd); //measured in audioFrames
	windowCount = pStrippedFeatures[0]->Row / windowStep;
  windowCount -= (unsigned long)ceil((double)(((windowCount-1)*windowStep + windowLength) - pStrippedFeatures[0]->Row) / (double)windowStep); //subtract all frames which overlap the amount of frames
  extra = pStrippedFeatures[0]->Row - (sclib::max(0, windowCount-1)*windowStep + windowLength) - 1; //the last window can be extra long, if there are some frames which are too less to make an extra window
  if (windowCount <= 0) {
		MFree_0D(pPriors);
		MFree_0D(pSaver);
		MFree_0D(pModelHandler);
    MFree_0D(pStrippedFeatures[0]);
		MFree_0D(pStrippedFeatures[1]);
		MFree_0D(pStrippedFeatures[2]);
    return SVLIB_Fail;
  }

  MArray_2D(distances, windowCount, 3, double, "SC_Segmentation_Changes_LZW.trainLuZhangWindow: distances"); //make it so big that it can hold the max. possible count of segments
  MArray_2D(segments, windowCount, 2, unsigned long int, "SC_Segmentation_Changes_LZW.trainLuZhangWindow: segments"); //make it so big that it can hold the max. possible count of segments
  MArray_2D(changeDist, windowCount, 3, double, "SC_Segmentation_Changes_LZW.trainLuZhangWindow: changeDist"); //make it so big that it can hold the max. possible count of segments
	MArray_2D(noChangeDist, windowCount, 3, double, "SC_Segmentation_Changes_LZW.trainLuZhangWindow: noChangeDist"); //make it so big that it can hold the max. possible count of segments
	MArray_1D(time2Change, windowCount, unsigned long int, "SC_Segmentation_Changes_LZW.trainLuZhangWindow: time2Change"); //make it so big that it can hold the max. possible count of segments

  //loop over all subsegments (x second windows with y second overlap), check if they're valid (long enough etc.)
  //compute divergence-shape distance for each two windows of lsps here
  for (y = 0; y < (unsigned long)windowCount; y++) {
		minLen = y * windowStep; //select a proper submatrix that contains the current window; omit the first column (contains only the frame-number)
		maxLen = minLen + windowLength + ((y == windowCount-1) ? extra : 0);
		minDim = 1;

		for (f = 0; f < 3; f++) {
			if (actCov[f] != NULL) { //save last valid cov
				MFree_2D(oldCov[f]); 
				oldCov[f] = actCov[f];
			}
			maxDim = pStrippedFeatures[f]->Col; 
			actCov[f] = this->pMatrixFunc->cov(pStrippedFeatures[f]->Mat, pStrippedFeatures[f]->Row, pStrippedFeatures[f]->Col, NULL, minLen, maxLen, minDim, maxDim); //compute new cov

			//fill distances
			segments[y][0] = (unsigned long)pStrippedFeatures[f]->Mat[minLen][0]; //original start-sample of first frame in the window
			segments[y][1] = (unsigned long)pStrippedFeatures[f]->Mat[maxLen][0] + pGT->getConverter()->audioFrame2sample(0, pStrippedFeatures[f]->Hdr.frameSize, pStrippedFeatures[f]->Hdr.frameStep, sclib::alignmentEnd); //original end-sample of last frame in the window
			if (y == 0) { //this is the first (valid) window in this scene, it has no predecessor to compute distance
				distances[y][f]	= -1.0;
				time2Change[changeCount] = 0;
			} else { //this window has a predecessor, so compute distance
				distances[y][f]	= this->pDist->divergenceShape(actCov[f], oldCov[f], maxDim-minDim);
				if (distances[y][f] <= 0.0) {
					//REPORT_ERROR(SVLIB_Fail, "Matrix Inversion failed... maybe too short window?");
					//distances[y][f] = 0;
				}
			}
		} //for f

		//create raw data to estimate probabilities
		includesChangePoint = (pGT->isSpeakerHomogenious(y > 0 ? segments[y-1][0] : segments[y][0], segments[y][1], sclib::modeGroundtruth) == false) ? true : false;
		if (includesChangePoint == true) {
			for (f = 0; f < 3; f++) {
				changeDist[changeCount][f] = distances[y][f];
			}
			changeCount++;
			time2Change[changeCount] = 0;
		} else {
			time2Change[changeCount]++; //create a statistic with the number of successive segments without change-point
			for (f = 0; f < 3; f++) {
				noChangeDist[y-changeCount][f] = distances[y][f];
			}
		}

	} //for all windows

	//free some memory
	for (f = 0; f < 3; f++) {
    MFree_0D(pStrippedFeatures[f]);
		MFree_2D(oldCov[f]);
		MFree_2D(actCov[f]);
	}
	MFree_2D(distances);
  MFree_2D(segments);
	MFree_2D(distances);
  MFree_2D(segments);
	
	//approximate distance-distribution per feature under each hypothesis via a single gaussian
	changeMean = this->pMatrixFunc->mean(changeDist, changeCount, 3);
	changeVariance = this->pMatrixFunc->variance(changeDist, changeCount, 3, changeMean);
	noChangeMean = this->pMatrixFunc->mean(noChangeDist, windowCount-changeCount, 3);
	noChangeVariance = this->pMatrixFunc->variance(noChangeDist, windowCount-changeCount, 3, noChangeMean);
	
	//first approximation of the density of the number of windows before a change occurs
	time2ChangeMean = this->pMatrixFunc->mean(time2Change, changeCount+1);
	time2ChangeVar = this->pMatrixFunc->variance(time2Change, changeCount+1, &time2ChangeMean);

	//fill the prior-matrix and save it
	pPriors->Mat[0][0] = (float)((double)(changeCount+1) / (double)(windowCount+1)); //laplace estimator: adding a small constant (e.g. 1) to both parts of the fraction to avoid zero probability
	pPriors->Mat[0][1] = (float)(time2ChangeMean);
	pPriors->Mat[0][2] = (float)(time2ChangeVar);
	pPriors->Mat[1][0] = 1 - pPriors->Mat[0][0];
	pPriors->Mat[1][1] = 0.0;
	pPriors->Mat[1][2] = 0.0;
	for (f = 0; f < 3; f++) {
		pPriors->Mat[2][f] = (float)(changeMean[f]);
		pPriors->Mat[3][f] = (float)(changeVariance[f]);
		pPriors->Mat[4][f] = (float)(noChangeMean[f]);
		pPriors->Mat[5][f] = (float)(noChangeVariance[f]);
	}
	pSaver->saveFeature(this->pTweak->segmentationChangesLz.priorsFileName, pPriors);

	//sclib::matrixOut("priors.txt", pPriors->Mat, 6, 3, this->pTweak);
	//sclib::vectorOut("timesToChange.txt", time2Change, changeCount+1, false, this->pTweak);

	//better approximation of the time2change density by a gmm of somehow optimal size
	pTime2ChangeData = new SV_Data(changeCount+1, 1);
	for (int l = 0; l < pTime2ChangeData->Row; l++) {
		pTime2ChangeData->Mat[l][0] = (float)(time2Change[l]);
	}
	bestOrder = pModelHandler->guessModelOrder(pTime2ChangeData, NULL, sclib::mtGMM_new, 1, 32);
	pTime2ChangeModel = pModelHandler->buildModel(pTime2ChangeData, NULL, bestOrder, sclib::mtGMM_new);
	pModelHandler->saveModel(this->pTweak->segmentationChangesLz.time2ChangeModelFileName, pTime2ChangeModel, false);

	//clean up
	MFree_0D(pTime2ChangeData);
	MFree_0D(pTime2ChangeModel);
	MFree_0D(pPriors);
	MFree_0D(pSaver);
	MFree_2D(changeDist);
	MFree_2D(noChangeDist);
	MFree_1D(changeMean);
	MFree_1D(changeVariance);
	MFree_1D(noChangeMean);
	MFree_1D(noChangeVariance);
	MFree_1D(time2Change);
	MFree_0D(pModelHandler);

	return SVLIB_Ok;
}

//====================================================================================================================
//	detect changes in the characteristics of the audio-stream
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
//====================================================================================================================
int SC_Segmentation_Changes_LZW::detectChanges(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	int idxMFCC = sclib::bitPosition(sclib::featureMFCC), idxLSP = sclib::bitPosition(sclib::featureLSP), idxPitch = sclib::bitPosition(sclib::featurePitch);
	
	if (pFeatures[idxMFCC] == NULL || pFeatures[idxLSP] == NULL || pFeatures[idxPitch] == NULL ||
		  pFeatures[idxMFCC]->Hdr.frameSize != pFeatures[idxLSP]->Hdr.frameSize ||
			pFeatures[idxMFCC]->Hdr.frameStep != pFeatures[idxLSP]->Hdr.frameStep ||
			pFeatures[idxMFCC]->Hdr.sampleRate != pFeatures[idxLSP]->Hdr.sampleRate ||
			pFeatures[idxPitch]->Hdr.frameSize != pFeatures[idxLSP]->Hdr.frameSize ||
			pFeatures[idxPitch]->Hdr.frameStep != pFeatures[idxLSP]->Hdr.frameStep ||
			pFeatures[idxPitch]->Hdr.sampleRate != pFeatures[idxLSP]->Hdr.sampleRate) {
    return SVLIB_Fail; //pitch, mfccs and lsps must be extracted on the same frame-basis, otherwise return error
  } else {
		if (this->mode == sclib::modeSpeakerChange) {
			return LuZhangWindow(pGT, segmentStart, segmentEnd, pFeatures, sclib::atSpeech, sclib::atPause|sclib::atSilence, sclib::atSpeakerBoundary);
		} else {
			return LuZhangWindow(pGT, segmentStart, segmentEnd, pFeatures, sclib::atNoise, sclib::noType, sclib::atNoiseBoundary);
		}
  }
}

//====================================================================================================================
//  Algorithm (scd-part) from the paper "Unsupervised speaker segmentation and tracking in real-time audio content 
//  analysis", Lu, Zhang, 2005
//  SNR-dependant models are omitted. 
//
//	return-value: # of changes
//====================================================================================================================
unsigned int SC_Segmentation_Changes_LZW::LuZhangWindow(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, unsigned long int type, unsigned long int typesNot, unsigned long int boundaryType) {
	unsigned int changes = 0, f, featureIdx[3];
  long int windowCount;
  unsigned long int extra, start, end, oldEnd, changePoint, y, windowLength, windowStep, maxLen, minLen, maxDim, minDim;
	unsigned short int segmentsWithoutChange = 0;
	double adaptiveThresh, **oldCov[3], **actCov[3], dist[3];
  double **distances = NULL; //array to hold: distance to predecessor: [x][0] MFCC-based, [x][1] LSP-based, [x][2] Pitch-based
  unsigned long int **segments = NULL; //array to hold: [x][0] segStart; [x][1] segEnd (in correspondence with the same indices of distances[])
	SV_Data *pWindow[3], *pWindowHook[3], *pPriors;
	SC_Model *pCurrentModel[3], *pTime2ChangeModel;
  SV_Data *pStrippedFeatures[3]; //copied-together features, [0] MFCCs, [1] LSPs, [2] Pitch
	SC_ModelHandler *pBuilder = new SC_ModelHandler(this->pTweak, true); //TODO: verbose?
	SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);
	double bayesianThreshold = this->pTweak->segmentationChangesLz.bayesianThreshold, P_H0, P_H1;
	
	//some initializations
	for (f = 0; f < 3; f++) {
		oldCov[f] = NULL;
		actCov[f] = NULL;
		pWindow[f] = NULL;
		pWindowHook[f] = NULL;
		pCurrentModel[f] = NULL;
	}
	featureIdx[0] = sclib::bitPosition(sclib::featureMFCC);
	featureIdx[1] = sclib::bitPosition(sclib::featureLSP);
	featureIdx[2] = sclib::bitPosition(sclib::featurePitch);
	pPriors = pFeatureHandler->loadFeature(this->pTweak->segmentationChangesLz.priorsFileName);
	pTime2ChangeModel = pBuilder->loadModel(this->pTweak->segmentationChangesLz.time2ChangeModelFileName, sclib::mtGMM_new);

	//because we want do detect boundarys here, delete all previously set boundary-labels.
	pGT->setSegment(segmentStart, segmentEnd, boundaryType, false, sclib::noSpeaker, sclib::modeLabelRemove);
  
  //mark too short segments as too short (doesn't matter if it has been done before)
	pGT->markShortSpeechSegments(segmentStart, segmentEnd, sclib::max(1, pGT->getConverter()->ms2sample(this->pTweak->general.shortSpeechThreshold, sclib::alignmentEnd)));

  //collect all desired frames in this segment, with their frameNr in the zeroth col.
	//because all 3 feature-sets must have equal frame sizes and steps, the ones from mfccs are used here and later on without loss of generality
  segmentEnd = sclib::min(segmentEnd, pGT->getConverter()->audioFrame2sample(pGT->getLastAudioFrameNrInSegment(segmentStart, segmentEnd, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep), pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd));
  pStrippedFeatures[0] = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureMFCC)], segmentStart, segmentStart, segmentEnd, type, typesNot, 1, 1, true);
	pStrippedFeatures[1] = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureLSP)], segmentStart, segmentStart, segmentEnd, type, typesNot, 1, 1, true);
	pStrippedFeatures[2] = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featurePitch)], segmentStart, segmentStart, segmentEnd, type, typesNot, 1, 1, true);
	if (pStrippedFeatures[0] != NULL) { //no frames of desired type here

		//TODO
		//SC_SignalHandler sigWriter(this->pTweak, sclib::stWave);
		//sigWriter.storeSignal("speech.wav", segmentStart, segmentEnd, pGT->getSignalPrototype(), pGT, type, sclib::noType, sclib::modeHypothesized, sclib::atSpeakerBoundary, false, false);

		//determine window parameters
 		windowLength = pGT->getConverter()->ms2audioFrame(this->pTweak->segmentationChangesLz.detectorWindowLength, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd); //measured in audioFrames
		windowStep = pGT->getConverter()->ms2audioFrame(this->pTweak->segmentationChangesLz.detectorWindowStep, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd); //measured in audioFrames
		windowCount = pStrippedFeatures[0]->Row / windowStep;
		windowCount -= (unsigned long)ceil((double)(((windowCount-1)*windowStep + windowLength) - pStrippedFeatures[0]->Row) / (double)windowStep); //subtract all frames which overlap the amount of frames
		extra = pStrippedFeatures[0]->Row - (sclib::max(0, windowCount-1)*windowStep + windowLength) - 1; //the last window can be extra long, if there are some frames which are too less to make an extra window
		if (windowCount <= 0) {
			MFree_0D(pStrippedFeatures[0]);
			MFree_0D(pStrippedFeatures[1]);
			MFree_0D(pStrippedFeatures[2]);
			return SVLIB_Ok;
		}

		MArray_2D(distances, windowCount, 3, double, "SC_Segmentation_Changes_LZW.LuZhangWindow: distances"); //make it so big that it can hold the max. possible count of segments
		MArray_2D(segments, windowCount, 2, unsigned long int, "SC_Segmentation_Changes_LZW.LuZhangSegment: segments"); //make it so big that it can hold the max. possible count of segments

		//loop over all subsegments (x second windows with y second overlap), check if they're valid (long enough etc.)
		//compute divergence-shape distance for each two windows of lsps here
		for (y = 0; y < (unsigned long)windowCount; y++) {
			minLen = y * windowStep; //select a proper submatrix that contains the current window; omit the first column (contains only the frame-number)
			maxLen = minLen + windowLength + ((y == windowCount-1) ? extra : 0);
			minDim = 1;

			for (f = 1; f < 2; f++) {
				if (actCov[f] != NULL) { //save last valid cov
					MFree_2D(oldCov[f]); 
					oldCov[f] = actCov[f];
				}
				maxDim = pStrippedFeatures[f]->Col; 
				actCov[f] = this->pMatrixFunc->cov(pStrippedFeatures[f]->Mat, pStrippedFeatures[f]->Row, pStrippedFeatures[f]->Col, NULL, minLen, maxLen, minDim, maxDim); //compute new cov

				//fill distances
				segments[y][0] = (unsigned long)pStrippedFeatures[f]->Mat[minLen][0]; //original start-sample of first frame in the window
				segments[y][1] = (unsigned long)pStrippedFeatures[f]->Mat[maxLen][0] + pGT->getConverter()->audioFrame2sample(0, pStrippedFeatures[f]->Hdr.frameSize, pStrippedFeatures[f]->Hdr.frameStep, sclib::alignmentEnd); //original end-sample of last frame in the window
				if (y == 0) { //this is the first (valid) window in this scene, it has no predecessor to compute distance
					distances[y][f]	= 0.0;
				} else { //this window has a predecessor, so compute distance
					distances[y][f]	= this->pDist->divergenceShape(actCov[f], oldCov[f], maxDim-minDim);
					if (distances[y][f] <= 0.0) {
						//REPORT_ERROR(SVLIB_Fail, "Matrix Inversion failed... maybe too short window?");
						//distances[y][f] = 0;
						sclib::matrixOut("actCov.txt", actCov, maxDim-minDim, maxDim-minDim, this->pTweak);
						sclib::matrixOut("oldCov.txt", oldCov, maxDim-minDim, maxDim-minDim, this->pTweak);
					}
				}
			}

			///*
			//if (changePoints2 > 0) {
				//SC_SignalHandler sigWriter(this->pTweak, sclib::stWave);
				//char buff[sclib::bufferSize];
				//sprintf(buff, "toSee_%04d_%d_%d_%d.wav", y, changePoints2, pGT->getConverter()->sample2videoFrame(segments[y][0]), pGT->getConverter()->sample2videoFrame(segments[y][1]));
				//sigWriter.storeSignal(buff, segments[y][0], segments[y][1], pGT->getSignalPrototype(), pGT, type, typesNot);
				//sprintf(buff, "real_%04d_%d_%d_%d.wav", y, changePoints2, pGT->getConverter()->sample2videoFrame(segments[y][0]), pGT->getConverter()->sample2videoFrame(segments[y][1]));
				//sigWriter.storeSignal(buff, segments[y][0], segments[y][1], pGT->getSignalPrototype(), pGT, sclib::noType, sclib::noType);
			//}
			//*/

		} //for all windows
		for (f = 0; f < 3; f++) {
			MFree_0D(pStrippedFeatures[f]);
			MFree_2D(oldCov[f]);
			MFree_2D(actCov[f]);
		}

		if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
			sclib::stringOut("scd_debug.txt", "videoFramStart; videoFrameEnd; LSP-distance; adaptiveThreshold; trueChangePoint; likelihoodRatio; hypothesizedThreshold", this->pTweak, "\n");
		}

		//step trough the distances-array to detect speaker-changes
		for (y = 0; y < (unsigned long)windowCount; y++) {
			start = segments[y][0];
			end = segments[y][1];
			oldEnd = (y == 0) ? 0 : segments[y-1][1];
			adaptiveThresh = adaptiveThreshold(distances, 1, y, sclib::min(segmentsWithoutChange, this->pTweak->segmentationChangesLz.lastNdistances));

			if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
				changePoint = (pGT->isSpeakerHomogenious(y > 0 ? segments[y-1][0] : segments[y][0], segments[y][1], sclib::modeGroundtruth) == false) ? 1 : 0; //gt-information if there is a true changepoint in this window; just for debug-output purposes
				sclib::scalarOut("scd_debug.txt", pGT->getConverter()->sample2videoFrame(start), this->pTweak, false, "; ");
				sclib::scalarOut("scd_debug.txt", pGT->getConverter()->sample2videoFrame(end), this->pTweak, false, "; ");
				sclib::scalarOut("scd_debug.txt", distances[y][1], this->pTweak, false, "; ");
				sclib::scalarOut("scd_debug.txt", adaptiveThresh, this->pTweak, false, "; ");
				sclib::scalarOut("scd_debug.txt", changePoint, this->pTweak, false, "; ");
			}

			for (f = 0; f < 3; f++) {
				if (pWindow[f] == NULL) { //pWindows always points to the first window in an ongoing speaker-homogenious segment; pWindowHook alway points to the segment currently under consideration, with a possible change-point directly preceeding it
					pWindow[f] = pGT->copyFramesTogether(pFeatures[featureIdx[f]], segmentStart, start, end, type, typesNot);
					pWindowHook[f] = pWindow[f];
				} else {
					pWindowHook[f]->Next = pGT->copyFramesTogether(pFeatures[featureIdx[f]], segmentStart, start, end, type, typesNot);
					pWindowHook[f] = pWindowHook[f]->Next;
				}
			}
			
			if (y == 0) { //this is the first valid segment, so we must assume a speaker-change and create a new model here!
				for (f = 0; f < 3; f++) {
					pCurrentModel[f] = pBuilder->buildModel(pWindow[f], NULL, 32, sclib::mtQGMM, 1);
				}
				segmentsWithoutChange = 1;
				changePoint = refineChangePoint(pGT, start, end, oldEnd, type);
				pGT->setSegment(changePoint, changePoint, boundaryType|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //this boundary is artificial because it is set in any case (becasue we don't know any previous segment to which a distance could be computed)
				if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
					sclib::scalarOut("scd_debug.txt", 0, this->pTweak, false, "; "); //ratio 
					sclib::scalarOut("scd_debug.txt", 0, this->pTweak, false, "\n"); //hypothesized boundary
				}

			}	else { //this not the first segment
				
				//do we have a potential speaker change at the current position?
				//if ((((long)y < windowCount-1) && (distances[y][1] > adaptiveThresh) && (distances[y][1] > distances[y-1][1]) && (distances[y][1] > distances[y+1][1])) || //the standart case
				//		(((long)y == windowCount-1) && (distances[y][1] > adaptiveThresh) && (distances[y][1] > distances[y-1][1]))) {	//we examine the last segment, so it can't be compared with the next one
				if (distances[y][1] > adaptiveThresh) {
					
					//Bayesian fusion scheme using prior probabilities from training data
					P_H0 = pPriors->Mat[0][0];
					P_H1 = pPriors->Mat[1][0];
					for (f = 0; f < 3; f++) {
						dist[f] = pBuilder->testModel(pCurrentModel[f], pWindowHook[f], 1);
						P_H0 *= ((SC_MixtureModel_GMM*)pTime2ChangeModel)->gaussSolver.gaussian(dist[f], pPriors->Mat[2][f], pPriors->Mat[3][f]);
						P_H1 *= ((SC_MixtureModel_GMM*)pTime2ChangeModel)->gaussSolver.gaussian(dist[f], pPriors->Mat[4][f], pPriors->Mat[5][f]);
					}
					if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
						sclib::scalarOut("scd_debug.txt", P_H0/P_H1, this->pTweak, false, "; "); //ratio
					}

					if (P_H0/P_H1 >= bayesianThreshold) {
						for (f = 0; f < 3; f++) {
							pWindow[f] = sclib::destructLinkedList(pWindow[f], segmentsWithoutChange); //kill data of last speaker's segment
							MFree_0D(pCurrentModel[f]);
							pCurrentModel[f] = pBuilder->buildModel(pWindow[f], NULL, 32, sclib::mtQGMM, 1);
						}
						changes++;
						segmentsWithoutChange = 1;
						changePoint = refineChangePoint(pGT, start, end, oldEnd, type);
						pGT->setSegment(changePoint, changePoint, boundaryType, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
						if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
							sclib::scalarOut("scd_debug.txt", 1, this->pTweak, false, "\n");
						}

					}	else { //no real change, update current model
						for (f = 0; f < 3; f++) {
							pCurrentModel[f]->TrainModel(pWindow[f], 1);
						}
						segmentsWithoutChange++;
						if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
							sclib::scalarOut("scd_debug.txt", 0, this->pTweak, false, "\n");
						}
					}

				} else { //no potential change
					for (f = 0; f < 3; f++) {
						pCurrentModel[f]->TrainModel(pWindowHook[f], 1);
					}
					segmentsWithoutChange++;
					if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true && sclib::bitTest(type, sclib::atSpeech) == true) {
						sclib::scalarOut("scd_debug.txt", 0, this->pTweak, false, "; "); //ratio
						sclib::scalarOut("scd_debug.txt", 0, this->pTweak, false, "\n");
					}
				} //not potential change
			} //not the first segment
		} //for all distances

	} else {//pStrippedFeatures[0] != NULL 
		for (f = 0; f < 3; f++) {
			MFree_0D(pStrippedFeatures[f]);
			MFree_2D(oldCov[f]);
			MFree_2D(actCov[f]);
		}
	}
	
	for (f = 0; f < 3; f++) {
		sclib::destructLinkedList(pWindow[f], segmentsWithoutChange);
		MFree_0D(pCurrentModel[f]);
	}
	MFree_2D(distances);
  MFree_2D(segments);
	MFree_0D(pBuilder);
	MFree_0D(pFeatureHandler);
	MFree_0D(pPriors);
	MFree_0D(pTime2ChangeModel);

	return changes;
}

//====================================================================================================================
//	calculates the adaptive treshold according to the paper
//	'Real-Time Unsupervised Speaker Change Detection', L.Lu, H.J.Zhang, 2002 (IEEE)
//
//	it takes into account the last min(maxN, currSeg) distances to compute the threshold, alpha is a amplifying factor
//
//	formula: Th = alpha * 1/min(maxN, currSeg) * SUM(distances[n]),
//					 SUM from n=0 to min(maxN, currSeg)
//====================================================================================================================
double SC_Segmentation_Changes_LZW::adaptiveThreshold(double** distances, unsigned int column, unsigned int currSeg, unsigned int maxN) {
	double alpha = this->pTweak->segmentationChangesLz.adaptiveThresholdAlpha;
	double threshold = 0.0;
	unsigned int N = maxN;

	if (currSeg < N) {
		N = currSeg;
	}

	if (N > 0) {
		for (unsigned int n = 1; n <= N; n++) {
			threshold += distances[currSeg - n + 1][column];
		}
		threshold *= alpha / (double)N;
	}

	return threshold;
}

//====================================================================================================================
// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
// are hadnled therein)
//====================================================================================================================
unsigned long int SC_Segmentation_Changes_LZW::getUncertaintyRegionWidth(void) {
	unsigned long int width = 0;

	//factors due to frame-based analysis
	width += 2 * this->pTweak->featureMfcc.frameStep;

	//factors due to windowed analysis
	//width += 2 * this->pTweak->segmentationChangesLz.detectorWindowStep;
	width += 2 * this->pTweak->segmentationChangesLz.detectorWindowLength; //TODO

	return width;
}
