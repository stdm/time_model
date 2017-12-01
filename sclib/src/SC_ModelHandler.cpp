/**************************************************************************/
/*    Responsibility:																											*/
/*		  - provides the possbility to build speaker- or noise models of    */
/*        all implemented types                                           */
/*      - helps to save/load models to/from a file                        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.02.2006																								*/
/**************************************************************************/

#include <limits.h>
#include "SC_ModelHandler.h"
#include "SC_Model_qGMM.h"
#include "SC_MixtureModel_MIXMAX.h"
#include "SC_MixtureModel_bGMM.h"
#include "SC_MixtureModel_GMM.h"
#include "SC_MixtureModel_GMM_UBM.h"
#include "SC_MixtureModel_MIX2MAX.h"
#include "SC_MixtureModel_MIX2MAXex.h"
#include "SC_Model_Pareto.h"
#include "SC_Model_SVM.h"
#include "SC_Model_HMM.h"
#include "SC_Model_FullGauss.h"
#include "SC_Model_VQ.h"
#include "SC_Model_Time.h"
#include "SC_Model_MetaGMM.h"
#include "SC_SignalHandler.h"
#include "SC_FeatureHandler.h"


//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_ModelHandler::SC_ModelHandler(SC_TweakableParameters *pTweak, bool verbose) {
  this->pTweak = pTweak;
  this->verbose = verbose;
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_ModelHandler::~SC_ModelHandler() {

}

//TODO: very video-task like, should be moved to a corpus-class?!? -> also usable for speaker-id on TIMIT data, maybe just rename...
//====================================================================================================================
//  Speaker Modelling with Noise-Compensation
//
//  ExplicitModels is a table whichs rows have the form 'sceneNr|segmentNr|fileName of explicit background model'.
//  SceneNr thereby is the number provided by the calling function, and segmentNr is the value of validSegmentCount in
//  this function. The first index gives the row, the second one the columns, the third one is because each entry is 
//  a string of 255 char's. The sceneNr and SegmentNr are only stored in the first byte of this characterfield, so 
//  they must be between 0 and 127. If explicitModelCount >0, explicitModels is parsed, and for a segment specified by 
//  scene- and  segmentNr (0 is a wildcard in  this case: 0 as sceneNr means all scenes, dito for segmentNr), 
//  the background-model will be loaded from the specified file instead of the following (standard) procedure:
//
//  To Compensate for variable background sounds, explicit background-GMM's will be build up from before and after
//	speech as well as from speech-pauses. These will be used to estimate GMM-IB's
//
//  Returned are the speaker-models (together with their respective features and segments-boundarys packed into 
//  cluster-objects) of all valid (long enough) segments as the linked list 'pClusters'; the features and segments-
//  boundarys of all invalid segments are returned in the linked list 'pInvalidClusters' so that overall each single
//  speech-segment of the given segment is found in one of the clusters afterwards
//====================================================================================================================
int	SC_ModelHandler::buildSpeakerModels(SC_Corpus* pCorpus, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data* pFeatures, SC_Cluster* &pClusters, SC_Cluster* &pInvalidClusters, char*** explicitModels, unsigned short int explicitModelCount) {
	unsigned long int x, y, sceneNr, speakerID, validSegmentCount = 0;
	unsigned short int dim;
	unsigned long int lastFrame = pCorpus->getGT()->getLastAudioFrameNrInSegment(segmentStart, segmentEnd, pFeatures->Hdr.frameSize, pFeatures->Hdr.frameStep), lastSample = pCorpus->getGT()->getConverter()->audioFrame2sample(lastFrame, pFeatures->Hdr.frameSize, pFeatures->Hdr.frameStep, sclib::alignmentEnd);
	long int speechStart, speechEnd = segmentStart, oldSpeechEnd, newSpeechEnd, noiseStart, noiseEnd, overallStart, overallEnd;
  SV_Data *pNoise, *pSpeech;
  SC_Model *pBackgroundModel = NULL;
	SC_Cluster *pFirstCluster = NULL, *pActualCluster, *pLastCluster, *pFirstInvalidCluster = NULL, *pLastInvalidCluster;
  SC_SignalHandler *pSaver = new SC_SignalHandler(this->pTweak, pCorpus->getGT()->getSignalPrototype()->getSignalType());
	SC_Signal *pSig2Save = NULL;
	SC_FeatureHandler extractor(this->pTweak, this->verbose);
	char speakerName[sclib::bufferSize], *fileName;
  double snr;
  int res = 0;

  assert(pFeatures != NULL);
  assert(pCorpus->getGT() != NULL);
  assert(segmentStart <= segmentEnd);
  //assert(this->pTweak->modelHandler.foregroundModelType == sclib::mtGMM_UBM || this->pTweak->modelHandler.foregroundModelType == sclib::mtGMM_new || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIXMAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtBGMM || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX_ex || this->pTweak->modelHandler.foregroundModelType == sclib::mtPareto); 
  assert(this->pTweak->modelHandler.backgroundModelType == sclib::mtGMM_UBM || this->pTweak->modelHandler.backgroundModelType == sclib::mtGMM_new || this->pTweak->modelHandler.backgroundModelType == sclib::mtPareto); 

  dim = pFeatures->Col;
  lastSample = min(lastSample, segmentEnd);

  //mark too short segments as too short (doesn't matter if it has been done before)
	pCorpus->getGT()->markShortSpeechSegments(segmentStart, segmentEnd, sclib::max(1, pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->general.shortSpeechThreshold)));

	for (y = segmentStart; y <= lastSample; y++) {
    oldSpeechEnd = speechEnd;

		pCorpus->getGT()->getNextBoundary(y, speechStart, speechEnd, sclib::atSpeakerBoundary, sclib::searchForward);
    if ((speechStart != sclib::noSegment) && (speechEnd != sclib::noSegment) && (speechStart < (long)(lastSample))) {				
      snr = 120.0;
			sceneNr = pCorpus->getGT()->sample2scene(speechStart);
      
      //if there was no further speaker-boundary in this scene, we must manually set the speech-end
			if (speechEnd > (long)(lastSample)) {speechEnd = (long)(lastSample);}

			//find the real end of this speech-segment
			speechEnd = pCorpus->getGT()->getSpeechEnd(speechEnd);
			if (speechEnd < speechStart) {speechEnd = speechStart;}

      //shall the segment of this speaker be evaluated?
      if (strncmp(this->pTweak->modelHandler.onlyThisSpeaker, "", sclib::bufferSize) == 0 || pCorpus->getGT()->isSpeakerHomogenious(speechStart, speechEnd, this->pTweak->modelHandler.onlyThisSpeaker) == true) {

			  //get startpoint of noise before the actual speech-segment
				noiseStart = (y == segmentStart) ? segmentStart : oldSpeechEnd;
			  if (noiseStart < (long)(speechStart - pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->modelHandler.msPerGaussian))) {
          noiseStart = speechStart - pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->modelHandler.msPerGaussian);
        }

        //get endpoint of noise after the actual speech-segment
			  pCorpus->getGT()->getNextBoundary(speechEnd+1, noiseEnd, newSpeechEnd, sclib::atSpeakerBoundary, sclib::searchForward);
				noiseEnd = (noiseEnd != sclib::noSegment) ? noiseEnd-1 : lastSample;
			  if (noiseEnd > (long)(lastSample)) {
          noiseEnd = (long)(lastSample);
        }
			  if (noiseEnd > (long)(speechEnd + pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->modelHandler.msPerGaussian))) {
          noiseEnd = (long)(speechEnd + pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->modelHandler.msPerGaussian));
        }

			  //copy all feature-vectors of the before-, between- and after-speech-noise-part together; TODO: types ok???
        pNoise  = pCorpus->getGT()->copyFramesTogether(pFeatures, segmentStart, noiseStart, noiseEnd, sclib::atNoise|sclib::atSilence); //, sclib::atUnvoiced|sclib::atPause);
			  pSpeech = pCorpus->getGT()->copyFramesTogether(pFeatures, segmentStart, speechStart, speechEnd, sclib::atSpeech); //, sclib::atShort|sclib::atPause|sclib::atSilence|sclib::atNoise|sclib::atUnvoiced);
  		
        if (pSpeech != NULL) {
          if (pSpeech->Row > (long)(pCorpus->getGT()->getConverter()->ms2audioFrame(this->pTweak->modelHandler.msPerGaussian, pFeatures->Hdr.frameSize, pFeatures->Hdr.frameStep))) { //minimum duration for models
						//sclib::matrixOutEx("speech.txt", pSpeech->Mat, pSpeech->Row, pSpeech->Col, this->pTweak); 
					  validSegmentCount++; //call the first valid segment '1', not '0'

            //only the gmm-ib/wesley-model/pareto-model needs a background-model here
            if (this->pTweak->modelHandler.foregroundModelType == sclib::mtMIXMAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX_ex || this->pTweak->modelHandler.foregroundModelType == sclib::mtPareto) {
              
              //is there an explicit model to load?
					    if (explicitModels != NULL) {
						    for (x = 0; x < explicitModelCount; x++) {
						      if (((explicitModels[x][0][0] == sceneNr)           || (explicitModels[x][0][0] == 0)) && 
							        ((explicitModels[x][1][0] == validSegmentCount) || (explicitModels[x][1][0] == 0))) {
							      pBackgroundModel = loadModel(explicitModels[x][2], this->pTweak->modelHandler.backgroundModelType);
								    break;
							    }
						    } 
					    }

					    //build background-model, if possible
					    if (pBackgroundModel == NULL) {
						    if (pNoise != NULL) {
							    snr = getSNR(pSpeech, pNoise, ((pSpeech->Hdr.ID != sclib::featureMFCC && pSpeech->Hdr.ID != sclib::featureSTE) ? pSpeech->Col-1 : 0)); //for all but mfcc/energy features, all cols must be summed up in order to get energy
							    if (snr <= pTweak->modelHandler.SNRthreshold) {
										if (this->pTweak->modelHandler.outlierRemovalMode != sclib::outlierRemoveNone) {
											printf("\n(%#6.5f%% outliers removed from noise data for model %d)", extractor.removeOutliers(pNoise, this->pTweak->modelHandler.outlierRemovalMode), speakerID);
										}
                    pBackgroundModel = buildModel(pNoise, NULL, sclib::modeBackground, 1);
								    if (pBackgroundModel != NULL && sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbWav)) {
									    fileName = new char[sclib::bufferSize];
									    sprintf(fileName, "%s_%04d%s", "noise", sceneNr*100 + validSegmentCount, ".wav");
											pSig2Save = pCorpus->loadSignal((unsigned long&)noiseStart, (unsigned long&)noiseEnd, true);
                      pSaver->storeSignal(fileName, noiseStart, noiseEnd, pSig2Save, pCorpus->getGT(), sclib::atNoise|sclib::atSilence, sclib::atUnvoiced|sclib::atPause);
											MFree_0D(pSig2Save);
									    MFree_1D(fileName);
								    }
							    }
						    }
					    }
            } //modelType == MIXMAX/MIX2MAX/Pareto

            //build a dummy model if none has been build before
            if (pBackgroundModel == NULL) {
              pBackgroundModel = createRawModel(this->pTweak->modelHandler.backgroundModelType, NULL, 0, dim);
            }

						//remove outliers, if wished
            speakerID = sceneNr*100 + validSegmentCount;
						if (this->pTweak->modelHandler.outlierRemovalMode != sclib::outlierRemoveNone) {
							printf("\n(%#6.5f%% outliers removed from speech data for model %d)", extractor.removeOutliers(pSpeech, this->pTweak->modelHandler.outlierRemovalMode), speakerID);
						}
					  sprintf(speakerName, "%4d\0", sceneNr*100 + validSegmentCount);
            pActualCluster = buildModel(pSpeech, speechStart, speechEnd, pBackgroundModel, sclib::modeForeground, speakerName, speakerID);
						if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSpeechTrainData) == true) {
							sclib::matrixOut("speechTrainData.txt", pSpeech->Mat, pSpeech->Row, pSpeech->Col, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
						}
            res++;

            MFree_0D(pBackgroundModel);
            if (this->verbose == true) {printf(".");}

            //print the SNR, if desired
            if (this->pTweak->debug.debugMode & sclib::dbSNR) {sclib::scalarOutEx("SNR.txt", snr, speakerName, this->pTweak);}

            //make the collected sounds hearable, if desired
            if (this->pTweak->debug.debugMode & sclib::dbWav) {
						  fileName = new char[sclib::bufferSize];
						  sprintf(fileName, "%s_%04d%s", "complete", sceneNr*100 + validSegmentCount, ".wav");
							overallStart = sclib::min(noiseStart, speechStart);
							overallEnd = sclib::max(noiseEnd, speechEnd);
							pSig2Save = pCorpus->loadSignal((unsigned long&)overallStart, (unsigned long&)overallEnd, true);
              pSaver->storeSignal(fileName, overallStart, overallEnd, pSig2Save, pCorpus->getGT());
							MFree_0D(pSig2Save);
						  MFree_1D(fileName);			
						  fileName = new char[sclib::bufferSize];
						  sprintf(fileName, "%s_%04d%s", "speech", sceneNr*100 + validSegmentCount, ".wav");
							pSig2Save = pCorpus->loadSignal((unsigned long&)speechStart, (unsigned long&)speechEnd, true);
              pSaver->storeSignal(fileName, speechStart, speechEnd, pSig2Save, pCorpus->getGT(), sclib::atSpeech, sclib::atPause|sclib::atUnvoiced|sclib::atSilence|sclib::atNoise);
							MFree_0D(pSig2Save);
						  MFree_1D(fileName);
					  }

            //build a linked list of all clusters of this scene			
            if (pActualCluster != NULL) {
              if (validSegmentCount == 1) {
						    pFirstCluster = pActualCluster;
						    pLastCluster = pFirstCluster;
					    } else {
						    pLastCluster->Next = pActualCluster;
						    pLastCluster = pLastCluster->Next;
					    }
            } else { //actual-cluster==NULL => speakermodel-order==0
						  validSegmentCount--;
						  MFree_0D(pSpeech);
						  pCorpus->getGT()->setSegment(speechStart, speechEnd, sclib::atShort);
						  pCorpus->getGT()->setSegment(speechStart, speechStart, sclib::atSpeechSegmentStart);
						  pCorpus->getGT()->setSegment(speechEnd, speechEnd, sclib::atSpeechSegmentEnd);
					  }
				  }	else { //speech not long enough to estimate models and build an initial cluster
            MFree_0D(pSpeech);
					  pCorpus->getGT()->setSegment(speechStart, speechEnd, sclib::atShort);
					  pCorpus->getGT()->setSegment(speechStart, speechStart, sclib::atSpeechSegmentStart);
					  pCorpus->getGT()->setSegment(speechEnd, speechEnd, sclib::atSpeechSegmentEnd);
          }
			  } else { //speech does not exist
				  pCorpus->getGT()->setSegment(speechStart, speechEnd, sclib::atShort);
				  pCorpus->getGT()->setSegment(speechStart, speechStart, sclib::atSpeechSegmentStart);
				  pCorpus->getGT()->setSegment(speechEnd, speechEnd, sclib::atSpeechSegmentEnd);
			  }
        
        MFree_0D(pSpeech);
        MFree_0D(pNoise);
      } else { //this speaker-name shall be regarded
        if (this->verbose == true) {printf("\n      Segment %i-%i skipped (not from '%s')", speechStart, speechEnd, this->pTweak->modelHandler.onlyThisSpeaker);}
      }

      y = speechEnd;
		} else {break;}
	}

  //collect the information of the too-short segments
  //this can't be done in the above ...else... branches, because there are short segments labeled as such in previous algorithms, so it has to be done here to collect all such segments
  for (unsigned long int y = segmentStart; y <= segmentEnd; y++) {
		pCorpus->getGT()->getNextSegment(y, speechStart, speechEnd, sclib::atShort);
		
    if ((speechStart != sclib::noSegment) && (speechEnd != sclib::noSegment) && (speechStart <= (long)(segmentEnd))) {
			sceneNr = pCorpus->getGT()->sample2scene(speechStart);

      if (speechEnd > (long)(segmentEnd)) {
        speechEnd = (long)(segmentEnd);
      }

      //shall the segment of this speaker be evaluated?
      if (strncmp(this->pTweak->modelHandler.onlyThisSpeaker, "", sclib::bufferSize) == 0 || pCorpus->getGT()->isSpeakerHomogenious(speechStart, speechEnd, this->pTweak->modelHandler.onlyThisSpeaker) == true) {
        pSpeech = pCorpus->getGT()->copyFramesTogether(pFeatures, segmentStart, speechStart, speechEnd, sclib::atSpeech, sclib::atPause|sclib::atSilence|sclib::atNoise|sclib::atUnvoiced);
        if (pSpeech != NULL) {
          validSegmentCount++;
          speakerID = sceneNr*100 + validSegmentCount;
				  sprintf(speakerName, "%4d\0", sceneNr*100 + validSegmentCount);
          pActualCluster = new SC_Cluster(this->pTweak, pSpeech, (unsigned long*)&speechStart, (unsigned long*)&speechEnd, NULL, NULL, NULL, 1, speakerID, speakerName, false); //a cluster just to pack segment-boundarys and features together, no models!
      		
          //build a linked list of all clusters of this scene			
          if (pActualCluster != NULL) {
            if (pFirstInvalidCluster == NULL) {
					    pFirstInvalidCluster = pActualCluster;
					    pLastInvalidCluster = pFirstInvalidCluster;
				    } else {
					    pLastInvalidCluster->Next = pActualCluster;
					    pLastInvalidCluster = pLastInvalidCluster->Next;
				    }
          }

          MFree_0D(pSpeech);      
        }
      }

      y = speechEnd; //normally y=segmentEnd+1, but the next loop increments y for us automatically		
    } else {
      break;
    }
  }

	//create the segment-statistics
	if (this->pTweak->debug.debugMode & sclib::dbSegmentStatistics) {
    pCorpus->getGT()->segmentStatisticsOut("segStat.txt", segmentStart, segmentEnd);
	}

  MFree_0D(pSaver);

	pClusters = pFirstCluster;
  pInvalidClusters = pFirstInvalidCluster;

	return res;
}

//====================================================================================================================
//	Used by buildSpeakerModels(): gets the features, a background-model and a name/id for this segment and constructs 
//  the model, which is returned in a cluster-object (for this reason, the additionalInfo is needed); debug-outputting 
//  is also handled; layer can be sclib::modeForeground or sclib::modeBackground to determine whether a speech- or 
//  background model should be built; segment-start and -end are just for the purpose of storing the start- and 
//  endsamples of the features in the cluster object.
//====================================================================================================================
SC_Cluster* SC_ModelHandler::buildModel(SV_Data *pFeatures, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Model* pBackgroundModel, unsigned long int layer, char* name, unsigned long int ID) {
  SC_Model *pModel = NULL, *pBgModel;
  SC_Cluster *pActualCluster = NULL;
  
  pModel = buildModel(pFeatures, pBackgroundModel, layer, 1);
  if (pModel != NULL) {
    //TODO: background-dummy for pareto?
    //pack all the collected infomation of this segment in a cluster-object
    pBgModel = (pBackgroundModel == NULL) ?  (new SC_MixtureModel_GMM(this->pTweak, 0, pFeatures->Col)) : pBackgroundModel;
	  pActualCluster = new SC_Cluster(this->pTweak, pFeatures, (unsigned long *)(&segmentStart), (unsigned long *)(&segmentEnd), NULL, pBgModel, pModel, 1, ID, name, false);
    MFree_0D(pModel);
    if (pBackgroundModel == NULL) {
      MFree_0D(pBgModel);
    }
  }
  
  return pActualCluster;
}

//====================================================================================================================
//	Used by buildSpeakerModels(): gets the features and maybe a background-model for this segment and constructs 
//  the model, which is returned; debug-outputting is also handled
//  layer can be sclib::modeForeground or sclib::modeBackground to determine whether a speech- or background model 
//  should be built
//====================================================================================================================
SC_Model* SC_ModelHandler::buildModel(SV_Data *pFeatures, SC_Model* pBackgroundModel, unsigned long int layer, unsigned long int segmentsToMerge) {
  unsigned long int res, modelType;
  unsigned short int modelOrder, maxOrder;
  SC_Model *pModel = NULL;
  
  if (layer == sclib::modeForeground) {
    if (this->pTweak->modelHandler.foregroundModelType == sclib::mtMIXMAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX_ex) {
      //TODO: bg-model for pareto
      //don't build MIXMAX (or MIX2MAX) model here if no background-model is available
      if (pBackgroundModel == NULL || (pBackgroundModel != NULL && ((SC_MixtureModel*)pBackgroundModel)->getDim() < 1)) {
        modelType = sclib::mtGMM_new;
      } else {
        modelType = this->pTweak->modelHandler.foregroundModelType;
      }
    } else {
      modelType = this->pTweak->modelHandler.foregroundModelType;
    }
    maxOrder = this->pTweak->modelHandler.maxSpeakerModelOrder;
  } else {
    modelType = this->pTweak->modelHandler.backgroundModelType;
    maxOrder = this->pTweak->modelHandler.maxNoiseModelOrder;
  }

  modelOrder = guessModelOrder(pFeatures, pBackgroundModel, modelType, 1, maxOrder, segmentsToMerge);
  if (modelOrder > 0) {
    pModel = createRawModel(modelType, pBackgroundModel, modelOrder, pFeatures->Col);

    if (pModel != NULL) {
			if (pFeatures != NULL) {
				res = pModel->TrainModel(pFeatures, segmentsToMerge);
			}
			pModel->setBackground(NULL); //the setted background-model (in case of the mixmax) get's destructed, what remains is a copy in the cluster, but with a different address
      
      if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSpeakerModels) && layer == sclib::modeForeground) {
				sclib::classOut("speakerModels.txt", pModel, this->pTweak);
      } else if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbNoiseModels) && layer != sclib::modeForeground) {
				sclib::classOut("noiseModels.txt", pModel, this->pTweak);
      }
    }
  }

  return pModel;
}

//====================================================================================================================
//	Generates a Model for a given segment of feature-vectors of specific type
//  segmentStart end segmentEnd (in samples) must be start- and endsample of the given features!!!
//====================================================================================================================
SC_Model* SC_ModelHandler::buildModel(SC_GroundTruth* pGT, SV_Data* pFeatures, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Model* pBackground, unsigned long int types, unsigned long int typesNot, unsigned short int modelOrder, unsigned short int modelType) {
	unsigned long int end = segmentStart + pGT->getConverter()->audioFrame2sample(pFeatures->Row, pFeatures->Hdr.frameSize, pFeatures->Hdr.frameStep, sclib::alignmentEnd);
  SV_Data* pFrames = pGT->copyFramesTogether(pFeatures, segmentStart, segmentStart, end, types, typesNot);
	SC_Model* pModel = NULL;
	
  if (pFrames != NULL) {
    pModel = createRawModel(modelType, pBackground, modelOrder, pFrames->Col);
  
    if (pModel != NULL) {
	    pModel->TrainModel(pFrames);
    }

	  MFree_0D(pFrames);
  }
	
  return pModel;
}

//====================================================================================================================
//	Generates a Model and returns it
//  in contrast to the abovementioned buildModel()-method, this one can merge all the frames in the linked list 
//  pFeature together and builds the model instead of copying frames of desired type together; if the given 
//  modelOrder==0, a search for the best order is conducted
//====================================================================================================================
SC_Model* SC_ModelHandler::buildModel(SV_Data* pFeatures, SC_Model* pBackground, unsigned short int modelOrder, unsigned long int modelType, unsigned long int segmentsToMerge) {
	int res = SVLIB_Fail;
  SC_Model* pModel = NULL;
	unsigned short int order = modelOrder;

	if (order == 0) {
		order = guessModelOrder(pFeatures, pBackground, modelType, 1, this->pTweak->modelHandler.maxSpeakerModelOrder, segmentsToMerge);
	}

  pModel = createRawModel(modelType, pBackground, order, pFeatures->Col);
  
  if (pModel != NULL) {
    pModel->TrainModel(pFeatures, segmentsToMerge);
  }

	return pModel;
}

//====================================================================================================================
//	Generates and saves a Model for a given segment of feature-vectors of specific type; returns nr. of bytes written
//====================================================================================================================
int SC_ModelHandler::saveModel(const char* fileName, SC_GroundTruth* pGT, SV_Data* pFeatures, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Model* pBackground, unsigned long int types, unsigned long int typesNot, unsigned short int modelOrder, unsigned short int modelType) {
  SC_Model *pModel = buildModel(pGT, pFeatures, segmentStart, segmentEnd, pBackground, types, typesNot, modelOrder, modelType);
  int res = saveModel(fileName, pModel);

  return res;
}

//====================================================================================================================
//	Generates a Model and saves it while returning nr. of bytes written
//  in contrast to the abovementioned buildModel()-method, this one just merges all the frames in the linked list 
//  pFeature together  and  builds the model instead of copying frames of desired type together
//====================================================================================================================
int SC_ModelHandler::saveModel(const char* fileName, SV_Data* pFeatures, SC_Model* pBackground, unsigned short int modelOrder, unsigned short int modelType) {
  SC_Model *pModel = buildModel(pFeatures, pBackground, modelOrder, modelType, 0);
  int res = saveModel(fileName, pModel);
  MFree_0D(pModel);

  return res;
}

//====================================================================================================================
//	save a model to a file; handles report-printing as well; uses debug-dir as a prefix to the filename only if wanted
//====================================================================================================================
int SC_ModelHandler::saveModel(const char *fileName, SC_Model *pModel, bool useDebugDir) {
	int res = SVLIB_Fail;
  char *fNameModel = NULL;
  char *fNameAscii = NULL, *fNameReport = NULL;
	char buffer[sclib::bufferSize];

  if (pModel != NULL) {
		if (useDebugDir == false) { //save debug-dir and set it to "" (necessary because SC_Model::saveModel() uses the debug-dir internally and not all these classes should be altered...); don't use setDebugDir() method because the internal caching isn't wanted here
			sprintf(buffer, "%s", this->pTweak->debug.debugDir);
			sprintf(this->pTweak->debug.debugDir, "");
		}

		fNameModel = new char[strlen(this->pTweak->debug.debugDir) + strlen(fileName) + 1];
		sprintf(fNameModel, "%s%s", this->pTweak->debug.debugDir, fileName);

		if (this->pTweak->debug.debugMode & sclib::dbModelGenerationReport) {
      fNameAscii = sclib::exchangeFileExtension(fileName, ".txt");
      fNameReport = sclib::exchangeFileExtension(fileName, ".report");

      this->pTweak->toggleVerboseMode(true);
			sclib::classOut(fNameReport, this->pTweak, this->pTweak);
			sclib::classOut(fNameAscii, pModel, this->pTweak);

		  MFree_1D(fNameAscii);
		  MFree_1D(fNameReport);
		}

	  pModel->OpenFile(fNameModel, WRITE_MODEL);
	  res = pModel->SaveModel();
	  pModel->CloseFile();

	  MFree_1D(fNameModel);

		if (useDebugDir == false) { //restore previous debug-dir from saved one
			sprintf(buffer, "%s", this->pTweak->debug.debugDir);
			sprintf(this->pTweak->debug.debugDir, "");
		}
  }

	return res;
}

//====================================================================================================================
//	load a model from a file according to the given model type
//====================================================================================================================
SC_Model* SC_ModelHandler::loadModel(const char* fileName, unsigned long int modelType) {
  SC_Model *pModel = NULL, *pBackgroundDummy = NULL;
  
  if (sclib::fileExists(fileName) == true) {
    pModel = createRawModel(modelType, NULL, 0, 1);

    if (pModel != NULL) {
      pModel->OpenFile(fileName, READ_MODEL);
	    pModel->LoadModel();
	    pModel->CloseFile();
    }
  }

  return pModel;
}

//====================================================================================================================
//	Creates a new object of the desired type
//====================================================================================================================
SC_Model* SC_ModelHandler::createRawModel(unsigned long int modelType, SC_Model *pBackgroundModel, unsigned short int modelOrder, unsigned short int dim) {
  SC_Model *pModel = NULL;

  switch (modelType) {
    case (sclib::mtMIXMAX): {
      if (!(pBackgroundModel == NULL || pBackgroundModel->Hdr.ModelType == sclib::mtMIXMAX || pBackgroundModel->Hdr.ModelType == sclib::mtMIX2MAX || pBackgroundModel->Hdr.ModelType == sclib::mtGMM_new || pBackgroundModel->Hdr.ModelType == sclib::mtBGMM || pBackgroundModel->Hdr.ModelType == sclib::mtMIX2MAX_ex)) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model must be a mixture-model");
      } else {
				pModel = new SC_MixtureModel_MIXMAX(this->pTweak, (SC_MixtureModel*)pBackgroundModel, modelOrder, dim); 
      }
      break;
    }
    case (sclib::mtMIX2MAX): {
      if (!(pBackgroundModel == NULL || pBackgroundModel->Hdr.ModelType == sclib::mtMIXMAX || pBackgroundModel->Hdr.ModelType == sclib::mtMIX2MAX || pBackgroundModel->Hdr.ModelType == sclib::mtGMM_new || pBackgroundModel->Hdr.ModelType == sclib::mtBGMM || pBackgroundModel->Hdr.ModelType == sclib::mtMIX2MAX_ex)) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model must be a mixture-model");
      } else {
				pModel = new SC_MixtureModel_MIX2MAX(this->pTweak, (SC_MixtureModel*)pBackgroundModel, modelOrder, dim); 
      }
      break;
    }
    case (sclib::mtMIX2MAX_ex): {
      if (!(pBackgroundModel == NULL || pBackgroundModel->Hdr.ModelType == sclib::mtMIXMAX || pBackgroundModel->Hdr.ModelType == sclib::mtMIX2MAX || pBackgroundModel->Hdr.ModelType == sclib::mtGMM_new || pBackgroundModel->Hdr.ModelType == sclib::mtBGMM || pBackgroundModel->Hdr.ModelType == sclib::mtMIX2MAX_ex)) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model must be a mixture-model");
      } else {
				pModel = new SC_MixtureModel_MIX2MAXex(this->pTweak, (SC_MixtureModel*)pBackgroundModel, modelOrder, dim); 
      }
      break;
    }
    case (sclib::mtBGMM): {
      pModel = new SC_MixtureModel_bGMM(this->pTweak, modelOrder, dim, this->pTweak->mixtureModelBgmm.fullCovariance, verbose);
      break;
    }
    case (sclib::mtGMM_new): {
      pModel = new SC_MixtureModel_GMM(this->pTweak, modelOrder, dim); 
      break;
    }
    case (sclib::mtGMM_UBM): {
      SC_MixtureModel_GMM *pUBM = (SC_MixtureModel_GMM*)loadModel(this->pTweak->mixtureModelGmmubm.ubmFileName, sclib::mtGMM_new);
      if (pUBM != NULL) {
        pModel = new SC_MixtureModel_GMM_UBM(this->pTweak, pUBM);
        MFree_0D(pUBM);
      } else {
        REPORT_ERROR(SVLIB_BadArg, "Can't load specified UBM from file!");
      }
      break;
    }
    case (sclib::mtPareto): {
      if (!(pBackgroundModel == NULL || pBackgroundModel->Hdr.ModelType == sclib::mtPareto)) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model must be a Pareto-model");
      } else {
        pModel = new SC_Model_Pareto(this->pTweak, (SC_Model_Pareto*)pBackgroundModel, this->pTweak->modelPareto.useMarginalDistributions); 
      }
      break;
    }
		case (sclib::mtQGMM): {
			pModel = new SC_Model_qGMM(this->pTweak, dim, modelOrder);
			break;
		}
		case (sclib::mtSVM): {
			pModel = new SC_Model_SVM(this->pTweak, this->pTweak->modelSvm.distanceBasedTesting, this->pTweak->modelSvm.doParameterSearch, this->verbose);
			break;
		}
		case (sclib::mtHMM): {
			pModel = new SC_Model_HMM(this->pTweak, this->pTweak->modelHmm.stateCount, this->pTweak->modelHmm.transitionStructure, this->pTweak->modelHmm.mixturesPerState, this->pTweak->modelHmm.useOrthogonalTransform, this->pTweak->modelHmm.leftToRight, this->pTweak->modelHmm.maxIterations, this->pTweak->modelHmm.verbose);
			break;
		}
		case (sclib::mtFullGauss): {
			pModel = new SC_Model_FullGauss(this->pTweak);
			break;
		}
		case (sclib::mtVQ): {
			pModel = new SC_Model_VQ(this->pTweak, this->pTweak->modelVq.codebookSize, this->pTweak->modelVq.splitMethod, this->pTweak->modelVq.maxIterations, this->verbose);
			break;
		}
		case (sclib::mtTime): {
			pModel = new SC_Model_Time(this->pTweak, dim, this->pTweak->modelTime.syllableLength, this->pTweak->modelTime.trajectoryStep, this->pTweak->modelTime.subModelType, this->pTweak->modelTime.removeTiming, this->pTweak->modelTime.templateCount, this->pTweak->modelTime.clusteringIterations, this->pTweak->modelTime.replaceTrainingData, this->pTweak->modelTime.checkForTrajectorization, this->pTweak->modelTime.worldModelFile, this->pTweak->modelTime.normalizationFile, true, this->verbose);
			break;
		}
		case (sclib::mtMetaGMM): {
			pModel = new SC_Model_MetaGMM(this->pTweak, modelOrder, dim);
			break;
		}
		default: {
      REPORT_ERROR(SVLIB_BadArg, "Unknown model-type");
      break;
    }
  }

  return pModel;
}

//====================================================================================================================
//	Creates a copy of the given parent-model of its own kind
//====================================================================================================================
SC_Model* SC_ModelHandler::copyModel(SC_Model *pParentModel, bool justLink) {
  unsigned long int mt = pParentModel->Hdr.ModelType;
  SC_Model *pModel = NULL;

  switch (mt) {
    case (sclib::mtMIXMAX): {
      pModel = new SC_MixtureModel_MIXMAX((SC_MixtureModel_MIXMAX&)(*pParentModel)); 
      break;
    }
    case (sclib::mtMIX2MAX): {
      pModel = new SC_MixtureModel_MIX2MAX((SC_MixtureModel_MIX2MAX&)(*pParentModel)); 
      break;
    }
    case (sclib::mtMIX2MAX_ex): {
      pModel = new SC_MixtureModel_MIX2MAXex((SC_MixtureModel_MIX2MAXex&)(*pParentModel)); 
      break;
    }
    case (sclib::mtBGMM): {
      pModel = new SC_MixtureModel_bGMM((SC_MixtureModel_bGMM&)(*pParentModel));
      break;
    }
    case (sclib::mtGMM_new): {
      pModel = new SC_MixtureModel_GMM((SC_MixtureModel_GMM&)(*pParentModel));
      break;
    }
    case (sclib::mtGMM_UBM): {
      pModel = new SC_MixtureModel_GMM_UBM((SC_MixtureModel_GMM_UBM&)(*pParentModel));
      break;
    }
    case (sclib::mtPareto): {
      pModel = new SC_Model_Pareto((SC_Model_Pareto&)(*pParentModel)); 
      break;
    }
    case (sclib::mtSVM): {
      pModel = new SC_Model_SVM((SC_Model_SVM&)(*pParentModel), justLink);
      break;
    }
    case (sclib::mtHMM): {
      pModel = new SC_Model_HMM((SC_Model_HMM&)(*pParentModel)); 
      break;
    }
    case (sclib::mtFullGauss): {
      pModel = new SC_Model_FullGauss((SC_Model_FullGauss&)(*pParentModel)); 
      break;
    }
    case (sclib::mtVQ): {
      pModel = new SC_Model_VQ((SC_Model_VQ&)(*pParentModel)); 
      break;
    }
		case (sclib::mtTime): {
      pModel = new SC_Model_Time((SC_Model_Time&)(*pParentModel), justLink);
      break;
    }
		case (sclib::mtMetaGMM): {
			pModel = new SC_Model_MetaGMM((SC_Model_MetaGMM&)(*pParentModel));
			break;
		}
		default: {
      REPORT_ERROR(SVLIB_BadArg, "Unknown model-type");
      break;
    }
  }

  return pModel;
}

//====================================================================================================================
//	Combine 2 models by adding their "components"
//  In contrats to the methods in the different models themselfes, this function "knows" which model can be combined
//  with which one
//====================================================================================================================
SC_Model* SC_ModelHandler::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
  SC_Model *pModel = NULL, *p1 = pFirst, *p2 = pSecond;
  bool exchangePositions = false, ready = true;
  
  //1.: decide which one of the two models should determine the typeof the combined model
  //if the two model have different type, it is important to ecide which one should be the first, 
  //because it determines the type of the combined model, and ome modeltypes incorporated more knowledge/information
  //about the modeled process than others. this first steps therefore decides which model should be the first one
  //according to the "preserve as much information as possible"-criterion
  if (pFirst->Hdr.ModelType == sclib::mtBGMM && pSecond->Hdr.ModelType != sclib::mtBGMM) { 
    //a bgmm can be combined only with a bgmm to form a new bgmm, but it can be combined with any other mixture-model 
    //to form another mixture-model; also, it does not incorporate important new knowledge about the process, just other 
    //trainings-algorithms
    exchangePositions = true;
  } else if (pSecond->Hdr.ModelType == sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
    //a gmm-ib incorporates more knowlege as a standard mixture-model (gmm/bgmm), 
    //so it should be the first one deciding the resultant model's type
    exchangePositions = true;
  } else if (pSecond->Hdr.ModelType == sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
    //the same as gmm-ib (a wesley-model is essentially a gm-ib with different algorithms, but identical meaning)
    exchangePositions = true;
  }  else if (pSecond->Hdr.ModelType == sclib::mtMIX2MAX_ex && pFirst->Hdr.ModelType != sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
    //the same as gmm-ib (a wesleyex-model is essentially a gm-ib with different algorithms, but identical meaning)
    exchangePositions = true;
  }

  //2.: maybe exchange the order of the two models
  //this means that if we excahnge the models, the keepFirstsNext-parameter can't be regarded anymore!!!
  if (exchangePositions == true) {
    p1 = pSecond;
    p2 = pFirst;
  }

  //3.: check if it is possible to combine the first model with one of the second's type to a new instance of its own type
  switch (p1->Hdr.ModelType) {
    case sclib::mtGMM_new: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM) { //pSecond can't be a GMM-IB/MIX2MAX-model here; if one model is such a model, it is the pFirst
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a GMM with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtMIXMAX: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtMIXMAX && p2->Hdr.ModelType != sclib::mtMIX2MAX && p2->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a GMM-IB with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtMIX2MAX: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtMIXMAX && p2->Hdr.ModelType != sclib::mtMIX2MAX && p2->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a MIX2MAX-Model with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtMIX2MAX_ex: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtMIXMAX && p2->Hdr.ModelType != sclib::mtMIX2MAX && p2->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a MIX2MAXex-Model with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtBGMM: {
      if (p2->Hdr.ModelType != sclib::mtBGMM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Baggenstoss-GMM with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtGMM_UBM: {
      REPORT_ERROR(SVLIB_BadArg, "Cannot combine two GMM-UBMs without reestimation");
      ready = false;
      break;
    }

    case sclib::mtPareto: {
      if (p2->Hdr.ModelType != sclib::mtPareto) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Pareto-Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtSVM: {
      if (p2->Hdr.ModelType != sclib::mtSVM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a one-class SVM Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtHMM: {
      if (p2->Hdr.ModelType != sclib::mtHMM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Hidden Markov Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtFullGauss: {
      if (p2->Hdr.ModelType != sclib::mtFullGauss) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Full-Covariance Gaussian Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtVQ: {
      if (p2->Hdr.ModelType != sclib::mtVQ) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Vector-Quantisation Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtTime: {
      if (p2->Hdr.ModelType != sclib::mtTime) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Time Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtMetaGMM: {
      if (p2->Hdr.ModelType != sclib::mtMetaGMM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a MetaGMM Model with any other model");
        ready = false;
      }
      break;
    }

    default: {
      REPORT_ERROR(SVLIB_BadArg, "Unknown modeltype or combining these models is not posible");
      ready = false;
      break;
    }
  } //switch

  //4.: create a combined model of the fist's type
  if (ready == true) {
    pModel = p1->combineModels(p1, p2, keepFirstsNext);
  }

  return pModel;
}

//====================================================================================================================
//	Combine 2 models by retraining a new one on the concatenated training data of both parants
//  In contrats to the methods in the different models themselfes, this function "knows" which model can be combined
//  with which one
//====================================================================================================================
SC_Model* SC_ModelHandler::combineModels(SC_Model* pFirst, SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  SC_Model *pModel = NULL, *p1 = pFirst, *p2 = pSecond;
  bool exchangePositions = false, ready = true;
  
  //1.: decide which one of the two models should determine the typeof the combined model
  //if the two model have different type, it is important to ecide which one should be the first, 
  //because it determines the type of the combined model, and ome modeltypes incorporated more knowledge/information
  //about the modeled process than others. this first steps therefore decides which model should be the first one
  //according to the "preserve as much information as possible"-criterion
  if (pSecond->Hdr.ModelType == sclib::mtGMM_UBM && pFirst->Hdr.ModelType != sclib::mtGMM_UBM) {
    //a bgmm has "better" training-algorithms so if we have one in the line, use this as the first one
    exchangePositions = true;
  } else if (pSecond->Hdr.ModelType == sclib::mtBGMM && pFirst->Hdr.ModelType == sclib::mtGMM_new) { 
    //a bgmm has "better" training-algorithms so if we have one in the line, use this as the first one
    //except that the first one incorporates more knowledge (as a gmmib/wesely-model does)
    exchangePositions = true;
  } else if (pSecond->Hdr.ModelType == sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
    //a gmm-ib incorporates more knowlege as a standard mixture-model (gmm/bgmm), 
    //so it should be the first one deciding the resultant model's type
    exchangePositions = true;
  } else if (pSecond->Hdr.ModelType == sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
    //the same as gmm-ib (a wesley-model is essentially a gmm-ib with different algorithms, but identical meaning)
    exchangePositions = true;
  } else if (pSecond->Hdr.ModelType == sclib::mtMIX2MAX_ex && pFirst->Hdr.ModelType != sclib::mtMIXMAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX && pFirst->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
    //the same as gmm-ib (a wesleyex-model is essentially a gmm-ib with different algorithms, but identical meaning)
    exchangePositions = true;
  }

  //2.: maybe exchange the order of the two models
  //this means that if we excahnge the models, the keepFirstsNext-parameter can't be regarded anymore!!!
  if (exchangePositions == true) {
    p1 = pSecond;
    p2 = pFirst;
  }

  //3.: check if it is possible to combine the first model with one of the second's type to a new instance of its own type
  switch (p1->Hdr.ModelType) {
    case sclib::mtGMM_new: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM) { //pSecond can't be a GMM-IB/MIX2MAX-model here; if one model is such a model, it is the pFirst
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a GMM with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtGMM_UBM: {
      if (p2->Hdr.ModelType != sclib::mtGMM_UBM && p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM) { //pSecond can't be a GMM-IB/MIX2MAX-model here; if one model is such a model, it is the pFirst
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a GMM-UBM with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtMIXMAX: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtMIXMAX && p2->Hdr.ModelType != sclib::mtMIX2MAX && p2->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a GMM-IB with a non-mixture-model");
        ready = false;
      }
      if (pBackgroundModels->Hdr.ModelType != sclib::mtBGMM && pBackgroundModels->Hdr.ModelType != sclib::mtGMM_new && pBackgroundModels->Hdr.ModelType != sclib::mtMIXMAX && pBackgroundModels->Hdr.ModelType != sclib::mtMIX2MAX && pBackgroundModels->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model is not a mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtMIX2MAX: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtMIXMAX && p2->Hdr.ModelType != sclib::mtMIX2MAX && p2->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a MIX2MAX-Model with a non-mixture-model");
        ready = false;
      }
      if (pBackgroundModels->Hdr.ModelType != sclib::mtBGMM && pBackgroundModels->Hdr.ModelType != sclib::mtGMM_new && pBackgroundModels->Hdr.ModelType != sclib::mtMIXMAX && pBackgroundModels->Hdr.ModelType != sclib::mtMIX2MAX && pBackgroundModels->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model is not a mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtMIX2MAX_ex: {
      if (p2->Hdr.ModelType != sclib::mtGMM_new && p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtMIXMAX && p2->Hdr.ModelType != sclib::mtMIX2MAX && p2->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a MIX2MAXex-Model with a non-mixture-model");
        ready = false;
      }
      if (pBackgroundModels->Hdr.ModelType != sclib::mtBGMM && pBackgroundModels->Hdr.ModelType != sclib::mtGMM_new && pBackgroundModels->Hdr.ModelType != sclib::mtMIXMAX && pBackgroundModels->Hdr.ModelType != sclib::mtMIX2MAX && pBackgroundModels->Hdr.ModelType != sclib::mtMIX2MAX_ex) {
        REPORT_ERROR(SVLIB_BadArg, "Background-model is not a mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtBGMM: {
      if (p2->Hdr.ModelType != sclib::mtBGMM && p2->Hdr.ModelType != sclib::mtGMM_new) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Baggenstoss-GMM with a non-mixture-model");
        ready = false;
      }
      break;
    }

    case sclib::mtPareto: {
      if (p2->Hdr.ModelType != sclib::mtPareto) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Pareto-Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtSVM: {
      if (p2->Hdr.ModelType != sclib::mtSVM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a one-class SVM Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtHMM: {
      if (p2->Hdr.ModelType != sclib::mtHMM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Hidden Markov Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtFullGauss: {
      if (p2->Hdr.ModelType != sclib::mtFullGauss) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Full-Covariance Gaussian Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtVQ: {
      if (p2->Hdr.ModelType != sclib::mtVQ) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Vector-Quantisation Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtTime: {
      if (p2->Hdr.ModelType != sclib::mtTime) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a Time Model with any other model");
        ready = false;
      }
      break;
    }

    case sclib::mtMetaGMM: {
      if (p2->Hdr.ModelType != sclib::mtMetaGMM) {
        REPORT_ERROR(SVLIB_BadArg, "Cannot combine a MetaGMM Model with any other model");
        ready = false;
      }
      break;
    }
		
		default: {
      REPORT_ERROR(SVLIB_BadArg, "Unknown modeltype or combining these models is not posible");
      ready = false;
      break;
    }
  } //switch

  //4.: create a combined model of the fist's type
  if (ready == true) {
    pModel = p1->combineModels(p2, pSpeechFrames, segmentsToMerge, pBackgroundModels);
  }

  return pModel;
}

//====================================================================================================================
// For mixture models: Guess the optimal model-order depending on the guessMode in the tweakable parameters:
//  - Do a search over the space of available modelorders, building a rough model (only few EM steps)  
//    for each and returning the order which maximized BIC criterion
//  ...or...
//  - Following a heuristic, which maps the length of the segment to the # of gaussians
//====================================================================================================================
unsigned short int SC_ModelHandler::guessModelOrder(SV_Data* pFeatures, SC_Model* pBackground, unsigned long int modelType, unsigned short int minOrder, unsigned short int maxOrder, unsigned long int segmentsToMerge) {
  unsigned long int segmentLength, samplesPerGaussian, featureCount = pFeatures->Row, count = 0;
  unsigned short int adaptedMaxOrder, currentOrder, bestOrder = 0;
  double **scores; //order|likelihood|BIC
  double weightLimit;
  SC_MixtureModel *pCandidate = NULL;
  SV_Data *pRes = NULL;
  SC_DistanceMeasures *pDist = NULL;
	SC_Conversion converter(pFeatures->Hdr.sampleRate);

  //count actual number of used features
  if (segmentsToMerge != 1) {
    pRes = pFeatures->Next;
    while (pRes != NULL && (count < segmentsToMerge || segmentsToMerge == 0)) {
      featureCount += pRes->Row;
      pRes = pRes->Next;
      count++;
    }
    pRes = NULL;
  }

  if (isMixtureModel(modelType) == true && modelType != sclib::mtGMM_UBM) {//only mixture models != GMM-UBM have tweakable orders
    if (this->pTweak->modelHandler.orderGuessMode == sclib::guessParameterSearch) { //conduct a search for the best model-size

      MArray_2D(scores, maxOrder-minOrder+1, 3, double, "SC_ModelHandler.guessModelOrder: scores");
      pDist = new SC_DistanceMeasures(this->pTweak);
      
      if (this->verbose == true) {printf("\nMixtureModel-Parameter-Search:\n");}

      for (currentOrder = minOrder; currentOrder <= maxOrder; currentOrder++) {
        pCandidate = (SC_MixtureModel*)createRawModel(modelType, pBackground, currentOrder, pFeatures->Col);
        
        weightLimit = pCandidate->getWeightLimit();
        if (weightLimit > 0.0) {
          adaptedMaxOrder = sclib::min(pFeatures->Row, (unsigned short)floor(1.0 / weightLimit)); //maximum number of gaussians, so that each gaussian models (under the assumption that each gaussian models an equal count of feature vectors) as much feature vectors that the weightlimit (=percent feature vectors to be modeled by one component) isn't violated
        } else {
          adaptedMaxOrder = sclib::min(pFeatures->Row, maxOrder); //avoid div-by-0
        }
        if (currentOrder <= adaptedMaxOrder) { //only build models if there is the possibility that the weightlimit can hold
          pCandidate->setMaxEMsteps(this->pTweak->modelHandler.orderGuessEMsteps); //restrict number of EM steps to save time, because the general trend (low/high likelihood) can be seen very early
          pCandidate->TrainModel(pFeatures, segmentsToMerge);
          pRes = pCandidate->TestModel(pFeatures, segmentsToMerge);
					if (currentOrder != pCandidate->getMixtureCount()) { //happens if there is too less data to create that many mixtures
						if (this->verbose == true) {printf(" Order: %d\timpossible to create, skipped this and remaining orders!", currentOrder);}
						maxOrder = currentOrder-1; //only regard smaller mixtureCounts in finding the best size
						break;
					}
          scores[currentOrder-minOrder][0] = (double)(currentOrder);
          scores[currentOrder-minOrder][1] = (double)(pRes->Mat[0][0]);
          scores[currentOrder-minOrder][2] = pDist->BIC(pCandidate, scores[currentOrder-minOrder][1]*featureCount, featureCount); //the logLikelihood is divided by T to compensate for differing feature-length; because the BIC does do this also, re-multiply it here!
          if (this->verbose == true) {printf(" Order: %d\tLikelihood: %f\tBIC: %f\n", currentOrder, scores[currentOrder-minOrder][1], scores[currentOrder-minOrder][2]);}
        } else { //otherwise set the values this way that they can't be the minimum
          scores[currentOrder-minOrder][0] = (double)(currentOrder);
          scores[currentOrder-minOrder][1] = -1.0 * numeric_limits<double>::max();
          scores[currentOrder-minOrder][2] = numeric_limits<double>::max();
        }
        
        MFree_0D(pRes);
        MFree_0D(pCandidate);
      }

      //pick order of model with smalles BIC-value
      sclib::quickSort(scores, 0, maxOrder-minOrder, 3, 2); //sort scores ascending by bic-value, so the row with the smallest BIC will be the one with index 0
      bestOrder = (unsigned short int)(scores[0][0]);

      if (this->verbose == true) {printf(" Order %d choosen to be best fitting", bestOrder);}

      MFree_2D(scores);
      MFree_0D(pDist);  

    } else { //use a heuristic that tells how many samples must be available to create one component
  		segmentLength = converter.audioFrame2sample(featureCount, pFeatures->Hdr.frameSize, pFeatures->Hdr.frameStep);
			samplesPerGaussian = converter.ms2sample(this->pTweak->modelHandler.msPerGaussian);
      bestOrder = sclib::min(maxOrder, segmentLength/samplesPerGaussian);
    }
  } else {
    bestOrder = 1; //so that the function resturns "true" for non-mixture models
  }
  
  return bestOrder;
}

//====================================================================================================================
// computes the SNR between given signal- and noise features. The features are meant to contain frame-Energys in 
// col 0 (e.g. EnergyZCR features, or MFCC). If they contain other features, this function can also calculate a SNR 
// over 0 to maxColIdx cols, but maybe the result isn't meaningful then.
//====================================================================================================================
double SC_ModelHandler::getSNR(SV_Data* pSignalEnergy, SV_Data* pNoiseEnergy, unsigned int maxColIdx) {
  double snr, signalEnergy = 0.0, noiseEnergy = 0.0;
  unsigned long int i, j;

  if (pNoiseEnergy != NULL) {
    //sum up the feature-values in the desired cols
    for (i = 0; i < (unsigned long)(sclib::max(pSignalEnergy->Row, pNoiseEnergy->Row)); i++) {
      for (j = 0; j <= maxColIdx; j++) {
        if (i < (unsigned long)(pNoiseEnergy->Row)) {
          noiseEnergy += pNoiseEnergy->Mat[i][j] * pNoiseEnergy->Mat[i][j];
        }
        if (i < (unsigned long)(pSignalEnergy->Row)) {
          signalEnergy += pSignalEnergy->Mat[i][j] * pSignalEnergy->Mat[i][j];
        }
      }
    }

    snr = 10.0 * log10(signalEnergy / noiseEnergy); //the dB-formula according to http://en.wikipedia.org/wiki/Decibel
  } else {
    //there was no noise
    snr = 120.0;
  }

  return snr;
}

//====================================================================================================================
//	Returns true if the given model is a child of SC_MixtureModel
//====================================================================================================================
bool SC_ModelHandler::isMixtureModel(SC_Model* pModel) {
  return isMixtureModel(pModel->Hdr.ModelType);
}

//====================================================================================================================
//	Returns true if the given model-type belongs to a child of SC_MixtureModel
//====================================================================================================================
bool SC_ModelHandler::isMixtureModel(long int modelType) {
  bool res = false;

  if (modelType == sclib::mtBGMM || 
      modelType == sclib::mtMIXMAX || 
      modelType == sclib::mtGMM_new || 
      modelType == sclib::mtGMM_UBM ||
      modelType == sclib::mtMIX2MAX ||
      modelType == sclib::mtMIX2MAX_ex) {
    res = true;
  }

  return res;
}

//====================================================================================================================
//	This methods returns the score for the given data and the given model; the speciality is the last parameter can be
//  used to force a model to compute un-normalized (i.e. no likelihood ratios as in the GMM-UBM) scores without 
//  knowledge of the particular abilities of the model; this way, all special knowlege about the model-internas 
//  remains within this class; as a sugar, the score is retuned as a double instead of an SV_Data object; if 
//  averageScores is false, the result will not be divided by the number of testpatterns as is usually (but not al-
//  ways) done by the models.
//====================================================================================================================
double SC_ModelHandler::testModel(SC_Model* pModel, SV_Data* pData, unsigned long int segmentsToMerge, bool forceNoNormalization, SC_Model* pBackgroundModel, bool averageScores) {
  SV_Data *pScore = NULL;
  double res = 0.0;
	SC_FeatureHandler handler(this->pTweak, this->verbose);
	
  //force no normalization, if wished
	unsigned int oldSm = pModel->getTweak()->mixtureModelGmmubm.scoringMethod;
  if (pModel->Hdr.ModelType == sclib::mtGMM_UBM && forceNoNormalization == true) {
		pModel->getTweak()->mixtureModelGmmubm.scoringMethod = (unsigned int)(sclib::scoringGMM);
  }

  pModel->setBackground(pBackgroundModel);
  pScore = pModel->TestModel(pData, segmentsToMerge);
  if (pScore != NULL) {
    res = pScore->Mat[0][0];
    MFree_0D(pScore);

		if (averageScores == false && pModel->scoreIsAverage() == true) {
			res *= (double)(handler.countFeaturesInList(pData, segmentsToMerge));
		} else if (averageScores == true && pModel->scoreIsAverage() == false) {
			res /= (double)(handler.countFeaturesInList(pData, segmentsToMerge));
		}
  } else {
    REPORT_ERROR(SVLIB_Fail, "Failure during model-test: No result was returned");
  }

  //restore previous state
  if (pModel->Hdr.ModelType == sclib::mtGMM_UBM && forceNoNormalization == true) {
		pModel->getTweak()->mixtureModelGmmubm.scoringMethod = oldSm;
  }

  return res;
}

//====================================================================================================================
//	Returns true if the given model type returns averaged scores (i.e. scores divided by the number of test patterns)
//====================================================================================================================
bool SC_ModelHandler::averagesItsScores(long int modelType) {
	SC_Model *pModel = createRawModel(modelType);
	bool res = pModel->scoreIsAverage();

	MFree_0D(pModel);

	return res;
}
