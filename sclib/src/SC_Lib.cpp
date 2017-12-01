/**************************************************************************/
/*	This is the general include file for the SC_Lib library, which should */
/*  be included by all it's users instead of the individual SC_*.h files  */
/*                                                                        */
/*  It also manages (in the correspondig .cpp file) dll-management for    */
/*  windows and provides the "main" function for calls via JNI (probably  */
/*  only important for videana)																						*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 05.03.2004																								*/
/**************************************************************************/

void libConstruct(void);
void libDestruct(void);

//====================================================================================================================
//  Windows DLL (and linux .so) management stuff...
//====================================================================================================================
#ifdef _WIN32
	#include <windows.h>
	#include "SC_Aux.h"
	#include <SV_Error.h>

	BOOL APIENTRY DllMain(HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved) {
		switch (ul_reason_for_call) {
			case DLL_PROCESS_ATTACH:
			case DLL_THREAD_ATTACH:
    		libConstruct();
				break;

			case DLL_THREAD_DETACH:
			case DLL_PROCESS_DETACH:
				libDestruct();
  			break;
		}
	  
		return TRUE;
	}
#else
	void __attribute__ ((constructor)) libConstruct(void);
	void __attribute__ ((destructor)) libDestruct(void);
#endif //_WIN32

//====================================================================================================================
//  Includes - must be preceeded by DllMain()!
//====================================================================================================================
#include "SC_Lib.h"
#ifdef SC_USE_JNI
	#include <jni_SC_Wrapper.h>
#endif
#ifdef SC_USE_MATLAB
	#include "SC_Matlab.h"
#endif

//====================================================================================================================
//  Library constructor/destructor: sets error handler and destructs static members of classes
//====================================================================================================================
void libConstruct(void) {
	//install new error handler
	sclib::errorHandlerThrows(false, false);
	SV_SetErrorHandler(sclib::errorHandler);
	
	//initialize random seeds to have repeatable results
	sclib::getRandomizer()->srandom(sclib::randomSeed);
	srand(sclib::randomSeed);
}

void libDestruct(void) {
	//destructors of global/static objects
	SC_MixtureModel::gaussSolver.~SC_Gauss();
	SC_MixtureModel_GMM_UBM::topMixtureCache.~SC_MixtureCache();

//#ifdef SC_USE_MATLAB
//	SC_Matlab::terminate();
//#endif
}

//====================================================================================================================
//  Global variable (grrr) to hold the progress (in percent) of the current stage (segmentation or speaker 
//  recognition) of the audioSegmentation() function; to be returned by getProgress()
//====================================================================================================================
double progress;

//====================================================================================================================
//  Global variable (grrr) to hold the name of the stage of the audioSegmentation() function; to be returned by 
//  getProgressStage()
//====================================================================================================================
char progressStage[sclib::bufferSize];

//====================================================================================================================
//  Checked (and returned) by doProceed() and doFinish(); if true, the audioSegmentation() function will immeditely 
//  return (can be set on caller's side, e.g. in the VideoAnalyzer)
//====================================================================================================================
bool immediateReturn;

//====================================================================================================================
//  Increases currentProgress by one and sets the progress member
//====================================================================================================================
bool doProceed(long int &currentProgress, long int maxProgress, unsigned int count = 1) {
	if (currentProgress < maxProgress) {
		currentProgress += count;
		progress = (double)(currentProgress) / (double)(maxProgress);
	}

	//printf(" %f%%", progress*100.0);

	return immediateReturn;
}

//====================================================================================================================
//  Increases currentProgress by one and sets the progress member
//====================================================================================================================
bool doFinish(void) {
	progress = 1.0;

	return immediateReturn;
}

//====================================================================================================================
//  Called from outside to tell the audioSegmentation() function to immediately return
//====================================================================================================================
void doAbort(void) {
	immediateReturn = true;

	return;
}

//====================================================================================================================
//  Called by audioSegmentation() function to know if immediate return is wished
//====================================================================================================================
bool isAborted(void) {
	return immediateReturn;
}

//====================================================================================================================
//  Returns the value of the "progress" global variable, the progress (in percent) of the current stage (quite coarse: 
//  segmentation or speaker-id) of the audioSegmentation() function
//====================================================================================================================
double getProgress(void) {
	return progress;
}

//====================================================================================================================
//  Returns the value of the "progressStage" global variable, the name of the stage of the audioSegmentation() 
//  function; memory for this variable is allocated and freed inside the library, so don't touch it!
//====================================================================================================================
char* getProcessingStage(void) {
	return progressStage;
}

//====================================================================================================================
//  "Main" function to be called by Videana via JNI (does in principle what SCiVo does...)
//  Returned are lists with the following information per video-frame:
//    - segmentation-results with bitflags giving the labels belonging to each frame (as the return value)
//    - speaker-ids (as parameter)
//    - probabilities of segmentation results in columns 0-31 and of speaker-id in column 32 (as parameter), if not 
//      switched off in the tweakable parameters
//    - probabilityCols is the number of columns in the probability-matrix
//    - length of each of the three lists (=nr. of video-frames -1)
//  the shotList contains videoFrames-numbers at which shot transition occur
//====================================================================================================================
unsigned long int audioSegmentation(const char *fileName, SC_TweakableParameters *pTweak, double videoFrameRate, long int* segmentationResults, long int* speakerIDs, long int* probabilities, int* probabilityCols, int *shotList, int shotListLength) {
  unsigned long int sceneNr = 1, resultlength, enhancedSignalLength, speakerModelFeature;
  long int sceneStart, sceneEnd, res;
	char previousResults[sclib::bufferSize];
	long int *resultSegmentation = NULL, *resultSpeakerIDs = NULL;
	double **resultProbabilties = NULL;
	long int progressCount, progressMax;
	short *enhancedSignalPointer;
  SC_FeatureHandler *pExtractor = new SC_FeatureHandler(pTweak, false);
  SC_ModelHandler *pModeller = new SC_ModelHandler(pTweak);
  SC_SegmentationHandler *pSegmenter = new SC_SegmentationHandler(pTweak, false);
  SC_SignalHandler *pSigHandler = new SC_SignalHandler(pTweak);
  SC_SpeakerIdentificator *pIdentificator = NULL;
  SC_SpeakerClusterer *pClusterer = NULL;
  SC_Model *pUBM = NULL;      
	SC_Corpus *pCorpus = NULL;
  SC_Signal *pSignal = NULL;
  SV_Data **pFeatures = NULL, *pPitch = NULL, *pConcat;//, *pNorm;
	SC_Cluster *pAllClusters = NULL, *pClusters	= NULL, *pAllShortClusters = NULL, *pShortClusters = NULL, *pHook = NULL;
	SC_SpeakerScore_Clustering *pScore = NULL;
	SC_SpeakerScore_Identification *pScoreID = NULL;
  SC_Enhancement* pEnhancer = NULL;
	bool featuresOk, readAllScenes = true;
	bool segmentationAvailable = false;

	immediateReturn = false;
	
	setbuf(stdout, NULL); //get instant printf()'s
  setbuf(stderr, NULL); //get instant printf()'s

  //decide which corpus to use:
	pCorpus = new SC_Corpus_MPEG7(pTweak, videoFrameRate, fileName, shotList, shotListLength);
	if (pTweak->modelHandler.onlyThisSpeaker != NULL && strlen(pTweak->modelHandler.onlyThisSpeaker) > 0) {
    pTweak->setDebugPrefix(pTweak->modelHandler.onlyThisSpeaker);
  }
	pCorpus->getGT()->checkConsistency(0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::modeGroundtruth);

	sprintf(progressStage, "Audio segmentation\0");

	//try to load previous segmentation results
	if (pTweak->general.preClusteringResultsPrefix != NULL && strncmp(pTweak->general.preClusteringResultsPrefix, "", 1) != 0) {
		sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, ".gt");
	}
	if (sclib::fileExists(previousResults) == true) { 
		//progressMax = 3;
		//progressCount = 0;

		if (pCorpus->getGT()->load(previousResults) == NULL) {
			REPORT_ERROR(SVLIB_FileErr, "Error loading ground truth class!!!");
		}
		segmentationAvailable = true;
		//immediateReturn = doProceed(progressCount, progressMax); //TODO: handle "abort"

		sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, "_speakers.cluster_1");
		if (sclib::fileExists(previousResults)) {
			pAllClusters = new SC_Cluster(pTweak);
			sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, "_speakers");
			if (pAllClusters->load(previousResults) == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "Error loading speaker clusters!!!");
			}
		}
		//immediateReturn = doProceed(progressCount, progressMax);

		sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, "_short.cluster_1");
		if (sclib::fileExists(previousResults)) {
			pShortClusters = new SC_Cluster(pTweak);
			sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, "_short");
			if (pShortClusters->load(previousResults) == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "Error loading short clusters!!!");
			}
		}
		//immediateReturn = doProceed(progressCount, progressMax);
		if (pAllClusters!=NULL && segmentationAvailable==true) {
			progressMax = 3;
			progressCount = 0;
			immediateReturn = doProceed(progressCount, progressMax, progressMax);
		}
	} 

	if (pAllClusters == NULL) { //calculate new segmentation results if there are no clusters loaded (then, presumably, only segmentation has been carried out and this is the subsequent clustering run
		//TODO: make feature selection more sophisticated
		speakerModelFeature = (pTweak->modelHandler.foregroundModelType == sclib::mtMIXMAX) ? sclib::featureFbE : pTweak->modelHandler.speakerModelFeature; //MIXMAX models only work on FBE features

		progressMax = pCorpus->getGT()->getRealSceneCount() * 10;
		progressCount = 0;
		if (pTweak->enhancement.doEnhancement == true) {
			pEnhancer = new SC_Enhancement(pCorpus->getGT(), pTweak, pTweak->enhancement.speechModelFile);
		}

		//processing scene-wise...
		for (unsigned long int y = 0; y < pCorpus->getGT()->getAudioSampleCount(); y++) {
	    
			pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
			if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
			
				//skip scene?
				if ((sceneNr < pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & pTweak->general.sceneSelection) && (pTweak->general.sceneSelection != 0))) {
					sceneNr++; 
					y = sceneEnd; 
					continue;
				}

				//load signal, extract features (2 proceeds)
				//------------------------------------------
				//load previous work if wished and existant
				if (pTweak->general.featurePrefix!=NULL && strncmp(pTweak->general.featurePrefix, "", 1)!=0) {
					sprintf(previousResults, "%s_%d%s\0", pTweak->general.featurePrefix, sceneNr, ".dat");
				} else {
					sprintf(previousResults, "");
				}
				if (sclib::fileExists(previousResults) == true) {
					pFeatures = pExtractor->loadFeatures(previousResults);
					featuresOk = pExtractor->checkFeatureSet(pFeatures, pSegmenter->getUsedFeatures(sclib::algorithm_nothing, 2, -1), NULL); 
					immediateReturn = doProceed(progressCount, progressMax, 2);
				} else {
 					pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd); //load signal
					immediateReturn = doProceed(progressCount, progressMax);
					pExtractor->setTweak(pTweak); //extract features for standard algorithms
					if (pTweak->enhancement.doEnhancement == true) { //in case of speech enhancement to be done, these features are only needed for non-speech-specific segmentation and are extracted once more on the enhanced signal later on
						pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(sclib::algorithm_nothing, 2, -1));
						featuresOk = pExtractor->checkFeatureSet(pFeatures, pSegmenter->getUsedFeatures(sclib::algorithm_nothing, 2, -1), NULL); 
					} else { //if no enhancement is to be carried out, we don't need the signal any more
						pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures()|speakerModelFeature, pSegmenter->getSpecialFeatureParameters());
						featuresOk = pExtractor->checkFeatureSet(pFeatures, pSegmenter->getUsedFeatures()|speakerModelFeature, pSegmenter->getSpecialFeatureParameters());
						MFree_0D(pSignal);
					}

					//save previous work for future runs if wished
					if (pTweak->enhancement.doEnhancement==false && pTweak->general.featurePrefix!=NULL && strncmp(pTweak->general.featurePrefix, "", 1)!=0) {
						sprintf(previousResults, "%s_%d%s\0",pTweak->general.featurePrefix, sceneNr, ".dat");
						pExtractor->saveFeatures(previousResults, pFeatures);
					}
					immediateReturn = doProceed(progressCount, progressMax);
				}

				if (featuresOk == true) {
					//non-speech-specific segmentation (3 proceeds)
					//---------------------------------------------
					if (segmentationAvailable == false) {
						pSegmenter->silenceDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, pTweak->segmentationHandler.silenceDetectorMode); //silence detection
						immediateReturn = doProceed(progressCount, progressMax);
						pSegmenter->audioClassification(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, pTweak->segmentationHandler.audioTypeMode); //audio type classification
						immediateReturn = doProceed(progressCount, progressMax);
						res = pSegmenter->acousticChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, pTweak->segmentationHandler.changeDetectorMode); //acd
						immediateReturn = doProceed(progressCount, progressMax);
					} else {
						immediateReturn = doProceed(progressCount, progressMax, 3);
					}

					//speech enhancement (2 proceeds)
					//------------------------------
					if (pTweak->enhancement.doEnhancement == true) {
						pEnhancer->setSegmentStart(sceneStart);
						pEnhancer->setOriginalSignalPointer(pSignal->GetBuf_L(), pSignal->GetLen());
						if (pEnhancer->enhance(false) == true) {
							pEnhancer->getEnhancedSignalPointer(enhancedSignalPointer, enhancedSignalLength);
							pEnhancer->forgetEnhancedSignalPointer();
							pSignal->setBuf_L(enhancedSignalPointer, enhancedSignalLength);
							immediateReturn = doProceed(progressCount, progressMax);
							pExtractor->freeFeatures(pFeatures); //re-extract features for standard algorithms on the enhanced signal for the now following speech-specific algorithms
							pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(sclib::algorithm_nothing, -1, 3)|speakerModelFeature, pSegmenter->getSpecialFeatureParameters());
							immediateReturn = doProceed(progressCount, progressMax);
						} else {
							immediateReturn = doProceed(progressCount, progressMax, 2);
							REPORT_ERROR(0, "Unknown Error during enhancement!");
						}
						MFree_0D(pSignal);
					} else {
						immediateReturn = doProceed(progressCount, progressMax, 2);
					}

					//speech-specific segmentation (2 proceeds)
					//-----------------------------------------
					if (segmentationAvailable == false) {
						pPitch = pSegmenter->unvoicedDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, pTweak->segmentationHandler.vUvDetectorMode); //label unvoiced speech
						if (pTweak->featurePitch.method == sclib::modeAS && (sclib::bitTest(speakerModelFeature, sclib::featurePitch) == true || sclib::bitTest(pSegmenter->getUsedFeatures(), sclib::featurePitch) == true)) {
							sclib::destructLinkedList(pFeatures[sclib::bitPosition(sclib::featurePitch)]);
							double invalidValue[1] = {0.0};
							SV_Data *pConvertedPitch = pExtractor->convertFrameRate(pPitch, sceneEnd-sceneStart+1, pCorpus->getGT()->getConverter()->ms2sample(pTweak->featurePitch.frameSize), pCorpus->getGT()->getConverter()->ms2sample(pTweak->featurePitch.frameStep), invalidValue);
							pFeatures[sclib::bitPosition(sclib::featurePitch)] = pConvertedPitch;
							MFree_0D(pPitch);
						} else {
							MFree_0D(pPitch);
						}
						immediateReturn = doProceed(progressCount, progressMax);
						res = pSegmenter->speakerChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, pTweak->segmentationHandler.changeDetectorMode); //scd
						immediateReturn = doProceed(progressCount, progressMax);
					} else {
						immediateReturn = doProceed(progressCount, progressMax, 2);
					}

					//build initial clusters with speaker models (1 proceed)
					//------------------------------------------------------
					pExtractor->prepareFeatureSet(pFeatures, NULL); //bring features for standrad parameters (i.e. for speaker clustering) "up" in the linked list
					pExtractor->equalizeFrameParameters(pFeatures, pExtractor->getFeatureCount(), speakerModelFeature);
					pConcat = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount(), speakerModelFeature);
					pExtractor->freeFeatures(pFeatures); //free features
					//if (pCorpus->getGT()->getGTtype() == sclib::gtTIMIT && pTweak->modelHandler.foregroundModelType == sclib::mtGroup) {
					//	((SC_GroundTruth_TIMIT*)pCorpus->getGT())->addPhoneLabels(pConcat, sceneStart, sceneEnd, sclib::modeHypothesized);
					//}
					//if (sclib::isPowerOfTwo(speakerModelFeature) == false) {
					//	pNorm = pExtractor->createNormalizationMatrix(pConcat);
					//	pExtractor->normalize(pConcat, pNorm);
					//	MFree_0D(pNorm);
					//}
					if (pTweak->speakerClusterer.doClustering == true || pTweak->speakerIdentification.doIdentification == true) {
						res = pModeller->buildSpeakerModels(pCorpus, sceneStart, sceneEnd, pConcat, pClusters, pShortClusters);
					}
					//print the speaker model features (regardless of if modeling was done) to files
					if (sceneNr > pTweak->general.firstScene) {
						sclib::stringOut("audio.llf", "", pTweak); //separate the features of different shots by a blank line
					}
					((SC_Corpus_MPEG7*)pCorpus)->featureOut("audio.llf", sceneStart, sceneEnd, pConcat, (sceneNr == pTweak->general.firstScene), false, NULL, true);
					MFree_0D(pConcat);
					immediateReturn = doProceed(progressCount, progressMax);
		      
					//link normal clusters together
					pHook = sclib::getLastInList(pAllClusters);
					if (pHook != NULL) {
						pHook->Next	= pClusters;
					} else {
						pAllClusters = pClusters;
					}

					//link short "clusters" together (they are not really clustes because they lack a speaker model [in fact the segments therein where too short to build one] but just hold the features and segment-borders)
					pHook = sclib::getLastInList(pAllShortClusters);
					if (pHook != NULL) {
						pHook->Next	= pShortClusters;
					} else {
						pAllShortClusters = pShortClusters;
					}
				} else { //features not ok => shot/scene typically to short, so write a dummy llf file entry
					if (sceneNr > pTweak->general.firstScene) {
						sclib::stringOut("audio.llf", "", pTweak); //separate the features of different shots by a blank line
					}
					sclib::stringOut("audio.llf", "0", pTweak); //a line with only one "0" tells that this shot has no features!
					immediateReturn = doProceed(progressCount, progressMax, 8);
				}

				sceneNr++;
				y = sceneEnd;

				pCorpus->getGT()->checkConsistency(sceneStart, sceneEnd, sclib::modeHypothesized);

				if (sceneNr > pTweak->general.lastScene || immediateReturn == true) {
					readAllScenes = false;
					break;
				}
			} else {
				break;	//scene not valid
			} 
		}

		//care for the fact that the audio decoder most ofton report less frames than the video decoder and there might be shots in the cutlist after the end of the file as we know it that also need features!
		if (pCorpus->getGT()->getGTtype() == sclib::gtMPEG7) {
			if (readAllScenes == true) { //only consider those shots if the end of the file should be evluated
				for (int i = 0; i < ((SC_GroundTruth_MPEG7*)pCorpus->getGT())->getPostEndCutCount(); i++) {
					sclib::stringOut("audio.llf", "", pTweak); //separate the features of different shots by a blank line
					sclib::stringOut("audio.llf", "0", pTweak); //a line with only one "0" tells that this shot has no features!
				}				
			}
		}

		//try saving segmentation results
		immediateReturn = doFinish();
		if (immediateReturn == false) {
			if (pTweak->general.preClusteringResultsPrefix != NULL && strncmp(pTweak->general.preClusteringResultsPrefix, "", 1) != 0) {
				sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, ".gt");
				if (pCorpus->getGT()->save(previousResults) != true) {
					REPORT_ERROR(SVLIB_FileErr, "Error saving ground truth class!");
				}
				if (pAllClusters != NULL) {
					sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, "_speakers");
					if (pAllClusters != NULL && pAllClusters->save(previousResults) != true) {
						REPORT_ERROR(SVLIB_FileErr, "Error saving speaker clusters!");
					}
				}
				if (pShortClusters != NULL) {
					sprintf(previousResults, "%s%s\0", pTweak->general.preClusteringResultsPrefix, "_short");
					if (pShortClusters != NULL && pShortClusters->save(previousResults) != true) {
						REPORT_ERROR(SVLIB_FileErr, "Error saving short clusters!");
					}
				}
			}
		}
	}

	sprintf(progressStage, "Speaker recognition\0");
	progressMax = sclib::getListCount(pAllClusters) + sclib::getListCount(pShortClusters);
	progressCount = 0;

  if (pAllClusters != NULL && pTweak->speakerClusterer.doClustering == true && immediateReturn == false) {
    //Clustering Speakers...
		//after this, pAllClusters is linked in the first partition of pScore->pPartitionList; it get's destructed when pScore gets destructed
		//pShortClusters is again a linked list of unhandled (due to short duration) speeech segments
		pClusterer = new SC_SpeakerClusterer(pTweak);
	  pScore = pClusterer->clusterSpeakers(pCorpus->getGT(), pAllClusters, pShortClusters, progressCount, progressMax, &doProceed); 
		if (pAllShortClusters == NULL) {
			pAllShortClusters = pShortClusters;
		} else {
			if (pShortClusters != NULL) { //set the new short clusters on top of the complete list because it may contain the bigger segments
				pHook = sclib::getLastInList(pShortClusters);
				pHook->Next = pAllShortClusters;
				pAllShortClusters = pShortClusters;
			}
		}
	  MFree_0D(pClusterer);
		immediateReturn = isAborted();

    //Identifying Speakers of small segments
		if (pAllShortClusters != NULL && pTweak->speakerIdentification.doIdentification == true && immediateReturn == false) { //try speaker-identification for those segments too short to build models and/or clusters
      if (pScore->getFinalPartition() != NULL) {
				pUBM = pModeller->loadModel(pTweak->mixtureModelGmmubm.ubmFileName, sclib::mtGMM_new);
				pIdentificator = new SC_SpeakerIdentificator(pTweak, (SC_MixtureModel*)pUBM, true);
				pScoreID = pIdentificator->identifySpeakers(pCorpus->getGT(), pAllShortClusters, pScore->getFinalPartition()->getClusters(), 0.0, true, false, progressCount, progressMax, &doProceed); 
				MFree_0D(pUBM);
				MFree_0D(pIdentificator);
				immediateReturn = isAborted();
				if (immediateReturn == false) {
					pTweak->speakerClusterer.speechSegLengthThreshold = pTweak->general.shortSpeechThreshold; //now include previously too short segments into clustering-scoring
					pScore->calcScores(); //the partition was altered with the result of the identification, so recalc the cluster-scores!
				}
      }
    }
	}
	
	if (immediateReturn == false) {
		resultlength = ((SC_GroundTruth_MPEG7*)pCorpus->getGT())->getResults(resultSegmentation, resultSpeakerIDs, resultProbabilties);
		*probabilityCols = (int)(pCorpus->getGT()->getProbabilityListDim());
		*segmentationResults = (long int)(&resultSegmentation[0]);
		*speakerIDs = (long int)(&resultSpeakerIDs[0]);
		*probabilities = (long int)(&resultProbabilties[0]);
		immediateReturn = doFinish();
	} else {
		resultlength = 0;
		*probabilityCols = 0;
		*segmentationResults = 0;
		*speakerIDs = 0;
		*probabilities = 0;
	}
	
	//Cleaning up...
  sclib::destructLinkedList(pAllShortClusters);
  MFree_0D(pScoreID);
  MFree_0D(pAllShortClusters);
  MFree_0D(pScore); //this releases all the partitions created by the clustering and, within, all the clusters, too
  MFree_0D(pCorpus);
  MFree_0D(pModeller);
  MFree_0D(pExtractor);
  MFree_0D(pSegmenter);
  MFree_0D(pSigHandler);
	MFree_0D(pEnhancer);

	return resultlength;
}

//====================================================================================================================
//  some wrappers to be called from java via JNI
//====================================================================================================================
#ifdef SC_USE_JNI
//====================================================================================================================
//  For mathew's GSM phone project: gets one frame of samples, returns its MFCC representation according to standard
//  parameters
//====================================================================================================================
#include "SC_Feature_MFCC.h"
JNIEXPORT jdoubleArray JNICALL Java_jni_SC_1Wrapper_frame2feature(JNIEnv *env, jobject caller, jshortArray samples, jint sampleRate) {
	sclib::errorHandlerThrows(false, true);
	try	 {
		int sampleCount = env->GetArrayLength(samples);
		SC_Feature_MFCC extractor(sampleRate, sampleCount, sampleCount, 20, (sampleRate>8000)?24:20, sclib::checkFftSize(512, sampleCount), sclib::wndHamming, 0.97, NULL, false, false, false, false, sclib::modeSClib, sclib::scaleMel, sclib::smoothNone, 0.0, (double)(sampleRate)/2.0);
		short *cs = env->GetShortArrayElements(samples, NULL);

		extractor.CopySignal(cs, sampleCount);
		env->ReleaseShortArrayElements(samples, cs, 0);
		SV_Data *pFrame = extractor.ExtractFeature();
		
		jdoubleArray jFrame = env->NewDoubleArray(pFrame->Col);
		for (int d = 0; d < sampleCount; d++) {
			double coeff = pFrame->Mat[0][d];
			env->SetDoubleArrayRegion(jFrame, d, 1, &coeff);
		}
		MFree_0D(pFrame);

		return jFrame;
	} catch (const char *errStr) {
		printf("%s%s", "ERROR: ", errStr);
		return NULL;
	}
}

//====================================================================================================================
//  For mathew's GSM phone project: gets 2 MFCCs, returns their euclidean distance
//====================================================================================================================
JNIEXPORT jdouble JNICALL Java_jni_SC_1Wrapper_frameDistance(JNIEnv *env, jobject caller, jdoubleArray x, jdoubleArray y) {
	sclib::errorHandlerThrows(false, true);
	try	 {
		unsigned short int dim = (unsigned short int)(env->GetArrayLength(x));
		double *cx = env->GetDoubleArrayElements(x, NULL);
		double *cy = env->GetDoubleArrayElements(y, NULL);
		
		jdouble res = SC_DistanceMeasures::euclid(cx, cy, dim);

		env->ReleaseDoubleArrayElements(x, cx, 0);
		env->ReleaseDoubleArrayElements(x, cy, 0);
		
		return res;
	} catch (const char *errStr) {
		printf("%s%s", "ERROR: ", errStr);
		return std::numeric_limits<double>::max();
	}
}

//====================================================================================================================
//  4 wrappers to be called from java via JNI that get the audio samples via a java stream; for webVoice webservice
//====================================================================================================================
JNIEXPORT jstring JNICALL Java_jni_SC_1Wrapper_wav2splice(JNIEnv *env, jobject caller, jobject jStream, jstring iniFile, jstring resultPath, jboolean usePhaseSynthesis, jint preservedBlockSize, jint intermediateFrameCount, jdouble steepness, jdouble olaErrorTarget, jint olaMaxIterationCount, jint spectrumFrameSize, jint spectrumFrameStep, jint spectrumFftSize) {
	char error[sclib::bufferSize];
	sclib::errorHandlerThrows(false, true);

	try	 {
		SC_Corpus *pCorpus;
		SC_Signal *pSignal;
		long int sceneStart, sceneEnd;
		SV_Data **pFeatures, *pTmp;
		const char *iniFileName = env->GetStringUTFChars(iniFile, 0);
		SC_TweakableParameters *pTweak = new SC_TweakableParameters(iniFileName);
		SC_FeatureHandler extractor(pTweak, true);
		SC_Synthesis synthesizer(pTweak, pTweak->featureMfcc.fftSize, 1.0, pTweak->featureMfcc.window);
		char audioFile[sclib::bufferSize];
		const char *resultPathName = env->GetStringUTFChars(resultPath, 0);
		char *path = sclib::makePath(resultPathName);
		unsigned long int featureType = (usePhaseSynthesis != false) ? sclib::featureSpectrum : sclib::featureSamples;
		bool res;

		setbuf(stdout, NULL); //get instant printf()s
		setbuf(stderr, NULL); //get instant printf()s
		std::ios::sync_with_stdio(); //synchronize couts and printf()s

		//release the strings previously copied from java
		env->ReleaseStringUTFChars(iniFile, iniFileName);
		env->ReleaseStringUTFChars(resultPath, resultPathName);

		pCorpus = new SC_Corpus_MPEG7(pTweak, 1.0, env, jStream);

		//set some parameters
		pTweak->featureSpectrum.frameSize = sclib::getBetween(1, spectrumFrameSize, 3000);
		pTweak->featureSpectrum.frameStep = sclib::getBetween(1, spectrumFrameStep, pTweak->featureSpectrum.frameStep);
		pTweak->featureSpectrum.FFTsize = sclib::checkFftSize(spectrumFftSize, pCorpus->getGT()->getConverter()->ms2sample(pTweak->featureSpectrum.frameSize));
		pTweak->featureSpectrum.preEmphasizeFactor = 0.0;
		pTweak->featureSpectrum.lowCut = 0.0;
		pTweak->featureSpectrum.highCut = 0.0;
		pTweak->featureSpectrum.createPhase = false;
		pTweak->featureSpectrum.logarithmize = true;
		pTweak->featureSpectrum.window = sclib::wndHamming;
		pTweak->featureSamples.frameSize = pTweak->featureSpectrum.frameSize;
		pTweak->featureSamples.frameStep = pTweak->featureSpectrum.frameStep;
		pTweak->featureSamples.highCut = pTweak->featureSpectrum.highCut;
		pTweak->featureSamples.lowCut = pTweak->featureSpectrum.lowCut;
		synthesizer.setFftLength(pTweak->featureSpectrum.FFTsize);
		synthesizer.setOlaMaxIterations(olaMaxIterationCount);
		synthesizer.setOlaErrorTarget(olaErrorTarget);
		synthesizer.setTaperingLength(pTweak->transform.taperingLength);
		synthesizer.setTaperingMode(sclib::wndHamming);

		//load the signal
		sceneStart = 0;
		sceneEnd = pCorpus->getGT()->getAudioSampleCount()-1;
		pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, featureType);
		extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list

		//splice & blend & resynthesize
		extractor.splice(pFeatures[sclib::bitPosition(featureType)], preservedBlockSize); //randomize time-order of frames
		pTmp = extractor.blend(pFeatures[sclib::bitPosition(featureType)], intermediateFrameCount, steepness, preservedBlockSize);
		sprintf(audioFile, "%s%s%s", path, tmpnam(NULL)+1, "wav"); //get a new temporary fileName
		res = synthesizer.feature2wav(audioFile, pTmp, NULL, true);
		if (res == false) {
			sprintf(error, "%s%s", "ERROR: no result returned by resynthesis, maybe no features given/extracted?");
		}
		MFree_0D(pTmp);
		MFree_1D(path);

		extractor.freeFeatures(pFeatures);
		MFree_0D(pSignal);
		MFree_0D(pCorpus);

		return (res==true) ? env->NewStringUTF(audioFile) : env->NewStringUTF(error);	
	} catch (const char *errStr) {
		sprintf(error, "%s%s", "ERROR: ", errStr);

		return env->NewStringUTF(error);
	}
}

JNIEXPORT jstring JNICALL Java_jni_SC_1Wrapper_wav2features(JNIEnv *env, jobject caller, jobject jStream, jstring iniFile, jstring resultPath, jint frameSize, jint frameStep, jboolean useMFCC, jboolean addPitch, jint mfccOrder, jint lpcOrder, jdouble preemphasis, jint window, jint fftSize, jint frequencyScale, jint filterbankSize, jdouble minFilterbankFrequency, jdouble maxFilterbankFrequency, jdouble olaErrorTarget, jint olaMaxIterationCount, jint pitchFrameSize, jint pitchFrameStep, jdouble pitchCandidateThresh, jdouble pitchLagWeight, jdouble pitchFrequencyWeight, jdouble pitchTransitionCost, jdouble pitchTransitionAmpModCost, jdouble pitchTransitionSpecModCost, jdouble pitchVoiceBias, jdouble pitchDoublingCost, jdouble pitchMinF0, jdouble pitchMaxF0, jint pitchCandidateCount, jdouble pitchCandidateWindowSize) {
	char error[sclib::bufferSize];
	sclib::errorHandlerThrows(false, true);

	try {
		SC_Corpus *pCorpus;
		SC_Signal *pSignal;
		long int sceneStart, sceneEnd;
		SV_Data **pFeatures, *pFeature;
		const char *iniFileName = env->GetStringUTFChars(iniFile, 0);
		SC_TweakableParameters *pTweak = new SC_TweakableParameters(iniFileName);
		SC_FeatureHandler extractor(pTweak, true);
		SC_Synthesis synthesizer(pTweak, pTweak->featureMfcc.fftSize, 1.0, pTweak->featureMfcc.window);
		char audioFile[sclib::bufferSize];
		const char *resultPathName = env->GetStringUTFChars(resultPath, 0);
		char *path = sclib::makePath(resultPathName);
		bool res;

		setbuf(stdout, NULL); //get instant printf()s
		setbuf(stderr, NULL); //get instant printf()s
		std::ios::sync_with_stdio(); //syncronize couts and printf()s

		//release the strings previously copied from java
		env->ReleaseStringUTFChars(iniFile, iniFileName);
		env->ReleaseStringUTFChars(resultPath, resultPathName);

		pCorpus = new SC_Corpus_MPEG7(pTweak, 1.0, env, jStream);

		//set some parameters
		if (addPitch != 0) {
			pTweak->modelHandler.speakerModelFeature = sclib::featurePitch;
		}
		if (mfccOrder>0 && useMFCC!=0) {
			pTweak->modelHandler.speakerModelFeature |= sclib::featureMFCC;
		} else if (lpcOrder>0 && useMFCC==0) {
			pTweak->modelHandler.speakerModelFeature |= sclib::featureLPC;
		}
		pTweak->featureMfcc.frameSize = sclib::getBetween(1, frameSize, 3000);
		pTweak->featureMfcc.frameStep = sclib::getBetween(1, frameStep, pTweak->featureMfcc.frameSize);
		pTweak->featureMfcc.fftSize = sclib::checkFftSize(fftSize, pCorpus->getGT()->getConverter()->ms2sample(pTweak->featureMfcc.frameSize));
		pTweak->featureMfcc.filterBankSize = sclib::getBetween(1, filterbankSize, pTweak->featureMfcc.filterBankSize/2+1);
		pTweak->featureMfcc.MFCCorder = sclib::getBetween((unsigned short int)(0), mfccOrder, pTweak->featureMfcc.filterBankSize);
		pTweak->featureMfcc.preEmphasizeFactor = sclib::getBetween(0.0, preemphasis, 1.0);
		pTweak->featureMfcc.sclib_minFilterBankFrequency = sclib::getBetween(0.0, minFilterbankFrequency, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0);
		pTweak->featureMfcc.sclib_maxFilterBankFrequency = sclib::getBetween(pTweak->featureMfcc.sclib_minFilterBankFrequency, maxFilterbankFrequency, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0);
		pTweak->featureMfcc.addDeltaDeltas = false; //the folowwing parameters can't be made audible and thus are not changeable by the user
		pTweak->featureMfcc.addDeltas = false;
		pTweak->featureMfcc.coeffSelection[0] = 0xFF;
		pTweak->featureMfcc.coeffSelection[1] = 0xFF;
		pTweak->featureMfcc.coeffSelection[2] = 0xFF;
		pTweak->featureMfcc.coeffSelection[3] = 0xFF;
		pTweak->featureMfcc.coeffSelection[4] = 0xFF;
		pTweak->featureMfcc.coeffSelection[5] = 0xFF;
		pTweak->featureMfcc.coeffSelection[6] = 0xFF;
		pTweak->featureMfcc.coeffSelection[7] = 0xFF;
		pTweak->featureMfcc.CMN = false;
		pTweak->featureMfcc.dEnergy = false;
		pTweak->featureMfcc.highCut = 0.0;
		pTweak->featureMfcc.lowCut = 0.0;
		pTweak->featureMfcc.method = sclib::modeSClib;
		pTweak->featureMfcc.sclib_frequencyScale = sclib::getBetween(0, frequencyScale, 4);
		pTweak->featureMfcc.sclib_smoothing = sclib::smoothNone;
		pTweak->featureMfcc.window = sclib::getBetween(0, window, 3);
		synthesizer.setFftLength(pTweak->featureMfcc.fftSize);
		synthesizer.setOlaMaxIterations(olaMaxIterationCount);
		synthesizer.setOlaErrorTarget(olaErrorTarget);
		synthesizer.setTaperingLength(pTweak->transform.taperingLength);
		synthesizer.setTaperingMode(pTweak->featureMfcc.window);
		pTweak->featureLpc.lowCut = 0.0;
		pTweak->featureLpc.highCut = 0.0;
		pTweak->featureLpc.frameSize = sclib::getBetween(1, frameSize, 3000);;
		pTweak->featureLpc.frameStep = sclib::getBetween(1, frameStep, pTweak->featureLpc.frameSize);
		pTweak->featureLpc.LPCorder = sclib::getBetween((unsigned short int)(0), lpcOrder, (unsigned short int)(pCorpus->getGT()->getConverter()->ms2sample(frameSize)));
		pTweak->featureLpc.preEmphasizeFactor = sclib::getBetween(0.0, preemphasis, 1.0);
		pTweak->featureLpc.window = sclib::getBetween(0, window, 3);
		pTweak->featureLpc.computeGain = false;
		pTweak->featurePitch.frameSize = sclib::getBetween(1, pitchFrameSize, 3000);
		pTweak->featurePitch.frameStep = sclib::getBetween(1, pitchFrameStep, pTweak->featurePitch.frameSize);
		pTweak->featurePitch.method = sclib::modeESPS;
		pTweak->featurePitch.esps_cand_thresh = (float)(sclib::getBetween(0.01, pitchCandidateThresh,	0.99));
		pTweak->featurePitch.esps_lag_weight = (float)(sclib::getBetween(0.0, pitchLagWeight, 1.0));
		pTweak->featurePitch.esps_freq_weight = (float)(sclib::getBetween(0.0, pitchFrequencyWeight, 1.0));
		pTweak->featurePitch.esps_trans_cost = (float)(sclib::getBetween(0.0, pitchTransitionCost, 1.0));
		pTweak->featurePitch.esps_trans_amp = (float)(sclib::getBetween(0.0, pitchTransitionAmpModCost, 100.0));
		pTweak->featurePitch.esps_trans_spec = (float)(sclib::getBetween(0.0, pitchTransitionSpecModCost, 100.0));
		pTweak->featurePitch.esps_voice_bias = (float)(sclib::getBetween(-1.0, pitchVoiceBias, 1.0));
		pTweak->featurePitch.esps_double_cost = (float)(sclib::getBetween(0.0, pitchDoublingCost, 10.0));
		pTweak->featurePitch.esps_min_f0 = (float)(sclib::getBetween(10.0, pitchMinF0, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0));
		pTweak->featurePitch.esps_max_f0 = (float)(sclib::getBetween((double)(pTweak->featurePitch.esps_min_f0), pitchMaxF0, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0));
		pTweak->featurePitch.esps_n_cands = sclib::getBetween(3, pitchCandidateCount, 100);
		pTweak->featurePitch.esps_wind_dur = (float)(sclib::getBetween(10.0/(double)(pCorpus->getGT()->getAudioSampleRate()), pitchCandidateWindowSize, 0.1));

		//load the signal
		sceneStart = 0;
		sceneEnd = pCorpus->getGT()->getAudioSampleCount()-1;
		pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pTweak->modelHandler.speakerModelFeature);
		extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list

		//resynthesize features
		sprintf(audioFile, "%s%s%s", path, tmpnam(NULL)+1, "wav"); //get a new temporary fileName
		if (pFeatures[sclib::bitPosition(sclib::featureMFCC)] != NULL) {
			pFeature = pFeatures[sclib::bitPosition(sclib::featureMFCC)];
		} else if (pFeatures[sclib::bitPosition(sclib::featureLPC)] != NULL) {
			pFeature = pFeatures[sclib::bitPosition(sclib::featureLPC)];
		} else {
			pFeature = pFeatures[sclib::bitPosition(sclib::featurePitch)];
		}
		res = synthesizer.feature2wav(audioFile, pFeature, pFeatures[sclib::bitPosition(sclib::featurePitch)], true);
		if (res == false) {
			sprintf(error, "%s%s", "ERROR: no result returned by resynthesis, maybe no features given/extracted?");
		}
		MFree_1D(path);

		extractor.freeFeatures(pFeatures);
		MFree_0D(pSignal);
		MFree_0D(pCorpus);

		return (res==true) ? env->NewStringUTF(audioFile) : env->NewStringUTF(error);	
	} catch (const char *errStr) {
		sprintf(error, "%s%s", "ERROR: ", errStr);

		return env->NewStringUTF(error);
	}
}

JNIEXPORT jstring JNICALL Java_jni_SC_1Wrapper_wav2gmm(JNIEnv *env, jobject caller, jobject jStream, jstring iniFile, jstring resultPath, jint mixtureCount, jboolean fullCovar, jdouble varianceLimit, jint maxEmIterations, jdouble emThreshold, jint frameSize, jint frameStep, jboolean useMFCC, jboolean addPitch, jint mfccOrder, jint lpcOrder, jdouble preemphasis, jint window, jint fftSize, jint frequencyScale, jint filterbankSize, jdouble minFilterbankFrequency, jdouble maxFilterbankFrequency, jdouble olaErrorTarget, jint olaMaxIterationCount, jint pitchFrameSize, jint pitchFrameStep, jdouble pitchCandidateThresh, jdouble pitchLagWeight, jdouble pitchFrequencyWeight, jdouble pitchTransitionCost, jdouble pitchTransitionAmpModCost, jdouble pitchTransitionSpecModCost, jdouble pitchVoiceBias, jdouble pitchDoublingCost, jdouble pitchMinF0, jdouble pitchMaxF0, jint pitchCandidateCount, jdouble pitchCandidateWindowSize) {
	char error[sclib::bufferSize];
	sclib::errorHandlerThrows(false, true);

	try {
		SC_Corpus *pCorpus;
		SC_Signal *pSignal;
		long int sceneStart, sceneEnd;
		SV_Data **pFeatures, *pConcat, *pTmp, *pFeature, *pPitch;
		const char *iniFileName = env->GetStringUTFChars(iniFile, 0);
		SC_TweakableParameters *pTweak = new SC_TweakableParameters(iniFileName);
		SC_FeatureHandler extractor(pTweak, true);
		SC_Synthesis synthesizer(pTweak, pTweak->featureMfcc.fftSize, 1.0, pTweak->featureMfcc.window);
		SC_Model *pModel;
		std::map<unsigned long int, std::pair<unsigned int, unsigned int> > featureMap;
		char audioFile[sclib::bufferSize];
		const char *resultPathName = env->GetStringUTFChars(resultPath, 0);
		char *path = sclib::makePath(resultPathName);
		bool res;

		setbuf(stdout, NULL); //get instant printf()s
		setbuf(stderr, NULL); //get instant printf()s
		std::ios::sync_with_stdio(); //syncronize couts and printf()s

		//release the strings previously copied from java
		env->ReleaseStringUTFChars(iniFile, iniFileName);
		env->ReleaseStringUTFChars(resultPath, resultPathName);

		pCorpus = new SC_Corpus_MPEG7(pTweak, 1.0, env, jStream);

		//set some parameters
		pTweak->modelHandler.foregroundModelType = (fullCovar==0) ? sclib::mtGMM_new : sclib::mtBGMM;
		pTweak->modelHandler.maxSpeakerModelOrder = sclib::getBetween((short int)(1), mixtureCount, (short int)(2048));
		if (pTweak->modelHandler.foregroundModelType == sclib::mtGMM_new) {
			pTweak->mixtureModelGmm.varianceLimit = sclib::max(0.0, varianceLimit);
			pTweak->mixtureModelGmm.maxEMiterations = sclib::getBetween(1, maxEmIterations, 1000);
			pTweak->mixtureModelGmm.EMthreshold = sclib::max(0.0, emThreshold);
			pTweak->mixtureModelGmm.weightLimit = 0.0;
		} else {
			pTweak->mixtureModelBgmm.varianceLimit = sclib::max(0.0, varianceLimit);
			pTweak->mixtureModelBgmm.maxEMiterations = sclib::getBetween(1, maxEmIterations, 1000);
			pTweak->mixtureModelBgmm.EMthreshold = sclib::max(0.0, emThreshold);
			pTweak->mixtureModelBgmm.weightLimit = 0.0;
			pTweak->mixtureModelBgmm.fullCovariance = true;
		}
		if (addPitch != 0) {
			pTweak->modelHandler.speakerModelFeature = sclib::featurePitch;
		}
		if (mfccOrder>0 && useMFCC!=0) {
			pTweak->modelHandler.speakerModelFeature |= sclib::featureMFCC;
		} else if (lpcOrder>0 && useMFCC==0) {
			pTweak->modelHandler.speakerModelFeature |= sclib::featureLPC;
		}
		pTweak->featureMfcc.frameSize = sclib::getBetween(1, frameSize, 3000);
		pTweak->featureMfcc.frameStep = sclib::getBetween(1, frameStep, pTweak->featureMfcc.frameSize);
		pTweak->featureMfcc.fftSize = sclib::checkFftSize(fftSize, pCorpus->getGT()->getConverter()->ms2sample(pTweak->featureMfcc.frameSize));
		pTweak->featureMfcc.filterBankSize = sclib::getBetween(1, filterbankSize, pTweak->featureMfcc.filterBankSize/2+1);
		pTweak->featureMfcc.MFCCorder = sclib::getBetween((unsigned short int)(0), mfccOrder, pTweak->featureMfcc.filterBankSize);
		pTweak->featureMfcc.preEmphasizeFactor = sclib::getBetween(0.0, preemphasis, 1.0);
		pTweak->featureMfcc.sclib_minFilterBankFrequency = sclib::getBetween(0.0, minFilterbankFrequency, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0);
		pTweak->featureMfcc.sclib_maxFilterBankFrequency = sclib::getBetween(pTweak->featureMfcc.sclib_minFilterBankFrequency, maxFilterbankFrequency, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0);
		pTweak->featureMfcc.addDeltaDeltas = false; //the folowwing parameters can't be made audible and thus are not changeable by the user
		pTweak->featureMfcc.addDeltas = false;
		pTweak->featureMfcc.coeffSelection[0] = 0xFF;
		pTweak->featureMfcc.coeffSelection[1] = 0xFF;
		pTweak->featureMfcc.coeffSelection[2] = 0xFF;
		pTweak->featureMfcc.coeffSelection[3] = 0xFF;
		pTweak->featureMfcc.coeffSelection[4] = 0xFF;
		pTweak->featureMfcc.coeffSelection[5] = 0xFF;
		pTweak->featureMfcc.coeffSelection[6] = 0xFF;
		pTweak->featureMfcc.coeffSelection[7] = 0xFF;
		pTweak->featureMfcc.CMN = false;
		pTweak->featureMfcc.dEnergy = false;
		pTweak->featureMfcc.highCut = 0.0;
		pTweak->featureMfcc.lowCut = 0.0;
		pTweak->featureMfcc.method = sclib::modeSClib;
		pTweak->featureMfcc.sclib_frequencyScale = sclib::getBetween(0, frequencyScale, 4);
		pTweak->featureMfcc.sclib_smoothing = sclib::smoothNone;
		pTweak->featureMfcc.window = sclib::getBetween(0, window, 3);
		synthesizer.setFftLength(pTweak->featureMfcc.fftSize);
		synthesizer.setOlaMaxIterations(olaMaxIterationCount);
		synthesizer.setOlaErrorTarget(olaErrorTarget);
		synthesizer.setTaperingLength(pTweak->transform.taperingLength);
		synthesizer.setTaperingMode(pTweak->featureMfcc.window);
		pTweak->featureLpc.lowCut = 0.0;
		pTweak->featureLpc.highCut = 0.0;
		pTweak->featureLpc.frameSize = sclib::getBetween(1, frameSize, 3000);;
		pTweak->featureLpc.frameStep = sclib::getBetween(1, frameStep, pTweak->featureLpc.frameSize);
		pTweak->featureLpc.LPCorder = sclib::getBetween((unsigned short int)(0), lpcOrder, (unsigned short int)(pCorpus->getGT()->getConverter()->ms2sample(frameSize)));
		pTweak->featureLpc.preEmphasizeFactor = sclib::getBetween(0.0, preemphasis, 1.0);
		pTweak->featureLpc.window = sclib::getBetween(0, window, 3);
		pTweak->featureLpc.computeGain = false;
		pTweak->featurePitch.frameSize = sclib::getBetween(1, pitchFrameSize, 3000);
		pTweak->featurePitch.frameStep = sclib::getBetween(1, pitchFrameStep, pTweak->featurePitch.frameSize);
		pTweak->featurePitch.method = sclib::modeESPS;
		pTweak->featurePitch.esps_cand_thresh = (float)(sclib::getBetween(0.01, pitchCandidateThresh,	0.99));
		pTweak->featurePitch.esps_lag_weight = (float)(sclib::getBetween(0.0, pitchLagWeight, 1.0));
		pTweak->featurePitch.esps_freq_weight = (float)(sclib::getBetween(0.0, pitchFrequencyWeight, 1.0));
		pTweak->featurePitch.esps_trans_cost = (float)(sclib::getBetween(0.0, pitchTransitionCost, 1.0));
		pTweak->featurePitch.esps_trans_amp = (float)(sclib::getBetween(0.0, pitchTransitionAmpModCost, 100.0));
		pTweak->featurePitch.esps_trans_spec = (float)(sclib::getBetween(0.0, pitchTransitionSpecModCost, 100.0));
		pTweak->featurePitch.esps_voice_bias = (float)(sclib::getBetween(-1.0, pitchVoiceBias, 1.0));
		pTweak->featurePitch.esps_double_cost = (float)(sclib::getBetween(0.0, pitchDoublingCost, 10.0));
		pTweak->featurePitch.esps_min_f0 = (float)(sclib::getBetween(10.0, pitchMinF0, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0));
		pTweak->featurePitch.esps_max_f0 = (float)(sclib::getBetween((double)(pTweak->featurePitch.esps_min_f0), pitchMaxF0, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0));
		pTweak->featurePitch.esps_n_cands = sclib::getBetween(3, pitchCandidateCount, 100);
		pTweak->featurePitch.esps_wind_dur = (float)(sclib::getBetween(10.0/(double)(pCorpus->getGT()->getAudioSampleRate()), pitchCandidateWindowSize, 0.1));

		//load the signal
		sceneStart = 0;
		sceneEnd = pCorpus->getGT()->getAudioSampleCount()-1;
		pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pTweak->modelHandler.speakerModelFeature);
		extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list

		//do modeling
		extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), pTweak->modelHandler.speakerModelFeature);
		pConcat = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), featureMap, pTweak->modelHandler.speakerModelFeature);
		//printf("\n(%#6.5f%% outliers removed from data)", extractor.removeOutliers(pConcat, pTweak->modelHandler.outlierRemovalMode));
		if (pTweak->modelHandler.foregroundModelType == sclib::mtGMM_new) {
			pModel = new SC_MixtureModel_GMM(pTweak, pTweak->modelHandler.maxSpeakerModelOrder, pConcat->Col);
			pModel->TrainModel(pConcat, 0);
		} else {
			pModel = new SC_MixtureModel_bGMM(pTweak, pTweak->modelHandler.maxSpeakerModelOrder, pConcat->Col, true, true);
			((SC_MixtureModel_bGMM*)pModel)->TrainModel(pConcat, 0, -1, true, 1.0, false, -1.0);
		}
		
		//resynthesize modelled speech
		pTmp = pModel->drawSamplesFromDistribution(pConcat->Row);
		pTmp->Hdr = pConcat->Hdr;
		MFree_0D(pModel);
		MFree_0D(pConcat);
		if (featureMap.count(sclib::featurePitch) > 0) {
			pPitch = extractor.getSubSet(pTmp, featureMap[sclib::featurePitch].first, featureMap[sclib::featurePitch].second);
			pPitch->Hdr = pFeatures[sclib::bitPosition(sclib::featurePitch)]->Hdr;
			extractor.cleanPitch(pPitch, 0, pTweak->featurePitch.esps_min_f0);
		} else {
			pPitch = NULL;
		}
		if (featureMap.count(sclib::featureMFCC) > 0) {
			pFeature = extractor.getSubSet(pTmp, featureMap[sclib::featureMFCC].first, featureMap[sclib::featureMFCC].second);
			pFeature->Hdr = pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr;
		} else if (featureMap.count(sclib::featureLPC) > 0) {
			pFeature = extractor.getSubSet(pTmp, featureMap[sclib::featureLPC].first, featureMap[sclib::featureLPC].second);
			pFeature->Hdr = pFeatures[sclib::bitPosition(sclib::featureLPC)]->Hdr;
		} else {
			pFeature = pPitch;
		}
		MFree_0D(pTmp);
		sprintf(audioFile, "%s%s%s", path, tmpnam(NULL)+1, "wav"); //get a new temporary fileName
		res = synthesizer.feature2wav(audioFile, pFeature, pPitch, true);
		if (res == false) {
			sprintf(error, "%s%s", "ERROR: no result returned by resynthesis, maybe no features given/extracted?");
		}

		//clean up
		if (pPitch != pFeature) {
			MFree_0D(pPitch);
		}
		MFree_0D(pFeature);
		MFree_1D(path);
		extractor.freeFeatures(pFeatures);
		MFree_0D(pSignal);
		MFree_0D(pCorpus);

		return (res==true) ? env->NewStringUTF(audioFile) : env->NewStringUTF(error);	
	} catch (const char *errStr) {
		sprintf(error, "%s%s", "ERROR: ", errStr);

		return env->NewStringUTF(error);
	}
}

JNIEXPORT jstring JNICALL Java_jni_SC_1Wrapper_wav2hmm(JNIEnv *env, jobject caller, jobject jStream, jstring iniFile, jstring resultPath, jint frameCountForSplit, jint frameCountForOverlap, jint stateCount, jint mixturesPerState, jstring transitionStructure, jboolean isLeft2Right, jint maxIterations, jint frameSize, jint frameStep, jboolean useMFCC, jboolean addPitch, jint mfccOrder, jint lpcOrder, jdouble preemphasis, jint window, jint fftSize, jint frequencyScale, jint filterbankSize, jdouble minFilterbankFrequency, jdouble maxFilterbankFrequency, jdouble olaErrorTarget, jint olaMaxIterationCount, jint pitchFrameSize, jint pitchFrameStep, jdouble pitchCandidateThresh, jdouble pitchLagWeight, jdouble pitchFrequencyWeight, jdouble pitchTransitionCost, jdouble pitchTransitionAmpModCost, jdouble pitchTransitionSpecModCost, jdouble pitchVoiceBias, jdouble pitchDoublingCost, jdouble pitchMinF0, jdouble pitchMaxF0, jint pitchCandidateCount, jdouble pitchCandidateWindowSize) {
	sclib::errorHandlerThrows(false, true);
	char error[sclib::bufferSize];

	try {
		SC_Corpus *pCorpus;
		SC_Signal *pSignal;
		long int sceneStart, sceneEnd;
		SV_Data **pFeatures, *pConcat, *pTmp, *pFeature, *pPitch;
		const char *iniFileName = env->GetStringUTFChars(iniFile, 0);
		const char *transitionStructureString = env->GetStringUTFChars(transitionStructure, 0);
		SC_TweakableParameters *pTweak = new SC_TweakableParameters(iniFileName);
		SC_FeatureHandler extractor(pTweak, true);
		SC_Synthesis synthesizer(pTweak, pTweak->featureMfcc.fftSize, 1.0, pTweak->featureMfcc.window);
		SC_Model *pModel;
		std::map<unsigned long int, std::pair<unsigned int, unsigned int> > featureMap;
		char audioFile[sclib::bufferSize];
		const char *resultPathName = env->GetStringUTFChars(resultPath, 0);
		char *path = sclib::makePath(resultPathName);
		bool res;

		setbuf(stdout, NULL); //get instant printf()s
		setbuf(stderr, NULL); //get instant printf()s
		std::ios::sync_with_stdio(); //syncronize couts and printf()s

		//release the strings previously copied from java
		env->ReleaseStringUTFChars(iniFile, iniFileName);
		env->ReleaseStringUTFChars(resultPath, resultPathName);

		pCorpus = new SC_Corpus_MPEG7(pTweak, 1.0, env, jStream);

		//set some parameters
		pTweak->modelHandler.foregroundModelType = sclib::mtHMM;
		pTweak->modelHmm.leftToRight = (isLeft2Right!=0)?true:false;
		pTweak->modelHmm.maxIterations = sclib::getBetween(1, maxIterations, 1000);
		pTweak->modelHmm.mixturesPerState = sclib::getBetween(1, mixturesPerState, 2048);
		pTweak->modelHmm.stateCount = sclib::getBetween(1, stateCount, 2048);
		MFree_1D(pTweak->modelHmm.transitionStructure);
		MArray_1D(pTweak->modelHmm.transitionStructure, sclib::bufferSize, char, "Java_jni_SC_1Wrapper_wav2hmm: pTweak->modelHmm.transitionStructure");
		printf(pTweak->modelHmm.transitionStructure, "%s", transitionStructureString);
		env->ReleaseStringUTFChars(transitionStructure, transitionStructureString);
		if (addPitch != 0) {
			pTweak->modelHandler.speakerModelFeature = sclib::featurePitch;
		}
		if (mfccOrder>0 && useMFCC!=0) {
			pTweak->modelHandler.speakerModelFeature |= sclib::featureMFCC;
		} else if (lpcOrder>0 && useMFCC==0) {
			pTweak->modelHandler.speakerModelFeature |= sclib::featureLPC;
		}
		pTweak->featureMfcc.frameSize = sclib::getBetween(1, frameSize, 3000);
		pTweak->featureMfcc.frameStep = sclib::getBetween(1, frameStep, pTweak->featureMfcc.frameSize);
		pTweak->featureMfcc.fftSize = sclib::checkFftSize(fftSize, pCorpus->getGT()->getConverter()->ms2sample(pTweak->featureMfcc.frameSize));
		pTweak->featureMfcc.filterBankSize = sclib::getBetween(1, filterbankSize, pTweak->featureMfcc.filterBankSize/2+1);
		pTweak->featureMfcc.MFCCorder = sclib::getBetween((unsigned short int)(0), mfccOrder, pTweak->featureMfcc.filterBankSize);
		pTweak->featureMfcc.preEmphasizeFactor = sclib::getBetween(0.0, preemphasis, 1.0);
		pTweak->featureMfcc.sclib_minFilterBankFrequency = sclib::getBetween(0.0, minFilterbankFrequency, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0);
		pTweak->featureMfcc.sclib_maxFilterBankFrequency = sclib::getBetween(pTweak->featureMfcc.sclib_minFilterBankFrequency, maxFilterbankFrequency, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0);
		pTweak->featureMfcc.addDeltaDeltas = false; //the folowwing parameters can't be made audible and thus are not changeable by the user
		pTweak->featureMfcc.addDeltas = false;
		pTweak->featureMfcc.coeffSelection[0] = 0xFF;
		pTweak->featureMfcc.coeffSelection[1] = 0xFF;
		pTweak->featureMfcc.coeffSelection[2] = 0xFF;
		pTweak->featureMfcc.coeffSelection[3] = 0xFF;
		pTweak->featureMfcc.coeffSelection[4] = 0xFF;
		pTweak->featureMfcc.coeffSelection[5] = 0xFF;
		pTweak->featureMfcc.coeffSelection[6] = 0xFF;
		pTweak->featureMfcc.coeffSelection[7] = 0xFF;
		pTweak->featureMfcc.CMN = false;
		pTweak->featureMfcc.dEnergy = false;
		pTweak->featureMfcc.highCut = 0.0;
		pTweak->featureMfcc.lowCut = 0.0;
		pTweak->featureMfcc.method = sclib::modeSClib;
		pTweak->featureMfcc.sclib_frequencyScale = sclib::getBetween(0, frequencyScale, 4);
		pTweak->featureMfcc.sclib_smoothing = sclib::smoothNone;
		pTweak->featureMfcc.window = sclib::getBetween(0, window, 3);
		synthesizer.setFftLength(pTweak->featureMfcc.fftSize);
		synthesizer.setOlaMaxIterations(olaMaxIterationCount);
		synthesizer.setOlaErrorTarget(olaErrorTarget);
		synthesizer.setTaperingLength(pTweak->transform.taperingLength);
		synthesizer.setTaperingMode(pTweak->featureMfcc.window);
		pTweak->featureLpc.lowCut = 0.0;
		pTweak->featureLpc.highCut = 0.0;
		pTweak->featureLpc.frameSize = sclib::getBetween(1, frameSize, 3000);;
		pTweak->featureLpc.frameStep = sclib::getBetween(1, frameStep, pTweak->featureLpc.frameSize);
		pTweak->featureLpc.LPCorder = sclib::getBetween((unsigned short int)(0), lpcOrder, (unsigned short int)(pCorpus->getGT()->getConverter()->ms2sample(frameSize)));
		pTweak->featureLpc.preEmphasizeFactor = sclib::getBetween(0.0, preemphasis, 1.0);
		pTweak->featureLpc.window = sclib::getBetween(0, window, 3);
		pTweak->featureLpc.computeGain = false;
		pTweak->featurePitch.frameSize = sclib::getBetween(1, pitchFrameSize, 3000);
		pTweak->featurePitch.frameStep = sclib::getBetween(1, pitchFrameStep, pTweak->featurePitch.frameSize);
		pTweak->featurePitch.method = sclib::modeESPS;
		pTweak->featurePitch.esps_cand_thresh = (float)(sclib::getBetween(0.01, pitchCandidateThresh,	0.99));
		pTweak->featurePitch.esps_lag_weight = (float)(sclib::getBetween(0.0, pitchLagWeight, 1.0));
		pTweak->featurePitch.esps_freq_weight = (float)(sclib::getBetween(0.0, pitchFrequencyWeight, 1.0));
		pTweak->featurePitch.esps_trans_cost = (float)(sclib::getBetween(0.0, pitchTransitionCost, 1.0));
		pTweak->featurePitch.esps_trans_amp = (float)(sclib::getBetween(0.0, pitchTransitionAmpModCost, 100.0));
		pTweak->featurePitch.esps_trans_spec = (float)(sclib::getBetween(0.0, pitchTransitionSpecModCost, 100.0));
		pTweak->featurePitch.esps_voice_bias = (float)(sclib::getBetween(-1.0, pitchVoiceBias, 1.0));
		pTweak->featurePitch.esps_double_cost = (float)(sclib::getBetween(0.0, pitchDoublingCost, 10.0));
		pTweak->featurePitch.esps_min_f0 = (float)(sclib::getBetween(10.0, pitchMinF0, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0));
		pTweak->featurePitch.esps_max_f0 = (float)(sclib::getBetween((double)(pTweak->featurePitch.esps_min_f0), pitchMaxF0, (double)(pCorpus->getGT()->getAudioSampleRate())/2.0));
		pTweak->featurePitch.esps_n_cands = sclib::getBetween(3, pitchCandidateCount, 100);
		pTweak->featurePitch.esps_wind_dur = (float)(sclib::getBetween(10.0/(double)(pCorpus->getGT()->getAudioSampleRate()), pitchCandidateWindowSize, 0.1));

		//load the signal
		sceneStart = 0;
		sceneEnd = pCorpus->getGT()->getAudioSampleCount()-1;
		pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pTweak->modelHandler.speakerModelFeature);
		extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list

		//do modeling
		extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), pTweak->modelHandler.speakerModelFeature);
		pConcat = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), featureMap, pTweak->modelHandler.speakerModelFeature);
		//printf("\n(%#6.5f%% outliers removed from data)", extractor.removeOutliers(pConcat, pTweak->modelHandler.outlierRemovalMode));
		int splitCount = sclib::getBetween(1, frameCountForSplit, pConcat->Row/2);
		int overlap = sclib::getBetween(1, frameCountForOverlap, splitCount);
		pTmp = extractor.splitFeatureSet(pConcat, splitCount, overlap);
		pTmp->Hdr = pConcat->Hdr;
		pModel = new SC_Model_HMM(pTweak, pTweak->modelHmm.stateCount, pTweak->modelHmm.transitionStructure, pTweak->modelHmm.mixturesPerState, false, pTweak->modelHmm.leftToRight, pTweak->modelHmm.maxIterations, false);
		pModel->TrainModel(pTmp, 0);
		sclib::destructLinkedList(pTmp);

		//resynthesize modelled speech
		pTmp = pModel->drawSamplesFromDistribution(pConcat->Row);
		pTmp->Hdr = pConcat->Hdr;
		MFree_0D(pModel);
		MFree_0D(pConcat);
		if (featureMap.count(sclib::featurePitch) > 0) {
			pPitch = extractor.getSubSet(pTmp, featureMap[sclib::featurePitch].first, featureMap[sclib::featurePitch].second);
			pPitch->Hdr = pFeatures[sclib::bitPosition(sclib::featurePitch)]->Hdr;
			extractor.cleanPitch(pPitch, 0, pTweak->featurePitch.esps_min_f0);
		} else {
			pPitch = NULL;
		}
		if (featureMap.count(sclib::featureMFCC) > 0) {
			pFeature = extractor.getSubSet(pTmp, featureMap[sclib::featureMFCC].first, featureMap[sclib::featureMFCC].second);
			pFeature->Hdr = pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr;
		} else if (featureMap.count(sclib::featureLPC) > 0) {
			pFeature = extractor.getSubSet(pTmp, featureMap[sclib::featureLPC].first, featureMap[sclib::featureLPC].second);
			pFeature->Hdr = pFeatures[sclib::bitPosition(sclib::featureLPC)]->Hdr;
		} else {
			pFeature = pPitch;
		}
		MFree_0D(pTmp);
		sprintf(audioFile, "%s%s%s", path, tmpnam(NULL)+1, "wav"); //get a new temporary fileName
		res = synthesizer.feature2wav(audioFile, pFeature, pPitch, true);
		if (res == false) {
			sprintf(error, "%s%s", "ERROR: no result returned by resynthesis, maybe no features given/extracted?");
		}

		//clean up
		if (pPitch != pFeature) {
			MFree_0D(pPitch);
		}
		MFree_0D(pFeature);
		MFree_1D(path);
		extractor.freeFeatures(pFeatures);
		MFree_0D(pSignal);
		MFree_0D(pCorpus);

		return (res==true) ? env->NewStringUTF(audioFile) : env->NewStringUTF(error);	
	} catch (const char *errStr) {
		sprintf(error, "%s%s", "ERROR: ", errStr);
		jstring jError = env->NewStringUTF(error);
		return env->NewStringUTF(error);
	}
}
#endif

//====================================================================================================================
//  A bunch of memory-deallocation methods to be called from a (java) application and an associated enumeration to 
//  tell the methods which type the given object has; at the moment, only those types currently passed to videana are
//  recognized in the enumkeration; TODO: enhance it with all exported types
//====================================================================================================================
union BaseType {
  long *pL;
	long **ppL;
	double *pD;
	double **ppD;
	SC_TweakableParameters *pTweak;
	SC_TweakableParameters **ppTweak;
	SV_Data *pData;
	SV_Data **ppData;
	SC_ModelHandler *pModelHandler;
	SC_ModelHandler **ppModelHandler;
	SC_Model *pModel;
	SC_Model **ppModel;
	SC_Centroid *pCentroid;
	SC_Centroid **ppCentroid;
	SC_Signature *pSignature;
	SC_Signature **ppSignature;
	SC_DistanceMeasures *pDist;
	SC_DistanceMeasures **ppDist;
};

void deleteScalar(void *scalar, SC_Types baseType) {
	BaseType bt;

	switch (baseType) {
		case scLong:
			bt.pL = (long*)(scalar);
			MFree_0D(bt.pL);
			break;
		case scDouble:
			bt.pD = (double*)(scalar);
			MFree_0D(bt.pD);
			break;
		case scTweakableParameters:
			bt.pTweak = (SC_TweakableParameters*)(scalar);
			MFree_0D(bt.pTweak);
			break;
		case svData:
			bt.pData = (SV_Data*)(scalar);
			MFree_0D(bt.pData);
			break;
		case scModelHandler:
			bt.pModelHandler = (SC_ModelHandler*)(scalar);
			MFree_0D(bt.pModelHandler);
			break;
		case scModel:
			bt.pModel = (SC_Model*)(scalar);
			MFree_0D(bt.pModel);
			break;
		case scCentroid:
			bt.pCentroid = (SC_Centroid*)(scalar);
			MFree_0D(bt.pCentroid);
			break;
		case scSignature:
			bt.pSignature = (SC_Signature*)(scalar);
			MFree_0D(bt.pSignature);
			break;
		case scDistanceMeasure:
			bt.pDist = (SC_DistanceMeasures*)(scalar);
			MFree_0D(bt.pDist);
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "Unknown base type");
			break;
	}

	return;
}

void deleteArray(void *array, SC_Types baseType) {
	BaseType bt;

	switch (baseType) {
		case scLong:
			bt.pL = (long*)(array);
			MFree_1D(bt.pL);
			break;
		case scDouble:
			bt.pD = (double*)(array);
			MFree_1D(bt.pD);
			break;
		case scTweakableParameters:
			bt.pTweak = (SC_TweakableParameters*)(array);
			MFree_1D(bt.pTweak);
			break;
		case svData:
			bt.pData = (SV_Data*)(array);
			MFree_1D(bt.pData);
			break;
		case scModelHandler:
			bt.pModelHandler = (SC_ModelHandler*)(array);
			MFree_1D(bt.pModelHandler);
			break;
		case scModel:
			bt.pModel = (SC_Model*)(array);
			MFree_1D(bt.pModel);
			break;
		case scCentroid:
			bt.pCentroid = (SC_Centroid*)(array);
			MFree_1D(bt.pCentroid);
			break;
		case scSignature:
			bt.pSignature = (SC_Signature*)(array);
			MFree_1D(bt.pSignature);
			break;
		case scDistanceMeasure:
			bt.pDist = (SC_DistanceMeasures*)(array);
			MFree_1D(bt.pDist);
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "Unknown base type");
			break;
	}

	return;
}

void deleteMatrix(void **matrix, SC_Types baseType) {
	BaseType bt;

	switch (baseType) {
		case scLong:
			bt.ppL = (long**)(matrix);
			MFree_2D(bt.ppL);
			break;
		case scDouble:
			bt.ppD = (double**)(matrix);
			MFree_2D(bt.ppD);
			break;
		case scTweakableParameters:
			bt.ppTweak = (SC_TweakableParameters**)(matrix);
			MFree_2D(bt.ppTweak);
			break;
		case svData:
			bt.ppData = (SV_Data**)(matrix);
			MFree_2D(bt.ppData);
			break;
		case scModelHandler:
			bt.ppModelHandler = (SC_ModelHandler**)(matrix);
			MFree_2D(bt.ppModelHandler);
			break;
		case scModel:
			bt.ppModel = (SC_Model**)(matrix);
			MFree_2D(bt.ppModel);
			break;
		case scCentroid:
			bt.ppCentroid = (SC_Centroid**)(matrix);
			MFree_2D(bt.ppCentroid);
			break;
		case scSignature:
			bt.ppSignature = (SC_Signature**)(matrix);
			MFree_2D(bt.ppSignature);
			break;
		case scDistanceMeasure:
			bt.ppDist = (SC_DistanceMeasures**)(matrix);
			MFree_2D(bt.ppDist);
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "Unknown base type");
			break;
	}

	return;
}

//====================================================================================================================
//  Methods to work with array of pointers (to objects) from within java applications
//====================================================================================================================
void* constructPointerArray(SC_Types baseType, int dim) {
	BaseType bt;

	switch (baseType) {
		case scTweakableParameters:
			MArray_1D(bt.ppTweak, dim, SC_TweakableParameters*, "constructPointerArray: ppTweak");
			return bt.ppTweak;
		case svData:
			MArray_1D(bt.ppData, dim, SV_Data*, "constructPointerArray: ppData");
			return bt.ppData;
		case scModelHandler:
			MArray_1D(bt.ppModelHandler, dim, SC_ModelHandler*, "constructPointerArray: ppModelHandler");
			return bt.ppModelHandler;
		case scModel:
			MArray_1D(bt.ppModel, dim, SC_Model*, "constructPointerArray: ppModel");
			return bt.ppModel;
		case scCentroid:
			MArray_1D(bt.ppCentroid, dim, SC_Centroid*, "constructPointerArray: ppCentroid");
			return bt.ppCentroid;
		case scSignature:
			MArray_1D(bt.ppSignature, dim, SC_Signature*, "constructPointerArray: ppSignature");
			return bt.ppSignature;
		case scDistanceMeasure:
			MArray_1D(bt.ppDist, dim, SC_DistanceMeasures*, "constructPointerArray: ppDist");
			return bt.ppDist;
		default:
			REPORT_ERROR(SVLIB_BadArg, "Unknown base type");
			break;
	}

	return NULL;
}

void* getPointerArrayElement(void **array, int idx) {
	return array[idx];
}

int setPointerArrayElement(void **array, int idx, void *newPointer) {
	array[idx] = newPointer;
	return 1;
}

void deletePointerArray(void **array, SC_Types baseType) {
	BaseType bt;

	switch (baseType) {
		case scTweakableParameters:
			bt.ppTweak = (SC_TweakableParameters**)(array);
			MFree_1D(bt.ppTweak);
			break;
		case svData:
			bt.ppData = (SV_Data**)(array);
			MFree_1D(bt.ppData);
			break;
		case scModelHandler:
			bt.ppModelHandler = (SC_ModelHandler**)(array);
			MFree_1D(bt.ppModelHandler);
			break;
		case scModel:
			bt.ppModel = (SC_Model**)(array);
			MFree_1D(bt.ppModel);
			break;
		case scCentroid:
			bt.ppCentroid = (SC_Centroid**)(array);
			MFree_1D(bt.ppCentroid);
			break;
		case scSignature:
			bt.ppSignature = (SC_Signature**)(array);
			MFree_1D(bt.ppSignature);
			break;
		case scDistanceMeasure:
			bt.ppDist = (SC_DistanceMeasures**)(array);
			MFree_1D(bt.ppDist);
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "Unknown base type");
			break;
	}

	return;
}

//====================================================================================================================
//  Get instant printf()'s from e.g. Java
//====================================================================================================================
void getInstantPrintfs(void) {
	setbuf(stdout, NULL);
	setbuf(stderr, NULL);

	return;
}

//====================================================================================================================
//  Interface to the SC_TweakableParameters class
//====================================================================================================================
SC_TweakableParameters* scTweakableParameters_construct(const char *fileName, SC_Bool verbose) {
	return new SC_TweakableParameters(fileName, (verbose!=0)?true:false);
}

SC_Bool scTweakableParameters_setByName(SC_TweakableParameters *pTweak, const char *parameterName, const char *value) {
	return (pTweak->setByName(parameterName, value) == true) ? 1 : 0;
}

SC_TweakableParameters* createTweak(const char *fileName) { //deprectated, just exists for compatibility reasons
	return scTweakableParameters_construct(fileName, true);
}

SC_Bool setParameterByName(SC_TweakableParameters *pTweak, char *parameterName, char *value) { //deprectated, just exists for compatibility reasons
	return scTweakableParameters_setByName(pTweak, parameterName, value);
}

//====================================================================================================================
//  Interface to the SV_Data class
//====================================================================================================================
SV_Data* svData_construct(int rows, int cols) {
	if (rows > 0 && cols > 0) {
		return new SV_Data(rows, cols);
	} else {
		return new SV_Data();
	}
}

SC_Bool svData_setCol(SV_Data *pData, int col) {
	pData->Col = col;
	return 1;
}

int svData_getCol(SV_Data *pData) {
	return pData->Col;
}

SC_Bool svData_setRow(SV_Data *pData, int row) {
	pData->Row = row;
	return 1;
}

int svData_getRow(SV_Data *pData) {
	return pData->Row;
}

SC_Bool svData_alloc(SV_Data *pData) {
	pData->Alloc();
	return 1;
}

SC_Bool svData_setMat(SV_Data *pData, int row, int col, float value) {
	pData->Mat[row][col] = value;
	return 1;
}

float svData_getMat(SV_Data *pData, int row, int col) {
	return pData->Mat[row][col];
}

SC_Bool svData_setNext(SV_Data *pData, SV_Data *pNext) {
	pData->Next = pNext;
	return 1;
}

SV_Data* svData_getNext(SV_Data *pData) {
	return pData->Next;
}

SV_Data* svData_mergeData(SV_Data *pData, int maxSegments) {
	return pData->MergeData(maxSegments);
}

//====================================================================================================================
//  Interface to the SC_ModelHandler class
//====================================================================================================================
SC_ModelHandler* scModelHandler_construct(SC_TweakableParameters *pTweak, SC_Bool verbose) {
	return new SC_ModelHandler(pTweak, (verbose!=0)?true:false);
}

int scModelHandler_guessModelOrder(SC_ModelHandler *pModelHandler, SV_Data *pFeatures, SC_Model *pBackground, int modelType, int minOrder, int maxOrder, int segmentsToMerge) {
	return pModelHandler->guessModelOrder(pFeatures, pBackground, modelType, minOrder, maxOrder, segmentsToMerge);
}

SC_Model* scModelHandler_buildModel(SC_ModelHandler *pModelHandler, SV_Data *pFeatures, SC_Model *pBackground, int modelOrder, int modelType, int segmentsToMerge) {
	return pModelHandler->buildModel(pFeatures, pBackground, modelOrder, modelType, segmentsToMerge);
}

double scModelHandler_testModel(SC_ModelHandler *pModelHandler, SC_Model *pModel, SV_Data *pFeatures, int segmentsToMerge, SC_Bool forceNoNormalization, SC_Model *pBackground) {
	return pModelHandler->testModel(pModel, pFeatures, segmentsToMerge, (forceNoNormalization!=0)?true:false, pBackground);
}

int scModelHandler_saveModel(SC_ModelHandler *pModelHandler, const char *fileName, SC_Model *pModel, SC_Bool useDebugDir) {
	return pModelHandler->saveModel(fileName, pModel, (useDebugDir!=0)?true:false);
}

SC_Model* scModelHandler_loadModel(SC_ModelHandler *pModelHandler, const char *fileName, int modelType) {
	return pModelHandler->loadModel(fileName, modelType);
}

SC_Bool scModelHandler_setTweak(SC_ModelHandler *pModelHandler, SC_TweakableParameters *pTweak) {
	pModelHandler->setTweak(pTweak);
	return 1;
}

//====================================================================================================================
//  Interface to SC_Centroid_Point class
//====================================================================================================================
SC_Centroid_Point* scCentroidPoint_construct(SC_TweakableParameters *pTweak, int dim, double *coordinate, SC_Bool justLink) {
	return new SC_Centroid_Point(pTweak, dim, coordinate, (justLink!=0)?true:false);
}

double* scCentroidPoint_getCoordinate(SC_Centroid_Point* pCentroid) {
	return pCentroid->getCoordinate();
}

int scCentroidPoint_getDim(SC_Centroid_Point* pCentroid) {
	return pCentroid->getDim();
}

int scCentroidPoint_setCoordinate(SC_Centroid_Point* pCentroid, double *newCoordinate) {
	pCentroid->setCoordinate(newCoordinate);
	return 1;
}

int scCentroidPoint_setDim(SC_Centroid_Point* pCentroid, int newDim) {
	pCentroid->setDim(newDim);
	return 1;
}

double scCentroidPoint_getDistance(SC_Centroid_Point* pCentroid, SC_Centroid *secondCentroid) {
	return pCentroid->getDistance(secondCentroid);
}

//====================================================================================================================
//  Interface to SC_Signature class
//====================================================================================================================
SC_Signature* scSignature_construct(SC_Centroid **centroids, double *weights, int n, SC_Bool justLinkWeights, double smallestUnnormalizedWeight) {
	return new SC_Signature(centroids, weights, n, (justLinkWeights!=0)?true:false, smallestUnnormalizedWeight);
}

//====================================================================================================================
//  Interface to SC_DistanceMeasures class
//====================================================================================================================
SC_DistanceMeasures* sc_DistanceMeasures_construct(SC_TweakableParameters* pTweak, SC_MatrixFunctions* pMatrixFunc, SC_Bool verbose) {
	return new SC_DistanceMeasures(pTweak, pMatrixFunc, (verbose!=0)?true:false);
}

double scDistanceMeasures_EMD(SC_Signature* pSignature1, SC_Signature* pSignature2, SC_TweakableParameters *pTweak) {
	double dist = SC_DistanceMeasures::EMD(pSignature1, pSignature2, pTweak);

	return dist;
}
