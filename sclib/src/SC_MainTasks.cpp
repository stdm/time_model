/**************************************************************************/
/*	This class is used to implement some common tasks that would normally */
/*  by solved by writing a small program that uses the sclib to do some   */
/*  work, such as training a specific algorithm or building specified     */
/*  models. In order to minimize the amount of needed projects of this    */
/*  type, the same functionality can be coded in one single "main"        */
/*  function here that can be feeded with parameters from a uniform       */
/*  external script that can such accomplish several tasks. In the long   */
/*  run this class could exchange programs like trainer, modeller, scivo, */
/*  ...                                                                   */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 04.03.2008																								*/
/**************************************************************************/

#include <list>
#include <vector>
#include "SC_MainTasks.h"
#include "SC_Lib.h"
#include "SC_Matlab.h"


//====================================================================================================================
//	constructor
//====================================================================================================================
SC_MainTasks::SC_MainTasks(SC_TweakableParameters* pTweak, const char *audioFileName, const char *segmentationFileName, const char *sceneFileName, double videoFrameRate) {
	initParameters(pTweak, audioFile, segmentationFile, sceneFile, videoFrameRate);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_MainTasks::~SC_MainTasks() {

}

//====================================================================================================================
//	(re)-initialize the class members
//====================================================================================================================
void SC_MainTasks::initParameters(SC_TweakableParameters* pTweak, const char *audioFileName, const char *segmentationFileName, const char *sceneFileName, double videoFrameRate) {
	this->pTweak = pTweak;
	this->videoFrameRate = videoFrameRate;

	sprintf(this->audioFile, "%s", (audioFileName!=NULL)?audioFileName:"");
	sprintf(this->sceneFile, "%s", (sceneFileName!=NULL)?sceneFileName:"");
	sprintf(this->segmentationFile, "%s", (segmentationFileName!=NULL)?segmentationFileName:"");
	
	return;
}

//====================================================================================================================
//	training script for the LZL audio type classifier
//====================================================================================================================
bool SC_MainTasks::trainAudioTypeLZL() {
	bool finalResult = true, existsFeatureFiles = true;
  unsigned long int segmentStart, segmentEnd;
  long int res;
	char fileNamePrefix[sclib::bufferSize], ext[sclib::bufferSize], *temp = NULL;
  time_t startTime	= time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Signal *pSignal = NULL;
  SV_Data **pFeatures = NULL;
  SC_FeatureHandler *pExtractor = NULL;
  SC_SegmentationHandler *pSegmenter = NULL;
  SC_Segmentation_AudioType_LZL *pATClassificator = NULL;

  printf("\nThis is Trainer: Learning algorithm training via the SC_Lib");
  printf("\n============================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
	this->pTweak->groundTruth.storeProbabilityInformation = false; //otherwise the memory-consumption would be to high
	this->pTweak->segmentationHandler.audioTypeMode = sclib::algorithm_ac_LZL; //so that the correct features get extracted
  pExtractor = new SC_FeatureHandler(this->pTweak);
  pSegmenter = new SC_SegmentationHandler(this->pTweak);
  pATClassificator = new SC_Segmentation_AudioType_LZL(this->pTweak);
  pCorpus = new SC_Corpus_MAC(this->pTweak, this->audioFile);
	printf("done!\n");

  printf("\nPhase 1:\tLoading and preprocessing data...");
  unsigned long int approxLoadLength = pCorpus->getGT()->getConverter()->ms2sample(10*60*1000); //10 minutes
	long int types[5] = {sclib::atPureSpeech, sclib::atNoisySpeech, sclib::atMusic, sclib::atBackground, sclib::atAction};

	//check if the feature files are already written
	sprintf(fileNamePrefix, "%s", this->pTweak->segmentationAudioTypeLzl.featureFileName);
	for (int x = 1; x <= 5; x++) {
		sprintf(ext, ".%d\0", x);
		temp = sclib::exchangeFileExtension(fileNamePrefix, ext);
		if (!sclib::fileExists(temp)) {
			existsFeatureFiles = false;
			MFree_1D(temp);
			break; 
		}
		MFree_1D(temp);
	}

  if (existsFeatureFiles == false) {
    //loop over all desired types
    for (int x = 0; x < 5; x++) {
      printf("\n Processing type %d...", types[x]);
      segmentStart = 0;
      
      //loop over all samples in the corpus
      while (segmentStart < pCorpus->getGT()->getAudioSampleCount()) {
        //load signal
        printf("\n    Phase I:\tLoading Signal...");
 	      pSignal = ((SC_Corpus_MAC*)pCorpus)->loadSignal(segmentStart, segmentEnd, approxLoadLength, types[x]);
        if (pSignal != NULL) {
          printf("done! We've reached %#5.1f%%", (double)(segmentEnd) / (double)(pCorpus->getGT()->getAudioSampleCount()) * 100.0);

          //extract features
          printf("\n    Phase II:\tExtracting Features...");
          pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, pSegmenter->getUsedFeatures(), pSegmenter->getSpecialFeatureParameters(), true);
          //if (!pExtractor->checkFeatures(pFeatures)) {
          //  printf("\n Error in the features!!!");
          //}
          MFree_0D(pSignal);
	        printf("done!");

          //silence detection
          printf("\n    Phase III:\tRemoving Silence...");
	        pSegmenter->silenceDetection(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
	        printf("done!");

          //feed learning algorithm and free features
					pExtractor->prepareFeatureSet(pFeatures, pSegmenter->getSpecialFeatureParameters(sclib::algorithm_ac_LZL));
          pATClassificator->partiallyLoadTrainingData(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures, types[x]);
          pExtractor->freeFeatures(pFeatures);

          segmentStart = segmentEnd + 1;
        } else {
          printf("failed! No more samples of this type to load, proceeding to the next one");
          segmentStart = pCorpus->getGT()->getAudioSampleCount(); //to finish the while-loop
        }
      }
      printf("\n done!\n");
    }
    printf("\ndone!\n");
  } else {
    printf("done! (will load cached features...)\n");
  }

  printf("\nPhase 2:\tTraining the Algorithm...");
  res = pATClassificator->trainClassifier(pCorpus->getGT(), 0, pCorpus->getGT()->getAudioSampleCount()-1, NULL);  
  printf("\ndone! Result: %d\n", res);

	printf("\nPhase 3:\tCleaning up...");
  MFree_0D(pCorpus);
  MFree_0D(pExtractor);
  MFree_0D(pSegmenter);
  MFree_0D(pATClassificator);

	printf("done!");
	printf("\n\nThe SC_Lib Trainer says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	the enhancer script to enhance (cleanse) the speech in the audio file (SCiVo corpus)
//====================================================================================================================
bool SC_MainTasks::enhancer() {
	bool finalResult = true;
  unsigned long int sceneNr = 1, res;
  long int sceneStart, sceneEnd;
	char *audioFile = NULL, *sceneFile = NULL, *segmentFile = NULL, *speechModelFile = NULL, *noiseModelFile = NULL;
	double videoFrameRate = 0;
  time_t startTime	= time(NULL);
  SC_TweakableParameters* pTweak = NULL;
  SC_Corpus *pCorpus = NULL;
  SC_FeatureHandler *pExtractor = NULL;
  SC_SegmentationHandler *pSegmenter = NULL;
  SC_SignalHandler *pSignalHandler = NULL;
  SC_Enhancement* pEnhance = NULL;
  SC_Signal *pSignal = NULL, *pNewSignal = NULL;
  SV_Data **pFeatures = NULL, *pFeatureHook = NULL, *pFeatureFirst = NULL, *pFeature = NULL;

  printf("\nThis is Enhancer: Speech enhancement via the SC_Lib");
  printf("\n=====================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile);
  pExtractor = new SC_FeatureHandler(this->pTweak);
  pSegmenter = new SC_SegmentationHandler(this->pTweak);
  pSignalHandler = new SC_SignalHandler(this->pTweak, sclib::stGuess);
	printf("done!\n");

  //load the complete audio file to insert the enhanced parts into it
  sceneStart = 0;
  sceneEnd = pCorpus->getGT()->getAudioSampleCount();
  pNewSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
	
  printf("\nPhase 1:\tProcessing Scenes...\n");
  for (unsigned long int y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
    
		pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
    if (sceneStart != -1 && sceneEnd != -1) {
		
      printf("\n  Processing Scene %d...", sceneNr);
			
      //skip scene?
      if ((sceneNr < pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & pTweak->general.sceneSelection) && (pTweak->general.sceneSelection != 0))) {
        printf(" skipped!"); 
        sceneNr++; 
        y = sceneEnd; 
        continue;
      }
      
      //load signal
      printf("\n    Phase I:\tLoading Signal...");
	    pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
      printf("done!");

      //extract features
      // - pFeatures[0]: MFCC's or NULL
      // - pFeatures[1]: Energy&ZCR or NULL
      // - pFeatures[2]: FbE's or NULL
      // - pFeatures[3]: Spectrum or NULL
      printf("\n    Phase II:\tExtracting Features...");
			pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(), pSegmenter->getSpecialFeatureParameters()); //sclib::featureEnergyZCR|sclib::featureFbE|sclib::featureMFCC|SCLIB_FEATURE_SPECTRUM
	    printf("done!");

      //silence detection
      printf("\n    Phase IIIa:\tRemoving Silence...");
	    pSegmenter->silenceDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
	    printf("done!");
    
      //label unvoiced speech
	    printf("\n    Phase IIIb:\tLabeling Unvoiced Speech...");
	    SV_Data *pPitch = pSegmenter->unvoicedDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);
      MFree_0D(pPitch);
	    printf("done!");

      //scd
      printf("\n    Phase IIIc:\tDetecting Speaker Changes...");
      res = pSegmenter->speakerChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.changeDetectorMode);
      printf("done! (%d changes found)", res);
      
      //acd
      printf("\n    Phase IIId:\tDetecting Acoustic Changes...");
	    res = pSegmenter->acousticChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.changeDetectorMode);
      printf("done! (%d changes found)", res);

      //do the enhancement
      printf("\n    Phase IV:\tEnhancing speech...");
      pEnhance = new SC_Enhancement(pCorpus->getGT(), this->pTweak, speechModelFile, noiseModelFile, sceneStart, pSignal->GetBuf_L(), pSignal->GetLen());
      if (!pEnhance->enhance(false)) {
        REPORT_ERROR(0, "Unknown Error during enhancement!");
				finalResult = false;
      }
      MFree_0D(pSignal);

      //implant the enhanced scene into the complete original samples
      if (pSignalHandler->implantSamples(pNewSignal, sceneStart, pEnhance->getEnhancedSignalPointer(), pEnhance->getEnhancedSignalLength()) != true) {
        REPORT_ERROR(0, "Unknown Error during sample-implantation!");
				finalResult = false;
      }
      MFree_0D(pEnhance);
      
      printf("done!\n");
			pExtractor->freeFeatures(pFeatures);

			sceneNr++;
      y = sceneEnd;

      if (sceneNr > pTweak->general.lastScene) {printf("\nSkipping remaining scenes!\n"); break;}
    } else {
			break;	//scene not valid
		} 
  }

  printf("\ndone!\n");
	printf("\nPhase 3:\tCleaning up...");

  pSignalHandler->saveSignal("new.wav", pNewSignal, true);
  MFree_0D(pNewSignal);

  MFree_0D(pSignalHandler);
  MFree_0D(pCorpus);
	MFree_0D(pSegmenter);
  MFree_0D(pExtractor);
	MFree_1D(speechModelFile);
  MFree_1D(noiseModelFile);

	printf("done!");
	printf("\n\nEnhancer says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	the modeller script to build individual models of given type, feature and audioportion
//====================================================================================================================
bool SC_MainTasks::modeller(const char *modelFile, const char *featureFile, const char *normalizationFile, unsigned long int audioTypes, unsigned long int audioTypesNot, unsigned long int featureType, unsigned long int modelType, unsigned short int modelOrder) {
	bool finalResult = true;
  unsigned long int y, sceneNr = 1, res;
  long int sceneStart, sceneEnd;
  time_t startTime	= time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Model *pModel = NULL;
  SC_SegmentationHandler *pSegmenter;
  SC_FeatureHandler *pExtractor;
  SC_ModelHandler *pModeller;
  SC_Signal *pSignal;
  SV_Data **pFeatures, *pFeatureHook = NULL, *pFeatureFirst = NULL, *pFeature = NULL, *pTmpFeatures = NULL, *pNorm = NULL;
	SC_SVMproblem *pProblem;

  printf("\nThis is Modeller: Building explicit modells via the SC_Lib");
  printf("\n===========================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  pSegmenter =	new SC_SegmentationHandler(this->pTweak);
  pModeller = new SC_ModelHandler(this->pTweak);
  pExtractor = new SC_FeatureHandler(this->pTweak);

  //decide which corpus to use:
  char *extension = sclib::extractExtension(this->audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) == 0) { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile);
	} else if (strncmp(extension, "mac", sclib::bufferSize) == 0) { //a MAC corpus file
		pCorpus = new SC_Corpus_MAC(this->pTweak, this->audioFile);
	} else if (strncmp(this->segmentationFile, "", sclib::bufferSize) == 0) { //no gt => MPEG7 corpus
    pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	} else { //no corpus-file, so assume SCiVo type data
    pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile);
  }
  MFree_1D(extension);
	printf("done!\n");

	//load normalization-matrix, if exists
	if (sclib::fileExists(normalizationFile) == true) {
		pNorm = pExtractor->loadFeature(normalizationFile);
	}

  if (sclib::fileExists(featureFile) == false) {
    
		printf("\nPhase 1:\tProcessing Audio-Data piece-wise...\n");
		y = 0;
    while (y < pCorpus->getGT()->getAudioSampleCount()) {
    
			//use different methods to chop audio-data into smaller pieces according to the corpus used: scene-wise or file-wise
			if (pCorpus->getGT()->getGTtype() == sclib::gtMAC) { //file-wise
        printf("\n    Phase I:\tLoading Signal...");
				sceneStart = y;
 	      pSignal = ((SC_Corpus_MAC*)pCorpus)->loadSignal((unsigned long&)(sceneStart), (unsigned long&)(sceneEnd), 0, audioTypes);
				if (pSignal != NULL) {
					printf("done! We've reached %#5.1f%%", (double)(sceneEnd) / (double)(pCorpus->getGT()->getAudioSampleCount()) * 100.0);
				} else {
          printf("failed! No more samples of this type to load.");
          sceneEnd = pCorpus->getGT()->getAudioSampleCount(); //to finish the while-loop
				}
			} else { //scene-wise
				pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
				if (sceneStart == -1 || sceneEnd == -1) {
					break; //scene not valid
				} else {
					printf("\n  Processing Scene %d...", sceneNr);
				}
        //skip scene?
        if ((sceneNr < this->pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & pTweak->general.sceneSelection) && (this->pTweak->general.sceneSelection != 0))) {
          printf(" skipped!"); 
          sceneNr++; 
          y = sceneEnd + 1; 
          continue;
        }
				//load signal
				printf("\n    Phase I:\tLoading Signal...");
				pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd); //sceneStart and -end may get altered here...
				printf("done!");
			}

			if (pSignal != NULL && pSignal->GetLen() > 0) {
				//extract features
				printf("\n    Phase II:\tExtracting Features...");
				pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures()|featureType, pSegmenter->getSpecialFeatureParameters());
				MFree_0D(pSignal);
				printf("done!");

				//silence detection
				printf("\n    Phase IIIa:\tRemoving Silence...");
				pSegmenter->silenceDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
				printf("done!");

				//label unvoiced speech
				printf("\n    Phase IIIb:\tLabeling Unvoiced Speech...");
				SV_Data *pPitch = pSegmenter->unvoicedDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);
				MFree_0D(pPitch);
				printf("done!");
	      
				//scd
				printf("\n    Phase IIIc:\tDetecting Speaker Changes...");
				res = pSegmenter->speakerChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.changeDetectorMode);
				printf("done! (%d changes found)", res);

				//acd
				printf("\n    Phase IIId:\tDetecting Acoustic Changes...");
				res = pSegmenter->acousticChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.changeDetectorMode);
				printf("done! (%d changes found)", res);

				//maybe concatenate multiple feature-sets (a normlization-matrix must be given, then)
				if (sclib::isPowerOfTwo(featureType) == false) {
					pExtractor->prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
					pExtractor->equalizeFrameParameters(pFeatures, pExtractor->getFeatureCount(), featureType);
					pTmpFeatures = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount(), featureType);
					if (pNorm != NULL)  {
						pExtractor->normalize(pTmpFeatures, pNorm);
					}
				} else {
					pTmpFeatures = pFeatures[sclib::bitPosition(featureType)];
				}

				//copy desired frames together and free original feature-sets
				pFeature = pCorpus->getGT()->copyFramesTogether(pTmpFeatures, sceneStart, sceneStart, sceneEnd, audioTypes, audioTypesNot);
				pExtractor->freeFeatures(pFeatures);
				if (sclib::isPowerOfTwo(featureType) == false) {
					MFree_0D(pTmpFeatures); //only free it if it hasn't been done before, i.e. if it is a new construction, not just a pointer into pFeatures
				}

				//build trajectroeis a priori if the time model is choosen so that the loadAndMerge mechanism can be used 
				if (modelType == sclib::mtTime || modelType == this->pTweak->modelTime.subModelType) {
					pTmpFeatures = pExtractor->createTrajectories(pFeature, pCorpus->getGT()->getConverter()->ms2audioFrame(this->pTweak->modelTime.syllableLength, pFeature->Hdr.frameSize, pFeature->Hdr.frameStep, sclib::alignmentStart), this->pTweak->modelTime.trajectoryStep, this->pTweak->modelTime.removeTiming, this->pTweak->modelTime.templateCount, this->pTweak->modelTime.clusteringIterations);
					MFree_0D(pFeature);
					pFeature = pTmpFeatures;
				}
	  
				//build linked list of all desired features of this video
				//pFeatureFirst is the pointer to the first element of the linked list, pFeatureHook to the current last element
				if (pFeature != NULL) {
					if (pFeatureFirst == NULL) {
						pFeatureFirst = pFeature;
						pFeatureHook = pFeatureFirst; 
					} else {
						pFeatureHook->Next = pFeature;
						pFeatureHook = pFeatureHook->Next;
					}
				}
			} //pSignal valid

			//special treatment for chopping-schemes
			if (pCorpus->getGT()->getGTtype() == sclib::gtMAC) { 
				y = sceneEnd + 1;
			} else {
				printf("\n");
				sceneNr++;
				y = sceneEnd + 1;
				if (sceneNr > pTweak->general.lastScene) {
					printf("\nSkipping remaining scenes!\n"); 
					break;
				}
			}

    }
    printf("\ndone!\n");

    printf("\nPhase 2:\tSaving Features...");
    pExtractor->saveFeatureList(featureFile, pFeatureFirst);
    sclib::destructLinkedList(pFeatureFirst);
    pFeatureFirst = NULL;
    pFeatureHook = NULL;
    printf("done!\n");

  } //!fileExists(featureFile)

  printf("\nPhase 3:\tBuilding models...");
  if (sclib::fileExists(featureFile) == true) {
    //save and load data to/from disk, to save one copy in memory...
		if (modelType==sclib::mtSVM && modelOrder>0 && (strncmp(normalizationFile, "", sclib::bufferSize)==0 || sclib::fileExists(normalizationFile)==true)) {
			pProblem = pExtractor->loadAndMergeSvmFeatures(featureFile);
			pFeature = NULL;
		} else {
			pFeature = pExtractor->loadAndMergeFeatures(featureFile);
			pProblem = NULL;
		}

		//normalization-filename given but file doesn't exist => build it!
		if (strncmp(normalizationFile, "", sclib::bufferSize) != 0 && sclib::fileExists(normalizationFile) == false) {
			printf("\nNormalization matrix was needed but not existent => I now build it, then please restart the program and extract features again (delete old unnormalized feature file first!)");
			MFree_0D(pNorm); //should already be NULL normally...
			pNorm = pExtractor->createNormalizationMatrix(pFeature);
			if (pExtractor->saveFeature(normalizationFile, pNorm) != true) {
				printf("\nError while saving normalization matrix...");
				finalResult = false;
			}
		}

		//finally, build the model
		if (modelOrder > 0) { //use given modelOrder
			if (pProblem != NULL) {
				pModel = new SC_Model_SVM(this->pTweak, this->pTweak->modelSvm.distanceBasedTesting, this->pTweak->modelSvm.doParameterSearch, true);
				((SC_Model_SVM*)pModel)->TrainModel(pProblem);
			} else {
				pModel = pModeller->buildModel(pFeature, NULL, modelOrder, modelType, 0);
			}
		} else { //guess model order
			pTweak->modelHandler.foregroundModelType = modelType;
			pModel = pModeller->buildModel(pFeature, NULL, sclib::modeForeground, 0);
		}
		res = pModeller->saveModel(modelFile, pModel);
		if (res <= 0) {
			printf("\nSorry, the model could not be saved\n");
			finalResult = false;
		}
		MFree_0D(pModel); //kills also the SVM-Problem, because it is its training data
    MFree_0D(pFeature);
		//SC_Classifier_SVM::killSVMproblem(pProblem);
		//MFree_0D(pProblem);

	  printf("done!\n");
  } else {
    printf("no data!\n");
		finalResult = false;
  }

	printf("\nPhase 3:\tCleaning up...");

	MFree_0D(pNorm);
  MFree_0D(pCorpus);
	MFree_0D(pSegmenter);
  MFree_0D(pModeller);
  MFree_0D(pExtractor);

	printf("done!");
	printf("\n\nModeller says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	the combineModels script from modeller
//====================================================================================================================
bool SC_MainTasks::combineModels(const char *modelFile1, const char *modelFile2, const char *combinedModelFile, unsigned long int modelType) {
	bool finalResult = true;
	unsigned long int oldC1, oldC2;
	int res;
	SC_Model *p1, *p2, *pC;
	SC_ModelHandler *pModeller = new SC_ModelHandler(this->pTweak, true);

	//pTweak->mixtureModelGmm.weightLimit = -1.0;
	p1 = pModeller->loadModel(modelFile1, modelType);
  oldC1 = p1->getTrainingDataCount();
  p1->setTrainindDataCount(1); //to give both parts equal weights regardless their training-data

	p2 = pModeller->loadModel(modelFile2, modelType);
  oldC2 = p2->getTrainingDataCount();
  p2->setTrainindDataCount(1); //to give both parts equal weights regardless their training-data

	pC = pModeller->combineModels(p1, p2);
  pC->setTrainindDataCount(oldC1 + oldC2);
  MFree_0D(p1);
  MFree_0D(p2);

	res = pModeller->saveModel(combinedModelFile, pC);
  MFree_0D(pC);

	MFree_0D(pModeller);

	finalResult = (res > 0);
	return finalResult;
}

//====================================================================================================================
//	the statistician script; prints densities of speaker duration distributions and detectability of change points
//====================================================================================================================
bool SC_MainTasks::statistician() {
	bool finalResult = true;
	char *extension;
	SC_Corpus *pCorpus;
	SC_Segmentation_Changes_LZW *pCD;
	SC_Model_Pareto *pModel;

  printf("\nThis is the Statisticion: Information about the ground truth of a media file");
  printf("\n=============================================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
	
  //decide which corpus to use:
	extension = sclib::extractExtension(audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) != 0) { //no corpus-file, so assume SCiVo or MPEG7 type data
		if (strncmp(this->segmentationFile, "", sclib::bufferSize) != 0) {
			pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile); //SCiVo-files do have ground-truth
		} else {
			pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile); //...MPEG7-files don't (maybe cutlists...)
		}
  } else { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile);
  }
	MFree_1D(extension);

	//check groundtruth
	if (true != pCorpus->getGT()->checkConsistency(0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::modeGroundtruth)) {
		printf("%s", "Consistency errors exist in the ground truth -> see transcript!");
	}
	printf("done!\n");

	printf("\nPhase 1:\tInfering segment-durations and detectability...");
	pCD = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeSpeakerChange); //which detector to take doesn't matter, the methods belong to the base class anyway...
	pModel = pCD->getSegmentDurationDistribution(pCorpus->getGT(), 0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atSpeech);
	sclib::classOut("durations.txt", pModel, pTweak);
	MFree_0D(pModel);
	pModel = pCD->getDetectability(pCorpus->getGT(), 0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atSpeech);
	sclib::classOut("detectability.txt", pModel, pTweak);
	MFree_0D(pModel);
	printf("done!\n");

	printf("\nPhase 4:\tCleaning up...");
  delete pCorpus;
	printf("done!");
	printf("\n\nThe SC_Lib Testsuite says: \"Good bye\".\n");

	return finalResult;
}

//====================================================================================================================
//	Bing's testscript for the v/uv stuff
//====================================================================================================================
bool SC_MainTasks::vuvDetection() {
	bool finalResult = true;
  unsigned long int sceneNr = 1;
  long int sceneStart, sceneEnd;
  int segmentLength = 0;
  time_t startTime	= time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Signal *pSignal = NULL;
  SV_Data **pFeatures = NULL;
  SC_FeatureHandler *pExtractor = NULL;
  SC_ModelHandler *pModeller = NULL;
  SC_SegmentationHandler *pSegmenter = NULL;
  SC_SignalHandler *pSigHandler = NULL;

  printf("\nThis is VUVDetection: Speaker Classification via the SC_Lib");
  printf("\n============================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  pModeller = new SC_ModelHandler(this->pTweak);
  pExtractor = new SC_FeatureHandler(this->pTweak);
  pSegmenter = new SC_SegmentationHandler(this->pTweak);
  pSigHandler = new SC_SignalHandler(this->pTweak);
  
  //decide which corpus to use:
	char *extension = sclib::extractExtension(audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) != 0) { //no corpus-file, so assume VUVDetection type data
   pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile);
  } else { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile);
  }
	MFree_1D(extension);
	printf("done!\n");
  printf("\nPhase 1:\tProcessing Scenes...\n");
	 
	// run the alg
	for (unsigned long int y = 0; y < pCorpus->getGT()->getAudioSampleCount(); y++) {
		pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
		if (sceneStart != -1 && sceneEnd != -1) {
			
			printf("\n  Processing Scene %d...", sceneNr);
			//skip scene?
			if ((sceneNr < this->pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & this->pTweak->general.sceneSelection) && (this->pTweak->general.sceneSelection != 0))) {
				printf(" skipped!"); 
				sceneNr++; 
				y = sceneEnd; 
				continue;
			}
			
			//load signal
			printf("\n    Phase Ia:\tLoading Signal...");
 			pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
			printf("done!");

			//extract features
			printf("\n    Phase II:\tExtracting Features...");
			pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(this->pTweak->segmentationHandler.vUvDetectorMode), pSegmenter->getSpecialFeatureParameters(this->pTweak->segmentationHandler.vUvDetectorMode), true);
			MFree_0D(pSignal);
			printf("done!");

			//label unvoiced speech
			printf("\n Phase III: Labeling Unvoiced Speech...");
			SV_Data *pPitch = pSegmenter->unvoicedDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);
			MFree_0D(pPitch);
			printf("done!");
		 
			sceneNr++;
			y = sceneEnd;
			pExtractor->freeFeatures(pFeatures);
			if (sceneNr > this->pTweak->general.lastScene) {
				printf("\nSkipping remaining scenes!\n"); 
				break;
			}

		} else {
			break;	//scene not valid
		} 
	}
  printf("\ndone!\n");
	
	//Scoring
	printf("\n Phase IV: Scoring for V/Uv-Segmentation!\n");
	//pCorpus->getGT()->frameListOut("f.txt",0,pCorpus->getGT()->getAudioSampleCount()-1);
	SC_Score_VUv  *pScoreVUv = new SC_Score_VUv(this->pTweak, pCorpus->getGT());
  pScoreVUv->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1);
  pScoreVUv->printReport("ReportVUv.txt", true, (unsigned long)(time(NULL) - startTime));

	//pCorpus->getGT()->frameListOut("frames.txt", 0, pCorpus->getGT()->getAudioSampleCount()-1);

  printf("done!");

	//exit
	printf("\nPhase 3:\tCleaning up...");
	MFree_0D(pScoreVUv);
  MFree_0D(pCorpus);
  MFree_0D(pModeller);
  MFree_0D(pExtractor);
  MFree_0D(pSegmenter);
	MFree_0D(pSigHandler);

	printf("done!");
	printf("\n\nThe SC_Lib Testsuite says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	the Wesley script
//====================================================================================================================
bool SC_MainTasks::wesley(long int whichExperiment, unsigned long int noiseRobustModelType, char *rawDataDir) {
	bool finalResult = true;
  unsigned long int sampleRate;
  SC_Corpus_Wesley *pCorpus = NULL;
	SC_Cluster *pTrainData = NULL;
	char tmp1[sclib::bufferSize], tmp2[sclib::bufferSize];

	//some constants to choose which experiment to do
	const int SCLIB_ME_MFC_GMM = 1;
	const int SCLIB_ME_MFC_SGMM = 2;
	const int SCLIB_ME_FBK_GMM = 3;
	const int SCLIB_ME_FBK_SGMM = 4;
	const int SCLIB_WES_MFC_GMM = 5;
	const int SCLIB_WES_MFC_SGMM = 6;
	const int SCLIB_WES_FBK_GMM = 7;
	const int SCLIB_WES_FBK_SGMM = 8;
	const int SCLIB_NOIZ_MFC_GMM = 9;
	const int SCLIB_NOIZ_MFC_SGMM = 10;
	const int SCLIB_NOIZ_FBK_GMM = 11;
	const int SCLIB_NOIZ_FBK_SGMM = 12;
	//const int SCLIB_WES_MFC_WESMOD_SGMM = 9;

  printf("\nThis is a Testsuit for Wesleys HTK-files with the SC_Lib algorithms");
  printf("\n====================================================================\n");
  
  printf("\nloading parameters...");
	printf("\n    Experiment: %d, Voice-Model: %d, Raw-Data: %s", whichExperiment, noiseRobustModelType, rawDataDir);
	//pTweak->modelHandler.maxSpeakerModelOrder = 40;
	//pTweak->modelHandler.maxSpeakerModelOrder = 32;
	//pTweak->debug.debugMode |= SCLIB_DB_SPEAKER_MODELS;
	switch (whichExperiment) { //wesley's files are 22kHz, mine are 16kHz, NOIZEUS is 8kHz; also use filesystem-knowledge to set the directory for model-storing correctly here
		case 1: case 2: case 3: case 4: 
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/gmm");
			this->pTweak->setDebugDir(tmp1);
			sampleRate = 16000; 
			break;
		case 5: case 6: case 7: case 8: 
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/gmm");
			this->pTweak->setDebugDir(tmp1);
			sampleRate = 22000; 
			break;
		case 9: case 10: case 11: case 12: 
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/gmm");
			this->pTweak->setDebugDir(tmp1);
			sampleRate = 8000; 
			break;
	}
  pCorpus = new SC_Corpus_Wesley(this->pTweak, sampleRate, false, sclib::modeSClib);

  printf("\ndone!\n");
  
  switch (whichExperiment) {
    case SCLIB_ME_FBK_GMM: {
      printf("starting in mode: SCLIB_ME_FBK_GMM\n");
      this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
      //build voice-models
      printf("build & save voice models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_me_fbk_gmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_ME_FBK_SGMM: {
      printf("starting in mode: SCLIB_ME_FBK_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = noiseRobustModelType;
      //build noise-models needed for training with MIXMAX models
      printf("build noise-models needed for training with MIXMAX models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_train-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild voice-models using the prebuild noise-models if MIXMAX models are wished...");

      //build voice-models using the prebuild noise-models if MIXMAX models are wished
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild noise-models needed for testing with MIXMAX models...");

      //build noise-models needed for testing with MIXMAX models
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_test-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_test-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_fbk_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_me_fbk_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_ME_MFC_GMM: {
      printf("starting in mode: SCLIB_ME_MFC_GMM\n");
      this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
      //build voice-models
      printf("build & save voice models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_me_mfc_gmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_ME_MFC_SGMM: {
      printf("starting in mode: SCLIB_ME_MFC_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = noiseRobustModelType;
      //build noise-models needed for training with MIXMAX models
      printf("build noise-models needed for training with MIXMAX models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_train-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild voice-models using the prebuild noise-models if MIXMAX models are wished...");

      //build voice-models using the prebuild noise-models if MIXMAX models are wished
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild noise-models needed for testing with MIXMAX models...");

      //build noise-models needed for testing with MIXMAX models
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_test-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_test-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/me_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/me_mfc_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_me_mfc_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_WES_FBK_GMM: {
      printf("starting in mode: SCLIB_WES_FBK_GMM\n");
      this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
      //build voice-models
      printf("build & save voice models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_wes_fbk_gmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_WES_FBK_SGMM: {
      printf("starting in mode: SCLIB_WES_FBK_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = noiseRobustModelType;
      //build noise-models needed for training with MIXMAX models
      printf("build noise-models needed for training with MIXMAX models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild voice-models using the prebuild noise-models if MIXMAX models are wished...");

      //build voice-models using the prebuild noise-models if MIXMAX models are wished
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild noise-models needed for testing with MIXMAX models...");

      //build noise-models needed for testing with MIXMAX models
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_test-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_test-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_fbk_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_wes_fbk_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_WES_MFC_GMM: {
      printf("starting in mode: SCLIB_WES_MFC_GMM\n");
      this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
      //build voice-models
      printf("build & save voice models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_wes_mfc_gmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_WES_MFC_SGMM: {
      printf("starting in mode: SCLIB_WES_MFC_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = noiseRobustModelType;
      //build noise-models needed for training with MIXMAX models
      printf("build noise-models needed for training with MIXMAX models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild voice-models using the prebuild noise-models if MIXMAX models are wished...");

      //build voice-models using the prebuild noise-models if MIXMAX models are wished
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild noise-models needed for testing with MIXMAX models...");

      //build noise-models needed for testing with MIXMAX models
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_test-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_test-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "wesley_corpus/wes_mfc_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_wes_mfc_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_NOIZ_FBK_GMM: {
      printf("starting in mode: SCLIB_NOIZ_FBK_GMM\n");
      this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
      //build voice-models
      printf("build & save voice models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_noiz_fbk_gmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_NOIZ_FBK_SGMM: {
      printf("starting in mode: SCLIB_NOIZ_FBK_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = noiseRobustModelType;
      //build noise-models needed for training with MIXMAX models
      printf("build noise-models needed for training with MIXMAX models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_train-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild voice-models using the prebuild noise-models if MIXMAX models are wished...");

      //build voice-models using the prebuild noise-models if MIXMAX models are wished
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild noise-models needed for testing with MIXMAX models...");

      //build noise-models needed for testing with MIXMAX models
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/test-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_test-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/fbk_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_noiz_fbk_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_NOIZ_MFC_GMM: {
      printf("starting in mode: SCLIB_NOIZ_MFC_GMM\n");
      this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;
      //build voice-models
      printf("build & save voice models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_noiz_mfc_gmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    case SCLIB_NOIZ_MFC_SGMM: {
      printf("starting in mode: SCLIB_NOIZ_MFC_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = noiseRobustModelType;
      //build noise-models needed for training with MIXMAX models
      printf("build noise-models needed for training with MIXMAX models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_train-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild voice-models using the prebuild noise-models if MIXMAX models are wished...");

      //build voice-models using the prebuild noise-models if MIXMAX models are wished
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_train-vocal.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeForeground, "-v.gmm");
      //res = pCorpus->saveModels(pTrainData, "-v.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\nbuild noise-models needed for testing with MIXMAX models...");

      //build noise-models needed for testing with MIXMAX models
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/test-noise.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_test-noise.crp");
      pTrainData = pCorpus->trainModels(tmp1, tmp2, sclib::modeBackground, "-m.gmm");
      //res = pCorpus->saveModels(pTrainData, "-m.gmm");
      pTrainData->killAllLinkedData(pTrainData);
      printf("done!\ntest the testcorpus against the prebuild models...");

      //test the testcorpus against the prebuild models
      printf("\n load models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/train-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_train-vocal.crp");
      pTrainData = pCorpus->loadModels(tmp1, tmp2, sclib::modeForeground);
      printf("done!\n test models...");
			sprintf(tmp1, "%s%s\0", rawDataDir, "NOIZEUS/test-vocal.lst");
			sprintf(tmp2, "%s%s\0", rawDataDir, "NOIZEUS/mfc_test-vocal.crp");
      pCorpus->testModels(tmp1, tmp2, pTrainData, "result_noiz_mfc_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }

    /*
    case SCLIB_WES_MFC_WESMOD_SGMM: {
      printf("starting in mode: SCLIB_WES_MFC_WESMOD_SGMM\n");
      this->pTweak->modelHandler.foregroundModelType = SCLIB_MT_GMM_IB;
      //pCorpus->setModelOrigin(SCLIB_MT_WESLEY);

      //build voice models for use in weselys system by using his noise-models
      //has nothing to do with the below training, so comment one of both blocks!!!
      if (argc > 1) {
        pTrainData = pCorpus->trainModels(argv[2], argv[3], SCLIB_AT_SPEECH);
      } else {
        pTrainData = pCorpus->trainModels("../../../Rohdaten/wesley_corpus/wes_mfc_train_male-vocal.lst", "../../../Rohdaten/wesley_corpus/wes_mfc_train-vocal.crp", sclib::modeForeground);
      }
      res = pCorpus->saveModels(pTrainData, ".32-8");
      pTrainData->killAllLinkedData(pTrainData);
      
      //test the testcorpus against the prebuild models
      printf("\n load Wesley's models...");
      pTrainData = pCorpus->loadModels("../../../Rohdaten/wesley_corpus/wes_mfc_WESMOD_train-vocal.lst", "../../../Rohdaten/wesley_corpus/wes_mfc_train-vocal.crp", sclib::modeForeground);
      printf("done!\n test models...");
      pCorpus->testModels("../../../Rohdaten/wesley_corpus/wes_mfc_WESMOD_test-vocal.lst", "../../../Rohdaten/wesley_corpus/wes_mfc_test-vocal.crp", pTrainData, "result_wes_mfc_WESMOD_sgmm.txt");
	    pTrainData->killAllLinkedData(pTrainData);
      break;
    }
    */
  }

  printf("done!\n");
  printf("cleaning up...");

	MFree_0D(pCorpus);

	printf("done!");
	printf("\n\nThe Wesley's Testsuite says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	the SCiVo script
//====================================================================================================================
bool SC_MainTasks::scivo(char *explicitModelFile, char ***explicitModels) {
	bool finalResult = true;
  unsigned long int sceneNr = 1, atcUncertainty, scdUncertainty, vuvUncertainty, enhancedSignalLength, speakerModelFeature;
  long int sceneStart, sceneEnd, res;
  int segmentLength = 0;
	char previousResults[sclib::bufferSize];
	double videoFrameRate = 0;
	unsigned short int explicitModelCount = 0;
	short *enhancedSignalPointer;
	char *extension;
  time_t startTime = time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Signal *pSignal = NULL;
  SV_Data **pFeatures = NULL, *pPitch = NULL, *pConcat; //, *pNorm;
	SC_Cluster *pAllClusters = NULL, *pClusters	= NULL, *pAllShortClusters = NULL, *pShortClusters = NULL, *pHook = NULL;
  SC_SpeakerScore_Clustering *pScore = NULL;
  SC_SpeakerScore_Identification *pScoreID = NULL;
	SC_SpeakerScore_Classification *pFinalScore = NULL;
	SC_Score_AudioTypeClassification *pScoreATC = NULL;
  SC_Score_ChangeDetection *pScoreSCD = NULL;
	SC_Score_VUv *pScoreVUv = NULL;
  SC_FeatureHandler *pExtractor = NULL;
  SC_ModelHandler *pModeller = NULL;
  SC_SegmentationHandler *pSegmenter = NULL;
  SC_SpeakerClusterer *pClusterer = NULL;
  SC_SpeakerIdentificator *pIdentificator = NULL;
  SC_Model *pUBM = NULL;
  SC_SignalHandler *pSigHandler = NULL;
	SC_Enhancement *pEnhancer = NULL;

	printf("\nThis is SCiVo: Speaker Classification via the SC_Lib");
  printf("\n=====================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
	pModeller = new SC_ModelHandler(this->pTweak);
  pExtractor = new SC_FeatureHandler(this->pTweak);
  pSegmenter = new SC_SegmentationHandler(this->pTweak);
  pSigHandler = new SC_SignalHandler(this->pTweak);
  //decide which corpus to use:
	extension = sclib::extractExtension(audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) != 0) { //no corpus-file, so assume SCiVo or MPEG7 type data
		if (strncmp(this->segmentationFile, "-", sclib::bufferSize) != 0) {
			pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile); //SCiVo-files do have ground-truth
		} else {
			pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile); //...MPEG7-files don't (maybe cutlists...)
		}
  } else { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, true, false, 0, true); //TODO //true, false, false, 0); //
  }
	pScoreATC = new SC_Score_AudioTypeClassification(this->pTweak, pCorpus->getGT());
	pScoreSCD = new SC_Score_ChangeDetection(this->pTweak, pCorpus->getGT());
	pScoreVUv = new SC_Score_VUv(this->pTweak, pCorpus->getGT());
	if (this->pTweak->enhancement.doEnhancement == true) {
		pEnhancer = new SC_Enhancement(pCorpus->getGT(), this->pTweak, this->pTweak->enhancement.speechModelFile);
	}
	MFree_1D(extension);
  explicitModelCount = pCorpus->getGT()->readExplicitModelList(explicitModelFile, explicitModels);
  if (this->pTweak->modelHandler.onlyThisSpeaker != NULL && strlen(this->pTweak->modelHandler.onlyThisSpeaker) > 0) {
    this->pTweak->setDebugPrefix(this->pTweak->modelHandler.onlyThisSpeaker);
  }
	if (true != pCorpus->getGT()->checkConsistency(0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::modeGroundtruth)) {
		printf("%s", "Consistency errors exist in the ground truth -> see transcript!");
	}
	atcUncertainty = pCorpus->getGT()->getConverter()->ms2sample(pSegmenter->getUncertaintyRegionWidth(this->pTweak->segmentationHandler.audioTypeMode));
	vuvUncertainty = pCorpus->getGT()->getConverter()->ms2sample(pSegmenter->getUncertaintyRegionWidth(this->pTweak->segmentationHandler.vUvDetectorMode));
	scdUncertainty = pCorpus->getGT()->getConverter()->ms2sample(pSegmenter->getUncertaintyRegionWidth(this->pTweak->segmentationHandler.changeDetectorMode));
	printf("done!\n");

	//try to load previous segmentation results
	if (this->pTweak->general.preClusteringResultsPrefix != NULL && strncmp(this->pTweak->general.preClusteringResultsPrefix, "", 1) != 0) {
		sprintf(previousResults, "%s%s", this->pTweak->general.preClusteringResultsPrefix, ".gt");
	}
	if (sclib::fileExists(previousResults) == true) { 
	  printf("\nPhase 1:\tLoading previous processing results...\n");
		if (pCorpus->getGT()->load(previousResults) == NULL) {
			REPORT_ERROR(SVLIB_FileErr, "Error loading ground truth class!!!");
		}
		sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, "_speakers.cluster_1");
		if (sclib::fileExists(previousResults)) {
			pAllClusters = new SC_Cluster(this->pTweak);
			sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, "_speakers");
			if (pAllClusters->load(previousResults) == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "Error loading speaker clusters!!!");
			}
		}
		sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, "_short.cluster_1");
		if (sclib::fileExists(previousResults)) {
			pShortClusters = new SC_Cluster(this->pTweak);
			sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, "_short");
			if (pShortClusters->load(previousResults) == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "Error loading short clusters!!!");
			}
		}
		printf("done!\n");
	} 
	//calculate new segmentation results
	else { 
		//TODO: make feature selection more sophisticated
		speakerModelFeature = (this->pTweak->modelHandler.foregroundModelType == sclib::mtMIXMAX) ? sclib::featureFbE : this->pTweak->modelHandler.speakerModelFeature; //MIXMAX models only work on FBE features

		printf("\nPhase 1:\tProcessing Scenes...\n");
		for (unsigned long int y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
	    
			pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
			if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
			
				printf("\n  Processing Scene %d...", sceneNr);
			
				//skip scene?
				if ((sceneNr < this->pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & this->pTweak->general.sceneSelection) && (this->pTweak->general.sceneSelection != 0))) {
					printf(" skipped!"); 
					sceneNr++; 
					y = sceneEnd; 
					continue;
				}

				//load previous work if wished and existant
				sprintf(previousResults, "%s_%d%s\0", this->pTweak->general.featurePrefix, sceneNr, ".dat");
				if (sclib::fileExists(previousResults) == true) {
					pFeatures = pExtractor->loadFeatures(previousResults);
				} else {

					//load signal
					printf("\n    Phase Ia:\tLoading Signal...");
 					pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
					printf("done!");

					//extract features
					printf("\n    Phase II:\tExtracting Features...");
					if (this->pTweak->enhancement.doEnhancement == true) { //in case of speech enhancement to be done, these features are only needed for non-speech-specific segmentation and are extracted once more on the enhanced signal later on
						pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(sclib::algorithm_nothing, 2, -1));
					} else { //if no enhancement is to be carried out, we don't need the signal any more
						pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures()|speakerModelFeature, pSegmenter->getSpecialFeatureParameters());
						MFree_0D(pSignal); 
					}
					printf("done!");
		
					//save previous work for future runs if wished
					if (this->pTweak->enhancement.doEnhancement==false && this->pTweak->general.featurePrefix!=NULL && strncmp(this->pTweak->general.featurePrefix, "", 1)!=0) {
						sprintf(previousResults, "%s_%d%s\0", this->pTweak->general.featurePrefix, sceneNr, ".dat");
						pExtractor->saveFeatures(previousResults, pFeatures);
					}
				} //nothing to load, extract features from scratch

				//silence detection
				printf("\n    Phase IIIa:\tRemoving Silence...");
				pSegmenter->silenceDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
				printf("done!");

				//audio type classification
				printf("\n    Phase IIIb:\tAudio Type Classification...");
				pSegmenter->audioClassification(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.audioTypeMode);
				//pScoreATC->calcScores(sceneStart, sceneEnd, atcUncertainty);
				//pScoreATC->printReport("atc_score_per_scene.txt", false, 0, false, true);
				printf("done!");

				//acd
				printf("\n    Phase IIIc:\tDetecting Acoustic Changes...");
				//TODO: comment out for scd tests
				res = pSegmenter->acousticChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.changeDetectorMode);
				printf("done! (%d changes found)", res);

				printf("\n    Phase IIId:\tDoing Speech Enhancement...");
				if (this->pTweak->enhancement.doEnhancement == true) {
					pEnhancer->setSegmentStart(sceneStart);
					pEnhancer->setOriginalSignalPointer(pSignal->GetBuf_L(), pSignal->GetLen());
					if (pEnhancer->enhance(false) == true) {
						pEnhancer->getEnhancedSignalPointer(enhancedSignalPointer, enhancedSignalLength);
						pEnhancer->forgetEnhancedSignalPointer();
						pSigHandler->exchangeSignalBuffer(pSignal, enhancedSignalPointer, enhancedSignalLength);
						pExtractor->freeFeatures(pFeatures); //re-extract features for standard algorithms on the enhanced signal for the now following speech-specific algorithms
						pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(sclib::algorithm_nothing, -1, 3)|speakerModelFeature, pSegmenter->getSpecialFeatureParameters());
						if (this->pTweak->general.featurePrefix != NULL && strncmp(this->pTweak->general.featurePrefix, "", 1) != 0) {
							sprintf(previousResults, "%s_%d%s\0", this->pTweak->general.featurePrefix, sceneNr, ".dat");
							pExtractor->saveFeatures(previousResults, pFeatures);
						}
						//pSignal->SaveSignal("enhanced.wav");
					} else {
						REPORT_ERROR(0, "Unknown Error during enhancement!");
					}
					MFree_0D(pSignal);
					printf("done!");
				} else {
					printf("skipped!");
				}

				//label unvoiced speech
				printf("\n    Phase IIIe:\tLabeling Unvoiced Speech...");
				pPitch = pSegmenter->unvoicedDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);
				if (this->pTweak->featurePitch.method == sclib::modeAS && (sclib::bitTest(speakerModelFeature, sclib::featurePitch) == true || sclib::bitTest(pSegmenter->getUsedFeatures(), sclib::featurePitch) == true)) {
					sclib::destructLinkedList(pFeatures[sclib::bitPosition(sclib::featurePitch)]);
					double invalidValue[1] = {0.0};
					SV_Data *pConvertedPitch = pExtractor->convertFrameRate(pPitch, sceneEnd-sceneStart+1, pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featurePitch.frameSize), pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featurePitch.frameStep), invalidValue);
					pFeatures[sclib::bitPosition(sclib::featurePitch)] = pConvertedPitch;
					MFree_0D(pPitch);
				} else {
					MFree_0D(pPitch);
				}
				if (pCorpus->getGT()->getGTtype() == sclib::gtMPEG7) { //do it here to get new pitch information...
					((SC_Corpus_MPEG7*)pCorpus)->lowLevelFeaturesOut("llFeatures.csv", sceneStart, sceneEnd, pFeatures);
				}
				//pScoreVUv->calcScores(sceneStart, sceneEnd, vuvUncertainty);
				//pScoreVUv->printReport("vuv_score_per_scene.txt", false, 0, false, true);
				printf("done!");

				//artificially shorten all speech segments? 
				if (segmentLength > 0) {
					pCorpus->getGT()->setSpeechSegLength(sceneStart, sceneEnd, pCorpus->getGT()->getConverter()->ms2sample(segmentLength, sclib::alignmentEnd), pCorpus->getGT()->getConverter()->ms2sample(segmentLength, sclib::alignmentEnd));
				}

				//pCorpus->speakerFeaturesOut("speakers.txt", sceneStart, sceneEnd, pFeatures, sclib::featureMFCC|SCLIB_FEATURE_LSP|sclib::featurePitch);

				/*
				SC_Segmentation_Changes_LZW *pSCD = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeSpeakerChange);
				pSCD->trainClassifier(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures);
				exit(0);
				*/

				//TODO: test feature normalization
				/*printf("\n    Phase xxx:\tNormalizing features...");
				for (unsigned int f = 0; f < pExtractor->getFeatureCount(); f++) {
					if (pFeatures[f] != NULL) {
						SV_Data *pNorm = pExtractor->createNormalizationMatrix(pFeatures[f]);
						pExtractor->normalize(pFeatures[f], pNorm);
						MFree_0D(pNorm);
						printf(".");
					}
				}*/

				//scd
				printf("\n    Phase IIIf:\tDetecting Speaker Changes...");
				res = pSegmenter->speakerChangeDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.changeDetectorMode);
				//pScoreSCD->calcScores(sceneStart, sceneEnd, scdUncertainty);
				//pScoreSCD->printReport("scd_score_per_scene.txt", false, 0, false, true);
				//pCorpus->getGT()->frameListOut("frames.txt", sceneStart, sceneEnd);
				printf("done! (%d changes found)", res);

				printf("\n    Phase IV:\tBuilding Speaker-Models...");
				pExtractor->prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
				pExtractor->equalizeFrameParameters(pFeatures, pExtractor->getFeatureCount(), speakerModelFeature);
				pConcat = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount(), speakerModelFeature);
				pExtractor->freeFeatures(pFeatures); //free features
				//if (pCorpus->getGT()->getGTtype() == sclib::gtTIMIT && this->pTweak->modelHandler.foregroundModelType == sclib::mtGroup) {
				//	((SC_GroundTruth_TIMIT*)pCorpus->getGT())->addPhoneLabels(pConcat, sceneStart, sceneEnd, sclib::modeHypothesized);
				//}
				//if (sclib::isPowerOfTwo(speakerModelFeature) == false) {
				//	pNorm = pExtractor->createNormalizationMatrix(pConcat);
				//	pExtractor->normalize(pConcat, pNorm);
				//	MFree_0D(pNorm);
				//}	
				res = pModeller->buildSpeakerModels(pCorpus, sceneStart, sceneEnd, pConcat, pClusters, pShortClusters, explicitModels, explicitModelCount);
				MFree_0D(pConcat);
				
				//link normal clusters together
				pHook = sclib::getLastInList(pAllClusters);
				if (pHook != NULL) {
					pHook->Next	= pClusters;
				} else {
					pAllClusters = pClusters;
				}

				//link short "clusters" together (they are not really clusters because they lack a speaker model [in fact the segments therein where too short to build one] but just hold the features and segment-borders)
				pHook = sclib::getLastInList(pAllShortClusters);
				if (pHook != NULL) {
					pHook->Next	= pShortClusters;
				} else {
					pAllShortClusters = pShortClusters;
				}

				printf("done! (%d models built, %d short segments)\n", res, sclib::getListCount(pShortClusters));
				sceneNr++;
				y = sceneEnd;

				if (true != pCorpus->getGT()->checkConsistency(sceneStart, sceneEnd, sclib::modeHypothesized)) {
					printf("%s", "Consistency errors exist in the algorithmic results -> see transcript!");
				}

				if (sceneNr > this->pTweak->general.lastScene) {printf("\nSkipping remaining scenes!\n"); break;}
			} else {
				break;	//scene not valid
			} 
		}

		//try saving segmentation results
		if (this->pTweak->general.preClusteringResultsPrefix != NULL && strncmp(this->pTweak->general.preClusteringResultsPrefix, "", 1) != 0) {
			sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, ".gt");
			if (pCorpus->getGT()->save(previousResults) != true) {
				REPORT_ERROR(SVLIB_FileErr, "Error saving ground truth class!");
				finalResult = false;
			}
			sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, "_speakers");
			if (pAllClusters != NULL && pAllClusters->save(previousResults) != true) {
				REPORT_ERROR(SVLIB_FileErr, "Error saving speaker clusters!");
				finalResult = false;
			}
			sprintf(previousResults, "%s%s\0", this->pTweak->general.preClusteringResultsPrefix, "_short");
			if (pShortClusters != NULL && pShortClusters->save(previousResults) != true) {
				REPORT_ERROR(SVLIB_FileErr, "Error saving short clusters!");
				finalResult = false;
			}
		}

		printf("\ndone!\n");
	}

	printf("\nPhase 1b:\tScoring ChangePoints before Clustering...");
	if (this->pTweak->segmentationHandler.changeDetectorMode != sclib::algorithm_nothing && pCorpus->getGT()->getGTtype() != sclib::gtStandard && pCorpus->getGT()->getGTtype() != sclib::gtWesley && pCorpus->getGT()->getGTtype() != sclib::gtMPEG7) {
		pScoreSCD->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1, scdUncertainty);
		pScoreSCD->printReport("scd_score_overall.txt", true, (unsigned long)(time(NULL) - startTime));
		printf("done!\n");
	} else {
		printf("skipped!\n");
	}

  if (pAllClusters != NULL && this->pTweak->speakerClusterer.doClustering == true) {
    printf("\nPhase 2a:\tClustering Speakers...");
		//after this, pAllClusters is linked in the first partition of pScore->pPartitionList; it get's destructed when pScore gets destructed
		//pShortClusters is again a linked list of unhandled (due to short duration) speeech segments
    pClusterer = new SC_SpeakerClusterer(this->pTweak);
		
		//TODO!!!
		//double distr[] = {0.0, 40.0/121.0, 38.0/121.0, 22.0/121.0, 15.0/121.0, 4.0/121.0, 0.0/121.0, 0.0/121.0, 1.0/121.0, 0.0/121.0, 1.0/121.0};
		//double distr[] = {0.0, 31.0/128.0, 62.0/128.0, 21.0/128.0, 10.0/128.0, 2.0/128.0, 0.0/128.0, 2.0/128.0, 0.0/128.0, 0.0/128.0, 0.0/128.0};
		//double distr[] = {0.0, 11.0/143.0, 130.0/143.0, 1.0/143.0, 1.0/143.0, 0.0/143.0, 0.0/143.0, 0.0/143.0, 0.0/143.0, 0.0/143.0, 0.0/143.0};
		//pClusterer->randomClustering("random.txt", pCorpus->getGT(), pAllClusters, distr, 10, 10000);
		//end TODO
		
		pScore = pClusterer->clusterSpeakers(pCorpus->getGT(), pAllClusters, pShortClusters); 
		if (pAllShortClusters == NULL) {
			pAllShortClusters = pShortClusters;
		} else {
			if (pShortClusters != NULL) { //set the new short clusters on top of the complete list because it may contain the bigger segments
				pHook = sclib::getLastInList(pShortClusters);
				pHook->Next = pAllShortClusters;
				pAllShortClusters = pShortClusters;
			}
		}
    pScore->printReport("clustering.txt", true, (unsigned long)(time(NULL) - startTime));
	  MFree_0D(pClusterer);
	  printf("done! (%d distinct speakers found)\n", pScore->getClusterCount());

		printf("\nPhase 2b:\tIdentifying Speakers of small segments...");
    if (pAllShortClusters != NULL && this->pTweak->speakerIdentification.doIdentification == true) { //try speaker-identification for those segments too short to build models and/or clusters
      if (pScore->getFinalPartition() != NULL) {
				pUBM = (this->pTweak->speakerIdentification.useUBMs == true) ? pModeller->loadModel(this->pTweak->mixtureModelGmmubm.ubmFileName, sclib::mtGMM_new) : NULL;
        pIdentificator = new SC_SpeakerIdentificator(this->pTweak, (SC_MixtureModel*)pUBM, true);
				pScoreID = pIdentificator->identifySpeakers(pCorpus->getGT(), pAllShortClusters, pScore->getFinalPartition()->getClusters(), std::numeric_limits<double>::max()*-1.0, true, false); //TODO: find proper impostor threshold 
				pScoreID->printReport("id.txt", true, (unsigned long)(time(NULL) - startTime));
        MFree_0D(pUBM);
        MFree_0D(pIdentificator);
				this->pTweak->speakerClusterer.speechSegLengthThreshold = this->pTweak->general.shortSpeechThreshold; //now include previously too short segments into clustering-scoring
				pScore->calcScores(); //the partition was altered with the result of the identification, so recalc the cluster-scores!
        pScore->printReport("clusteringAfterID.txt", true, (unsigned long)(time(NULL) - startTime));
      }
			printf("done! (%d distinct speakers found)\n", pScore->getClusterCount());
		} else {
			printf("skipped!\n");
		}
	}

	//do outputting of results, but only if ground-truth data is available
	printf("\nPhase 3:\tScoring...");
	if (pCorpus->getGT()->getGTtype() == sclib::gtMPEG7) {
		((SC_Corpus_MPEG7*)pCorpus)->highLevelFeaturesOut("hlFeatures.csv");
	} else if ( pCorpus->getGT()->getGTtype() != sclib::gtStandard && pCorpus->getGT()->getGTtype() != sclib::gtWesley) {
		if (this->pTweak->segmentationHandler.audioTypeMode != sclib::algorithm_nothing) {
			pScoreATC->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1, atcUncertainty);
			pScoreATC->printReport("atc_score_overall.txt", true, (unsigned long)(time(NULL) - startTime));
		}
		//no scd scoring here because storeIDs() within clusterSpeakers() removes all old boundaries
		if (this->pTweak->segmentationHandler.vUvDetectorMode != sclib::algorithm_nothing && pCorpus->getGT()->getGTtype() == sclib::gtTIMIT) {
			pScoreVUv->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1, vuvUncertainty);
			pScoreVUv->printReport("vuv_score_overall.txt", true, (unsigned long)(time(NULL) - startTime));
		}
		pFinalScore = new SC_SpeakerScore_Classification(this->pTweak, pCorpus->getGT());
		pFinalScore->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1, scdUncertainty); //normally, scd and speaker classification use the same (or at least comparable) features
		pFinalScore->printReport("finalResult.txt", true, (unsigned long)(time(NULL) - startTime), true);
	}
	//pCorpus->getGT()->frameListOut("frames.txt", 0, pCorpus->getGT()->getAudioSampleCount()-1);
	printf("done!");

	printf("\nPhase 4:\tCleaning up...");
  sclib::destructLinkedList(pAllShortClusters);
  MFree_0D(pScoreID);
  MFree_0D(pAllShortClusters);
	MFree_0D(pScoreATC);
	MFree_0D(pScoreSCD);
	MFree_0D(pScoreVUv);
	MFree_0D(pFinalScore);
  MFree_0D(pScore); //this releases all the partitions created by the clustering and, within, all the clusters, too
	sclib::destructLinkedList(pAllClusters);
	MFree_0D(pAllClusters);
  MFree_0D(pCorpus);
  MFree_0D(pModeller);
  MFree_0D(pExtractor);
  MFree_0D(pSegmenter);
  MFree_0D(pSigHandler);
	MFree_0D(pEnhancer);

	MFree_3D(explicitModels);

	printf("done!");
	printf("\n\nThe SC_Lib Testsuite says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	simple testscript for the ESPS pitch tracker
//====================================================================================================================
bool SC_MainTasks::pitchTest() {
	bool finalResult = true;
	long int segmentStart, segmentEnd;
  SC_Corpus *pCorpus;
  SC_FeatureHandler extractor(this->pTweak);
  SC_Signal *pSignal;
	SV_Data **pFeatures;
	char fileName[sclib::bufferSize];
	SC_Synthesis synthesizer(this->pTweak);

	printf("\nThis is the ESPS pitch tracker tester");
  printf("\n=======================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	//this are the (visible) standard parameters in wavesurfer 1.8.5 that is used for comparison
	/*this->pTweak->featurePitch.method = sclib::modeESPS;
	this->pTweak->featurePitch.esps_max_f0 = 400.0f;
	this->pTweak->featurePitch.esps_min_f0 = 60.0f;
	this->pTweak->featurePitch.esps_wind_dur = 0.0075f;
	this->pTweak->featurePitch.frameSize = 10;
	this->pTweak->featurePitch.frameStep = 10;
	this->pTweak->featurePitch.esps_freq_weight = 0.02f;*/
	printf("done!\n");

  printf("\n    Phase I:\tLoading Signal...");
	segmentStart = 0;
	segmentEnd = pCorpus->getGT()->getAudioSampleCount()-1;
	pSignal = pCorpus->loadSignal((unsigned long int&)segmentStart, (unsigned long int&)segmentEnd);
	printf("done!");

	if (pSignal != NULL && pSignal->GetLen() > 0) {
		//extract features
		printf("\n    Phase II:\tExtracting Features...\n");
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, sclib::featurePitch|sclib::featureMFCC);
		
		//output result in textual and audible form
		sclib::classOut("esps_pitch.txt", pFeatures[sclib::bitPosition(sclib::featurePitch)], this->pTweak);
		sprintf(fileName, "%shumm.wav", this->pTweak->debug.debugDir);
		bool harmonics[] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
		synthesizer.generateHumFromPitch(fileName, pFeatures[sclib::bitPosition(sclib::featurePitch)], harmonics, 64, true);

		extractor.freeFeatures(pFeatures);
		printf("done!");
	}
  
	printf("\nPhase 3:\tCleaning up...");
  MFree_0D(pCorpus);

	printf("done! => g'bye\n");

	return finalResult;
}

//====================================================================================================================
//	"training script" (prints the 2 thresholds to the console) for the LZL silence detector
//====================================================================================================================
bool SC_MainTasks::trainSilenceLZL() {
	bool finalResult = true;
	SC_Corpus *pCorpus;
	SV_Data **pFeatures = NULL;
	SC_FeatureHandler extractor(this->pTweak, true);
	SC_Segmentation_Silence_LZL silencer(this->pTweak);
	SC_Signal *pSignal = NULL;
	unsigned long int segmentStart, segmentEnd;

  printf("\nThis is the training script for the LZL silence detector");
  printf("\n=========================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  //a file with pure silence from start to end is expected, hence no groun-truth is necessary, thus use the MPEG7-corpus (it is gt-free)
	pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	segmentStart = 0;
	segmentEnd = pCorpus->getGT()->getAudioSampleCount() - 1;
	printf("done!\n");

	printf("\nPhase 1:\tInfering proper thresholds on ZCR and STE to detect silence");
	pSignal = pCorpus->loadSignal(segmentStart, segmentEnd);
	pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, silencer.getUsedFeatures(), silencer.getSpecialFeatureParameters());
	MFree_0D(pSignal);
	extractor.prepareFeatureSet(pFeatures, silencer.getSpecialFeatureParameters());
	silencer.trainClassifier(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures);
	extractor.freeFeatures(pFeatures);
	printf("done!\n");

	printf("\nPhase 4:\tCleaning up...");
  MFree_0D(pCorpus);
	printf("done!");
	printf("\n\nThe training script says: \"Good bye\".\n");

	return finalResult;
}

//====================================================================================================================
//	test script for a audio type classifier
//====================================================================================================================
bool SC_MainTasks::testAudioType() {
	bool finalResult = true;
  unsigned long int sceneNr = 1, atcUncertainty;
  long int sceneStart, sceneEnd;
	char *extension, previousResults[sclib::bufferSize];
  time_t startTime = time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Signal *pSignal = NULL;
  SV_Data **pFeatures = NULL, *pPitch = NULL;
	SC_Score_AudioTypeClassification *pScoreATC = NULL;
  SC_FeatureHandler *pExtractor = NULL;
  SC_SegmentationHandler *pSegmenter = NULL;
  SC_SignalHandler *pSigHandler = NULL;

	printf("\nThis is a testscript for audio type classifiers in the sclib");
  printf("\n=============================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  pExtractor = new SC_FeatureHandler(this->pTweak);
  pSegmenter = new SC_SegmentationHandler(this->pTweak);
  pSigHandler = new SC_SignalHandler(this->pTweak);
  //decide which corpus to use:
	extension = sclib::extractExtension(this->audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) != 0) { //no corpus-file, so assume SCiVo or MPEG7 type data
		if (strncmp(this->segmentationFile, "", sclib::bufferSize) != 0) {
			pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile); //SCiVo-files do have ground-truth
		} else {
			pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile); //...MPEG7-files don't (maybe cutlists...)
		}
  } else { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, false); //TODO: create scenes normally
  }
	pScoreATC = new SC_Score_AudioTypeClassification(this->pTweak, pCorpus->getGT());
	MFree_1D(extension);
  if (this->pTweak->modelHandler.onlyThisSpeaker != NULL && strlen(this->pTweak->modelHandler.onlyThisSpeaker) > 0) {
    this->pTweak->setDebugPrefix(this->pTweak->modelHandler.onlyThisSpeaker);
  }
	if (true != pCorpus->getGT()->checkConsistency(0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::modeGroundtruth)) {
		printf("%s", "Consistency errors exist in the ground truth -> see transcript!");
	}
	atcUncertainty = pCorpus->getGT()->getConverter()->ms2sample(pSegmenter->getUncertaintyRegionWidth(this->pTweak->segmentationHandler.audioTypeMode));
	printf("done!\n");

	printf("\nPhase 1:\tProcessing Scenes...\n");
	for (unsigned long int y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
	  
		pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
		if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
		
			printf("\n  Processing Scene %d...", sceneNr);
		
			//skip scene?
			if ((sceneNr < this->pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & this->pTweak->general.sceneSelection) && (this->pTweak->general.sceneSelection != 0))) {
				printf(" skipped!"); 
				sceneNr++; 
				y = sceneEnd; 
				continue;
			}

			//load previous work if wished and existant
			sprintf(previousResults, "%s_%d%s\0", this->pTweak->general.featurePrefix, sceneNr, ".dat");
			if (sclib::fileExists(previousResults) == true) {
				pFeatures = pExtractor->loadFeatures(previousResults);
			} else {

				//load signal
				printf("\n    Phase Ia:\tLoading Signal...");
 				pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
				printf("done!");

				//extract features
				printf("\n    Phase II:\tExtracting Features...");
				pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pSegmenter->getUsedFeatures(), pSegmenter->getSpecialFeatureParameters());
				MFree_0D(pSignal); 
				printf("done!");
	
				//save previous work for future runs if wished
				if (this->pTweak->general.featurePrefix!=NULL && strncmp(this->pTweak->general.featurePrefix, "", 1)!=0) {
					sprintf(previousResults, "%s_%d%s\0", this->pTweak->general.featurePrefix, sceneNr, ".dat");
					pExtractor->saveFeatures(previousResults, pFeatures);
				}
			} //nothing to load, extract features from scratch

			//silence detection
			printf("\n    Phase IIIa:\tRemoving Silence...");
			pSegmenter->silenceDetection(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
			printf("done!");

			//audio type classification
			printf("\n    Phase IIIb:\tAudio Type Classification...");
			pSegmenter->audioClassification(pCorpus->getGT(), sceneStart, sceneEnd, pFeatures, this->pTweak->segmentationHandler.audioTypeMode);
			pScoreATC->calcScores(sceneStart, sceneEnd, atcUncertainty);
			pScoreATC->printReport("atc_score_per_scene.txt", false, 0, false, true);
			printf("done!");

			//pSigHandler->storeSignal("noise.wav", sceneStart, sceneEnd, pCorpus->getGT()->getSignalPrototype(), pCorpus->getGT(), sclib::atNoise);

			sceneNr++;
			y = sceneEnd;

			if (true != pCorpus->getGT()->checkConsistency(sceneStart, sceneEnd, sclib::modeHypothesized)) {
				printf("%s", "Consistency errors exist in the algorithmic results -> see transcript!");
			}

			if (sceneNr > this->pTweak->general.lastScene) {printf("\nSkipping remaining scenes!\n"); break;}
		} else {
			break;	//scene not valid
		} 
	}

	printf("\ndone!\n");

	//do outputting of results, but only if ground-truth data is available
	printf("\nPhase 3:\tScoring...");
	if (pCorpus->getGT()->getGTtype() != sclib::gtStandard && pCorpus->getGT()->getGTtype() != sclib::gtWesley) {
		pScoreATC->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1, atcUncertainty);
		pScoreATC->printReport("atc_score_overall.txt", true, (unsigned long)(time(NULL) - startTime));
	}
	//pCorpus->getGT()->frameListOut("frames.txt", 0, pCorpus->getGT()->getAudioSampleCount()-1);
	printf("done!");

	printf("\nPhase 4:\tCleaning up...");
	MFree_0D(pScoreATC);
  MFree_0D(pCorpus);
  MFree_0D(pExtractor);
  MFree_0D(pSegmenter);
  MFree_0D(pSigHandler);

	printf("done!");
	printf("\n\nThe SC_Lib Testsuite says: \"Good bye\".\n");
	return finalResult;
}

//====================================================================================================================
//	tests various feature-re-synthesis approaches like mfcc2wav etc.
//====================================================================================================================
bool SC_MainTasks::synthesisTest() {
	bool finalResult = true;
	long int segmentStart, segmentEnd;
  SC_Corpus *pCorpus;
  SC_FeatureHandler extractor(this->pTweak);
	SC_SegmentationHandler segmenter(this->pTweak);
	SC_Synthesis synthesizer(this->pTweak);
	SC_Conversion converter;
  SC_Signal *pSignal;
	SV_Data **pFeatures, *pFeature, *pTemp;
	char fileName[sclib::bufferSize], *extension;
	SC_Model **pModelSpec, **pModelPitch;
	int i;

	printf("\nThis is the Synthesizer");
  printf("\n=========================\n");

	printf("\nPhase 0:\tLoading Parameters...");
	extension = sclib::extractExtension(this->audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) != 0) { //no corpus-file, so assume SCiVo or MPEG7 type data
		if (strncmp(this->segmentationFile, "-", sclib::bufferSize) != 0) {
			pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile); //SCiVo-files do have ground-truth
		} else {
			pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile); //...MPEG7-files don't (maybe cutlists...)
		}
  } else { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, true);
  }
	MFree_1D(extension);
	pCorpus->getGT()->setSegment(0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atSpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
	pCorpus->getGT()->setSegment(0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atSpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
	this->pTweak->segmentationHandler.audioTypeMode = 0; //switch off unnecessary algorithms here (we need only to detect silence and voicing in speech)
	this->pTweak->segmentationHandler.changeDetectorMode = 0;
	this->pTweak->featureSpectrum.frameSize = 32;
	this->pTweak->featureSpectrum.frameStep = 16;
	this->pTweak->featureSpectrum.FFTsize = 512;
	this->pTweak->featureSpectrum.preEmphasizeFactor = 0.0;
	this->pTweak->featureSpectrum.lowCut = 0.0;
	this->pTweak->featureSpectrum.highCut = pCorpus->getGT()->getAudioSampleRate() / 2.0;
	this->pTweak->featureSpectrum.createPhase = true;
	this->pTweak->featureSpectrum.logarithmize = true;
	this->pTweak->featureLpc.LPCorder = 14;
	this->pTweak->featureLpc.computeGain = false;
	this->pTweak->featureLpc.frameSize = this->pTweak->featureFormant.frameSize;
	this->pTweak->featureLpc.frameStep = this->pTweak->featureFormant.frameStep;
	synthesizer.setFftLength(this->pTweak->featureMfcc.fftSize);
	synthesizer.setOlaMaxIterations(1000);
	synthesizer.setOlaErrorTarget(4.0);
	synthesizer.setTaperingLength(this->pTweak->transform.taperingLength);
	synthesizer.setTaperingMode(this->pTweak->featureMfcc.window!=0 ? sclib::wndHamming : sclib::wndRectangle);
	printf("done!\n");

  printf("\n    Phase I:\tLoading Signal...");
	segmentStart = 0;
	segmentEnd = pCorpus->getGT()->getAudioSampleCount()-1;
	pSignal = pCorpus->loadSignal((unsigned long int&)segmentStart, (unsigned long int&)segmentEnd);
	converter.setAudioSampleRate(pSignal->SigPar.SRate);
	printf("done!");

	if (pSignal != NULL && pSignal->GetLen() > 0) {
		//extract features
		printf("\n    Phase II:\tExtracting Features...\n");
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, sclib::featureFormant|sclib::featureSpectrum|sclib::featurePitch|sclib::featureLPC|sclib::featureMFCC|segmenter.getUsedFeatures(), segmenter.getSpecialFeatureParameters(), false);
		//signal is still needed down there...

		//do necessary segmentation
		printf("\n    Phase III:\tDetecting Silence and Voicing of Speech...\n");
		segmenter.silenceDetection(pCorpus->getGT(), 0, pCorpus->getGT()->getAudioSampleCount()-1, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
		segmenter.unvoicedDetection(pCorpus->getGT(), 0, pCorpus->getGT()->getAudioSampleCount()-1, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);

		//splice & blend test
		printf("\n    Phase IVa:\tSplicing/Blending Test...\n");
		pFeature = pCorpus->getGT()->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureSpectrum)], 0, 0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atVoiced, sclib::atSilence); //get all voiced speech
		extractor.splice(pFeature, 12); //randomize time-order of frames
		for (i = 1; i <= 10; i += 2) {
			sprintf(fileName, "%sspliced_%d.wav", this->pTweak->debug.debugDir, i);
			pTemp = extractor.blend(pFeature, 5, i*0.25, 12);
			synthesizer.spectrum2wav(fileName, pTemp, false, this->pTweak->featureSpectrum.logarithmize, true);		
			MFree_0D(pTemp);
		}
		MFree_0D(pFeature);

		//general mfcc2wav test
		printf("\n    Phase IVb:\tMFCC2wav Test...\n");
		sprintf(fileName, "%smfcc2wav.wav", this->pTweak->debug.debugDir);
		synthesizer.mfcc2wav(fileName, pFeatures[sclib::bitPosition(sclib::featureMFCC)], this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 0, 1.0, pFeatures[sclib::bitPosition(sclib::featurePitch)]);
		
		//general lpc2wav test
		printf("\n    Phase IVc:\tLPC2wav Test...\n");
		sprintf(fileName, "%slpc2wav.wav", this->pTweak->debug.debugDir);
		synthesizer.lpc2wav(fileName, pFeatures[sclib::bitPosition(sclib::featureLPC)], this->pTweak->featureLpc.preEmphasizeFactor, this->pTweak->featureLpc.window, 0, 1.0, pFeatures[sclib::bitPosition(sclib::featurePitch)]);

    //tests for bark band averaging
		printf("\n    Phase IVd:\tBark Band Averaging Test...\n");
		sprintf(fileName, "%sbarkvoice.wav", this->pTweak->debug.debugDir);
		synthesizer.generateVoiceFromCB(fileName, pSignal->GetBuf_L(), pSignal->GetLen(), pSignal->SigPar.SRate, converter.ms2sample(this->pTweak->featureFbe.frameSize), converter.ms2sample(this->pTweak->featureFbe.frameStep));

		//pitch2hum test
		printf("\n    Phase IVe:\tPitch2Hum Test...\n");
		sprintf(fileName, "%shumm.wav", this->pTweak->debug.debugDir);
		bool harmonics[] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
		synthesizer.generateHumFromPitch(fileName, pFeatures[sclib::bitPosition(sclib::featurePitch)], harmonics, 64, true);

		//model inversion test
		printf("\n    Phase IVf:\tModel Inversion Test...\n");
		pFeature = pCorpus->getGT()->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureMFCC)], 0, 0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atVoiced, sclib::atSilence); //get all voiced speech
		SV_Data *pTempPitch, *pPitch = pCorpus->getGT()->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featurePitch)], 0, 0, pCorpus->getGT()->getAudioSampleCount()-1, sclib::atVoiced, sclib::atSilence); //get all voiced speech
		MArray_1D(pModelSpec, 4, SC_Model*, "SC_MainTasks.synthesisTest: pModelSpec");
		MArray_1D(pModelPitch, 4, SC_Model*, "SC_MainTasks.synthesisTest: pModelPitch");
		pModelSpec[0] = new SC_Model_FullGauss(this->pTweak);
		pModelSpec[1] = new SC_MixtureModel_GMM(this->pTweak, 1, pFeature->Col);
		pModelSpec[2] = new SC_Model_HMM(this->pTweak, 5, "1 1 1 0 0;0 1 1 1 0;0 0 1 1 1;0 0 0 1 1;0 0 0 0 1", 5, false, true, 100, true);
		pModelSpec[3] = new SC_MixtureModel_bGMM(this->pTweak, 10, pFeature->Col, true, true);
		pModelPitch[0] = new SC_Model_FullGauss(this->pTweak);
		pModelPitch[1] = new SC_MixtureModel_GMM(this->pTweak, 1, pPitch->Col);
		pModelPitch[2] = new SC_Model_HMM(this->pTweak, 5, "1 1 1 0 0;0 1 1 1 0;0 0 1 1 1;0 0 0 1 1;0 0 0 0 1", 5, false, true, 100, true);
		pModelPitch[3] = new SC_MixtureModel_bGMM(this->pTweak, 10, pPitch->Col, true, true);
		for (i = 0; i < 4; i++) {
			if (i == 2) { //prepeare features: special treatment for the hmm features to provide several examples of temporal order
				pTemp = extractor.splitFeatureSet(pFeature, 10, 10);
				pTemp->Hdr = pFeature->Hdr;
				pTempPitch = extractor.splitFeatureSet(pPitch, 10, 10);
				pTempPitch->Hdr = pPitch->Hdr;
			}	else {
				pTemp = pFeature;
				pTempPitch = pPitch;
			}
			if (i != 3) { //train model
				pModelSpec[i]->TrainModel(pTemp, 0);
				pModelPitch[i]->TrainModel(pTempPitch, 0);
			} else {
				((SC_MixtureModel_bGMM*)pModelSpec[i])->TrainModel(pTemp, 0, -1, true, 1.0, true, -1.0);
				((SC_MixtureModel_bGMM*)pModelPitch[i])->TrainModel(pTempPitch, 0, -1, true, 1.0, true, -1.0);
			}
			if (pTemp != pFeature) { //sample from the model
				sclib::destructLinkedList(pTemp);
				sclib::destructLinkedList(pTempPitch);
			}
			pTempPitch = pModelPitch[i]->drawSamplesFromDistribution((int)(ceil(200.0 * (double)(pPitch->Row)/(double)(pFeature->Row))));
			pTempPitch->Hdr = pPitch->Hdr;
			extractor.cleanPitch(pTempPitch, 0, this->pTweak->featurePitch.esps_min_f0);
			//sclib::quickSort(pTempPitch->Mat, 0, pTempPitch->Row-1, pTempPitch->Col, 0); //TODO: remove pitch sorting from low to high?!
			pTemp = pModelSpec[i]->drawSamplesFromDistribution(200);
			pTemp->Hdr = pFeature->Hdr;
			sprintf(fileName, "%smodelinversion_%d.wav", this->pTweak->debug.debugDir, i); //make the samples audible
			//pTemp = sclib::directOutput(extractor.blend(pTemp, 10, 1.0), pTemp, false);
			//synthesizer.spectrum2wav(fileName, pTemp, false, this->pTweak->featureSpectrum.logarithmize, true);
			//synthesizer.mfcc2wav(fileName, pTemp, this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 10, 1.0, NULL);
			//sprintf(fileName, "%smodelinversion_pitch_%d.wav", this->pTweak->debug.debugDir, i); //make the samples audible
			synthesizer.mfcc2wav(fileName, pTemp, this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 10, 1.0, pTempPitch);
			MFree_1D(pTemp);
			MFree_1D(pTempPitch);
			MFree_0D(pModelSpec[i]);
			MFree_0D(pModelPitch[i]);
		}
		MFree_0D(pFeature);
		MFree_0D(pPitch);
		MFree_1D(pModelSpec);
		MFree_1D(pModelPitch);

		//formant removal test
		printf("\n    Phase IVg:\tFormant Removal Test...\n");
		sprintf(fileName, "%sformantremoval.wav", this->pTweak->debug.debugDir);
		synthesizer.removeFormantInfluence(pSignal->GetBuf_L(), pSignal->GetLen(), pSignal->SigPar.SRate, pFeatures[sclib::bitPosition(sclib::featureFormant)], 2, 3.0);
		extractor.freeFeatures(pFeatures);
		pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, sclib::featureSpectrum, NULL, false);
		synthesizer.spectrum2wav(fileName, pFeatures[sclib::bitPosition(sclib::featureSpectrum)], false, this->pTweak->featureSpectrum.logarithmize, true);

		extractor.freeFeatures(pFeatures);
		printf("done!");
	}
  
	printf("\nPhase 3:\tCleaning up...");
	MFree_0D(pSignal);
  MFree_0D(pCorpus);
	printf("done! => g'bye\n");

	return finalResult;
}

//====================================================================================================================
//	tests the method called by java/videana
//====================================================================================================================
bool SC_MainTasks::javaTest(int *shotList, int shotListLength) {
	bool finalResult = true;
	unsigned long int res;
	long segRes = 0, ids = 0, probs = 0;
	int cols = 0;
	double **prob;
	long int *seg, *id;
	
	printf("\nThis is javaTest");
  printf("\n=================\n");

	printf("\nPhase 1:\tCalling audioSegmentation...");
	res = audioSegmentation(this->audioFile, this->pTweak, this->videoFrameRate, &segRes, &ids, &probs, &cols, shotList, shotListLength);
	seg = (long int*)(segRes);
	id = (long int*)(ids);
	prob = (double**)(probs);
	finalResult = (res > 0) ? true : false;
	printf("done! => %d rows in result, %d cols in prob. matrix\n", res, cols);

	printf("\nPhase 2:\tCleaning up...");
	MFree_1D(seg);
	MFree_1D(id);
	MFree_2D(prob);
	printf("done! => g'bye\n");

	return finalResult;
}

//====================================================================================================================
//	testscript for the branch&bound feature selection on conTIMIT data
//====================================================================================================================
bool SC_MainTasks::kottiBICtest(char *bicTrainFile, char *durationModelFile, char *fsClass1File, char *fsClass2File) {
	bool finalResult = true;
	SC_Corpus *pCorpus = NULL;
	unsigned long int y, i;
	long int start, end;
	SC_Signal *pSignal;
	SC_FeatureHandler *pExtractor = new SC_FeatureHandler(this->pTweak, true);
	SV_Data **pFeatures, *pFeatureStart = NULL, *pFeatureHook;
	bool *labels = NULL;
	unsigned char *selected = NULL;
	unsigned int resultingColumnCount;
	char *fsFiles[2] = {fsClass1File, fsClass2File};
	unsigned int classCount[2] = {0, 0};
	std::list<unsigned long int> durationList;
	unsigned long int *durations;
	SC_MatrixFunctions matFunc;
	double meanDuration, lambdaDuration, sampleDispersionDuration, lambda;
	SC_SegmentationHandler segmenter(this->pTweak, false);
	bool alreadyTweaked = false;
  SC_Score_ChangeDetection *pScoreSCD = NULL;
	int changes;

	printf("\nThis is kottiBICtest");
  printf("\n=====================\n");



	//step 0: prepare some parameters
	this->pTweak->segmentationHandler.changeDetectorMode = sclib::algorithm_cd_KBK;
	resultingColumnCount = 24;
	this->pTweak->segmentationChangesKbk.lambda = 1.3; //1.7625; //0.586291; //0.982; //0.586291; //0.75;
	this->pTweak->segmentationChangesKbk.r = 1.5584415; //3.116883 / 2.0; //3.122581;
	this->pTweak->segmentationChangesKbk.tolerance = 1.0;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[0] = 0xFF; //0xE7;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[1] = 0xFF; //0x7F; //0xFE; //0x7F;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[2] = 0xFF; //0x78; //0x1E; //0x78;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[3] = 0xFF; //0xE8; //0x17; //0xE8;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[4] = 0xFF; //0x0E; //0x70; //0x0E;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[5] = 0x00;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[6] = 0x00;
	this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[7] = 0x00;
	this->pTweak->segmentationChangesKbk.mfccParameters.MFCCorder = 13; //36;
	this->pTweak->segmentationChangesKbk.mfccParameters.addDeltaDeltas = false; //false
	this->pTweak->segmentationChangesKbk.mfccParameters.addDeltas = false; //false
	this->pTweak->segmentationChangesKbk.mfccParameters.frameSize = 20; //16;
	this->pTweak->segmentationChangesKbk.mfccParameters.frameStep = 10;

	this->pTweak->segmentationChangesKbk.mfccParameters.method = sclib::modeSClib; //TODO
	this->pTweak->segmentationChangesKbk.mfccParameters.dEnergy = false;
	this->pTweak->segmentationChangesKbk.mfccParameters.CMN = false;
	this->pTweak->segmentationChangesKbk.mfccParameters.fftSize = 512;
	this->pTweak->segmentationChangesKbk.mfccParameters.filterBankSize = 24;
	this->pTweak->segmentationChangesKbk.mfccParameters.preEmphasizeFactor = 0.97;
	this->pTweak->segmentationChangesKbk.mfccParameters.sclib_frequencyScale = sclib::scaleMel;
	this->pTweak->segmentationChangesKbk.mfccParameters.sclib_smoothing = sclib::smoothNone;
	this->pTweak->segmentationChangesKbk.mfccParameters.window = sclib::wndHamming;
	this->pTweak->segmentationChangesKbk.mfccParameters.sclib_maxFilterBankFrequency = 0.0;
	this->pTweak->segmentationChangesKbk.mfccParameters.sclib_minFilterBankFrequency = 0.0;

	alreadyTweaked = true; //set true if the above parameters are already tweaked so that training can be discarded

	/*this->pTweak->segmentationChangesKbk.mfccParameters.addDeltaDeltas = false;
	this->pTweak->segmentationChangesKbk.mfccParameters.addDeltas = false;
	this->pTweak->segmentationChangesKbk.mfccParameters.MFCCorder = 12;
	this->pTweak->segmentationChangesKbk.mfccParameters.filterBankSize = 24;*/


	//step 1: find optimal feature set
	if (alreadyTweaked == false) {
		printf("\nPhase 1:\tFinding best feature set:\n");
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[0] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[1] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[2] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[3] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[4] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[5] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[6] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[7] = 0xFF;
		this->pTweak->segmentationChangesKbk.mfccParameters.addDeltas = false;
		this->pTweak->segmentationChangesKbk.mfccParameters.addDeltaDeltas = false;
		for (i = 0; i < 2; i++) {
			printf("\n\tReading data of class %d", i);
			pCorpus = new SC_Corpus_MAC(this->pTweak, fsFiles[i], false);
			for (y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
				pCorpus->getGT()->getNextBoundary(y, start, end, sclib::atSceneBoundary, sclib::searchForward);
				if (start != sclib::noSegment && end != sclib::noSegment) {
					pSignal = pCorpus->loadSignal((unsigned long &)start, (unsigned long &)end);
					pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, start, end, segmenter.getUsedFeatures(sclib::algorithm_cd_KBK), segmenter.getSpecialFeatureParameters(sclib::algorithm_cd_KBK), true);
					MFree_0D(pSignal); 

					//unite and link features to a list
					if (pFeatureStart == NULL) {
						pFeatureStart = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount());
						pFeatureHook = pFeatureStart;
					} else {
						pFeatureHook->Next = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount());
						pFeatureHook = pFeatureHook->Next;
					}
					pExtractor->freeFeatures(pFeatures);
					classCount[i]++;

					printf(".");

					y = end;
				} else {
					break;
				}
			}
			MFree_0D(pCorpus);
		} //for the two classes
		if (classCount[0] > 0 && classCount[1] > 0) { //find best feature set
			//create labels for the data elements that correspond with the linked list of features
			MArray_1D(labels, classCount[0]+classCount[1], bool, "SC_MainTasks.kottiBICtest: labels");
			for (y = 0; y < classCount[0]+classCount[1]; y++) {
				labels[y] = (y < classCount[0]) ? true : false;
			}
			printf("\n\tDo Feature Selection: ", resultingColumnCount, pFeatureStart->Col);
			selected = pExtractor->branchBoundFeatureSelection(pFeatureStart, labels, resultingColumnCount);
			printf("\n\tFinally selected columns (%d out of %d): ", resultingColumnCount, pFeatureStart->Col);
			for (y = 0; y < (unsigned long int)(pFeatureStart->Col); y++) {
				if (sclib::bitTest(selected[y/8], sclib::bit(y%8)) == true) { //the y'th col in the feature set is the y%8'th bit in selected[y/8]
					printf("%d, ", y);
				}
			}
			printf("\r\r!\n");
			for (y = 0; y < 8; y++) { //store the result
				if (y < (unsigned long int)(ceil(sclib::ld((double)(pFeatureStart->Col))))) {
					this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[y] = selected[y];
				} else {
					this->pTweak->segmentationChangesKbk.mfccParameters.coeffSelection[y] = 0x00;
				}
			}
			MFree_0D(labels);
			MFree_0D(selected);
			sclib::destructLinkedList(pFeatureStart);
		}



		//step 2: find distribution of speaker change points in conTIMIT test data
		printf("\nPhase 2:\tFind utterence length distribution parameters:\n");
		pCorpus = new SC_Corpus_TIMIT(this->pTweak, durationModelFile, false);
		for (y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
			pCorpus->getGT()->getNextBoundary(y, start, end, sclib::atSpeakerBoundary, sclib::searchForward);
			if (start != sclib::noSegment && end != sclib::noSegment) {
				durationList.push_back(sclib::round(pCorpus->getGT()->getConverter()->sample2ms(end-start+1)/1000.0)); //model in seconds, not samples
				y = end;
			} else {
				break;
			}
		}
		y = (unsigned long int)(durationList.size());
		MArray_1D(durations, y, unsigned long int, "SC_MainTasks.kottiBICtest: durations");
		for (i = 0; i < y; i++) {
			durations[i] = durationList.front();
			durationList.pop_front();
		}
		meanDuration = matFunc.mean(durations, y);
		lambdaDuration = matFunc.lambda(durations, y, &meanDuration);
		sampleDispersionDuration = sqrt(matFunc.variance(durations, y, &meanDuration));
		MFree_1D(durations);
		//SC_Segmentation_Changes_LZW segmenter(this->pTweak, sclib::modeSpeakerChange);
		//SC_Model_Pareto *pModel = segmenter.getSegmentDurationDistribution(pCorpus->getGT(), 0, pCorpus->getGT()->getAudioSampleCount()-1);
		//sclib::classOut("duration.txt", pModel, this->pTweak);
		MFree_0D(pCorpus);
		printf("Found mean duration under inverse gaussian model with lambda=%f: %fs, standard deviation: %fs\n", lambdaDuration, meanDuration, matFunc.igStd(meanDuration, lambdaDuration));
		this->pTweak->segmentationChangesKbk.r = meanDuration / 2.0; //store the result


		//step 3: find a suitable lambda
		printf("\nPhase 4:\tFind a suitable lambda:\n");
		pCorpus = new SC_Corpus_TIMIT(this->pTweak, bicTrainFile, false);
		start = 0;
		end = pCorpus->getGT()->getAudioSampleCount()-1;
		pSignal = pCorpus->loadSignal((unsigned long &)start, (unsigned long &)end);
		pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, start, end, segmenter.getUsedFeatures(sclib::algorithm_cd_KBK), segmenter.getSpecialFeatureParameters(sclib::algorithm_cd_KBK), true);
		MFree_0D(pSignal);
		SC_Segmentation_Changes_KBK kbk(this->pTweak, this->pTweak->segmentationChangesKbk.tolerance);
		pScoreSCD = new SC_Score_ChangeDetection(this->pTweak, pCorpus->getGT());
		this->pTweak->segmentationChangesKbk.lambda = 0.75684551082902707; //TODO: something big, but not too big to yield infinite bic-values
		double bestLambda = this->pTweak->segmentationChangesKbk.lambda, bestF1 = 0.0, f1;
		while (true) {
			//TODO
			lambda = kbk.findLambdaKotti(pCorpus->getGT(), start, end, pFeatures, this->pTweak->segmentationChangesKbk.r, this->pTweak->segmentationChangesKbk.tolerance);
			pScoreSCD->calcScores(start, end, pCorpus->getGT()->getConverter()->ms2sample(segmenter.getUncertaintyRegionWidth(sclib::algorithm_cd_KBK)));
			char buf[sclib::bufferSize];
			f1 = pScoreSCD->getF1(pScoreSCD->class2idx(sclib::atSpeakerBoundary));
			sprintf(buf, "lambda=%f REC=%f PRC=%f F1=%f MDR=%f FAR=%f", this->pTweak->segmentationChangesKbk.lambda, pScoreSCD->getRecall(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getPrecision(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), f1, pScoreSCD->getMDR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getFAR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)));
			printf("\n%s", buf);
			sclib::stringOut("l.txt", buf, this->pTweak);
			if (f1 > bestF1) {
				bestF1 = f1;
				bestLambda = this->pTweak->segmentationChangesKbk.lambda;
			}
			//pScoreSCD->printReport("report.txt", false);
			if (lambda < this->pTweak->segmentationChangesKbk.lambda) {
				this->pTweak->segmentationChangesKbk.lambda = lambda;
			} else {
				break;
			}
			//if (pScoreSCD->getMDR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)) > 0.3) {
			//	break;
			//}
		}
		this->pTweak->segmentationChangesKbk.lambda = bestLambda;
		//this->pTweak->segmentationChangesKbk.lambda = kbk.findLambda(pCorpus->getGT(), start, end, pFeatures, this->pTweak->segmentationChangesKbk.r, this->pTweak->segmentationChangesKbk.resolutionFactor, this->pTweak->segmentationChangesKbk.tolerance);
		pExtractor->freeFeatures(pFeatures);
		MFree_0D(pScoreSCD);
		MFree_0D(pCorpus);
	} //alreadyTweaked==false

	//step 4: test the algorithm
	printf("\nPhase 4:\tTest the algorithm:\n");	
	pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, false);
	pScoreSCD = new SC_Score_ChangeDetection(this->pTweak, pCorpus->getGT());
	for (y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
		pCorpus->getGT()->getNextBoundary(y, start, end, sclib::atSceneBoundary, sclib::searchForward);
		if (start != sclib::noSegment && end != sclib::noSegment) {
			pSignal = pCorpus->loadSignal((unsigned long &)start, (unsigned long &)end, true);
			pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, start, end, segmenter.getUsedFeatures(sclib::algorithm_cd_KBK), segmenter.getSpecialFeatureParameters(sclib::algorithm_cd_KBK), true);
			MFree_0D(pSignal);

			printf("\n\t\tSCD: ");
			changes = segmenter.speakerChangeDetection(pCorpus->getGT(), start, end, pFeatures, sclib::algorithm_cd_KBK);
			printf("%d changes found.", changes);
			pExtractor->freeFeatures(pFeatures);
			pScoreSCD->calcScores(start, end, pCorpus->getGT()->getConverter()->ms2sample(segmenter.getUncertaintyRegionWidth(sclib::algorithm_cd_KBK)));
			char *buffer = ((SC_Corpus_TIMIT*)pCorpus)->getCurrentAudioFileName(y);
			buffer = sclib::directOutput(sclib::extractFileName(buffer), buffer, true);
			printf("\n%s\tREC=%f\tPRC=%f\tF1=%f\tMDR=%f\tFAR=%f\n", buffer, pScoreSCD->getRecall(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getPrecision(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getF1(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getMDR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getFAR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)));
			char buf[sclib::bufferSize];
			sprintf(buf, "%f\t%f\t%f\t%f\t%f", pScoreSCD->getRecall(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getPrecision(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getF1(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getMDR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getFAR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)));
			sclib::stringOut("res.txt", buf, this->pTweak);
			MFree_1D(buffer);
			pScoreSCD->printReport("scdres.txt", y==0);

			y = end;
		} else {
			break;
		}
	}
	pScoreSCD->calcScores(0, pCorpus->getGT()->getAudioSampleCount()-1, pCorpus->getGT()->getConverter()->ms2sample(segmenter.getUncertaintyRegionWidth(sclib::algorithm_cd_KBK)));
	printf("\nOverall REC=%f\tPRC=%f\tF1=%f\tMDR=%f\tFAR=%f\n", pScoreSCD->getRecall(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getPrecision(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getF1(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getMDR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)), pScoreSCD->getFAR(pScoreSCD->class2idx(sclib::atSpeakerBoundary)));
	pScoreSCD->printReport("scdres_overall.txt");
	MFree_0D(pCorpus);
	MFree_0D(pScoreSCD);

	printf("\nPhase 5:\tCleaning up...");
	MFree_0D(pExtractor);
	printf("done! => g'bye\n");

	return finalResult;
}

//====================================================================================================================
//	model building part of the experiments from Reynolds, "Speaker Identification and Verification using Gaussian 
//  Mixture Speaker Models", 1995; Reynolds, Rose, "Robust Text-Independent Speaker Identification using Gaussian
//  Mixture Speaker Models", 1995
//====================================================================================================================
bool SC_MainTasks::buildReynoldsSpeakerModels(double reduceByPercent) {
	SC_Timer timer(3, this->pTweak);
	bool failure = false;
	char modelName[sclib::bufferSize];
	long int speakerId, segmentStart, segmentEnd , spkCount = 0;
	int res;
  SC_Corpus *pCorpus;
  SC_FeatureHandler extractor(this->pTweak);
	SC_SegmentationHandler segmenter(this->pTweak, false);
  SC_ModelHandler modeller(this->pTweak);
  SC_Signal *pSignal;
	SV_Data **pFeatures, *pFeature, *pFeatureFirst = NULL, *pFeatureHook;
	SC_Model *pModel; 

	printf("\nThis builds speaker models for the experiments in Reynolds, 1995");
  printf("\n=================================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  //decide which corpus to use:
  pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile);
	this->pTweak->unsetDebugPrefix(); //so that the filename of the model doesn't get messed up!
	this->pTweak->debug.useDebugPrefix = false;
	/*
	this->pTweak->modelHandler.maxSpeakerModelOrder = 32; //or: 16
	this->pTweak->mixtureModelGmm.EMthreshold = 0.01; //to always do maxEmIterations
	this->pTweak->mixtureModelGmm.kMeansIterations = 10;
	this->pTweak->mixtureModelGmm.maxEMiterations = 10; //or: 100
	this->pTweak->mixtureModelGmm.varianceLimit = 0.05; //between 0.01 .. 0.1
	this->pTweak->mixtureModelGmm.weightLimit = 0.0; //no weight limit
	this->pTweak->featureMfcc.addDeltaDeltas = false;
	this->pTweak->featureMfcc.addDeltas = false;
	this->pTweak->featureMfcc.CMN = false;
	this->pTweak->featureMfcc.coeffSelection[0] = 0xFE; //discard c0, retain all others
	this->pTweak->featureMfcc.coeffSelection[1] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[2] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[3] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[4] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[5] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[6] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[7] = 0xFF;
	this->pTweak->featureMfcc.dEnergy = false; //discarded anyway...
	this->pTweak->featureMfcc.fftSize = 512; //256 for 8kHz
	this->pTweak->featureMfcc.filterBankSize = 20;
	this->pTweak->featureMfcc.frameSize = 20;
	this->pTweak->featureMfcc.frameStep = 10;
	this->pTweak->featureMfcc.highCut = 0.0;
	this->pTweak->featureMfcc.lowCut = 0.0;
	this->pTweak->featureMfcc.method = sclib::modeSClib;
	this->pTweak->featureMfcc.MFCCorder = 20;
	this->pTweak->featureMfcc.preEmphasizeFactor = 0.97;
	this->pTweak->featureMfcc.window = sclib::wndHamming;
	this->pTweak->featureMfcc.sclib_smoothing = sclib::smoothNone;
	this->pTweak->featureMfcc.sclib_minFilterBankFrequency = 0.0;
	this->pTweak->featureMfcc.sclib_maxFilterBankFrequency = 4600.0; //end of davis & mermelsteins filterbank
	this->pTweak->featureMfcc.sclib_frequencyScale = sclib::scaleMel;
	this->pTweak->featureMfcc.sampleRate = 0.0;
	*/
	printf("done!\n");

	//SC_Signal *pSignalHook, *pSignalStart = NULL; //TODO wave output
	timer.startTimer(1);
  while (pCorpus->getGT()->getSpeakersCallback(speakerId) == true) {
		printf("\nPhase 1:\tProcessing Audio-Data speaker-wise for speaker %d\n", speakerId);
		spkCount++;

		while (pCorpus->getGT()->getSpeakersSegmentsCallback(speakerId, segmentStart, segmentEnd) == true) {
      printf("\n    Phase I:\tLoading Signal of speaker nr. %d...", spkCount);
			pSignal = pCorpus->loadSignal((unsigned long int&)segmentStart, (unsigned long int&)segmentEnd);
			printf("done!");

			if (pSignal != NULL && pSignal->GetLen() > 0) {
				//if (pSignalStart == NULL) { //TODO wave output
				//	pSignalStart = pSignal;
				//	pSignalHook = pSignalStart;
				//} else {
				//	pSignalHook->Next = pSignal;
				//	pSignalHook = pSignalHook->Next;
				//}
				
				//extract features
				printf("\n    Phase II:\tExtracting Features...");
				pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, sclib::featureMFCC);
				MFree_0D(pSignal); //TODO wave output
				printf("done!");

				/*
				//silence detection
				printf("\n    Phase IIIa:\tRemoving Silence...");
				segmenter.silenceDetection(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
				printf("done!");

				//label unvoiced speech
				printf("\n    Phase IIIb:\tLabeling Unvoiced Speech...");
				SV_Data *pPitch = segmenter.unvoicedDetection(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);
				MFree_0D(pPitch);
				printf("done!");
				*/

				//copy all speech frames together
				pFeature = pCorpus->getGT()->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureMFCC)], segmentStart, segmentStart, segmentEnd, sclib::atSpeech); //reynolds used v+uv!
				extractor.freeFeatures(pFeatures);

				//build linked list of all desired features of this video
				//pFeatureFirst is the pointer to the first element of the linked list, pFeatureHook to the current last element
				if (pFeature != NULL) {
					if (pFeatureFirst == NULL) {
						pFeatureFirst = pFeature;
						pFeatureHook = pFeatureFirst; 
					} else {
						pFeatureHook->Next = pFeature;
						pFeatureHook = pFeatureHook->Next;
					}
				}
			} //pSignal valid

		} //for all segments of one speaker
		printf("\ndone!\n");
		//SC_SignalHandler sigHandler(this->pTweak); //TODO wave output
		//sprintf(modelName, "%s_train.wav", pCorpus->getGT()->getSpeakerName(speakerId));
		//sigHandler.saveSignal(modelName, pSignalStart, true);
		//sclib::destructLinkedList(pSignalStart);

		printf("\nPhase 2:\tBuilding model...");
		if (pFeatureFirst->Next != NULL) {
			//ok = extractor.mergeOnDisk(pFeatureFirst); //save and load data to/from disk, to save one copy in memory...
			pFeature = pFeatureFirst->MergeData(); //merge in memory because we assume that one speaker doesn't have this much speech...
			sclib::destructLinkedList(pFeatureFirst);
		}
		if (pFeature != NULL) {
			if (reduceByPercent > 0.0) { //TODO: ensure reduceByPercent==0.0 for original Reyndolds experiments!
				extractor.reduceByPercent(pFeature, reduceByPercent);
				sclib::scalarOut("audioLength_train.txt", (double)(pCorpus->getGT()->getConverter()->audioFrame2ms(pFeature->Row, pFeature->Hdr.frameSize, pFeature->Hdr.frameStep))/1000.0, this->pTweak); //in [s]
			}
			//pModel = new SC_MixtureModel_GMM(this->pTweak, this->pTweak->modelHandler.maxSpeakerModelOrder, pFeature->Col);
			SC_ModelHandler handler(this->pTweak, false);
			int modelOrder = this->pTweak->modelHandler.maxSpeakerModelOrder;
			timer.startTimer(2);
			if (this->pTweak->modelHandler.orderGuessMode == sclib::guessParameterSearch &&
	  			this->pTweak->modelHandler.foregroundModelType == sclib::mtGMM_new) {
				timer.startTimer(3);
				modelOrder = handler.guessModelOrder(pFeature, NULL, sclib::mtGMM_new, 1, modelOrder, 1);
				timer.stopTimer(3);
			}
			pModel = handler.createRawModel(this->pTweak->modelHandler.foregroundModelType, NULL, modelOrder, pFeature->Col); //TODO: ensure the modeltype choosen is mtGMM_new with 32 mixtures for original reynolds experiment!

			res = pModel->TrainModel(pFeatureFirst);
			timer.stopTimer(2);
			if (handler.isMixtureModel(pModel) == true) {
				sclib::scalarOut("mixtureCount.txt", ((SC_MixtureModel*)pModel)->getMixtureCount(), this->pTweak);
			}
			SC_Signature *pSig = pModel->toSignature();
			if (pSig != NULL) {
				sclib::classOut("signatures.txt", pSig, this->pTweak, ios_base::out|ios_base::app, pCorpus->getGT()->getSpeakerName(speakerId));
			}
			//sprintf(modelName, "%s_train.txt", pCorpus->getGT()->getSpeakerName(speakerId)); //TODO model matlab output
			//sclib::classOut(modelName, pModel, this->pTweak, ios_base::out|ios_base::app, pCorpus->getGT()->getSpeakerName(speakerId));
			MFree_0D(pSig);
			if (res != SVLIB_Fail) {
				sprintf(modelName, "%s.gmm", pCorpus->getGT()->getSpeakerName(speakerId));
				res = modeller.saveModel(modelName, pModel, true);
				if (res <= 0) {
					printf("\nSorry, the model could not be saved\n");
					failure = true;
				} else {
					sprintf(modelName, "%s%s.gmm", this->pTweak->debug.debugDir, pCorpus->getGT()->getSpeakerName(speakerId));
					sclib::stringOut("reynoldsModelList.txt", modelName, this->pTweak, "\n");
				}
			} else {
				printf("\nSorry, the model could not be built\n");
				failure = true;
			}
			MFree_0D(pModel);
			MFree_0D(pFeatureFirst);
		} else {
			printf("\nSorry, data couldn't be loaded\n");
			failure = true;
		}

	  printf("done!\n");		
	} //for all speakers
	timer.stopTimer(1);
  
	printf("\nPhase 3:\tCleaning up...");

  MFree_0D(pCorpus);
	timer.logIt("timer_train.txt");

	printf("done!");
	printf("\nSimon says: \"Good bye\".\n");

	return !failure;
}

//====================================================================================================================
//	model test part of the experiments from Reynolds, "Speaker Identification and Verification using Gaussian 
//  Mixture Speaker Models", 1995; Reynolds, Rose, "Robust Text-Independent Speaker Identification using Gaussian
//  Mixture Speaker Models", 1995
//====================================================================================================================
bool SC_MainTasks::testReynoldsSpeakerModels(char *modelListFile, double reduceByPercent) {
	SC_Timer timer(2, this->pTweak);
	bool failure = false;
	char *tmp, line[sclib::bufferSize], *lastDot;
	long int speakerId, segmentStart, segmentEnd;
	int modelCnt = 0, segmentCnt, spkCnt = 0;
  SC_Corpus *pCorpus;
  SC_ModelHandler modeller(this->pTweak);
  SC_FeatureHandler extractor(this->pTweak);
	SC_SegmentationHandler segmenter(this->pTweak, false);
  SC_Signal *pSignal;
	SV_Data **pFeatures, *pFeature, *pTemp;
	SC_Cluster *pFirst = NULL, *pHook;
	SC_Model *pModel;
	FILE *modelList;
	double maxScore;
	long int maxId;
	unsigned int tested = 0, correct = 0, fail = 0;

	printf("\nThis tests speaker models like in the experiments in Reynolds, 1995");
  printf("\n====================================================================\n");

	printf("\nPhase 0a:\tLoading Parameters...");
  pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, true, true);
	this->pTweak->unsetDebugPrefix(); //so that the filename of the model doesn't get messed up!
	this->pTweak->debug.useDebugPrefix = false;
	/*
	this->pTweak->modelHandler.maxSpeakerModelOrder = 32; //or: 16
	this->pTweak->mixtureModelGmm.EMthreshold = 0.01; //to always do maxEmiterations
	this->pTweak->mixtureModelGmm.kMeansIterations = 10;
	this->pTweak->mixtureModelGmm.maxEMiterations = 10; //or: 100
	this->pTweak->mixtureModelGmm.varianceLimit = 0.05; //between 0.01 .. 0.1
	this->pTweak->mixtureModelGmm.weightLimit = 0.0; //no weight limit
	this->pTweak->featureMfcc.addDeltaDeltas = false;
	this->pTweak->featureMfcc.addDeltas = false;
	this->pTweak->featureMfcc.CMN = false;
	this->pTweak->featureMfcc.coeffSelection[0] = 0xFE; //discard c0, retain all others
	this->pTweak->featureMfcc.coeffSelection[1] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[2] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[3] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[4] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[5] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[6] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[7] = 0xFF;
	this->pTweak->featureMfcc.dEnergy = false; //discarded anyway...
	this->pTweak->featureMfcc.fftSize = 512; //256 for 8kHz
	this->pTweak->featureMfcc.filterBankSize = 20;
	this->pTweak->featureMfcc.frameSize = 20;
	this->pTweak->featureMfcc.frameStep = 10;
	this->pTweak->featureMfcc.highCut = 0.0;
	this->pTweak->featureMfcc.lowCut = 0.0;
	this->pTweak->featureMfcc.method = sclib::modeSClib;
	this->pTweak->featureMfcc.MFCCorder = 20;
	this->pTweak->featureMfcc.preEmphasizeFactor = 0.97;
	this->pTweak->featureMfcc.window = sclib::wndHamming;
	this->pTweak->featureMfcc.sclib_smoothing = sclib::smoothNone;
	this->pTweak->featureMfcc.sclib_minFilterBankFrequency = 0.0;
	this->pTweak->featureMfcc.sclib_maxFilterBankFrequency = 4600.0; //end of davis & mermelsteins filterbank
	this->pTweak->featureMfcc.sclib_frequencyScale = sclib::scaleMel;
	this->pTweak->featureMfcc.sampleRate = 0.0;
	*/
	printf("done!\n");

	printf("\nPhase 0b:\tLoading Models...");
	modelList = fopen(modelListFile, "r");
	if (modelList != NULL) {
		while (!feof(modelList)) {
			sclib::readline(modelList, line, sclib::bufferSize);
			if (sclib::fileExists(line) == true) {
				pModel = modeller.loadModel(line, this->pTweak->modelHandler.foregroundModelType); //TODO: ensure mtGMM_new for original reynolds experiment //sclib::mtGMM_new);
				if (pModel != NULL) { //build a linked list of all the models
					tmp = sclib::extractFileName(line);
					lastDot = strrchr(tmp, '.');
					if (lastDot != NULL) {
						*lastDot = '\0'; //remove file extension to really get the speaker name if the file was named e.g. "some_path/speakername.gmm"
					}
					if (pFirst == NULL) {
						pFirst = new SC_Cluster(this->pTweak, NULL, NULL, NULL, NULL, pModel, 1, pCorpus->getGT()->getSpeakerIDfromName(tmp, (int)(strlen(tmp))), tmp, true);
						pHook = pFirst;
					} else {
						pHook->Next = new SC_Cluster(this->pTweak, NULL, NULL, NULL, NULL, pModel, 1, pCorpus->getGT()->getSpeakerIDfromName(tmp, (int)(strlen(tmp))), tmp, true);
						pHook = pHook->Next;
					}
					MFree_1D(tmp);
					modelCnt++;
					MFree_0D(pModel);
				} else {
					REPORT_ERROR(SVLIB_FileErr, "Problem loading one of the models");
					failure = true;
				}
			}
		}	
		fclose(modelList);
	}
	if (pFirst == NULL) {
		REPORT_ERROR(SVLIB_FileErr, "No models could be loaded.\n");
		MFree_0D(pCorpus);
		return false;
	}
	printf("done: %d models loaded!\n", modelCnt);

	timer.startTimer(1);
	MArray_1D(tmp, sclib::bufferSize, char, "SC_MainTasks.testReynoldsSpeakerModels: tmp");
  while (pCorpus->getGT()->getSpeakersCallback(speakerId) == true) {
		printf("\nPhase 1:\tProcessing Audio-Data speaker-wise for speaker nr. %d\n", ++spkCnt);
		segmentCnt = 0;

		while (pCorpus->getGT()->getSpeakersSegmentsCallback(speakerId, segmentStart, segmentEnd) == true) {
      printf("\n    Phase I:\tLoading Signal...");
			pSignal = pCorpus->loadSignal((unsigned long int&)segmentStart, (unsigned long int&)segmentEnd);
			printf("done!");

			if (pSignal != NULL && pSignal->GetLen() > 0) {
				//extract features
				printf("\n    Phase II:\tExtracting Features...");
				pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, sclib::featureMFCC);
				MFree_0D(pSignal);
				printf("done!");

				/*
				//silence detection
				printf("\n    Phase IIIa:\tRemoving Silence...");
				segmenter.silenceDetection(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationHandler.silenceDetectorMode);
				printf("done!");

				//label unvoiced speech
				printf("\n    Phase IIIb:\tLabeling Unvoiced Speech...");
				SV_Data *pPitch = segmenter.unvoicedDetection(pCorpus->getGT(), segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationHandler.vUvDetectorMode);
				MFree_0D(pPitch);
				printf("done!");
				*/

				//copy all speech frames together
				pFeature = pCorpus->getGT()->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureMFCC)], segmentStart, segmentStart, segmentEnd, sclib::atSpeech);
				extractor.freeFeatures(pFeatures);
				if (reduceByPercent > 0.0) { //TODO: ensure reduceByPercent==0.0 for original Reyndolds experiments!
					extractor.reduceByPercent(pFeature, reduceByPercent);
					sclib::scalarOut("audioLength_test.txt", (double)(pCorpus->getGT()->getConverter()->audioFrame2ms(pFeature->Row, pFeature->Hdr.frameSize, pFeature->Hdr.frameStep))/1000.0, this->pTweak); //in [s]
				}

				//test against all the models
				printf("\n    Phase IV:\tTesting models vs. speaker...\n");
				modelCnt = 0;
				pHook = pFirst;
				maxScore = std::numeric_limits<double>::max() * -1.0;

				if (sclib::bitTest(this->pTweak->speakerClusterer.distanceMeasure, sclib::dmGLR) == true) { //TODO: block exists only for signature output
					SC_ModelHandler handler(this->pTweak, false);
					SC_Model *pTestModel = NULL;
					int modelOrder = this->pTweak->modelHandler.maxSpeakerModelOrder;
					if (this->pTweak->modelHandler.orderGuessMode == sclib::guessParameterSearch &&
	  					this->pTweak->modelHandler.foregroundModelType == sclib::mtGMM_new) {
						modelOrder = handler.guessModelOrder(pFeature, NULL, sclib::mtGMM_new, 1, modelOrder, 1);
					}
					pTestModel = handler.createRawModel(this->pTweak->modelHandler.foregroundModelType, NULL, modelOrder, pFeature->Col);
					pTestModel->TrainModel(pFeature);
					SC_Signature *pSig = pTestModel->toSignature();
					if (pSig != NULL) {
						sclib::classOut("signatures.txt", pSig, this->pTweak, ios_base::out|ios_base::app, pCorpus->getGT()->getSpeakerName(speakerId));
					}
					MFree_0D(pSig);
					MFree_0D(pTestModel);
				}

				while (pHook != NULL) {
					timer.startTimer(2);
					pTemp = pHook->getMergedModel()->TestModel(pFeature);
					timer.stopTimer(2);
					//pTemp->Mat[0][0] *= (float)(pFeature->Row); //to undo normalization
					sprintf(tmp, "Speaker: %s\tModel: %s => %f", pCorpus->getGT()->getSpeakerName(speakerId), pHook->getSpeakerName(), pTemp->Mat[0][0]);
					sclib::stringOut("details.txt", tmp, this->pTweak);
					if (pTemp->Mat[0][0] > maxScore) {
						maxScore = pTemp->Mat[0][0];
						maxId = pHook->getSpeakerID();
					}
					MFree_0D(pTemp);
					modelCnt++;
					pHook = pHook->Next;
				}
				tested++;
				if (maxId == speakerId) {
					correct++;
				} else {
					fail++;
				}
				printf("%d recognised as %d; current accuracy: %f", speakerId, maxId, (double)(correct)/(double)(tested));
				MFree_0D(pFeature);
				printf("done!\n");

			} //pSignal valid

			segmentCnt++;
		} //for all segments of one speaker
		printf("\ndone!\n");
	} //for all speakers
	MFree_1D(tmp);
	timer.stopTimer(1);

	printf("\n\nResult: %d correct, %d total, %f accuracy!", correct, tested, (double)(correct)/(double)(tested));
	
	char buf[sclib::bufferSize];
	sprintf(buf, "\n\nResult: %d correct, %d total, %f accuracy (percent: %f%% off)!", correct, tested, (double)(correct)/(double)(tested), reduceByPercent);
	sclib::stringOut("res.txt", buf, this->pTweak);
  
	printf("\nPhase 3:\tCleaning up...");
  MFree_0D(pCorpus);
	pHook = pFirst;
	while (pHook != NULL) {
		pFirst = pHook;
		pHook = pHook->Next;
		MFree_0D(pFirst);
	}

	timer.logIt("timer_test.txt");
	printf("done!");
	printf("\nSimon says: \"Good bye\".\n");

	return !failure;
}

//====================================================================================================================
//	testcontainer for Jun's HHT-code
//====================================================================================================================
//bool SC_MainTasks::HHTtest(unsigned long int segStart, unsigned long int segEnd) {
bool SC_MainTasks::HHTtest(double SampleRate){
  bool finalResult = true;
  unsigned long int segStart = 0;
  SC_SignalHandler *pSighandler = new SC_SignalHandler(pTweak,sclib::stWave);
  SC_HHT *pHHT = new SC_HHT(pTweak);

  SC_Corpus_MPEG7 *pCorpus = new SC_Corpus_MPEG7(this->pTweak, SampleRate, this->audioFile);
  unsigned long int segEnd = pCorpus->getGT()->getAudioSampleCount()-1;

  //load the Signal for analysis
  SC_Signal *pSignal = pCorpus->loadSignal(segStart,segEnd , true);

  //get all IMFs and Residue from original Signal after EMD
  SV_Data *imf = pHHT->emd(pSignal->GetBuf_L(), pSignal->GetLen());

//------------------
//use this code fragment to merge IMFs to 2 parts(high-frequecy und low-frequency)
//------------------
/*
  SV_Data *n = new SV_Data(1,imf->Col);
  SV_Data *v = new SV_Data(1,imf->Col);

  for (long int j=0; j<pSignal->GetLen();j++){
	n->Mat[0][j] = 0.0;
	v->Mat[0][j] = 0.0;}

  for (long int j=0; j<imf->Col;j++)
	for (long int i=0;i<2;i++)
	  n->Mat[0][j] += imf->Mat[i][j];

  for (long int j=0; j<imf->Col;j++)
	for (long int i=2;i<imf->Row;i++)
	  v->Mat[0][j] += imf->Mat[i][j];

  pSighandler->storeSignal("musik.wav", n, pCorpus->getGT()->getSignalPrototype(), -1);
  MFree_0D(n);

  pSighandler->storeSignal("voice.wav", v, pCorpus->getGT()->getSignalPrototype(), -1);
  MFree_0D(v);
 */

  
  //copy the IMFs and store to several wave Data
  SV_Data *imf_wav = new SV_Data(*imf, false); //really create a copy!!!
  pSighandler->storeSignal("imf.wav", imf_wav, pCorpus->getGT()->getSignalPrototype(), -1);
  MFree_0D(imf_wav);
  

  //*fre is instantaneous frequency based on Hz(per second)
	SV_Data *F=new SV_Data(imf->Row,imf->Col);
	SV_Data *A=new SV_Data(imf->Row,imf->Col);

	SV_Data *f_frame=pHHT->insten_frequency(imf,imf->Row,imf->Col,videoFrameRate,F,A);

	SC_MatrixFunctions mFunc;
	
	float **newMat = mFunc.transpose(F->Mat,F->Row, F->Col, false);
	sclib::matrixOut("f_sample.txt",newMat, F->Col,F->Row, this->pTweak);
	MFree_2D(newMat);

	
	float **newMat2 = mFunc.transpose(f_frame->Mat,f_frame->Row, f_frame->Col, false);
	sclib::matrixOut("f_frame.txt",newMat2, f_frame->Col,f_frame->Row, this->pTweak);
	MFree_2D(newMat2);

  //spec is a N*M matrix, where N is the maximal instantaneous frequency every 10Hz , M is the time interval(samples)
  //
  //for example
  //the value of spec[0][j] shows the amplitude of the time point j(sample) with 0 to 10Hz Frequency
  //the value of spec[5][k] shows the amplitude of the time point k(sample) with 50 to 60Hz Frequency
  //...
  //the value of spec[100][l] shows the amplitude of the time point l(sample) with 1000Hz and higher Frequency
	SV_Data *spec = pHHT->hilbert_spectrum(F,A,imf,F->Row,F->Col);
	sclib::matrixOut("spec.txt",spec->Mat, spec->Row,spec->Col, this->pTweak);

	MFree_0D(imf);
	MFree_0D(f_frame);
	MFree_0D(spec);
	MFree_0D(pSighandler);
	MFree_0D(pHHT);
	MFree_0D(pSignal);
	MFree_0D(F);
	MFree_0D(A);

	return finalResult;
}

//====================================================================================================================
//	resynthesizes the speaker-clustering-relevant audio (feature2wav & model2wav) in the given corpus
//====================================================================================================================
bool SC_MainTasks::synthesizeCorpus() {
  bool finalResult = true;
	SC_Corpus *pCorpus;
	SC_Signal *pSignal;
	long int sceneStart, sceneEnd;
	SV_Data **pFeatures, *pConcat, *pTmp, *pMFCC, *pPitch;
	SC_FeatureHandler extractor(this->pTweak, true);
	SC_ModelHandler modeller(this->pTweak, true);
	SC_Synthesis synthesizer(this->pTweak, this->pTweak->featureMfcc.fftSize, 1.0, this->pTweak->featureMfcc.window);
	SC_Model *pModel;
	std::map<unsigned long int, std::pair<unsigned int, unsigned int> > featureMap;
	char *fileName;
	SC_MatrixFunctions matFunc;

	//this works only with MFCCs, possibly accompanied by pitch
	if (!((this->pTweak->modelHandler.speakerModelFeature == sclib::featureMFCC) ||
			  (this->pTweak->modelHandler.speakerModelFeature == (sclib::featureMFCC|sclib::featurePitch)))) {
		REPORT_ERROR(SVLIB_BadArg, "corpus synthesizing only possible for MFCC and/or pitch features");
	}
	this->pTweak->featureMfcc.addDeltaDeltas = false;
	this->pTweak->featureMfcc.addDeltas = false;
	this->pTweak->featureMfcc.coeffSelection[0] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[1] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[2] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[3] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[4] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[5] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[6] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[7] = 0xFF;

	pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, true, false, 0, false); //only set scene-bounds @ file-bounds

	for (unsigned long int y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
    pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
		if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
			pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
			pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, this->pTweak->modelHandler.speakerModelFeature);
			extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list

			//resynthesize extracted features (MFCC & pitch)
			fileName = sclib::addPostfixToFilename(pCorpus->getCurrentAudioFileName(sceneStart), "_sf");
			synthesizer.mfcc2wav(fileName, pFeatures[sclib::bitPosition(sclib::featureMFCC)], this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 0, 1.0, pFeatures[sclib::bitPosition(sclib::featurePitch)]);
			MFree_1D(fileName);

			//do modeling
			extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
			pConcat = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), featureMap, this->pTweak->modelHandler.speakerModelFeature);
			pTmp = pCorpus->getGT()->copyFramesTogether(pConcat, sceneStart, sceneStart, sceneEnd, sclib::atSpeech);
			printf("\n(%#6.5f%% outliers removed from data)", extractor.removeOutliers(pConcat, this->pTweak->modelHandler.outlierRemovalMode));
			pModel = modeller.buildModel(pTmp, NULL, sclib::modeForeground, 1);
			MFree_0D(pTmp);

			//resynthesize modelled speech
			pTmp = pModel->drawSamplesFromDistribution(pConcat->Row);
			pTmp->Hdr = pConcat->Hdr;
			MFree_0D(pModel);
			MFree_0D(pConcat);
			if (featureMap.count(sclib::featureMFCC) > 0) {
				pMFCC = extractor.getSubSet(pTmp, featureMap[sclib::featureMFCC].first, featureMap[sclib::featureMFCC].second);
				pMFCC->Hdr = pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr;
			} else {
				pMFCC = NULL;
			}
			if (featureMap.count(sclib::featurePitch) > 0) {
				pPitch = extractor.getSubSet(pTmp, featureMap[sclib::featurePitch].first, featureMap[sclib::featurePitch].second);
				pPitch->Hdr = pFeatures[sclib::bitPosition(sclib::featurePitch)]->Hdr;
			} else {
				pPitch = NULL;
			}
			MFree_0D(pTmp);
			fileName = sclib::addPostfixToFilename(pCorpus->getCurrentAudioFileName(sceneStart), "_sm");
			if (pMFCC != NULL) {
 				synthesizer.mfcc2wav(fileName, pMFCC, this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 0, 1.0, pPitch);
			} else {
				synthesizer.generateHumFromPitch(fileName, pPitch);
			}
			MFree_1D(fileName);
			MFree_0D(pMFCC);
			MFree_0D(pPitch);

			extractor.freeFeatures(pFeatures);
			MFree_0D(pSignal);

			y = sceneEnd;
		} else {
			break;
		}
	}

	MFree_0D(pCorpus);

	return finalResult;
}

//====================================================================================================================
//	test task for a webservice: gets a wav-file (name) as input, builds a GMM from all samples inside and returns a
//  new wav file name of the resynthesized, inverted GMM
//====================================================================================================================
bool SC_MainTasks::wav2gmm(char* &resultFileName) {
  bool finalResult = true;
	SC_Corpus *pCorpus;
	SC_Signal *pSignal;
	long int sceneStart, sceneEnd;
	SV_Data **pFeatures, *pConcat, *pTmp, *pMFCC, *pPitch;
	SC_FeatureHandler extractor(this->pTweak, true);
	SC_ModelHandler modeller(this->pTweak, true);
	SC_Synthesis synthesizer(this->pTweak, this->pTweak->featureMfcc.fftSize, 1.0, this->pTweak->featureMfcc.window);
	SC_Model *pModel;
	std::map<unsigned long int, std::pair<unsigned int, unsigned int> > featureMap;
	char *fileName;
	SC_MatrixFunctions matFunc;

	//set some parameters
	if (!((this->pTweak->modelHandler.speakerModelFeature == sclib::featureMFCC) || //this works only with MFCCs, possibly accompanied by pitch
			  (this->pTweak->modelHandler.speakerModelFeature == (sclib::featureMFCC|sclib::featurePitch)))) {
		REPORT_ERROR(SVLIB_BadArg, "wav2gmm is only possible for MFCC and/or pitch features");
		this->pTweak->modelHandler.speakerModelFeature = sclib::featureMFCC;
	}
	this->pTweak->featureMfcc.addDeltaDeltas = false;
	this->pTweak->featureMfcc.addDeltas = false;
	this->pTweak->featureMfcc.coeffSelection[0] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[1] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[2] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[3] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[4] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[5] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[6] = 0xFF;
	this->pTweak->featureMfcc.coeffSelection[7] = 0xFF;
	//this->pTweak->modelHandler.foregroundModelType = sclib::mtGMM_new;

	//load the signal
	pCorpus = new SC_Corpus_MPEG7(this->pTweak, 1.0, this->audioFile);
	sceneStart = 0;
	sceneEnd = pCorpus->getGT()->getAudioSampleCount()-1;
	pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
	pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, this->pTweak->modelHandler.speakerModelFeature);
	extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list

	//do modeling
	extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
	pConcat = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), featureMap, this->pTweak->modelHandler.speakerModelFeature);
	printf("\n(%#6.5f%% outliers removed from data)", extractor.removeOutliers(pConcat, this->pTweak->modelHandler.outlierRemovalMode));
	pModel = modeller.buildModel(pConcat, NULL, sclib::modeForeground, 1);

	//resynthesize modelled speech
	pTmp = pModel->drawSamplesFromDistribution(pConcat->Row);
	pTmp->Hdr = pConcat->Hdr;
	MFree_0D(pModel);
	MFree_0D(pConcat);
	if (featureMap.count(sclib::featureMFCC) > 0) {
		pMFCC = extractor.getSubSet(pTmp, featureMap[sclib::featureMFCC].first, featureMap[sclib::featureMFCC].second);
		pMFCC->Hdr = pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr;
	} else {
		pMFCC = NULL;
	}
	if (featureMap.count(sclib::featurePitch) > 0) {
		pPitch = extractor.getSubSet(pTmp, featureMap[sclib::featurePitch].first, featureMap[sclib::featurePitch].second);
		pPitch->Hdr = pFeatures[sclib::bitPosition(sclib::featurePitch)]->Hdr;
	} else {
		pPitch = NULL;
	}
	MFree_0D(pTmp);
	fileName = sclib::addPostfixToFilename(this->audioFile, "_sm");
	if (pMFCC != NULL) {
		synthesizer.mfcc2wav(fileName, pMFCC, this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 0, 1.0, pPitch);
	} else {
		synthesizer.generateHumFromPitch(fileName, pPitch);
	}
	MFree_1D(resultFileName);
	resultFileName = fileName;
	MFree_0D(pMFCC);
	MFree_0D(pPitch);

	extractor.freeFeatures(pFeatures);
	MFree_0D(pSignal);
	MFree_0D(pCorpus);

	return finalResult;	
}

//====================================================================================================================
//	all test stuff for the time model
//====================================================================================================================
bool SC_MainTasks::timeModelTest() {
	bool finalResult = true;
  long int sceneStart, sceneEnd;
  time_t startTime	= time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Model *pModel = NULL;
  SC_FeatureHandler *pExtractor;
  SC_ModelHandler *pModeller;
  SC_Signal *pSignal;
  SV_Data **pFeatures, *pFeatureHook = NULL, *pFeatureFirst = NULL, *pFeature = NULL, *pTmpFeatures = NULL, *pNorm = NULL;

  printf("\nTestframework for time models");
  printf("\n===============================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  pModeller = new SC_ModelHandler(this->pTweak);
  pExtractor = new SC_FeatureHandler(this->pTweak);
	
  //decide which corpus to use:
  char *extension = sclib::extractExtension(this->audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) == 0) { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile);
	} else if (strncmp(this->segmentationFile, "", sclib::bufferSize) == 0) { //no gt => MPEG7 corpus
    pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	} else { //no corpus-file, so assume SCiVo type data
    pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile);
  }
  MFree_1D(extension);
	printf("done!\n");

	//load normalization-matrix, if exists
	pNorm = pExtractor->loadFeature(this->pTweak->modelTime.normalizationFile);

	printf("\nPhase 1:\tProcessing Audio-Data ...\n");
	//load signal
	printf("\n    Phase I:\tLoading Signal...");
	sceneStart = 0;
	sceneEnd = pCorpus->getGT()->getAudioSampleCount() - 1;
	pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd); //sceneStart and -end may get altered here...
	printf("done!");

	if (pSignal != NULL && pSignal->GetLen() > 0) {
		//extract features
		printf("\n    Phase II:\tExtracting Features...");
		pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, this->pTweak->modelHandler.speakerModelFeature);
		MFree_0D(pSignal);
		printf("done!");

		//maybe concatenate multiple feature-sets (a normlization-matrix must be given, then)
		pExtractor->prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
		pExtractor->equalizeFrameParameters(pFeatures, pExtractor->getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
		pTmpFeatures = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
		if (pNorm != NULL)  {
			pExtractor->normalize(pTmpFeatures, pNorm);
		}
		pExtractor->freeFeatures(pFeatures);

		//copy desired frames together and free original feature-sets
		pFeature = pCorpus->getGT()->copyFramesTogether(pTmpFeatures, sceneStart, sceneStart, sceneEnd, sclib::atSpeech, sclib::noType);
		MFree_0D(pTmpFeatures);
		
		//this->pTweak->classifierSvm.gamma = 0.001;
		//this->pTweak->classifierSvm.nu = 0.2;
		//this->pTweak->modelSvm.doParameterSearch = false;

		//pModel = new SC_Model_Time(this->pTweak, pFeature->Col, this->pT, 1, sclib::mtSVM, NULL, true); 
		//pModel = new SC_Model_SVM(this->pTweak, false, false, true);
		//pModel->TrainModel(pFeature);
		
		MFree_0D(pFeature);
		pFeature = pExtractor->loadAndMergeFeatures("../data/test/timit_train.mfcc_pitch");
		pModel = pModeller->loadModel(this->pTweak->modelTime.worldModelFile, sclib::mtSVM);
		((SC_Model_SVM*)pModel)->setDistanceBasedTesting(false);

		/*int *labels = NULL;
		MFree_0D(pFeature);
		pFeature = pExtractor->loadLibSVMfeatures("../data/test/svmtoy.txt", labels);
		MFree_1D(labels);
		pModel = new SC_Model_SVM(this->pTweak, false, true, true);
		pModel->TrainModel(pFeature);*/

		SV_Data *pScore = pModel->TestModel(pFeature);
		printf("result: %f", pScore->Mat[0][0]);
		MFree_0D(pModel);
		MFree_0D(pFeature);
		MFree_0D(pScore);
	}

	printf("\nPhase 3:\tCleaning up...");

	MFree_0D(pNorm);
  MFree_0D(pCorpus);
  MFree_0D(pModeller);
  MFree_0D(pExtractor);

	printf("done!");
	return finalResult;
}

//====================================================================================================================
//	test container to determine length of voice-dependent speech "sounds"
//====================================================================================================================
bool SC_MainTasks::splicer(int blockSizeMs, int olaMaxIterations, double olaErrorTarget) {
	SC_Corpus *pCorpus;
	SC_Signal *pSignal;
	long int sceneStart, sceneEnd;
	SV_Data **pFeatures, *pFeature;
	SC_FeatureHandler extractor(this->pTweak, true);
	SC_Synthesis synthesizer(this->pTweak, this->pTweak->featureMfcc.fftSize, 1.0, pTweak->featureMfcc.window);
	char *fileName, *newFileName, buffer[sclib::bufferSize];
	bool res = true;
	int preservedBlockSize;

	//set some parameters
	//this->pTweak->featureSpectrum.frameSize = 32;
	//this->pTweak->featureSpectrum.frameStep = 16;
	this->pTweak->featureSpectrum.FFTsize = 512;
	this->pTweak->featureSpectrum.preEmphasizeFactor = 0.0;
	this->pTweak->featureSpectrum.lowCut = 0.0;
	this->pTweak->featureSpectrum.highCut = 0.0;
	this->pTweak->featureSpectrum.createPhase = false;
	this->pTweak->featureSpectrum.logarithmize = true;
	this->pTweak->featureSamples.frameSize = this->pTweak->featureSpectrum.frameSize;
	this->pTweak->featureSamples.frameStep = this->pTweak->featureSpectrum.frameStep;
	this->pTweak->featureSamples.highCut = this->pTweak->featureSpectrum.highCut;
	this->pTweak->featureSamples.lowCut = this->pTweak->featureSpectrum.lowCut;
	synthesizer.setFftLength(this->pTweak->featureSpectrum.FFTsize);
	synthesizer.setOlaMaxIterations(olaMaxIterations);
	synthesizer.setOlaErrorTarget(olaErrorTarget);
	synthesizer.setTaperingLength(this->pTweak->transform.taperingLength);
	synthesizer.setTaperingMode(this->pTweak->featureSpectrum.window!=0 ? sclib::wndHamming : sclib::wndRectangle);

	//load the signal
	pCorpus = new SC_Corpus_MPEG7(this->pTweak, 1.0, this->audioFile, "");
	sceneStart = 0;
	sceneEnd = pCorpus->getGT()->getAudioSampleCount()-1;
	pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd);
	pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, sclib::featureSamples);
	extractor.prepareFeatureSet(pFeatures, NULL);
	pFeature = pFeatures[sclib::bitPosition(sclib::featureSamples)];

	//splice & blend
	preservedBlockSize = pCorpus->getGT()->getConverter()->ms2audioFrame(blockSizeMs, pFeatures[sclib::bitPosition(sclib::featureSamples)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureSamples)]->Hdr.frameStep);
	extractor.splice(pFeature, preservedBlockSize); //randomize time-order of frames
	//pTmp = extractor.blend(pFeatures[sclib::bitPosition(featureType)], intermediateFrameCount, steepness, preservedBlockSize);

	//filter out low-energy blocks
	/*
	int blockCount = sclib::getRowCount(pFeature->Row, preservedBlockSize, preservedBlockSize);
	int t, d, blockNr = 0;
	SC_MatrixFunctions matFunc;
	SV_Data *pNrg = new SV_Data(blockCount, 1);
	matFunc.fillWithValue(pNrg->Mat, pNrg->Row, pNrg->Col, 0.0f);
	for (t = 0; t < pFeature->Row; t++) { //get avg. energy of each block
		if (t%preservedBlockSize==0 && t>0 && blockNr<blockCount-1) {
			blockNr++;
		}
		for (d = 0; d < pFeature->Col; d++) {
			pNrg->Mat[blockNr][0] += (float)(fabs(pFeature->Mat[t][d]) / (double)(preservedBlockSize*pFeature->Col));
		}
	}
	SC_MixtureModel_GMM *pGMM = new SC_MixtureModel_GMM(this->pTweak, 2, 1); //cluster the avg. nrg's via an gmm
	pGMM->TrainModel(pNrg);
	SV_Data *pScore = pGMM->TestModel(pNrg);
	int highNrgMixture = (pGMM->getMean(0, 0) > pGMM->getMean(1, 0)) ? 0 : 1; //which of the 2 mixtures represents the high energy cluster?
	int highNrgCount = 0;
	for (t = 1; t <= blockCount; t++) { //how many high energy frames are there
		if (pScore->Mat[t][0] == (float)(highNrgMixture)) {
			highNrgCount++;
		}
	}
	SV_Data *pTmp = new SV_Data(highNrgCount*preservedBlockSize, pFeature->Col);
	int tt = 0;
	blockNr = 0;
	for (t = 0; t < pFeature->Row; t++) {
		if (t%preservedBlockSize==0 && t>0 && blockNr<blockCount-1) {
			blockNr++;
		}
		if (pScore->Mat[blockNr+1][0]==(float)(highNrgMixture) && tt<pTmp->Row) {
			for (d = 0; d < pFeature->Col; d++) {
				pTmp->Mat[tt][d] = pFeature->Mat[t][d];
			}
			tt++;
		}
	}
	MFree_0D(pScore);
	MFree_0D(pGMM);
	pTmp->Hdr = pFeature->Hdr;
	pFeature = pTmp;
	*/

	int blockCount=sclib::getRowCount(pFeature->Row, preservedBlockSize, preservedBlockSize);
	for (int b = 0; b < blockCount; b++) {
		SV_Data *pTmp = new SV_Data(preservedBlockSize, pFeature->Col);
		pTmp->Hdr = pFeature->Hdr;
		int tt = 0;
		for (int t = b*preservedBlockSize; t < (b+1)*preservedBlockSize; t++) {
			for (int d = 0; d < pFeature->Col; d++) {
				pTmp->Mat[tt][d] = pFeature->Mat[t][d];
			}
			tt++;
		}		
		fileName = sclib::extractFileName(this->audioFile);
		sprintf(buffer, "_%dms_%d_resyn.wav", blockSizeMs, b);
		newFileName = sclib::exchangeFileExtension(fileName, buffer);
		sprintf(buffer, "%s%s", this->pTweak->debug.debugDir, newFileName);
		synthesizer.samples2wav(buffer, pTmp);
		MFree_1D(fileName);
		MFree_1D(newFileName);
		MFree_0D(pTmp);
	}

	/*
	//resynthesize
	fileName = sclib::extractFileName(this->audioFile);
	sprintf(buffer, "_%dms_resyn.wav", blockSizeMs);
	newFileName = sclib::exchangeFileExtension(fileName, buffer);
	sprintf(buffer, "%s%s", this->pTweak->debug.debugDir, newFileName);
	synthesizer.samples2wav(buffer, pFeature);
	MFree_1D(fileName);
	MFree_1D(newFileName);
	*/

	extractor.freeFeatures(pFeatures);
	//MFree_0D(pTmp);
	MFree_0D(pSignal);
	MFree_0D(pCorpus);

	return res;
}

//====================================================================================================================
//	extracts feature(s) for given file; skips sampleFeed samples at the beginning of the signal in order to evaluate 
//  what effect small displacements in frame position may have on the feature vectors.
//  switch off pre-emphasis, resampling and low-/highpass-filtering for the choosen (signle) feature-set (via the 
//  speaker-model-feature in the modelHandlers tweakable parameters) in order to get comparable results 
//====================================================================================================================
bool SC_MainTasks::featureDisplacemetTest(int sampleFeed, bool filePerFrame) {
	bool finalResult = true;
	SC_Corpus *pCorpus;
	SV_Data **pFeatures = NULL, *pFeature;
	SC_FeatureHandler extractor(this->pTweak, true);
	SC_Signal *pSignal = NULL;
	unsigned long int segmentStart, segmentEnd;
	char *fileName, *withoutPath, tmp[sclib::bufferSize];

  printf("\nThis is a script to analyze the effect of frame position displacement");
  printf("\n======================================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  //a file with pure silence from start to end is expected, hence no groun-truth is necessary, thus use the MPEG7-corpus (it is gt-free)
	pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	segmentStart = 0 + sampleFeed;
	segmentEnd = pCorpus->getGT()->getAudioSampleCount() - 1;
	printf("done!\n");

	printf("\nPhase 1:\tExtracting features...");
	pSignal = pCorpus->loadSignal(segmentStart, segmentEnd, true);
	pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, this->pTweak->modelHandler.speakerModelFeature);
	pFeature = pFeatures[sclib::bitPosition(this->pTweak->modelHandler.speakerModelFeature)];
	printf("done!\n");

	printf("\nPhase 1:\tWriting them to one file per frame...");
	withoutPath = sclib::extractFileName(this->audioFile);
	if (filePerFrame == true) {
		for (int t = 0; t < pFeature->Row; t++) {
			sprintf(tmp, "_frame_%d.txt", t);
			fileName = sclib::exchangeFileExtension(withoutPath, tmp);
			sclib::vectorOut(fileName, pFeature->Mat[t], pFeature->Col, false, this->pTweak, " ");
			MFree_1D(fileName);
		}
	} else {
		fileName = sclib::exchangeFileExtension(withoutPath, "_features.txt");
		sprintf(tmp, "displacement of %d samples\n", sampleFeed);
		sclib::classOut(fileName, pFeature, this->pTweak, ios_base::out|ios_base::app, tmp);
		MFree_1D(fileName);
	}
	MFree_1D(withoutPath)
	printf("done!\n");

	printf("\nPhase 2:\tCleaning up...");
	extractor.freeFeatures(pFeatures);
	MFree_0D(pSignal);
  MFree_0D(pCorpus);
	printf("done!");
	printf("\n\nThe test script says: \"Good bye\".\n");

	return finalResult;
}

//====================================================================================================================
//	receives the filename of a saved feature set as "audioFile", loads it and then performs k-means clustering and 
//  saves the found centroids
//====================================================================================================================
bool SC_MainTasks::createTemplates(int k) {
	bool finalResult = true;
  SV_Data *pFeature = NULL, *pCentroids = NULL;
	SC_Clusterer clusterer(this->pTweak);
  SC_FeatureHandler extractor(this->pTweak, true);
	char fileName[sclib::bufferSize];

  printf("\nCreating VQ templates for a given feature set");
  printf("\n==============================================\n");

  if (sclib::fileExists(this->audioFile) == true) {
		printf("\nLoading features...");
    pFeature = extractor.loadAndMergeFeatures(this->audioFile);

		printf("\nPerforming clustering...");
		pCentroids = clusterer.kMeans(k, pFeature, 10, false);

		sprintf(fileName, "%s.vq_%d", this->audioFile, k);
		extractor.saveFeature(fileName, pCentroids);

		MFree_0D(pFeature);
		MFree_0D(pCentroids);
  } else {
    printf("\nno data!\n");
		finalResult = false;
  }

	printf("\nCleaning up...");
	printf("done!");

	printf("\n\ncreateTemplates says: \"Good bye\".\n");

	return finalResult;
}

//====================================================================================================================
//	creates the context vectors (trajectories) for each speaker in the given training set and creates a matrix of 
//  those vectors (along with a file specifiying the speaker id for each column)
//====================================================================================================================
bool SC_MainTasks::createContextVectorMatrix(char *matrixFileName, char *idFileName) {
	bool finalResult = true;
	long int t, sceneStart, sceneEnd, speakerId;
  time_t startTime	= time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Model *pModel = NULL;
  SC_FeatureHandler extractor(this->pTweak, true);
  SC_ModelHandler handler(this->pTweak, true);
  SC_Signal *pSignal;
	char buffer[sclib::bufferSize];
  SV_Data **pFeatures, *pFeatureHook = NULL, *pFeatureFirst = NULL, *pFeature = NULL, *pTmpFeatures = NULL, *pNorm = NULL;
	unsigned long int trajectoryLength;
	std::vector<long int> idMap;

  printf("\nTestframework for time models");
  printf("\n===============================\n");

	printf("\nPhase 0:\tLoading Parameters...");
	pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, true, false, 0, true);
	pNorm = extractor.loadFeature(this->pTweak->modelTime.normalizationFile);
	sprintf(buffer, "%stmpTrajectories.dat", this->pTweak->debug.debugDir);
	printf("done!\n");
		
	while (pCorpus->getGT()->getSpeakersCallback(speakerId) == true) {
		printf("\nPhase 1:\tProcessing Audio-Data speaker-wise for speaker %d\n", speakerId);

		while (pCorpus->getGT()->getSpeakersSegmentsCallback(speakerId, sceneStart, sceneEnd) == true) {
      printf("\n    Phase I:\tLoading signal...");
			pSignal = pCorpus->loadSignal((unsigned long int&)sceneStart, (unsigned long int&)sceneEnd);

			if (pSignal != NULL && pSignal->GetLen() > 0) {
	      printf("\n    Phase II:\tExtracting features...");
				pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, this->pTweak->modelHandler.speakerModelFeature);
				extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
				SC_SignalHandler sigHandler(this->pTweak);
				sigHandler.saveSignal("train.wav", pSignal, true);
				MFree_0D(pSignal);

	      printf("\n    Phase II:\tForming features...");
				if (this->pTweak->modelHandler.speakerModelFeature == sclib::featureSDP) {
					pFeature = pFeatures[sclib::bitPosition(sclib::featureSDP)];
					SC_MatrixFunctions matFunc;
					//sclib::classOut("training.txt", pFeature, this->pTweak);
					matFunc.addNoise(pFeature->Mat, pFeature->Row, pFeature->Col, 1.0f);
				} else {
					//maybe concatenate multiple feature-sets (a normlization-matrix must be given, then)
					extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
					pTmpFeatures = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
					if (pNorm != NULL)  {
						extractor.normalize(pTmpFeatures, pNorm);
					}
					trajectoryLength = pCorpus->getGT()->getConverter()->ms2audioFrame(this->pTweak->modelTime.syllableLength, pTmpFeatures->Hdr.frameSize, pTmpFeatures->Hdr.frameStep, sclib::alignmentStart);
					pFeature = extractor.createTrajectories(pTmpFeatures, trajectoryLength, this->pTweak->modelTime.trajectoryStep, this->pTweak->modelTime.removeTiming, this->pTweak->modelTime.templateCount, this->pTweak->modelTime.clusteringIterations);
					MFree_0D(pTmpFeatures);
				}

	      printf("\n    Phase III:\tSaving temporaries...");
				if (!extractor.saveFeature(buffer, pFeature)) {
					printf("failed!");
				}
				for (t = 0; t < pFeature->Row; t++) { //fill the dictionary of which column later contains samples from which speaker
					idMap.push_back(speakerId);
				}				
				MFree_0D(pFeature);
				extractor.freeFeatures(pFeatures);
			}
		}
	}

	printf("\nPhase 3:\tCreating ground truth...");
	remove(idFileName);
	sclib::mapOut(idFileName, &idMap);
	idMap.clear();

	printf("\nPhase 4:\tCreating final matrix...");
	pFeature = extractor.loadAndMergeFeatures(buffer, -1, -1, true);
	remove(buffer);
	extractor.normalizeColumns(pFeature);

	printf("\nPhase 5:\tSaving final matrix...");
	remove(matrixFileName);
	sclib::matrixOut(matrixFileName, pFeature->Mat, pFeature->Row, pFeature->Col);
	//if (!extractor.saveFeature(matrixFileName, pFeature)) {
	//	printf("failed!");
	//	finalResult = false;
	//}
	
	printf("\nPhase 6:\tCleaning up...");
	MFree_0D(pNorm);
  MFree_0D(pCorpus);

	printf("done!");
	return finalResult;
}

//====================================================================================================================
//	finds a sparse representation of the audioFile's context vectors against the trained matrix via Matlab's l1-Magic
//  and evaluates it
//====================================================================================================================
bool SC_MainTasks::testContextVectorMatrix(char *matrixFileName, char *idFileName) {
	bool finalResult = true;
	long int t, sceneStart, sceneEnd;
  time_t startTime	= time(NULL);
  SC_Corpus *pCorpus = NULL;
  SC_Model *pModel = NULL;
  SC_FeatureHandler extractor(this->pTweak, true);
  SC_ModelHandler handler(this->pTweak, true);
  SC_Signal *pSignal;
  SV_Data **pFeatures, *pFeatureHook = NULL, *pFeatureFirst = NULL, *pFeature = NULL, *pTmpFeatures = NULL, *pNorm = NULL, *pTrain = NULL;
	unsigned long int trajectoryLength;
	int *idMap = NULL;
	SC_Matlab mHandler;
	SC_MatrixFunctions matFunc;

  printf("\nTestframework for time models");
  printf("\n===============================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  //decide which corpus to use:
  char *extension = sclib::extractExtension(this->audioFile);
  if (strncmp(extension, "crp", sclib::bufferSize) == 0) { //a corpus-file (extension ".crp") as audioFile, so assume the TIMIT corpus as data source
    pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, false, false, 0, true);
	} else if (strncmp(this->segmentationFile, "", sclib::bufferSize) == 0) { //no gt => MPEG7 corpus
    pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	} else { //no corpus-file, so assume SCiVo type data
    pCorpus = new SC_Corpus_SCiVo(this->pTweak, this->videoFrameRate, this->audioFile, this->sceneFile, this->segmentationFile);
  }
  MFree_1D(extension);
	pNorm = extractor.loadFeature(this->pTweak->modelTime.normalizationFile);
	pTrain = extractor.loadFeature(matrixFileName);
	mHandler.initialize();
	MatlabRealMatrix mTrain = mHandler.matrix2matlab(pTrain->Mat, pTrain->Row, pTrain->Col);
	int trainCols = pTrain->Col;
	MFree_0D(pTrain);
	MArray_1D(idMap, trainCols, int, "SC_MainTasks.testContextVectorMatrix: idMap");
	ifstream idFile;
	idFile.open(idFileName);
	for (t = 0; t < trainCols; t++) {
		idFile >> idMap[t];
	}
	idFile.close();
	printf("done!\n");

	printf("\nPhase 1:\tProcessing Audio-Data ...\n");
	for (unsigned long int y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
		pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
		if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {

			printf("\n    Phase I:\tLoading signal...");
			pSignal = pCorpus->loadSignal((unsigned long int&)sceneStart, (unsigned long int&)sceneEnd);
			if (pSignal != NULL && pSignal->GetLen() > 0) {
				printf("\n    Phase II:\tExtracting features...");
				pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, this->pTweak->modelHandler.speakerModelFeature);
				extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
				//SC_SignalHandler sigHandler(this->pTweak);
				//sigHandler.saveSignal("test.wav", pSignal, true);
				MFree_0D(pSignal);

				printf("\n    Phase II:\tForming features...");
				if (this->pTweak->modelHandler.speakerModelFeature == sclib::featureSDP) {
					pFeature = pFeatures[sclib::bitPosition(sclib::featureSDP)];
					sclib::classOut("testing.txt", pFeature, this->pTweak);
					//matFunc.addNoise(pFeature->Mat, pFeature->Row, pFeature->Col, 0.3f);
				} else {
					//maybe concatenate multiple feature-sets (a normalization-matrix must be given, then)
					extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
					extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
					pTmpFeatures = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
					if (pNorm != NULL)  {
						extractor.normalize(pTmpFeatures, pNorm);
					}
					trajectoryLength = pCorpus->getGT()->getConverter()->ms2audioFrame(this->pTweak->modelTime.syllableLength, pTmpFeatures->Hdr.frameSize, pTmpFeatures->Hdr.frameStep, sclib::alignmentStart);
					pFeature = extractor.createTrajectories(pTmpFeatures, trajectoryLength, this->pTweak->modelTime.trajectoryStep, this->pTweak->modelTime.removeTiming, this->pTweak->modelTime.templateCount, this->pTweak->modelTime.clusteringIterations);
					MFree_0D(pTmpFeatures);
				}

				printf("\n    Phase II:\tEvaluating test patterns...");
				extractor.normalizeRows(pFeature);
				double *z = matFunc.zeros(trainCols);
				MatlabRealMatrix epsilon = mHandler.scalar2matlab(0.01);
				MatlabRealMatrix lbtol = mHandler.scalar2matlab(1e-3);
				MatlabRealMatrix mu = mHandler.scalar2matlab(10);
				for (t = 0; t < pFeature->Row; t++) {
					if (t > 0) {
						z[t-1] = 0.0;
					}
					if (t < trainCols) {
						z[t] = 1.0;
					}
					MatlabRealMatrix x0 = mHandler.vector2matlab(z, trainCols, false);
					MatlabRealMatrix example = mHandler.vector2matlab(pFeature->Mat[t], pFeature->Col, false);
					MatlabRealMatrix xp = mHandler.sc_mlfL1qc_logbarrier(x0, mTrain, example, epsilon, lbtol, mu);
					//MatlabRealMatrix xp = mHandler.sc_mlfL1eq_pd(x0, mTrain, example);
					mHandler.freeMatrix(example);
					double *weights = mHandler.matlab2array(xp);
					sclib::vectorOut("spk.txt", weights, trainCols, false, this->pTweak);
					mHandler.freeMatrix(xp);
					mHandler.freeMatrix(x0);
					std::map<int, int> speakers;
					int cnt = 0;
					for (int i = 0; i < trainCols; i++) {
						if (weights[i] > 0.0) {
							cnt++;
							speakers[idMap[i]]++;
						}
					}
					std::map<int, int>::const_iterator j;
					for (j = speakers.begin(); j != speakers.end(); ++j) {
						sclib::scalarOut("speakers.txt", j->first, this->pTweak, false, " ");
						sclib::scalarOut("speakers.txt", j->second, this->pTweak, false, ",");
					}
					sclib::scalarOut("speakers.txt", cnt, this->pTweak, false, "\n");
					MFree_1D(weights);
				}
				MFree_1D(z);
				mHandler.freeMatrix(epsilon);
				mHandler.freeMatrix(lbtol);
				mHandler.freeMatrix(mu);

				extractor.freeFeatures(pFeatures);
				if (this->pTweak->modelHandler.speakerModelFeature != sclib::featureSDP) {
					MFree_0D(pFeature);
				}
			}

			y = sceneEnd;
		} else {
			break;
		}
	}

	printf("\nPhase 3:\tCleaning up...");
	MFree_0D(pNorm);
	mHandler.freeMatrix(mTrain);
  MFree_0D(pCorpus);
	MFree_1D(idMap);

	printf("done!");
	return finalResult;
}

//====================================================================================================================
//	produces synthetic data to build a GMM from and visualize it
//====================================================================================================================
bool SC_MainTasks::gmmVisualizationTest() {
	bool res = true;
	double *w, **m, **v;

	MArray_1D(w, 3, double, "gmmVisualizationTest.w");
	MArray_2D(m, 3, 2, double, "gmmVisualizationTest.w");
	MArray_2D(v, 3, 2, double, "gmmVisualizationTest.w");

	w[0] = 0.2; w[1] = 0.3; w[2]= 0.5;

	m[0][0] = 3.0; m[0][1] = 2.0;
	m[1][0] = 5.0; m[1][1] = 6.5;
	m[2][0] = 9.0; m[2][1] = 4.0;

	v[0][0] = 2.0; v[0][1] = 1.0;
	v[1][0] = 1.0; v[1][1] = 4.0;
	v[2][0] = 0.5; v[2][1] = 0.5;

	SC_MixtureModel_GMM *pGMM = new SC_MixtureModel_GMM(this->pTweak, 3, 2, w, m, v);
	SV_Data *pTrainData = pGMM->drawSamplesFromDistribution(10000);


	pGMM->TrainModel(pTrainData);
	sclib::classOut("gmm.txt", pGMM, this->pTweak);

	MFree_1D(w);
	MFree_2D(m);
	MFree_2D(v);
	MFree_0D(pTrainData);
	MFree_0D(pGMM);

	return res;
}

//should rather be called timeModel2wac()... (by stdm 2017)
bool SC_MainTasks::wav2timeModel(char *resultPath, int syllableLength, int trajectoryStep, char subModelType, bool removeTiming, int templateCount, int clusteringIterations, int frameSize, int frameStep, bool useMFCC, bool addPitch, int mfccOrder, int lpcOrder, double preemphasis, int window, int fftSize, int frequencyScale, int filterbankSize, double minFilterbankFrequency, double maxFilterbankFrequency, double olaErrorTarget, int olaMaxIterationCount, int pitchFrameSize, int pitchFrameStep, double pitchCandidateThresh, double pitchLagWeight, double pitchFrequencyWeight, double pitchTransitionCost, double pitchTransitionAmpModCost, double pitchTransitionSpecModCost, double pitchVoiceBias, double pitchDoublingCost, double pitchMinF0, double pitchMaxF0, int pitchCandidateCount, double pitchCandidateWindowSize) {
	char error[sclib::bufferSize];
	//sclib::errorHandlerThrows(false, true);

	//try {
		SC_Corpus *pCorpus;
		SC_Signal *pSignal;
		long int sceneStart, sceneEnd;
		SV_Data **pFeatures, *pConcat, *pTmp, *pFeature, *pPitch;
		//const char *iniFileName = env->GetStringUTFChars(iniFile, 0);
		//SC_TweakableParameters *pTweak = new SC_TweakableParameters(iniFileName);
		SC_FeatureHandler extractor(pTweak, true);
		SC_Synthesis synthesizer(pTweak, pTweak->featureMfcc.fftSize, 1.0, pTweak->featureMfcc.window);
		SC_Model *pModel;
		std::map<unsigned long int, std::pair<unsigned int, unsigned int> > featureMap;
		char audioFile[sclib::bufferSize];
		const char *resultPathName = resultPath; //env->GetStringUTFChars(resultPath, 0);
		char *path = sclib::makePath(resultPathName);
		bool res;

		//setbuf(stdout, NULL); //get instant printf()s
		//setbuf(stderr, NULL); //get instant printf()s
		//std::ios::sync_with_stdio(); //syncronize couts and printf()s

		//release the strings previously copied from java
		//env->ReleaseStringUTFChars(iniFile, iniFileName);
		//env->ReleaseStringUTFChars(resultPath, resultPathName);

		pCorpus = new SC_Corpus_MPEG7(pTweak, 1.0, this->audioFile, this->sceneFile); //, env, jStream);

		//set some parameters
		pTweak->modelHandler.foregroundModelType = sclib::mtTime;
		pTweak->modelTime.syllableLength = sclib::max(0, syllableLength);
		pTweak->modelTime.trajectoryStep = sclib::getBetween(1, trajectoryStep, 50);
		pTweak->modelTime.subModelType = subModelType;
		pTweak->modelTime.removeTiming = removeTiming;
		pTweak->modelTime.templateCount = sclib::getBetween(0, templateCount, 10000);
		pTweak->modelTime.clusteringIterations = sclib::getBetween(0, clusteringIterations, 1000);
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
		pTweak->featureMfcc.addDeltaDeltas = false; //the following parameters can't be made audible and thus are not changeable by the user
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
		pModel = new SC_Model_Time(pTweak, pConcat->Col, pTweak->modelTime.syllableLength, pTweak->modelTime.trajectoryStep, pTweak->modelTime.subModelType, pTweak->modelTime.removeTiming, pTweak->modelTime.templateCount, pTweak->modelTime.clusteringIterations, pTweak->modelTime.replaceTrainingData, pTweak->modelTime.checkForTrajectorization, pTweak->modelTime.worldModelFile, pTweak->modelTime.normalizationFile, false, true);
		pModel->TrainModel(pConcat, 0);
		
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
			return false; //TODO change in java version if this will ever be built...
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

		return true;

		//return (res==true) ? env->NewStringUTF(audioFile) : env->NewStringUTF(error);	
	//} catch (const char *errStr) {
	//	sprintf(error, "%s%s", "ERROR: ", errStr);

	//	return env->NewStringUTF(error);
	//}
}

//====================================================================================================================
//	audio feature extraction for markus (audio-"sift"-idea)
//====================================================================================================================
bool SC_MainTasks::videoFeatureExtractor(char *fileName, unsigned long int featureTypes) {
	time_t startTime = time(NULL);
	printf("\nVideo feature extraction for Markus ('audio-SIFT-test')");
  printf("\n=========================================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
	SC_Corpus *pCorpus = new SC_Corpus_MPEG7(this->pTweak, this->videoFrameRate, this->audioFile);
	unsigned long int segmentStart = 0;
	unsigned long int segmentEnd = pCorpus->getGT()->getAudioSampleCount() - 1;
	printf("done!\n");

	printf("\nPhase 1:\tExtracting features...");
	SC_Signal *pSignal = pCorpus->loadSignal(segmentStart, segmentEnd, true);
	SC_FeatureHandler extractor(this->pTweak, true);
	SV_Data **pFeatures = extractor.extractFeatures(pCorpus, pSignal, segmentStart, segmentEnd, featureTypes);
	printf("done!\n");

	printf("\nPhase 2:\tEqualizing and output...");
	extractor.prepareFeatureSet(pFeatures, NULL);
	extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), featureTypes);
	SV_Data *pConcat = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), featureTypes);
	char *featureNames = extractor.getFeatureName(featureTypes);
	pCorpus->featureOut(fileName, featureNames, segmentStart, segmentEnd, pConcat, true);
	MFree_0D(featureNames);
	//sclib::classOut(fileName, pConcat, this->pTweak);
	printf("done!\n");

	printf("\nPhase 2:\tCleaning up...");
	extractor.freeFeatures(pFeatures);
	MFree_0D(pSignal);
  MFree_0D(pCorpus);
	printf("done!");

	unsigned long int elapsedSeconds = (unsigned long)(time(NULL) - startTime);
	cout << "Elapsed time: " << setw(5) << (elapsedSeconds / 60) << ":" << setw(2) << (elapsedSeconds % 60) << "\n";

	printf("\n\nThe test script says: \"Good bye\".\n");

	return true;
}

//====================================================================================================================
//	creates the context vectors (trajectories) for each speaker in the given TIMIT corpus file, and saves them in a 
//  new text file under the name fileName
//  uses the timeModel normalization matrix and speakerModelFeature
//====================================================================================================================
bool SC_MainTasks::createTIMITcontextVectors(char *fileName) {
  printf("\ncreateTIMITcontextVectors()");
  printf("\n============================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  SC_FeatureHandler extractor(this->pTweak, true);
	SC_Corpus* pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile, true, true, false, 0, false); //mind the last "false" -> file-wise scene-boudaries!
	SV_Data* pNorm = extractor.loadFeature(this->pTweak->modelTime.normalizationFile);
	char* fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
	sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	remove(fileName);
	fstream	fileOut;
	fileOut.open(fName, ios_base::out|ios_base::app);
	fileOut << "utterance?;speaker-id;file-name;trajectory-vector" << "\n";
	printf("done!\n");

	printf("\nPhase 1:\tProcessing Files...\n");
	long int sceneStart, sceneEnd, sceneNr = 0;
	for (unsigned long int y = 0; y <= pCorpus->getGT()->getAudioSampleCount(); y++) {
		pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
		if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
			printf("\n  Processing Scene %d...", ++sceneNr);

	//while (pCorpus->getGT()->getSpeakersCallback(speakerId) == true) {
		//printf("\nPhase 1:\tProcessing Audio-Data speaker-wise for speaker %d\n", speakerId);
		//while (pCorpus->getGT()->getSpeakersSegmentsCallback(speakerId, sceneStart, sceneEnd) == true) {
      printf("\n    Phase I:\tLoading signal...");
			SC_Signal* pSignal = pCorpus->loadSignal((unsigned long int&)sceneStart, (unsigned long int&)sceneEnd);

			if (pSignal != NULL && pSignal->GetLen() > 0) {
	      printf("\n    Phase II:\tExtracting features...");
				SV_Data** pFeatures = extractor.extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, this->pTweak->modelHandler.speakerModelFeature);
				extractor.prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
				MFree_0D(pSignal);

	      printf("\n    Phase II:\tForming features...");
				//maybe concatenate multiple feature-sets (a normlization-matrix must be given, then)
				extractor.equalizeFrameParameters(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
				SV_Data* pTmpFeatures = extractor.combineFeatureVectors(pFeatures, extractor.getFeatureCount(), this->pTweak->modelHandler.speakerModelFeature);
				if (pNorm != NULL)  {
					extractor.normalize(pTmpFeatures, pNorm);
				}
				//no silence removal takes place here.
				unsigned long int trajectoryLength = pCorpus->getGT()->getConverter()->ms2audioFrame(this->pTweak->modelTime.syllableLength, pTmpFeatures->Hdr.frameSize, pTmpFeatures->Hdr.frameStep, sclib::alignmentStart);
				SV_Data* pFeature = extractor.createTrajectories(pTmpFeatures, trajectoryLength, this->pTweak->modelTime.trajectoryStep, this->pTweak->modelTime.removeTiming, this->pTweak->modelTime.templateCount, this->pTweak->modelTime.clusteringIterations);
				MFree_0D(pTmpFeatures);

	      printf("\n    Phase III:\tSaving features...");
				unsigned long int offset, sampleCount;
				char *wavFileName = ((SC_GroundTruth_TIMIT*)pCorpus->getGT())->whichFile(sceneStart, offset, sampleCount); 
				const char *speakerName = pCorpus->getGT()->getSpeakerName(pCorpus->getGT()->getMajorSpeakerID(sceneStart, sceneEnd, sclib::modeGroundtruth));
				char sex = ((SC_GroundTruth_TIMIT*)pCorpus->getGT())->getGenderFromFilePath(wavFileName);
				unsigned int utteranceNr = ((SC_GroundTruth_TIMIT*)pCorpus->getGT())->getUtteranceNumber(wavFileName, speakerName);
				for (long int t = 0; t < pFeature->Row; t++) {
					fileOut << utteranceNr << ";" << sex << speakerName << ";" << wavFileName; //next ";" comes in the loop over dimensions
					for (long int d = 0; d < pFeature->Col; d++) {
						fileOut << ";" << pFeature->Mat[t][d];
					}
					fileOut << "\n";
				}

				y = sceneEnd;
				MFree_0D(pFeature);
				extractor.freeFeatures(pFeatures);
			}
		}
	}

	printf("\nPhase 6:\tCleaning up...");
	MFree_0D(pNorm);
  MFree_0D(pCorpus);
	fileOut.close();
	MFree_1D(fName);

	printf("done!");
	return true;
}

//====================================================================================================================
//	training the time world model & normalization matrix for ACM MM 2009 experiments, stripped by anything not needed 
//  to work on the used TIMIT data; 
//  uses timeModel normalization matrix and worldModel filenames
//  based on modeller() script
//====================================================================================================================
bool SC_MainTasks::trainTimeWorldModelAndNormMatrix(const char *tmpFeaturesFile) {
	bool finalResult = true;
  time_t startTime = time(NULL);
  SV_Data *pFeature = NULL, *pFeatureHook = NULL, *pFeatureFirst = NULL, *pNorm = NULL;

  printf("\ntrainTimeWorldModelAndNormMatrix()");
  printf("\n===================================\n");

	printf("\nPhase 0:\tLoading Parameters...");
  SC_ModelHandler* pModeller = new SC_ModelHandler(this->pTweak);
  SC_FeatureHandler* pExtractor = new SC_FeatureHandler(this->pTweak);
  SC_Corpus* pCorpus = new SC_Corpus_TIMIT(this->pTweak, this->audioFile);
	printf("done!\n");

	//load normalization-matrix, if exists
	if (sclib::fileExists(pTweak->modelTime.normalizationFile) == true) {
		pNorm = pExtractor->loadFeature(pTweak->modelTime.normalizationFile);
	}

  if (sclib::fileExists(tmpFeaturesFile) == false) {
		printf("\nPhase 1:\tProcessing Audio-Data piece-wise...\n");
		unsigned long int y = 0, sceneNr = 1;
    while (y < pCorpus->getGT()->getAudioSampleCount()) {
    
			//chop audio-data scene-wise into smaller pieces
		  long int sceneStart, sceneEnd;
			pCorpus->getGT()->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward);
			if (sceneStart == -1 || sceneEnd == -1) {
				break; //scene not valid
			} else {
				printf("\n  Processing Scene %d...", sceneNr);
			}

			//load signal
			printf("\n    Phase I:\tLoading Signal...");
			SC_Signal* pSignal = pCorpus->loadSignal((unsigned long &)sceneStart, (unsigned long &)sceneEnd); //sceneStart and -end may get altered here...
			printf("done!");

			if (pSignal != NULL && pSignal->GetLen() > 0) {
				//extract features
				printf("\n    Phase II:\tExtracting Features...");
				SV_Data** pFeatures = pExtractor->extractFeatures(pCorpus, pSignal, sceneStart, sceneEnd, pTweak->modelHandler.speakerModelFeature);
				MFree_0D(pSignal);
				printf("done!");

				//concatenate multiple feature-sets (a normlization-matrix must be given, then)
				pExtractor->prepareFeatureSet(pFeatures, NULL); //bring features for standard parameters (i.e. for speaker clustering) "up" in the linked list
				pExtractor->equalizeFrameParameters(pFeatures, pExtractor->getFeatureCount(), pTweak->modelHandler.speakerModelFeature);
				SV_Data* pTmpFeature = pExtractor->combineFeatureVectors(pFeatures, pExtractor->getFeatureCount(), pTweak->modelHandler.speakerModelFeature);
				pExtractor->freeFeatures(pFeatures);
				if (pNorm != NULL)  {
					pExtractor->normalize(pTmpFeature, pNorm);
				}

				//build trajectories a priori so that the loadAndMerge mechanism can be used 
				pFeature = pExtractor->createTrajectories(pTmpFeature, pCorpus->getGT()->getConverter()->ms2audioFrame(this->pTweak->modelTime.syllableLength, pTmpFeature->Hdr.frameSize, pTmpFeature->Hdr.frameStep, sclib::alignmentStart), this->pTweak->modelTime.trajectoryStep, this->pTweak->modelTime.removeTiming, this->pTweak->modelTime.templateCount, this->pTweak->modelTime.clusteringIterations);
				MFree_0D(pTmpFeature);
	  
				//build linked list of all desired features of this corpus
				//pFeatureFirst is the pointer to the first element of the linked list, pFeatureHook to the current last element
				if (pFeature != NULL) {
					if (pFeatureFirst == NULL) {
						pFeatureFirst = pFeature;
						pFeatureHook = pFeatureFirst; 
					} else {
						pFeatureHook->Next = pFeature;
						pFeatureHook = pFeatureHook->Next;
					}
				}
			} //pSignal valid

			printf("\n");
			sceneNr++;
			y = sceneEnd + 1;
    }
    printf("\ndone!\n");

    printf("\nPhase 2:\tSaving Features...");
    pExtractor->saveFeatureList(tmpFeaturesFile, pFeatureFirst);
    sclib::destructLinkedList(pFeatureFirst);
    pFeatureFirst = NULL;
    pFeatureHook = NULL;
    printf("done!\n");

  } //!fileExists(featureFile)

  printf("\nPhase 3:\tBuilding models...");
  if (sclib::fileExists(tmpFeaturesFile) == true) {
    //save and load data to/from disk, to save one copy in memory...
		pFeature = pExtractor->loadAndMergeFeatures(tmpFeaturesFile); 
		//somehow, allocating the space for whole TIMIT train seems to crash if attempted directly after processing the data with the code in phase 2; however, after a restart of the programm (and then without phase 2, becasue the file already exists), it works

		//normalization-filename given but file doesn't exist => build it!
		if (strncmp(pTweak->modelTime.normalizationFile, "", sclib::bufferSize) != 0 && sclib::fileExists(pTweak->modelTime.normalizationFile) == false) {
			printf("\nNormalization matrix was needed but not existent => I now build it, then please restart the program and extract features again (I delete the old unnormalized feature file below!)");
			MFree_0D(pNorm); //should already be NULL normally...
			pNorm = pExtractor->createNormalizationMatrix(pFeature);
			if (pExtractor->saveFeature(pTweak->modelTime.normalizationFile, pNorm) != true) {
				printf("\nError while saving normalization matrix...");
				finalResult = false;
			}
			remove(tmpFeaturesFile);
			return true;
		}

		//finally, build the model
		pTweak->modelTime.checkForTrajectorization = false; //features have already been trajectorized above, and bevause somehow the relevant bit in the header is lost in one of the copying steps above, re-trajectorization would take place which makes no sense
		SC_Model* pModel = pModeller->buildModel(pFeature, NULL, 0, sclib::mtTime, 0);
		unsigned long int res = pModeller->saveModel(pTweak->modelTime.worldModelFile, pModel);
		if (res <= 0) {
			printf("\nSorry, the model could not be saved\n");
			finalResult = false;
		}
		MFree_0D(pModel);
    MFree_0D(pFeature);

	  printf("done!\n");
  } else {
    printf("no data!\n");
		finalResult = false;
  }

	printf("\nPhase 3:\tCleaning up...");

	MFree_0D(pNorm);
  MFree_0D(pCorpus);
  MFree_0D(pModeller);
  MFree_0D(pExtractor);

	printf("done!");
	printf("\n\"Good bye\".\n");
	return finalResult;
}
