/**************************************************************************/
/*    Responsibility:																											*/
/*      - Implements the audio-segmenter published in "Content-based Audio*/
/*        Classification and Segmentation by Using Support Vector         */
/*        Machines", Lu/Zhang/Li 2003                                     */
/*      - Provides possibilities to train & use the classifier            */
/*      - Capable to separate speech/non-pure speech/music/noise from     */
/*        audiofeatures of a non-silent signal                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.03.2006																								*/
/**************************************************************************/

#include "SC_Segmentation_AudioType_LZL.h"
#include "SC_Aux.h"
#include "SC_ClassifierTree.h"
#include "SC_FeatureHandler.h"
#include "SC_MatrixFunctions.h"

//experimental code
#include "SC_ModelHandler.h"
//end experimental

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_AudioType_LZL::SC_Segmentation_AudioType_LZL(SC_TweakableParameters* pTweak, bool cacheFeatures2Disk) : SC_Segmentation_AudioType(pTweak) {
  unsigned int x;

  this->cacheFeatures = cacheFeatures2Disk;
	this->cacheFilePrefix = sclib::exchangeFileExtension(this->pTweak->segmentationAudioTypeLzl.featureFileName, "");

  //deletion has do be done by hand by the user...
  /*
  //delete all old cache-files
  for (x = 0; x < 32; x++) {
    sprintf(fileName, "%s.%d\0", this->cacheFilePrefix, x);
    remove(fileName);
  }
  */

  for (x = 0; x < 32; x++) {
    this->pTrainingData[x] = NULL;
  }
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_AudioType_LZL::~SC_Segmentation_AudioType_LZL() {
  MFree_1D(this->cacheFilePrefix);

  for (unsigned int x = 0; x < 32; x++) {
    MFree_0D(this->pTrainingData[x]);
  }
}

//====================================================================================================================
// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
// implementation
// The features needed for training are meant to reside in an array of feature-containers as returned by 
// SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
//====================================================================================================================
int SC_Segmentation_AudioType_LZL::trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
  int res = SVLIB_Fail;
  
  //assure that all needed features and feature-types are provided
  if (pFeatures == NULL || checkFeatures(pFeatures) == true) {
    res = trainLzlAlgorithm(pGT, segmentStart, segmentEnd, pFeatures);
  }

  return res;
}

//====================================================================================================================
//	Check if the provided features are suitable for the purpose of this class
//====================================================================================================================
bool SC_Segmentation_AudioType_LZL::checkFeatures(SV_Data **pFeatures) {
  bool res = false;
  
  //assure that all needed features and feature-types are provided
  if (pFeatures[sclib::bitPosition(sclib::featureMFCC)] != NULL && pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Col >= 8) {
    if (pFeatures[sclib::bitPosition(sclib::featureZCR)] != NULL) {
      if (pFeatures[sclib::bitPosition(sclib::featureSubbandPower)] != NULL) {
        if (pFeatures[sclib::bitPosition(sclib::featureBrightnessBandwidth)] != NULL) {
          if (pFeatures[sclib::bitPosition(sclib::featureSpectrumFlux)] != NULL) {
            if (pFeatures[sclib::bitPosition(sclib::featureBandPeriodicity)] != NULL) {
              if (pFeatures[sclib::bitPosition(sclib::featureNFR)] != NULL) {
                res = true;
              } else {
                REPORT_ERROR(SVLIB_BadArg, "NFR feature missing, but necessary for audio-type training");
              }
            } else {
              REPORT_ERROR(SVLIB_BadArg, "Band-periodicity feature missing, but necessary for audio-type training");
            }
          } else {
            REPORT_ERROR(SVLIB_BadArg, "Spectrum-flux feature missing, but necessary for audio-type training");
          }
        } else {
          REPORT_ERROR(SVLIB_BadArg, "Brightness & bandwidth feature missing, but necessary for audio-type training");
        }
      } else {
        REPORT_ERROR(SVLIB_BadArg, "Subband-power feature missing, but necessary for audio-type training");
      }
    } else {
      REPORT_ERROR(SVLIB_BadArg, "ZCR feature missing, but necessary for audio-type training");
    }
  } else {
    REPORT_ERROR(SVLIB_BadArg, "MFCC feature missing, but necessary for audio-type training");
  }

  return res;
}

//====================================================================================================================
//	This is the real training algorithm as described in the paper "Content-based Audio Classification and Segmentation 
//  by Using Support Vector Machines", Lu/Zhang/Li 2003. 
//  The feature-vectors should be classified into silence/non-silence beforehand because this information is used 
//  inside
//====================================================================================================================
int SC_Segmentation_AudioType_LZL::trainLzlAlgorithm(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
  int res = SVLIB_Ok, *labels = NULL;
  int **treeTable, pos_0[2], pos_1[1], pos_2[1], pos_3[1], neg_0[3], neg_1[1], neg_2[2], neg_3[1], *posLabels[4], *negLabels[4], posLabelCount[4], negLabelCount[4], posResults[4], negResults[4];
  unsigned long int x, subClipLength = pGT->getConverter()->ms2sample(this->pTweak->segmentationAudioTypeLzl.subClipLength);
  long int neededTypes[5] = {sclib::atPureSpeech, sclib::atNoisySpeech, sclib::atBackground, sclib::atAction, sclib::atMusic}; //this MUST be single (not or-concatenated) types!
	double errors;
  SV_Data *pNorm, *pFeatureList, *pFeatureHook,*pFinalFeatures = NULL;
  SC_Segmentation_AudioType_LZL::SC_Labels *pLabelsList, *pLabelsHook;
  SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);
  SC_ClassifierTree *pTree = NULL;
	SC_TweakableParameters::SC_ClassifierSVMPar oldSvmParameters;

  //fill the traingData-members (maybe in addition to what has already been done)
  for (x = 0; x < 5; x++) {
    partiallyLoadTrainingData(pGT, segmentStart, segmentEnd, pFeatures, neededTypes[x]);
  }

  if (this->cacheFeatures == true) { //load from disk-cache
    res = (loadFromCache(neededTypes, 5, pFinalFeatures, labels) == true) ? SVLIB_Ok : SVLIB_Fail;
  } else {//construct labels and features from memory-cache
    for (x = 0; x < 5; x++) {
      //check traningData-members and construct linked list of features and labels, if ok
      if (this->pTrainingData[neededTypes[x]] != NULL) {
        if (x == 0) { //remember roots of linked lists
          pLabelsList= new SC_Segmentation_AudioType_LZL::SC_Labels(neededTypes[x], this->pTrainingData[neededTypes[x]]->Row);
          pLabelsHook = pLabelsList;
          pFeatureList = this->pTrainingData[neededTypes[x]];
          pFeatureHook = pFeatureList;
        } else { //add elemets to the tail of the linked lists
          pLabelsHook = pLabelsHook->addList(neededTypes[x], this->pTrainingData[neededTypes[x]]->Row);
          pFeatureHook->Next = this->pTrainingData[neededTypes[x]];
          pFeatureHook = pFeatureHook->Next;
        }
      } else {
        REPORT_ERROR(SVLIB_BadArg, "Needed Training-data for audio-type classification was partially not provided");
        res = SVLIB_Fail;
        break;
      }
    }

    //free partially loaded traing data
    for (x = 0; x < 32; x++) {
      MFree_0D(this->pTrainingData[x]);
    }

    //construct final features and labels
    pFinalFeatures = pFeatureList->MergeData();
    sclib::destructLinkedList(pFeatureList);
    labels = pLabelsList->merge(x);
    assert(x == pFinalFeatures->Row); //each feature-vector should have one label
    sclib::destructLinkedList(pLabelsList);
  }

  if (res != SVLIB_Fail && pFinalFeatures != NULL && labels != NULL) { //don't do anything if feature-types where missing
    //create normalization-matrix
		pNorm = pFeatureHandler->createStandardizationMatrix(pFinalFeatures);

    //normalize features
    pFeatureHandler->standardize(pFinalFeatures, pNorm);

    //save normalization matrix for use in testing
    pFeatureHandler->saveFeature(this->pTweak->segmentationAudioTypeLzl.normalizationFileName, pNorm);

    //construction of binary classification tree with node-numbers and semantic meaning:
    //
    //                    0 (speech/non-speech)
    //             +      |      -
    //     -------------------------------
    //    /                               \
    //   1 (pure-speech/noisy-speech)      2 (music/background-sound)
    //                                +    |    -
    //                           -------------------- 
    //                          /                    \
    //                       music                    3 (background/action)
    //
    MArray_2D(treeTable, 4, 2, int, "SC_Segmentation_AudioTYpe_LZL.trainLzlAlgorithm: treeTable");
    treeTable[0][0] = 1;  treeTable[0][1] = 2; //left and right (pos/neg) sibling for 0th node
    treeTable[1][0] = -1; treeTable[1][1] = -1;
    treeTable[2][0] = -1; treeTable[2][1] = 3;
    treeTable[3][0] = -1; treeTable[3][1] = -1;

    //create label-lists for the nodes: which node has to be trainied with which labels?
    posLabelCount[0] = 2; posLabels[0] = pos_0;
    posLabelCount[1] = 1; posLabels[1] = pos_1;
    posLabelCount[2] = 1; posLabels[2] = pos_2;
    posLabelCount[3] = 1; posLabels[3] = pos_3;
    negLabelCount[0] = 3; negLabels[0] = neg_0;
    negLabelCount[1] = 1; negLabels[1] = neg_1;
    negLabelCount[2] = 2; negLabels[2] = neg_2;
    negLabelCount[3] = 1; negLabels[3] = neg_3;
    pos_0[0] = sclib::atPureSpeech; pos_0[1] = sclib::atNoisySpeech;
    neg_0[0] = sclib::atBackground; neg_0[1] = sclib::atMusic; neg_0[2] = sclib::atAction;
    pos_1[0] = sclib::atPureSpeech;
    neg_1[0] = sclib::atNoisySpeech;
    pos_2[0] = sclib::atMusic;
    neg_2[0] = sclib::atBackground; neg_2[1] = sclib::atAction;
    pos_3[0] = sclib::atBackground;
    neg_3[0] = sclib::atAction;

    //result to be returned from each node (only relevant for the leaves, though)
    posResults[0] = sclib::atSpeech;                      negResults[0] = sclib::atNoise;
    posResults[1] = sclib::atSpeech|sclib::atPureSpeech;  negResults[1] = sclib::atSpeech|sclib::atNoisySpeech;
    posResults[2] = sclib::atNoise|sclib::atMusic;        negResults[2] = sclib::atNoise;
    posResults[3] = sclib::atNoise|sclib::atBackground;   negResults[3] = sclib::atNoise|sclib::atAction;

    //create the classifier and free temporary data (the tree construction data has been copied inside the tree)
		oldSvmParameters = this->pTweak->classifierSvm; //save general parameters...
		this->pTweak->classifierSvm = this->pTweak->segmentationAudioTypeLzl.svmParameters; //...and use special ones
    pTree = new SC_ClassifierTree(this->pTweak, false, true, false, 4, treeTable, posLabels, posLabelCount, negLabels, negLabelCount, posResults, negResults, sclib::ctSVM);

    //train & save classifier
    //TODO: doScaling?!?
    errors = pTree->trainMultiClass(pFinalFeatures, labels);
    if (errors != (double)(SVLIB_Fail)) {
      res = pTree->saveClassifier( this->pTweak->segmentationAudioTypeLzl.classifierFileName);
		} else {
			res = SVLIB_Fail;
		}

		this->pTweak->classifierSvm = oldSvmParameters; //restore original parameters
  }

  MFree_0D(pTree);
  MFree_0D(pFinalFeatures);
  MFree_1D(labels);
  MFree_0D(pFeatureHandler);
  MFree_2D(treeTable);

  return res;
}

//====================================================================================================================
//	Loads all features from the disk-cache in ready-to-use form (no more merging needed, all done on disk) togehter
//  with the respective labels; returns true on success, false otherwise
//====================================================================================================================
bool SC_Segmentation_AudioType_LZL::loadFromCache(long int *neededTypes, int typeCount, SV_Data* &pFeatures, int* &labels) {
  bool res;
  int z, cols;
  char fileName[sclib::bufferSize];
  long int rows, x, y, overallRows = 0;
  SV_DataIO IO;
  SV_Data *pTemp;

  //prepare variables
  MFree_0D(pFeatures);
  MFree_1D(labels);
  
  //measure overall size
  res = true;
  for (z = 0; z < typeCount; z++) {
    sprintf(fileName, "%s.%d\0", this->cacheFilePrefix, sclib::bitPosition(neededTypes[z])); //construct the final fileName with log(audioType) as the extension
    if (sclib::fileExists(fileName) == true) {
      IO.OpenFile(fileName, READ_REC);
      rows = IO.getAllRecRowCount(cols);
      if (cols > 0 && rows > 0) {
        overallRows += rows;
      } else {
        res = false;
      }
      IO.CloseFile();
    } else {
      res = false;
    }
  }

  //really load data now
  if (res != false) {
    pFeatures = new SV_Data(overallRows, cols);
    labels = new int[overallRows];
    rows = 0;

    for (z = 0; z < typeCount; z++) {
      sprintf(fileName, "%s.%d\0", this->cacheFilePrefix, sclib::bitPosition(neededTypes[z]));
      IO.OpenFile(fileName, READ_REC);
      pTemp = IO.GetDataRec();
      while (pTemp != NULL) {
        for (y = 0; y < pTemp->Row; y++) {
          for (x = 0; x < pTemp->Col; x++) {
            pFeatures->Mat[rows + y][x] = pTemp->Mat[y][x];
          }
          labels[rows + y] = neededTypes[z];
        }
        rows += pTemp->Row;
        MFree_0D(pTemp);
        pTemp = IO.GetDataRec();
      }
      IO.CloseFile();
    }
  }
  
  return res;
}

//====================================================================================================================
//  If the training-database is too big to be hold in memory completely (not forgetting the number of copies for 
//  copying-together, aggregation, normalization, and inside the classifier), this function can be called several 
//  times prior to call the trainClassifier(): Each time, a different subpart, representing one needed audio-type,
//  can be provided (even the same audioType can be provided in parts); it's feature-vectors get aggregated per 
//  subclip and afterwards stored for training, so the original features can be destroyed. After providing 
//  feature-sets for all relevant audio-types, the training-algorithm can be called.
//  pFeatures is meant to be an array as returned by SC_FeatureHandler.extractFeatures(), and must contain frames 
//  of the specified (single, not or-concatenated) audio-type (but can also contain others, though this gives no 
//  sense...)
//====================================================================================================================
int SC_Segmentation_AudioType_LZL::partiallyLoadTrainingData(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned long int audioType) {
  int res = SVLIB_Fail, index = sclib::bitPosition(audioType), x;
  int idx[7] = {sclib::bitPosition(sclib::featureMFCC), sclib::bitPosition(sclib::featureSubbandPower), sclib::bitPosition(sclib::featureZCR), sclib::bitPosition(sclib::featureBrightnessBandwidth), sclib::bitPosition(sclib::featureSpectrumFlux), sclib::bitPosition(sclib::featureBandPeriodicity), sclib::bitPosition(sclib::featureNFR)};
  int startCol[7] = {0, 0, 0, 0, 0, 0, 0};
  int endCol[7] = {7, 0, 0, 0, 0, 1, 0};
  bool addSd[7] = {true, true, true, true, false, false, false};
  bool divByLengthMinusOne[7] = {false, false, false, false, true, false, false};
  unsigned long int subClipLength = pGT->getConverter()->ms2sample(this->pTweak->segmentationAudioTypeLzl.subClipLength);
  char buffer[sclib::bufferSize];
	unsigned long int row, column;
  SV_Data *pTempAggregated[7], *pTempWOS, *pTempCombined;
  SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);

  if (pFeatures != NULL && checkFeatures(pFeatures) == true) {
    if (pGT->existsSegmentType(segmentStart, segmentEnd, audioType) == true) {
      //only allow "single" audioTypes, no combinations here!
      if (sclib::isPowerOfTwo(audioType) == true) {
        //aggregate each feature per subclip in the appropriate manner and store it in the final matrix
        //omit all feature-vectors previously classified as silence from this step
        for (x = 0; x < 7; x++) {
          pTempWOS = pGT->copyFramesTogether(pFeatures[idx[x]], segmentStart, segmentStart, segmentEnd, audioType, sclib::atSilence);
          pTempAggregated[x] = aggregatePerSubClip(pTempWOS, subClipLength, startCol[x], endCol[x], addSd[x], divByLengthMinusOne[x], false);
          MFree_0D(pTempWOS);
        }
        pTempCombined = pFeatureHandler->combineFeatureVectors(pTempAggregated, 7);
        if (!pFeatureHandler->checkFeatures(pTempCombined, row, column)) {
          sprintf(buffer, "Corrupted features after aggregation in row %d and column %d!", row, column);
          REPORT_ERROR(SVLIB_BadData, buffer);
        }

        for (x = 0; x < 7; x++) {
          MFree_0D(pTempAggregated[x]);
        }

        if (cacheFeatures == true) { //save features to disk if caching was wished during construction of this object; this saves lots of memory!
          sprintf(buffer, "%s.%d\0", this->cacheFilePrefix, index);
          pFeatureHandler->saveFeature(buffer, pTempCombined); //construct the final fileName with log(audioType) as the extension
          MFree_0D(pTempCombined);
        } else {
          if (this->pTrainingData[index] == NULL) { //this is the first feature-vector-portion of this type => just store it
            this->pTrainingData[index] = pTempCombined;
          } else { //there are already some feature-vector of this type => combine them with the new ones
            this->pTrainingData[index]->Next = pTempCombined;
            pTempWOS = this->pTrainingData[index]->MergeData();
            MFree_0D(this->pTrainingData[index]);
            MFree_0D(pTempCombined);
            this->pTrainingData[index] = pTempWOS;
          }
        }

        res = SVLIB_Ok;
      }
    } else {
      if (this->pTrainingData[index] == NULL && this->cacheFeatures == false) { //only complain if the respective training-data is still empty, but shouldn't be
        REPORT_ERROR(SVLIB_BadArg, "Given feature-type missing, but necessary for audio-type training");
      }
    }
  }

  MFree_0D(pFeatureHandler);

  return res;
}

//====================================================================================================================
// classifiy the given audio-segment according to the underlying audio-type (speech/noise/music/...)
// pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
// the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
//====================================================================================================================
int SC_Segmentation_AudioType_LZL::classifyAudioType(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
  int res = SVLIB_Fail;
  
  //assure that all needed features are provided
  if (checkFeatures(pFeatures) == true) {
    res = lzlAlgorithm(pGT, segmentStart, segmentEnd, pFeatures);

    //this threshold (in ms) defines how silence-regions within speech should be handled: shorter silences will be regarded 
    //as pauses which means that they are not taken to estimate speaker-models, but they don't fragment a 
    //speech-segment in 2; longer segments are not considered and divide the speech-segment into 2 distinct ones 
    //(which could be spoken by different speakers)
	  pGT->silence2pause(segmentStart, segmentEnd, pGT->getConverter()->ms2sample(pTweak->general.pauseSilenceThreshold));
  }

  return res;
}

//====================================================================================================================
//	This is the realy classification/segmentation algorithm as described in the paper "Content-based Audio 
//  Classification and Segmentation by Using Support Vector Machines", Lu/Zhang/Li 2003. Though most if it's work is
//  done in the feature-extraction/feature-handler/classifier-classes, it is rather short & simple
//  Because the classifier is trained on signals of a specific sampleRate and the testdata needs to look alike, this
//  Algorithm possibly receives features which have a different SR as stated in the groundTruth (but match the 
//  classifier). This is compensated for here (and only here throughout the lib except feature extraction, at the 
//  moment).
//====================================================================================================================
int SC_Segmentation_AudioType_LZL::lzlAlgorithm(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
  int *results = NULL, res = SVLIB_Fail;
  unsigned long int subClipLength, start, end, x, y; //, z;
  double srRatio = (double)(pGT->getAudioSampleRate()) / (double)(pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.sampleRate);
	double pos_0, pos_1, pos_2, pos_3, neg_0, neg_1, neg_2, neg_3;
  SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);
  SC_Classifier *pClassifier;
  SV_Data *pNorm = NULL, *pFinalFeatures = NULL, **pTemp, *pProbabilities = NULL;
	SC_TweakableParameters::SC_ClassifierSVMPar oldSvmParameters;

	oldSvmParameters = this->pTweak->classifierSvm; //save general parameters...
	this->pTweak->classifierSvm = this->pTweak->segmentationAudioTypeLzl.svmParameters; //...and use special ones
	pClassifier = new SC_ClassifierTree(this->pTweak, false, true, true);

	//experimental code following...
	/*
	unsigned long int startFrame, endFrame, frameSize, frameStep, idx = sclib::bitPosition(sclib::featureMFCC);
	double actionScore;
	SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak, false);
	SC_Model *pActionModel = pModelHandler->loadModel(this->pTweak->segmentationAudioTypeLzl.actionModelFileName, sclib::mtGMM_new);
	SV_Data *pFrame = new SV_Data(1, pFeatures[idx]->Col);
	*/
	//end experimental

	//another experimental snippet...
	/*
	MFree_0D(pClassifier);
	MFree_0D(pFeatureHandler);
	return actionThreshClassifyer(pGT, pFeatures[sclib::bitPosition(sclib::featureEnergyZCR)], segmentStart, segmentEnd);
	*/
	//end snippet

  //try to load the specified classifier (don't call sclib::fileExists() because maybe the specific loading operation first has to process the given fileName in a way before it really fits a physical file...)
  if (this->pTweak->segmentationAudioTypeLzl.classifierFileName != NULL && 
      strcmp("", this->pTweak->segmentationAudioTypeLzl.classifierFileName) != 0 && 
      pClassifier->loadClassifier(this->pTweak->segmentationAudioTypeLzl.classifierFileName) != SVLIB_Fail) {
    
    //try to load normalization parameters
    if (this->pTweak->segmentationAudioTypeLzl.normalizationFileName != NULL && 
        strcmp("", this->pTweak->segmentationAudioTypeLzl.normalizationFileName) != 0 && 
        sclib::fileExists(this->pTweak->segmentationAudioTypeLzl.normalizationFileName) == true && 
        (pNorm = pFeatureHandler->loadFeature(this->pTweak->segmentationAudioTypeLzl.normalizationFileName)) != NULL) {

      //delete old (groundTruth-inherited audio-type labels in the hypothesized part of GT's frameList
			pGT->setSegment(segmentStart, segmentEnd, SC_GroundTruth::getDetailedAudioTypes(), false, sclib::noSpeaker, sclib::modeLabelRemove);
      
      //aggregate each features per subclip in the appropriate manner and store it in the final matrix
      subClipLength = pGT->getConverter()->ms2sample(this->pTweak->segmentationAudioTypeLzl.subClipLength, sclib::alignmentStart, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.sampleRate);
      MArray_1D(pTemp, 7, SV_Data*, "SC_Segmentation_AudioType_LZL.trainLzlAlgorithm: pTemp");
      pTemp[0] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureMFCC)], subClipLength, 0, 7, true, false, true); //MFCC, only 1-8
      pTemp[1] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureSubbandPower)], subClipLength, 0, 0, true, false, true); //SubBandPower
      pTemp[2] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureZCR)], subClipLength, 0, 0, true, false, true); //ZCR
      pTemp[3] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureBrightnessBandwidth)], subClipLength, 0, 0, true, false, true); //Brightness&Bandwidth
      pTemp[4] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureSpectrumFlux)], subClipLength, 0, 0, false, true, true); //SpectrumFlux
      pTemp[5] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureBandPeriodicity)], subClipLength, 0, 1, false, false, true); //BandPeriodicity
      pTemp[6] = aggregatePerSubClip(pFeatures[sclib::bitPosition(sclib::featureNFR)], subClipLength, 0, 0, false, false, true); //NFR
      pFinalFeatures = pFeatureHandler->combineFeatureVectors(pTemp, 7);
      for (x = 0; x < 7; x++) {
        MFree_0D(pTemp[x]);
      }
      MFree_1D(pTemp);

      //normalize final matrix
      pFeatureHandler->standardize(pFinalFeatures, pNorm);

      //classify final feature-vectors
      results = pClassifier->classify(pFinalFeatures, pProbabilities);
      
      //sclib::vectorOut("results.txt", results, pFinalFeatures->Row, true, this->pTweak);
      //sclib::classOut("probs.txt", pProbabilities, this->pTweak);

      //save results in the frameList of the ground-truth
      for (y = 0; y < (unsigned long int)(pFinalFeatures->Row); y++) {
        //do smoothing
        if (y > 0 && y < (unsigned long int)(pFinalFeatures->Row-1)) {
		      if (results[y] != results[y-1] && results[y] != sclib::atSilence && sclib::bitTest(results[y], sclib::atAction) == false && results[y] != results[y+1]) { //new: no smoothing for ACTION
  		      results[y] = results[y-1];
						for (x = 0; x < (unsigned long int)(pProbabilities->Col); x++) { //smooth also probabilities
							pProbabilities->Mat[y][x] = pProbabilities->Mat[y-1][x];
						}
		      }
        }

        //segmentStart is based on global sampleRate as stored in the GT
        //the feature's frameSize is based on a possibly different sampleRate
        //this must be compensated for when writing to the sample-based (global SR dependant) frameList
        start = segmentStart + (y * (long)sclib::round(pFinalFeatures->Hdr.frameStep * srRatio));
        end = sclib::min((start + (long)sclib::round(pFinalFeatures->Hdr.frameSize * srRatio) - 1), segmentEnd); //because the last subClip might be underpopulated...

				/*
				//experimental code
				//special treatment for BACKGROUND subclips to possibly find ACTION within
				//if (sclib::bitTest(results[y], sclib::atBackground) == true) {
					frameSize = (unsigned long)(pFeatures[idx]->Hdr.frameSize * srRatio);
					frameStep = (unsigned long)(pFeatures[idx]->Hdr.frameStep * srRatio);
					startFrame = pGT->sample2audioFrame(start-segmentStart, frameSize, frameStep);
					endFrame = pGT->sample2audioFrame(end-segmentStart, frameSize, frameStep);
					for (z = startFrame; z <= sclib::min(endFrame, pFeatures[idx]->Row-1); z++) {
						for (x = 0; x < (unsigned long)(pFeatures[idx]->Col); x++) {
							pFrame->Mat[0][x] = pFeatures[idx]->Mat[z][x];
						}
						actionScore = pModelHandler->testModel(pActionModel, pFrame, 0);
						pGT->setProbability(segmentStart + pGT->audioFrame2sample(z, frameSize, frameStep, sclib::alignmentStart), segmentStart + pGT->audioFrame2sample(z, frameSize, frameStep, sclib::alignmentEnd), sclib::atUndefined, actionScore);
					}
					//pSubClip = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(sclib::featureMFCC)], segmentStart, start, end);
					//actionScore = pModelHandler->testModel(pActionModel, pSubClip, 0);
					//MFree_0D(pSubClip);
				//}
				//end experimental
				*/



				//get raw probability-information
				//ATTENTION: this uses internal knowledge of this class about how the tree-classifier is constructed to access the correct cols in pProbabilities as well as that the probabilities are real probabilities in [0..1]!
				//           for rememberance (tree as constructed by the training-algorithm above): 
				//
				//                    0 (speech/non-speech)
				//             +      |      -
				//     -------------------------------
				//    /                               \
				//   1 (pure-speech/noisy-speech)      2 (music/background-sound)
				//                                +    |    -
				//                           -------------------- 
				//                          /                    \
				//                       music                    3 (background/action)
				//
				//first, compute the probabilities of turning left and right (pos/neg) at each branch
				pos_0 = (sclib::round(pProbabilities->Mat[y][0]) == sclib::labelPositive) ? pProbabilities->Mat[y][1] : 1.0 - pProbabilities->Mat[y][1];
				neg_0 = (sclib::round(pProbabilities->Mat[y][0]) == sclib::labelNegative) ? pProbabilities->Mat[y][1] : 1.0 - pProbabilities->Mat[y][1];
				pos_1 = (sclib::round(pProbabilities->Mat[y][2]) == sclib::labelPositive) ? pProbabilities->Mat[y][3] : 1.0 - pProbabilities->Mat[y][3];
				neg_1 = (sclib::round(pProbabilities->Mat[y][2]) == sclib::labelNegative) ? pProbabilities->Mat[y][3] : 1.0 - pProbabilities->Mat[y][3];
				pos_2 = (sclib::round(pProbabilities->Mat[y][4]) == sclib::labelPositive) ? pProbabilities->Mat[y][5] : 1.0 - pProbabilities->Mat[y][5];
				neg_2 = (sclib::round(pProbabilities->Mat[y][4]) == sclib::labelNegative) ? pProbabilities->Mat[y][5] : 1.0 - pProbabilities->Mat[y][5];
				pos_3 = (sclib::round(pProbabilities->Mat[y][6]) == sclib::labelPositive) ? pProbabilities->Mat[y][7] : 1.0 - pProbabilities->Mat[y][7];
				neg_3 = (sclib::round(pProbabilities->Mat[y][6]) == sclib::labelNegative) ? pProbabilities->Mat[y][7] : 1.0 - pProbabilities->Mat[y][7];

				//do specificity corrections according to user's wishes; again, this uses internal knowledge!
				if (this->pTweak->segmentationAudioTypeLzl.speechSpecificity!=0.5 && this->pTweak->segmentationAudioTypeLzl.musicSpecificity!=0.5) {
					pos_0 = sclib::getBetween(0.0, sclib::invSigmoid(this->pTweak->segmentationAudioTypeLzl.speechSpecificity, 0.4, 1.0, 2.5, 0.000001, 0.0), 1.0); //TODO: use sigmoidal instead of linar scaling? (p=0.3)
					neg_0 = 1.0 - pos_0;
					pos_2 = sclib::getBetween(0.0, sclib::invSigmoid(this->pTweak->segmentationAudioTypeLzl.musicSpecificity, 0.4, 1.0, 2.5, 0.000001, 0.0), 1.0);
					neg_2 = 1.0 - pos_2;
					if (pos_0 > neg_0) { //refine classification based on specificity computations
						if (pos_1 > pos_2) {
							results[y] = sclib::atSpeech|sclib::atPureSpeech;
						} else {
							results[y] = sclib::atSpeech|sclib::atNoisySpeech;
						}
					} else {
						if (pos_2 > neg_2) {
							results[y] = sclib::atNoise|sclib::atMusic;
						} else {
							if (pos_3 > neg_3) {
								results[y] = sclib::atNoise;
							} else {
								results[y] = sclib::atNoise|sclib::atAction;
							}
						}
					}
				}

        //write classification result to the frameList
        pGT->setSegmentIf(start, end, sclib::noType, false, sclib::atSilence, true, results[y], false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //don't overwrite the (maybe previously written) SILENCE tag
				//TODO: this helps to improve the detection of speech in music for Meine_Eltern.mpg
				/*if (pos_0 > 0.05) {
	        pGT->setSegmentIf(start, end, sclib::noType, false, sclib::atSilence, true, sclib::atSpeech|sclib::atPureSpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //don't overwrite the (maybe previously written) SILENCE tag
				}*/

				//then, compute the probability of each audio-type by multiplicating all probabilities along the path from the tree's root to the type itself
				pGT->setProbability(start, end, sclib::atSpeech, pos_0);
				pGT->setProbability(start, end, sclib::atNoise, neg_0);
				pGT->setProbability(start, end, sclib::atPureSpeech, pos_0*pos_1);
				pGT->setProbability(start, end, sclib::atNoisySpeech, pos_0*neg_1);
				pGT->setProbability(start, end, sclib::atMusic, neg_0*pos_2);
				pGT->setProbability(start, end, sclib::atBackground, neg_0*neg_2*neg_3);
				pGT->setProbability(start, end, sclib::atAction, neg_0*neg_2*pos_3);
      }

  		//give all passages that have not been processed (if any) the UNDEF label
      pGT->setSegmentIf(segmentStart, segmentEnd, sclib::noType, false, sclib::atSilence|sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atMusic|sclib::atBackground|sclib::atAction|sclib::atBreath, false, sclib::atUndefined);

      //clean up
      MFree_0D(pFinalFeatures);
      MFree_0D(pNorm);
      MFree_1D(results);
      MFree_0D(pProbabilities);

      res = SVLIB_Ok;
    } else {
      REPORT_ERROR(SVLIB_FileErr, "Specified normalization file could not be loaded from the file-system");
    }
  } else {
    REPORT_ERROR(SVLIB_FileErr, "Specified classifier could not be loaded from the file-system");
  }

	//experimental code following...
	/*
	MFree_0D(pModelHandler);
	MFree_0D(pActionModel);
	MFree_0D(pFrame);
	*/
	//end experimental

	this->pTweak->classifierSvm = oldSvmParameters; //restore original parameters
  MFree_0D(pClassifier);
  MFree_0D(pFeatureHandler);

  return res;
}

//TODO: testtesttest
/*
int SC_Segmentation_AudioType_LZL::actionThreshClassifyer(SC_GroundTruth *pGT, SV_Data *pEnergy, unsigned long int segmentStart, unsigned long int segmentEnd) {
	int clipSize = 250; //ms
	int clipCount = (int)ceil((clipSize - (pGT->sample2ms(pEnergy->Hdr.frameSize-pEnergy->Hdr.frameStep, pEnergy->Hdr.sampleRate))) / (double)(pGT->sample2ms(pEnergy->Hdr.frameStep, pEnergy->Hdr.sampleRate))); //number of frames needed to fill a clip of wanted size
	int start = 0, end;
	double *mean, srRatio = (double)(pGT->getAudioSampleRate()) / (double)(pEnergy->Hdr.sampleRate);
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	unsigned long int frameSize, frameStep;
	double thresh = 63.0;
	int action = 0, no_action = 0;
	double ratio;
	
	while (start < pEnergy->Row) {
		end = sclib::min(pEnergy->Row-1, start + clipCount);
		mean = pFunc->getMeanV(pEnergy->Mat, pEnergy->Row, pEnergy->Col, start, end);
		
		//sclib::scalarOut("energy4.txt", mean[0], this->pTweak);

		frameSize = (unsigned long)(pEnergy->Hdr.frameSize * srRatio);
		frameStep = (unsigned long)(pEnergy->Hdr.frameStep * srRatio);
		pGT->setProbability(segmentStart + pGT->audioFrame2sample(start, frameSize, frameStep, sclib::alignmentStart), segmentStart + pGT->audioFrame2sample(end, frameSize, frameStep, sclib::alignmentEnd), sclib::atAction, mean[0]-thresh);

		if (mean[0] >= thresh) {
			action++;
			pGT->setSegment(segmentStart + pGT->audioFrame2sample(start, frameSize, frameStep, sclib::alignmentStart), segmentStart + pGT->audioFrame2sample(end, frameSize, frameStep, sclib::alignmentEnd), sclib::atAction);
		} else {
			no_action++;
		}

		MFree_1D(mean);
		start = end + 1;
	}
	
	ratio = action / (double)(no_action);
	MFree_0D(pFunc);

	return SVLIB_Ok;
}
*/

//====================================================================================================================
//	Convert the per-frame-features to per-subclip-features and return a new SV_Data object containing:
//   - The following values for the columns starting with startCol and ending with endCol (or the last available col,
//     if endCol=0)
//   - The mean (if divByLengthMinusOne=false; otherwise mean*(T/(T-1)) with T number of frames in a subclip)
//   - The standard deviation of the features in the subclip, if addSd=true
//  Subclip-length has to be given in samples
//  At the end there may be some frames missing that didn't fully fit into a subclip; we don't care about them if 
//  createUnderpopulatedFrames==false
//====================================================================================================================
SV_Data* SC_Segmentation_AudioType_LZL::aggregatePerSubClip(SV_Data *pFeatures, unsigned long int subClipLength, unsigned int startCol, unsigned int endCol, bool addSd, bool divByLengthMinusOne, bool createUnderpopulatedFrames) {
  long int firstFrame, lastFrame, subClip, subClipCount, subClipLengthInFrames, d, D = (endCol == 0 || endCol > (unsigned int)pFeatures->Col) ? pFeatures->Col : (long int)endCol+1;
  double *mean = NULL, *sd = NULL;
  SC_MatrixFunctions *pMFunc = new SC_MatrixFunctions();
  SV_Data *pResult = NULL;

  //this is the length of a subClip in frames; the result has min. subClipLength samples, max. subClipLength+frameLength-1 samples
  subClipLengthInFrames = ((subClipLength - pFeatures->Hdr.frameSize) / pFeatures->Hdr.frameStep) + 1;

  //don't allow underpolpulated frames normally, they will deteriorate the statistics and in really odd cases (only 1 frame) give -1.#IND results in case of sd (division by length-1=0!)
	if (createUnderpopulatedFrames == false) {
		subClipCount = sclib::getRowCount(pFeatures->Row, subClipLengthInFrames, subClipLengthInFrames); //for training
	} else {
		subClipCount = pFeatures->Row / subClipLengthInFrames + (pFeatures->Row % subClipLengthInFrames > 0 ? 1 : 0); //for testing
	}

  if (subClipCount > 0) {
    pResult = new SV_Data(subClipCount, D * ((addSd == true) ? 2 : 1));
    pResult->Hdr.frameSize = subClipLength; //((subClipLengthInFrames - 1) * pFeatures->Hdr.frameStep) + pFeatures->Hdr.frameSize;
    pResult->Hdr.frameStep = pResult->Hdr.frameSize;
    pResult->Hdr.sampleRate = pFeatures->Hdr.sampleRate;

    for (subClip = 0; subClip < subClipCount; subClip++) {
      firstFrame = subClip * subClipLengthInFrames;
      lastFrame = sclib::min(firstFrame+subClipLengthInFrames-1, pFeatures->Row-1);

      mean = pMFunc->mean(pFeatures->Mat, pFeatures->Row, pFeatures->Col, firstFrame, lastFrame, startCol, D);
      for (d = 0; d < D; d++) {
        pResult->Mat[subClip][d] = (float)(mean[d]);
        if (divByLengthMinusOne == true) {
          pResult->Mat[subClip][d] *= (float)(subClipLengthInFrames) / (float)(subClipLengthInFrames-1);
        }
      }
      
      if (addSd == true) {
        sd = pMFunc->std(pFeatures->Mat, pFeatures->Row, pFeatures->Col, mean, firstFrame, lastFrame, startCol, D);
        for (d = 0; d < D; d++) {
					pResult->Mat[subClip][D+d] = (sclib::isFinite(sd[d]) == true) ? (float)(sd[d]) : 0.0f; //in case of underpopulated subclips with length 1, a div-by-0 error might have occured that we correct here...
        }
        MFree_1D(sd);
      }

      MFree_1D(mean);
    }
  }

  MFree_0D(pMFunc);

  return pResult;
}

//====================================================================================================================
// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
// are handled therein)
//====================================================================================================================
unsigned long int SC_Segmentation_AudioType_LZL::getUncertaintyRegionWidth(void) {
	unsigned long int width = 0;
	int max = 0;

	//factors due to frame-based analysis
	if (max < this->pTweak->featureBandPeriodicity.frameStep) {
		max = this->pTweak->featureBandPeriodicity.frameStep;
	}
	if (max < this->pTweak->featureBrightnessBandwidth.frameStep) {
		max = this->pTweak->featureBrightnessBandwidth.frameStep;
	}
	if (max < this->pTweak->featureMfcc.frameStep) {
		max = this->pTweak->featureMfcc.frameStep;
	}
	if (max < this->pTweak->featureNfr.frameStep) {
		max = this->pTweak->featureNfr.frameStep;
	}
	if (max < this->pTweak->featureSpectrumFlux.frameStep) {
		max = this->pTweak->featureSpectrumFlux.frameStep;
	}
	if (max < this->pTweak->featureSubBandPower.frameStep) {
		max = this->pTweak->featureSubBandPower.frameStep;
	}
	if (max < this->pTweak->featureZcr.frameStep) {
		max = this->pTweak->featureZcr.frameStep;
	}
	width += (unsigned long int)(2 * max);

	//factors due to windowed (sub-clip) analysis
	width += 2 * this->pTweak->segmentationAudioTypeLzl.subClipLength;

	return width;
}

//====================================================================================================================
// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
// parameters
//====================================================================================================================
SC_TweakableParameters::SC_FeaturePar* SC_Segmentation_AudioType_LZL::getSpecialFeatureParameters(void) {
	//properly link the needed parameter-sets
	this->pTweak->segmentationAudioTypeLzl.bbParameters.Next = &(this->pTweak->segmentationAudioTypeLzl.bpParameters);
	this->pTweak->segmentationAudioTypeLzl.bpParameters.Next = &(this->pTweak->segmentationAudioTypeLzl.mfccParameters);
	this->pTweak->segmentationAudioTypeLzl.mfccParameters.Next = &(this->pTweak->segmentationAudioTypeLzl.sbpParameters);
	this->pTweak->segmentationAudioTypeLzl.sbpParameters.Next = &(this->pTweak->segmentationAudioTypeLzl.sfParameters);
	this->pTweak->segmentationAudioTypeLzl.sfParameters.Next = &(this->pTweak->segmentationAudioTypeLzl.zcrParameters);
	this->pTweak->segmentationAudioTypeLzl.zcrParameters.Next = &(this->pTweak->segmentationAudioTypeLzl.nfrParameters);
	this->pTweak->segmentationAudioTypeLzl.nfrParameters.Next = NULL;

	return &(this->pTweak->segmentationAudioTypeLzl.bbParameters);
}
