/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle Wei-Ho (Wesley) Tsais test    */
/*        data                                                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 31.08.2005																								*/
/**************************************************************************/

#include <limits.h>
#include "SC_GroundTruth.h"
#include "SC_FeatureHandler.h"
#include "SC_Aux.h"
#include "SC_MixtureModel_GMM.h"
#include "SC_MixtureModel_MIXMAX.h"
#include "SC_MixtureModel_MIX2MAX.h"
#include "SC_MixtureModel_MIX2MAXex.h"
#include "SC_Corpus_Wesley.h"
#include "SC_Transform.h"
#include "SC_Signal_WAVE.h"
#include "SC_GroundTruth_Wesley.h"
#include <SV_General.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Corpus_Wesley::SC_Corpus_Wesley(SC_TweakableParameters* pTweak, unsigned long int sampleRate, bool mfcc2lfbe, unsigned long int modelOrigin) : SC_Corpus(pTweak) {
  this->pModeller = new SC_ModelHandler(this->pTweak);
  this->convertFeatures = mfcc2lfbe;
  this->modelOrigin = modelOrigin;
  this->pGT = new SC_GroundTruth_Wesley(this->pTweak, sampleRate);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Corpus_Wesley::~SC_Corpus_Wesley() {
  MFree_0D(this->pModeller);
}

//====================================================================================================================
//	load or train models, return linked list of clusters
//  if the postfix is given, the trained models get saved immediately as if saved by savedModel()
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::constructModels(const char* targetListFileName, const char* corpusFileName, unsigned long int layer, const char* savePostFix) {
  unsigned long int count = 0;
  FILE *listFile;
  SC_Model *pBackgroundModel = NULL;
  SC_Cluster *pFirstCluster = NULL, *pClusterHook = NULL, *pActualCluster = NULL;
  char ID[30], modelFileName[sclib::bufferSize], bgModelFileName[sclib::bufferSize];

  if ((listFile = fopen(targetListFileName, "r")) != NULL) {
    while (fscanf(listFile, "%s %s %s", ID, modelFileName, bgModelFileName) != EOF) {
      if (this->modelOrigin == sclib::modeSClib) {
        pBackgroundModel = this->pModeller->loadModel(bgModelFileName, sclib::mtGMM_new);
      } else {
        pBackgroundModel = loadWesleyGMM(bgModelFileName);
      }
      pActualCluster = constructModel(corpusFileName, ID, pBackgroundModel, modelFileName, layer, savePostFix);
      if (pActualCluster == NULL) {
        REPORT_ERROR(SVLIB_FileErr, "Error while building model")
      }
      MFree_0D(pBackgroundModel);

      if (count == 0) {
        pFirstCluster = pActualCluster;
        pClusterHook = pActualCluster;
      } else {
        pClusterHook->Next = pActualCluster;
        pClusterHook = pClusterHook->Next;
      }

      printf(".");
      count++;
    }
    fclose(listFile);
  } else {
		REPORT_ERROR(SVLIB_FileErr, "Corpus-file not found or could not be opened")
  }

  return pFirstCluster; 
}

//====================================================================================================================
//	load or train a single model, return a cluster-object
//  if the postfix is given, the trained models get saved immediately as if saved by savedModel()
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::constructModel(const char* corpusFileName, const char* targetID, SC_Model* pBackgroundModel, const char* modelFileName, unsigned long int layer, const char* savePostFix) {
  SC_Model *pModel = NULL, *pOriginalBackgroundModel = NULL;
  SC_Cluster *pCluster = constructFeatures(corpusFileName, targetID);

  if (pCluster != NULL) {
    //try to load the model to save the time that the EM algorithm consumes
    if (this->modelOrigin == sclib::modeSClib) {
      pModel = this->pModeller->loadModel(modelFileName, (layer == sclib::modeForeground) ? this->pTweak->modelHandler.foregroundModelType : this->pTweak->modelHandler.backgroundModelType);
    } else {
      if (layer == sclib::modeForeground && (this->pTweak->modelHandler.foregroundModelType == sclib::mtMIXMAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX || this->pTweak->modelHandler.foregroundModelType == sclib::mtMIX2MAX_ex)) {
        pModel = loadWesleyMIXMAX (modelFileName);
      } else {
        pModel = loadWesleyGMM(modelFileName);
      }
    }
    if (pModel != NULL && (pModel->Hdr.ModelType == sclib::mtMIXMAX || pModel->Hdr.ModelType == sclib::mtMIX2MAX || pModel->Hdr.ModelType == sclib::mtMIX2MAX_ex)) {
      pOriginalBackgroundModel = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pBackgroundModel));
      switch (pModel->Hdr.ModelType) {
        case sclib::mtMIXMAX: {
          ((SC_MixtureModel_MIXMAX*)pModel)->setOriginalBackground(((SC_MixtureModel*)pOriginalBackgroundModel));
          break;
        }
        case sclib::mtMIX2MAX: {
          ((SC_MixtureModel_MIX2MAX*)pModel)->setOriginalBackground(((SC_MixtureModel*)pOriginalBackgroundModel));
          break;
        }
        case sclib::mtMIX2MAX_ex: {
          ((SC_MixtureModel_MIX2MAXex*)pModel)->setOriginalBackground(((SC_MixtureModel*)pOriginalBackgroundModel));
          break;
        }
      }
    }

    //if loading the model failed or if there was no filename given to load, build the model from scratch
    if (pModel == NULL) { 
      pModel = this->pModeller->buildModel(pCluster->getSpeechFrames(), pBackgroundModel, layer, 1);
    }

    pCluster->setMergedModel(pModel); //creates a copy of the model!!!
    pCluster->setBackgroundModel(pBackgroundModel);

    if (((this->pTweak->debug.debugMode && sclib::dbSpeakerModels != 0) && layer == sclib::modeForeground) || ((this->pTweak->debug.debugMode && sclib::dbNoiseModels != 0) && layer != sclib::modeForeground)) {
      char *fileName = sclib::exchangeFileExtension(pCluster->getSpeakerName(), ".txt");
			sclib::classOut(fileName, pModel, this->pTweak);
      MFree_1D(fileName);
    }

		//try to save immediately
		if (strncmp(savePostFix, "", sclib::bufferSize) != 0) {
			saveModel(pCluster, savePostFix, false);
		}

    MFree_0D(pModel);
  } else {
    REPORT_ERROR(SVLIB_BadData, "Data couldn't be loaded");
  }

  return pCluster;
}

//====================================================================================================================
//	collect all the data for one singer and return it as a cluster-object
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::constructFeatures(const char* corpusFileName, const char* targetID) {
  SC_FeatureHandler *pHTKloader = new SC_FeatureHandler(this->pTweak);
  SV_Data *pMFCCstart = NULL, *pMFCChook = NULL, *pMFCC = NULL;
  SC_Cluster *pCluster = NULL;
  FILE *corpusFile;
  char fileName[sclib::bufferSize], identity[30], oldIdentity[30], *strippedFileName;
  long ID = sclib::noSpeaker;
  unsigned long corpusSize = 0, segStart = 0, segEnd = 0;
  bool finished = false;
  int scanStatus;

	if ((corpusFile = fopen(corpusFileName, "r")) != NULL) {
		while (finished == false) {
      corpusSize++;
      segEnd++;
      scanStatus = fscanf(corpusFile, "%s %s\n", fileName, identity);

      //id-change, so all the data of one singer is collected 
      if (((strcmp(identity, oldIdentity) != 0) && (corpusSize > 1)) || scanStatus == EOF) {
        ID++; //the id for the cluster-object
        
        //target idetified => merge the data, build the model, store it in the linked list of models (a cluster)
        if (strncmp(oldIdentity, targetID, 30) == 0) {
          pMFCC = pMFCCstart->MergeData();
          sclib::destructLinkedList(pMFCCstart);
          pCluster = new SC_Cluster(this->pTweak, pMFCC, &segStart, &segEnd, NULL, NULL, NULL, 1, ID, oldIdentity);
          MFree_0D(pMFCC);
          finished = true;
        } //id == target
        
        if (scanStatus == EOF) {
          finished = true;
        }

        segStart = segEnd + 1; 
      } //if id-change       

      if ((scanStatus != EOF) && (strncmp(identity, targetID, 30) == 0)) {
        if (pMFCCstart == NULL) {
          pMFCCstart =  pHTKloader->loadHTKfeatures(fileName);
          if (pMFCCstart == NULL) {
            continue;
          } else {
            pMFCChook = pMFCCstart;
          }
        } else {
          pMFCChook->Next = pHTKloader->loadHTKfeatures(fileName);
          if (pMFCChook->Next == NULL) {
            continue;
          } else {
            pMFCChook = pMFCChook->Next;
          }
        } //mfccstart == NULL

        //fill in the frame-parameters based on prior knowledge... ugly thing ;-)
        pMFCChook->Hdr.frameSize = this->pGT->getConverter()->ms2sample(32);
        pMFCChook->Hdr.frameStep = this->pGT->getConverter()->ms2sample(10);
        pMFCChook->Hdr.sampleRate = this->pGT->getAudioSampleRate();

        if (this->convertFeatures == true) { //convert to log-filterbank features, if wished
          mfcc2lfbe(pMFCChook);
        }

        if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbFeatures) == true) {
          sprintf(fileName, "%s.txt", fileName);
          strippedFileName = sclib::extractFileName(fileName);
          sclib::classOut(strippedFileName, pMFCChook, this->pTweak, ios_base::out);
          MFree_1D(strippedFileName);
        }

        strcpy(oldIdentity, identity);
      } //scanStatus != EOF
      
 		} //while fscanf()
  	fclose(corpusFile);
	} else {
		REPORT_ERROR(SVLIB_FileErr, "Corpus-file not found or could not be opened")
	}

  MFree_0D(pHTKloader);

  return pCluster;
}

//====================================================================================================================
//	builds a model for each ID mentioned in the targetListFile and found in the corpusFile; returns a linked list of
//  clusters containing the models and the data
//  if savePostFix != '', the models are, not only trained, but also saved immediately after creation (sving via 
//  saveModels() needs all models trained before saving the first, so there is lot of unsafed work if training is 
//  slow...)
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::trainModels(const char* targetListFileName, const char* corpusFileName, unsigned long int layer, const char* savePostFix) {
  return constructModels(targetListFileName, corpusFileName, layer, savePostFix);
}

//====================================================================================================================
//	reads the corpus-file and builds the model for the specified target; returns a cluster-object containing the model
//  and the training-data
//  if modelFileName != NULL, it is tried to load the model from this file instead of buildig it from scratch
//  if savePostFix != '', the models are not only trained, but also saved immediately after creation 
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::trainModel(const char* corpusFileName, const char* targetID, SC_Model* pBackgroundModel, unsigned long int layer, const char* savePostFix) {
  return constructModel(corpusFileName, targetID, pBackgroundModel, NULL, layer, savePostFix);
}

//====================================================================================================================
//	just a shortcut to constructModel(), because nearly the same functionality is needed in both functions
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::loadModel(const char* corpusFileName, const char* targetID, const char* modelFileName, SC_Model* pBackgroundModel, unsigned long int layer) {
  return constructModel(corpusFileName, targetID, pBackgroundModel, modelFileName, layer);
}

//====================================================================================================================
//	fill the cluster-objects as in trainModels(), but load the models from saved files instead of rebuilding them
//  if loading a model fails, it is rebuilded from the feature vectors as in trainModels()
//====================================================================================================================
SC_Cluster* SC_Corpus_Wesley::loadModels(const char* targetListFileName, const char* corpusFileName, unsigned long int layer) {
  return constructModels(targetListFileName, corpusFileName, layer);
}

//====================================================================================================================
//	write the models to binary files; the filenames get the given postfix (should include the '.')
//====================================================================================================================
long SC_Corpus_Wesley::saveModels(SC_Cluster* pModels, const char* postfix) {
  SC_Cluster *pHook = pModels;
  long bytesWritten = 0;
  
  while (pHook != NULL) {
		bytesWritten += saveModel(pHook, postfix, true);
    pHook = pHook->Next;
  }

  return bytesWritten;
}

//====================================================================================================================
//	write the model to a binary file; the filename gets the given postfix (should include the '.')
//====================================================================================================================
long SC_Corpus_Wesley::saveModel(SC_Cluster* pModelWithInfo, const char* postfix, bool verbose) {
  SC_Model *pModel = pModelWithInfo->getMergedModel();
  long res;
  char fileName[200];

	sprintf(fileName, "%s%s%s\0", this->pTweak->debug.debugDir, pModelWithInfo->getSpeakerName(), postfix);
  
  pModel->OpenFile(fileName, WRITE_MODEL);
	res = pModel->SaveModel();
  if (res <= 0) {
    REPORT_ERROR(SVLIB_FileErr, "Model could not be written");
  }
	pModel->CloseFile();

	if (verbose == true) {
		printf(".");
	}

	return res;
}

//====================================================================================================================
//	test the prebuild models against the files test-corpus
//====================================================================================================================
void SC_Corpus_Wesley::testModels(const char* targetListFileName, const char* corpusFileName, SC_Cluster* pTrainedModels, const char* resultFileName) {
  SV_Data *pMFCC = NULL;
  SC_Cluster *pCluster, *pBestModel, *pHook;
  SC_Model *pModel = NULL, *pBackgroundModel = NULL;
  FILE *listFile, *corpusFile, *resultFile;
  char identity[30], modelFileName[sclib::bufferSize], bgModelFileName[sclib::bufferSize], resFileName[sclib::bufferSize];
	char otherResults[10*sclib::bufferSize], tmp[sclib::bufferSize];
  double res, bestResult;
  unsigned int correct = 0, overall = 0;
  
  //bool printed = false;
  
  sprintf(resFileName, "%s%s\0", this->pTweak->debug.debugDir, resultFileName);
  if ((resultFile = fopen(resFileName, "w")) == NULL) {
    REPORT_ERROR(SVLIB_FileErr, "Result file could not be written");
    return;
  }
  fclose(resultFile);

	if ((listFile = fopen(targetListFileName, "r")) != NULL) {
	  if ((corpusFile = fopen(corpusFileName, "r")) != NULL) {
		
      while (fscanf(listFile, "%s %s %s\n", identity, modelFileName, bgModelFileName) != EOF) {
        pCluster = constructFeatures(corpusFileName, identity);
        if (this->modelOrigin == sclib::modeSClib) {
          pBackgroundModel = this->pModeller->loadModel(bgModelFileName, sclib::mtGMM_new);
        } else {
          pBackgroundModel = loadWesleyGMM(bgModelFileName);
        }
        if (pCluster != NULL) {
          pHook = pTrainedModels;
          pBestModel = NULL;
          bestResult = (-1.0) * std::numeric_limits<double>::max();
					sprintf(otherResults, "");

          while (pHook != NULL) {
            pModel = pHook->getMergedModel();
            res = this->pModeller->testModel(pModel, pCluster->getSpeechFrames(), pCluster->getSegmentCount(), false, pBackgroundModel);
            if (res > bestResult) {
              bestResult = res;
              pBestModel = pHook;
            }
						sprintf(tmp, "[%s %f]", pHook->getSpeakerName(), res);
						strcat(otherResults, tmp);
						//sprintf(otherResults, "[%s %f] %s", pHook->getSpeakerName(), res, otherResults);
            pHook = pHook->Next;
            printf(".");
          }

          if (pBestModel != NULL) {
            resultFile = fopen(resFileName, "a");
            fprintf(resultFile, "%s %s %s %f (%s)\n", modelFileName, identity, pBestModel->getSpeakerName(), bestResult, otherResults);
            overall++;
            if (strncmp(identity, pBestModel->getSpeakerName(), 2) == 0) {
              correct++;
            }
          } else {
            fprintf(resultFile, "%s %s %s %f\n", modelFileName, identity, "!!!ERROR!!!", bestResult);
          }
          fclose(resultFile);
        }
        
        MFree_0D(pBackgroundModel);
        MFree_0D(pCluster);
        printf(":");

      } //while fscanf()

  	  fclose(corpusFile);
    } else {
      REPORT_ERROR(SVLIB_FileErr, "Corpus-file not found or could not be opened")
    }
    fclose(listFile);
	} else {
		REPORT_ERROR(SVLIB_FileErr, "Modellist-file not found or could not be opened")
	}

  resultFile = fopen(resFileName, "a");
  fprintf(resultFile, "=> Result: %i / %i = %f%", correct, overall, (double)(correct) / (double)(overall) * 100.0);
  fclose(resultFile);

  return;
}

//====================================================================================================================
//	convert MFCC features to log-FbE features by invoking an inverse dct on the MFCCs
//====================================================================================================================
void SC_Corpus_Wesley::mfcc2lfbe(SV_Data* pMFCC) {
  long int y, x;
  double *idctVec = NULL;
  double *mfccVec = new double[pMFCC->Col];
  SC_Transform *pTrans = new SC_Transform(this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.window, this->pTweak->transform.taperingLength);
  
  for (y = 0; y < pMFCC->Row; y++) {
    for (x = 0; x < pMFCC->Col; x++) {
      mfccVec[x] = pMFCC->Mat[y][x];
    }
  
    //sclib::vectorOut("mfcc.txt", mfccVec, pMFCC->Col);
    idctVec = pTrans->idct(mfccVec, pMFCC->Col);
    //sclib::vectorOut("lfbe.txt", idctVec, pMFCC->Col);    

    for (x = 0; x < pMFCC->Col; x++) {
      pMFCC->Mat[y][x] = (float)(idctVec[x]);
    }
    MFree_1D(idctVec);
  }

  MFree_0D(pTrans);
  MFree_1D(mfccVec);

  return;
}

//====================================================================================================================
//	load a GMM model saved by Wesleys algorithms
//  (is always of type sclib::mtGMM_new, not of type sclib::mtBGMM or what else could be specified in 
//  this->pTweak->modelType)
//  TODO: make it machine-independant (endianness, 64/32bit), also in Wesleys programms then
//====================================================================================================================
SC_Model* SC_Corpus_Wesley::loadWesleyGMM(const char* fileName) {
  double *newMean, *newVar, *newWeight;
  int i, j, newDim, newMixtureCount;
  SC_MixtureModel_GMM *pModel = NULL;
  FILE *modelFile = fopen(fileName,"rb");

  if (modelFile != NULL) {
		fread(&newDim, sizeof(int), 1, modelFile);
		fread(&newMixtureCount, sizeof(int), 1, modelFile);

    pModel = new SC_MixtureModel_GMM(this->pTweak, newMixtureCount, newDim);

    // set new mean
    for (i = 0; i < newMixtureCount; i++) {
      newMean = pModel->getMean(i);
      fread(newMean, sizeof(double), newDim, modelFile);
    }

    //set new variance and sd
    for (i = 0; i < newMixtureCount; i++) {
      newVar = pModel->getVariance(i);
      fread(newVar, sizeof(double), newDim, modelFile);
      for (j = 0; j < newDim; j++) {
        pModel->setSd(i, j, sqrt(newVar[j]));
      }
    }

    //set new weight
    newWeight = pModel->getWeight();
		fread(newWeight, sizeof(double), newMixtureCount, modelFile);

    //set additional members
    pModel->setTrainindDataCount(1); //because we don't know...

    fclose(modelFile);
	}
  
  return pModel;
}

//====================================================================================================================
//	load a MIXMAX model saved by Wesleys algorithms
//  TODO: make it machine-independant (endianness, 64/32bit), also in Wesleys programms then
//====================================================================================================================
SC_Model* SC_Corpus_Wesley::loadWesleyMIXMAX(const char* fileName) {
  double *newMean, *newVar, *newWeight;
  int i, j, newDim, newMixtureCount;
  SC_MixtureModel *pModel = NULL;
  SC_MixtureModel_GMM *pDummy = NULL;
  FILE *modelFile = fopen(fileName,"rb");
  SC_ModelHandler *pHandler = new SC_ModelHandler(this->pTweak);

  if (modelFile != NULL) {
		fread(&newDim, sizeof(int), 1, modelFile);
		fread(&newMixtureCount, sizeof(int), 1, modelFile);

    pModel = (SC_MixtureModel*)pHandler->createRawModel(this->pTweak->modelHandler.foregroundModelType, NULL, newMixtureCount, newDim);

    //set new mean and maskLevel
    for (i = 0; i < newMixtureCount; i++) {
      newMean = pModel->getMean(i);
      fread(newMean, sizeof(double), newDim, modelFile);
    }

    //set new variance and sd
    for (i = 0; i < newMixtureCount; i++) {
      newVar = pModel->getVariance(i);
      fread(newVar, sizeof(double), newDim, modelFile);
      for (j = 0; j < newDim; j++) {
        pModel->setSd(i, j, sqrt(newVar[j]));
      }
    }

    //set new weight
    newWeight = pModel->getWeight();
		fread(newWeight, sizeof(double), newMixtureCount, modelFile);

    //set additional members
    pModel->setTrainindDataCount(1); //because we don't know...

    //Set a dummy as the original background
    pDummy = new SC_MixtureModel_GMM(this->pTweak, 0, newDim);
    switch (this->pTweak->modelHandler.foregroundModelType) {
      case sclib::mtMIXMAX: {
        for (i = 0; i < newMixtureCount; i++) {
          ((SC_MixtureModel_MIXMAX*)pModel)->setMaskLevel(i, 0.0); 
        }
        ((SC_MixtureModel_MIXMAX*)pModel)->setNoiseCorruptionType(sclib::nctMax);
        ((SC_MixtureModel_MIXMAX*)pModel)->setOriginalBackground(pDummy);
        break;
      }
      case sclib::mtMIX2MAX: {
        for (i = 0; i < newMixtureCount; i++) {
          ((SC_MixtureModel_MIX2MAX*)pModel)->setMaskLevel(i, 0.0); 
        }
        ((SC_MixtureModel_MIX2MAX*)pModel)->setOriginalBackground(pDummy);
        break;
      }
      case sclib::mtMIX2MAX_ex: {
        for (i = 0; i < newMixtureCount; i++) {
          ((SC_MixtureModel_MIX2MAXex*)pModel)->setMaskLevel(i, 0.0); 
        }
        ((SC_MixtureModel_MIX2MAXex*)pModel)->setOriginalBackground(pDummy);
        break;
      }
    }
    
    fclose(modelFile);
	}

  MFree_0D(pHandler);
  
  return pModel;
}

//====================================================================================================================
//	alter the modelOrigin-member; return true, if new value is valid, otherwise return false and don't change anything
//====================================================================================================================
bool SC_Corpus_Wesley::setModelOrigin(unsigned long int newOrigin) {
  if (newOrigin == sclib::modeForeign || newOrigin == sclib::modeSClib) {
    this->modelOrigin = newOrigin;
    return true;
  } else {
    return false;
  }
}

//====================================================================================================================
// change the sampleRate
//====================================================================================================================
void SC_Corpus_Wesley::setSampleRate(unsigned long int sampleRate) {
  ((SC_GroundTruth_Wesley*)this->pGT)->setSampleRate(sampleRate);
  return;
}

