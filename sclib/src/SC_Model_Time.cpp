/**************************************************************************/
/*    This is the testbed for a new kind of model that also regards the   */
/*    time order of feature vectors in the given feature sets             */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.02.2009																								*/
/**************************************************************************/

#include <map>
#include "SC_Model_Time.h"
#include "SC_Aux.h"
#include "SC_ModelHandler.h"
#include "SC_FeatureHandler.h"
#include "SC_Model_SVM.h"
#include "SC_Synthesis.h"

SC_Model* SC_Model_Time::pWorldModel = NULL; //to save memory, all instances share one world model; TODO: change if usecase changes!
int SC_Model_Time::instanceCounter = 0;

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_Time::SC_Model_Time(SC_TweakableParameters *pTweak, unsigned int dim, unsigned int syllableLength, unsigned int trajectoryStep, char subModelType, bool removeTiming, unsigned int templateCount, unsigned int clusteringIterations, bool replaceTrainingData, bool checkForTrajectorization, const char *worldModelFile, const char *normalizationFile, bool cacheResults, bool verbose) : SC_Model(pTweak) {
	this->Hdr.ModelType = sclib::mtTime;
	this->subModelType = subModelType;
	this->verbose = verbose;
	this->dim = dim;
	this->syllableLength = syllableLength;
	this->trajectoryStep = trajectoryStep;
	this->pSubModel = NULL;
	this->removeTiming = removeTiming;
	this->templateCount = templateCount;
	this->clusteringIterations = clusteringIterations;
	this->replaceTrainingData = replaceTrainingData; 
	this->checkForTrajectorization = checkForTrajectorization;
	this->cacheResults = cacheResults;
	SC_Model_Time::instanceCounter++;
	if (SC_Model_Time::pWorldModel==NULL && sclib::fileExists(worldModelFile)==true) {
		SC_ModelHandler loader(this->pTweak, false);
		SC_Model_Time::pWorldModel = loader.loadModel(worldModelFile, this->subModelType);
	}// else {
	//	this->pWorldModel = NULL;
	//}
	if (sclib::fileExists(normalizationFile) == true) {
		SC_FeatureHandler loader(this->pTweak, false);
		this->pNorm = loader.loadFeature(normalizationFile, 0);
	} else {
		this->pNorm = NULL;
	}
}

//====================================================================================================================
// copy-constructor
// the normalization matrix pNorm (as created by the feature handler's createNormalizationMatrix() method) is stored
// and freed inside this class, although (maybe) created externally!
//====================================================================================================================
SC_Model_Time::SC_Model_Time(SC_Model_Time& pParent, bool justLink) : SC_Model(pParent) {
	SC_ModelHandler handler(this->pTweak, false);
	
	this->subModelType = pParent.subModelType;
	this->verbose = pParent.verbose;
	this->dim = pParent.dim;
	this->syllableLength = pParent.syllableLength;
	this->trajectoryStep = pParent.trajectoryStep;
	this->removeTiming = pParent.removeTiming;
	this->templateCount = pParent.templateCount;
	this->clusteringIterations = pParent.clusteringIterations;
	this->replaceTrainingData = pParent.replaceTrainingData; 
	this->checkForTrajectorization = pParent.checkForTrajectorization;
	this->cacheResults = pParent.cacheResults;
	if (pParent.pNorm != NULL) {
		this->pNorm = new SV_Data(*pParent.pNorm, false);
	} else {
		this->pNorm = NULL;
	}

	if (pParent.pSubModel != NULL) {
		this->pSubModel = handler.copyModel(pParent.pSubModel, justLink);
	} else {
		this->pSubModel = NULL;
	}
	SC_Model_Time::instanceCounter++;
	//if (pParent.pWorldModel != NULL) {
	//	this->pWorldModel = handler.copyModel(pParent.pWorldModel, justLink);
	//} else {
	//	this->pWorldModel = NULL;
	//}
	std::map<std::string, double>::const_iterator i;
	for (i = pParent.checksumCache.map.begin(); i != pParent.checksumCache.map.end(); ++i) {
		this->checksumCache.map[i->first] = i->second;
	}
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_Time::~SC_Model_Time() {
	SC_Model_Time::instanceCounter--;
	if (SC_Model_Time::instanceCounter == 0) {
		MFree_0D(SC_Model_Time::pWorldModel);
	}
	MFree_0D(this->pSubModel);
	MFree_0D(this->pNorm);
	this->checksumCache.map.clear();
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_Time& SC_Model_Time::operator=(SC_Model_Time& pParent) {
	SC_ModelHandler handler(this->pTweak, false);

	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		this->subModelType = pParent.subModelType;
		this->verbose = pParent.verbose;
		this->dim = pParent.dim;
		this->syllableLength = pParent.syllableLength;
		this->trajectoryStep = pParent.trajectoryStep;
		this->removeTiming = pParent.removeTiming;
		this->templateCount = pParent.templateCount;
		this->clusteringIterations = pParent.clusteringIterations;
		this->replaceTrainingData = pParent.replaceTrainingData; 
		this->checkForTrajectorization = pParent.checkForTrajectorization;
		this->cacheResults = pParent.cacheResults;
		if (pParent.pNorm != NULL) {
			this->pNorm = new SV_Data(*pParent.pNorm, false);
		} else {
			this->pNorm = NULL;
		}

		if (pParent.pSubModel != NULL) {
			this->pSubModel = handler.copyModel(pParent.pSubModel);
		} else {
			this->pSubModel = NULL;
		}
		SC_Model_Time::instanceCounter++;
		//if (pParent.pWorldModel != NULL) {
		//	this->pWorldModel = handler.copyModel(pParent.pWorldModel);
		//} else {
		//	this->pWorldModel = NULL;
		//}

		std::map<std::string, double>::const_iterator i;
		for (i = pParent.checksumCache.map.begin(); i != pParent.checksumCache.map.end(); ++i) {
			this->checksumCache.map[i->first] = i->second;
		}
	}

  return *this;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_Time::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_Time*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_Time*)pBackgroundModels));
}
SC_Model_Time* SC_Model_Time::combineModels(SC_Model_Time* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_Time* pBackgroundModels) {
	SC_Model_Time *pNewModel = new SC_Model_Time(this->pTweak, this->dim);

	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  return pNewModel;
}

//====================================================================================================================
// Train a model; here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_Time::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Train a model; here, one can specify how many segments of a linked list of feature-vectors should be used
// (0 means all)
//====================================================================================================================
int	SC_Model_Time::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	int t, d, res = SVLIB_Ok, segmentCounter = 0, T;
	SC_ModelHandler modeler(this->pTweak, this->verbose);
	SC_Conversion converter(pData->Hdr.sampleRate);
	unsigned long int trajectoryLength = converter.ms2audioFrame(this->syllableLength, pData->Hdr.frameSize, pData->Hdr.frameStep, sclib::alignmentStart);
  SV_Data *pCompleteData, *pHook = pData, *pTrajectories = NULL, *pTrajectoriesHook, *pTmp;
	SC_FeatureHandler extractor(this->pTweak, false);

	MFree_0D(this->pSubModel);
	this->trainingDataCount = 0;
	this->dim = pData->Col;
	this->checksumCache.map.clear();

	if (this->checkForTrajectorization==false || pData->Hdr.Signature[2]>0) {
		pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
		this->trainingDataCount = pCompleteData->Row;
	} else {
		do { //create trajectories per segment, then merge the trajectories so that no trajectorie is built that didn't exist in the features (because of artificial concatenation)
			this->trainingDataCount += pHook->Row;

			//create trajectory vectors of given syllable-length
			pTmp = pHook;
			if (this->pNorm != NULL) {
				pTmp = extractor.normalize(pTmp, this->pNorm, this->replaceTrainingData);
			}
			if (pTrajectories == NULL) {
				pTrajectories = extractor.createTrajectories(pTmp, trajectoryLength, this->trajectoryStep, this->removeTiming, this->templateCount, this->clusteringIterations);
				pTrajectoriesHook = pTrajectories;
			} else {
				pTrajectoriesHook->Next = extractor.createTrajectories(pTmp, trajectoryLength, this->trajectoryStep, this->removeTiming, this->templateCount, this->clusteringIterations);
				pTrajectoriesHook = pTrajectoriesHook->Next;
			}

			if (this->replaceTrainingData == true) {
				pTmp->Hdr = pTrajectoriesHook->Hdr;
				pTmp->Row = pTrajectoriesHook->Row;
				pTmp->Col = pTrajectoriesHook->Col;
				if (pTmp->getJustLinked() == false) {
					MFree_2D(pTmp->Mat);
				}
				pTmp->Mat = pTrajectoriesHook->Mat;
				pTmp->setJustLinked(false);
				pTrajectoriesHook->setJustLinked(true);
			}
			if (pHook != pTmp) {
				MFree_0D(pTmp);
			}

			pHook = pHook->Next;
			segmentCounter++;
		} while (pHook!=NULL && (unsigned long int)(segmentCounter)<segmentsToMerge);
		pCompleteData = pTrajectories->MergeData();
		sclib::destructLinkedList(pTrajectories);
	}

	//remove common parts of speech if a world model is given
	if (SC_Model_Time::pWorldModel!=NULL && this->subModelType==sclib::mtSVM) { //TODO: enhance implementation to fit for other model types
		((SC_Model_SVM*)SC_Model_Time::pWorldModel)->setDistanceBasedTesting(false); //make sure we get a result per data row, not just an overall one
		pTmp = SC_Model_Time::pWorldModel->TestModel(pCompleteData);
		T = 0;
		for (t = 1; t < pTmp->Row; t++) { //count the number of uncommon trajectories
			if (sclib::round(pTmp->Mat[t][0]) != sclib::labelPositive) {
				T++;
			}
		}
		if (this->verbose == true) {
			printf("\n      Time-Model: %f%% (%d/%d) training data removed by world model", (double)((pCompleteData->Row-T)*100.0)/(double)(pCompleteData->Row), pCompleteData->Row-T, pCompleteData->Row);
		}
		pHook = new SV_Data(T, pCompleteData->Col);
		pHook->Hdr = pCompleteData->Hdr;
		T = 0;
		for (t = 1; t < pTmp->Row; t++) { //copy the non-common trajectories to a new set
			if (sclib::round(pTmp->Mat[t][0]) != sclib::labelPositive) {
				for (d = 0; d < pCompleteData->Col; d++) {
					pHook->Mat[T][d] = pCompleteData->Mat[t-1][d];
				}
				T++;
			}
		}
		MFree_0D(pTmp);

		/*sclib::toggle(this->pTweak->featureMfcc.dEnergy, false, false);
		SC_Synthesis synthesizer(this->pTweak);
		SC_FeatureHandler handler(this->pTweak);
		pTmp = handler.unwindTrajectories(pCompleteData, trajectoryLength, this->trajectoryStep);
		handler.unNormalize(pTmp, this->pNorm, true);
		SV_Data *pTmp2 = new SV_Data(pTmp->Row, pTmp->Col);
		pTmp2->Hdr = pTmp->Hdr;
		pTmp2->Hdr.ID = sclib::featureMFCC;
		for (t = 0; t < pTmp->Row; t++) {
			pTmp2->Mat[t][0] = 10.0;
			for (d = 0; d < pTmp->Col-1; d++) {
				pTmp2->Mat[t][d+1] = pTmp->Mat[t][d];
			}
		}
		synthesizer.feature2wav("orig.wav", pTmp2);
		MFree_0D(pTmp);
		MFree_0D(pTmp2);
		pTmp = handler.unwindTrajectories(pHook, trajectoryLength, this->trajectoryStep);
		handler.unNormalize(pTmp, this->pNorm, true);
		pTmp2 = new SV_Data(pTmp->Row, pTmp->Col);
		pTmp2->Hdr = pTmp->Hdr;
		pTmp2->Hdr.ID = sclib::featureMFCC;
		for (t = 0; t < pTmp->Row; t++) {
			pTmp2->Mat[t][0] = 10.0;
			for (d = 0; d < pTmp->Col-1; d++) {
				pTmp2->Mat[t][d+1] = pTmp->Mat[t][d];
			}
		}
		synthesizer.feature2wav("worldremoved.wav", pTmp2);
		MFree_0D(pTmp);
		MFree_0D(pTmp2);
		sclib::toggle(this->pTweak->featureMfcc.dEnergy, true);*/

		if (pCompleteData != pData) {
			MFree_0D(pCompleteData);
		}
		pCompleteData = pHook; //make the new set the complete set
	}

	//create and train submodel with training trajectories
	//sclib::toggle(this->pTweak->modelSvm.doParameterSearch, false, true);
	bool oldDbt = this->pTweak->modelSvm.distanceBasedTesting;
	this->pTweak->modelSvm.distanceBasedTesting = false;
	this->pSubModel = modeler.createRawModel(this->subModelType, NULL, 0, this->dim);
	res = this->pSubModel->TrainModel(pCompleteData);
	//sclib::toggle(this->pTweak->modelSvm.doParameterSearch, true);
	this->pTweak->modelSvm.distanceBasedTesting = oldDbt;
	MFree_0D(pTrajectories);

	if (pCompleteData != pData) {
		MFree_0D(pCompleteData);
	}

  return res;
}

//====================================================================================================================
// Test the model while merging no segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_Time::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;
}

//====================================================================================================================
// Test the model with specified nr. of segments in the linked list of feature-vectors.
//====================================================================================================================
SV_Data* SC_Model_Time::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
	SV_Data *pHook = TestData, *pCompleteData, *pTrajectories = NULL, *pTrajectoriesHook, *pScore = NULL, *pTmp;
	SC_Conversion converter(TestData->Hdr.sampleRate);
	SC_FeatureHandler extractor(this->pTweak, false);
	unsigned long int trajectoryLength = converter.ms2audioFrame(this->syllableLength, TestData->Hdr.frameSize, TestData->Hdr.frameStep, sclib::alignmentStart);
	int t, T, d, segmentCounter = 0;
	std::string checksum;

	if (this->cacheResults == true) {
		char *tmp = extractor.getChecksum(TestData, segmentsToMerge);
		checksum = tmp;
		MFree_1D(tmp);
		std::map<std::string, double>::const_iterator i = this->checksumCache.map.find(checksum);
		if (i != this->checksumCache.map.end()) {
			pScore = new SV_Data(1, 1);
			pScore->Mat[0][0] = (float)(i->second);
			return pScore;
		}
	}

	if (this->pSubModel != NULL) {
		if (this->checkForTrajectorization==true && TestData->Hdr.Signature[2]>0) {
			pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge);
		} else {
			do { //create trajectories per segment, then merge the trajectories so that no trajectorie is built that didn't exist in the features (because of artificial concatenation)
				//create trajectory vectors of given syllable-length
				pTmp = pHook;
				if (this->pNorm != NULL) {
					pTmp = extractor.normalize(pTmp, pNorm, false);
				}
				if (pTrajectories == NULL) {
					pTrajectories = extractor.createTrajectories(pTmp, trajectoryLength, this->trajectoryStep, this->removeTiming, this->templateCount, this->clusteringIterations);
					pTrajectoriesHook = pTrajectories;
				} else {
					pTrajectoriesHook->Next = extractor.createTrajectories(pTmp, trajectoryLength, this->trajectoryStep, this->removeTiming, this->templateCount, this->clusteringIterations);
					pTrajectoriesHook = pTrajectoriesHook->Next;
				}

				if (pHook != pTmp) {
					MFree_0D(pTmp);
				}

				pHook = pHook->Next;
				segmentCounter++;
			} while (pHook!=NULL && (unsigned long int)(segmentCounter)<segmentsToMerge);
			pCompleteData = pTrajectories->MergeData();
			sclib::destructLinkedList(pTrajectories);
		}

		//remove common parts of speech if a world model is given
		if (SC_Model_Time::pWorldModel!=NULL && this->subModelType==sclib::mtSVM) { //TODO: enhance implementation to fit for other model types
			((SC_Model_SVM*)SC_Model_Time::pWorldModel)->setDistanceBasedTesting(false); //make sure we get a result per data row, not just an overall one
			pScore = SC_Model_Time::pWorldModel->TestModel(pCompleteData);
			T = 0;
			for (t = 1; t < pScore->Row; t++) { //count the number of uncommon trajectories
				if (sclib::round(pScore->Mat[t][0]) != sclib::labelPositive) {
					T++;
				}
			}
			if (this->verbose == true) {
				printf("\n      Time-Model: %f%% (%d/%d) test data removed by world model", (double)((pCompleteData->Row-T)*100.0)/(double)(pCompleteData->Row), pCompleteData->Row-T, pCompleteData->Row);
			}
			pHook = new SV_Data(T, pCompleteData->Col);
			pHook->Hdr = pCompleteData->Hdr;
			T = 0;
			for (t = 1; t < pScore->Row; t++) { //copy the non-common trajectories to a new set
				if (sclib::round(pScore->Mat[t][0]) != sclib::labelPositive) {
					for (d = 0; d < pCompleteData->Col; d++) {
						pHook->Mat[T][d] = pCompleteData->Mat[t-1][d];
					}
					T++;
				}
			}
			MFree_0D(pScore);
			if (pCompleteData != TestData) {
				MFree_0D(pCompleteData);
			}
			pCompleteData = pHook; //make the new set the complete set
		}

		//finally, score the data
		pScore = this->pSubModel->TestModel(pCompleteData);
		MFree_0D(pTrajectories);

		if (this->cacheResults == true) {
			this->checksumCache.map[checksum] = pScore->Mat[0][0];
		}
	}

	if (pCompleteData != TestData) {
		MFree_0D(pCompleteData);
	}

	return pScore;
}

//====================================================================================================================
// load a model from file
//====================================================================================================================
SV_Model*	SC_Model_Time::LoadModel(void) {
	int res;
	bool modelAvailable;
	char *buffer = NULL, extension[sclib::bufferSize];
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);
	SC_ModelHandler handler(this->pTweak, false);

	this->checksumCache.map.clear();
	
  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Failed!");}
	
	//read own members
	io.readScalar(&(this->DFile), this->subModelType, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->verbose, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->dim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->syllableLength, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->trajectoryStep, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->removeTiming, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->templateCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->clusteringIterations, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->replaceTrainingData, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->checkForTrajectorization, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), this->cacheResults, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time Model Failed!");}
	io.readScalar(&(this->DFile), modelAvailable, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_ANN Model Failed!");}

	//read sub-models, if available
	if (modelAvailable == true) {
		MFree_0D(this->pSubModel);
		sprintf(extension, ".time_model");
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, extension);
		this->pSubModel = handler.createRawModel(this->subModelType, NULL, 0, 1);
		this->pSubModel->OpenFile(buffer, READ_MODEL);
		this->pSubModel = (SC_Model*)(this->pSubModel->LoadModel());
		if (this->pSubModel != NULL) {
			this->pSubModel->CloseFile();
		} else {
			REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Time sub-model Failed!");
		}
		MFree_1D(buffer);
	}

	return(this);
}

//====================================================================================================================
// save a model to file
//====================================================================================================================
int	SC_Model_Time::SaveModel(void) {
 	int res, bytes, cnt = 0;
	SV_DataIO io;
	char *buffer, extension[sclib::bufferSize];

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time model Failed!");}

  //write trainingDataCount
	bytes = io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}

  //write own stuff
	bytes += io.writeScalar(&(this->DFile), this->subModelType);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->verbose);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->dim);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->syllableLength);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->trajectoryStep);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->removeTiming);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->templateCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->clusteringIterations);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->replaceTrainingData);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->checkForTrajectorization);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->cacheResults);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), (this->pSubModel!=NULL)); //store if a sub-model file (.time_model) is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time Model Failed!");}

	//write sub-model file
	if (this->pSubModel != NULL) {
		sprintf(extension, ".time_model");
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, extension);
		this->pSubModel->OpenFile(buffer, WRITE_MODEL);
		res = this->pSubModel->SaveModel();
		if (res <= 0) {
			REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Time sub-model Failed!");
		} else {
			bytes += res;
		}
		this->pSubModel->CloseFile();
		MFree_1D(buffer);
	}

  return bytes + MHLen; //MHLen may be incorrect...
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_Time::modelOut(ostream& os) {
	os << "Syllable length [ms]:\t" << this->syllableLength << endl;
	os << "Trajectory step [#frames]:\t" << this->trajectoryStep << endl;
	os << "Remove timing:\t" << this->removeTiming << endl;
	os << "Nr. of Templates:\t" << this->templateCount << endl;
	os << "Clustering iterations:\t" << this->clusteringIterations << endl;
	os << "Replace trainign data:\t" << this->replaceTrainingData << endl;
	os << "Check for trajectorization:\t" << this->checkForTrajectorization << endl;
	os << "World-Model exists:\t" << (SC_Model_Time::pWorldModel != NULL) << endl;
	os << "Sub-Model exists:\t" << (this->pSubModel != NULL) << endl;
	os << "Sub-Model type:\t" << this->subModelType << endl;
	os << "Traing-Data count:\t" << this->trainingDataCount << endl;
	os << "Traing-Data dimensionality:\t" << this->dim << endl;

	os << endl;
	os << this->pSubModel << endl;
	os << endl;

	return(os);
}

//====================================================================================================================
// for computing BIC etc.
//====================================================================================================================
unsigned int SC_Model_Time::getFreeParameterCount(void) {
	return this->pSubModel->getFreeParameterCount();
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_Time::drawSamplesFromDistribution(unsigned long int count) {
	SV_Data *pSamples = NULL;

	if (this->pSubModel != NULL) {
		pSamples = new SV_Data(count, this->dim);

		SV_Data *pTrajectories = this->pSubModel->drawSamplesFromDistribution(count); //that'll be too much trajectories, but we don't know the number of frames per trajectory yet

		int tt = 0;
		for (int t = 0; t < pTrajectories->Row; t++) {
			for (unsigned int x = 0; x < pTrajectories->Col/this->dim; x++) {
				if (tt < pSamples->Row) {
					for (unsigned int d = 0; d < this->dim; d++) {
						pSamples->Mat[tt][d] = pTrajectories->Mat[t][x*this->dim +d];
					}
					tt++;
				}
			}
		}

		MFree_0D(pTrajectories);

		if (this->pNorm != NULL) {
			SC_FeatureHandler handler(this->pTweak);
			handler.unNormalize(pSamples, this->pNorm, true);
		}
	}
	
	return pSamples; 
}

//====================================================================================================================
// generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
// method needs to return false;
//====================================================================================================================
bool SC_Model_Time::scoreIsAverage(void) {
	SC_ModelHandler handler(this->pTweak, false);

	return handler.averagesItsScores(this->subModelType);
}
