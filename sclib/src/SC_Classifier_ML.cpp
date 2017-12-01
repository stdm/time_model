/**************************************************************************/
/*    Responsibility:																											*/
/*      - conatins a maximum-likelihood classifier based on one of the    */
/*        implemented SC_Model classes                                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.03.2006																								*/
/**************************************************************************/

#include "SC_Classifier_ML.h"
#include "SC_ModelHandler.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Classifier_ML::SC_Classifier_ML(SC_TweakableParameters* pTweak, bool doScaling) : SC_Classifier(pTweak, doScaling) {
	this->classifierType = sclib::ctML;
  this->pModels = NULL;
  this->modelCount = 0;
  this->pMapping = NULL;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Classifier_ML::~SC_Classifier_ML() {
  freeModels();
  sclib::destructLinkedList(this->pMapping);
}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//====================================================================================================================
int SC_Classifier_ML::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative) {
  unsigned short int order;
  SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak);
  SV_Data *pData = NULL;

  //destroy previously trained classifier (and corrsponding data), if any
  MFree_0D(this->pScale);
  freeModels();
  sclib::destructLinkedList(this->pMapping);
  this->isTrained = false;
	this->classCount = 0;

  //find scaling parameters
  if (this->doScaling == true) {
    SV_Data *pHook = pPositive->Next; //save next-pointer
    SV_Data *pCompleteData;
    pPositive->Next = pNegative;
    pCompleteData = pPositive->MergeData(2);
    MFree_0D(this->pScale);
    this->pScale = findScalingParameters(pCompleteData);
    MFree_0D(pCompleteData);
    pPositive->Next = pHook; //restore next-pointer
  }

  this->modelCount = 2;

  MArray_1D(this->pModels, 2, SC_Model*, "SC_Classifier_ML.trainTwoClass: pModels");
  
  //build model for negative samples
  pData = scaleFeatures(pNegative, this->pScale, -1, true);
  this->pModels[0] = NULL;
  order = pModelHandler->guessModelOrder(pData, NULL, this->pTweak->classifierMl.modelType, this->pTweak->classifierMl.minOrder, this->pTweak->classifierMl.maxOrder);
 	this->pModels[0] = pModelHandler->buildModel(pData, NULL, order, this->pTweak->classifierMl.modelType);
  if (pData != pNegative) {
    MFree_0D(pData);
  }
  addMapping(sclib::labelNegative, 0); //negative (-1) labeled data is mapped to model with index 0

  //build model for positive samples
  pData = scaleFeatures(pPositive, this->pScale, -1, true);
  this->pModels[1] = NULL;

	order = pModelHandler->guessModelOrder(pData, NULL, this->pTweak->classifierMl.modelType, this->pTweak->classifierMl.minOrder, this->pTweak->classifierMl.maxOrder);
  this->pModels[1] = pModelHandler->buildModel(pData, NULL, order, this->pTweak->classifierMl.modelType); 
  if (pData != pPositive) {
    MFree_0D(pData);
  }
  addMapping(sclib::labelPositive, 1); //positive (1) labeled data is mapped to model with index 1

  MFree_0D(pModelHandler);

  if (this->pModels[0] == NULL || this->pModels[1] == NULL) {
    freeModels();
    MFree_0D(this->pScale);
    sclib::destructLinkedList(this->pMapping);
    return SVLIB_Fail;
  } else {
    this->isTrained = true;
		this->classCount = 2;
    return SVLIB_Ok;
  }
}

//====================================================================================================================
//	train a classifier for distinguishing between several classes
//  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
//  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
//  respective row of pData.
//  the return value is meant to be the training error [0..1] or -1.0 on failure
//====================================================================================================================
int SC_Classifier_ML::trainMultiClass(SV_Data *pData, int *classes) {
  int count, label, i, x, y, *sortedClasses = NULL;
  unsigned short int order;
  bool success = true;
  SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak);
  SV_Data *pScaledData = NULL;
  SC_MatrixFunctions *pMat = new SC_MatrixFunctions();

  //destroy previously trained classifier (and corrsponding data), if any
  MFree_0D(this->pScale);
  freeModels();
  sclib::destructLinkedList(this->pMapping);
  this->isTrained = false;
	this->classCount = 0;

  //find scaling parameters
  if (this->doScaling == true) {
    this->pScale = findScalingParameters(pData);
  }

  //find number of different classes and establish a mapping between model-indices and class-labels
  this->modelCount = 0;
  sortedClasses = pMat->copy(classes, pData->Row);
  sclib::quickSort(sortedClasses, 0, pData->Row-1);
  for (i = 0; i < pData->Row-1; i++) {
    if (sortedClasses[i] != sortedClasses[i+1]) {
      addMapping(sortedClasses[i], this->modelCount);
      this->modelCount++;
    }
  }
  addMapping(sortedClasses[i], this->modelCount); //the last found class doesn't get added by the above loop
  this->modelCount++;
  MFree_1D(sortedClasses);

  //build models
  MArray_1D(this->pModels, this->modelCount, SC_Model*, "SC_Classifier_ML.trainMultiClass: pModels");
  for (i = 0; i < (int)(this->modelCount); i++) {
    label = getMappedLabel(i);

    //how much samples are in this class?
    count = 0;
    for (y = 0; y < pData->Row; y++) {
      if (classes[y] == label) {
        count++;
      }
    }
    pScaledData = new SV_Data(count, pData->Col);
    
    //copy samples together, perform scaling, if wished
    count = 0;
    for (y = 0; y < pData->Row; y++) {
      for (x = 0; x < pData->Col; x++) {
        if (classes[y] == label) {
          if (this->doScaling == true && this->pScale != NULL) { //do scaling here instead of with the explicit function to save time: copying has to be done here anyway
            pScaledData->Mat[count][x] = (pData->Mat[y][x] - this->pScale->Mat[0][x]) / (this->pScale->Mat[1][x] - this->pScale->Mat[0][x]); //(x-min)/(max-min)
          } else {
            pScaledData->Mat[count][x] = pData->Mat[y][x];
          }
        }
      }
    }

    //train model
    this->pModels[i] = NULL;
    order = pModelHandler->guessModelOrder(pScaledData, NULL, this->pTweak->classifierMl.modelType, this->pTweak->classifierMl.minOrder, this->pTweak->classifierMl.maxOrder);
    this->pModels[i] = pModelHandler->buildModel(pScaledData, NULL, order, this->pTweak->classifierMl.modelType);
    MFree_0D(pScaledData);
  
    if (this->pModels[i] == NULL) {
      success = false;
      break;
    }
  }

  MFree_0D(pModelHandler);
  MFree_0D(pMat);

  if (success == false) {
    freeModels();
    MFree_0D(this->pScale);
    sclib::destructLinkedList(this->pMapping);
    return SVLIB_Fail;
  } else {
    this->isTrained = true;
		this->classCount = this->modelCount;
    return SVLIB_Ok;
  }
}

//====================================================================================================================
//	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
//  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
//
//  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
//====================================================================================================================
int* SC_Classifier_ML::classify(SV_Data *pData, SV_Data* &pProbabilities) {
  double maxLikelihood, score;
  int *classes = NULL, bestIdx = -1;
  SV_Data *pScaledData = NULL;
  SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak, false);

  if (this->doScaling == true && this->pScale == NULL) {
    REPORT_ERROR(SVLIB_BadData, "Can't do scaling if no scaling parameters where found");
  }

  if (this->isTrained == true) {
    if (this->pModels != NULL && this->modelCount > 0) {
      
			MFree_0D(pProbabilities);
			pProbabilities = new SV_Data(pData->Row, this->modelCount);
      MArray_1D(classes, pData->Row, int, "SC_Classifier_ML.classify: classes");
      
      for (long int y = 0; y < pData->Row; y++) { //classify each feature vector separately
        pScaledData = scaleFeatures(pData, this->pScale, ((pData->Row > 1) ? y : -1), true); //save time for copying in case of just one row
        maxLikelihood = -1.0 * numeric_limits<double>::max();
        bestIdx = -1;

        for (int i = 0; i < (int)(this->modelCount); i++) {
          score = pModelHandler->testModel(this->pModels[i], pScaledData, 1, false);
          if (score > maxLikelihood) {
            maxLikelihood = score;
            bestIdx = i; //pick the model with the maximum likelihood value
          }
					pProbabilities->Mat[y][i] = (float)(score);
        }

        if (pScaledData != pData) {
          MFree_0D(pScaledData);
        }
        classes[y] = getMappedLabel(bestIdx);
      }

    } else {
      this->isTrained = false;
    }
  }

  MFree_0D(pModelHandler);
  return classes;
}

//====================================================================================================================
//	save a trained classifier to a file
//  write one summaryfile conatining the filenames of the single models
//====================================================================================================================
int SC_Classifier_ML::saveClassifier(const char *fileName) {
  int i, res = SVLIB_Fail;
  char *scaleFileName = NULL, *modelFileName = NULL, temp[sclib::bufferSize];
  SV_DataIO io;
	SC_Classifier_ML::SC_ClassMapping *pHook = this->pMapping;

  if (strlen(fileName) > 0) {
    if (this->isTrained == true) {
      if (this->pModels != NULL && this->modelCount > 0) {

        sprintf(temp, "%d\0", this->modelCount);
        sclib::stringOut(fileName, temp); //write the model-count in the summary-file

				//write class mapping to the summary-file
        sprintf(temp, "%d\0", sclib::getListCount(this->pMapping));
        sclib::stringOut(fileName, temp); //nr. of class-label-mappings to come
				while (pHook != NULL) {
					sprintf(temp, "%d %d\0", pHook->classLabel, pHook->modelIndex);
					sclib::stringOut(fileName, temp); //class-label/model-index pair
					pHook = pHook->Next;
				}

				for (i = 0; i < (int)(this->modelCount); i++) {
          sprintf(temp, ".model%d\0", i);
          modelFileName = sclib::exchangeFileExtension(fileName, temp);
	        this->pModels[i]->OpenFile(modelFileName, WRITE_MODEL);
	        res = this->pModels[i]->SaveModel(); //write each model in a separate file
	        this->pModels[i]->CloseFile();
          sprintf(temp, "%s %d\0", modelFileName, getMappedLabel(i)); //write the modelFileName and the label mapped to this model into the summary
          sclib::stringOut(fileName, temp);
          MFree_1D(modelFileName);
          if (res <= 0) {
            REPORT_ERROR(res, "Error during classifier saving");
          }
        }

        if (this->doScaling == true) { //save also the scaling parameters of this training set to use it with the test-data to come
          scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
          io.OpenFile(scaleFileName, WRITE_REC);
		      res = io.PutDataRec(*this->pScale);
          res = (res == 0) ? SVLIB_Fail : SVLIB_Ok; //translate between different traditions to report errrors or success...
          io.CloseFile();
          MFree_1D(scaleFileName);
        }
      } else {
        this->isTrained = false;
      }
    }
  }

  return res;
}

//====================================================================================================================
//	load a trained classifier from a file
//====================================================================================================================
int SC_Classifier_ML::loadClassifier(const char *fileName) {
  int modelIndex, classLabel, mappingLength, res = SVLIB_Ok, label, count = 0;
  char *modelFileName = NULL, *scaleFileName = NULL, buffer[sclib::bufferSize];
  SV_DataIO IO;
  SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak);
	FILE* inFile = NULL; 

  //destroy previously trained classifier, if any
  freeModels();
  sclib::destructLinkedList(this->pMapping);
  this->isTrained = false;

  //load models
  if ((inFile = fopen(fileName, "r")) != NULL) {
    sclib::readline(inFile, buffer, sclib::bufferSize); //read model-count
    this->modelCount = sclib::getNextIntFromString(buffer, sclib::bufferSize);
		this->classCount = this->modelCount;

		//read class mapping from the summary-file
    sclib::readline(inFile, buffer, sclib::bufferSize); //read model-count
    mappingLength = sclib::getNextIntFromString(buffer, sclib::bufferSize);
		for (int i = 0; i < mappingLength; i++) {
			if (sclib::readline(inFile, buffer, sclib::bufferSize) > 0) { //read data
				classLabel = sclib::getNextIntFromString(buffer, sclib::bufferSize);
				modelIndex = sclib::getNextIntFromString(buffer, sclib::bufferSize);
				addMapping(classLabel, modelIndex);
			} else {
				res = SVLIB_Fail;
				break;
			}
		}

    if (res !=  SVLIB_Fail && this->modelCount != SVLIB_Fail) {
      MArray_1D(this->pModels, this->modelCount, SC_Model*, "SC_Classifier_ML.loadClassifier: pModels");

	    while (!feof(inFile)) {
        sclib::readline(inFile, buffer, sclib::bufferSize);
        if (strlen(buffer) > 2) { //ignore blank lines
          modelFileName = sclib::getNextStringFromString(buffer, sclib::bufferSize);
          label = sclib::getNextIntFromString(buffer, sclib::bufferSize);

          if (sclib::fileExists(modelFileName) == true) {
            this->pModels[count] = pModelHandler->loadModel(modelFileName, this->pTweak->classifierMl.modelType);
            if (this->pModels[count] == NULL) {
              res = SVLIB_Fail;
              break;
            }
          } else {
            res = SVLIB_Fail;
            break;
          }
          if (addMapping(label, count) == false) {
            res = SVLIB_Fail;
            break;
          }

          MFree_1D(modelFileName);
          count++;
        }
      }
    } else {
      res = SVLIB_Fail;
    }

		fclose(inFile);
  }

  if (res == SVLIB_Ok) {
    this->isTrained = true;
  } else {
    freeModels();
    sclib::destructLinkedList(this->pMapping);
    this->isTrained = false;
  }

  //also try to load (if exists) the scaling parameters:
  MFree_0D(this->pScale);
  if (res == SVLIB_Ok && this->doScaling == true) {
    scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
    if (sclib::fileExists(scaleFileName) == true) {
      IO.OpenFile(scaleFileName, READ_REC);
      this->pScale = IO.GetAllRec();
      IO.CloseFile();
    }
    MFree_1D(scaleFileName);
  }

  MFree_0D(pModelHandler);

  return res;
}

//====================================================================================================================
//	establishes a mapping between class-labels and model-indices; it is not checked for collisions
//====================================================================================================================
bool SC_Classifier_ML::addMapping(int classLabel, int modelIndex) {
  SC_Classifier_ML::SC_ClassMapping *mapping = new SC_Classifier_ML::SC_ClassMapping(), *pHook;

  mapping->classLabel = classLabel;
  mapping->modelIndex = modelIndex;
  mapping->Next = NULL;
  
  if (this->pMapping == NULL) { //this is the first entry
    this->pMapping = mapping;
  } else {
    pHook = this->pMapping;
    while (pHook->Next != NULL) {
      pHook = pHook->Next;
    }
    pHook->Next = mapping;
  }
  
  return true;
}

//====================================================================================================================
//	return the first fitting label; return SVLIB_fail (-1, could also be a valid label... arghh) if no fitting label 
//  is found
//====================================================================================================================
int SC_Classifier_ML::getMappedLabel(int modelIndex) {
  SC_Classifier_ML::SC_ClassMapping *pHook = this->pMapping; 
  int classLabel = SVLIB_Fail;

  while (pHook != NULL) {
    if (pHook->modelIndex == modelIndex) {
      classLabel = pHook->classLabel;
      break;
    }
    pHook = pHook->Next;
  }

  return classLabel;
}

//====================================================================================================================
//	return the first fitting index; return SVLIB_fail (-1) if no fitting index is found
//====================================================================================================================
int SC_Classifier_ML::getMappedIndex(int classLabel) {
  SC_Classifier_ML::SC_ClassMapping *pHook = this->pMapping; 
  int modelIndex = SVLIB_Fail;

  while (pHook != NULL) {
    if (pHook->classLabel == classLabel) {
      modelIndex = pHook->modelIndex;
      break;
    }
    pHook = pHook->Next;
  }

  return modelIndex;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns -1 if something goes wrong
//====================================================================================================================
long int SC_Classifier_ML::label2idx(long int label) {
	long int res = -1;
	
	if (this->isTrained == true) {
		res = getMappedIndex(label);
	}

	return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns sclib::noType if something goes wrong
//====================================================================================================================
long int SC_Classifier_ML::idx2label(long int idx) {
	long int res = sclib::noType;
	
	if (this->isTrained == true) {
		if (idx >= 0 && idx < (long)(this->modelCount)) {
			res = getMappedLabel(idx);
		}
	}

	return res;
}

//====================================================================================================================
//	free the memory of the models
//====================================================================================================================
void SC_Classifier_ML::freeModels(void) {
  if (this->pModels != NULL) {
    for (unsigned int x = 0; x < this->modelCount; x++) {
      MFree_0D(this->pModels[x]);
    }
    MFree_1D(this->pModels);
  }

  return;
}
