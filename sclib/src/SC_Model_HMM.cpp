/**************************************************************************/
/*    A wrapper around SV_Lib's SV_Model_CHMM to implement a continous    */
/*    density Hidden Markov Model for kinds of layouts (i.e. ergodic/     */
/*    fully-connected, left-to-right/Bakis, you name it). The code in     */
/*    SV_Model_CHMM was slightly changed to incorporate the initial state */
/*    probability to allow for this flexibility. Now, the implementation  */
/*    closely follows Rabiner, "A Tutorial on Hidden Markov Models and    */
/*    Selected Applications in Speech Recognition", 1989. All additional  */
/*    sugar is implemented here: compliance with the SC_Model interface,  */
/*    including sampling from the model, and a useable interface to best- */
/*    path estimation via the Viterbi algorithm.                          */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 30.04.2008																								*/
/**************************************************************************/

#include <limits>
#include <time.h>
#include "SC_Model_HMM.h"
#include "SC_MatrixFunctions.h"
#include "SC_Aux.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_HMM::SC_Model_HMM(SC_TweakableParameters *pTweak, unsigned int stateCount, const char *transitionStructure, unsigned int mixturesPerState, bool useOrthogonalTransform, bool leftToRight, unsigned int maxIterations, bool verbose) : SC_Model(pTweak) {
	this->Hdr.ModelType = sclib::mtHMM;

	//all the parameters are actually stored in the encapsulated sv-lib chmm object
	this->pHMM = new SV_Model_CHMM();
	this->pHMM->Verbose = verbose;
	this->pHMM->Mixtures = mixturesPerState;
	this->pHMM->WithOrth = useOrthogonalTransform;
	this->pHMM->MaxIter = maxIterations;
	this->pHMM->RandSeed = (int)(time(NULL)); //change to a constant value to have deterministic behaviour for debugging purposes
	this->pHMM->NState = stateCount;
	this->pHMM->ConfMat = parseTransitionStructureString(transitionStructure, stateCount);
	this->pHMM->isLeftRightModel = leftToRight;
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model_HMM::SC_Model_HMM(SC_Model_HMM& pParent) : SC_Model(pParent) {
	int res;
	char fileName[sclib::bufferSize];
	SC_MatrixFunctions matFunc;

	//create a copy by saving/loading to/from a temp file to copy the private members of the svlib-hmm
	sprintf(fileName, "%s%s%s", this->pTweak->debug.debugDir, tmpnam(NULL)+1, "tmp"); //+1 to overread the heading "\"
	pParent.OpenFile(fileName, WRITE_MODEL);
	res = pParent.SaveModel();
	pParent.CloseFile();
	if (res > 0) {
		this->OpenFile(fileName, READ_MODEL);
		this->LoadModel();
		this->CloseFile();
		if (this->pHMM == NULL) {
			REPORT_ERROR(SVLIB_FileErr, "couldn't load temporary hmm model");
		}
	} else {
		REPORT_ERROR(SVLIB_FileErr, "couldn't save temporary hmm model");
	}
	remove(fileName);

	//the public parameters wheren't saved in the SV_Model_CHMM class, so fill it now
	this->pHMM->Verbose = pParent.pHMM->Verbose;
	this->pHMM->Mixtures = pParent.pHMM->Mixtures;
	this->pHMM->WithOrth = pParent.pHMM->WithOrth;
	this->pHMM->MaxIter = pParent.pHMM->MaxIter;
	this->pHMM->RandSeed = pParent.pHMM->RandSeed; //TODO: ok?
	this->pHMM->NState =  pParent.pHMM->NState;
	this->pHMM->ConfMat = (pParent.pHMM->ConfMat != NULL) ? matFunc.copy(pParent.pHMM->ConfMat, this->pHMM->NState, this->pHMM->NState) : NULL;
	this->pHMM->isLeftRightModel = pParent.pHMM->isLeftRightModel;
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_HMM::~SC_Model_HMM() {
	MFree_1D(this->pHMM->ConfMat);
	MFree_0D(this->pHMM);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_HMM& SC_Model_HMM::operator=(SC_Model_HMM& pParent) {
	int res;
	char fileName[sclib::bufferSize];
	SC_MatrixFunctions matFunc;

	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		//create a copy by saving/loading to/from a temp file to copy the private members of the svlib-hmm
		sprintf(fileName, "%s%s%s", this->pTweak->debug.debugDir, tmpnam(NULL)+1, "tmp"); //+1 to overread the heading "\"
		pParent.OpenFile(fileName, WRITE_MODEL);
		res = pParent.SaveModel();
		pParent.CloseFile();
		if (res > 0) {
			this->OpenFile(fileName, READ_MODEL);
			this->LoadModel();
			this->CloseFile();
			if (this->pHMM == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "couldn't load temporary hmm model");
			}
		} else {
			REPORT_ERROR(SVLIB_FileErr, "couldn't save temporary hmm model");
		}
		remove(fileName);

		//the public parameters wheren't saved in the SV_Model_CHMM class, so fill it now
		this->pHMM->Verbose = pParent.pHMM->Verbose;
		this->pHMM->Mixtures = pParent.pHMM->Mixtures;
		this->pHMM->WithOrth = pParent.pHMM->WithOrth;
		this->pHMM->MaxIter = pParent.pHMM->MaxIter;
		this->pHMM->RandSeed = pParent.pHMM->RandSeed; //TODO: ok?
		this->pHMM->NState =  pParent.pHMM->NState;
		this->pHMM->ConfMat = (pParent.pHMM->ConfMat != NULL) ? matFunc.copy(pParent.pHMM->ConfMat, this->pHMM->NState, this->pHMM->NState) : NULL;
		this->pHMM->isLeftRightModel = pParent.pHMM->isLeftRightModel;
	}

	return *this;
}

//====================================================================================================================
// Parses a string that represents a boolean matrix (columns delimited by a space ' ', rows delimited by a semicolon
// ';', true by '1', false by '0'). The entries in row i, column j in this matrix mean if the hmm-state i should have
// a link (state transition probability possibly > 0) to state j or not. Returned is the reconstructed matrix with
// integer datatype as needed by the SV_Model_CHMM class for the ConfMat member. if the string is invalid or "" or 
// NULL, a fully-connected model is assumed or too short or missing rows are filled up with '1'.
//====================================================================================================================
int** SC_Model_HMM::parseTransitionStructureString(const char *transitionStructure, unsigned int stateCount) {
	SC_MatrixFunctions matFunc;
	int **confMat = matFunc.initMatrix(stateCount, stateCount, (int)(1));
	unsigned int i = 0, j = 0;
	char *c = const_cast<char*>(transitionStructure);

	//parse the string character by character and set elements in the matrix to zero 
	//attention: this is not robust against departures from the specified format: exactly one ' ' for column-change, only one ';' for row-change!
	while (c != NULL && *c != '\0') {
		switch (*c) {
			case '0': //a zero is found, mark this in the confMat
				if (i < stateCount && j < stateCount) {
					confMat[i][j] = 0;
				} else {
					REPORT_ERROR(SVLIB_BadArg, "Given transition structure string implies a bigger matrix than actually allowed");
				}
				break;
			case ' ': //increase column-count
				j++;
				break;
			case ';': //increase row-count
				i++;
				j = 0;
				break;
			default:
				//do nothing
				break;
		}
		c++;
	}

	return confMat;
}

//====================================================================================================================
// The other way round
//====================================================================================================================
char* SC_Model_HMM::constructTransitionStructureString(int **transitionMatrix, unsigned int stateCount) {
	char *transitionString = new char[stateCount*2*stateCount +1];
	int c = 0;
	
	for (unsigned int i = 0; i < stateCount; i++) {
		for (unsigned int j = 0; j < stateCount; j++) {
			transitionString[c++] = (transitionMatrix[i][j] != 0) ? '1' : '0';
			transitionString[c++] = (j < stateCount-1) ? ' ' : ';';
		}
	}
	transitionString[c] = '\0';

	return transitionString;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_HMM::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_HMM*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_HMM*)pBackgroundModels));
}
SC_Model_HMM* SC_Model_HMM::combineModels(SC_Model_HMM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_HMM* pBackgroundModels) {
	char *transitionStructure = constructTransitionStructureString(this->pHMM->ConfMat, this->pHMM->NState);
	SC_Model_HMM *pNewModel = new SC_Model_HMM(this->pTweak, this->pHMM->NState, transitionStructure, this->pHMM->Mixtures, (this->pHMM->WithOrth!=0)?true:false, this->pHMM->isLeftRightModel, this->pHMM->MaxIter,(this->pHMM->Verbose!=0)?true:false);

	MFree_1D(transitionStructure);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  return pNewModel;
}

//====================================================================================================================
// Train a model; here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_HMM::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Train a model; here, one can specify how many segments of a linked list of feature-vectors should be used
// (0 means all)
//====================================================================================================================
int	SC_Model_HMM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	int res, featureCnt = 0;
  SV_Data *pHook, *pSave, *pTmp;

	//we don't need to copy things together here (the svlib-chmm is able to learn from linked lists), just in case cut the list
	if (segmentsToMerge > 0) {
		pHook = sclib::getListWithIndex(pData, segmentsToMerge-1);
		pSave = pHook->Next;
		pHook->Next = NULL;
	} else {
		pHook = pData;
		pSave = pData->Next;
	}

	res = this->pHMM->TrainModel(pData);
	
	//get training datat count
	pTmp = pData;
	while (pTmp != NULL) {
		featureCnt += pTmp->Row;
		pTmp = pTmp->Next;
	}
	this->trainingDataCount = featureCnt;

	pHook->Next = pSave; //restore a possibly cutted list

  return res;
}

//====================================================================================================================
// Test the model while merging no segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_HMM::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;
}

//====================================================================================================================
// Test the model with specified nr. of segments in the linked list of feature-vectors.
//====================================================================================================================
SV_Data* SC_Model_HMM::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
	double sumScores = 0.0, sumFeatures = 0.0;
  SV_Data *pHook, *pSave, *pRes = NULL, *pTemp;

	if (this->pHMM->GetStaNum() > 0) { //the model is trained
		//we must not to copy things together here, just in case cut the list; this is becaue a hmm models a temporal structure,
		//and copying together would change it; instead, we score each list item individually, and return a later on combined
		//score along with the individual segment's scores following form row-id 1.
		if (segmentsToMerge > 0) {
			pHook = sclib::getListWithIndex(TestData, segmentsToMerge-1);
			pSave = pHook->Next;
			pHook->Next = NULL;
		} else {
			pHook = TestData;
			pSave = TestData->Next;
		}

		pTemp = this->pHMM->TestModel(TestData);
		pHook->Next = pSave; //restore a possibly cutted list

		//create a new zeroth row in the resultset that contains the combined and properly weighted scores
		pHook = TestData;
		if (pTemp != NULL) {
			pRes = new SV_Data(pTemp->Row+1, 1);
			for (int i = 0; i < pTemp->Row; i++) {
				pRes->Mat[i+1][0] = pTemp->Mat[i][0];
				sumScores += pTemp->Mat[i][0] * (double)(pHook->Row);
				sumFeatures += (double)(pHook->Row);
				pHook = pHook->Next;
			}
			pRes->Mat[0][0] = (float)(sumScores / sumFeatures);
		}
	}

	return pRes;
}

//====================================================================================================================
// load a model from file
//====================================================================================================================
SV_Model*	SC_Model_HMM::LoadModel(void) {
	int res;
	bool modelFileExists;
	char *buffer;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_HMM Failed!");}

	//read model
	io.readScalar(&(this->DFile), modelFileExists, codeSizes, fileSizes); //read if a model file is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_HMM Model Failed!");}
	if (modelFileExists == true) {
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".hmm_model");
		this->pHMM->LoadModel();
		MFree_1D(buffer);
	}

	return(this);
}

//====================================================================================================================
// save a model to file
//====================================================================================================================
int	SC_Model_HMM::SaveModel(void) {
	int bytes;
	char *buffer;
	SV_DataIO io;

  //write trainingDataCount
	bytes = io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_HMM Model Failed!");}

	//write HMM
	bytes += io.writeScalar(&(this->DFile), this->pHMM!=NULL); //store if a model file is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_HMM Model Failed!");}
	if (this->pHMM != NULL) {
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".hmm_model");
		bytes += this->pHMM->SaveModel();
		MFree_1D(buffer);
	}

  return bytes;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_HMM::modelOut(ostream& os) {
	os << this->pHMM;

	return(os);
}

//====================================================================================================================
// for computing BIC etc.
//====================================================================================================================
unsigned int SC_Model_HMM::getFreeParameterCount(void) {
	unsigned int cnt = this->pHMM->GetStaNum() + this->pHMM->GetStaNum()*this->pHMM->GetStaNum() + this->pHMM->GetStaNum()*this->pHMM->GetMixNum() + this->pHMM->GetStaNum()*this->pHMM->GetMixNum()*this->pHMM->GetVecDim()*2;
	//                 initial state probs.    + transition probs.                               + mixture weights                                 + mixture-means and -variances
	return cnt;
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_HMM::drawSamplesFromDistribution(unsigned long int count) {
	unsigned long int t;
	unsigned int state, lastState, mixture, d;
  SV_Data *pSamples = NULL;
	double **mixtureWeight;

	if (this->pHMM->GetStaNum() > 0) { //the model is trained
    pSamples = new SV_Data(count, this->pHMM->GetVecDim());

		MArray_2D(mixtureWeight, this->pHMM->GetStaNum(), this->pHMM->GetMixNum(), double, "SC_Model_HMM.drawSamplesFromDistribution: mixtureWeights");
		for (state = 0; state < (unsigned int)(this->pHMM->GetStaNum()); state++) {
			for (mixture = 0; mixture < (unsigned int)(this->pHMM->GetMixNum()); mixture++) {
				mixtureWeight[state][mixture] = this->pHMM->GetEmit()[state][mixture].wgt; //create a copy that can directly be passed to sclib::drawIndexFromDistribution()
			}
		}

	  for (t = 0; t < count; t++) {
			//choose a state
			if (t == 0) {
				state = sclib::drawIndexFromDistribution(this->pHMM->getInitialStateProbability(), this->pHMM->GetStaNum());
			} else {
				state = sclib::drawIndexFromDistribution(this->pHMM->GetTran()[lastState], this->pHMM->GetStaNum());
			}

			//choose a mixture in the state's gmm
			mixture = sclib::drawIndexFromDistribution(mixtureWeight[state], this->pHMM->GetMixNum());

			//draw a sample from the choosen gmm mixture
      for (d = 0; d < (unsigned int)(this->pHMM->GetVecDim()); d++) {
				pSamples->Mat[t][d] = (float)(sclib::getRandomizer()->rand_gaus(this->pHMM->GetEmit()[state][mixture].mean[d], sqrt(this->pHMM->GetEmit()[state][mixture].vari[d]))); //this uses the Box-Muller algorithm inside as in the paper
      }

			lastState = state;
		}

		MFree_2D(mixtureWeight);
	}

  return pSamples;
}

//====================================================================================================================
// uses the viterbi-algorithm to find the best state sequence through the previously trained model. an int-array with 
// rows corresponding to the pFeatures-rows is returned, where each entry gives a state-number between 0 and 
// stateCount-1
//====================================================================================================================
int* SC_Model_HMM::getBestStateSequence(SV_Data *pFeatures, double &logLikelihood) {
	int *stateSequence = NULL;
	double **mixtureOutputs;

	logLikelihood = -1.0 * std::numeric_limits<double>::max();
	if (this->pHMM->GetStaNum() > 0) { //the model is trained
		MArray_1D(stateSequence, pFeatures->Row, int, "SC_Model_HMM.getBestStateSequence: stateSequence");
		MArray_2D(mixtureOutputs, this->pHMM->GetStaNum(), pFeatures->Row, double, "SC_Model_HMM.getBestStateSequence: mixtureOutputs");

		this->pHMM->getMixtureEmissions(mixtureOutputs, pFeatures->Mat, pFeatures->Row);
		logLikelihood = this->pHMM->doViterbi(mixtureOutputs, stateSequence, pFeatures->Row);

		MFree_2D(mixtureOutputs);
	}

	return stateSequence;
}
