/**************************************************************************/
/*    A model that doesn't solve the (hard) problem of estimating the pdf */
/*    of the data but instead the simpler problem of estimating the level */
/*    set (support or center of mass of the pdf) via a one-class nu-SVM   */
/*    as described in                                                     */
/*      "An Online Kernel Change Detection Algorithm", Desobry, Davy,     */
/*      Doncarli, IEEE Trans. Sig. Proc., Vol 53, No. 8, 2005,            */
/*      pp. 2961-2974                                                     */
/*    and                                                                 */
/*      "Unsupervised Speaker Indexing Using One-Class Support Vector     */
/*      Machines", Fergani, Davy, Houacine, 2006                          */
/*                                                                        */
/*    The model is typically evaluated via computing a distance (proposed */
/*    in the same papers) to a second model of same type; such, when      */
/*    evaluated giving only a dataset, this is first modeled.             */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 17.04.2007																								*/
/**************************************************************************/

#include "SC_Model_SVM.h"
#include "SC_DistanceMeasures.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_SVM::SC_Model_SVM(SC_TweakableParameters *pTweak, bool distanceBasedTesting, bool doParameterSearch, bool verbose) : SC_Model(pTweak) {
	this->verbose = verbose;
	this->doParameterSearch = doParameterSearch;
	this->distanceBasedTesting = distanceBasedTesting;
	this->pSVM = new SC_Classifier_SVM(this->pTweak, this->doParameterSearch, false, this->verbose);
	this->Hdr.ModelType = sclib::mtSVM;
  this->dim = 0;
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model_SVM::SC_Model_SVM(const SC_Model_SVM& pParent, bool justLink) : SC_Model(pParent) {
	this->doParameterSearch = pParent.doParameterSearch;
	this->distanceBasedTesting = pParent.distanceBasedTesting;
  this->pSVM = new SC_Classifier_SVM(*pParent.pSVM, justLink);
  this->dim = pParent.dim;
	this->verbose = pParent.verbose;
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_SVM::~SC_Model_SVM() {
	MFree_0D(this->pSVM);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_SVM& SC_Model_SVM::operator=(const SC_Model_SVM& pParent) {
	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		//Destruct old parameters
		if (this->pSVM != NULL) {
			MFree_0D(this->pSVM);
		}	

		this->verbose = pParent.verbose;
		this->dim = pParent.dim;
		this->doParameterSearch = pParent.doParameterSearch;
		this->distanceBasedTesting = pParent.distanceBasedTesting;
		this->pSVM = new SC_Classifier_SVM(*pParent.pSVM);
	}

  return *this;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_SVM::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_SVM*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_SVM*)pBackgroundModels));
}
SC_Model_SVM* SC_Model_SVM::combineModels(SC_Model_SVM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_SVM* pBackgroundModels) {
  SC_Model_SVM *pNewModel = new SC_Model_SVM(this->pTweak);

	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  return pNewModel;
}

//====================================================================================================================
// Train a model; here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_SVM::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Train a model; here, one can specify how many segments of a linked list of feature-vectors should be merged
// (0 means all)
//====================================================================================================================
int	SC_Model_SVM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	int res = SVLIB_Ok;
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);

  //create new one
  this->trainingDataCount = pCompleteData->Row;
  this->dim = pCompleteData->Col;

  //train one-class svm
	this->pSVM->setDoParameterSearch(this->doParameterSearch);
	res = this->pSVM->trainOneClass(pCompleteData);
	if (res == SVLIB_Fail) {
		this->trainingDataCount = 0;
		this->dim = 0;
	}

	//clean up
  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  return res;
}

//====================================================================================================================
// Train a model using SVM-Problem data structure directly as input
//====================================================================================================================
int	SC_Model_SVM::TrainModel(SC_SVMproblem *pProblem) {
	int res = SVLIB_Ok, i, j;

  //create new one
  this->trainingDataCount = pProblem->l;
	this->dim = 0; //find max. dimensionality of data in sparse matrix representation of the svm-problem data structure
	for (i = 0; i < pProblem->l; i++) {
		j = 0;
		while (pProblem->x[i][j].index != -1) {
			j++;
		}
		if (j > 0) {
			if (pProblem->x[i][j-1].index > (int)(this->dim)) {
				this->dim = pProblem->x[i][j-1].index;
			}
		}
	}

  //train one-class svm
	this->pSVM->setDoParameterSearch(this->doParameterSearch);
	res = this->pSVM->trainOneClass(pProblem);
	if (res == SVLIB_Fail) {
		this->trainingDataCount = 0;
		this->dim = 0;
	}

  return res;
}

//====================================================================================================================
// Test the model while merging no segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_SVM::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;
}

//====================================================================================================================
// Test the model with specified nr. of segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_SVM::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge); //merge all feature-vectors in a linked list together
  SV_Data *pScore = NULL, *pProbabilities = NULL;
	SC_Model_SVM *pSecondModel = NULL;
	int *labels;

	if (this->dim > 0) {
		if (this->distanceBasedTesting == true) {
			pSecondModel = new SC_Model_SVM(this->pTweak);
			pSecondModel->TrainModel(TestData, segmentsToMerge);
			pScore = TestModel(pSecondModel); //TODO: make score (more) like an average log-likelihood to be usabel in BIC etc.
			MFree_0D(pSecondModel);
		} else {
			int oldP = this->pTweak->classifierSvm.probability;
			this->pTweak->classifierSvm.probability = 0;
			labels = this->pSVM->classify(pCompleteData, pProbabilities);
			this->pTweak->classifierSvm.probability = oldP;
			MFree_0D(pProbabilities);
			if (labels != NULL) {
				pScore = new SV_Data(pCompleteData->Row+1, 1); //0th row contains overall result, following rows contain classification label of corresponding test data row+1
				pScore->Mat[0][0] = 0.0f;
				for (int t = 0; t < pCompleteData->Row; t++) {
					pScore->Mat[0][0] += (labels[t]==sclib::labelPositive) ? 1.0f : 0.0f;
					pScore->Mat[t+1][0] = (float)(labels[t]);
				}
				pScore->Mat[0][0] = (float)(sclib::sLog(pScore->Mat[0][0] / (double)(pCompleteData->Row)));
			} else {
				pScore = new SV_Data(1, 1); 
				pScore->Mat[0][0] = (float)(sclib::sLog(0.0));
			}
			MFree_1D(labels);
		}
	} else {
		REPORT_ERROR(SVLIB_Fail, "Can't test on an untrained model");
	}
  
  if (TestData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  return pScore;
}

//====================================================================================================================
// scoring is always accomplished here by giving a distance between this model and a second model (possibly temporary 
// created from a TestData dataset)
//====================================================================================================================
SV_Data* SC_Model_SVM::TestModel(SC_Model_SVM *pModel) {
	SV_Data *pScore = NULL;
	SC_DistanceMeasures *pDist = NULL;

	if (this->dim > 0) {
		pScore = new SV_Data(1, 1);
		pDist = new SC_DistanceMeasures(this->pTweak);
		pScore->Mat[0][0] = (float)(pDist->svmArcDistance(this->pSVM->getSVM(), this->pSVM->getClassifier(), pModel->getSVMmodel()));
		MFree_0D(pDist);
	} else {
    REPORT_ERROR(SVLIB_Fail, "Can't test on an untrained model");
	}

	return pScore; //TODO: make score (more) like an average log-likelihood to be usabel in BIC etc.
}

//====================================================================================================================
// load a model from file
//====================================================================================================================
SV_Model*	SC_Model_SVM::LoadModel(void) {
	int res, newDim = 0;
	char *buffer;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

	//read scalar parameters
	io.readScalar(&(this->DFile), this->dim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {return(NULL);}
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_SVM Failed!");}
	io.readScalar(&(this->DFile), this->verbose, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	io.readScalar(&(this->DFile), this->doParameterSearch, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	io.readScalar(&(this->DFile), this->distanceBasedTesting, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}

	//read SVM via methods by SC_Classifier_SVM
	buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".svm_model");
	res = this->pSVM->loadClassifier(buffer);
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_SVM model Failed!");}
	MFree_1D(buffer);

	return(this);
}

//====================================================================================================================
// save a model to file
//====================================================================================================================
int	SC_Model_SVM::SaveModel(void) {
 	int res, bytes;
	SV_DataIO io;
	char *buffer;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}

  //write scalar parameters
	bytes = io.writeScalar(&(this->DFile), this->dim);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->verbose);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->doParameterSearch);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->distanceBasedTesting);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}

	//write SVM via methods provided by SC_Classifier_SVM
	buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".svm_model");
	res = this->pSVM->saveClassifier(buffer);
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_SVM Failed!");}
	MFree_1D(buffer);

  return bytes + MHLen; //MHLen may be incorrect...
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_SVM::modelOut(ostream& os) {
	os.setf(ios::fixed|ios::basefield);
  os.precision(5);

  if (this->dim > 0) { 
		os << "doParameterSearch: " << this->doParameterSearch << endl;
		os << "distanceBasedScoring: " << this->distanceBasedTesting << endl;
		os << "Dimensions: " << this->dim << endl << endl;
		//TODO: put out SVM's SVs and alphas...
  } else {
    os << "The model is still untrained..." << endl;
  }

	os << endl;

	return(os);
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_SVM::drawSamplesFromDistribution(unsigned long int count) {
	SV_Data *pSamples = NULL;

	if (this->pSVM != NULL) {
		SV_Data *pScalingParameters = this->pSVM->getScalingParameters();
		pSamples = new SV_Data(count, this->dim);

		for (unsigned long int t = 0; t < count; t++) {
			int randomIdx = sclib::rand(this->pSVM->getClassifier()->l - 1); //select a support vector at random
			int d = 0, dd = 0;
			while (this->pSVM->getClassifier()->SV[randomIdx][d].index > -1) {
				for (int c = dd; c < this->pSVM->getClassifier()->SV[randomIdx][d].index-1; c++) { //svmnodes are sparse; fill all the columns in the svdata-representation with zeros that do not occur in the sparse representation of the svmnode
					pSamples->Mat[t][c] = 0.0f;
					dd++;
				}
				pSamples->Mat[t][dd] = (float)(this->pSVM->getClassifier()->SV[randomIdx][d].value);
				if (this->pSVM->doesScaling() == true) { //remove scaling, i.e. undo (x-min)/(max-min)// TODO: include this in SC_Classifier for sake of encapsulation
					pSamples->Mat[t][dd] = pSamples->Mat[t][dd] * (pScalingParameters->Mat[1][dd]-pScalingParameters->Mat[0][dd]) + pScalingParameters->Mat[0][dd];
				}
				dd++;
				d++;
			}
		}
	}

	return pSamples;
}
