/**************************************************************************/
/*    A wrapper around SV_Lib's SV_Model_Gaus to implement a full-        */
/*    covariance multivariate gaussian model (no mixture-model!).         */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.05.2008																								*/
/**************************************************************************/

#include "SC_Model_FullGauss.h"
#include "SC_Aux.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_FullGauss::SC_Model_FullGauss(SC_TweakableParameters *pTweak) : SC_Model(pTweak) {
	this->Hdr.ModelType = sclib::mtFullGauss;
	this->pGauss = new SV_Model_Gaus();
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model_FullGauss::SC_Model_FullGauss(SC_Model_FullGauss& pParent) : SC_Model(pParent) {
	int res;
	char fileName[sclib::bufferSize];

	this->pGauss = new SV_Model_Gaus();

	//create a copy by saving/loading to/from a temp file to copy the private members of the svlib-model
	sprintf(fileName, "%s%s%s", this->pTweak->debug.debugDir, tmpnam(NULL)+1, "tmp"); //+1 to overread the heading "\"
	pParent.pGauss->OpenFile(fileName, WRITE_MODEL);
	res = pParent.pGauss->SaveModel();
	pParent.pGauss->CloseFile();
	if (res > 0) {
		this->pGauss->OpenFile(fileName, READ_MODEL);
		this->pGauss->LoadModel();
		this->pGauss->CloseFile();
		if (this->pGauss == NULL) {
			REPORT_ERROR(SVLIB_FileErr, "couldn't load temporary full-gauss model");
		}
	} else {
		REPORT_ERROR(SVLIB_FileErr, "couldn't save temporary full-gauss model");
	}
	remove(fileName);
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_FullGauss::~SC_Model_FullGauss() {
	MFree_0D(this->pGauss);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_FullGauss& SC_Model_FullGauss::operator=(SC_Model_FullGauss& pParent) {
	int res;
	char fileName[sclib::bufferSize];

	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		//create a copy by saving/loading to/from a temp file to copy the private members of the svlib-model
		sprintf(fileName, "%s%s%s", this->pTweak->debug.debugDir, tmpnam(NULL)+1, "tmp"); //+1 to overread the heading "\"
		pParent.pGauss->OpenFile(fileName, WRITE_MODEL);
		res = pParent.pGauss->SaveModel();
		pParent.pGauss->CloseFile();
		if (res > 0) {
			this->pGauss->OpenFile(fileName, READ_MODEL);
			this->pGauss->LoadModel();
			this->pGauss->CloseFile();
			if (this->pGauss == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "couldn't load temporary full-gauss model");
			}
		} else {
			REPORT_ERROR(SVLIB_FileErr, "couldn't save temporary full-gauss model");
		}
		remove(fileName);
	}

  return *this;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_FullGauss::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_FullGauss*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_FullGauss*)pBackgroundModels));
}
SC_Model_FullGauss* SC_Model_FullGauss::combineModels(SC_Model_FullGauss* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_FullGauss* pBackgroundModels) {
	SC_Model_FullGauss *pNewModel = new SC_Model_FullGauss(this->pTweak);

	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  return pNewModel;
}

//====================================================================================================================
// Train a model; here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_FullGauss::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Train a model; here, one can specify how many segments of a linked list of feature-vectors should be used
// (0 means all)
//====================================================================================================================
int	SC_Model_FullGauss::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	int res = this->pGauss->TrainModel(pData);

	this->trainingDataCount = pCompleteData->Row;

	if (pCompleteData != pData) {
		MFree_0D(pCompleteData);
	}

  return res;
}

//====================================================================================================================
// Test the model while merging no segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_FullGauss::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;
}

//====================================================================================================================
// Test the model with specified nr. of segments in the linked list of feature-vectors.
//====================================================================================================================
SV_Data* SC_Model_FullGauss::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge);
	SV_Data *pScore = this->pGauss->TestModel(pCompleteData);

	pScore->Mat[0][0] /= (float)(pCompleteData->Row); //log-likelihood => average log-likelihood as proposed in "Speaker Verification Using Adapted Gaussian Mixture Models"

	if (pCompleteData != TestData) {
		MFree_0D(pCompleteData);
	}

	return pScore;
}

//====================================================================================================================
// load a model from file
//====================================================================================================================
SV_Model*	SC_Model_FullGauss::LoadModel(void) {
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
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_FullGauss Failed!");}

	//read model
	io.readScalar(&(this->DFile), modelFileExists, codeSizes, fileSizes); //read if a model file is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_FullGauss Model Failed!");}
	if (modelFileExists == true) {
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".gauss_model");
		this->pGauss->LoadModel();
		MFree_1D(buffer);
	}

	return(this);
}

//====================================================================================================================
// save a model to file
//====================================================================================================================
int	SC_Model_FullGauss::SaveModel(void) {
	int bytes;
	char *buffer;
	SV_DataIO io;

  //write trainingDataCount
	bytes = io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_FullGauss Model Failed!");}

	//write model
	bytes += io.writeScalar(&(this->DFile), this->pGauss!=NULL); //store if a model file is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_FullGauss Model Failed!");}
	if (this->pGauss != NULL) {
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".gauss_model");
		bytes += this->pGauss->SaveModel();
		MFree_1D(buffer);
	}

  return bytes;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_FullGauss::modelOut(ostream& os) {
	os << *this->pGauss;

	return(os);
}

//====================================================================================================================
// for computing BIC etc.
//====================================================================================================================
unsigned int SC_Model_FullGauss::getFreeParameterCount(void) {
	unsigned int cnt = this->pGauss->getDim()*this->pGauss->getDim() + this->pGauss->getDim();
	//                 covariance                                    + mean
	return cnt;
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_FullGauss::drawSamplesFromDistribution(unsigned long int count) {
	int t, d, res, D = this->pGauss->getDim(); 
	double **randoms = NULL, **covar;
	SV_Data *pSamples = NULL;
	SC_MatrixFunctions matFunc;

	covar = matFunc.copy(this->pGauss->getInvertedCov(), D, D);
	res = matFunc.inv(covar, D);

	if (res != SVLIB_Fail) {
		randoms = sclib::randN(count, D, this->pGauss->getMean(), covar);
	}
	MFree_2D(covar);

	if (randoms != NULL) {
		pSamples = new SV_Data(count, D);
		for (t = 0; t < pSamples->Row; t++) {
			for (d = 0; d < pSamples->Col; d++) {
				pSamples->Mat[t][d] = (float)(randoms[t][d]);
			}
		}
		MFree_2D(randoms);
	}

  return pSamples;
}
