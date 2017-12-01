/**************************************************************************/
/*    A wrapper around SV_Lib's SV_Model_VQ to implement a vector-        */
/*    quatisation model trained by LBG algorithm.                         */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 01.08.2008																								*/
/**************************************************************************/

#include "SC_Model_VQ.h"
#include "SC_Aux.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_VQ::SC_Model_VQ(SC_TweakableParameters *pTweak, unsigned int codebookSize, unsigned int splitMethod, unsigned int maxIterations, bool verbose) : SC_Model(pTweak) {
	this->Hdr.ModelType = sclib::mtVQ;
	this->pVQ = new SV_Model_VQ();

	this->pVQ->Verbose = 0; //(verbose == true) ? 1 : 0;
	this->pVQ->MaxIter = maxIterations;
	this->pVQ->SplitMethod = (splitMethod = sclib::modeLBG) ? 0 : 1;
	this->pVQ->CBSize = codebookSize;
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model_VQ::SC_Model_VQ(SC_Model_VQ& pParent) : SC_Model(pParent) {
	int res;
	char fileName[sclib::bufferSize];

	this->pVQ = new SV_Model_VQ();

	//create a copy by saving/loading to/from a temp file to copy the private members of the svlib-model
	sprintf(fileName, "%s%s%s", this->pTweak->debug.debugDir, tmpnam(NULL)+1, "tmp"); //+1 to overread the heading "\"
	pParent.pVQ->OpenFile(fileName, WRITE_MODEL);
	res = pParent.pVQ->SaveModel();
	pParent.pVQ->CloseFile();
	if (res > 0) {
		this->pVQ->OpenFile(fileName, READ_MODEL);
		this->pVQ->LoadModel();
		this->pVQ->CloseFile();
		if (this->pVQ == NULL) {
			REPORT_ERROR(SVLIB_FileErr, "couldn't load temporary vq model");
		}
	} else {
		REPORT_ERROR(SVLIB_FileErr, "couldn't save temporary vq model");
	}
	remove(fileName);
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_VQ::~SC_Model_VQ() {
	MFree_0D(this->pVQ);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_VQ& SC_Model_VQ::operator=(SC_Model_VQ& pParent) {
	int res;
	char fileName[sclib::bufferSize];

	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		//create a copy by saving/loading to/from a temp file to copy the private members of the svlib-model
		sprintf(fileName, "%s%s%s", this->pTweak->debug.debugDir, tmpnam(NULL)+1, "tmp"); //+1 to overread the heading "\"
		pParent.pVQ->OpenFile(fileName, WRITE_MODEL);
		res = pParent.pVQ->SaveModel();
		pParent.pVQ->CloseFile();
		if (res > 0) {
			this->pVQ->OpenFile(fileName, READ_MODEL);
			this->pVQ->LoadModel();
			this->pVQ->CloseFile();
			if (this->pVQ == NULL) {
				REPORT_ERROR(SVLIB_FileErr, "couldn't load temporary vq model");
			}
		} else {
			REPORT_ERROR(SVLIB_FileErr, "couldn't save temporary vq model");
		}
		remove(fileName);
	}

  return *this;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_VQ::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_VQ*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_VQ*)pBackgroundModels));
}
SC_Model_VQ* SC_Model_VQ::combineModels(SC_Model_VQ* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_VQ* pBackgroundModels) {
	SC_Model_VQ *pNewModel = new SC_Model_VQ(this->pTweak);

	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  return pNewModel;
}

//====================================================================================================================
// Train a model; here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_VQ::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Train a model; here, one can specify how many segments of a linked list of feature-vectors should be used
// (0 means all)
//====================================================================================================================
int	SC_Model_VQ::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	int res = this->pVQ->TrainModel(pData);

	this->trainingDataCount = pCompleteData->Row;

	if (pCompleteData != pData) {
		MFree_0D(pCompleteData);
	}

  return res;
}

//====================================================================================================================
// Test the model while merging no segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_VQ::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;
}

//====================================================================================================================
// Test the model with specified nr. of segments in the linked list of feature-vectors.
//====================================================================================================================
SV_Data* SC_Model_VQ::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge);
	SV_Data *pScore = this->pVQ->TestModel(pCompleteData);

	pScore->Mat[0][0] /= (float)(pCompleteData->Row); //log-likelihood => average log-likelihood as proposed in "Speaker Verification Using Adapted Gaussian Mixture Models"

	if (pCompleteData != TestData) {
		MFree_0D(pCompleteData);
	}

	return pScore;
}

//====================================================================================================================
// load a model from file
//====================================================================================================================
SV_Model*	SC_Model_VQ::LoadModel(void) {
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
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_VQ Failed!");}

	//read model
	io.readScalar(&(this->DFile), modelFileExists, codeSizes, fileSizes); //read if a model file is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_VQ Model Failed!");}
	if (modelFileExists == true) {
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".vq_model");
		this->pVQ->LoadModel();
		MFree_1D(buffer);
	}

	return(this);
}

//====================================================================================================================
// save a model to file
//====================================================================================================================
int	SC_Model_VQ::SaveModel(void) {
	int bytes;
	char *buffer;
	SV_DataIO io;

  //write trainingDataCount
	bytes = io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_VQ Model Failed!");}

	//write model
	bytes += io.writeScalar(&(this->DFile), this->pVQ!=NULL); //store if a model file is expected to exist
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_VQ Model Failed!");}
	if (this->pVQ != NULL) {
		buffer = sclib::exchangeFileExtension(this->lastUsedFileName, ".vq_model");
		bytes += this->pVQ->SaveModel();
		MFree_1D(buffer);
	}

  return bytes;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_VQ::modelOut(ostream& os) {
	os << this->pVQ;

	return(os);
}

//====================================================================================================================
// for computing BIC etc.
//====================================================================================================================
unsigned int SC_Model_VQ::getFreeParameterCount(void) {
	unsigned int cnt = this->pVQ->GetCNum() * this->pVQ->GetCDim();
	//                 codebook-size          feature-dimension
	return cnt;
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_VQ::drawSamplesFromDistribution(unsigned long int count) {
	int  d, D = this->pVQ->GetCDim(), T = this->pVQ->GetCNum();
	long int idx;
	SV_Data *pSamples = NULL;

	if (count > 0) {
		pSamples = new SV_Data(count, D);
		for (unsigned long int t = 0; t < count; t++) {
			idx = sclib::round(sclib::rand(0, T)); //each codebook-entry is equally likely
			for (d = 0; d < D; d++) {
				pSamples->Mat[t][d] = this->pVQ->GetCBook()[idx][d]; //it is copied directly because the model does not contain any statistical uncertainty that would require randomness here (just like a gaussian kernel with zero variance: always the mean will be returned)
			}
		}
	}

	return pSamples; 
}
