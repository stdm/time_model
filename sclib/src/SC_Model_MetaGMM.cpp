/**************************************************************************/
/*    This is a "meta-GMM": it holds one special GMM per dimension, so    */
/*    that each dimension can have a speratae mixture count. Maybe this   */
/*    helps to save some parameters.                                      */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 07.09.2009																								*/
/**************************************************************************/

#include "SC_Model_MetaGMM.h"
#include "SC_ModelHandler.h"
#include "SC_Centroid_Signature.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_MetaGMM::SC_Model_MetaGMM(SC_TweakableParameters* pTweak, unsigned short int maxMixtureCount, unsigned short int dim) : SC_Model(pTweak) {
	this->maxMixtureCount = maxMixtureCount; 
	this->dim = dim;
	this->Hdr.ModelType = sclib::mtMetaGMM;
  this->mixtureCount = NULL;
	this->pGMM = NULL;
	this->orthTransMat = NULL;
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model_MetaGMM::SC_Model_MetaGMM(const SC_Model_MetaGMM& pParent) : SC_Model(pParent) {
	this->dim = pParent.dim;
	this->maxMixtureCount = pParent.maxMixtureCount;

	if (pParent.pGMM != NULL) {
		MArray_1D(this->mixtureCount, this->dim, unsigned short int, "SC_Model_MetaGMM: mixtureCount");
		MArray_1D(this->pGMM, this->dim, SC_MixtureModel_GMM*, "SC_Model_MetaGMM: pGMM");
		MArray_2D(this->orthTransMat, this->dim, this->dim, double, "SC_Model_MetaGMM: orthTransMat");

		SC_ModelHandler handler(this->pTweak, false);
		for (int d = 0; d < this->dim; d++) {
			this->pGMM[d] = (SC_MixtureModel_GMM*)(handler.copyModel(pParent.pGMM[d], false));
			for (int dd = 0; dd < this->dim; dd++) {
				this->orthTransMat[d][dd] = pParent.orthTransMat[d][dd];
			}
		}
	} else {
		this->mixtureCount = NULL;
		this->pGMM = NULL;
		this->orthTransMat = NULL;
	}
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_MetaGMM::~SC_Model_MetaGMM() {
	MFree_1D(this->mixtureCount);
	MFree_2D(this->orthTransMat);
	if (this->pGMM != NULL) {
		for (int d = 0; d < this->dim; d++) {
			MFree_0D(this->pGMM[d]);
		}
		MFree_1D(this->pGMM);
	}
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_MetaGMM& SC_Model_MetaGMM::operator=(const SC_Model_MetaGMM& pParent) {
	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		this->dim = pParent.dim;
		this->maxMixtureCount = pParent.maxMixtureCount;
		MFree_2D(this->orthTransMat);
		MFree_1D(this->mixtureCount);
		if (this->pGMM != NULL) {
			for (int d = 0; d < this->dim; d++) {
				MFree_0D(this->pGMM[d]);
			}
			MFree_1D(this->pGMM);
		}

		if (pParent.pGMM != NULL) {
			MArray_1D(this->mixtureCount, this->dim, unsigned short int, "SC_Model_MetaGMM: mixtureCount");
			MArray_1D(this->pGMM, this->dim, SC_MixtureModel_GMM*, "SC_Model_MetaGMM: pGMM");
			MArray_2D(this->orthTransMat, this->dim, this->dim, double, "SC_Model_MetaGMM: orthTransMat");

			SC_ModelHandler handler(this->pTweak, false);
			for (int d = 0; d < this->dim; d++) {
				this->pGMM[d] = (SC_MixtureModel_GMM*)(handler.copyModel(pParent.pGMM[d], false));
				for (int dd = 0; dd < this->dim; dd++) {
					this->orthTransMat[d][dd] = pParent.orthTransMat[d][dd];
				}
			}  
		}
	}
  
  return *this;
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by rose/reynolds.
//
// here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_MetaGMM::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by rose/reynolds.
//
// here, one can specify how many segments of a linked list of feature-vectors should be merged
//====================================================================================================================
int SC_Model_MetaGMM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	MFree_0D(this->mixtureCount);
	MFree_2D(this->orthTransMat);
	MArray_1D(this->mixtureCount, this->dim, unsigned short int, "SC_Model_MetaGMM.TrainModel: mixtureCount");
	if (this->pGMM != NULL) {
		for (int d = 0; d < this->dim; d++) {
			MFree_0D(this->pGMM[d]);
		}
		MFree_1D(this->pGMM);
	}
	MArray_1D(this->pGMM, this->dim, SC_MixtureModel_GMM*, "SC_Model_MetaGMM.TrainModel: pGMM");

	// Orthogonal transform data as in SV_Model_GMM to make dimensions nearly independant (see Jialong He's ICASSP'99 paper:
	//"On the use of orthogonal GMM in speaker recognition", L.Liu and J.He, ICASSP 1999)
	MArray_2D(this->orthTransMat, this->dim, this->dim, double, "SC_Model_MetaGMM.TrainModel: orthTransMat");
	SV_Data *pCompleteData = pData->MergeData(segmentsToMerge);
	double **CMat = CovMatrix(pCompleteData); // overall covariance matrix
	double *EigVal;
	MArray_1D(EigVal, this->dim, double, "SC_Model_MetaGMM.TrainModel: EigVal");
	GN_Matrix MEng;
	MEng.Eigen(CMat, this->orthTransMat, EigVal, this->dim);  
	OrthTrans(pCompleteData, this->orthTransMat);
	MFree_1D(EigVal);
	MFree_2D(CMat);

	SC_ModelHandler handler(this->pTweak, false);
	int res = SVLIB_Ok;
	for (int d = 0; d < this->dim; d++) {
		SV_Data *pOneDim = pCompleteData->MergeData(segmentsToMerge, d); //just to extract dimension d, no merging necessary
		if (d == 0) {
			this->trainingDataCount = pOneDim->Row;
		}

		this->mixtureCount[d] = handler.guessModelOrder(pOneDim, NULL, sclib::mtGMM_new, 1, this->maxMixtureCount, 1);
		this->pGMM[d] = (SC_MixtureModel_GMM*)(handler.buildModel(pOneDim, NULL, this->mixtureCount[d], sclib::mtGMM_new, 1));

		if (this->pGMM[d] == NULL) {
			res = SVLIB_Fail;
			break;
		}
		MFree_0D(pOneDim);
	}

	MFree_0D(pCompleteData);
	
	//TODO
	sclib::vectorOut("mixtureCounts.txt", this->mixtureCount, this->dim, false, this->pTweak);

	return res;
}

//====================================================================================================================
// Test GMM, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_Model_MetaGMM::TestModel(SV_Data *pData) {
  SV_Data	*pScore = TestModel(pData, 1);

  return pScore;					
}

//====================================================================================================================
// Test GMM with a specified nr. of segments out of the linked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_Model_MetaGMM::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
	SV_Data *pScore = NULL;
	int T;

	if (this->pGMM != NULL) {
		pScore = new SV_Data(1, 1);
		pScore->Mat[0][0] = 0.0f;
		SV_Data *pCompleteData = TestData->MergeData(segmentsToMerge);
		this->OrthTrans(pCompleteData, this->orthTransMat);

		for (int d = 0; d < this->dim; d++) {
			SV_Data *pOneDim = pCompleteData->MergeData(segmentsToMerge, d); //just extract dimension d, no merging necessary
			if (d == 0) {
				T = pOneDim->Row;
			}
			SV_Data *pTmpScore = this->pGMM[d]->TestModel(pOneDim);
			pScore->Mat[0][0] += pTmpScore->Mat[0][0]*(float)(pOneDim->Row);
			MFree_0D(pOneDim);
			MFree_0D(pTmpScore);
		}

		MFree_0D(pCompleteData);
		pScore->Mat[0][0] /= (float)(T);
	}

  return pScore;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_MetaGMM::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_MetaGMM*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_MetaGMM*)pBackgroundModels));
}
SC_Model_MetaGMM* SC_Model_MetaGMM::combineModels(SC_Model_MetaGMM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_MetaGMM* pBackgroundModels) {
	SC_Model_MetaGMM *pNewModel = new SC_Model_MetaGMM(this->pTweak, this->maxMixtureCount, this->dim);

	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  return pNewModel;
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes written, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_Model_MetaGMM::SaveModel(void) {
 	int res, bytes, cnt = 0;
	SV_DataIO io;
	char *buffer, extension[sclib::bufferSize];

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM model Failed!");}

  //write trainingDataCount
	bytes = io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM Model Failed!");}

  //write own stuff
	bytes += io.writeScalar(&(this->DFile), this->maxMixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->dim);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM Model Failed!");}
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount!=NULL); //is trained?
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM Model Failed!");}
	if (this->mixtureCount != NULL) {
		bytes += io.writeArray(&(this->DFile), this->mixtureCount, this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM Model Failed!");}
		bytes += io.writeMatrix(&(this->DFile), this->orthTransMat, this->dim, this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM Model Failed!");}
	}

	//write individual sub-model files
	if (this->pGMM != NULL) {
		for (int d = 0; d < this->dim; d++) {
			sprintf(extension, ".dimension_model_%d", d);
			buffer = sclib::exchangeFileExtension(this->lastUsedFileName, extension);
			this->pGMM[d]->OpenFile(buffer, WRITE_MODEL);
			res = this->pGMM[d]->SaveModel();
			if (res <= 0) {
				REPORT_ERROR(SVLIB_Fail, "Save SC_Model_MetaGMM sub-model Failed!");
			} else {
				bytes += res;
			}
			this->pGMM[d]->CloseFile();
			MFree_1D(buffer);
		}
	}

  return bytes + MHLen; //MHLen may be incorrect...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model* SC_Model_MetaGMM::LoadModel(void) {
	int res;
	char *buffer = NULL, extension[sclib::bufferSize];
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);
	SC_ModelHandler handler(this->pTweak, false);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM Failed!");}
	
	//read own members
	io.readScalar(&(this->DFile), this->maxMixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM Model Failed!");}
	io.readScalar(&(this->DFile), this->dim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM Model Failed!");}

	bool wasTrained;
	io.readScalar(&(this->DFile), wasTrained, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM Model Failed!");}

	MFree_2D(this->orthTransMat);
	MFree_1D(this->mixtureCount);
	if (this->pGMM != NULL) {
		for (int d = 0; d < this->dim; d++) {
			MFree_0D(this->pGMM[d]);
		}
		MFree_1D(this->pGMM);
	}
	if (wasTrained == true) {
		MArray_1D(this->mixtureCount, this->dim, unsigned short int, "SC_Model_MetaGMM.LoadModel: mixtureCount");
		MArray_2D(this->orthTransMat, this->dim, this->dim, double, "SC_Model_MetaGMM: orthTransMat");

		io.readArray(&(this->DFile), this->mixtureCount, this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM Model Failed!");}
		io.readMatrix(&(this->DFile), this->orthTransMat, this->dim, this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM Model Failed!");}

		MArray_1D(this->pGMM, this->dim, SC_MixtureModel_GMM*, "SC_Model_MetaGMM.LoadModel: pGMM");
		SC_ModelHandler handler(this->pTweak, false);
		for (int d = 0; d < this->dim; d++) {
			sprintf(extension, ".dimension_model_%d", d);
			buffer = sclib::exchangeFileExtension(this->lastUsedFileName, extension);
			this->pGMM[d] = (SC_MixtureModel_GMM*)(handler.createRawModel(sclib::mtGMM_new, NULL, 0, 1));
			this->pGMM[d]->OpenFile(buffer, READ_MODEL);
			this->pGMM[d] = (SC_MixtureModel_GMM*)(this->pGMM[d]->LoadModel());
			if (this->pGMM[d] != NULL) {
				this->pGMM[d]->CloseFile();
			} else {
				REPORT_ERROR(SVLIB_Fail, "Load SC_Model_MetaGMM dimension-model Failed!");
			}
			MFree_1D(buffer);
		}
	}
		
	return(this);
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_MetaGMM::modelOut(ostream& os) {
	os << "Traing-Data count:\t" << this->trainingDataCount << endl;
	os << "Traing-Data dimensionality:\t" << this->dim << endl;
	if (this->pGMM != NULL) {
		os << "Mixture-Counts per dimension:\t";
		for (int d = 0; d < this->dim; d++) {
			os << this->mixtureCount[d] << "\t";
		}
		os << endl;

		os << "Orthogonal transform matrix:" << endl;
		for (int y = 0; y < this->dim; y++) {
			for (int x = 0; x < this->dim; x++) {
				os << setw(8) << this->orthTransMat[y][x] << "\t";
			}
			os << endl;
		}
		os << endl;

		os << endl;
		for (int d = 0; d < this->dim; d++) {
			os << "Dimension:\t" << d << endl;
			os << this->pGMM[d] << endl;
		}
		os << endl;
	}

	return(os);
}

//====================================================================================================================
// for computing BIC etc.
//====================================================================================================================
unsigned int SC_Model_MetaGMM::getFreeParameterCount(void) {
	unsigned int cnt = 0;

	for (int d = 0; d < this->dim; d++) {
		cnt += this->pGMM[d]->getFreeParameterCount();
	}

	return cnt;
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_MetaGMM::drawSamplesFromDistribution(unsigned long int count) {
	SV_Data *pSamples = NULL;

	if (this->pGMM != NULL) {
		pSamples = new SV_Data(count, this->dim);
	
		for (int d = 0; d < this->dim; d++) { //draw each dimension independantly
			SV_Data *pOneDim = this->pGMM[d]->drawSamplesFromDistribution(count);
			for (unsigned long int t = 0; t < count; t++) {
				pSamples->Mat[t][d] = pOneDim->Mat[t][0];
			}
			MFree_0D(pOneDim);
		}

		SC_MatrixFunctions matFunc;
		double **invTransMat = matFunc.copy(this->orthTransMat, this->dim, this->dim);
		int res = matFunc.inv(this->orthTransMat, this->dim);
		this->OrthTrans(pSamples, invTransMat); //transform back to original coordinate system
		MFree_2D(invTransMat);
	}

	return pSamples; 
}

//====================================================================================================================
// generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
// method needs to return false;
//====================================================================================================================
bool SC_Model_MetaGMM::scoreIsAverage(void) {
	return true;
}

//====================================================================================================================
// create a (linked, not copied) SC_Signature-view on this model (for distance computation )
//====================================================================================================================
SC_Signature* SC_Model_MetaGMM::toSignature(void) {
  SC_Centroid **centroids = NULL;
  SC_Signature *pSignature;
  
  MArray_1D(centroids, this->dim, SC_Centroid*, "SC_Model_MetaGMM.toSignature: centroids");
  for (unsigned short int i = 0; i < this->dim; i++)   {
		centroids[i] = new SC_Centroid_Signature(this->pTweak, this->pGMM[i]->toSignature(), false); //do not justLink, so that the created signature gets destroyed when the centroid dies
  }
	pSignature = new SC_Signature(centroids, this->dim);

  return pSignature;
}

//====================================================================================================================
// destruct a signature created by this model's toSignature()-method
//====================================================================================================================
void SC_Model_MetaGMM::killSignature(SC_Signature *pSignature) {
  MFree_0D(pSignature);

  return;
}