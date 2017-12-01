/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a Baggenstoss'-GMM as described in	*/
/*				'Robust Text-Independent Speaker Identification Using Gaussian  */
/*				'Mixture Speaker Models', D.Reynolds, R.C.Rose, 1995 (IEEE)			*/
/*				'Statistical Modeling Using Gaussian Mixtures and HMM with      */
/*         Matlab', P.M.Baggenstoss, 2002 (web)														*/
/*																																				*/
/*		Some Issues:																												*/
/*		  -	This model represents a GMM using P.M.Baggenstoss' algorithms		*/
/*				for training etc. in the full-covariance- and altered version   */
/*        for diagonal covariances. It is capable of deleting, splitting  */
/*        and merging mixtures during training/combining :-)							*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 29.03.2005																								*/
/**************************************************************************/

#include "SC_Aux.h"
#include "SC_MixtureModel_bGMM.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel_bGMM::SC_MixtureModel_bGMM(SC_TweakableParameters* pTweak, unsigned short int mixtureCount, unsigned short int dim, bool fullCovariance, bool verbose) : SC_MixtureModel(pTweak) {
	this->mixtureCount = sclib::max(mixtureCount, 0); 
	this->dim = dim;
	this->fullCovariance = fullCovariance;
	this->Hdr.ModelType	= sclib::mtBGMM;
  this->minVariance = NULL;
  this->maxMixtureCount = 100;
  this->maxEMsteps = this->pTweak->mixtureModelBgmm.maxEMiterations;
	this->verbose = verbose;

	this->weight = NULL; //unused in this class
	this->variance = NULL; //unused in this class
	this->sd = NULL; //unused in this class
	this->mean = NULL; //unused in this class

	if (this->fullCovariance == true) {
		this->pEM = new SC_BaggenstossEM(this->pTweak, this->verbose);
	} else {
		this->pEM = new SC_BaggenstossEMex(this->pTweak, this->verbose);
	}
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_MixtureModel_bGMM::SC_MixtureModel_bGMM(const SC_MixtureModel_bGMM& pParent) : SC_MixtureModel(pParent) {
  this->maxMixtureCount = pParent.maxMixtureCount;
	this->verbose = pParent.verbose;
	this->fullCovariance = pParent.fullCovariance;
	this->minVariance = NULL;
	setMinVariance(pParent.minVariance);

	if (this->fullCovariance == true) {
		this->pEM = new SC_BaggenstossEM((SC_BaggenstossEM &)*(pParent.pEM));
	} else {
		this->pEM = new SC_BaggenstossEMex((SC_BaggenstossEMex &)*(pParent.pEM));
	}
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_MixtureModel_bGMM::~SC_MixtureModel_bGMM() {
  MFree_1D(this->minVariance);
	MFree_0D(this->pEM);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_MixtureModel_bGMM& SC_MixtureModel_bGMM::operator=(const SC_MixtureModel_bGMM& pParent) {
	if (this != &pParent) {
		this->SC_MixtureModel::operator=(pParent);

		this->maxMixtureCount = pParent.maxMixtureCount;
		this->verbose = pParent.verbose;
		this->fullCovariance = pParent.fullCovariance;
		setMinVariance(pParent.minVariance);

		MFree_0D(this->pEM);
		if (this->fullCovariance == true) {
			this->pEM = new SC_BaggenstossEM((SC_BaggenstossEM &)*(pParent.pEM));
		} else {
			this->pEM = new SC_BaggenstossEMex((SC_BaggenstossEMex &)*(pParent.pEM));
		}
	}
  
  return *this;
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by baggenstoss.
//
// here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_MixtureModel_bGMM::TrainModel(SV_Data *pData) {
	return TrainModel(pData, 1, -1, false, 1.0, true, -1.0);
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by baggenstoss.
//
// here, the specified number of segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_MixtureModel_bGMM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	return TrainModel(pData, segmentsToMerge, -1, false, 1.0, true, -1.0);
}

//====================================================================================================================
//  calls the Baggenstoss' training procedure
//====================================================================================================================
int SC_MixtureModel_bGMM::TrainModel(SV_Data* pData, unsigned long int segmentsToMerge, long int samplesPerMode, bool bias, double maxCloseness, bool addModes, double splitThresh) {
  unsigned int iterationCount;
	SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
  
	if (this->fullCovariance == false) {
		this->pEM->init(pCompleteData, this->mixtureCount, this->pTweak->mixtureModelBgmm.varianceLimit, this->pTweak->mixtureModelBgmm.randomInitialization);
	} else {
		this->pEM->init(pCompleteData, this->mixtureCount, sqrt(this->pTweak->mixtureModelBgmm.varianceLimit), this->pTweak->mixtureModelBgmm.randomInitialization); //here, the minimal standard deviation is needed instead of the min. variance!
	}

	this->pEM->train(pCompleteData, iterationCount, this->maxEMsteps, samplesPerMode, bias, maxCloseness, addModes, splitThresh);
	this->trainingDataCount = pCompleteData->Row;
	this->mixtureCount = this->pEM->getMixtureCount(); //might change due to merging/splitting/deflation

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

	return iterationCount;
}

//====================================================================================================================
// Test GMM, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_bGMM::TestModel(SV_Data *pData) {
	SV_Data	*pScore = TestModel(pData, 1);

	return pScore;					
}

//====================================================================================================================
// Test GMM with a specified nr. of segments out of the lonked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_bGMM::TestModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  SV_Data	*pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(0); //merge all patterns in this linked list together
	SV_Data	*pScore = new SV_Data(1, 1);

	pScore->Mat[0][0] = (float)(this->pEM->test(pCompleteData) / (double)(pCompleteData->Row)); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...

	return (pScore);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_bGMM::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_MixtureModel_bGMM*)pSecond, false));
}
SC_MixtureModel_bGMM* SC_MixtureModel_bGMM::combineModels(SC_MixtureModel_bGMM* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_bGMM::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel_bGMM*)pFirst, (SC_MixtureModel_bGMM*)pSecond, keepFirstsNext));
}
SC_MixtureModel_bGMM* SC_MixtureModel_bGMM::combineModels(SC_MixtureModel_bGMM* pFirst, SC_MixtureModel_bGMM* pSecond, bool keepFirstsNext) {
	return NULL; //TODO: could be done better...
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_MixtureModel_bGMM::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_MixtureModel_bGMM*)pSecond, pSpeechFrames, segmentsToMerge, pBackgroundModels));
}
SC_MixtureModel_bGMM* SC_MixtureModel_bGMM::combineModels(SC_MixtureModel_bGMM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
	SC_MixtureModel_bGMM* pNewModel;

  assert(this->dim == pSecond->getDim());
	
	pNewModel = new SC_MixtureModel_bGMM(this->pTweak, sclib::min(this->mixtureCount+pSecond->getMixtureCount(), this->pTweak->modelHandler.maxSpeakerModelOrder), this->dim, this->fullCovariance, this->verbose);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

	return pNewModel;
}

//====================================================================================================================
// Set the minimum-variance from outside
//====================================================================================================================
void SC_MixtureModel_bGMM::setMinVariance(double *newMinVariance) {
  MFree_1D(this->minVariance); 

	if (newMinVariance != NULL) {
		MArray_1D(this->minVariance, this->dim, double, "SC_MixtureModel_bGMM.setMinVariance: minVariance"); 
		for (int x = 0; x < this->dim; x++) {
			this->minVariance[x] = newMinVariance[x];
		} 
	}

  return;
}

//====================================================================================================================
// Set the minimum-variance from outside
//====================================================================================================================
void SC_MixtureModel_bGMM::setMinVariance(double newMinVariance) {
  MFree_1D(this->minVariance); 

  MArray_1D(this->minVariance, this->dim, double, "SC_MixtureModel_bGMM.setMinVariance: minVariance"); 
  for (int x = 0; x < this->dim; x++) {
    this->minVariance[x] = newMinVariance;
  } 

  return;
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_MixtureModel_bGMM::SaveModel(void) {
 	int res, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_bGMM model Failed!");}

  //write scalars
	bytes = io.writeScalar(&(this->DFile), this->maxMixtureCount);
	bytes += io.writeScalar(&(this->DFile), this->verbose);
	bytes += io.writeScalar(&(this->DFile), this->fullCovariance);
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount);
	bytes += io.writeScalar(&(this->DFile), this->maxEMsteps);
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
	bytes += io.writeScalar(&(this->DFile), this->dim);
	if (this->DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_bGMM Model Failed!");}

  //write minVariance-vector
	bytes += io.writeArray(&(this->DFile), this->minVariance, this->pEM->getDim());
	if (this->DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_bGMM Model Failed!");}

	//write Baggenstoss-stuff
	bytes += this->pEM->write(&(this->DFile));

  return bytes + MHLen; //MHLen is just an estimate, not necessarily true...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model* SC_MixtureModel_bGMM::LoadModel(void) {
	int res;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

	//read scalars
	io.readScalar(&(this->DFile), this->maxMixtureCount, codeSizes, fileSizes);
	io.readScalar(&(this->DFile), this->verbose, codeSizes, fileSizes);
	io.readScalar(&(this->DFile), this->fullCovariance, codeSizes, fileSizes);
	io.readScalar(&(this->DFile), this->mixtureCount, codeSizes, fileSizes);
	io.readScalar(&(this->DFile), this->maxEMsteps, codeSizes, fileSizes);
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
	io.readScalar(&(this->DFile), this->dim, codeSizes, fileSizes);
	if (this->DFile.good() != TRUE) {return(NULL);}

	//read min-Variance vector
	io.readArray(&(this->DFile), this->minVariance, dim, codeSizes, fileSizes);
	if (this->DFile.good() != TRUE) {return(NULL);}

	//read Baggenstoss-stuff
	res = this->pEM->read(&(this->DFile), &fileSizes);
	if (this->DFile.good() != TRUE || res <= 0) {return(NULL);}

	return(this);
}

//====================================================================================================================
// for computing BIC etc.
//====================================================================================================================
unsigned int SC_MixtureModel_bGMM::getFreeParameterCount(void) {
	return (this->fullCovariance == true) ? this->pEM->getDim()*this->pEM->getDim()*this->pEM->getMixtureCount()+this->pEM->getDim()*this->pEM->getMixtureCount()+this->pEM->getMixtureCount() : 2*this->pEM->getDim()*this->pEM->getMixtureCount()+this->pEM->getMixtureCount();
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
// this method follows the instructions in "A Monte-Carlo Method for Score Normalization in Automatic Speaker 
// Verification using Kullback-Leibler Distances" by Ben, Blouet, Bimbot on ICASSP 2002
//====================================================================================================================
SV_Data* SC_MixtureModel_bGMM::drawSamplesFromDistribution(unsigned long int count) {
  unsigned long int t;
  unsigned short int i, d, dd;
	double *sample = NULL, value;
  SV_Data *pSamples = NULL;

  if (this->mixtureCount > 0 && this->dim > 0) {
    pSamples = new SV_Data(count, this->dim);

		if (this->fullCovariance == true) {
			MArray_1D(sample, this->dim, double, "SC_MixtureModel_bGMM.drawSamplesFromDistribution: sample");
		}
    
		for (t = 0; t < count; t++) {
      //1. randomly draw a mixture-index respecting the mixture weights
			i = sclib::drawIndexFromDistribution(this->pEM->getWeight(), this->mixtureCount);

      //2. generate a gaussian random variable according to the parameters of the i'th mixture selected above
			if (this->fullCovariance == false) {
				for (d = 0; d < this->dim; d++) {
					pSamples->Mat[t][d] = (float)(sclib::getRandomizer()->rand_gaus(this->pEM->getMean()[i][d], sqrt(((SC_BaggenstossEMex*)this->pEM)->getVariance()[i][d]))); //this uses the Box-Muller algorithm inside as in the paper
				}
			} else {
				for (d = 0; d < this->dim; d++) {
					sample[d] = sclib::getRandomizer()->rand_gaus(0.0, 1.0); //normally, i.e. N(0,1), distributed random numbers
				}
				for (d = 0; d < this->dim; d++) {
					value = this->pEM->getMean()[i][d]; //this adds the desired mean => we have pseudo random numbers according to N(mean, covar)!
					for (dd = 0; dd < this->dim; dd++) {
						value += sample[dd] * this->pEM->getCholeskyCovar()[i][dd][d]; //this is tmp * cholesky, which models the desired correlation (no transpose on cholesky, because it is already an upper triangular matrix here, see gaussian-faq by John D'Errico)
					}
					pSamples->Mat[t][d] = (float)(value);
				}
			} //if full-covar
		} //for t

		if (this->fullCovariance == true) {
			MFree_1D(sample);
		}
  }

  return pSamples;
}

//====================================================================================================================
// create a (linked, not copied) SC_Signature-view on this model (for distance computation )
//====================================================================================================================
SC_Signature* SC_MixtureModel_bGMM::toSignature(void) {
  SC_Centroid **centroids = NULL;
  SC_Signature *pSignature;
  
  MArray_1D(centroids, this->mixtureCount, SC_Centroid*, "SC_MixtureModel_bGMM.toSignature: centroids");

  for (unsigned short int i = 0; i < this->mixtureCount; i++)   {
		centroids[i] = new SC_Centroid_Gaussian(this->pTweak, this->dim, this->getMean(i), this->getVariance(i)); //TODO: getVariance() works only for the diagonal variance case
  }

	pSignature = new SC_Signature(centroids, this->getWeight(), this->mixtureCount);

  return pSignature;
}

//====================================================================================================================
// destruct a signature created by this model's toSignature()-method
//====================================================================================================================
void SC_MixtureModel_bGMM::killSignature(SC_Signature *pSignature) {
  MFree_0D(pSignature);

  return;
}
