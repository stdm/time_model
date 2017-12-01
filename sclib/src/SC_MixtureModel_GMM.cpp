/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a GMM as described in							*/
/*				'Robust Text-Independent Speaker Identification Using Gaussian  */
/*				'Mixture Speaker Models', D.Reynolds, R.C.Rose, 1995 (IEEE)			*/
/*																																				*/
/*		Some Issues:																												*/
/*			- The Model Uses only diagonal Covariance-Matrices, because this	*/
/*				has proven to yield best results in the above paper. To draw		*/
/*				maximum reward from this simplification, we use only a vector		*/
/*				to store the variances. This simplifies some calculations as		*/
/*				well as the memory-usage																				*/
/*			- Variance-limiting is applied.																		*/
/*			- Initialization of the model prior to estimating parameters			*/
/*			  isn't well researched; here, the means are build randomly, and  */
/*				some iterations of k-means algorithm privides somehow				    */
/*				meaningfull distribution of the data along the mixtures 				*/
/*				accordning to the paper.																				*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 06.04.2004																								*/
/**************************************************************************/

#include <math.h>
#include <float.h>
#include <assert.h>
#include "SC_Aux.h"
#include "SC_MixtureModel_GMM.h"
#include <SV_Error.h>

#include "SC_Timer.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel_GMM::SC_MixtureModel_GMM(SC_TweakableParameters* pTweak, unsigned short int mixtureCount, unsigned short int dim) : SC_MixtureModel(pTweak) {
	this->mixtureCount = mixtureCount; 
	this->dim = dim;
	this->Hdr.ModelType = sclib::mtGMM_new;
  this->maxEMsteps = this->pTweak->mixtureModelGmm.maxEMiterations;

	//normally, the mixtureCount must be >0, but there is one exception from this rule:
	//maybe we need only a "link" in the linked list of models, but have no data to model.
	//then, a model is created with only a "dummy" mixture: the mixture-count remains 0 
  //in this case to indicate the dummy model, but there is indeed one mixture component
  //with DUMMY_* parameters.
  if (this->mixtureCount > 0) {
    MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: sd");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: mean");
    MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_GMM: weight");
  	for (int i = 0; i < this->mixtureCount; i++) {
			this->weight[i] = 0.0;
		}
	} else {
    MArray_2D(this->variance, 1, this->dim, double, "SC_MixtureModel_GMM: variance");
    MArray_2D(this->sd, 1, this->dim, double, "SC_MixtureModel_GMM: varsdiance");
    MArray_2D(this->mean, 1, this->dim, double, "SC_MixtureModel_GMM: mean");
    MArray_1D(this->weight, 1, double, "SC_MixtureModel_GMM: weight");
    for (int d = 0; d < this->dim; d++) {
      this->variance[0][d] = this->dummyVariance;
			this->mean[0][d] = this->dummyMean;
      this->sd[0][d] = sqrt(this->dummyVariance);
    }
		this->weight[0] = this->dummyWeight;
	}
}

//====================================================================================================================
// constructor to simply fed already existing parameters (for synthesis, e.g.)
//====================================================================================================================
SC_MixtureModel_GMM::SC_MixtureModel_GMM(SC_TweakableParameters *pTweak, unsigned short int mixtureCount, unsigned short int dim, double *weight, double **mean, double **variance) : SC_MixtureModel(pTweak) {
	this->mixtureCount = mixtureCount; 
	this->dim = dim;
	this->Hdr.ModelType = sclib::mtGMM_new;
  this->maxEMsteps = this->pTweak->mixtureModelGmm.maxEMiterations;

  MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: variance");
  MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: varsdiance");
  MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: mean");
  MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_GMM: weight");
	for (int i = 0; i < this->mixtureCount; i++) {
		for (int d = 0; d < this->dim; d++) {
			this->variance[i][d] = variance[i][d];
			this->mean[i][d] = mean[i][d];
			this->sd[i][d] = sqrt(this->variance[i][d]);
		}
		this->weight[i] = weight[i];
	}
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_MixtureModel_GMM::SC_MixtureModel_GMM(const SC_MixtureModel_GMM& pParent) : SC_MixtureModel(pParent) {
	if (this->mixtureCount <= 0) { //see explanation in the constructor for these weird things...
    MArray_2D(this->variance, 1, this->dim, double, "SC_MixtureModel_GMM: variance");
    MArray_2D(this->sd, 1, this->dim, double, "SC_MixtureModel_GMM: sd");
    MArray_2D(this->mean, 1, this->dim, double, "SC_MixtureModel_GMM: mean");
    MArray_1D(this->weight, 1, double, "SC_MixtureModel_GMM: weight");
    for (int d = 0; d < this->dim; d++) {
      this->variance[0][d] = this->dummyVariance;
      this->sd[0][d] = sqrt(this->dummyVariance);
			this->mean[0][d] = this->dummyMean;
    }
		this->weight[0] = this->dummyWeight;
	}
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_MixtureModel_GMM::~SC_MixtureModel_GMM() {

}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_MixtureModel_GMM& SC_MixtureModel_GMM::operator=(const SC_MixtureModel_GMM& pParent) {
	if (this != &pParent) {
		this->SC_MixtureModel::operator=(pParent);

		if (this->mixtureCount <= 0) {
			MArray_2D(this->variance, 1, this->dim, double, "SC_MixtureModel_GMM: variance");
			MArray_2D(this->sd, 1, this->dim, double, "SC_MixtureModel_GMM: sd");
			MArray_2D(this->mean, 1, this->dim, double, "SC_MixtureModel_GMM: mean");
			MArray_1D(this->weight, 1, double, "SC_MixtureModel_GMM: weight");
			for (int d = 0; d < this->dim; d++) {
				this->variance[0][d] = this->dummyVariance;
				this->sd[0][d] = sqrt(this->dummyVariance);
				this->mean[0][d] = this->dummyMean;
			}
			this->weight[0] = this->dummyWeight;
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
int SC_MixtureModel_GMM::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by rose/reynolds.
//
// here, one can specify how many segments of a linked list of feature-vectors should be merged
//====================================================================================================================
int SC_MixtureModel_GMM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	double Log_p_X = 0.0, Log_p_X_old; //the log-likelihood of the data given the actual/last model: p(X|lambda)
	double *p_i; //class-probability for acoustic classes (mixture-component) i: p(i|x_t, lambda) for specific t and all i
	double *Sum_p_i; //SUM_t=1..T p(i|x_t, lambda) for all i 
	double **Sum_p_i_x; //SUM_t=1..T p(i|x_t, lambda)*x_t for all i 
	double **Sum_p_i_xx; //SUM_t=1..T p(i|x_t, lambda)*(x_t)^2 for all i 
	double *pb; //values of gauss-densities multiplied with the class-weights for a specific feature-vector and all mixtures: p_i * b_i(x_t)
	double Sum_pb; //SUM_i=1..M p_i * b_i(x_t) for a specific t
	double Log_Prod_pb; //log(PROD_t=1..T SUM_i=1..M p_i * b_i(x_t))
	unsigned long int t; //0..T over all feature-vectors
	unsigned short int i; //0..M over all mixture-components
	unsigned short int d; //0..D over all components of the feature-vectors
	unsigned long int iterationCount = 0;
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
  double **scaledVariance;

	MArray_1D(p_i, this->mixtureCount, double, "SC_MixtureModel_GMM.TrainModel: p_i");
	MArray_1D(pb, this->mixtureCount, double, "SC_MixtureModel_GMM.TrainModel: pb");
	MArray_1D(Sum_p_i, this->mixtureCount, double, "SC_MixtureModel_GMM.TrainModel: Sum_p_i");
	MArray_2D(Sum_p_i_x, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM.TrainModel: Sum_p_i_x");
	MArray_2D(Sum_p_i_xx, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM.TrainModel: Sum_p_i_xx");
  MArray_2D(scaledVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM.TrainModel: scaledVariance");
	
	//The model-parameters must be initialized some way prior to the em-algorithm
	initParameters(pCompleteData, this->pTweak->mixtureModelGmm.varianceLimit, this->pTweak->mixtureModelGmm.kMeansIterations);
	if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbModelCreation) == true) {sclib::classOut("gmm.txt", this, this->pTweak);}

  //initialize the scaled variance to speed up gaussian computation
	for (i = 0; i < this->mixtureCount; i++) {
  	for (d = 0; d < this->dim; d++) {
      scaledVariance[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
		}			
  }

	//now we can start... ladys 'n gentleman: the EM-Algorithm ;-)
	do {
		Log_p_X_old = Log_p_X;
		iterationCount++;

		//calculate all SUM/PROD terms of the formulas a-priori in some loops:
		Log_Prod_pb = 0.0;
		for (i = 0; i < this->mixtureCount; i++) { //1. initialize with zero... needed by the third loop
			Sum_p_i[i] = 0.0;
			for (d = 0; d < this->dim; d++) {
				Sum_p_i_x[i][d] = 0.0;
				Sum_p_i_xx[i][d] = 0.0;
			}
		}

		for (t = 0; t < (unsigned long int)pCompleteData->Row; t++) {
			for (i = 0; i < this->mixtureCount; i++) { //2. compute the sum_pb over all mixtures i needed for p_X and the third loop
        pb[i] = log(this->weight[i]);
				for (d = 0; d < this->dim; d++) { //the multivariate gaussian pdf may be evalutated this way (component-wise), because we have diagonal covariances (=> features are assumed uncorrelated)!
					pb[i] += this->gaussSolver.fastLogGaussian(pCompleteData->Mat[t][d], this->mean[i][d], this->variance[i][d], scaledVariance[i][d]);
          //pb[i] += scaledVariance[i][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[i][d]) * (pCompleteData->Mat[t][d] - this->mean[i][d]) / this->variance[i][d]);
        }
        Sum_pb = (i == 0) ? pb[i] : sclib::sLogAdd(Sum_pb, pb[i]);
			}//i=0..M

      Log_Prod_pb += Sum_pb;

			for (i = 0; i < this->mixtureCount; i++) { //3. compute the p_i's needed by the reestimation formulas
        p_i[i] = sclib::sExp(pb[i] - Sum_pb);
        Sum_p_i[i] += p_i[i];
        for (d = 0; d < this->dim; d++) {
					Sum_p_i_x[i][d] += p_i[i] * pCompleteData->Mat[t][d];
          Sum_p_i_xx[i][d] += p_i[i] * pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d];
				}
			} //i=0..M
		} //t=0..T

		Log_p_X = Log_Prod_pb;

		//let's start with the reestimation process:
		for (i = 0; i < this->mixtureCount; i++) {
			this->weight[i] = Sum_p_i[i] / pCompleteData->Row; //p_i* = 1/T SUM_t=1..T p(i|x_t, lambda)
      //assert(this->weight[i] > 0.0);
      //assert(Sum_p_i[i] != 0.0);
      if (this->weight[i] > 0.0) {
  		  for (d = 0; d < this->dim; d++) {
				  this->mean[i][d] = Sum_p_i_x[i][d] / Sum_p_i[i]; //mu_i* = (SUM_t=1..T p(i|x_t, lambda)*x_t) / (SUM_t=1..T p(i|x_t, lambda))
				  this->variance[i][d] = (Sum_p_i_xx[i][d] / Sum_p_i[i]) - (this->mean[i][d] * this->mean[i][d]); //sigma_i* = ((SUM_t=1..T p(i|x_t, lambda)*(x_t)^2) / (SUM_t=1..T p(i|x_t, lambda))) - (mu_i*)^2
				  if (this->variance[i][d] < this->pTweak->mixtureModelGmm.varianceLimit) {this->variance[i][d] = this->pTweak->mixtureModelGmm.varianceLimit;} //do variance-limiting
          this->sd[i][d] = sqrt(this->variance[i][d]);
          scaledVariance[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
			  }			
      } else {
        killMixture(i);
        printf("\n     mixture %i killed in run %i, %i remaining", i, iterationCount, this->mixtureCount);
      }
  	}
  	
		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("gmm.txt", this, this->pTweak);}
		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::scalarOut("gmm_likelihood.txt", Log_p_X, this->pTweak);}

    //assert(Log_p_X >= Log_p_X_old || iterationCount == 1);
    //if (!(Log_p_X >= Log_p_X_old || iterationCount == 1)) {
    //  printf("\n     abnormal log-likelihood: %f (old: %f)", Log_p_X, Log_p_X_old);
    //}
    
    if (iterationCount == this->maxEMsteps) {
      break;
    }

	} while (fabs(Log_p_X - Log_p_X_old) > fabs(this->pTweak->mixtureModelGmm.EMthreshold));

	this->trainingDataCount = pCompleteData->Row;

  //kill too small mixtures
  this->sortParameters();
  for (i = 0; i < this->mixtureCount; i++) {
    if (this->weight[i] < this->pTweak->mixtureModelGmm.weightLimit) { //kill too small mixtures
      this->killMixture(i);
      i--;
    }
  }

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }
	MFree_1D(p_i);
	MFree_1D(pb);
	MFree_1D(Sum_p_i);
	MFree_2D(Sum_p_i_x);
	MFree_2D(Sum_p_i_xx);
  MFree_2D(scaledVariance);
	
	return iterationCount;
}

//====================================================================================================================
// Test GMM, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_GMM::TestModel(SV_Data *pData) {
  SV_Data	*pScore = TestModel(pData, 1);

  return pScore;					
}

//====================================================================================================================
// Test GMM with a specified nr. of segments out of the linked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_GMM::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge); //merge all patterns in this linked list together
	unsigned long int t, T = (unsigned long int)pCompleteData->Row; //0..T over all feature-vectors
	unsigned short int i, I = this->mixtureCount, maxMixture; //0..M over all mixture-components
	unsigned short int d, D = this->dim; //0..D over all components of the feature-vectors
	double Log_Prod_pb; //PROD_t=1..T SUM_i=1..M p_i * b_i(x_t)
	double Sum_pb; //SUM_i=1..M p_i * b_i(x_t) for a specific t
	double pb; //value of gauss-densities multiplied with the class-weights for a specific feature-vector and all mixtures: p_i * b_i(x_t)
	SV_Data *pScore = new SV_Data(pCompleteData->Row+1, 1);
  double **scaledVariance, maxLL;

  //initialize the scaled variance to speed up gaussian computation
  MArray_2D(scaledVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM.TrainModel: scaledVariance");
	for (i = 0; i < this->mixtureCount; i++) {
  	for (d = 0; d < this->dim; d++) {
      scaledVariance[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
		}			
  }  

  //compute the likelihood
  Log_Prod_pb = 0.0;
	for (t = 0; t < T; t++) {
		maxLL = -1.0 * std::numeric_limits<double>::max();
		for (i = 0; i < I;  i++) {
      pb = log(this->weight[i]);
			for (d = 0; d < D; d++) {
				pb += this->gaussSolver.fastLogGaussian(pCompleteData->Mat[t][d], this->mean[i][d], this->variance[i][d], scaledVariance[i][d]);
        //pb += scaledVariance[i][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[i][d]) * (pCompleteData->Mat[t][d] - this->mean[i][d]) / this->variance[i][d]);
      }
			if (pb > maxLL) {
				maxLL = pb;
				maxMixture = i;
			}
      Sum_pb = (i == 0) ? pb : sclib::sLogAdd(Sum_pb, pb);
		}
		pScore->Mat[t+1][0] = (float)(maxMixture); //remember the mixture component that was closest to vector t in the result; (i.e. the clustering result)
    Log_Prod_pb += Sum_pb;
	}
	
  if (TestData != pCompleteData) {
	  MFree_0D(pCompleteData);
  }
  MFree_2D(scaledVariance);
	pScore->Mat[0][0] = (float)Log_Prod_pb / (float)T; //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...;

  return (pScore);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_GMM::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_MixtureModel*)pSecond, false));
}
SC_MixtureModel_GMM* SC_MixtureModel_GMM::combineModels(SC_MixtureModel* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_GMM::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_GMM* SC_MixtureModel_GMM::combineModels(SC_MixtureModel* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext) {
  unsigned short int i, d, newMixtureCount;
	unsigned long int newTrainingDataCount;
	double *newWeight, **newMean, **newVariance, **newSd;
	SC_MixtureModel_GMM *pNewModel;
  bool tooSmall = false;
	
  //maybe both of the models are dummys. then, just return a single dummy.
  if (pSecond == NULL || (pFirst->getMixtureCount() == 0 && pSecond->getMixtureCount() == 0)) {
    pNewModel = new SC_MixtureModel_GMM(*(SC_MixtureModel_GMM*)pFirst);
  } else {
    newMixtureCount = sclib::max(1, pFirst->getMixtureCount()) + sclib::max(1, pSecond->getMixtureCount());

    //if none of both is a dummy, they must have equal feature-dim
    assert(pFirst->getDim() == pSecond->getDim());

    //determine new traingdata count
    if (pFirst->getMixtureCount() == 0) {pFirst->setTrainindDataCount(pSecond->getTrainingDataCount());}
    if (pSecond->getMixtureCount() == 0) {pSecond->setTrainindDataCount(pFirst->getTrainingDataCount());}
    newTrainingDataCount = pFirst->getTrainingDataCount() + pSecond->getTrainingDataCount();

	  MArray_1D(newWeight, newMixtureCount, double, "SC_MixtureModel_GMM.combineModels: newWeight");
	  MArray_2D(newMean, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_GMM.combineModels: newMean");
	  MArray_2D(newVariance, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_GMM.combineModels: newVariance");
    MArray_2D(newSd, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_GMM.combineModels: newSd");
    pNewModel = new SC_MixtureModel_GMM(pFirst->getTweak(), newMixtureCount, pFirst->getDim());

	  //copy the mixture-components together
	  for (i = 0; i < newMixtureCount; i++) {
		  if (i < sclib::max(1, pFirst->getMixtureCount())) {
			  newWeight[i] = (pFirst->getWeight(i) * pFirst->getTrainingDataCount()) / newTrainingDataCount; //re-normalize the weights so that their sum = 1.0
			  for (d = 0; d < pFirst->getDim(); d++) {
				  newMean[i][d] = pFirst->getMean(i, d);
				  newVariance[i][d] = pFirst->getVariance(i, d);
          newSd[i][d] = pFirst->getSd(i, d);
			  }
		  } else {
			  newWeight[i] = (pSecond->getWeight(i - sclib::max(1, pFirst->getMixtureCount())) * pSecond->getTrainingDataCount()) / newTrainingDataCount; //re-normalize the weights so that their sum = 1.0
			  for (d = 0; d < pFirst->getDim(); d++) {
				  newMean[i][d] = pSecond->getMean(i - sclib::max(1, pFirst->getMixtureCount()), d);
				  newVariance[i][d] = pSecond->getVariance(i - sclib::max(1, pFirst->getMixtureCount()), d);
          newSd[i][d] = pSecond->getSd(i - sclib::max(1, pFirst->getMixtureCount()), d);
			  }
		  }
      if (newWeight[i] < pFirst->getTweak()->mixtureModelGmm.weightLimit) {tooSmall = true;}
	  }

	  //make the new mixture-components the new parameters of the new model
    pNewModel->setWeight(newWeight);
    pNewModel->setMean(newMean);
    pNewModel->setVariance(newVariance);
    pNewModel->setSd(newSd);
    pNewModel->setTrainindDataCount(newTrainingDataCount);
    pNewModel->setMixtureCount(newMixtureCount);

    //kill too small mixtures
    if (tooSmall == true) {
      pNewModel->sortParameters();
      for (i = 0; i < pNewModel->getMixtureCount(); i++) {
        if (pNewModel->getWeight(i) < pFirst->getTweak()->mixtureModelGmm.weightLimit) {
          pNewModel->killMixture(i);
          i--;
        }
      }
    }
  
  } //min. 1 model was no dummy

  if (keepFirstsNext == true) {
    pNewModel->Next = pFirst->Next;
  } else {
    pNewModel->Next = NULL;
  }

	return pNewModel;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_MixtureModel_GMM::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pSecond, pSpeechFrames, segmentsToMerge, pBackgroundModels));
}
SC_MixtureModel_GMM* SC_MixtureModel_GMM::combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
	SC_MixtureModel_GMM* pNewModel;

  assert(this->dim == pSecond->getDim());
	
	pNewModel = new SC_MixtureModel_GMM(this->pTweak, sclib::min(this->mixtureCount+pSecond->getMixtureCount(), this->pTweak->modelHandler.maxSpeakerModelOrder), this->dim);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

	return pNewModel;
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes written, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_MixtureModel_GMM::SaveModel(void) {
 	int res, x, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_GMM model Failed!");}

  //write mixtureCount and featureDim
	bytes = io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_GMM Model Failed!");}

  //write weight-vector
	bytes += io.writeArray(&(this->DFile), this->weight, sclib::max(1, this->mixtureCount));
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_GMM Model Failed!");}

	//write mean and variance matrices
	for (x = 0; x < sclib::max(1, this->mixtureCount); x++) {
		bytes += io.writeArray(&(this->DFile), this->mean[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->variance[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->sd[x], this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_GMM Model Failed!");}
	}

  //write trainingDataCount
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_GMM Model Failed!");}

  return bytes + MHLen; //MHLen may only be an estimate...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model* SC_MixtureModel_GMM::LoadModel(void) {
	int res, x; 
	unsigned short int newDim = 0, newMixtureCount = 0;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

	//read mixtureCount and featureDim
	io.readScalar(&(this->DFile), newDim, codeSizes, fileSizes);
	if ((DFile.good() != TRUE) || (newDim == 0)) {return(NULL);}

	io.readScalar(&(this->DFile), newMixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {return(NULL);}

	//prepare this model for changing it's parameters
	if ((newDim != this->dim) || (newMixtureCount != this->mixtureCount)) {
		MFree_1D(this->weight);
		MFree_2D(this->mean);
		MFree_2D(this->variance);
    MFree_2D(this->sd);

		this->dim = newDim;
		this->mixtureCount = newMixtureCount;

 		MArray_1D(this->weight, sclib::max(1, this->mixtureCount), double, "SC_MixtureModel_GMM.loadModel: weight");
    MArray_2D(this->mean, sclib::max(1, this->mixtureCount), this->dim, double, "SC_MixtureModel_GMM.loadModel: mean");
		MArray_2D(this->variance, sclib::max(1, this->mixtureCount), this->dim, double, "SC_MixtureModel_GMM.loadModel: variance");
    MArray_2D(this->sd, sclib::max(1, this->mixtureCount), this->dim, double, "SC_MixtureModel_GMM.loadModel: sd");
	};

	//read weight, mean and variance
	io.readArray(&(this->DFile), this->weight, sclib::max(1, this->mixtureCount), codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_GMM Failed!");}
	for (x = 0; x < sclib::max(1, this->mixtureCount); x++) {
		io.readArray(&(this->DFile), this->mean[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->variance[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->sd[x], this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_GMM Failed!");}
	}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_GMM Failed!");}

  this->maxEMsteps = 150;

	return(this);
}
