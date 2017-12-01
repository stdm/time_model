/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a GMM-UBM as described in				  */
/*				'Speaker Verification Using Adapted Gaussian Mixture Models',   */
/*				D.A.Reynolds, T.F.Quatieri, R.B.Dunn, 2000 (Academic Press)     */
/*																																				*/
/*		Some Issues:																												*/
/*			- All instances share a common pointer to a background model      */
/*        to speed up testing by caching                                  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.04.2006																								*/
/**************************************************************************/

#include <math.h>
#include <float.h>
#include <assert.h>
#include "SC_Aux.h"
#include "SC_MixtureModel_GMM_UBM.h"
#include <SV_Error.h>

SC_MixtureModel_GMM_UBM::SC_MixtureCache SC_MixtureModel_GMM_UBM::topMixtureCache(0);

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel_GMM_UBM::SC_MixtureModel_GMM_UBM(SC_TweakableParameters* pTweak, SC_MixtureModel_GMM *pUBM) : SC_MixtureModel(pTweak) {
  this->pUBM = new SC_MixtureModel_GMM((SC_MixtureModel_GMM &)*pUBM);

  if (this->pUBM == NULL) {
    REPORT_ERROR(SVLIB_BadArg, "Need a UBM to adapt a new GMM");
  }

	this->mixtureCount = this->pUBM->getMixtureCount(); 
	this->dim = this->pUBM->getDim();
	this->Hdr.ModelType = sclib::mtGMM_UBM;
  this->pTopMixtureList = NULL;
  this->maxEMsteps = 0; //not needed in this class

  if (this->mixtureCount > 0) {
    MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM: sd");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM: mean");
    MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_GMM_UBM: weight");
		for (int x = 0; x < this->mixtureCount; x++) {
      this->weight[x] = this->pUBM->getWeight(x);
      for (int y = 0; y < this->dim; y++) {
        this->variance[x][y] = this->pUBM->getVariance(x, y);
        this->sd[x][y] = this->pUBM->getSd(x, y);
			  this->mean[x][y] = this->pUBM->getMean(x, y);
      }
		}
	}
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_MixtureModel_GMM_UBM::SC_MixtureModel_GMM_UBM(const SC_MixtureModel_GMM_UBM& pParent) : SC_MixtureModel(pParent) {
  this->pUBM = new SC_MixtureModel_GMM((SC_MixtureModel_GMM &)*pParent.pUBM);
  this->pTopMixtureList = NULL;
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_MixtureModel_GMM_UBM::~SC_MixtureModel_GMM_UBM() {
  MFree_2D(this->pTopMixtureList);
  MFree_0D(this->pUBM);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_MixtureModel_GMM_UBM& SC_MixtureModel_GMM_UBM::operator=(const SC_MixtureModel_GMM_UBM& pParent) {
	if (this != &pParent) {
		this->SC_MixtureModel::operator=(pParent);

		MFree_2D(this->pTopMixtureList);
		MFree_0D(this->pUBM);
	  
		this->pUBM = new SC_MixtureModel_GMM((SC_MixtureModel_GMM &)*pParent.pUBM);
	}

  return *this;
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by rose/reynolds.
//
// here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_MixtureModel_GMM_UBM::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// the gmm is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by rose/reynolds.
//
// here, one can specify how many segments of a linked list of feature-vectors should be merged
//====================================================================================================================
int SC_MixtureModel_GMM_UBM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	double Log_p_X = 0.0; //the log-likelihood of the data given the actual/last model: p(X|lambda)
	double *p, Sum_p, *Pr, *n, **Ex, **Exx, adaptationCoeff, **oldMean, weightSum, **scaledVariance;
	unsigned long int t; //0..T over all feature-vectors
	unsigned short int i; //0..M over all mixture-components
	unsigned short int d; //0..D over all components of the feature-vectors
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);

  MArray_1D(p, this->mixtureCount, double, "SC_MixtureModel_GMM_UBM.TrainModel: p");
  MArray_1D(Pr, this->mixtureCount, double, "SC_MixtureModel_GMM_UBM.TrainModel: Pr");
  MArray_1D(n, this->mixtureCount, double, "SC_MixtureModel_GMM_UBM.TrainModel: n");
  MArray_2D(Ex, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM.TrainModel: Ex");
  MArray_2D(Exx, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM.TrainModel: Exx");
  MArray_2D(oldMean, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM.TrainModel: oldMean");
  MArray_2D(scaledVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM.TrainModel: scaledVariance");
	
	//initialize with zeros...
  for (i = 0; i < this->mixtureCount; i++) { 
		n[i] = 0.0;
    for (d = 0; d < this->dim; d++) {
      Ex[i][d] = 0.0;
      Exx[i][d] = 0.0;
    }
	}

  //initialize the scaled variance to speed up gaussian computation
	for (i = 0; i < this->mixtureCount; i++) {
  	for (d = 0; d < this->dim; d++) {
      scaledVariance[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
		}			
  }

	for (t = 0; t < (unsigned long int)pCompleteData->Row; t++) {
		//calculate log-gaussian density of the current feature vector with the complete model (Sum_pb) and each single mixture (pb)
    for (i = 0; i < this->mixtureCount; i++) {
      p[i] = log(this->weight[i]);
			for (d = 0; d < this->dim; d++) { //the multivariate gaussian pdf may be evalutated this way (component-wise), because we have diagonal covariances (=> features are assumed uncorrelated)!
        p[i] += scaledVariance[i][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[i][d]) * (pCompleteData->Mat[t][d] - this->mean[i][d]) / this->variance[i][d]);
      }
      Sum_p = (i == 0) ? p[i] : sclib::sLogAdd(Sum_p, p[i]);
		}//i=0..M
    Log_p_X += Sum_p;

    //calculate the sufficient statistics: count, first and second moment
		for (i = 0; i < this->mixtureCount; i++) {
      Pr[i] = sclib::sExp(p[i] - Sum_p);
      n[i] += Pr[i];
      for (d = 0; d < this->dim; d++) {
        Ex[i][d] += Pr[i] * pCompleteData->Mat[t][d];
        Exx[i][d] += Pr[i] * pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d]; //x^2 = diag(x*x')
			}
		} //i=0..M
	} //t=0..T

  //finish calculateing statistics
	for (i = 0; i < this->mixtureCount; i++) {
    for (d = 0; d < this->dim; d++) {
      Ex[i][d] /= n[i];
      Exx[i][d] /= n[i];
		}
	} //i=0..M

	//let's start with the reestimation process:
  weightSum = 0.0;
	for (i = 0; i < this->mixtureCount; i++) {
    adaptationCoeff = n[i] / (n[i] + this->pTweak->mixtureModelGmmubm.relevanceFactor);
    if (this->pTweak->mixtureModelGmmubm.adaptWeights == true) {
      this->weight[i] = (adaptationCoeff*n[i]/(double)(pCompleteData->Row)) + ((1.0-adaptationCoeff)*this->weight[i]);
      weightSum += this->weight[i];
    }
    for (d = 0; d < this->dim; d++) {
      oldMean[i][d] = this->mean[i][d]; //make a copy that is used during variance-adaptation
			if (this->pTweak->mixtureModelGmmubm.adaptMeans == true) {
        this->mean[i][d] = (adaptationCoeff*Ex[i][d]) + ((1.0-adaptationCoeff)*this->mean[i][d]);
      }
      if (this->pTweak->mixtureModelGmmubm.adaptVariances == true) {
			  this->variance[i][d] = (adaptationCoeff*Exx[i][d]) + ((1.0-adaptationCoeff)*(this->variance[i][d]+(oldMean[i][d]*oldMean[i][d])-(this->mean[i][d]*this->mean[i][d])));
			  if (this->variance[i][d] < this->pTweak->mixtureModelGmmubm.varianceLimit) {this->variance[i][d] = this->pTweak->mixtureModelGmmubm.varianceLimit;} //do variance-limiting
        this->sd[i][d] = sqrt(this->variance[i][d]);
      }
		}			
  }
  if (this->pTweak->mixtureModelGmmubm.adaptWeights == true) {
  	for (i = 0; i < this->mixtureCount; i++) {
      this->weight[i] /= weightSum; //make them adding up to 1.0
    }
  }

	if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("gmmUbm.txt", this, this->pTweak);}
	if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::scalarOut("gmmUbm_likelihood.txt", Log_p_X, this->pTweak);}
  
	this->trainingDataCount = pCompleteData->Row;

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }
	MFree_1D(p);
	MFree_1D(Pr);
	MFree_1D(n);
	MFree_2D(Ex);
	MFree_2D(Exx);
  MFree_2D(oldMean);
  MFree_2D(scaledVariance);
	
	return 1;
}

//====================================================================================================================
// Test GMM, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_GMM_UBM::TestModel(SV_Data *pData) {
	SV_Data	*pScore = TestModel(pData, 1);

	return pScore;					
}

//====================================================================================================================
// Test GMM with a specified nr. of segments out of the linked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_GMM_UBM::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(0); //merge all patterns in this linked list together
	unsigned long int t, T = (unsigned long int)pCompleteData->Row; //0..T over all feature-vectors
	unsigned short int idx, i, I = this->mixtureCount; //0..M over all mixture-components
	unsigned short int d, D = this->dim; //0..D over all components of the feature-vectors
	double ubmLikelihood = 0.0, gmmLikelihood = 0.0, l, sumUBM, sumGMM;
	SV_Data *pScore = new SV_Data(1+pCompleteData->Row, 1);
  double **scaledVarianceUBM = NULL, **scaledVarianceGMM = NULL;

  //initialize the scaled variance of the GMM (needed in any case) to speed up gaussian computation
  MArray_2D(scaledVarianceGMM, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM_UBM.TrainModel: scaledVarianceGMM");
	for (i = 0; i < this->mixtureCount; i++) {
  	for (d = 0; d < this->dim; d++) {
      scaledVarianceGMM[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi));
		}			
  }
  
  if (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM_UBM || (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM_UBM_CACHE && this->topMixtureCache.cacheHit(pCompleteData) == false)) { //just use standard GMM-UBM scoring

    //initialize the scaled variance of the UBM (needed only if cache isn't used ) to speed up gaussian computation
    MArray_2D(scaledVarianceUBM, this->pUBM->getMixtureCount(), this->pUBM->getDim(), double, "SC_MixtureModel_GMM_UBM.TrainModel: scaledVarianceUBM");
    for (i = 0; i < this->mixtureCount; i++) {
  	  for (d = 0; d < this->dim; d++) {
        scaledVarianceUBM[i][d] = log(1.0 / (this->pUBM->getSd(i, d) * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
		  }			
    }

    //init the top-mixture cache
    if (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM_UBM_CACHE) {
      this->topMixtureCache.initTopList(pCompleteData, this->pTweak->mixtureModelGmmubm.topCmixtures);
    }

	  for (t = 0; t < T; t++) {
      //determine the C top scoring mixtures for the UBM
		  for (i = 0; i < I;  i++) { 
        l = log(this->pUBM->getWeight(i));
			  for (d = 0; d < D; d++) {
          l += scaledVarianceUBM[i][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->pUBM->getMean(i, d)) * (pCompleteData->Mat[t][d] - this->pUBM->getMean(i, d)) / this->pUBM->getVariance(i, d));
        }
        insertTopMixture(i, l);
		  }

      //fill the top-mixture cache
      if (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM_UBM_CACHE) {
        this->topMixtureCache.copyTopList(t, this->pTopMixtureList);
      }
      
      //compute the UBM's and this model's scores only from the top scoring UBM mixtures
      for (i = 0; i < this->pTweak->mixtureModelGmmubm.topCmixtures; i++) {
        idx = getTopMixture(i+1, l);
        sumUBM = (i == 0) ? l : sclib::sLogAdd(sumUBM, l); //likelihood of the UBM for only the top-scoring mixtures

        l = log(this->weight[idx]);
			  for (d = 0; d < D; d++) {
          l += scaledVarianceGMM[idx][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[idx][d]) * (pCompleteData->Mat[t][d] - this->mean[idx][d]) / this->variance[idx][d]);
        }
        sumGMM = (i == 0) ? l : sclib::sLogAdd(sumGMM, l); //likelohood of the GMM (this model!) for only the top-scoring mixtures
      }

      ubmLikelihood += sumUBM;
      gmmLikelihood += sumGMM;
			pScore->Mat[t+1][0] = (float)(sumGMM-sumUBM);
	  }

    MFree_2D(scaledVarianceUBM);

  } else if (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM_UBM_CACHE && this->topMixtureCache.cacheHit(pCompleteData) == true) { //use the top-mixture cache

	  for (t = 0; t < T; t++) {
      //compute the UBM's and this model's scores only from the top scoring UBM mixtures
      for (i = 0; i < this->pTweak->mixtureModelGmmubm.topCmixtures; i++) {
        idx = this->topMixtureCache.getTopMixture(t, i+1, l); 
        sumUBM = (i == 0) ? l : sclib::sLogAdd(sumUBM, l); //likelihood of the UBM for only the top-scoring mixtures

        l = log(this->weight[idx]);
			  for (d = 0; d < D; d++) {
          l += scaledVarianceGMM[idx][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[idx][d]) * (pCompleteData->Mat[t][d] - this->mean[idx][d]) / this->variance[idx][d]);
        }
        sumGMM = (i == 0) ? l : sclib::sLogAdd(sumGMM, l); //likelohood of the GMM (this model!) for only the top-scoring mixtures
      }

      ubmLikelihood += sumUBM;
      gmmLikelihood += sumGMM;
			pScore->Mat[t+1][0] = (float)(sumGMM-sumUBM);
	  }

  } else if (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM) { //decide whether to use standard GMM scoring or special (normal case!) and fast GMM-UBM scoring (GMM scoring leaves out normalizing by the UBM's score...)

    //compute the score of this model as if it where just a standard GMM
	  for (t = 0; t < T; t++) {
		  for (i = 0; i < I;  i++) {
        l = log(this->weight[i]);
			  for (d = 0; d < D; d++) {
          l += scaledVarianceGMM[i][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[i][d]) * (pCompleteData->Mat[t][d] - this->mean[i][d]) / this->variance[i][d]);
        }
        sumGMM = (i == 0) ? l : sclib::sLogAdd(sumGMM, l);
		  }
      gmmLikelihood += sumGMM;
			pScore->Mat[t+1][0] = (float)(sumGMM);
	  }

  } else { //scoring method unknown
    REPORT_ERROR(SVLIB_BadArg, "Specified scoring method for a GMM-UBM unknown");
  }
	
  //clean up and store result
  if (TestData != pCompleteData) {
	  MFree_0D(pCompleteData);
  }
  MFree_2D(scaledVarianceGMM);
	pScore->Mat[0][0] = (float)(gmmLikelihood/(double)(T) - ubmLikelihood/(double)(T)); //divide by T seems better, because it fits to the other models... (there according to "Speaker Verification Using Adapted Gaussian Mixture Models")

  return (pScore);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_GMM_UBM::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_MixtureModel*)pSecond, false));
}
SC_MixtureModel_GMM_UBM* SC_MixtureModel_GMM_UBM::combineModels(SC_MixtureModel* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_GMM_UBM::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_GMM_UBM* SC_MixtureModel_GMM_UBM::combineModels(SC_MixtureModel* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext) {
  REPORT_ERROR(SVLIB_BadArg, "Two GMM-UBM's can't be combined by adding up their mixtures, they must be retrained!");

	return NULL; //TODO: any possibility to do this anyway?
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_MixtureModel_GMM_UBM::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pSecond, pSpeechFrames, segmentsToMerge, pBackgroundModels));
}
SC_MixtureModel_GMM_UBM* SC_MixtureModel_GMM_UBM::combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
	SC_MixtureModel_GMM_UBM* pNewModel;

  assert(this->dim == pSecond->getDim());
	
	pNewModel = new SC_MixtureModel_GMM_UBM(this->pTweak, this->pUBM);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

	return pNewModel;
}

//====================================================================================================================
// Add a mixture/likelihood-pair to the list of top-scoring ones; returns true if it was added (if it was 
// top-scoring), otherwise doesn't add and returns false
//====================================================================================================================
bool SC_MixtureModel_GMM_UBM::insertTopMixture(unsigned short int mixtureNum, double likelihood) {
  unsigned long int i, C = this->pTweak->mixtureModelGmmubm.topCmixtures;
  
  if (this->pTopMixtureList == NULL) { //initialize the top-list, if still uninitialized
    MArray_2D(this->pTopMixtureList, (long)(C), 2, double, "SC_MixtureModel_GMM_UBM.insertTopMixture: pTopMixtureList");
    for (i = 0; i < C; i++) {
      this->pTopMixtureList[i][0] = (double)(-1); //-1 indicates "uninitialized" here
      this->pTopMixtureList[i][1] = -1.0 * std::numeric_limits<double>::max();
    }
  }

  if (this->pTopMixtureList[0][1] < likelihood) { //the given mixture is top-scoring if it scores better than the last one in the list of top-scorers (idx=0 is the last one, idx=C-1 is the best one)
    for (i = 1; i < C; i++) { //sort the new element in (in ascending order, so the 0th elemt has the lowest rank/smallest likelihood)
                              //use some sort of incomplete bubbelsort here, because that is optimal in this case (prior knowledge, pre-sorted list)
      if (this->pTopMixtureList[i][1] < likelihood && i < C-1) {
        this->pTopMixtureList[i-1][0] = this->pTopMixtureList[i][0];
        this->pTopMixtureList[i-1][1] = this->pTopMixtureList[i][1];
      } else {
        this->pTopMixtureList[i-1][0] = this->pTopMixtureList[i][0];
        this->pTopMixtureList[i-1][1] = this->pTopMixtureList[i][1];
        this->pTopMixtureList[i][0] = (double)mixtureNum;
        this->pTopMixtureList[i][1] = likelihood;
        break;
      }
    }

    return true;
  } else {
    return false;
  }
}

//====================================================================================================================
// Returns the mixture-number of the top-scoring mixture with given rank (1..C)  or -1 in the case of error
//====================================================================================================================
short int SC_MixtureModel_GMM_UBM::getTopMixture(unsigned short int rank, double &likelihood) {
  unsigned long int C = this->pTweak->mixtureModelGmmubm.topCmixtures;
  unsigned int idx = C - rank;
  
  if (this->pTopMixtureList != NULL && rank <= C) {
    likelihood = this->pTopMixtureList[idx][1];
    return (short int)sclib::round(this->pTopMixtureList[idx][0]);
  } else {
    likelihood = -1.0 * std::numeric_limits<double>::max();
    return -1;
  }
}
