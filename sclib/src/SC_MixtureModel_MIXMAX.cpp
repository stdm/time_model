/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a GMM-IB as described in			  		*/
/*				'Integrated Models of Signal and Background with Application to */
/*				 Speaker Identification in Noise', R.C.Rose, E.M.Hofstetter,		*/
/*			   D.A.Reynolds, 1994 (IEEE)																			*/
/*				 &																															*/
/*				'Speech Recognition Using Noise-Adaptive Prototypes', A.Nadas,	*/
/*				 D.Nahamoo, M.A.Picheny, 1989 (IEEE)														*/
/*																																				*/
/*		Some Issues:																												*/
/*			- The model uses only diaginal Covariance-Matrices, because this	*/
/*				has proven to yield best results in the above paper. To draw		*/
/*				maximum reward from this simplification, we use only a vector		*/
/*				to store the variances. This simplifies some calculations as		*/
/*				well as the memory-usage																				*/
/*			- Variance-limiting is applied.																		*/
/*			- Initialization of the model prior to estimating parameters			*/
/*			  isn't well researched; here, the means are build randomly, and  */
/*				some iterations of k-means algorithm privide somehow			    	*/
/*				meaningfull distribution of the data along the mixtures 				*/
/*				according to the paper.																				  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 06.04.2004																								*/
/**************************************************************************/

#include <math.h>
#include <limits>
#include <assert.h> 
#include <float.h>
#include "SC_MixtureModel_MIXMAX.h"
#include "SC_Aux.h"
#include "SC_Gauss.h"
#include <SV_Error.h>

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel_MIXMAX::SC_MixtureModel_MIXMAX(SC_TweakableParameters* pTweak, SC_MixtureModel* pBackground, unsigned short int mixtureCount, unsigned short int dim) : SC_MixtureModel(pTweak) {
	this->pBackground = pBackground;
  this->pOriginalBackground	= NULL;
	this->mixtureCount = mixtureCount; 
	this->dim = dim;
  this->noiseCorruptionType = (this->pTweak != NULL) ? this->pTweak->mixtureModelMixMax.noiseCorruptionType : sclib::nctMax;
  this->maxEMsteps = this->pTweak->mixtureModelMixMax.maxEMiterations;

	if (this->pBackground != NULL) {
    this->pOriginalBackground	= new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pBackground));
	  assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
	}

	this->Hdr.ModelType = sclib::mtMIXMAX;

	//normally, the mixtureCount must be >0, but there is one exception from this rule:
	//maybe we need only a "link" in the linked list of models, but have no data to model.
	//then, a model without any mixtures represents such a link without modelling anything.
	if (this->mixtureCount > 0) {
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX: sd");
		MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX: mean");
		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_MIXMAX: weight");
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_MIXMAX: maskLevel");

		for (int x = 0; x < this->mixtureCount; x++) {
      for (int y = 0; y < this->dim; y++) {
			  this->variance[x][y] = 0.0;
        this->sd[x][y] = 0.0;
			  this->mean[x][y] = 0.0;
      }
			this->weight[x] = 0.0;
      this->maskLevel[x] = 0.0;
		}
	} else {
		this->variance = NULL;
    this->sd = NULL;
		this->mean = NULL;
		this->weight =	NULL;
    this->maskLevel = NULL;
	}
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_MixtureModel_MIXMAX::SC_MixtureModel_MIXMAX(const SC_MixtureModel_MIXMAX& pParent) : SC_MixtureModel(pParent) {
	this->pBackground = pParent.pBackground;
  this->pOriginalBackground	= (pParent.pOriginalBackground != NULL) ? new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pParent.pOriginalBackground)) : NULL;
  this->noiseCorruptionType = pParent.noiseCorruptionType;

	if (this->pBackground != NULL) {
    assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
	}

	if (this->mixtureCount > 0) {
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_GMM: maskLevel");
		for (int x = 0; x < this->mixtureCount; x++) {
      this->maskLevel[x] = pParent.maskLevel[x];
		}
	}
}

//====================================================================================================================
// destructor
// attention: doesn't delete the background-gmm, nore the linked gmm-ib's!
// but: does delete the original background gmm, because it's a copy, not a reference/pointer!
//====================================================================================================================
SC_MixtureModel_MIXMAX::~SC_MixtureModel_MIXMAX() {
  this->pBackground = NULL;

  MFree_0D(this->pOriginalBackground);
  MFree_1D(this->maskLevel);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_MixtureModel_MIXMAX& SC_MixtureModel_MIXMAX::operator=(const SC_MixtureModel_MIXMAX& pParent) {
	if (this != &pParent) {
		this->SC_MixtureModel::operator=(pParent);

		MFree_1D(this->maskLevel);
		MFree_0D(this->pOriginalBackground);
	  
		this->pBackground = pParent.pBackground;
		this->pOriginalBackground	= (pParent.pOriginalBackground != NULL) ? new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pParent.pOriginalBackground)) : NULL;
		this->noiseCorruptionType = pParent.noiseCorruptionType;

		if (this->pBackground != NULL) {
			assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
		}

		if (this->mixtureCount > 0) {
			MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_GMM: maskLevel");
			for (int x = 0; x < this->mixtureCount; x++) {
				this->maskLevel[x] = pParent.maskLevel[x];
			}
		}  
	}
  
  return *this;
}

//====================================================================================================================
// Inherited version of training-algorithm; just calls the new one (below), which is capable of handling many differnt
// noise-models for different parts of the feature-vectors
//====================================================================================================================
int SC_MixtureModel_MIXMAX::TrainModel(SV_Data* TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Trains the model, using segmentsToMerge segments from the linked list of features; evaluates the 
// noiseCorruptionType and chooses the suitable algorithm
//====================================================================================================================
int SC_MixtureModel_MIXMAX::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  int res;

  switch (this->noiseCorruptionType) {
    case sclib::nctAdditive: 
    case sclib::nctAdditiveLog:
      REPORT_ERROR(SVLIB_BadData, "Additive signal-noise-interaction is not maintained!");
      //res = TrainWithAdditiveNoise(pData, segmentsToMerge);
      res = -1;
      break;
    case sclib::nctMax:
      res = TrainWithMaxNoise(pData, segmentsToMerge);
      break;
    case sclib::nctMaxLog:
      res = TrainWithMaxNoise_LogArithmetic(pData, segmentsToMerge);
      break;
    default:
      REPORT_ERROR(SVLIB_BadData, "No such signal-noise-interaction defined!");
      res = -1;
  }

  return res;
}

//====================================================================================================================
// training with max signal-noise-relationship
//
// the MIXMAX is first initialized by a function and then the model parameters are estimated to best fit the data in a
// maximum likelihood sense. this is done by the em-algorithm for censored data as described by rose/reynolds.
//
// here's an addition to rose/reynolds algorithm:  This version can handle different background-models for different 
// parts of the given feature-vectors (if they and the background-models are organized in a linked list, for every new 
// list-element of the data, the next new list-element of the background-models is used). The parameter 
// "segmentsToMerge" specifies how many list-elements should be considered (1: standart-case, only first list-element; 
// 0: all list-elemts in the list)
//====================================================================================================================
int SC_MixtureModel_MIXMAX::TrainWithMaxNoise(SV_Data *pData, unsigned long int segmentsToMerge) {
  double **a, **A, **b, **B, bB; //gaussian and cumulative gaussian distributions of background (a) and signal (b)
  double **E_xlz, **E_xlz2; //E{x^n|x<z,i,t,lambda}, n={1,2}
  double p_xez; //p(x=z|i,j,lambda)
  double p_z; //p(z|i,j,lambda) for one component of the vector z
  double p_ij; //p(i,j|z,lambda)
  double ***E_x, ***E_x2; //E{x^n|z,i,j,lambda), n={1,2} for all i,j,d
  double **SumSum_pij_Ex, **SumSum_pij_Ex2; //SUM_t=1..T(SUM_j=1..N(p(i,j|z,lambda))) [* E{x|...}, *E{x^2|...}] for all i,d
  double *SumSumProd_pij; //SUM_t=1..T(SUM_j=1..N(PROD_d=1..D(p(i,j|z,lambda)))) for the whole vector z and all i
  double Prod_pz;
  double SumSumProd_pz_pq; //SUM_i=1..M(SUM_j=1..N(PROD_d=1..D(p(z|i,j,lambda))*p*q) for the whole vector instead of for each component
  double **Prod_pz_pq; //PROD_d=1..D(p(z|i,j,lambda) for all i,j
  double log_pZ = 0.0, log_pZ_old; //the log-likelihood of the complete observation given the model
  short i, j, d, I = this->mixtureCount, J = sclib::max(1, this->pBackground->getMaxMixturesInList()), D = this->dim;
  unsigned int iterationCount = 0; //count of loops of the em-algorithm
  double scaledX;
  unsigned long int idx;

  long actualT	= pData->Row; //count of feature-vectors in the actual elements of the linked list
	SV_Data	*pActualData = pData;	//pointer to actual element of linked list
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	SC_MixtureModel *pBackground	= this->pBackground; //don't touch the original pointer during training!
  long t, T = pCompleteData->Row; //actual number and complete count of feature vectors

  MArray_2D(a, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: a");
  MArray_2D(A, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: A");
  MArray_2D(b, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: b");
  MArray_2D(B, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: B");
  MArray_2D(E_xlz, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: E_xlz");
  MArray_2D(E_xlz2, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: E_xlz2");
  MArray_3D(E_x, I, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: E_x");
  MArray_3D(E_x2, I, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: E_x2");
  MArray_2D(SumSum_pij_Ex, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: SumSum_pij_Ex");
  MArray_2D(SumSum_pij_Ex2, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: SumSum_pij_Ex2");
  MArray_1D(SumSumProd_pij, I, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: SumSumProd_pij");
  MArray_2D(Prod_pz_pq, I, J, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise: Prod_pz_pq");

	//The model-parameters must be initialized some way prior to the em-algorithm
  initParameters(pCompleteData, this->pTweak->mixtureModelMixMax.varianceLimit, this->pTweak->mixtureModelMixMax.kMeansIterations);
	if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIXMAX.txt", this, this->pTweak);}

  do { //start with the em-algorithm for integrated signal-background-models 
       //with z = max(x, y) as the signal-noise-interaction function

    iterationCount++;
    log_pZ_old = log_pZ;
    log_pZ = 0.0;

    for (i = 0; i < I; i++) {
      for (d = 0; d < D; d++) {
        SumSum_pij_Ex[i][d] = 0.0;
        SumSum_pij_Ex2[i][d] = 0.0;
      }
      SumSumProd_pij[i] = 0.0;
      this->maskLevel[i] = 0.0;
    }

    actualT	= pData->Row;
	  pActualData = pData;
	  pBackground	= this->pBackground;
    J = sclib::max(1, pBackground->getMixtureCount());

    for (t = 0; t < T; t++) {
      
			//handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
		  if (t >= actualT) {
        assert(pActualData->Next != NULL);
        assert(pBackground->Next != NULL);
			  pActualData	= pActualData->Next;
			  actualT += pActualData->Row;
			  pBackground	= (SC_MixtureModel*)(pBackground->Next);
        J = sclib::max(1, pBackground->getMixtureCount());
		  }

      for (d = 0; d < D; d++) {
        for (i = 0; i < I; i++) { 
          scaledX = sclib::zTransform(pCompleteData->Mat[t][d], this->mean[i][d], this->sd[i][d]);
          idx = this->gaussSolver.getIdx(scaledX);
          this->gaussSolver.tabledGaussianAndErf(idx, this->sd[i][d], b[i][d], B[i][d]);
          bB = b[i][d] / B[i][d];
          if (sclib::isFinite(bB) == true) { //the standart case as reported by nadas et. al.
            E_xlz[i][d] = this->mean[i][d] - this->variance[i][d] * bB;
            E_xlz2[i][d] = (this->mean[i][d] *  this->mean[i][d] + this->variance[i][d]) - this->variance[i][d] * bB * (pCompleteData->Mat[t][d] + this->mean[i][d]);
          } else { //if b(z)/B(z)~0/0, then this expressions can be solved (using the rule of de l'hospital) to b'(z)/B'(z)=(mu-z)/sigma^2, which yields the following expression for the conditional expectatgions E{x^n|...}:
					  E_xlz[i][d] = pCompleteData->Mat[t][d]; //z;
					  E_xlz2[i][d] =	this->variance[i][d] + (pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d]); //s2 + (z * z);
          }
        }

        for (j = 0; j < J; j++) { //J >= 1 !!!
          scaledX = sclib::zTransform(pCompleteData->Mat[t][d], pBackground->getMean(j, d), pBackground->getSd(j, d));
          idx = this->gaussSolver.getIdx(scaledX);
          this->gaussSolver.tabledGaussianAndErf(idx, pBackground->getSd(j, d), a[j][d], A[j][d]);
        }
      }

      SumSumProd_pz_pq = 0.0;

      //1. do some preliminary caching so that this algorithm becomes efficient
      for (i = 0; i < I; i++) {
        for (j = 0; j < J; j++) {
          Prod_pz = 1.0;
          Prod_pz_pq[i][j] = this->weight[i] * pBackground->getWeight(j);

          for (d = 0; d < D; d++) {            
            p_z = (a[j][d] * B[i][d]) + (b[i][d] * A[j][d]);
            if (p_z > 0) {
              p_xez = (b[i][d] * A[j][d]) / p_z;
              E_x[i][j][d] = p_xez * pCompleteData->Mat[t][d] + (1.0 - p_xez) * E_xlz[i][d];
              E_x2[i][j][d] = p_xez * pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d] + (1.0 - p_xez) * E_xlz2[i][d];
              this->maskLevel[i] += 1.0 - p_xez;
            } else {
              E_x[i][j][d] = 0.0;  //the conditional expectations are not really =0 in the case that p(z|..)=0, 
							E_x2[i][j][d] = 0.0; //but because of p(z|...)=0 they don't play a role in the later summation, so they are set =0 for simplification purposes.
              this->maskLevel[i] += 1.0;
            }
            Prod_pz *= p_z;
          } //d
          
          Prod_pz_pq[i][j] *= Prod_pz;
          SumSumProd_pz_pq += Prod_pz_pq[i][j];
        } //j
      } //i

			//.2 calculate the sums needed to reestimate model-parameters
      for (i = 0; i < I; i++) {
        for (j = 0; j < J; j++) {
          //assert(SumSumProd_pz_pq > 0.0); //TODO: this and next line...
          p_ij = (SumSumProd_pz_pq != 0.0) ? Prod_pz_pq[i][j] / SumSumProd_pz_pq : 0.0; //HISTORICAL BIG BUG: p_ij = p_z_pq[i][j][d] / SumSum_pz_pq[d];
          SumSumProd_pij[i] += p_ij;
          for (d = 0; d < D; d++) {
            SumSum_pij_Ex[i][d] += p_ij * E_x[i][j][d];
            SumSum_pij_Ex2[i][d] += p_ij * E_x2[i][j][d];
          } //d
        } //j
      } //i

      //assert(SumSumProd_pz_pq > 0.0); //TODO
      log_pZ += sclib::sLog(SumSumProd_pz_pq);
    } //t

    //3. finally start with the reestimation process:
    for (i = 0; i < I; i++) {
      this->weight[i] = SumSumProd_pij[i] / (double)(T);
      //assert(this->weight[i] > 0.0);
      //assert(SumSum_pij[i][d] != 0.0);
      if (this->weight[i] > 0.0) {
        for (d = 0; d < D; d++) {
          this->mean[i][d] = SumSum_pij_Ex[i][d] / SumSumProd_pij[i];
          this->variance[i][d] = (SumSum_pij_Ex2[i][d] / SumSumProd_pij[i]) - (this->mean[i][d] * this->mean[i][d]);
          if (this->variance[i][d] < this->pTweak->mixtureModelMixMax.varianceLimit) {this->variance[i][d] = this->pTweak->mixtureModelMixMax.varianceLimit;} //do variance limiting
          this->sd[i][d] = sqrt(this->variance[i][d]);
        } //d
        this->maskLevel[i] /= (double)(T*D*J);
      } else {
        killMixture(i);
        I--;
        printf("\n     mixture %i killed in run %i, %i remaining", i, iterationCount, this->mixtureCount);
      }      
    } //i

		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIXMAX.txt", this, this->pTweak);}
		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::scalarOut("MIXMAX_likelihood.txt", log_pZ, this->pTweak);}

    //assert(log_pZ >= log_pZ_old || iterationCount == 1);

    if (iterationCount == this->maxEMsteps) {
      break;
    }

  } while (fabs(log_pZ_old - log_pZ) > this->pTweak->mixtureModelMixMax.EMthreshold);

  //kill too small mixtures
  this->sortParameters();
  for (i = 0; i < this->mixtureCount; i++) {
    if (this->weight[i] < this->pTweak->mixtureModelMixMax.weightLimit) { //kill too small mixtures
      this->killMixture(i);
      i--;
    }
  }

	this->trainingDataCount = pCompleteData->Row;
  MFree_0D(this->pOriginalBackground);
  this->pOriginalBackground =  new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)this->pBackground));

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }
  MFree_2D(a);
  MFree_2D(A);
  MFree_2D(b);
  MFree_2D(B);
  MFree_2D(E_xlz);
  MFree_2D(E_xlz2);
  MFree_3D(E_x);
  MFree_3D(E_x2);
  MFree_2D(Prod_pz_pq);
  MFree_1D(SumSumProd_pij);
  MFree_2D(SumSum_pij_Ex);
  MFree_2D(SumSum_pij_Ex2);
  
  return iterationCount;
}

//====================================================================================================================
// training with max signal-noise-relationship as above, but with logarithmic arithmetic
//====================================================================================================================
int SC_MixtureModel_MIXMAX::TrainWithMaxNoise_LogArithmetic(SV_Data *pData, unsigned long int segmentsToMerge) {
  double **a, **A, **b, **B, bB; //gaussian and cumulative gaussian distributions of background (a) and signal (b)
  double **E_xlz, **E_xlz2; //E{x^n|x<z,i,t,lambda}, n={1,2}
  double p_xez; //p(x=z|i,j,lambda)
  double p_z; //p(z|i,j,lambda) for one component of the vector z
  double p_ij; //p(i,j|z,lambda)
  double ***E_x, ***E_x2; //E{x^n|z,i,j,lambda), n={1,2} for all i,j,d
  double **SumSum_pij_Ex, **SumSum_pij_Ex2; //SUM_t=1..T(SUM_j=1..N(p(i,j|z,lambda))) [* E{x|...}, *E{x^2|...}] for all i,d
  double *SumSumProd_pij; //SUM_t=1..T(SUM_j=1..N(PROD_d=1..D(p(i,j|z,lambda)))) for the whole vector z and all i
  double Prod_pz;
  double SumSumProd_pz_pq; //SUM_i=1..M(SUM_j=1..N(PROD_d=1..D(p(z|i,j,lambda))*p*q) for the whole vector instead of for each component
  double **Prod_pz_pq; //PROD_d=1..D(p(z|i,j,lambda) for all i,j
  double log_pZ = 0.0, log_pZ_old; //the log-likelihood of the complete observation given the model
  short i, j, d, I = this->mixtureCount, J = sclib::max(1, this->pBackground->getMaxMixturesInList()), D = this->dim;
  unsigned int iterationCount = 0; //count of loop of the em-algorithm
  double scaledX, **logSd, **logBgSd;
  unsigned long int idx;

  long actualT	= pData->Row; //count of feature-vectors in the actual elements of the linked list
	SV_Data	*pActualData = pData;	//pointer to actual element of linked list
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	SC_MixtureModel *pBackground	= this->pBackground; //don't touch the original pointer during training!
  long t, T = pCompleteData->Row; //actual number and complete count of feature vectors

  MArray_2D(a, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: a");
  MArray_2D(A, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: A");
  MArray_2D(b, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: b");
  MArray_2D(B, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: B");
  MArray_2D(E_xlz, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: E_xlz");
  MArray_2D(E_xlz2, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: E_xlz2");
  MArray_3D(E_x, I, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: E_x");
  MArray_3D(E_x2, I, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: E_x2");
  MArray_2D(SumSum_pij_Ex, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: SumSum_pij_Ex");
  MArray_2D(SumSum_pij_Ex2, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: SumSum_pij_Ex2");
  MArray_1D(SumSumProd_pij, I, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: SumSumProd_pij");
  MArray_2D(Prod_pz_pq, I, J, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: Prod_pz_pq");
  MArray_2D(logSd, I, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: logSd");
  MArray_2D(logBgSd, J, D, double, "SC_MixtureModel_MIXMAX.TrainWithMaxNoise_LogArithmetic: logBgSd");

	//The model-parameters must be initialized some way prior to the em-algorithm
  initParameters(pCompleteData, this->pTweak->mixtureModelMixMax.varianceLimit, this->pTweak->mixtureModelMixMax.kMeansIterations);
	if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIXMAX.txt", this, this->pTweak);}
  for (d = 0; d < D; d++) {
    for (i = 0; i < I; i++) {
      logSd[i][d] = sclib::sLog(this->sd[i][d]);
    }
    for (j = 0; j < J; j++) {
      logBgSd[j][d] = sclib::sLog(pBackground->getSd(j, d));
    }
  }

  do { //start with the em-algorithm for integrated signal-background-models 
       //with z = max(x, y) as the signal-noise-interaction function

    iterationCount++;
    log_pZ_old = log_pZ;
    log_pZ = 0.0;

    for (i = 0; i < I; i++) {
      for (d = 0; d < D; d++) {
        SumSum_pij_Ex[i][d] = 0.0;
        SumSum_pij_Ex2[i][d] = 0.0;
      }
      SumSumProd_pij[i] = 0.0;
      this->maskLevel[i] = 0.0;
    }

    actualT	= pData->Row;
	  pActualData = pData;
	  pBackground	= this->pBackground;
    J = sclib::max(1, pBackground->getMixtureCount());

    for (t = 0; t < T; t++) {
      
			//handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
		  if (t >= actualT) {
        assert(pActualData->Next != NULL);
        assert(pBackground->Next != NULL);
			  pActualData	= pActualData->Next;
			  actualT += pActualData->Row;
			  pBackground	= (SC_MixtureModel*)(pBackground->Next);
        J = sclib::max(1, pBackground->getMixtureCount());
        for (j = 0; j < J; j++) {
          for (d = 0; d < D; d++) {
            logBgSd[j][d] = sclib::sLog(pBackground->getSd(j, d));
          }
        }
		  }

      for (d = 0; d < D; d++) {
        for (i = 0; i < I; i++) { 
          scaledX = sclib::zTransform(pCompleteData->Mat[t][d], this->mean[i][d], this->sd[i][d]);
          idx = this->gaussSolver.getIdx(scaledX);
          this->gaussSolver.tabledLogGaussianAndLogErf(idx, logSd[i][d], b[i][d], B[i][d]);
          bB = sclib::sExp(b[i][d] - B[i][d]);
          if (sclib::isFinite(bB) == true) { //the standart case as reported by nadas et. al.
            E_xlz[i][d] = this->mean[i][d] - this->variance[i][d] * bB;
            E_xlz2[i][d] = (this->mean[i][d] *  this->mean[i][d] + this->variance[i][d]) - this->variance[i][d] * bB * (pCompleteData->Mat[t][d] + this->mean[i][d]);
          } else { //if b(z)/B(z)~0/0, then this expressions can be solved (using the rule of de l'hospital) to b'(z)/B'(z)=(mu-z)/sigma^2, which yields the following expression for the conditional expectatgions E{x^n|...}:
					  E_xlz[i][d] = pCompleteData->Mat[t][d]; //z;
					  E_xlz2[i][d] =	this->variance[i][d] + (pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d]); //s2 + (z * z);
          }
        }

        for (j = 0; j < J; j++) { //J >= 1 !!!
          scaledX = sclib::zTransform(pCompleteData->Mat[t][d], pBackground->getMean(j, d), pBackground->getSd(j, d));
          idx = this->gaussSolver.getIdx(scaledX);
          this->gaussSolver.tabledLogGaussianAndLogErf(idx, logBgSd[j][d], a[j][d], A[j][d]);
        }
      }

      SumSumProd_pz_pq = 0.0;

      //1. do some preliminary caching so that this algorithm becomes efficient
      for (i = 0; i < I; i++) {
        for (j = 0; j < J; j++) {
          Prod_pz = 0.0;
          Prod_pz_pq[i][j] = sclib::sLog(this->weight[i] * pBackground->getWeight(j));

          for (d = 0; d < D; d++) {            
            p_z = sclib::sLogAdd(a[j][d] + B[i][d], b[i][d] + A[j][d]);
            p_xez = exp(b[i][d] + A[j][d] - p_z);
            if (sclib::isFinite(p_xez) == true) {
              E_x[i][j][d] = p_xez * pCompleteData->Mat[t][d] + (1.0 - p_xez) * E_xlz[i][d];
              E_x2[i][j][d] = p_xez * pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d] + (1.0 - p_xez) * E_xlz2[i][d];
              this->maskLevel[i] += 1.0 - p_xez;
            } else {
              E_x[i][j][d] = 0.0;  //the conditional expectations are not really =0 in the case that p(z|..)=0, 
							E_x2[i][j][d] = 0.0; //but because of p(z|...)=0 they don't play a role in the later summation, so they are set =0 for simplification purposes.
              this->maskLevel[i] += 1.0;
            }
            Prod_pz += p_z;
          } //d
          
          Prod_pz_pq[i][j] += Prod_pz;
          SumSumProd_pz_pq = (i == 0 && j == 0) ? Prod_pz_pq[i][j] : sclib::sLogAdd(SumSumProd_pz_pq, Prod_pz_pq[i][j]);
        } //j
      } //i

			//.2 calculate the sums needed to reestimate model-parameters
      for (i = 0; i < I; i++) {
        for (j = 0; j < J; j++) {
          p_ij = sclib::sExp(Prod_pz_pq[i][j] - SumSumProd_pz_pq); //HISTORICAL BIG BUG: p_ij = p_z_pq[i][j][d] / SumSum_pz_pq[d];
          SumSumProd_pij[i] += p_ij;
          for (d = 0; d < D; d++) {
            SumSum_pij_Ex[i][d] += p_ij * E_x[i][j][d];
            SumSum_pij_Ex2[i][d] += p_ij * E_x2[i][j][d];
          } //d
        } //j
      } //i

      log_pZ += SumSumProd_pz_pq;
    } //t

    //3. finally start with the reestimation process:
    for (i = 0; i < I; i++) {
      this->weight[i] = SumSumProd_pij[i] / (double)(T);
      //assert(this->weight[i] > 0.0);
      //assert(SumSum_pij[i][d] != 0.0);
      if (this->weight[i] > 0.0) {
        for (d = 0; d < D; d++) {
          this->mean[i][d] = SumSum_pij_Ex[i][d] / SumSumProd_pij[i];
          this->variance[i][d] = (SumSum_pij_Ex2[i][d] / SumSumProd_pij[i]) - (this->mean[i][d] * this->mean[i][d]);
          if (this->variance[i][d] < this->pTweak->mixtureModelMixMax.varianceLimit) {this->variance[i][d] = this->pTweak->mixtureModelMixMax.varianceLimit;} //do variance limiting
          this->sd[i][d] = sqrt(this->variance[i][d]);
          logSd[i][d] = sclib::sLog(this->sd[i][d]);
        } //d
        this->maskLevel[i] /= (double)(T*D*J);
      } else {
        killMixture(i);
        printf("\n     mixture %i killed in run %i, %i remaining", i, iterationCount, this->mixtureCount);
      }      
    } //i

		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIXMAX.txt", this, this->pTweak);}
		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::scalarOut("MIXMAX_likelihood.txt", log_pZ, this->pTweak);}

    //assert(log_pZ >= log_pZ_old || iterationCount == 1);

    if (iterationCount == this->maxEMsteps) {
      break;
    }

  } while (fabs(log_pZ_old - log_pZ) > this->pTweak->mixtureModelMixMax.EMthreshold);

  //kill too small mixtures
  this->sortParameters();
  for (i = 0; i < this->mixtureCount; i++) {
    if (this->weight[i] < this->pTweak->mixtureModelMixMax.weightLimit) { //kill too small mixtures
      this->killMixture(i);
      i--;
    }
  }

	this->trainingDataCount = pCompleteData->Row;
  MFree_0D(this->pOriginalBackground);
  this->pOriginalBackground =  new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)this->pBackground));

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }
  MFree_2D(a);
  MFree_2D(A);
  MFree_2D(b);
  MFree_2D(B);
  MFree_2D(E_xlz);
  MFree_2D(E_xlz2);
  MFree_3D(E_x);
  MFree_3D(E_x2);
  MFree_2D(Prod_pz_pq);
  MFree_1D(SumSumProd_pij);
  MFree_2D(SumSum_pij_Ex);
  MFree_2D(SumSum_pij_Ex2);
  MFree_2D(logSd);
  MFree_2D(logBgSd);
  
  return iterationCount;
}

//====================================================================================================================
// Test GMM-IB while considering only one element of the linked list 
// return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_MIXMAX::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;					
}

//====================================================================================================================
// Test GMM-IB with a specified nr. of segments out of the linked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_MIXMAX::TestModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  SV_Data *pScore = NULL;

  switch (this->noiseCorruptionType) {
    case sclib::nctAdditive:
    case sclib::nctAdditiveLog:
      REPORT_ERROR(SVLIB_BadArg, "MIXMAX with additive noise is not maintained, please check the code!");
      //pScore = TestWithAdditiveNoise(pData, segmentsToMerge);
      break;
    case sclib::nctMax:
      pScore = TestWithMaxNoise(pData, segmentsToMerge);
      break;
    case sclib::nctMaxLog:
      pScore = TestWithMaxNoise_LogArithmetic(pData, segmentsToMerge);
      break;
    default:
      REPORT_ERROR(SVLIB_BadData, "No such signal-noise-interaction defined!");
  }

  return pScore;
}

//====================================================================================================================
// Test under the assumption that z=max(x,y)
// Here, no logarithmic arithmetic is used because it's slower than normal arithmetic
//====================================================================================================================
SV_Data* SC_MixtureModel_MIXMAX::TestWithMaxNoise(SV_Data *pData, unsigned long int segmentsToMerge) {
  double Log_p_Z = 0.0; //the log-likelihood of the test-data given the model
	SV_Data *pScore	= new SV_Data(1, 1); //the return-value
	double **a, **A; //gaussian and cumulative gaussian distributions of background for all j,d
	double **b, **B; //gaussian and cumulative gaussian distributions of speech-signal for all d
  double scaledX;
  unsigned long int idx;

  unsigned short int	i, j, d, J, I = this->mixtureCount, D = this->dim;
	double Prod_pz_pq; //PROD_d=1..D(p(z|i,j,lambda)
	double SumSumProd_pz_pq; //SUM_i=1..M(SUM_j=1..N(PROD_d=1..D(p(z|i,j,lambda))*p*q) for the whole vector instead of for each component

	long int actualT = pData->Row; //count of feature-vectors in the actual elements of the linked list
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);	//merge selected list-elemts in a linked list together
	SV_Data *pActualData = pData; //pointer to actual element of linked list
  SC_MixtureModel *pBackground;
  if (this->pTweak->mixtureModelMixMax.bgModelCombination == true) {
    pBackground = (SC_MixtureModel*)(this->pBackground->combineModels(this->pBackground, this->pOriginalBackground, true)); //don't touch the original pointer during training!
  } else {
    pBackground = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)(this->pBackground)));
  }
  SC_MixtureModel *pBgHook = NULL;
  long int t, T = pCompleteData->Row;

  if (this->pTweak->mixtureModelMixMax.bgModelCombination == true) {
    J = (this->pBackground->getMaxMixturesInList() == 0  && this->pOriginalBackground->getMixtureCount() == 0) ? 1 : sclib::max(1, this->pBackground->getMaxMixturesInList()) + sclib::max(1, this->pOriginalBackground->getMixtureCount());
  } else {
    J =  sclib::max(1, pBackground->getMaxMixturesInList());
  }

  MArray_2D(a, J, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise: a");
  MArray_2D(A, J, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise: A");
  MArray_2D(b, I, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise: b");
  MArray_2D(B, I, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise: B");
  
  J = sclib::max(pBackground->getMixtureCount(), 1);
	for (t = 0; t < T; t++) {
		if (t >= actualT) { //handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
      assert(pActualData->Next != NULL);
      assert(pBackground->Next != NULL);
			pActualData	= pActualData->Next;
			actualT += pActualData->Row;
      pBgHook = (SC_MixtureModel*)(pBackground->Next);
      MFree_0D(pBackground);
      if (this->pTweak->mixtureModelMixMax.bgModelCombination == true) {
        pBackground = (SC_MixtureModel*)(this->pBackground->combineModels((SC_MixtureModel*)(pBgHook), this->pOriginalBackground, true));
      } else {
        pBackground = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)(pBgHook)));
      }
      J = sclib::max(1, pBackground->getMixtureCount());
		}

    for (d = 0; d < D; d++) {
      for (i = 0; i < I; i++) { 
        scaledX = sclib::zTransform(pCompleteData->Mat[t][d], this->mean[i][d], this->sd[i][d]);
        idx = this->gaussSolver.getIdx(scaledX);
        this->gaussSolver.tabledGaussianAndErf(idx, this->sd[i][d], b[i][d], B[i][d]);
      }
      for (j = 0; j < J; j++) { //J >= 1 !!!
        scaledX = sclib::zTransform(pCompleteData->Mat[t][d], pBackground->getMean(j, d), pBackground->getSd(j, d));
        idx = this->gaussSolver.getIdx(scaledX);
        this->gaussSolver.tabledGaussianAndErf(idx, pBackground->getSd(j, d), a[j][d], A[j][d]);
      }
    }

    SumSumProd_pz_pq = 0.0;
    for (i = 0; i < I; i++) { 
	    for (j = 0; j < J; j++) { //J >= 1 !!!
  	    Prod_pz_pq	= this->weight[i] * pBackground->getWeight(j);
		    for (d = 0; d < D; d++) { 
			    Prod_pz_pq *= (a[j][d] * B[i][d]) + (b[i][d] * A[j][d]);
		    }   
  	    SumSumProd_pz_pq += Prod_pz_pq;
      }
    }

    //accumulate results
    Log_p_Z += sclib::sLog(SumSumProd_pz_pq);
  } //t=0..T

  MFree_0D(pBackground);
  if (pData != pCompleteData) {
	  MFree_0D(pCompleteData);
  }
  MFree_2D(a);
  MFree_2D(A);
  MFree_2D(b);
  MFree_2D(B);

	pScore->Mat[0][0] = (float)(Log_p_Z) / (float)(T); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...
	return (pScore);
}

//====================================================================================================================
// Test under the assumption that z=max(x,y)
// Here, logarithmic arithmetic is used
//====================================================================================================================
SV_Data* SC_MixtureModel_MIXMAX::TestWithMaxNoise_LogArithmetic(SV_Data *pData, unsigned long int segmentsToMerge) {
  double Log_p_Z = 0.0; //the log-likelihood of the test-data given the model
	SV_Data *pScore	= new SV_Data(1, 1); //the return-value
	double **a, **A; //gaussian and cumulative gaussian distributions of background for all j,d
	double **b, **B; //gaussian and cumulative gaussian distributions of speech-signal for all d
  double scaledX, **logSd, **logBgSd;
  unsigned long int idx;

  unsigned short int	i, j, d, J, I = this->mixtureCount, D = this->dim;
	double Prod_pz_pq; //PROD_d=1..D(p(z|i,j,lambda)
	double SumSumProd_pz_pq; //SUM_i=1..M(SUM_j=1..N(PROD_d=1..D(p(z|i,j,lambda))*p*q) for the whole vector instead of for each component

	long int actualT = pData->Row; //count of feature-vectors in the actual elements of the linked list
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);	//merge selected list-elemts in a linked list together
	SV_Data *pActualData = pData; //pointer to actual element of linked list
  SC_MixtureModel *pBackground;
  if (this->pTweak->mixtureModelMixMax.bgModelCombination == true) {
    pBackground	= (SC_MixtureModel*)(this->pBackground->combineModels(this->pBackground, this->pOriginalBackground, true)); //don't touch the original pointer during training!
  } else {
    pBackground = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)(this->pBackground)));
  }
  SC_MixtureModel *pBgHook = NULL;
  long int t, T = pCompleteData->Row;

  if (this->pTweak->mixtureModelMixMax.bgModelCombination == true) {
    J = (this->pBackground->getMaxMixturesInList() == 0  && this->pOriginalBackground->getMixtureCount() == 0) ? 1 : sclib::max(1, this->pBackground->getMaxMixturesInList()) + sclib::max(1, this->pOriginalBackground->getMixtureCount());
  } else {
    J =  sclib::max(1, pBackground->getMaxMixturesInList());
  }

  MArray_2D(a, J, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise_LogArithmetic: a");
  MArray_2D(A, J, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise_LogArithmetic: A");
  MArray_2D(b, I, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise_LogArithmetic: b");
  MArray_2D(B, I, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise_LogArithmetic: B");
  MArray_2D(logSd, I, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise_LogArithmetic: logSd");
  MArray_2D(logBgSd, J, D, double, "SC_MixtureModel_MIXMAX.TestWithMaxNoise_LogArithmetic: logBgSd");
  
  //cache the log_sd
  for (d = 0; d < D; d++) {
    for (i = 0; i < I; i++) {
      logSd[i][d] = sclib::sLog(this->sd[i][d]);
    }
    for (j = 0; j < J; j++) {
      logBgSd[j][d] = sclib::sLog(pBackground->getSd(j, d));
    }
  }

  //main loop
  J = sclib::max(pBackground->getMixtureCount(), 1);
	for (t = 0; t < T; t++) {
		if (t >= actualT) { //handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
      assert(pActualData->Next != NULL);
      assert(pBackground->Next != NULL);
			pActualData	= pActualData->Next;
			actualT += pActualData->Row;
      pBgHook = (SC_MixtureModel*)(pBackground->Next);
      MFree_0D(pBackground);
      if (this->pTweak->mixtureModelMixMax.bgModelCombination == true) {
        pBackground = (SC_MixtureModel*)(this->pBackground->combineModels((SC_MixtureModel*)(pBgHook), this->pOriginalBackground, true));
      } else {
        pBackground = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)(pBgHook)));
      }
      J = sclib::max(1, pBackground->getMixtureCount());
      for (j = 0; j < J; j++) {
        for (d = 0; d < D; d++) {
          logBgSd[j][d] = sclib::sLog(pBackground->getSd(j, d));
        }
      }
		}

    for (d = 0; d < D; d++) {
      for (i = 0; i < I; i++) { 
        scaledX = sclib::zTransform(pCompleteData->Mat[t][d], this->mean[i][d], this->sd[i][d]);
        idx = this->gaussSolver.getIdx(scaledX);
        this->gaussSolver.tabledLogGaussianAndLogErf(idx, logSd[i][d], b[i][d], B[i][d]);
      }
      for (j = 0; j < J; j++) { //J >= 1 !!!
        scaledX = sclib::zTransform(pCompleteData->Mat[t][d], pBackground->getMean(j, d), pBackground->getSd(j, d));
        idx = this->gaussSolver.getIdx(scaledX);
        this->gaussSolver.tabledGaussianAndErf(idx, logBgSd[j][d], a[j][d], A[j][d]);
      }
    }

    SumSumProd_pz_pq = 0.0;
    for (i = 0; i < I; i++) { 
	    for (j = 0; j < J; j++) { //J >= 1 !!!
  	    Prod_pz_pq	= sclib::sLog(this->weight[i] * pBackground->getWeight(j));
		    for (d = 0; d < D; d++) { 
			    Prod_pz_pq += sclib::sLogAdd(a[j][d] + B[i][d], b[i][d] + A[j][d]);
		    }   
        SumSumProd_pz_pq = (i == 0 && j == 0) ? Prod_pz_pq : sclib::sLogAdd(SumSumProd_pz_pq, Prod_pz_pq);
      }
    }

    //accumulate results
    Log_p_Z += SumSumProd_pz_pq;
  } //t=0..T

  MFree_0D(pBackground);
  if (pData != pCompleteData) {
	  MFree_0D(pCompleteData);
  }
  MFree_2D(a);
  MFree_2D(A);
  MFree_2D(b);
  MFree_2D(B);
  MFree_2D(logSd);
  MFree_2D(logBgSd);

	pScore->Mat[0][0] = (float)(Log_p_Z) / (float)(T); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...;
	return (pScore);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_MIXMAX::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_MixtureModel*)pSecond, false));
}
SC_MixtureModel_MIXMAX* SC_MixtureModel_MIXMAX::combineModels(SC_MixtureModel* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model* SC_MixtureModel_MIXMAX::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel_MIXMAX*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_MIXMAX* SC_MixtureModel_MIXMAX::combineModels(SC_MixtureModel_MIXMAX* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext) {
  unsigned short int i, d, newMixtureCount = pFirst->getMixtureCount() + pSecond->getMixtureCount();
	unsigned long int newTrainingDataCount = pFirst->getTrainingDataCount() + pSecond->getTrainingDataCount();
	double *newWeight, *newMaskLevel, **newMean, **newVariance, **newSd;
	SC_MixtureModel_MIXMAX *pNewModel;
  SC_MixtureModel *pNewOriginalBackground;
  bool tooSmall = false;
  
	if (pFirst->getDim() != pSecond->getDim()) {REPORT_ERROR(SVLIB_BadData, "Can't combine Models trained on different feature-dim!");}

	MArray_1D(newWeight, newMixtureCount, double, "SC_MixtureModel_MIXMAX.combineModels: newWeight");
  MArray_1D(newMaskLevel, newMixtureCount, double, "SC_MixtureModel_MIXMAX.combineModels: newMaskLevel");
	MArray_2D(newMean, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIXMAX.combineModels: newMean");
	MArray_2D(newVariance, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIXMAX.combineModels: newVariance");
  MArray_2D(newSd, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIXMAX.combineModels: newSd");

  pNewModel = new SC_MixtureModel_MIXMAX(pFirst->getTweak(), NULL, newMixtureCount, pFirst->getDim());

	//copy the mixture-components together
	for (i = 0; i < newMixtureCount; i++) {
		if (i < pFirst->getMixtureCount()) {
			newWeight[i] = (pFirst->getWeight(i) * pFirst->getTrainingDataCount()) / newTrainingDataCount; //re-normalize the weights so that their sum = 1.0
      newMaskLevel[i] = pFirst->getMaskLevel(i);
      for (d = 0; d < pFirst->getDim(); d++) {
				newMean[i][d] = pFirst->getMean(i, d);
				newVariance[i][d] = pFirst->getVariance(i, d);
        newSd[i][d] = pFirst->getSd(i, d);
			}
		} else {
			newWeight[i] = (pSecond->getWeight(i - pFirst->getMixtureCount()) * pSecond->getTrainingDataCount()) / newTrainingDataCount; //re-normalize the weights so that their sum = 1.0
      if (pSecond->Hdr.ModelType == sclib::mtMIXMAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX_ex) {
        newMaskLevel[i] = ((SC_MixtureModel_MIXMAX*)pSecond)->getMaskLevel(i - pFirst->getMixtureCount());
      } else {
        newMaskLevel[i] = 0.0; //no gmm-ib => no noise-masking
      }
			for (d = 0; d < pFirst->getDim(); d++) {
				newMean[i][d] = pSecond->getMean(i - pFirst->getMixtureCount(), d);
				newVariance[i][d] = pSecond->getVariance(i - pFirst->getMixtureCount(), d);
        newSd[i][d] = pSecond->getSd(i - pFirst->getMixtureCount(), d);
			}
		}
    if (newWeight[i] < pFirst->getTweak()->mixtureModelMixMax.weightLimit) {tooSmall = true;}
	}

	//make the new mixture-components the new parameters of the new model
  pNewModel->setWeight(newWeight);
  pNewModel->setMaskLevel(newMaskLevel);
  pNewModel->setMean(newMean);
  pNewModel->setVariance(newVariance);
  pNewModel->setSd(newSd);
  pNewModel->setTrainindDataCount(newTrainingDataCount);
  pNewModel->setMixtureCount(newMixtureCount);
  if (pSecond->Hdr.ModelType == sclib::mtMIXMAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX_ex) {
    pNewOriginalBackground = (SC_MixtureModel*)(pFirst->getOriginalBackground()->combineModels(((SC_MixtureModel_MIXMAX*)pSecond)->getOriginalBackground()));
    pNewModel->setOriginalBackground(pNewOriginalBackground);
  } else {
    SC_MixtureModel_GMM *pNewDummy = new SC_MixtureModel_GMM(pFirst->getTweak(), 0, pFirst->getDim());
    pNewOriginalBackground = (SC_MixtureModel*)(pNewDummy->combineModels(pFirst->getOriginalBackground()));
    pNewModel->setOriginalBackground(pNewOriginalBackground);
    MFree_0D(pNewDummy);
  }

  //kill too small mixtures
  if (tooSmall == true) {
    pNewModel->sortParameters();
    for (i = 0; i < pNewModel->getMixtureCount(); i++) {
      if (pNewModel->getWeight(i) < pFirst->getTweak()->mixtureModelMixMax.weightLimit) {
        pNewModel->killMixture(i);
        i--;
      }
    }
  }

  if (keepFirstsNext == true) {
    pNewModel->Next = pFirst->Next;
  } else {
    pNewModel->Next = NULL;
  }

	return pNewModel;
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model; 
// Convert a dummy mixture (no mean, no weight, no variance) to a real mixture with equal result
//====================================================================================================================
SC_Model* SC_MixtureModel_MIXMAX::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepDummyMixture, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel_MIXMAX*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_MIXMAX* SC_MixtureModel_MIXMAX::combineModels(SC_MixtureModel_MIXMAX* pFirst, SC_MixtureModel* pSecond, bool keepDummyMixture, bool keepFirstsNext) {
  assert("Don't use this method!!!" == 0); //TODO: why???
  
  return combineModels(pFirst, pSecond, keepFirstsNext);
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_MixtureModel_MIXMAX::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pSecond, pSpeechFrames, segmentsToMerge, pBackgroundModels));
}
SC_MixtureModel_MIXMAX* SC_MixtureModel_MIXMAX::combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_MixtureModel* pBackgroundModels) {
	SC_MixtureModel_MIXMAX* pNewModel;
	
  assert("Don't use this method!!!" == 0); //TODO: why???
  assert(this->dim == pSecond->getDim());
	
	pNewModel = new SC_MixtureModel_MIXMAX(this->pTweak, pBackgroundModels, sclib::min(this->mixtureCount+pSecond->getMixtureCount(), this->pTweak->modelHandler.maxSpeakerModelOrder), this->dim);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

	return pNewModel;
}

//====================================================================================================================
// Sort the model-parameters according to the weights
//====================================================================================================================
void SC_MixtureModel_MIXMAX::sortParameters() {
	unsigned short int i, d;
	SV_Data *temp;
	double *oldWeight, *oldMaskLevel, **oldMean, **oldVariance, **oldSd;
	
	temp = new SV_Data(this->mixtureCount, 2);
	MArray_1D(oldWeight, this->mixtureCount, double, "SC_MixtureModel_MIXMAX.sortParameters: oldWeight");
	MArray_1D(oldMaskLevel, this->mixtureCount, double, "SC_MixtureModel_MIXMAX.sortParameters: oldMaskLevel");
	MArray_2D(oldMean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX.sortParameters: oldMean");
	MArray_2D(oldVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX.sortParameters: oldVariance");
  MArray_2D(oldSd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX.sortParameters: oldVariance");

	for (i = 0; i < this->mixtureCount; i++) {
		temp->Mat[i][0] = (float)(this->weight[i]);
		temp->Mat[i][1] = (float)(i); //remember the old index in the 2nd column of the data-matrix
		oldWeight[i] = this->weight[i];
    oldMaskLevel[i] = this->maskLevel[i];
		for (d = 0; d < this->dim; d++) {
			oldMean[i][d] = this->mean[i][d];
			oldVariance[i][d] = this->variance[i][d];
      oldSd[i][d] = this->sd[i][d];
		}
	}

	//sort the list of weights; remember the old index in the 2nd column of the data-matrix
	//sclib::quickSort(temp, 0, this->mixtureCount - 1, 2, 0);
  sclib::quickSort(temp->Mat, 0, this->mixtureCount - 1, 2, 0);

	for (i = 0; i < this->mixtureCount; i++) {
		this->weight[i] = oldWeight[(int)(temp->Mat[i][1])];
    this->maskLevel[i] = oldMaskLevel[(int)(temp->Mat[i][1])];
		for (d = 0; d < this->dim; d++) {
			this->mean[i][d] = oldMean[(int)(temp->Mat[i][1])][d];
			this->variance[i][d] = oldVariance[(int)(temp->Mat[i][1])][d];
      this->sd[i][d] = oldSd[(int)(temp->Mat[i][1])][d];
		}
	}

	MFree_0D(temp);
	MFree_1D(oldWeight);
  MFree_1D(oldMaskLevel);
	MFree_2D(oldMean);
	MFree_2D(oldVariance);
  MFree_2D(oldSd);

	return;
}

//====================================================================================================================
// Remove the indexed mixture-component from the models mean/variance/weight-vectors/matrices
//====================================================================================================================
void SC_MixtureModel_MIXMAX::killMixture(unsigned short int mixture) {
  unsigned short int i, i_idx = 0, d,  I = this->mixtureCount, D = this->dim;
  double *newWeight, *newMaskLevel, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I-1, double, "SC_MixtureModel_MIXMAX.killMixture: newWeight");
  MArray_1D(newMaskLevel, I-1, double, "SC_MixtureModel_MIXMAX.killMixture: newMaskLevel");
  MArray_2D(newVariance, I-1, D, double, "SC_MixtureModel_MIXMAX.killMixture: newVariance");
  MArray_2D(newMean, I-1, D, double, "SC_MixtureModel_MIXMAX.killMixture: newMean");
  MArray_2D(newSd, I-1, D, double, "SC_MixtureModel_MIXMAX.killMixture: newSd");

  for (i = 0; i < I; i++) { //copy mixtures, leave the one to kill out
    if (i != mixture) {
      newWeight[i_idx] = this->weight[i] * (1.0  / (1.0 - this->weight[mixture])); //renormalize weight
      newMaskLevel[i_idx] = this->maskLevel[i];
      for (d = 0; d < D; d++) {
        newVariance[i_idx][d] = this->variance[i][d];
        newMean[i_idx][d] = this->mean[i][d];
        newSd[i_idx][d] = this->sd[i][d];
      }
      i_idx++;
    }
  }

  MFree_1D(this->weight);
  MFree_1D(this->maskLevel);
  MFree_2D(this->variance);
  MFree_2D(this->mean);
  MFree_2D(this->sd);

  this->weight = newWeight;
  this->maskLevel = newMaskLevel;
  this->variance = newVariance;
  this->mean = newMean;
  this->sd = newSd;
  this->mixtureCount--;

  return;
}

//====================================================================================================================
// Add one new mixture component with all means=0, variances=varianceLimit, weight=1/mixtureCount
//====================================================================================================================
void SC_MixtureModel_MIXMAX::addMixture(void) {
  unsigned short int i, d, I = this->mixtureCount, D = this->dim;
  double *newWeight, *newMaskLevel, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I+1, double, "SC_MixtureModel_MIXMAX.addMixture: newWeight");
  MArray_1D(newMaskLevel, I+1, double, "SC_MixtureModel_MIXMAX.addMixture: newMaskLevel");
  MArray_2D(newVariance, I+1, D, double, "SC_MixtureModel_MIXMAX.addMixture: newVariance");
  MArray_2D(newMean, I+1, D, double, "SC_MixtureModel_MIXMAX.addMixture: newMean");
  MArray_2D(newSd, I+1, D, double, "SC_MixtureModel_MIXMAX.addMixture: newSd");

  //add new mixture
  newWeight[I] = 1.0 / (this->mixtureCount+1.0);
  newMaskLevel[I] = 0.0;
  for (d = 0; d < D; d++) {
    newVariance[I][d] = 1.0;
    newMean[I][d] = 0.0;
    newSd[I][d] = 1.0;
  }

  //recalc old mixtures
  for (i = 0; i < I; i++) { 
    newWeight[i] = this->weight[i] * (1.0 - newWeight[I]);
    for (d = 0; d < D; d++) {
      newVariance[i][d] = this->variance[i][d];
      newMean[i][d] = this->mean[i][d];
      newSd[i][d] = this->sd[i][d];
    }
  }

  MFree_1D(this->weight);
  MFree_1D(this->maskLevel);
  MFree_2D(this->variance);
  MFree_2D(this->mean);
  MFree_2D(this->sd);

  this->weight = newWeight;
  this->maskLevel = newMaskLevel;
  this->variance = newVariance;
  this->mean = newMean;
  this->sd = newSd;
  this->mixtureCount++;

  return;
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_MixtureModel_MIXMAX::SaveModel(void) {
 	int res, x, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX model Failed!");}

  //write mixtureCount and featureDim
	bytes = io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX Model Failed!");}

  //write weight-vector
	bytes += io.writeArray(&(this->DFile), this->weight, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX Model Failed!");}

  //write maskLevel-vector
	bytes += io.writeArray(&(this->DFile), this->maskLevel, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX Model Failed!");}

	//write mean and variance matrices
	for (x = 0; x < this->mixtureCount; x++) {
		bytes += io.writeArray(&(this->DFile), this->mean[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->variance[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->sd[x], this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX Model Failed!");}
	}

  //write trainingDataCount
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX Model Failed!");}

  //write noiseCorruptionType
	bytes += io.writeScalar(&(this->DFile), this->noiseCorruptionType);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIXMAX Model Failed!");}

  //write maxEMstep
	io.writeScalar(&(this->DFile), this->maxEMsteps);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}

  return bytes + MHLen; //MHLen is just an estimate of the written header-bytes...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model* SC_MixtureModel_MIXMAX::LoadModel(void) {
	int res, x;
	unsigned short int newDim = 0, newMixtureCount = 0;
  SC_MixtureModel_GMM *pDummy = NULL;
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
	if ((DFile.good() != TRUE) || (newMixtureCount == 0)) {return(NULL);}

	//prepare this model for changing it's parameters
	if ((newDim != this->dim) || (newMixtureCount != this->mixtureCount)) {
		MFree_1D(this->weight);
    MFree_1D(this->maskLevel);
		MFree_2D(this->mean);
		MFree_2D(this->variance);
    MFree_2D(this->sd);

		this->dim = newDim;
		this->mixtureCount = newMixtureCount;

 		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_MIXMAX: weight");
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_MIXMAX: maskLevel");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX: mean");
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIXMAX: sd");
	};

	//read weight, maskLevel, mean and variance
	io.readArray(&(this->DFile), this->weight, this->mixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}
	io.readArray(&(this->DFile), this->maskLevel, this->mixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}
	for (x = 0; x < this->mixtureCount; x++) {
		io.readArray(&(this->DFile), this->mean[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->variance[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->sd[x], this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}
	}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}

  //Set a dummy as the original background
  pDummy = new SC_MixtureModel_GMM(this->pTweak, 0, this->dim);
  this->pOriginalBackground = pDummy;

	//read maxEMstep
	io.readScalar(&(this->DFile), this->noiseCorruptionType, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}

  //read maxEMstep
	io.readScalar(&(this->DFile), this->maxEMsteps, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIXMAX Failed!");}

	return(this);
}

//====================================================================================================================
// create a (linked, not copied) SC_Signature-view on this model (for distance computation )
//====================================================================================================================
SC_Signature* SC_MixtureModel_MIXMAX::toSignature(void) {
  double *weights;
  SC_Centroid **centroids = NULL;
  SC_Signature *pSignature;
  
  MArray_1D(centroids, this->mixtureCount, SC_Centroid*, "SC_MixtureModel_MIXMAX.toSignature: centroids");
  MArray_1D(weights, this->mixtureCount, double, "SC_MixtureModel_MIXMAX.toSignature: weights");

  for (unsigned short int i = 0; i < this->mixtureCount; i++) {
    weights[i] = this->weight[i] * (1.0 - this->maskLevel[i]);
		centroids[i] = new SC_Centroid_Gaussian(this->pTweak, this->dim, this->mean[i], this->variance[i]);
  }

	pSignature = new SC_Signature(centroids, weights, this->mixtureCount);

  return pSignature;
}

//====================================================================================================================
// destruct a signature created by this model's toSignature()-method
//====================================================================================================================
void SC_MixtureModel_MIXMAX::killSignature(SC_Signature *pSignature) {
  if (pSignature != NULL) {
    MFree_0D(pSignature);
  }

  return;
}

//====================================================================================================================
//
//
//
// A T T E N T I O N :  Below here is untested, unrecommended and unmaintained code!!! 
//                      It is old and NOT up-to-date
//                      Don't use it without debugging/completing it first
//
//
//
//====================================================================================================================

/*
//====================================================================================================================
// training with additive signal-noise-relationship
//====================================================================================================================
int SC_MixtureModel_MIXMAX::TrainWithAdditiveNoise(SV_Data *pData, unsigned long int segmentsToMerge) {
  double p_z; //p(z|i,j,lambda) for one component of the vector z
  double p_ij; //p(i,j|z,lambda)
  double ***E_x, ***E_x2; //E{x^n|z,i,j,lambda), n={1,2} for all i,j,d
  double **SumSum_pij, **SumSum_pij_Ex, **SumSum_pij_Ex2; //SUM_t=1..T(SUM_j=1..N(p(i,j|z,lambda))) [* E{x|...}, *E{x^2|...}] for all i,d
  double *SumSumProd_pij; //SUM_t=1..T(SUM_j=1..N(PROD_d=1..D(p(i,j|z,lambda)))) for the whole vector z and all i
  double Prod_pz;
  double SumSumProd_pz_pq; //SUM_i=1..M(SUM_j=1..N(PROD_d=1..D(p(z|i,j,lambda))*p*q) for the whole vector instead of for each component
  double **Prod_pz_pq; //PROD_d=1..D(p(z|i,j,lambda) for all i,j
  double log_pZ = 0.0, log_pZ_old; //the log-likelihood of the complete observation given the model
  double z, m, s2, p, mb, s2b, q; //actual featire vector / mean, variance and weight of a specific signal/background mixture
  short i, j, d, I = this->mixtureCount, J = this->pBackground->getMaxMixturesInList(), D = this->dim;
  unsigned int iterationCount = 0; //count of loop of the em-algorithm

  long actualT	= pData->Row; //count of feature-vectors in the actual elements of the linked list
	SV_Data	*pActualData = pData;	//pointer to actual element of linked list
	SV_Data	*pCompleteData = pData->MergeData(segmentsToMerge);	//merge all feature-vectors in a linked list together
	SC_MixtureModel *pBackground	= this->pBackground; //don't touch the original pointer during training!
  long t, T = pCompleteData->Row; //actual number and complete count of feature vectors

  assert("Don't use this method!!!" == 0);

  MArray_3D(E_x, I, J, D, double, "SC_MixtureModel_MIXMAX.TrainModel: E_x");
  MArray_3D(E_x2, I, J, D, double, "SC_MixtureModel_MIXMAX.TrainModel: E_x2");
  MArray_2D(SumSum_pij, I, D, double, "SC_MixtureModel_MIXMAX.TrainModel: SumSum_pij");
  MArray_2D(SumSum_pij_Ex, I, D, double, "SC_MixtureModel_MIXMAX.TrainModel: SumSum_pij_Ex");
  MArray_2D(SumSum_pij_Ex2, I, D, double, "SC_MixtureModel_MIXMAX.TrainModel: SumSum_pij_Ex2");
  MArray_1D(SumSumProd_pij, I, double, "SC_MixtureModel_MIXMAX.TrainModel: SumSumProd_pij");
  MArray_2D(Prod_pz_pq, I, J, double, "SC_MixtureModel_MIXMAX.TrainModel: Prod_pz_pq");

	//The model-parameters must be initialized some way prior to the em-algorithm
  initParameters(pCompleteData, 0.01, 10);
	if (this->pTweak->debugMode & SCLIB_DB_MIXMAX) {sclib::classOut("MIXMAX.txt", this, this->pTweak);}
  
  do { //start with the em-algorithm for integrated signal-background-models 
       //with z = max(x, y) as the signal-noise-interaction function

    iterationCount++;
    log_pZ_old = log_pZ;
    log_pZ = 0.0;

    for (i = 0; i < I; i++) {
      for (d = 0; d < D; d++) {
        SumSum_pij_Ex[i][d] = 0.0;
        SumSum_pij[i][d] = 0.0;
        SumSum_pij_Ex2[i][d] = 0.0;
      }
      SumSumProd_pij[i] = 0.0;
    }

    actualT	= pData->Row;
	  pActualData = pData;
	  pBackground	= this->pBackground;

    for (t = 0; t < T; t++) {
      
			//handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
		  if (t >= actualT) {
        assert(pActualData->Next != NULL);
        assert(pBackground->Next != NULL);
			  pActualData	= pActualData->Next;
			  actualT += pActualData->Row;
			  pBackground	= (SC_MixtureModel*)(pBackground->Next);
        J = sclib::max(pBackground->getMixtureCount(), 1);
		  }

      SumSumProd_pz_pq = 0.0;

      //1. do some preliminary caching so that this algorithm becomes efficient
      for (i = 0; i < I; i++) {
        p = this->weight[i];

        for (j = 0; j < J; j++) {
          q = sclib::max(pBackground->getWeight(j), 1.0);  //if we need a virtual mixture due to that the actual background-model is just a dummy without any components
          Prod_pz = 1.0;
          Prod_pz_pq[i][j] = 1.0;

          for (d = 0; d < D; d++) {
            z = pCompleteData->Mat[t][d];
            m = this->mean[i][d];
            s2 = this->variance[i][d];
            mb = pBackground->getMean(j, d);
            s2b = pBackground->getVariance(j, d);
            
            p_z = this->func.gaussian(z, mb, s2b);
            
            if (p_z > 0) {
              p_z *= (sqrt(s2b) / sqrt(s2 + s2b)) * exp(((m*m) / (2*s2)) + (pow((m*s2b - mb*s2 + s2*z), 2) / (2*s2b*s2*(s2+s2b))));
              E_x[i][j][d] = (s2 / (s2+s2b)) * (z + ((s2b/s2)*m - mb));
              E_x2[i][j][d] = ((s2*s2b) / (s2+s2b)) * E_x[i][j][d] * E_x[i][j][d];
            } else {
              E_x[i][j][d] = 0.0;  //the conditional expectations are not really =0 in the case that p(z|..)=0, 
							E_x2[i][j][d] = 0.0; //but because of p(z|...)=0 they don't play a role in the later summation, so they are set =0 for simplification purposes.
            }

            Prod_pz *= p_z;
          } //d
          Prod_pz_pq[i][j] = Prod_pz * p * q;
          SumSumProd_pz_pq += Prod_pz_pq[i][j];
        } //j
      } //i

			//.2 calculate the sums needed to reestimate model-parameters
      for (i = 0; i < I; i++) {
        for (j = 0; j < J; j++) {
          p_ij = Prod_pz_pq[i][j] / SumSumProd_pz_pq; //HISTORICAL BIG BUG: p_ij = p_z_pq[i][j][d] / SumSum_pz_pq[d];
          SumSumProd_pij[i] += p_ij;
          for (d = 0; d < D; d++) {
            SumSum_pij[i][d] += p_ij;
            SumSum_pij_Ex[i][d] += p_ij * E_x[i][j][d];
            SumSum_pij_Ex2[i][d] += p_ij * E_x2[i][j][d];
          } //d
        } //j
      } //i

      assert(SumSumProd_pz_pq > 0);
      log_pZ += log(SumSumProd_pz_pq);
    } //t

    //3. finally start with the reestimation process:
    for (i = 0; i < I; i++) {
      this->weight[i] = SumSumProd_pij[i] / T;
      assert(this->weight[i] > 0.0);
      assert(SumSum_pij[i][d] != 0.0);
      for (d = 0; d < D; d++) {
        this->mean[i][d] = SumSum_pij_Ex[i][d] / SumSum_pij[i][d];
        this->variance[i][d] = (SumSum_pij_Ex2[i][d] / SumSum_pij[i][d]) - (this->mean[i][d] * this->mean[i][d]);
        if (this->variance[i][d] < this->pTweak->varianceLimit) {this->variance[i][d] = this->pTweak->varianceLimit;} //do variance limiting
      } //d
    } //i

 		if (this->pTweak->debugMode & SCLIB_DB_MIXMAX) {sclib::classOut("MIXMAX.txt", this, this->pTweak);}
		if (this->pTweak->debugMode & SCLIB_DB_MIXMAX) {sclib::scalarOut("MIXMAX_likelihood.txt", log_pZ, this->pTweak);}

    assert(log_pZ >= log_pZ_old || iterationCount == 1);

    if (iterationCount == this->maxEMsteps) {
      break;
    }

  } while (fabs(log_pZ_old - log_pZ) > this->pTweak->EMthreshold);

	this->trainingDataCount = pCompleteData->Row;

  //kill too small mixtures
  this->sortParameters();
  for (i = 0; i < this->mixtureCount; i++) {
    if (this->weight[i] < this->pTweak->weightLimit) { //kill too small mixtures
      this->killMixture(i);
      i--;
    }
  }

  MFree_3D(E_x);
  MFree_3D(E_x2);
  MFree_2D(Prod_pz_pq);
  MFree_1D(SumSumProd_pij);
  MFree_2D(SumSum_pij);
  MFree_2D(SumSum_pij_Ex);
  MFree_2D(SumSum_pij_Ex2);
  
  return iterationCount;
}

//====================================================================================================================
// Test under the assumption that z=x+y
//====================================================================================================================
SV_Data* SC_MixtureModel_MIXMAX::TestWithAdditiveNoise(SV_Data *pData, unsigned long int segmentsToMerge) {
  double              Log_p_Z = 0.0; //the log-likelihood of the test-data given the model
	SV_Data							*pScore	= new SV_Data(1, 1); //the return-value

  unsigned short int	i, j, d, I = this->mixtureCount, J = this->pBackground->getMaxMixturesInList(), N = this->pBackground->getMixtureCount(), D = this->dim;
  double              z, m, s2, mb, s2b;
  double							p_z; //p(z|i,j,lambda) for one component of the vector z
	double						  Prod_pz_pq; //PROD_d=1..D(p(z|i,j,lambda)
	double							SumSumProd_pz_pq; //SUM_i=1..M(SUM_j=1..N(PROD_d=1..D(p(z|i,j,lambda))*p*q) for the whole vector instead of for each component

	long int		        actualT				 = pData->Row; //count of feature-vectors in the actual elements of the linked list
	SV_Data							*pCompleteData = pData->MergeData(segmentsToMerge);	//merge selected list-elemts in a linked list together
	SV_Data							*pActualData	 = pData; //pointer to actual element of linked list
	SC_MixtureModel    				*pBackground	 = this->pBackground; //don't touch the original pointer during training!
  long int            t, T = pCompleteData->Row;

  assert("Don't use this method!!!" == 0);

  J = N;
	for (t = 0; t < T; t++) {
		//handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
		if (t >= actualT) {
      assert(pActualData->Next != NULL);
      assert(pBackground->Next != NULL);
			pActualData	= pActualData->Next;
			actualT		 += pActualData->Row;
			pBackground = (SC_MixtureModel*)(pBackground->Next);
      N = pBackground->getMixtureCount();
      J = sclib::max(N, 1);
		}

    SumSumProd_pz_pq = 0.0;
    for (i = 0; i < I; i++) { //these terms are only dependant on i, so they only need to be recalculated once for every t
	    for (j = 0; j < J; j++) {
  	    Prod_pz_pq	= this->weight[i] * ((N > 0) ? pBackground->getWeight(j) : 1.0);
		    for (d = 0; d < D; d++) { //these terms are only dependent on j, so they must only be recalculated if j changed
          z = pCompleteData->Mat[t][d];
          m = this->mean[i][d];
          s2 = this->variance[i][d];
          mb = pBackground->getMean(j, d);
          s2b = pBackground->getVariance(j, d);

          //ERROR: what about the virtual mixture in case of a dummy model?!? -> build it in case you want to use this function!!!
          
          p_z = this->func.gaussian(z, mb, s2b);
          if (p_z > 0) {
            p_z *= (sqrt(s2b) / sqrt(s2 + s2b)) * exp(((m*m) / (2*s2)) + (pow((m*s2b - mb*s2 + s2*z), 2) / (2*s2b*s2*(s2+s2b))));
          }
 			    Prod_pz_pq	*= p_z;
		    }   
  	    SumSumProd_pz_pq	+= Prod_pz_pq;
      }
    }

    //accumulate results
    assert(SumSumProd_pz_pq > 0.0);
    if (SumSumProd_pz_pq <= 0.0) {
      Log_p_Z -= std::numeric_limits<double>::max();
    } else {
      Log_p_Z += log(SumSumProd_pz_pq);
    }
  } //t=0..T

	MFree_0D(pCompleteData);

  pScore->Mat[0][0] = (float)(Log_p_Z) / (float)(T); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...
	return (pScore);
}
*/
