/**************************************************************************/
/*    Derived from:																												*/
/*      -                 																								*/
/*																																				*/
/*    Responsibility:																											*/
/*      - implementations of P.M.Baggenstoss' enhanced EM-Agorithm				*/
/*			- based on his own Matlab-implementation                          */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.02.2005																								*/
/**************************************************************************/

//#include <math.h>
#include <assert.h>
#include <time.h>
#include "SC_BaggenstossEM.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_BaggenstossEM::SC_BaggenstossEM(SC_TweakableParameters *pTweak, bool verbose) {
  this->pTweak = pTweak;
  this->verbose = verbose;
  this->minStd = NULL;
  this->weight = NULL;
  this->mean = NULL;
  this->choleskyCovar = NULL;
  this->maxMixtureCount = 100;
  this->pMatrixFunc = new SC_MatrixFunctions();
	this->U = NULL;
	this->V = NULL;
	this->S = NULL;
	this->logDet = NULL;
	this->transposedInvertedCholeskyCovar = NULL;
}

//====================================================================================================================
//	copy-constructor
//====================================================================================================================
SC_BaggenstossEM::SC_BaggenstossEM(const SC_BaggenstossEM& pParent) {
	this->pTweak = pParent.pTweak;
  this->pMatrixFunc = new SC_MatrixFunctions();
	this->verbose = pParent.verbose;
	this->dim = pParent.dim;
	this->maxMixtureCount = pParent.maxMixtureCount;
	this->mixtureCount = pParent.mixtureCount;
	this->trainingDataCount = pParent.trainingDataCount;

	if (pParent.minStd != NULL) {
		this->minStd = this->pMatrixFunc->copy(pParent.minStd, this->dim);
	} else {
		this->minStd = NULL;
	}
	this->mean = this->pMatrixFunc->copy(pParent.mean, this->mixtureCount, this->dim);
	this->weight = this->pMatrixFunc->copy(pParent.weight, this->dim);
	if (pParent.choleskyCovar != NULL) {
		MArray_1D(this->choleskyCovar, this->mixtureCount, double**, "SC_BaggenstossEM.SC_BaggenstossEM: this->choleskyCovar");
	} else {
		this->choleskyCovar = NULL;
	}
	if (pParent.U != NULL) {
		MArray_1D(this->U, this->mixtureCount, double**, "SC_BaggenstossEM.SC_BaggenstossEM: this->U");
	} else {
		this->U = NULL;
	}
	if (pParent.V != NULL) {
		MArray_1D(this->V, this->mixtureCount, double**, "SC_BaggenstossEM.SC_BaggenstossEM: this->V");
	} else {
		this->V = NULL;
	}
	if (pParent.S != NULL) {
		MArray_1D(this->S, this->mixtureCount, double**, "SC_BaggenstossEM.SC_BaggenstossEM: this->S");
	} else {
		this->S = NULL;
	}
	if (pParent.transposedInvertedCholeskyCovar != NULL) {
		MArray_1D(this->transposedInvertedCholeskyCovar, this->mixtureCount, double**, "SC_BaggenstossEM.SC_BaggenstossEM: this->transposedInvertedCholeskyCovar");
	} else {
		this->transposedInvertedCholeskyCovar = NULL;
	}
	if (pParent.logDet != NULL) {
		MArray_1D(this->logDet, this->mixtureCount, double, "SC_BaggenstossEM.SC_BaggenstossEM: this->logDet");
	} else {
		this->logDet = NULL;
	}
	for (int i = 0; i < this->mixtureCount; i++) {
		if (this->choleskyCovar != NULL) {
			this->choleskyCovar[i] = this->pMatrixFunc->copy(pParent.choleskyCovar[i], this->dim, this->dim);
		}
		if (this->U != NULL) {
			this->U[i] = this->pMatrixFunc->copy(pParent.U[i], this->dim, this->dim);
		}
		if (this->V != NULL) {
			this->V[i] = this->pMatrixFunc->copy(pParent.V[i], this->dim, this->dim);
		}
		if (this->S != NULL) {
			this->S[i] = this->pMatrixFunc->copy(pParent.S[i], this->dim, this->dim);
		}
		if (this->transposedInvertedCholeskyCovar != NULL) {
			this->transposedInvertedCholeskyCovar[i] = this->pMatrixFunc->copy(pParent.transposedInvertedCholeskyCovar[i], this->dim, this->dim);
		}
		if (this->logDet != NULL) {
			this->logDet[i] = pParent.logDet[i];
		}
	}
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_BaggenstossEM::~SC_BaggenstossEM() {
  MFree_0D(this->pMatrixFunc);
  MFree_1D(this->minStd);
  MFree_1D(this->weight);
  MFree_2D(this->mean);
	for (int i = 0; i < this->mixtureCount; i++) {
		if (this->choleskyCovar != NULL) {
			MFree_2D(this->choleskyCovar[i]);
		}
		if (this->U != NULL) {
			MFree_2D(this->U[i]);
		}
		if (this->V != NULL) {
			MFree_2D(this->V[i]);
		}
		if (this->S != NULL) {
			MFree_2D(this->S[i]);
		}
		if (this->transposedInvertedCholeskyCovar != NULL) {
			MFree_2D(this->transposedInvertedCholeskyCovar[i]);
		}
	}
	MFree_1D(this->choleskyCovar);
	MFree_1D(this->U);
	MFree_1D(this->V);
	MFree_1D(this->S);
	MFree_1D(this->transposedInvertedCholeskyCovar);
	MFree_1D(this->logDet);
}

//====================================================================================================================
//	C++ translation of P.M.Baggenstoss' function parm = init_gmix(data,n_modes,names,min_std,[rand_init]) for
//  Matlab (03.04.2001). All Matlab-Like functions are located in SC_MatrixFunctions
//
//  Changes to original version:
//    - minVariance can also be a scalar, not a vector, so that it can be equal for all dimensions
//====================================================================================================================
void SC_BaggenstossEM::init(SV_Data *pData, unsigned short int mixtureCount, double minStd, bool randInit) {
  double* minimumStd;
  MArray_1D(minimumStd, pData->Col, double, "SC_BaggenstossEM.init: minimumStd");

  for (short int d = 0; d < pData->Col; d++) {
    minimumStd[d] = minStd;
  }

  init(pData, mixtureCount, minimumStd, randInit);

  MFree_1D(minimumStd);
  return;
}

void SC_BaggenstossEM::init(SV_Data *pData, unsigned short int mixtureCount, double* minStd, bool randInit) {
  unsigned long int T = pData->Row;
  unsigned long int tryCount;
  unsigned long int *idx = NULL;
  unsigned short int i, d, D = pData->Col;
  unsigned short int minMixtures = 1;
  double maximum, minimum, *mean, *std, *xmin, *xmax, *tmpVec, *startStd, one = 1.0, tmp;
  bool illConditioned = false, idxHasDuplicates = true, idxIsUnique; 
 
  MFree_1D(this->minStd);
  MArray_1D(this->minStd, D, double, "SC_BaggenstossEM.init: this->minStd");
  this->minStd = this->pMatrixFunc->copy(minStd, D);  

  //test if data is ill-conditioned
  if (T > 1) {
    mean = this->pMatrixFunc->mean(pData->Mat, T, D);
    std = this->pMatrixFunc->std(pData->Mat, T, D, mean);
  } else {
    mean = *(this->pMatrixFunc->scopy(pData->Mat, 1, D, (double)(1)));
    this->pMatrixFunc->fillWithValue(&std, 1, D, 0.0);
  }
  maximum = numeric_limits<double>::min();
  minimum = numeric_limits<double>::max();
  for (d = 0; d < D; d++) {
    if (fabs(mean[d]) > 1000*std[d]) {
      illConditioned = true;
      break;
    }
    if (std[d] > maximum) {maximum = std[d];}
    if (std[d] < minimum) {minimum = std[d];}
  }
  MFree_1D(mean);
  MFree_1D(std);
  assert((illConditioned == false) && (maximum <= 100*minimum));
  //if ((illConditioned == true) || (max > 100*min)) {REPORT_ERROR(SVLIB_BadData, "SC_BaggenstossEM.init: Possible ill-conditioning. Your features are shit. You should scale them or remove means!");}

  //compute initial value of covariances
  if (T > 1) {
    xmin = this->pMatrixFunc->min(pData->Mat, T, D);
    xmax = this->pMatrixFunc->max(pData->Mat, T, D);
  } else {
    xmin = *(this->pMatrixFunc->scopy(pData->Mat, 1, D, (double)(1)));
    xmax = *(this->pMatrixFunc->scopy(pData->Mat, 1, D, (double)(1)));
  }
  tmpVec = *(this->pMatrixFunc->sub(&xmax, &xmin, 1, D));
  MFree_1D(xmin);
  MFree_1D(xmax);
  this->pMatrixFunc->divide(&tmpVec, 1, D, 4.0);
  startStd = *(this->pMatrixFunc->max(&tmpVec, &minStd, 1, D));
  MFree_1D(tmpVec);

  //select starting means from a random set of data, make sure they are none the same
  if (randInit == true) {
    while ((idxHasDuplicates == true) && (mixtureCount >= minMixtures)) {
     tryCount = 0; //try 1000 times
      
      while ((idxHasDuplicates == true) && tryCount < 1000) { //generate mixtureCount random numbers in 1..T
        MFree_1D(tmpVec);
        MFree_1D(idx);
        
        tmpVec = this->pMatrixFunc->zeros(mixtureCount);
				tmp = (double)(T);
        this->pMatrixFunc->execComponentWise(&tmpVec, 1, mixtureCount, &sclib::rand, &one, &tmp);
				sclib::quickSort(tmpVec, 0, mixtureCount-1);
        idx = this->pMatrixFunc->scopy(tmpVec, mixtureCount, (unsigned long int)(1));
        
        if (mixtureCount > 1) { //see if they are each unique
          idxIsUnique = true;
          for (i = 0; i < mixtureCount-1; i++) {
            if (idx[i] == idx[i+1]) {
              idxIsUnique = false; 
              break;
            }
          }
          if (idxIsUnique == true) {idxHasDuplicates = false;}
        } else {
	        idxHasDuplicates = false;
        }

        tryCount++;
      } //while tryCount < 1000
      
      if (idxHasDuplicates == true) { //if this still fails after 1000 tries, reduce number of modes and keep trying
        mixtureCount--;
      }

    } //while mixtureCount >= minMixtures

    assert(idxHasDuplicates == false);
    //number of modes has fallen below minimum, fail
    //if (idxHasDuplicates == true) {REPORT_ERROR(SVLIB_BadData, "SC_BaggenstossEM.init: mixtureCount has been reduced to less than allowed, you need more data!");}
  } else { //no random init
    mixtureCount = sclib::min(mixtureCount, T);
    MArray_1D(idx, mixtureCount, unsigned long int, "SC_BaggenstossEM.init: idx");
    for (i = 0; i < mixtureCount; i++) {
      idx[i] = i*(T/mixtureCount) + 1; //changed from i+1; this way, maybe the resulting means are more likely to be different
    }
  } 
	MFree_1D(tmpVec);

  //set up class-parameters with the computed starting-covars, random means and equal weights
  MFree_1D(this->weight);
  MFree_2D(this->mean);
	if (this->choleskyCovar != NULL) {
		for (i = 0; i < this->mixtureCount; i++) {
			MFree_2D(this->choleskyCovar[i]);
			MFree_2D(this->U[i]);
			MFree_2D(this->V[i]);
			MFree_2D(this->S[i]);
			MFree_2D(this->transposedInvertedCholeskyCovar[i]);
		}
		MFree_1D(this->choleskyCovar);
		MFree_1D(this->U);
		MFree_1D(this->V);
		MFree_1D(this->S);
		MFree_1D(this->transposedInvertedCholeskyCovar);
	}
	MFree_1D(this->logDet);
  this->mixtureCount = mixtureCount;
  this->trainingDataCount = T;
  this->dim = D;
  MArray_1D(this->weight, this->mixtureCount, double, "SC_BaggenstossEM.init: weight");
  MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_BaggenstossEM.init: mean");
  MArray_1D(this->choleskyCovar, this->mixtureCount, double**, "SC_BaggenstossEM.init: choleskyCovar");
  MArray_1D(this->U, this->mixtureCount, double**, "SC_BaggenstossEM.init: U");
  MArray_1D(this->V, this->mixtureCount, double**, "SC_BaggenstossEM.init: V");
  MArray_1D(this->S, this->mixtureCount, double**, "SC_BaggenstossEM.init: S");
  MArray_1D(this->transposedInvertedCholeskyCovar, this->mixtureCount, double**, "SC_BaggenstossEM.init: transposedInvertedCholeskyCovar");
  MArray_1D(this->logDet, this->mixtureCount, double, "SC_BaggenstossEM.init: logDet");
  for (i = 0; i < this->mixtureCount; i++) {
    this->weight[i] = 1.0 / (double)(this->mixtureCount);
    for (d = 0; d < this->dim; d++) {
      this->mean[i][d] = pData->Mat[idx[i]-1][d];
    }
    this->choleskyCovar[i] = this->pMatrixFunc->diag(startStd, D);
		this->U[i] = NULL; //filled later on demand
		this->V[i] = NULL;
		this->S[i] = NULL;
		this->transposedInvertedCholeskyCovar[i] = NULL;
		this->logDet[i] = 0.0; //filled later on demand
  }

  MFree_1D(idx);
  MFree_1D(startStd);

  return;
}

//====================================================================================================================
//	C++ translation of P.M.Baggenstoss' function 
//  parm=gmix_trainscript(parm, data, nit, [S_P_M], [BIAS], [maxclose],[addmodes],[kurt_thresh]) for Matlab 
//  (04.06.1999). All Matlab-Like functions are located in SC_MatrixFunctions
//
//  Changes to original version:
//    - the last value of the log-likelihood-function is returned instead of nothing
//====================================================================================================================
double SC_BaggenstossEM::train(SV_Data* pData, unsigned int &actualIterations, unsigned long int maxIterations, long int samplesPerMode, bool bias, double maxCloseness, bool addModes, double splitThresh) {
  unsigned short int qIdx = 0, convergenceCounter = 0;
  unsigned short int trainPeriod = 6; //How often to check for close modes and try mode splitting
  unsigned short int iMerge, iSplit, nModesOld;
  unsigned long int i;
  double stopCriterion = 0.0001 * pData->Row; //How little an improvement is needed to continue
  double Q, best;
  
  //set meaningful parameters in case of default values
  if (samplesPerMode < 0) {samplesPerMode = 4 * pData->Col;}
  if (maxCloseness > 0.0) {maxCloseness = -2.0 * pData->Col;} //the docu suggests a default value of about -0.5*dim
  if (splitThresh < 0.0) {splitThresh = 1.0;}

  //If small number of modes, do splitting first, else do merging first
  if (this->mixtureCount < 3) {
      iMerge = 0;
      iSplit = trainPeriod / 2;
  } else {
      iMerge = trainPeriod / 2;
      iSplit = 0;
  }

  for (i = 0; i < maxIterations; i++) {
    Q = emStep(pData, bias);
    if (i == 0) {best = Q;}

    nModesOld = this->mixtureCount; //remeber old nr of mixtures
    
    //kill weak modes
    deflate(sclib::min((double)samplesPerMode/pData->Row, 0.5), 1.5/pData->Row);
   
    //split or merge adequate modes all trainPeriod iterations
    if (i < maxIterations-trainPeriod) {
      if (((long)(i) > 2*trainPeriod) && ((i % trainPeriod) == iMerge)) {
        merge(maxCloseness);
      }
      if ((addModes == true) && ((i % trainPeriod) == iSplit)) {
        split(pData, splitThresh);
      }
    } 

    //If there has been a change in modes, allow more time to converge
    if (nModesOld != this->mixtureCount) {
			convergenceCounter = 0;
		}

    //has there been an improvement?
    if (Q-best > stopCriterion) {
        //reset counter if there's an improvement
        convergenceCounter = 0;
    } else {
        //Increment "no improvement" counter. 
        //If no improvement has been seen in trainPeriod steps, terminate
        convergenceCounter++;
        if (convergenceCounter > trainPeriod) {
          break;
        }
    }
    if (Q > best) {best = Q;}


		if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbModelCreation) == true) {
			sclib::classOut("bagg.txt", this, this->pTweak);
		}
  }
	actualIterations = i;

  return Q;
}

//====================================================================================================================
//	C++ translation of P.M.Baggenstoss' function [parm,Q] = gmix_step(parm,x, [bias],[data_wts]) for Matlab 
//  (27.09.2001). All Matlab-Like functions are located in SC_MatrixFunctions
//
//  Changes to original version:
//    - no individual data-weights for each input feature vector
//====================================================================================================================
double SC_BaggenstossEM::emStep(SV_Data *pData, bool bias) {
  unsigned long int t, T = pData->Row;
  unsigned short int d, dd, D = this->dim, i, I = this->mixtureCount;
  float **z = pData->Mat;
  double **tmpIdx, **q = NULL, **r = NULL, **S = NULL, **US = NULL, **tmpM, *tmpV, *Svec, **minStdMatrix;

  assert(D == pData->Col);
  MArray_2D(tmpIdx, (long)(T), (short)(D), double, "SC_BaggenstossEM.emStep: tmpIdx");
	  
  //E-Step
  double mx, px, norm, **w, *alpha, Q, tmpQ, **sumWz, uu;
  MArray_2D(w, (short)(I), (long)(T), double, "SC_BaggenstossEM.emStep: w");
  alpha = this->pMatrixFunc->zeros(I);
  sumWz = this->pMatrixFunc->zeros(I, D);
  
  Q = 0.0;
  for (t = 0; t < T; t++) {
    mx = -1.0 * numeric_limits<double>::max();
    norm = 0.0;
    tmpQ = 0.0;

    for (i = 0; i < I; i++) { //compute the ingredients for the final weights w and priors alpha
      //px = lqrEval(z[t], this->mean[i], this->choleskyCovar[i], this->dim);
			px = lqrEval(z[t], this->mean[i], this->choleskyCovar[i], this->transposedInvertedCholeskyCovar[i], this->dim, this->logDet[i]);
      w[i][t] = log(this->weight[i]) + px;
      //tmpQ += exp(w[i][t]);
			tmpQ = (i == 0) ? w[i][t] : sclib::sLogAdd(tmpQ, w[i][t]); //work entirely in the log-domain
      mx = sclib::max(mx, w[i][t]);
    }

    for (i = 0; i < I; i++) {
      w[i][t] = exp(w[i][t] - mx);
      norm += w[i][t];
    }

		Q += tmpQ; //sclib::sLog(tmpQ); //compute the log-likelihood Q

    for (i = 0; i < I; i++) { //compute the final (properly normalized) w's and alphas in the log-domain
      w[i][t] /= norm;
      alpha[i] += w[i][t];
      
      for (d = 0; d < D; d++) {
        sumWz[i][d] += w[i][t] * z[t][d];
      }
    }
	} //for t


  //M-Step

  //take minStd^2 as a matrix (needed below)
  MArray_2D(minStdMatrix, D, D, double, "SC_BaggenstossEM.emStep: minStdMatrix");
  for (d = 0; d < D; d++) {
    for (dd = 0; dd < D; dd++) {
      minStdMatrix[d][dd] = (d == dd) ? this->minStd[d] * this->minStd[d] : 0.0;
    }
  }

  double sumAlpha = this->pMatrixFunc->sum(alpha, I);

  //Update means, weights, covars and constrain covariance
  for (i = 0; i < I; i++) {
    //update mean
    for (d = 0; d < D; d++) {
      this->mean[i][d] = sumWz[i][d] / alpha[i]; 
    }
    
    //prepare covar-trimming
    for (t = 0; t < T; t++) {
      uu = sqrt(w[i][t] / alpha[i]);
      for (d = 0; d < D; d++) {
        tmpIdx[t][d] = (z[t][d] - this->mean[i][d]) * uu;
      }
    }
    this->pMatrixFunc->qrDecomposition(tmpIdx, T, D, q, r);
    MFree_2D(q);
    //this->pMatrixFunc->triu(r, D); //already done inside qrDecomposition!
    
    //do covar-trimming
    if (bias == false) { //constraint method
      //first determine eigen decomp of C=R' * R = V * S^2 * V'
      this->pMatrixFunc->SVD(r, this->U[i], US, this->S[i], this->V[i], D, D); //fill the class' SVD cache
      assert(U != NULL);
      MFree_2D(US);
      MFree_2D(r);

      //S = diag(sclib::max(diag(S),sqrt(diag( V' * diag(min_std.^2) * V ))));
      tmpM = this->pMatrixFunc->mult(minStdMatrix, this->V[i], D, D, D, D);
      this->pMatrixFunc->transpose(this->V[i], D, D, true); //must be reverted later
      this->pMatrixFunc->mult(tmpM, this->V[i], D, D, D, D, true);
      tmpV = this->pMatrixFunc->diag(tmpM, D);
      MFree_2D(tmpM);
      this->pMatrixFunc->execComponentWise(&tmpV, 1, D, &sqrt); //should be only the diag-vector
      Svec = this->pMatrixFunc->diag(this->S[i], D);
      this->pMatrixFunc->max(&Svec, &tmpV, 1, D, true); //should be only the diag-vector of the result and both matrices
      MFree_1D(tmpV);
      S = this->pMatrixFunc->diag(Svec, D);
      MFree_1D(Svec);

      //tmpM = U * S * V'
      tmpM = this->pMatrixFunc->mult(this->U[i], S, D, D, D, D);
      this->pMatrixFunc->mult(tmpM, this->V[i], D, D, D, D, true);
      MFree_2D(S);
      
			MFree_2D(this->choleskyCovar[i]);
      MFree_2D(this->U[i]);
			MFree_2D(this->V[i]); //this saves us from reverting the transpose above...
			MFree_2D(this->S[i]);
			MFree_2D(this->transposedInvertedCholeskyCovar[i]);
      this->pMatrixFunc->qrDecomposition(tmpM, D, D, q, this->choleskyCovar[i]);
      MFree_2D(q);
      MFree_2D(tmpM);
    } else { //bias method
      tmpM = this->pMatrixFunc->concat(r, minStdMatrix, D, D, D, D);
      MFree_2D(r);
      MFree_2D(this->choleskyCovar[i]);
      this->pMatrixFunc->qrDecomposition(tmpM, 2*D, D, q, this->choleskyCovar[i]);
      MFree_2D(q);
      MFree_2D(tmpM);
			MFree_2D(this->U[i]);
			MFree_2D(this->V[i]);
			MFree_2D(this->S[i]);
			MFree_2D(this->transposedInvertedCholeskyCovar[i]);
    }

    //make the new cholesky-covar upper triangular
    //this->pMatrixFunc->triu(this->choleskyCovar[i], D); //not neccessary; done inside of qrDecomposition()
  
    //update weights
    this->weight[i] = sclib::max(alpha[i] / sumAlpha, 0.0);
  }

  MFree_2D(sumWz);
  MFree_2D(tmpIdx);
  MFree_2D(minStdMatrix);
  //MFree_2D(z);
  MFree_1D(alpha);
  MFree_2D(w);

  return Q;
}

//====================================================================================================================
//	Returns the log-likelihood of the Data given this model; Inspired by 
//  Inspired by P.M.Baggenstoss' Matlab-version function lg=lqr_evp(parm,x,flag)
//
//  Changes to original version:
//    - it's not a direct copy...
//====================================================================================================================
double SC_BaggenstossEM::test(SV_Data *pData) {
  unsigned short int i;
	double px, q, Q = 0.0;

	//determine likelihood
  for (unsigned long int t = 0; t < (unsigned long)pData->Col; t++) {
    q = 0.0;

    for (i = 0; i < this->mixtureCount; i++) {
      //px = lqrEval(pData->Mat[t], this->mean[i], this->choleskyCovar[i], this->dim);
			px = lqrEval(pData->Mat[t], this->mean[i], this->choleskyCovar[i], this->transposedInvertedCholeskyCovar[i], this->dim, this->logDet[i]);
      q += this->weight[i] * exp(px);
    }

    if (q <= 0.0) {
      Q += std::numeric_limits<double>::min();
    } else {
		  Q += log(q); //compute the likelihood of the data given the model
    }
  }

  return Q;
}

//====================================================================================================================
//	Evaluates the log-gaussian PDF for the feature-vector x with mean-vector mean and cholesky-covar choleskyCovar
//  Adopted from P.M.Baggenstoss' Matlab-version function l = lqr_eval(x,means,cholesky_covar)
//
//  Changes to original version:
//    - runs for just one feature vector instead of all at once
//====================================================================================================================
double SC_BaggenstossEM::lqrEval(double* x, double* mean, double** choleskyCovar, unsigned short dim) {
  double logDet = 0.0;
  double *tmpIdx;
  double **tmp1, *tmp2;
  unsigned short int d;
  double l = 0.0;
  int res;

  MArray_1D(tmpIdx, dim, double, "SC_BaggenstossEM.lqrEval: tmpIdx");
  
  for (d = 0; d < dim; d++) {
    logDet += log(fabs(choleskyCovar[d][d])); //the determinant of a triangular matrix is just the sum of the diagonal elements (the "trace")
    tmpIdx[d] = x[d] - mean[d];
  }
  logDet *= 2.0;

  tmp1 = this->pMatrixFunc->transpose(choleskyCovar, dim, dim);
  res = this->pMatrixFunc->inv(tmp1, dim);
  assert(res == 0); //inverting tmp1 should work!
  tmp2 = this->pMatrixFunc->mult(tmp1, tmpIdx, dim, dim, dim);
  MFree_2D(tmp1);
  MFree_1D(tmpIdx);

  for (d = 0; d < dim; d++) {
    l += tmp2[d] * tmp2[d];
  }
  MFree_1D(tmp2);

	l = -0.5*l - ((double)(dim) / 2.0)*sclib::log_2pi - 0.5*logDet;

  return l;
}

//====================================================================================================================
// This is a altered version of the above function: It takes the data in a float-vector instead of a double-vector
//====================================================================================================================
double SC_BaggenstossEM::lqrEval(float* x, double* mean, double** choleskyCovar, unsigned short dim) {
  double logDet = 0.0;
  double *tmpIdx;
  double **tmp1, *tmp2;
  unsigned short int d;
  double l = 0.0;
  int res;

  MArray_1D(tmpIdx, dim, double, "SC_BaggenstossEM.lqrEval: tmpIdx");
  
  for (d = 0; d < dim; d++) {
    logDet += log(fabs(choleskyCovar[d][d])); //the determinant of a triangular matrix is just the sum of the diagonal elements (the "trace")
    tmpIdx[d] = x[d] - mean[d];
  }
  logDet *= 2.0;

  tmp1 = this->pMatrixFunc->transpose(choleskyCovar, dim, dim);
  res = this->pMatrixFunc->inv(tmp1, dim);
  assert(res == 0); //inverting tmp1 should work!
  tmp2 = this->pMatrixFunc->mult(tmp1, tmpIdx, dim, dim, dim);
  MFree_2D(tmp1);
  MFree_1D(tmpIdx);

  for (d = 0; d < dim; d++) {
    l += tmp2[d] * tmp2[d];
  }
  MFree_1D(tmp2);

	l = -0.5*l - ((double)(dim) / 2.0)*sclib::log_2pi - 0.5*logDet;

  return l;
}

//====================================================================================================================
// This is a altered version of the above function: It takes the precomputed lod-determinant of the cholesky-covar and
// transposed-inverted form to speed up things a little bit
//====================================================================================================================
double SC_BaggenstossEM::lqrEval(float *x, double *mean, double **choleskyCovar, double** &transposedInvertedCholeskyCovar, unsigned short dim, double &logDet) {
  double *tmpIdx;
  double *tmp;
  unsigned short int d;
  double l = 0.0;
	int res;

  MArray_1D(tmpIdx, dim, double, "SC_BaggenstossEM.lqrEval: tmpIdx");
	if (transposedInvertedCholeskyCovar == NULL) { //compute and return on use
		logDet = 0.0;
		for (d = 0; d < dim; d++) {
			logDet += log(fabs(choleskyCovar[d][d])); //the determinant of a triangular matrix is just the sum of the diagonal elements (the "trace")
			tmpIdx[d] = x[d] - mean[d];
		}
		logDet *= 2.0;

		transposedInvertedCholeskyCovar = this->pMatrixFunc->transpose(choleskyCovar, dim, dim);
		res = this->pMatrixFunc->inv(transposedInvertedCholeskyCovar, dim);
		assert(res == 0); //inverting the transposed cholesky covar should work!
	} else {
		for (d = 0; d < dim; d++) {
			tmpIdx[d] = x[d] - mean[d];
		}
	}
  tmp = this->pMatrixFunc->mult(transposedInvertedCholeskyCovar, tmpIdx, dim, dim, dim);
  MFree_1D(tmpIdx);

  for (d = 0; d < dim; d++) {
    l += tmp[d] * tmp[d];
  }
  MFree_1D(tmp);

	l = -0.5*l - ((double)(dim) / 2.0)*sclib::log_2pi - 0.5*logDet;

  return l;
}

//====================================================================================================================
// This is a altered version of the above function to accept a double input parameter as first argument
//====================================================================================================================
double SC_BaggenstossEM::lqrEval(double *x, double *mean, double **choleskyCovar, double** &transposedInvertedCholeskyCovar, unsigned short dim, double &logDet) {
  double *tmpIdx;
  double *tmp;
  unsigned short int d;
  double l = 0.0;
	int res;

  MArray_1D(tmpIdx, dim, double, "SC_BaggenstossEM.lqrEval: tmpIdx");
	if (transposedInvertedCholeskyCovar == NULL) { //compute and return on use
		logDet = 0.0;
		for (d = 0; d < dim; d++) {
			logDet += log(fabs(choleskyCovar[d][d])); //the determinant of a triangular matrix is just the sum of the diagonal elements (the "trace")
			tmpIdx[d] = x[d] - mean[d];
		}
		logDet *= 2.0;

		transposedInvertedCholeskyCovar = this->pMatrixFunc->transpose(choleskyCovar, dim, dim);
		res = this->pMatrixFunc->inv(transposedInvertedCholeskyCovar, dim);
		assert(res == 0); //inverting the transposed cholesky covar should work!
	} else {
		for (d = 0; d < dim; d++) {
			tmpIdx[d] = x[d] - mean[d];
		}
	}
  tmp = this->pMatrixFunc->mult(transposedInvertedCholeskyCovar, tmpIdx, dim, dim, dim);
  MFree_1D(tmpIdx);

  for (d = 0; d < dim; d++) {
    l += tmp[d] * tmp[d];
  }
  MFree_1D(tmp);

	l = -0.5*l - ((double)(dim) / 2.0)*sclib::log_2pi - 0.5*logDet;

  return l;
}

//====================================================================================================================
//	Kills weak mixtures according to their weight. All mixtures with weight<minWeightAll are killed, but only one 
//  with weight<minWeightOne.
//  Adopted from P.M.Baggenstoss' Matlab-version function parm = gmix_deflate(parm,min_weight_1,min_weight_all)
//
//  Changes to original version:
//    - returns the count of killed mixtures
//====================================================================================================================
unsigned short int SC_BaggenstossEM::deflate(double minWeightOne, double minWeightAll) {
  unsigned long int d, i, j, minWeightIdx;
  unsigned short int newMixtureCount, deflateCount = 0;
  double minWeight, sumWeights;
  double *newWeight, **newMean, ***newCholeskyCovar, ***newU, ***newV, ***newS, ***newTransposed, *newLogDet;
  bool proceed = true;

  if (this->mixtureCount > 1) {
    while (proceed == true) {
      minWeight = this->pMatrixFunc->min(this->weight, this->mixtureCount, &minWeightIdx);
      if ((minWeight < minWeightAll) || ((minWeight < minWeightOne) && (deflateCount == 0))) {
        if (minWeight >= minWeightAll) {deflateCount++;}
          
          //create new model parameters (by just leaving the deflated one out and renormalizing the weights)
          sumWeights = this->pMatrixFunc->sum(this->weight, this->mixtureCount) - this->weight[minWeightIdx];
          newMixtureCount = this->mixtureCount - 1;
          MArray_1D(newWeight, (long)newMixtureCount, double, "SC_BaggesnstossEM.deflate: newWeight");
          MArray_2D(newMean, (long)newMixtureCount, this->dim, double, "SC_BaggenstossEM.deflate: newMean");
          MArray_1D(newCholeskyCovar, (long)newMixtureCount, double**, "SC_BaggenstossEM.deflate: newCholeskyCovar");
					MArray_1D(newU, (long)newMixtureCount, double**, "SC_BaggenstossEM.deflate: newU");
					MArray_1D(newV, (long)newMixtureCount, double**, "SC_BaggenstossEM.deflate: newV");
					MArray_1D(newS, (long)newMixtureCount, double**, "SC_BaggenstossEM.deflate: newS");
					MArray_1D(newTransposed, (long)newMixtureCount, double**, "SC_BaggenstossEM.deflate: newTransposed");
					MArray_1D(newLogDet, (long)newMixtureCount, double, "SC_BaggenstossEM.deflate: newLogDet");
          for (i = 0, j = 0; i < this->mixtureCount; i++) {
            if (i != minWeightIdx) {
              newWeight[j] = this->weight[i] / sumWeights;
              newCholeskyCovar[j] = this->choleskyCovar[i];
							newU[j] = this->U[i];
							newV[j] = this->V[i];
							newS[j] = this->S[i];
							newTransposed[j] = this->transposedInvertedCholeskyCovar[i];
							newLogDet[j] = this->logDet[i];
              for (d = 0; d < this->dim; d++) {
                newMean[j][d] = this->mean[i][d];
              }
              j++;
						} else {
							MFree_2D(this->choleskyCovar[i]);
							MFree_2D(this->U[i]);
							MFree_2D(this->V[i]);
							MFree_2D(this->S[i]);
							MFree_2D(this->transposedInvertedCholeskyCovar[i]);
						}
          }

          //make the new parameters the actual ones
          MFree_1D(this->weight);
          MFree_2D(this->mean);
					MFree_1D(this->choleskyCovar);
					MFree_1D(this->U);
					MFree_1D(this->V);
					MFree_1D(this->S);
					MFree_1D(this->transposedInvertedCholeskyCovar);
					MFree_1D(this->logDet);
          this->weight = newWeight;
          this->mean = newMean;
          this->choleskyCovar = newCholeskyCovar;
					this->U = newU;
					this->V = newV;
					this->S = newS;
					this->transposedInvertedCholeskyCovar = newTransposed;
					this->logDet = newLogDet;
          this->mixtureCount = newMixtureCount;

      } else {
        proceed = false;
      }
    } //while proceed
  } //if mixtureCount>1

  return deflateCount;
}

//====================================================================================================================
//	Determines the closeness of two mixtures; zero means that the modes are identical.
//  Adopted from P.M.Baggenstoss' Matlab-version function d = mode_dist(m1,m2,v1,v2)
//
//  Changes to original version:
//    -
//====================================================================================================================
double SC_BaggenstossEM::mixtureCloseness(double* mean1, double* mean2, double** choleskyCovar1, double** choleskyCovar2, unsigned short int dim, double** &S1, double** &S2, double** &V1, double** &V2, double** &U1, double** &U2, double** transposedInvertedCholesky1, double** transposedInvertedCholesky2, double logDet1, double logDet2) {
  unsigned long int d, ev, i, t = 0, T = 2 * (2*dim + 1);
  double **x, **choleskyCovar, *mean, **U = NULL, **S = NULL, *s, **US = NULL, **V = NULL; 
  double p1x1 = 0.0, p1x2 = 0.0, p2x1 = 0.0, p2x2 = 0.0;
	bool svd = false;

  MArray_2D(x, (long)T, (long)dim, double, "SC_baggenstossEM.mixtureCloseness: x");

  //loop over the two mixtures and synthesize T samples
  for (i = 0; i < 2; i++) {

    //select the right values for this round
    if (i==0) {
      choleskyCovar = choleskyCovar1;
      mean = mean1; 
			V = V1;
			S = S1;
			U = U1;
    } else { 
      choleskyCovar = choleskyCovar2; 
      mean = mean2; 
			V = V2;
			S = S2;
			U = U2;
    }

    //get the eigenvectors (right singular vectors (columns of V) eig(tr)'
		if (S == NULL || V == NULL || U == NULL) {
			this->pMatrixFunc->SVD(choleskyCovar, U, US, S, V, dim, dim, false);
			MFree_2D(US);
		}
    s = this->pMatrixFunc->diag(S, dim);

    //add the center point 
    this->pMatrixFunc->setRow(x, mean, T, dim, t++);

    //loop over eigenvectors 
    for (ev = 0; ev < dim; ev++) {
      for (d = 0; d < dim; d++) {
        x[t][d] = mean[d] + (s[ev] * V[d][ev]);
        x[t+1][d] = mean[d] - (s[ev] * V[d][ev]);
      }
      t = t + 2;
    }
    MFree_1D(s);
  }
  
  //compute the cross-likelihoods
  for (t = 0; t < T; t++) {
    if (t > T/2) {
			p1x2 += (transposedInvertedCholesky1 == NULL) ? lqrEval(x[t], mean1, choleskyCovar1, dim) : lqrEval(x[t], mean1, choleskyCovar1, transposedInvertedCholesky1, dim, logDet1);
			p2x2 += (transposedInvertedCholesky2 == NULL) ? lqrEval(x[t], mean2, choleskyCovar2, dim) : lqrEval(x[t], mean2, choleskyCovar2, transposedInvertedCholesky2, dim, logDet2);
    } else {
			p2x1 += (transposedInvertedCholesky2 == NULL) ? lqrEval(x[t], mean2, choleskyCovar2, dim) : lqrEval(x[t], mean2, choleskyCovar2, transposedInvertedCholesky2, dim, logDet2);
			p1x1 += (transposedInvertedCholesky1 == NULL) ? lqrEval(x[t], mean1, choleskyCovar1, dim) : lqrEval(x[t], mean1, choleskyCovar1, transposedInvertedCholesky1, dim, logDet1);
    }
  }

  MFree_2D(x);

  return (p1x2 + p2x1) - (p2x2 + p1x1);
}

//====================================================================================================================
//	Merges mixtures closer than the threshold maxCloseness, which should be -1*dim. Identical mixtures have 
//  closeness zero.
//  Adopted from P.M.Baggenstoss' Matlab-version function parm = gmix_merge(parm,max_closeness)
//
//  Changes to original version:
//    - returns the count of merged mixtures
//====================================================================================================================
unsigned short int SC_BaggenstossEM::merge(double maxCloseness) {
  unsigned short int d, dd, i, j, mergeCounter = 0;
  double closeness;
  double newWeight, *newMean, **newCholeskyCovar;
  double w1, w2, **v1, **v2, **q = NULL, **tmp, sqrtDim = sqrt((double)this->dim);

  if (this->mixtureCount > 1) {

    MArray_1D(newMean, this->dim, double, "SC_BaggenstossEM.merge: newMean");
    MArray_2D(newCholeskyCovar, this->dim, this->dim, double, "SC_BaggenstossEM.merge: newCholeskyCovar");

    for (i = 0; i < this->mixtureCount; i++) {
      for (j = i+1; j < this->mixtureCount; j++) {
        if ((this->weight[i] > 0) && (this->weight[j] > 0)) {
          closeness = mixtureCloseness(this->mean[i], this->mean[j], this->choleskyCovar[i], this->choleskyCovar[j], this->dim, this->S[i], this->S[j], this->V[i], this->V[j], this->U[i], this->U[j], this->transposedInvertedCholeskyCovar[i], this->transposedInvertedCholeskyCovar[j], this->logDet[i], this->logDet[j]);
          if ((closeness > maxCloseness) && (mergeCounter < (this->mixtureCount/2 +1))) {
            
            mergeCounter++;
            
            //compute new weight
            newWeight = this->weight[i] + this->weight[j];
            w1 = this->weight[i] / newWeight;
            w2 = this->weight[j] / newWeight;
            
            //the central mean by weighted average
            for (d = 0; d < this->dim; d++) {
              newMean[d] = w1*this->mean[i][d] + w2*this->mean[j][d];
            }
            
            //the ROWS of v times the square root of DIM (v is the the QR of covariance) can be considered 
            //data samples about the means of each distribution.  
            v1 = this->pMatrixFunc->mult(this->choleskyCovar[i], sqrtDim, this->dim, this->dim);
            v2 = this->pMatrixFunc->mult(this->choleskyCovar[j], sqrtDim, this->dim, this->dim);

            //they need to be re-referenced to the new center
            for (d = 0; d < this->dim; d++) {
              for (dd = 0; dd < this->dim; dd++) {
	              v1[d][dd] = v1[d][dd] + (this->mean[i][dd] - newMean[dd]);
                v2[d][dd] = v2[d][dd] + (this->mean[j][dd] - newMean[dd]);
              }
            }

            //form a weighted augmented matrix 
            this->pMatrixFunc->mult(v1, sqrt(w1), this->dim, this->dim, true);
            this->pMatrixFunc->mult(v2, sqrt(w2), this->dim, this->dim, true);
            tmp = this->pMatrixFunc->concat(v1, v2, this->dim, this->dim, this->dim, this->dim);
            MFree_2D(v1);
            MFree_2D(v2);
            this->pMatrixFunc->qrDecomposition(tmp, 2*this->dim, this->dim, q, newCholeskyCovar);
            MFree_2D(tmp);
            MFree_2D(q);
            this->pMatrixFunc->mult(newCholeskyCovar, 1.0/sqrtDim, this->dim, this->dim, true);

            //update parameters
            if (this->weight[i] > this->weight[j]) { //merging mixture j into mixture i
        			this->weight[i] = newWeight;
              this->weight[j] = 0.0;
              MFree_2D(this->choleskyCovar[i]);
              this->choleskyCovar[i] = newCholeskyCovar;
              this->pMatrixFunc->setRow(this->mean, newMean, this->mixtureCount, this->dim, i);
							MFree_2D(this->U[i]);
							MFree_2D(this->V[i]);
							MFree_2D(this->S[i]);
							MFree_2D(this->transposedInvertedCholeskyCovar[i]);
            } else { //merging mixture i into mixture j
        			this->weight[j] = newWeight;
              this->weight[i] = 0.0;
              MFree_2D(this->choleskyCovar[j]);
              this->choleskyCovar[j] = newCholeskyCovar;
							MFree_2D(this->U[j]);
							MFree_2D(this->V[j]);
							MFree_2D(this->S[j]);
							MFree_2D(this->transposedInvertedCholeskyCovar[j]);
              this->pMatrixFunc->setRow(this->mean, newMean, this->mixtureCount, this->dim, j);
            }
            newCholeskyCovar = NULL;

          } //really merge this mixtures
        } //weights > 0
        if (this->weight[i] == 0) {break;}
      } //for j
    } //for i

    MFree_1D(newMean);
    MFree_2D(newCholeskyCovar);

	} //mixtureCount>1

  deflate(1.0e-10, 1.0e-10); //finally delete those modes which have been merged into another one

  return mergeCounter;
}

//====================================================================================================================
//	Kurtosis-based splitting of mixtures based on N. Vlassis A. Likas,  "The Kurtosis-EM algorithm for Gaussian 
//  mixture modelling", to be published IEEE SMC Transactions, in 1999
//  Adopted from P.M.Baggenstoss' Matlab-version function parm = gmix_kurt(parm,x,thresh,debug);
//
//  Changes to original version:
//    - returns the count of splitted mixtures
//====================================================================================================================
unsigned short int SC_BaggenstossEM::split(SV_Data *pData, double threshold) {
  unsigned short int dd, d, i, I = this->mixtureCount, D = this->dim, splitCounter = 0;
  unsigned long int t, T = pData->Row;
  float **z = pData->Mat;
  double mx, px, norm, **w;
  double **US = NULL, *sumW;

  MArray_2D(w, (short)(I), (long)(T), double, "SC_BaggenstossEM.emStep: w");
  sumW = this->pMatrixFunc->zeros(I);
  	
  //compute membership probabilities of modes for each sample, store it in w[I][T]
  for (t = 0; t < T; t++) {
    mx = -1.0 * numeric_limits<double>::max();
    norm = 0.0;

    for (i = 0; i < I; i++) { //compute the ingredients for the final weights w
      //px = lqrEval(z[t], this->mean[i], this->choleskyCovar[i], this->dim);
			px = lqrEval(z[t], this->mean[i], this->choleskyCovar[i], this->transposedInvertedCholeskyCovar[i], this->dim, this->logDet[i]);
      w[i][t] = log(this->weight[i]) + px;
      mx = sclib::max(mx, w[i][t]);
    }

    for (i = 0; i < I; i++) {
      w[i][t] = exp(w[i][t] - mx);
      norm += w[i][t];
    }

    for (i = 0; i < I; i++) { //compute the final (properly normalized) w's in the log-domain
      w[i][t] /= norm;
      sumW[i] += w[i][t];
    }
  }
  
  double k1, k2, k3, k4, tmp, **zt, fac, merit, mMax, *mAll;
  unsigned short int dMax, *dAll;
  MArray_2D(zt, (long)T, (short)D, double, "SC_BaggenstossEM.split: zt");
  MArray_1D(dAll, I, unsigned short int, "SC_BaggenstossEM.split: dAll");
  MArray_1D(mAll, I, double, "SC_BaggenstossEM.split: mAll");

  //check which mixtures need splitting
  for (i = 0; i < I; i++) {
    //get axes
		if (this->S[i] == NULL || this->V[i] == NULL) {
			this->pMatrixFunc->SVD(this->choleskyCovar[i], this->U[i], US, this->S[i], this->V[i], D, D, false);
			MFree_2D(US);
		}

    //remove mean
    for (t = 0; t < T; t++) {
      for (d = 0; d < this->dim; d++) {
        zt[t][d] = z[t][d] - this->mean[i][d];
      }
    }
    
    //for each dimension, project onto this vector
    mMax = -1.0 * numeric_limits<double>::max();
		dMax = 0;
    for (d = 0; d < D; d++) {
      k1 = 0.0; k2 = 0.0; k3 = 0.0; k4 = 0.0;

      for (t = 0; t < T; t++) {
        tmp = 0;
        for (dd = 0; dd < D; dd++) {
          tmp += zt[t][dd] * this->V[i][dd][d];
        }
        tmp /= this->S[i][d][d];      //tmp = zt * V(:,j) / S(j,j)

        k1 += w[i][t] * tmp; //pow(tmp, 1.0);
        k2 += w[i][t] * tmp*tmp; //pow(tmp, 2.0);
        k3 += w[i][t] * tmp*tmp*tmp; //pow(tmp, 3.0);
        k4 += w[i][t] * tmp*tmp*tmp*tmp; //pow(tmp, 4.0);
      }

			k1 /= (sumW[i] > 0) ? sumW[i] : 1.0; 
      k2 /= (sumW[i] > 0) ? sumW[i] : 1.0; 
      k3 /= (sumW[i] > 0) ? sumW[i] : 1.0; 
      k4 /= (sumW[i] > 0) ? sumW[i] : 1.0; 

      if (((sumW[i] / 30.0) / (double)D) > 10.0) {
        fac = 1.0;
      } else if (((sumW[i] / 30.9) / (double)D) > 2.0) {
        fac = 0.75;
      } else if (((sumW[i] / 30.0) / (double)D) > 1.0) {
        fac = 0.5;
      } else {
        fac = 0.0;
      }
      
      merit = 2.0*(fabs(k4-3.0)  + fabs(k3)) * fac * k2;
      if(merit > mMax) {
        mMax = merit; 
        dMax = d;
      }
    } //for d

    mAll[i] = mMax;
    dAll[i] = dMax;
  } //for i

  MFree_2D(zt);
	MFree_2D(w);
	MFree_1D(sumW);

  double *m1, *m2, *newWeight, **newMean, ***newCholeskyCovar, ***newU, ***newV, ***newS, ***newTransposed, *newLogDet;
  unsigned short int ii, newMixtureCount;

  //do the splitting
  for (i = 0; i < I; i++) { //yes, 'I' is the old mixtureCount!!!
    if (mAll[i] > threshold) {
      if (this->mixtureCount < this->maxMixtureCount) {
        dMax = dAll[i];
        
        //get axes
				if (this->S[i] == NULL || this->V[i] == NULL) {
					this->pMatrixFunc->SVD(this->choleskyCovar[i], this->U[i], US, this->S[i], this->V[i], D, D, false);
	        MFree_2D(US);
				}

        //compute new means;
        MArray_1D(m1, D, double, "SC_BaggenstossEM.split: m1");
        MArray_1D(m2, D, double, "SC_BaggenstossEM.split: m2");
        for (d = 0; d < D; d++) {
          m1[d] = this->mean[i][d] + (this->V[i][d][dMax] * (this->S[i][dMax][dMax] / 2.0));
          m2[d] = this->mean[i][d] - (this->V[i][d][dMax] * (this->S[i][dMax][dMax] / 2.0));
        }

        //compute parameters
        newMixtureCount = this->mixtureCount + 1;
        MArray_1D(newWeight, (long)newMixtureCount, double, "SC_BaggesnstossEM.split: newWeight");
        MArray_2D(newMean, (long)newMixtureCount, this->dim, double, "SC_BaggenstossEM.split: newMean");
        MArray_1D(newCholeskyCovar, (long)newMixtureCount, double**, "SC_BaggenstossEM.split: newCholeskyCovar");
				MArray_1D(newU, (long)newMixtureCount, double**, "SC_BaggenstossEM.split: newU");
				MArray_1D(newV, (long)newMixtureCount, double**, "SC_BaggenstossEM.split: newV");
				MArray_1D(newS, (long)newMixtureCount, double**, "SC_BaggenstossEM.split: newS");
				MArray_1D(newTransposed, (long)newMixtureCount, double**, "SC_BaggenstossEM.split: newTransposed");
				MArray_1D(newLogDet, (long)newMixtureCount, double, "SC_BaggenstossEM.split: newLogDet");
        for (ii = 0; ii < newMixtureCount; ii++) {
          if ((ii < this->mixtureCount) && (ii != i)) { //just copy old parameters
            newWeight[ii] = this->weight[ii];
            newCholeskyCovar[ii] = this->choleskyCovar[ii];
						newU[ii] = this->U[ii];
						newV[ii] = this->V[ii];
						newS[ii] = this->S[ii];
						newTransposed[ii] = this->transposedInvertedCholeskyCovar[ii];
						newLogDet[ii] = this->logDet[ii];
            for (d = 0; d < D; d++) {
              newMean[ii][d] = this->mean[ii][d];
            }
          } else if ((ii < this->mixtureCount) && (ii == i)) { //alter the splitted mixture
            newWeight[ii] = this->weight[i] / 2.0;
            newCholeskyCovar[ii] = this->choleskyCovar[ii];
						newU[ii] = this->U[ii];
						newV[ii] = this->V[ii];
						newS[ii] = this->S[ii];
						newTransposed[ii] = this->transposedInvertedCholeskyCovar[ii];
						newLogDet[ii] = this->logDet[ii];
            for (d = 0; d < D; d++) {
              newMean[ii][d] = m1[d];
            }
          } else if (ii >= this->mixtureCount) { //generate the new mixture
            newWeight[ii] = this->weight[i] / 2.0;
            newCholeskyCovar[ii] = this->pMatrixFunc->copy(this->choleskyCovar[i], D, D);
						newU[ii] = this->pMatrixFunc->copy(this->U[i], D, D);
						newV[ii] = this->pMatrixFunc->copy(this->V[i], D, D);
						newS[ii] = this->pMatrixFunc->copy(this->S[i], D, D);
						newTransposed[ii] = this->pMatrixFunc->copy(this->transposedInvertedCholeskyCovar[i], D, D);
						newLogDet[ii] = this->logDet[ii];
            for (d = 0; d < D; d++) {
              newMean[ii][d] = m2[d];
            }
          }
        }
        MFree_1D(m1);
        MFree_1D(m2);

        //set new parameters
        MFree_1D(this->weight);
        MFree_2D(this->mean);
				MFree_1D(this->choleskyCovar);
				MFree_1D(this->U);
				MFree_1D(this->V);
				MFree_1D(this->S);
				MFree_1D(this->transposedInvertedCholeskyCovar);
				MFree_1D(this->logDet);
        this->weight = newWeight;
        this->mean = newMean;
        this->choleskyCovar = newCholeskyCovar;
				this->U = newU;
				this->V = newV;
				this->S = newS;
				this->transposedInvertedCholeskyCovar = newTransposed;
				this->logDet = newLogDet;
        this->mixtureCount = newMixtureCount;
        newWeight = NULL;
        newMean = NULL;
        newCholeskyCovar = NULL;
				newU = NULL;
				newV = NULL;
				newS = NULL;
				newTransposed = NULL;
				newLogDet = NULL;

      } else { //maxMixtureCount reached
        return splitCounter;
      }
    } //if mAll[i] > threshold
  } //for i

  MFree_1D(mAll);
  MFree_1D(dAll);
 
  return splitCounter;
}

//====================================================================================================================
//  Computes and returns the covariance-matrix out of its stored cholesky-decomposition
//====================================================================================================================
double** SC_BaggenstossEM::getCovar(unsigned short int mixture) {
  double **coVar;

  coVar = this->pMatrixFunc->transpose(this->choleskyCovar[mixture], this->dim, this->dim);
  this->pMatrixFunc->mult(coVar, this->choleskyCovar[mixture], this->dim, this->dim, this->dim, this->dim, true); //result resides in covar
	//TODO: mutliply covar with 1/N?

  return coVar;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& operator<< (ostream& OutS, SC_BaggenstossEM& Data) {
	unsigned short int x, y, z;

	OutS.setf(ios::fixed|ios::basefield);
  OutS.precision(5);

  if (Data.verbose == true) {
	  OutS << "MixNr\tWeight\tMean" << endl;
  } else {
    OutS << "MixNr\tWeight\t|CholCvar|\tMean" << endl;
  }

	for (x = 0; x < Data.mixtureCount; x++) {
		if (Data.verbose == true) {
      OutS << x << "\t" << setw(7) << Data.weight[x] << "\t";
    } else {
      OutS << x << "\t" << setw(7) << Data.weight[x] << "\t" << setw(10) << Data.pMatrixFunc->det(Data.choleskyCovar[x], Data.dim) << "\t";
    }    
		for (y = 0; y < Data.dim; y++) {
			OutS << setw(10) << Data.mean[x][y] << ((y == Data.dim-1) ? "" : " ");
      //OutS << Data.mean[x][y] << ((y == Data.dim-1) ? "" : " ");
    }
		OutS << endl;
	}

  if (Data.verbose == true) {
	  for (x = 0; x < Data.mixtureCount; x++) {
  	  OutS << endl << "MixNr " << x << " Cholesky-Covariance" << endl;
		  for (y = 0; y < Data.dim; y++) {
        for (z = 0; z < Data.dim; z++) {
				  OutS << setw(10) << Data.choleskyCovar[x][y][z] << ((y == Data.dim-1) ? "" : " ");
          //OutS << Data.choleskyCovar[x][y][z] << ((y == Data.dim-1) ? "" : " ");
        }
        OutS << endl;
		  }	
		  OutS << endl;
	  }
  }

	OutS << endl;

	return(OutS);
}

//====================================================================================================================
// write parameters to an already opened, ready-to-write binary file
//====================================================================================================================
unsigned int SC_BaggenstossEM::write(fstream *file) {
	unsigned int bytes = 0;
	SV_DataIO io;

	bytes = io.writeScalar(file, this->verbose);
	bytes += io.writeScalar(file, this->maxMixtureCount);
	bytes += io.writeScalar(file, this->mixtureCount);
	bytes += io.writeScalar(file, this->dim);
	bytes += io.writeScalar(file, this->trainingDataCount);
	bytes += io.writeArray(file, this->minStd, this->dim);
	bytes += io.writeArray(file, this->weight, this->mixtureCount);
	bytes += io.writeMatrix(file, this->mean, this->mixtureCount, this->dim);
	for (int i = 0; i < this->mixtureCount; i++) {
		bytes += io.writeMatrix(file, this->choleskyCovar[i], this->dim, this->dim);
	}

	return bytes;
}

//====================================================================================================================
// read parameters from an already opened, correctly positioned binary file
//====================================================================================================================
unsigned int SC_BaggenstossEM::read(fstream *file, SV_DataIO::SV_DatatypeSizes *fileSizes) {
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes;
	io.getCurrentDatatypeSizes(codeSizes);

	io.readScalar(file, this->verbose, codeSizes, *fileSizes);
	io.readScalar(file, this->maxMixtureCount, codeSizes, *fileSizes);
	io.readScalar(file, this->mixtureCount, codeSizes, *fileSizes);
	io.readScalar(file, this->dim, codeSizes, *fileSizes);
	io.readScalar(file, this->trainingDataCount, codeSizes, *fileSizes);
	io.readArray(file, this->minStd, this->dim, codeSizes, *fileSizes);
	io.readArray(file, this->weight, this->mixtureCount, codeSizes, *fileSizes);
	io.readMatrix(file, this->mean, this->mixtureCount, this->dim, codeSizes, *fileSizes);
	for (int i = 0; i < this->mixtureCount; i++) {
			io.readMatrix(file, this->choleskyCovar[i], this->dim, this->dim, codeSizes, *fileSizes);
	}

	return (file->good() != TRUE) ? SVLIB_Fail : SVLIB_Ok;
}

//====================================================================================================================
// Below is some test-code for the SVD/qr/cholesky-functions in SC_MatrixFunctions; it can be inserted into the 
// emStep-function and should work there
//====================================================================================================================
/*
  //--------------------------- some tests

  //create a 4*4 matrix
  double **mat1;
  MArray_2D(mat1, 4, 4, double, "");
  mat1[0][0] = 2; mat1[0][1] = 3; mat1[0][2] = 7; mat1[0][3] = 1;
  mat1[1][0] = 4; mat1[1][1] = 0; mat1[1][2] = 8; mat1[1][3] = 6;
  mat1[2][0] = 9; mat1[2][1] = 8; mat1[2][2] = 9; mat1[2][3] = 4;
  mat1[3][0] = 0; mat1[3][1] = 4; mat1[3][2] = 3; mat1[3][3] = 1;

  //test svd
  pMatrixFunc->SVD(mat1, U, US, S, V, 4, 4);

  //calculate back from the SVD to mat1
  double **tmp1, **tmp2;
  pMatrixFunc->transpose(V, 4, 4, true);
  tmp1 = pMatrixFunc->mult(S, V, 4, 4, 4, 4);
  tmp2 = pMatrixFunc->mult(U, tmp1, 4, 4, 4, 4);

  //the error should be of the order of the machine-precision: 10^(-15) to 10^(-16)
  double **diff;
  diff = pMatrixFunc->sub(mat1, tmp2, 4, 4);

  MFree_2D(diff);
  MFree_2D(tmp1);
  MFree_2D(tmp2);
  MFree_2D(U);
  MFree_2D(US);
  MFree_2D(S);
  MFree_2D(V);

  //test qr
  pMatrixFunc->qrDecomposition(mat1, 4, 4, q, r);
  tmp2 = pMatrixFunc->mult(q, r, 4, 4, 4, 4);
  diff = pMatrixFunc->sub(mat1, tmp2, 4, 4);

  MFree_2D(tmp2);
  MFree_2D(diff);
  MFree_2D(q);
  MFree_2D(r);

  //test cholesky
  MFree_2D(mat1);
  MArray_2D(mat1, 9, 9, double, "");
  mat1[0][0] = 1; mat1[0][1] = 1; mat1[0][2] = 1; mat1[0][3] = 1; mat1[0][4] = 1; mat1[0][5] = 1; mat1[0][6] = 1; mat1[0][7] = 1; mat1[0][8] = 0; 
  mat1[1][0] = 1; mat1[1][1] = 2; mat1[1][2] = 2; mat1[1][3] = 2; mat1[1][4] = 2; mat1[1][5] = 2; mat1[1][6] = 2; mat1[1][7] = 2; mat1[1][8] = 1; 
  mat1[2][0] = 1; mat1[2][1] = 2; mat1[2][2] = 3; mat1[2][3] = 3; mat1[2][4] = 3; mat1[2][5] = 3; mat1[2][6] = 3; mat1[2][7] = 3; mat1[2][8] = 2; 
  mat1[3][0] = 1; mat1[3][1] = 2; mat1[3][2] = 3; mat1[3][3] = 4; mat1[3][4] = 4; mat1[3][5] = 4; mat1[3][6] = 4; mat1[3][7] = 4; mat1[3][8] = 3; 
  mat1[4][0] = 1; mat1[4][1] = 2; mat1[4][2] = 3; mat1[4][3] = 4; mat1[4][4] = 5; mat1[4][5] = 5; mat1[4][6] = 5; mat1[4][7] = 5; mat1[4][8] = 4; 
  mat1[5][0] = 1; mat1[5][1] = 2; mat1[5][2] = 3; mat1[5][3] = 4; mat1[5][4] = 5; mat1[5][5] = 6; mat1[5][6] = 6; mat1[5][7] = 6; mat1[5][8] = 5; 
  mat1[6][0] = 1; mat1[6][1] = 2; mat1[6][2] = 3; mat1[6][3] = 4; mat1[6][4] = 5; mat1[6][5] = 6; mat1[6][6] = 7; mat1[6][7] = 7; mat1[6][8] = 6; 
  mat1[7][0] = 1; mat1[7][1] = 2; mat1[7][2] = 3; mat1[7][3] = 4; mat1[7][4] = 5; mat1[7][5] = 6; mat1[7][6] = 7; mat1[7][7] = 8; mat1[7][8] = 7; 
  mat1[8][0] = 0; mat1[8][1] = 1; mat1[8][2] = 2; mat1[8][3] = 3; mat1[8][4] = 4; mat1[8][5] = 5; mat1[8][6] = 6; mat1[8][7] = 7; mat1[8][8] = 8; 
  
  //matrix-example (symmetric, positive definite) taken from http://www.ccr.jussieu.fr/ccr/Documentation/Calcul/pessl.html.en_US/html/pessl163.html
  //| 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0 |
  //| 1.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  1.0 |
  //| 1.0  2.0  3.0  3.0  3.0  3.0  3.0  3.0  2.0 |
  //| 1.0  2.0  3.0  4.0  4.0  4.0  4.0  4.0  3.0 |
  //| 1.0  2.0  3.0  4.0  5.0  5.0  5.0  5.0  4.0 |
  //| 1.0  2.0  3.0  4.0  5.0  6.0  6.0  6.0  5.0 |
  //| 1.0  2.0  3.0  4.0  5.0  6.0  7.0  7.0  6.0 |
  //| 1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  7.0 |
  //| 0.0  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0 |
  //

  tmp1 = pMatrixFunc->choleskyDecomposition(mat1, 9, true);
  if (tmp1 != NULL) {
    tmp2 = pMatrixFunc->transpose(tmp1, 9, 9);
    pMatrixFunc->mult(tmp1, tmp2, 9, 9, 9, 9, true);
    diff = pMatrixFunc->sub(mat1, tmp1, 9, 9);
  }

  MFree_2D(tmp1);
  MFree_2D(tmp2);
  MFree_2D(diff);

  //--------------------------- end of some tests
  */
