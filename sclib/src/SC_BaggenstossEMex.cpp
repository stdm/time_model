/**************************************************************************/
/*    Inspired by:																												*/
/*      - SC_BaggenstossEM 																								*/
/*																																				*/
/*    Responsibility:																											*/
/*      - this is a variant of SC_BaggenstossEM, which uses only          */
/*        diagonal covariances and is therefore somewhat optimized for    */
/*        not using so much matrix algebra                              	*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 23.02.2005																								*/
/**************************************************************************/

#include <math.h>
#include <assert.h>
#include <time.h>
#include "SC_BaggenstossEMex.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_BaggenstossEMex::SC_BaggenstossEMex(SC_TweakableParameters *pTweak, bool verbose) : SC_BaggenstossEM(pTweak, verbose) {
  this->minVariance = NULL;
  this->variance = NULL;
}

//====================================================================================================================
//	copy-constructor
//====================================================================================================================
SC_BaggenstossEMex::SC_BaggenstossEMex(const SC_BaggenstossEMex& pParent) : SC_BaggenstossEM(pParent) {
	this->pTweak = pParent.pTweak;
  this->pMatrixFunc = new SC_MatrixFunctions();
	this->verbose = pParent.verbose;
	this->dim = pParent.dim;
	this->maxMixtureCount = pParent.maxMixtureCount;
	this->mixtureCount = pParent.mixtureCount;
	this->trainingDataCount = pParent.trainingDataCount;

	this->minStd = NULL;
	this->minVariance = this->pMatrixFunc->copy(pParent.minVariance, this->dim);
	this->mean = this->pMatrixFunc->copy(pParent.mean, this->mixtureCount, this->dim);
	this->choleskyCovar = NULL;
	this->variance = this->pMatrixFunc->copy(pParent.variance, this->mixtureCount, this->dim);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_BaggenstossEMex::~SC_BaggenstossEMex() {
  MFree_0D(this->pMatrixFunc);
  MFree_1D(this->minVariance);
  MFree_1D(this->weight);
  MFree_2D(this->mean);
  MFree_2D(this->variance);
}

//====================================================================================================================
//	C++ translation of P.M.Baggenstoss' function parm = init_gmix(data,n_modes,names,min_std,[rand_init]) for
//  Matlab (03.04.2001). All Matlab-Like functions are located in SC_MatrixFunctions
//
//  Changes to original version:
//    - minStd is now minVariance
//    - minVariance can also be a scalar, not a vector, so that it can be equal for all dimensions
//====================================================================================================================
void SC_BaggenstossEMex::init(SV_Data *pData, unsigned short int mixtureCount, double minVariance, bool randInit) {
  double* minimumVariance;
  MArray_1D(minimumVariance, pData->Col, double, "SC_BaggenstossEMex.init: minimumVariance");

  for (short int d = 0; d < pData->Col; d++) {
    minimumVariance[d] = minVariance;
  }

  init(pData, mixtureCount, minimumVariance, randInit);

  MFree_1D(minimumVariance);
  return;
}

void SC_BaggenstossEMex::init(SV_Data *pData, unsigned short int mixtureCount, double* minVariance, bool randInit) {
  unsigned long int T = pData->Row;
  unsigned long int tryCount;
  unsigned long int *idx = NULL;
  unsigned short int i, d, D = pData->Col;
  unsigned short int minMixtures = 1;
  double maximum, minimum, *mean, *std, *xmin, *xmax, *tmpVec = NULL, *startStd, one = 1.0, tmp;
  bool illConditioned = false, idxHasDuplicates = true, idxIsUnique; 
 
  MFree_1D(this->minVariance);
  this->minVariance = this->pMatrixFunc->copy(minVariance, D);  

  //test if data is ill-conditioned
  if (T > 1) {
    mean = this->pMatrixFunc->mean(pData->Mat, T, D);
    std = this->pMatrixFunc->std(pData->Mat, T, D, mean);
  } else {
    mean = this->pMatrixFunc->scopy(pData->Mat[0], D, (double)(1));
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
  //if ((illConditioned == true) || (max > 100*min)) {REPORT_ERROR(SVLIB_BadData, "SC_BaggenstossEMex.init: Possible ill-conditioning. Your features are shit. You should scale them or remove means!");}

  //compute initial value of covariances
  if (T > 1) {
    xmin = this->pMatrixFunc->min(pData->Mat, T, D);
    xmax = this->pMatrixFunc->max(pData->Mat, T, D);
  } else {
    xmin = this->pMatrixFunc->scopy(pData->Mat[0], D, (double)(1));
    xmax = this->pMatrixFunc->scopy(pData->Mat[0], D, (double)(1));
  }
  MArray_1D(startStd, D, double, "SC_baggenstossEMex.init: startStd");
  for (d = 0; d < D; d++) {
    startStd[d] = sclib::max((xmax[d]-xmin[d])/4.0, sqrt(minVariance[d]));
  }
  MFree_1D(xmin);
  MFree_1D(xmax);

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
    //if (idxHasDuplicates == true) {REPORT_ERROR(SVLIB_BadData, "SC_BaggenstossEMex.init: mixtureCount has been reduced to less than allowed, you need more data!");}
  } else { //no random init
    mixtureCount = sclib::min(mixtureCount, T);
    MArray_1D(idx, mixtureCount, unsigned long int, "SC_BaggenstossEMex.init: idx");
    for (i = 0; i < mixtureCount; i++) {
      idx[i] = i + 1;
    }
  } 
  MFree_1D(tmpVec);

  //set up class-parameters with the computed starting-covars, random means and equal weights
  MFree_1D(this->weight);
  MFree_2D(this->mean);
  MFree_2D(this->variance);
  this->mixtureCount = mixtureCount;
  this->trainingDataCount = T;
  this->dim = D;
  MArray_1D(this->weight, this->mixtureCount, double, "SC_BaggenstossEMex.init: weight");
  MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_BaggenstossEMex.init: mean");
  MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_BaggenstossEMex.init: variance");

  for (i = 0; i < this->mixtureCount; i++) {
    this->weight[i] = 1.0 / (double)(this->mixtureCount);
    for (d = 0; d < this->dim; d++) {
      this->mean[i][d] = pData->Mat[idx[i]-1][d];
      this->variance[i][d] = startStd[d] * startStd[d];
    }
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
double SC_BaggenstossEMex::train(SV_Data* pData, unsigned int &actualIterations, unsigned long int maxIterations, long int samplesPerMode, bool bias, double maxCloseness, bool addModes, double splitThresh) {
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
    if (nModesOld != this->mixtureCount) {convergenceCounter = 0;}

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
			sclib::classOut("baggEx.txt", this, this->pTweak);
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
//    - uses a variance-vector (=diagonal covariance-matrix) instead of the cholesky-factor of the covar
//    - no individual data-weights for each input feature vector
//====================================================================================================================
double SC_BaggenstossEMex::emStep(SV_Data *pData, bool bias) {
  unsigned long int t, T = pData->Row;
  unsigned short int d, D = this->dim, i, I = this->mixtureCount;
  float **z = pData->Mat;
  double mx, px, norm, **w, *alpha, Q, tmpQ, **sumWz, **sumWz2;

  assert(D == pData->Col);
  
  //E-Step
  MArray_2D(w, (short)(I), (long)(T), double, "SC_BaggenstossEMex.emStep: w");
  alpha = this->pMatrixFunc->zeros(I);
  sumWz = this->pMatrixFunc->zeros(I, D);
  sumWz2 = this->pMatrixFunc->zeros(I, D);
  
  Q = 0.0;
  for (t = 0; t < T; t++) {
    mx = -1.0 * numeric_limits<double>::max();
    norm = 0.0;
    tmpQ = 0.0;

    for (i = 0; i < I; i++) { //compute the ingredients for the final weights w and priors alpha
      px = fullLogGaussian(z[t], this->mean[i], this->variance[i], this->dim);
      w[i][t] = log(this->weight[i]) + px;
      tmpQ += exp(w[i][t]);
      mx = sclib::max(mx, w[i][t]);
    }

    for (i = 0; i < I; i++) {
      w[i][t] = exp(w[i][t] - mx);
      norm += w[i][t];
    }

    Q += log(tmpQ); //compute the log-likelihood Q

    for (i = 0; i < I; i++) { //compute the final (properly normalized) w's and alphas in the log-domain
      w[i][t] /= norm;
      alpha[i] += w[i][t];
      
      for (d = 0; d < D; d++) {
        sumWz[i][d] += w[i][t] * z[t][d];
        sumWz2[i][d] += w[i][t] * z[t][d] * z[t][d];
      }
    }
  }

  //M-Step
  double sumAlpha = this->pMatrixFunc->sum(alpha, I);

  //Update means, weights, covars and constrain covariance
  for (i = 0; i < I; i++) {
    
    //update mean and variance
    for (d = 0; d < D; d++) {
      this->mean[i][d] = sumWz[i][d] / alpha[i]; 
      this->variance[i][d] = (sumWz2[i][d] / alpha[i]) - (this->mean[i][d] * this->mean[i][d]);
    
      //do covar-trimming
      if (bias == false) { //constraint method
        this->variance[i][d] = sclib::max(this->variance[i][d], this->minVariance[d]);
      } else {
        this->variance[i][d] += this->minVariance[d];
      }
    }

    //update weights
    this->weight[i] = sclib::max(alpha[i] / sumAlpha, 0.0);
  }

  MFree_1D(alpha);
  MFree_2D(sumWz);
  MFree_2D(sumWz2);
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
double SC_BaggenstossEMex::test(SV_Data *pData) {
  double px, q, Q = 0.0;

  for (unsigned long int t = 0; t < (unsigned long)pData->Col; t++) {
    q = 0.0;

    for (unsigned short int i = 0; i < this->mixtureCount; i++) {
      px = fullLogGaussian(pData->Mat[t], this->mean[i], this->variance[i], this->dim);
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
//	Evaluates the log-gaussian PDF for the feature-vector x with mean-vector mean and variance-vector variance
//  Adopted from P.M.Baggenstoss' Matlab-version function l = lqr_eval(x,means,cholesky_covar)
//
//  Changes to original version:
//    - uses a variance-vector (=diagonal covariance-matrix) instead of the cholesky-factor of the covar
//    - runs for just one feature vector instead of all at once
//====================================================================================================================
double SC_BaggenstossEMex::fullLogGaussian(double* x, double* mean, double* variance, unsigned short dim) {
  double sd, t, l = 0.0;

  for (unsigned short int d = 0; d < dim; d++) {
    sd = sqrt(variance[d]);
    t = (x[d] - mean[d]) / sd;
		l += log(1.0/(sclib::sqrt_2pi*sd)) + (-0.5 * t * t);
  }

  return l;
}

//====================================================================================================================
// This is a altered version of the above function: It takes the data in a float-vector instead of a double-vector
//====================================================================================================================
double SC_BaggenstossEMex::fullLogGaussian(float* x, double* mean, double* variance, unsigned short dim) {
  double sd, t, l = 0.0;

  for (unsigned short int d = 0; d < dim; d++) {
    sd = sqrt(variance[d]);
    t = (x[d] - mean[d]) / sd;
		l += log(1.0/(sclib::sqrt_2pi*sd)) + (-0.5 * t * t);
  }

  return l;
}

//====================================================================================================================
//	Kills weak mixtures according to their weight. All mixtures with weight<minWeightAll are killed, but only one 
//  with weight<minWeightOne.
//  Adopted from P.M.Baggenstoss' Matlab-version function parm = gmix_deflate(parm,min_weight_1,min_weight_all)
//
//  Changes to original version:
//    - uses a variance-vector (=diagonal covariance-matrix) instead of the cholesky-factor of the covar
//    - returns the count of killed mixtures
//====================================================================================================================
unsigned short int SC_BaggenstossEMex::deflate(double minWeightOne, double minWeightAll) {
  unsigned long int d, i, j, minWeightIdx;
  unsigned short int newMixtureCount, deflateCount = 0;
  double minWeight, sumWeights;
  double *newWeight, **newMean, **newVariance;
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
          MArray_2D(newMean, (long)newMixtureCount, this->dim, double, "SC_BaggenstossEMex.delate: newMean");
          MArray_2D(newVariance, (long)newMixtureCount, this->dim, double, "SC_BaggenstossEMex.delate: newVariance");
          for (i = 0, j = 0; i < this->mixtureCount; i++) {
            if (i != minWeightIdx) {
              newWeight[j] = this->weight[i] / sumWeights;
              for (d = 0; d < this->dim; d++) {
                newMean[j][d] = this->mean[i][d];
                newVariance[j][d] = this->variance[i][d];
              }
              j++;
            }
          }

          //make the new parameters the actual ones
          MFree_1D(this->weight);
          MFree_2D(this->mean);
          MFree_2D(this->variance);
          this->weight = newWeight;
          this->mean = newMean;
          this->variance = newVariance;
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
//    - uses a variance-vector (=diagonal covariance-matrix) instead of the cholesky-factor of the covar
//====================================================================================================================
double SC_BaggenstossEMex::mixtureCloseness(double* mean1, double* mean2, double* variance1, double* variance2, unsigned short int dim) {
  unsigned long int d, ev, i, t = 0, T = 2 * (2*dim + 1);
  double **x, *variance, *mean;
  double p1x1 = 0.0, p1x2 = 0.0, p2x1 = 0.0, p2x2 = 0.0;

  MArray_2D(x, (long)T, (long)dim, double, "SC_baggenstossEM.mixtureCloseness: x");

  //loop over the two mixtures and synthesize T samples
  for (i = 0; i < 2; i++) {

    //select the right means and variances
    if (i==0) {
      variance = variance1;
      mean = mean1; 
    } else { 
      variance = variance2; 
      mean = mean2; 
    }

    //add the center point 
    this->pMatrixFunc->setRow(x, mean, T, dim, t++);

    //loop over eigenvectors (which are, in the case of diagonal covar, just a vector full of zeros with the variance element in the right position)
    for (ev = 0; ev < dim; ev++) {
      for (d = 0; d < dim; d++) {
        x[t][d] = mean[d] + ((ev==d)?sqrt(variance[d]):0.0);
        x[t+1][d] = mean[d] - ((ev==d)?sqrt(variance[d]):0.0);
      }
      t = t + 2;
    }
  }
  
  //compute the cross-likelihoods
  for (t = 0; t < T; t++) {
    if (t > T/2) {
      p1x2 += fullLogGaussian(x[t], mean1, variance1, dim);
      p2x2 += fullLogGaussian(x[t], mean2, variance2, dim);
    } else {
      p2x1 += fullLogGaussian(x[t], mean2, variance2, dim);
      p1x1 += fullLogGaussian(x[t], mean1, variance1, dim);
    }
  }

  MFree_2D(x);

  return (p1x2 + p2x1) - (p2x2 + p1x1);
}

//====================================================================================================================
//	Merges mixtures closer than the threshold maxXloseness, which should be -1*dim. Identical mixtures have 
//  closeness zero.
//  Adopted from P.M.Baggenstoss' Matlab-version function parm = gmix_merge(parm,max_closeness)
//
//  Changes to original version:
//    - uses a variance-vector (=diagonal covariance-matrix) instead of the cholesky-factor of the covar
//    - returns the count of merged mixtures
//====================================================================================================================
unsigned short int SC_BaggenstossEMex::merge(double maxCloseness) {
  unsigned short int d, dd, i, j, mergeCounter = 0;
  double closeness;
  double newWeight, *newMean, *newVariance;
  double w1, w2, *m1, *m2, **C, sqrtDim = sqrt((double)this->dim);

  if (this->mixtureCount > 1) {            

    MArray_1D(newMean, this->dim, double, "SC_BaggenstossEMex.merge: newMean");
    MArray_1D(newVariance, this->dim, double, "SC_BaggenstossEMex.merge: newVariance");
    MArray_1D(m1, this->dim, double, "SC_BaggenstossEMex.merge: m1");
    MArray_1D(m2, this->dim, double, "SC_BaggenstossEMex.merge: m2");
    MArray_2D(C, this->dim*2, this->dim, double, "SC_BaggenstossEMex.merge: C");
            
    for (i = 0; i < this->mixtureCount; i++) {
      for (j = i+1; j < this->mixtureCount; j++) {
        if ((this->weight[i] > 0) && (this->weight[j] > 0)) {
          closeness = mixtureCloseness(this->mean[i], this->mean[j], this->variance[i], this->variance[j], this->dim);
          if ((closeness > maxCloseness) && (mergeCounter < (this->mixtureCount/2 +1))) {
            
            mergeCounter++;
            
            //compute new weight
            newWeight = this->weight[i] + this->weight[j];
            w1 = this->weight[i] / newWeight;
            w2 = this->weight[j] / newWeight;
            
            //the central mean by weighted average
            for (d = 0; d < this->dim; d++) {
              newMean[d] = w1*this->mean[i][d] + w2*this->mean[j][d];
              m1[d] = this->mean[i][d] - newMean[d];
              m2[d] = this->mean[j][d] - newMean[d];
            }

            //build the matrix C of samples of each distribution 1 and 2: it's variance is the new variance
            //each sample is build as sqrt(dim)*sqrt(variance), where variance is the diagonal covar-matrix
            //the samples are then re-referenced to the new mean and weighted by their weights
            for (dd = 0; dd < this->dim*2; dd++) {
              for (d = 0; d < this->dim; d++) {
                if (dd < this->dim) { //upper half with samples of distribution 1
                  if (dd == d) { //a diagonal element
                    C[dd][d] = ((sqrtDim * sqrt(this->variance[i][d])) + m1[d]) * w1;
                  } else { //a non-diagonal element
                    C[dd][d] = m1[d] * w1;
                  }
                } else { //lower half with samples of distribution 2
                  if (dd-this->dim == d) { //a diagonal element
                    C[dd][d] = ((sqrtDim * sqrt(this->variance[j][d])) + m2[d]) * w2;
                  } else { //a non-diagonal element
                    C[dd][d] = m2[d] * w2;
                  }
                }
              }
            }

            //compute the variance of the new matrix C
						MFree_1D(newVariance);
            newVariance = this->pMatrixFunc->variance(C, 2*this->dim, this->dim);
        
            //update parameters
            if (this->weight[i] > this->weight[j]) { //merging mixture j into mixture i
        			this->weight[i] = newWeight;
              this->weight[j] = 0.0;
              this->pMatrixFunc->setRow(this->variance, newVariance, this->maxMixtureCount, this->dim, i);
              this->pMatrixFunc->setRow(this->mean, newMean, this->mixtureCount, this->dim, i);
            } else { //merging mixture i into mixture j
        			this->weight[j] = newWeight;
              this->weight[i] = 0.0;
              this->pMatrixFunc->setRow(this->variance, newVariance, this->maxMixtureCount, this->dim, j);
              this->pMatrixFunc->setRow(this->mean, newMean, this->mixtureCount, this->dim, j);
            }

          } //really merge this mixtures
        } //weights > 0
        if (this->weight[i] == 0) {break;}
      } //for j
    } //for i

    MFree_1D(newMean);
    MFree_1D(newVariance);
    MFree_1D(m1);
    MFree_1D(m2);
    MFree_2D(C);

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
//    - uses a variance-vector (=diagonal covariance-matrix) instead of the cholesky-factor of the covar
//    - returns the count of splitted mixtures
//====================================================================================================================
unsigned short int SC_BaggenstossEMex::split(SV_Data *pData, double threshold) {
  unsigned short int dd, d, i, I = this->mixtureCount, D = this->dim, splitCounter = 0;
  unsigned long int t, T = pData->Row;
  float **z = pData->Mat;
  double mx, px, norm, **w;
  double *sumW;

  MArray_2D(w, (short)(I), (long)(T), double, "SC_BaggenstossEMex.emStep: w");
  sumW = this->pMatrixFunc->zeros(I);
  
  //compute membership probabilities of modes for each sample, store it in w[T][I]
  for (t = 0; t < T; t++) {
    mx = -1.0 * numeric_limits<double>::max();
    norm = 0.0;

    for (i = 0; i < I; i++) { //compute the ingredients for the final weights w
      px = fullLogGaussian(z[t], this->mean[i], this->variance[i], this->dim);
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
  
  double k1, k2, k3, k4, tmp, fac, merit, mMax, *mAll;
  unsigned short int dMax, *dAll;
  MArray_1D(dAll, I, unsigned short int, "SC_BaggenstossEMex.split: dAll");
  MArray_1D(mAll, I, double, "SC_BaggenstossEMex.split: mAll");

  //check which mixtures need splitting
  for (i = 0; i < I; i++) {
    
    //for each dimension, project onto this vector
    mMax = -1.0 * numeric_limits<double>::max();
    for (d = 0; d < D; d++) {
      k1 = 0.0; k2 = 0.0; k3 = 0.0; k4 = 0.0;

      for (t = 0; t < T; t++) {
        //tmp = zt * V(:,j) / S(j,j); zt=z-mean
        tmp = 0;
        for (dd = 0; dd < D; dd++) {
          tmp += (z[t][dd] - this->mean[i][dd]) * (dd==d)?sqrt(this->variance[i][dd]):0.0;
        }
     
        k1 += w[i][t] * pow(tmp, 1.0);
        k2 += w[i][t] * pow(tmp, 2.0);
        k3 += w[i][t] * pow(tmp, 3.0);
        k4 += w[i][t] * pow(tmp, 4.0);
      }

      k1 /= sumW[i]; 
      k2 /= sumW[i]; 
      k3 /= sumW[i]; 
      k4 /= sumW[i];

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

    mAll[i]= mMax;
    dAll[i] = dMax;
  } //for i

	MFree_1D(sumW);
	MFree_2D(w);
	

  double *m1, *m2, *newWeight, **newMean, **newVariance;
  unsigned short int ii, newMixtureCount;

  //do the splitting
  for (i = 0; i < I; i++) { //yes, 'I' is the old mixtureCount!!!
    if (mAll[i] > threshold) {
      if (this->mixtureCount < this->maxMixtureCount) {
        dMax = dAll[i];
        
        //compute new means;
        MArray_1D(m1, D, double, "SC_BaggenstossEMex.split: m1");
        MArray_1D(m2, D, double, "SC_BaggenstossEMex.split: m2");
        for (d = 0; d < D; d++) {
          m1[d] = this->mean[i][d] + ((d==dMax)?sqrt(this->variance[i][dMax]):0.0);//(V[d][dMax] * (S[dMax][dMax] / 2.0));
          m2[d] = this->mean[i][d] - ((d==dMax)?sqrt(this->variance[i][dMax]):0.0);//(V[d][dMax] * (S[dMax][dMax] / 2.0));
        }

        //compute parameters
        newMixtureCount = this->mixtureCount + 1;
        MArray_1D(newWeight, (long)newMixtureCount, double, "SC_BaggesnstossEM.deflate: newWeight");
        MArray_2D(newMean, (long)newMixtureCount, this->dim, double, "SC_BaggenstossEMex.delate: newMean");
        MArray_2D(newVariance, (long)newMixtureCount, this->dim, double, "SC_BaggenstossEMex.delate: newVariance");
        for (ii = 0; ii < newMixtureCount; ii++) {
          if ((ii < this->mixtureCount) && (ii != i)) { //just copy old parameters
            newWeight[ii] = this->weight[ii];
            for (d = 0; d < D; d++) {
              newMean[ii][d] = this->mean[ii][d];
              newVariance[ii][d] = this->variance[ii][d];
            }
          } else if ((ii < this->mixtureCount) && (ii == i)) { //alter the splitted mixture
            newWeight[ii] = this->weight[i] / 2.0;
            for (d = 0; d < D; d++) {
              newVariance[ii][d] = this->variance[ii][d];
              newMean[ii][d] = m1[d];
            }
          } else if (ii >= this->mixtureCount) { //generate the new mixture
            newWeight[ii] = this->weight[i] / 2.0;
            for (d = 0; d < D; d++) {
              newVariance[ii][d] = this->variance[i][d];
              newMean[ii][d] = m2[d];
            }
          }
        }
        MFree_1D(m1);
        MFree_1D(m2);

        //set new parameters
        MFree_1D(this->weight);
        MFree_2D(this->mean);
        MFree_2D(this->variance);
        this->weight = newWeight;
        this->mean = newMean;
        this->variance = newVariance;
        this->mixtureCount = newMixtureCount;
        newWeight = NULL;
        newMean = NULL;
        newVariance = NULL;

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
//  Computes and returns the covariance-matrix out of the stored variance-vector
//====================================================================================================================
double** SC_BaggenstossEMex::getCovar(unsigned short int mixture) {
  return this->pMatrixFunc->diag(this->variance[mixture], this->dim);
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& operator<< (ostream& OutS, SC_BaggenstossEMex& Data) {
	unsigned short int x, y;

	OutS.setf(ios::fixed|ios::basefield);
  OutS.precision(5);

  OutS << "MixNr\tWeight\tMean" << endl;
	for (x = 0; x < Data.mixtureCount; x++) {
    OutS << x << "\t" << setw(7) << Data.weight[x] << "\t";
		for (y = 0; y < Data.dim; y++) {
			OutS << setw(10) << Data.mean[x][y] << ((y == Data.dim-1) ? "" : " ");
    }
		OutS << endl;
	}

  OutS << endl << "MixNr\t\tVariance" << endl;
	for (x = 0; x < Data.mixtureCount; x++) {
  	OutS << x << "\t\t";
    for (y = 0; y < Data.dim; y++) {
      OutS << setw(10) << Data.variance[x][y] << ((y == Data.dim-1) ? "" : " ");
    }
		OutS << endl;
	}

	OutS << endl;

	return(OutS);
}

//====================================================================================================================
// write parameters to an already opened, ready-to-write binary file
//====================================================================================================================
unsigned int SC_BaggenstossEMex::write(fstream *file) {
	unsigned int bytes = 0;
	SV_DataIO io;

	bytes = io.writeScalar(file, this->verbose);
	bytes += io.writeScalar(file, this->maxMixtureCount);
	bytes += io.writeScalar(file, this->mixtureCount);
	bytes += io.writeScalar(file, this->dim);
	bytes += io.writeScalar(file, this->trainingDataCount);
	bytes += io.writeArray(file, this->minVariance, this->dim);
	bytes += io.writeArray(file, this->weight, this->mixtureCount);
	bytes += io.writeMatrix(file, this->mean, this->mixtureCount, this->dim);
	bytes += io.writeMatrix(file, this->variance, this->mixtureCount, this->dim);

	return bytes;
}

//====================================================================================================================
// read parameters from an already opened, correctly positioned binary file
//====================================================================================================================
unsigned int SC_BaggenstossEMex::read(fstream *file, SV_DataIO::SV_DatatypeSizes *fileSizes) {
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes;
	io.getCurrentDatatypeSizes(codeSizes);

	io.readScalar(file, this->verbose, codeSizes, *fileSizes);
	io.readScalar(file, this->maxMixtureCount, codeSizes, *fileSizes);
	io.readScalar(file, this->mixtureCount, codeSizes, *fileSizes);
	io.readScalar(file, this->dim, codeSizes, *fileSizes);
	io.readScalar(file, this->trainingDataCount, codeSizes, *fileSizes);
	io.readArray(file, this->minVariance, this->dim, codeSizes, *fileSizes);
	io.readArray(file, this->weight, this->mixtureCount, codeSizes, *fileSizes);
	io.readMatrix(file, this->mean, this->mixtureCount, this->dim, codeSizes, *fileSizes);
	io.readMatrix(file, this->variance, this->mixtureCount, this->dim, codeSizes, *fileSizes);

	return (file->good() != TRUE) ? SVLIB_Fail : SVLIB_Ok;
}
