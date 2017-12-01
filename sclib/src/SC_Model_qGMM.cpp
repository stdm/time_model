/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Model to represent a quasi-GMM as described in		  					*/
/*				'Real-Time Unsupervised Speaker Change Detection',							*/
/*				 L.Lu, H.J.Zhang, 2002 (IEEE)																		*/
/*																																				*/
/*    How it works:																												*/
/*			A GMM-1 is estimated from a initial chunk of trainig data; as			*/
/*			more training-data gets available, the estimation of the cov of		*/
/*			GMM-1 can be imporved; if after updating no improvement is found,	*/
/*			the model evolves to GMM-2, that means, a new cov for a new				*/
/*			mixture is created. this emerges till MAX_MIXTURES mixtures are		*/
/*			reached or no more trainingdata are available.										*/
/*																																				*/
/*			because the model is meant to be used for speaker-change-					*/
/*			detection, where few data is available, it yields methods to			*/
/*			compute the match-score between a model and a cov instead of			*/
/*			a model and some feature-vectors.																	*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 13.03.2004																								*/
/**************************************************************************/

#include <math.h>
#include "SC_Aux.h"
#include "SC_Model_qGMM.h"
#include "SC_TweakableParameters.h"
#include <SV_Error.h>

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_qGMM::SC_Model_qGMM(SC_TweakableParameters* pTweak, unsigned short int dim, unsigned short int maxMixtures) : SC_Model(pTweak) {
	this->complete = false;
	this->dim = dim;
	this->mixtureCount = 0; 
	this->pMatrixFunc = new SC_MatrixFunctions();
  this->pDist = new SC_DistanceMeasures(this->pTweak, this->pMatrixFunc);
  this->maxMixtures = maxMixtures;

	this->Hdr.ModelType = sclib::mtQGMM;

  MArray_1D(this->cov, this->maxMixtures, double**, "SC_Model_qGMM: cov");
  MArray_1D(this->mean, this->maxMixtures, double*, "SC_Model_qGMM: mean");
  MArray_1D(this->weight, this->maxMixtures, double, "SC_Model_qGMM: weights");
  MArray_1D(this->trainingDataCount, this->maxMixtures, unsigned long int, "SC_Model_qGMM: trainingDataCount");

	for (int x = 0; x < this->maxMixtures; x++) {
		this->cov[x] = NULL;
		this->mean[x] = NULL;
		this->weight[x] = 0.0;
    this->trainingDataCount[x] = 0;
	}
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_qGMM::~SC_Model_qGMM() {
  for (int x = 0; x < this->maxMixtures; x++) {
		MFree_2D(this->cov[x]);
		MFree_1D(this->mean[x]);
	}
 
  MFree_0D(this->pDist);
  MFree_0D(this->pMatrixFunc);
  MFree_1D(this->cov);
  MFree_1D(this->mean);
  MFree_1D(this->weight);
  MFree_1D(this->trainingDataCount);
}

//====================================================================================================================
// Generate and update quasi-GMM using only the first element of a linked list of feature-vectors
//====================================================================================================================
int SC_Model_qGMM::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Generate and update quasi-GMM using all linked feature-sets
//====================================================================================================================
int SC_Model_qGMM::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	double **trainingCov, **newModelCov;
	double *trainingMean, *newMean;
  double dist, equalDist;
  unsigned long int N;
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);

	if (pCompleteData->Col != this->dim) {REPORT_ERROR(SVLIB_BadArg, "Wrong dimension of features!");}

	//the easy case: the model is still new, without any mixture-component.
	//so all we have to do is estimate a cov-matrix for the training data and save it as the first mixture-component
	if (this->mixtureCount == 0) {
		this->mean[0] = this->pMatrixFunc->mean(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col); //MeanVec(pData);
    this->cov[0] = this->pMatrixFunc->cov(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col, this->mean[0]); //CovMatrix(pData);		
    this->trainingDataCount[0] = pCompleteData->Row;
		this->weight[0] = 1.0;
		this->mixtureCount = 1;
	}

	//the not-so-easy case: the model already has 1..MAX_MIXTURES components; we have to find out if we shall update
	//the last component or if this one is complete and we create the next one. 
	else if (this->complete != true) {
		trainingMean = this->pMatrixFunc->mean(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col); //MeanVec(pData);	
    trainingCov = this->pMatrixFunc->cov(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col, trainingMean); //CovMatrix(pData);
	
    //C' = N/(N+Nm)*C + Nm/(N+Nm)*Cm; m stands for the new training data 
    newModelCov = combineCov(this->cov[this->mixtureCount-1], trainingCov, this->trainingDataCount[this->mixtureCount-1], pCompleteData->Row, this->dim);

    dist = this->pDist->divergenceShape(this->cov[this->mixtureCount-1], newModelCov, this->dim);
		//sclib::scalarOut("qgmm_dist.txt", dist, this->pTweak);
	  //equalDist = getThreshold(this->cov[this->mixtureCount-1], this->pTweak->modelQgmm.percentDifference);
		//TODO:
		equalDist = 0.05; //std::numeric_limits<double>::max();

		//decide whether to update last mixture or create a new one
    //if dist==-1, a matrix inversion failed, which is an indicator for poorly estimated covs, so update the old one in this case
    //if dist>equalDist, there was a serious improvement in the cov-update, so update it
    if ((dist == -1.0) || (dist > equalDist)) {

			//sclib::scalarOut("updateDist.txt", dist, this->pTweak); //TODO

      MFree_2D(trainingCov);
	    MFree_2D(this->cov[this->mixtureCount-1]);
      this->cov[this->mixtureCount-1] = newModelCov;
      
			//compute new mean vector
			MArray_1D(newMean, this->dim, double, "TrainModel: newMean");
			for (unsigned short int k = 0; k < this->dim; k++) {
				newMean[k] = (double)(this->mean[this->mixtureCount-1][k]*this->trainingDataCount[this->mixtureCount-1] + trainingMean[k]*pCompleteData->Row) / (double)(this->trainingDataCount[this->mixtureCount-1]+pCompleteData->Row);
			}
      MFree_1D(trainingMean);
			MFree_1D(this->mean[this->mixtureCount-1]);
			this->mean[this->mixtureCount-1] = newMean;

      this->trainingDataCount[this->mixtureCount-1] += pCompleteData->Row;

    } else { //create new mixture if there is space left

      if (this->mixtureCount == this->maxMixtures) {
        MFree_1D(trainingMean);
        MFree_2D(trainingCov);
        MFree_2D(newModelCov);
        this->complete = true;
      } else {
        MFree_2D(newModelCov);
        this->cov[this->mixtureCount] = trainingCov;
        this->mean[this->mixtureCount] = trainingMean;
        this->trainingDataCount[this->mixtureCount] = pCompleteData->Row;
        this->mixtureCount++;
      }

    } //dist > threshold

		//renormalize weights
    N = getCompleteTrainingDataCount();
    for (unsigned short int i = 0; i < this->mixtureCount; i++) {
      this->weight[i] = this->trainingDataCount[i] / (double)N;
    }

  }  //complete != true

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  return 0;
}

//====================================================================================================================
// the model-traing-method needs a threshold to decide whether the actual mixture is estimated well enough or not;
// this method computes the divergenceShape-distance between the input cov and a another cov, which's values differ
// from the first one's about percentDiffernce percent. this can serve as the threshold.
// TODO: this is questionable, because it has never been evaluated
//====================================================================================================================
double SC_Model_qGMM::getThreshold(double** cov, double percentDifference) {
	double **differCov = this->pMatrixFunc->mult(cov, (double)(1.0 + percentDifference/100.0), this->dim, this->dim);
	double threshold = this->pDist->divergenceShape(cov, differCov, this->dim);
	
	MFree_2D(differCov);

	return threshold;
}

//====================================================================================================================
// compute the distance between this model and the cov/mean of the TestData
//====================================================================================================================
SV_Data* SC_Model_qGMM::TestModel(SV_Data *TestData) {
  SV_Data *pRes = TestModel(TestData, 1);

  return pRes;
}

//====================================================================================================================
// compute the distance between this model and the cov/mean of the (merged linked list of) TestData
//====================================================================================================================
SV_Data* SC_Model_qGMM::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pRes = new SV_Data(1, 1), *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge);
  double **cov, *mean;
  unsigned long int dataCount = pCompleteData->Row;

  mean = this->pMatrixFunc->mean(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col);
  cov = this->pMatrixFunc->cov(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col, mean);
  pRes->Mat[0][0] = (float)(TestModel(cov, mean, dataCount));

	MFree_2D(cov);
	MFree_1D(mean);
  if (TestData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  return pRes;
}

//====================================================================================================================
// Old: compute the distance between this model and a given cov-matrix formed by testN feature-vectors via deltaBIC
// New: the distance is the weighted divergence shape distance between the test-cov and all modell-component-covs
//      (Lu, Zhang, 2005: Unsupervised Speaker Segmentation and Tracking in Real-Time AudioComponents Analysis)
//====================================================================================================================
double SC_Model_qGMM::TestModel(double** testCov, double* testMean, unsigned long int testDataCount) {
  double dist = 0.0;
	//double **combinedCov;
	
  for (unsigned short int i = 0; i < this->mixtureCount; i++) {
		//Old:		
		//estimate combined cov (combination of all training & testing data)
    //combinedCov = combineCov(this->cov[i], testCov, this->mean[i], testMean, this->trainingDataCount[i], testDataCount, this->dim);
    //combinedCov = combineCov(this->cov[i], testCov, this->trainingDataCount[i], testDataCount, this->dim);
		//compute distance
    //dist += this->weight[i] * this->pDist->BIC(combinedCov, this->cov[i], this->trainingDataCount[i], testCov, testDataCount, this->dim, this->pTweak->modelQgmm.deltaBIClambda);
		//MFree_2D(combinedCov);

		//New:
		dist += this->weight[i] * this->pDist->divergenceShape(this->cov[i], testCov, this->dim);
  }

  return dist;
}

//====================================================================================================================
// return the sum of the n[i]
//====================================================================================================================
unsigned long int SC_Model_qGMM::getCompleteTrainingDataCount(void) {
  unsigned long int N = 0;

  for (unsigned short int i = 0; i < this->mixtureCount; i++) {
    N += this->trainingDataCount[i];
  }

  return N;
}

//====================================================================================================================
// combine a given cov-matrix with new data to the cov resulting from estimating a cov out of the combination of 
// the original and the new training data (using a method by p.m.baggenstoss)
//====================================================================================================================
double** SC_Model_qGMM::combineCov(double** cov1, double** cov2, double* mean1, double* mean2, double weight1, double weight2, unsigned short int dim) {
  unsigned short int d, dd;
  double **combinedCov, *combinedMean, combinedWeight;
  double w1, w2, **c1, **c2, **v1, **v2, **tmp, **q = NULL, **r = NULL, sqrtDim = sqrt((double)dim);

  //compute new weight
  combinedWeight = weight1 + weight2;
  w1 = weight1 / combinedWeight;
  w2 = weight2 / combinedWeight;
  
  //the central mean by weighted average
  MArray_1D(combinedMean, dim, double, "SC_Model_qGMM.combineCov: combinedMean");
  for (d = 0; d < dim; d++) {
    combinedMean[d] = w1*mean1[d] + w2*mean2[d];
  }
  
  //the ROWS of v times the square root of DIM (v is the QR of cov) can be considered 
  //data samples about the means of each distribution.
  c1 = this->pMatrixFunc->choleskyDecomposition(cov1, dim, false);
  c2 = this->pMatrixFunc->choleskyDecomposition(cov2, dim, false);
  v1 = this->pMatrixFunc->mult(cov1, sqrtDim, dim, dim);
  v2 = this->pMatrixFunc->mult(cov2, sqrtDim, dim, dim);
  MFree_2D(c1);
  MFree_2D(c2);

  //they need to be re-referenced to the new center
  for (d = 0; d < dim; d++) {
    for (dd = 0; dd < dim; dd++) {
	    v1[d][dd] = v1[d][dd] + mean1[dd] - combinedMean[dd];
      v2[d][dd] = v2[d][dd] + mean2[dd] - combinedMean[dd];
    }
  }

  //form a weighted augmented matrix 
  this->pMatrixFunc->mult(v1, sqrt(w1), dim, dim, true);
  this->pMatrixFunc->mult(v2, sqrt(w2), dim, dim, true);
  tmp = this->pMatrixFunc->concat(v1, v2, dim, dim, dim, dim);
  MFree_2D(v1);
  MFree_2D(v2);
  this->pMatrixFunc->qrDecomposition(tmp, 2*dim, dim, q, r);
  MFree_2D(tmp);
  MFree_2D(q);
  this->pMatrixFunc->mult(r, 1.0/sqrtDim, dim, dim, true);

  //r is the cholesky factor of the new combined cov
  tmp = this->pMatrixFunc->transpose(r, dim, dim);
  combinedCov = this->pMatrixFunc->mult(tmp, r, dim, dim, dim ,dim);
  MFree_2D(tmp);
  MFree_2D(r);

  return combinedCov;
}

//====================================================================================================================
// combine two covariance-matrices using the method of lu and zhang
//====================================================================================================================
double** SC_Model_qGMM::combineCov(double** cov_1, double** cov_2, unsigned long int n_1, unsigned long int n_2, unsigned short int dim) {
	double **tmp, **combinedCov;

  //C' = N/(N+Nm)*C + Nm/(N+Nm)*Cm; m stands for the new training data 
  tmp = this->pMatrixFunc->mult(cov_2, (double)(n_2)/(double)(n_1+n_2), this->dim, this->dim);
  combinedCov = this->pMatrixFunc->mult(cov_1, (double)(n_1)/(double)(n_1+n_2), this->dim, this->dim);
  combinedCov = this->pMatrixFunc->add(combinedCov, tmp, this->dim, this->dim, true);
  MFree_2D(tmp);

  return combinedCov;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_qGMM::modelOut(ostream& os) {
	unsigned short int x, y, z;

	os << "MixtureNr\tWeight\tn\tMean" << endl;
	for (x = 0; x < this->mixtureCount; x++) {
		os << x << "\t" << this->weight[x] << "\t" << this->trainingDataCount[x] << "\t";
		for (y = 0; y < this->dim; y++) {
			os << this->mean[x][y] << " ";
		}
		os << endl;
	}
	
	os << endl << "Cov's" << endl;

	for (x = 0; x < this->mixtureCount; x++) {
		for (y = 0; y < this->dim; y++) {
			for (z = 0; z < this->dim; z++) {
				os << this->cov[x][y][z];
			}
			os << endl;
		}	
		os << endl;
	}

	os << endl << endl;

	return(os);
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_Model_qGMM::SaveModel(void) {
 	int res, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_qGMM model Failed!");}

  //write mixtureCount, dim, complete
	bytes = io.writeScalar(&(this->DFile), this->mixtureCount);
	bytes += io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->complete);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_qGMM Model Failed!");}

  //write weight and count vectors
	bytes += io.writeArray(&(this->DFile), this->weight, this->mixtureCount);
	bytes += io.writeArray(&(this->DFile), this->trainingDataCount, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_qGMM Model Failed!");}

	//write mean and cov matrices
  bytes += io.writeMatrix(&(this->DFile), this->mean, this->mixtureCount, this->dim);
	for (unsigned int m = 0; m < this->mixtureCount; m++) {
		bytes += io.writeMatrix(&(this->DFile), this->cov[m], this->dim, this->dim);
	}
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_qGMM Model Failed!");}

  return bytes + MHLen; //MHLen may be incorrect...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model * SC_Model_qGMM::LoadModel(void) {
	int	res, x;
	unsigned short int newDim, newMixtureNr;
	bool newStatus;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}
	if (this->Hdr.ModelType != sclib::mtQGMM) {return(NULL);}

	//read mixtureCount, dim, complete
	io.readScalar(&(this->DFile), newMixtureNr, codeSizes, fileSizes);
	if ((DFile.good() != TRUE) || (newMixtureNr == 0)) {return(NULL);}
	io.readScalar(&(this->DFile), newDim, codeSizes, fileSizes);
	if ((DFile.good() != TRUE) || (newDim == 0)) {return(NULL);}
	io.readScalar(&(this->DFile), newStatus, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {return(NULL);}
  
	//prepare this model for changing it's parameters
	if ((newDim != this->dim) || (newMixtureNr != this->mixtureCount)) {
		for (x = 0; x < this->mixtureCount; x++) {
			MFree_1D(this->mean[x]);
			MFree_2D(this->cov[x]);
		}

		this->dim = newDim;
		this->mixtureCount = newMixtureNr;
		this->complete = newStatus;

		for (x = this->mixtureCount-1; x < this->maxMixtures; x++) {
			if (x > 0) {
				this->trainingDataCount[x] = 0;
				this->weight[x]  = 0.0;
			}
		}

		for (x = 0; x < this->mixtureCount; x++) {
			MArray_1D(this->mean[x], this->dim, double, "SC_Model_qGMM: mean");
			MArray_2D(this->cov[x], this->dim, this->dim, double, "SC_Model_qGMM: cov");
		}
	}

	//read weight, count, mean and cov
	io.readArray(&(this->DFile), this->weight, this->mixtureCount, codeSizes, fileSizes);
	io.readArray(&(this->DFile), this->trainingDataCount, this->mixtureCount, codeSizes, fileSizes);
	io.readMatrix(&(this->DFile), this->mean, this->mixtureCount, this->dim, codeSizes, fileSizes);
	for (unsigned int m = 0; m < this->mixtureCount; m++) {
		io.readMatrix(&(this->DFile), this->cov[m], this->dim, this->dim, codeSizes, fileSizes);
	}
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_qGMM Failed!");}

	return(this);
}
