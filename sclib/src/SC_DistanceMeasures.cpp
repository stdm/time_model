/**************************************************************************/
/*	A class to hold the implementations of all used diatnce-measures			*/
/*	throughout the SC_* project																						*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.04.2005																								*/
/**************************************************************************/

#include <math.h>
#include <assert.h>
#include "SC_DistanceMeasures.h"
#include "SC_FeatureHandler.h"
#include "SC_Aux.h"
#include <SV_Error.h>
#include "SC_Centroid_Point.h"
#include "SC_MixtureModel.h"

//====================================================================================================================
//	constructors
//====================================================================================================================
SC_DistanceMeasures::SC_DistanceMeasures(SC_TweakableParameters* pTweak, SC_MatrixFunctions* pMatrixFunc, bool verbose) {
	this->pTweak = pTweak;
  this->verbose = verbose;
  this->pMatrixFunc = pMatrixFunc;
  this->pModelHandler = new SC_ModelHandler(this->pTweak, this->verbose);
  this->linkedMatrixFunc = (this->pMatrixFunc != NULL) ? true : false;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_DistanceMeasures::~SC_DistanceMeasures() {
  MFree_0D(this->pModelHandler);
  if (this->linkedMatrixFunc == false) {
    MFree_0D(this->pMatrixFunc);
  }
}

//====================================================================================================================
//	access to protected member
//====================================================================================================================
SC_MatrixFunctions* SC_DistanceMeasures::getMatrixFunc(void) {
	if (this->pMatrixFunc == NULL) {
		this->pMatrixFunc = new SC_MatrixFunctions();
		this->linkedMatrixFunc = false;
	}

	return this->pMatrixFunc;
}

//====================================================================================================================
//	Cross-Likelihood-Ratio distance measure between this and the pSecond cluster
//  Also known as Cross-Entropy or Symmetric Kullback-Liebler Information Distance
//	According to: 'Clustering Speakers by Their Voices', Solomonoff, Mielke, Gish, 1998 (IEEE ICASSP)
//                'Unsupervised Speaker Indexing Using Speaker Model Selection Based On BIC', Nishida, Kawahara, 
//								 2003 (IEEE ICASSP)
//								'Blind Clustering of speech utterances based on speaker and language characteristics', Reynolds,
//								 Singer, Carlson, McLaughlin, O'Leary, Zissman, 1998 (IEEE ICSLP)
//====================================================================================================================
double SC_DistanceMeasures::CLR(SC_Model* pModel_1, SC_Model* pBkModels_1, SV_Data* pData_1, unsigned long int segmentCount_1, SC_Model* pModel_2, SC_Model* pBkModels_2, SV_Data* pData_2, unsigned long int segmentCount_2) {
  double log_P_X_lambdaX, log_P_X_lambdaY, log_P_Y_lambdaY, log_P_Y_lambdaX;
	double	dist;	//dist = log(P(X|lambdaX) / P(X|lambdaY)) + log(P(Y|lambdaY) / P(Y|lambdaX))

  log_P_X_lambdaX = this->pModelHandler->testModel(pModel_1, pData_1, segmentCount_1, false, pBkModels_1, false);
	log_P_X_lambdaY = this->pModelHandler->testModel(pModel_2, pData_1, segmentCount_1, false, pBkModels_1, false);
	log_P_Y_lambdaX = this->pModelHandler->testModel(pModel_1, pData_2, segmentCount_2, false, pBkModels_2, false);
	log_P_Y_lambdaY = this->pModelHandler->testModel(pModel_2, pData_2, segmentCount_2, false, pBkModels_2, false);

	dist = (log_P_X_lambdaX - log_P_X_lambdaY) + (log_P_Y_lambdaY - log_P_Y_lambdaX);
  assert(sclib::isFinite(dist));

	return dist;
}

//====================================================================================================================
//	Generalized Likelihood Ratio distance measure between this and the second cluster
//	According to: 'Clustering Speakers by Their Voices', Solomonoff, Mielke, Gish, 1998 (IEEE ICASSP)
//                'Online Speaker Clustering', Liu, Kubala, 2004 (IEEE ICASSP)
//								'Segregation Of Speakers For Speech Recognition And Speaker Identification', Gish, Schmidt, 1991
//								 (IEEE ICASSP)
//====================================================================================================================
double SC_DistanceMeasures::GLR(SC_Model* pModel_1, SC_Model* pBkModels_1, SV_Data* pData_1, unsigned long int segmentCount_1, SC_Model* pModel_2, SC_Model* pBkModels_2, SV_Data* pData_2, unsigned long int segmentCount_2) {
	double dist;	//dist = log(P(X) / (P(x1) * P(x2)))
	double log_P_X, log_P_x1, log_P_x2;
  SV_Data *pFeatureHook, *pFeatureSave;
	SC_Model *pCompleteModel, *pModelHook, *pModelSave;

  log_P_x1 = this->pModelHandler->testModel(pModel_1, pData_1, segmentCount_1, false, pBkModels_1, false);
  log_P_x2 = this->pModelHandler->testModel(pModel_2, pData_2, segmentCount_2, false, pBkModels_2, false);

	pFeatureHook = sclib::getListWithIndex(pData_1, segmentCount_1-1);
	pFeatureSave = pFeatureHook->Next; //store old link
	pFeatureHook->Next = pData_2;
	pModelHook = sclib::getListWithIndex(pBkModels_1, segmentCount_1-1);
	pModelSave = (SC_Model*)(pModelHook->Next);
	pModelHook->Next = pBkModels_2;

  switch (this->pTweak->distanceMeasure.mergeMode) {
    case sclib::mergeAddUp: {
      pCompleteModel = this->pModelHandler->combineModels(pModel_1, pModel_2);
      break;
    }
    case sclib::mergeRetrain: {
      pCompleteModel = this->pModelHandler->combineModels(pModel_1, pModel_2, pData_1, segmentCount_1+segmentCount_2, pBkModels_1);
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Specified model merge mode unknown");
      break;
    }
  }
  
  log_P_X = this->pModelHandler->testModel(pCompleteModel, pData_1, segmentCount_1+segmentCount_2, false, pBkModels_1, false);
  pFeatureHook->Next = pFeatureSave;
	pModelHook->Next = pModelSave;

  //this is -GLR as defined in 'Clustering Speakers by Their Voices', because standard GLR falls with rising similarity, which is the opposite (or negative) of a distance!
  dist = log_P_x1 + log_P_x2 - log_P_X;
  assert(sclib::isFinite(dist));

	MFree_0D(pCompleteModel);

	return dist;
}

//====================================================================================================================
//	if 'pComleteModel' comes from the merging of 'pModel_1' and 'pModel_2', the Bayesian Information Criterion shows 
//  if this merge is to prefer compared with the 2 single models (if deltaBIC > 0)
//	According to: 'A Robust Speaker Clustering Algorithm', Ajmera, Wooters, 2003
//								'Speaker, Environment And Channel Change Detection And Clustering Via The BIC', Chen, 
//								 Gopalakrishnan, 1998
//  Needs (as parameters) for all 3 models to test against each other:
//    - the model
//    - (linked list of) background-models
//    - (linked list of) training feature vectors
//    - count of list-elements to merge 
// alpha is to (de-)emphasize the penalty-factor for model-complexity; it should normally be =1.0 but is zero here by 
// default for historical reasons
//====================================================================================================================
double SC_DistanceMeasures::BIC(SC_Model* pCompleteModel, SC_Model* pCompleteBkModels, SV_Data* pCompleteData, unsigned long int completeSegmentCount, SC_Model* pModel_1, SC_Model* pBkModels_1, SV_Data* pData_1, unsigned long int segmentCount_1, SC_Model* pModel_2, SC_Model* pBkModels_2, SV_Data* pData_2, unsigned long int segmentCount_2, double alpha) {
  SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);
  unsigned long int completeT = pFeatureHandler->countFeaturesInList(pCompleteData, completeSegmentCount);
  double log_P_X_lambda, log_P_x1_lambda1, log_P_x2_lambda2;
  double penaltyFactor;
	double	res; //deltaBIC = log(L(X|M)) - (log(L(x1|m1)) + log(L(x2|m2))) - penaltyFactor
               //deltaBIC >= 0 => choose M instead of m1&m2

  //compute log-likelihoods; compensate for the division by the number of features so that the penalty-factor below works as it should work
  log_P_X_lambda = this->pModelHandler->testModel(pCompleteModel, pCompleteData, completeSegmentCount, false, pCompleteBkModels, false);
  log_P_x1_lambda1 = this->pModelHandler->testModel(pModel_1, pData_1, segmentCount_1, false, pBkModels_1, false);
  log_P_x2_lambda2 = this->pModelHandler->testModel(pModel_2, pData_2, segmentCount_2, false, pBkModels_2, false);

  penaltyFactor = alpha/2.0;
  penaltyFactor *= log((double)(completeT));
  penaltyFactor *= (double)(pCompleteModel->getFreeParameterCount()) - ((double)(pModel_1->getFreeParameterCount()) + (double)(pModel_2->getFreeParameterCount()));

  //this is similar to GLR (yes, the positive one, not like if GLR is used as a distance-measure!) when the penalty-factor is 0
  res = log_P_X_lambda - log_P_x1_lambda1 - log_P_x2_lambda2 - penaltyFactor;  //TODO -> penaltyFactor can get very large with GMM-UBM!!!

  assert(sclib::isFinite(res));
  MFree_0D(pFeatureHandler);

	return res; 
}

//====================================================================================================================
//	function to calcualte the delta-BIC between the cov-matrices of speech-segments of (maybe) different size
//
//	if return-value > 0, the segments differ significantly, so one can assume that they are uttered by different
//	speakers
//
//	according to the following paper:
//	'real-time unsupervised speaker change detection', l.lu, h.j.zhang, 2002
//
//	the formula in short:
//	BIC(i,j) =   0.5 * ( (n_1+n_2)*log|cov| - n_1*log|cov_1| - n_2*log|cov_2| )
//             - 0.5 * lambda * (dim + 0.5*dim*(dim+1)) * log(n_1+n_2)
//====================================================================================================================
double SC_DistanceMeasures::BIC(double** completeCov, double** cov_1, unsigned long int trainingDataCount_1, double** cov_2, unsigned long int trainingDataCount_2, unsigned short int dim, double lambda) {
	double deltaBIC = 0.0;
	double logDetC1, logDetC2, logDetC;
  unsigned long int n_1 = trainingDataCount_1, n_2 = trainingDataCount_2;
	unsigned long int n = n_1 + n_2;

  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

	logDetC1 = this->pMatrixFunc->logDet(cov_1, dim); //logDet() works better when cov is nearly singular so that det=0 would be returned
	logDetC2 = this->pMatrixFunc->logDet(cov_2, dim);
  logDetC = this->pMatrixFunc->logDet(completeCov, dim);

	//formula is nearly identical to "'Speaker, Environment And Channel Change Detection And Clustering 
	//Via The BIC', Chen, Gopalakrishnan, 1998" implemented below: there, the first term isn't multiplied with 0.5
	deltaBIC = 0.5 * ( (double)(n)*logDetC - (double)(n_1)*logDetC1 - (double)(n_2)*logDetC2 )
	         - 0.5 * lambda * ((double)(dim) + 0.5*(double)(dim)*(double)(dim+1)) * log((double)n);

	return deltaBIC;
}
double SC_DistanceMeasures::BIC(SV_Data *pFirst, SV_Data *pSecond) {
	double **cov1 = this->getMatrixFunc()->cov(pFirst->Mat, pFirst->Row, pFirst->Col);
	double **cov2 = this->getMatrixFunc()->cov(pSecond->Mat, pSecond->Row, pSecond->Col);

	SV_Data *pSave = pFirst->Next;
	pFirst->Next = pSecond;
	SV_Data *pMerged = pFirst->MergeData(2);
	double **cov = this->getMatrixFunc()->cov(pMerged->Mat, pMerged->Row, pMerged->Col);
	MFree_0D(pMerged);
	pFirst->Next = pSave;

	double bic = BIC(cov, cov1, pFirst->Row, cov2, pSecond->Row, pFirst->Col, 1.0);

	MFree_2D(cov);
	MFree_2D(cov1);
	MFree_2D(cov2);

	return bic;
}

//====================================================================================================================
//	The standard BIC for one concrete Model and it's test-result:
//  
//  BIC(model, data) = -2*log(likelihood(data|model)) + #(model)*log(#(data))
//
//  where #() is the count-operator (number of parameters/datapoints)
//  smaller values are indicate a better model here because of the factor of -2.0
//====================================================================================================================
double SC_DistanceMeasures::BIC(SC_Model *pModel, double logLikelihood, unsigned long int dataCount) {
  return (-2.0 * logLikelihood) + (pModel->getFreeParameterCount() * log((double)dataCount));
}

//====================================================================================================================
//	The Bayesian Information Criterion as defined in 'Speaker, Environment And Channel Change Detection And Clustering 
//  Via The BIC', Chen, Gopalakrishnan, 1998. It assumes gaussian clusters and doesn't use the cluster's models, but 
//  works on the covars of original features only. If x and y are non-NULL, they get filled with scalars representing 
//  the two terms in BIC=x-alpha*y, needed for tuning the alpha value outside this class
//====================================================================================================================
double SC_DistanceMeasures::BIC(SC_Cluster *pFirst, SC_Cluster *pSecond, double alpha, double *x, double *y) {
	double **cov1 = pFirst->getCovar(), **cov2 = pSecond->getCovar(), **cov, det1, det2, det, a, b;
	unsigned long int M = pFirst->getFeatureVectorCount(), N = pSecond->getFeatureVectorCount(), dim = pFirst->getSpeechFrames()->Col;
	SV_Data *pComplete, *pHook = pFirst->getSpeechFrames(pFirst->getSegmentCount()-1), *pSave;

  //get a matrix-function-object, if it doesn't already exist
  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

	//get the 3 determinats; its a little pain to compute the combined one, the rest is just plain ;-)
	pSave = pHook->Next;
	pHook->Next = pSecond->getSpeechFrames();
	pComplete = pFirst->getSpeechFrames()->MergeData(pFirst->getSegmentCount()+pSecond->getSegmentCount());
	pHook->Next = pSave;
	cov = this->pMatrixFunc->cov(pComplete->Mat, pComplete->Row, pComplete->Col);
	det = this->pMatrixFunc->logDet(cov, dim); //logDet() works better when cov is nearly singular so that det=0 would be returned
	MFree_0D(pComplete);
	det1 = this->pMatrixFunc->logDet(cov1, dim);
	det2 = this->pMatrixFunc->logDet(cov2, dim);
	MFree_2D(cov1);
	MFree_2D(cov2);
	MFree_2D(cov);

	//definition from 'Speaker, Environment And Channel Change Detection And Clustering Via The BIC', Chen, Gopalakrishnan, 1998
	//a = (double)(M+N)*log(det) - (double)(M)*log(det1) - (double)(N)*log(det2);
	a = (double)(M+N)*det - (double)(M)*det1 - (double)(N)*det2;
	b =  0.5 * ((double)(dim) + 0.5*(double)(dim)*(double)(dim+1)) * log((double)(M+N));

	if (x != NULL) {
		*x = a;
	}
	if (y != NULL) {
		*y = b;
	}

	return a - alpha*b;
}

//====================================================================================================================
//	Within Cluster Dispersion as defined in 'Automatic Speaker Clustering', Jin, Kubala, Schwartz, 1997 or
//  'Online Speaker Clustering', Liu, Kubala, IEEE ICASSP 2003
//
//  WCD = |W| + C*log(k)
//  with k = number of clusters
//       W = SUM_{all clusters}(#segmentsPerCluster * covarOfCluster)
//       C = pTweak->distanceMeasures.WCDpenaltyFactor; should be coosen such that the 2 terms have equal order of 
//           magnitude
//  if x and y are non-NULL, the two parts of WCD=x+penalty*y are returned for penalyt factor tuning purposes
//====================================================================================================================
double SC_DistanceMeasures::withinClusterDispersion(SC_Partition *pPartition, double *x, double *y) {
  SC_Cluster *pHook = pPartition->getClusters();
  unsigned int clusterCount = 0, dim = pHook->getFeatureDim();
  double res, det, **cov = NULL, **sum = NULL, penalty;

  //get a matrix-function-object, if it doesn't already exist
  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

  //initialize the sum-matrix
  sum = this->pMatrixFunc->zeros(dim, dim);

  //calc the sum and cluster-count
  while (pHook != NULL) {
    cov = pHook->getCovar();
    this->pMatrixFunc->mult(cov, (double)(pHook->getFeatureVectorCount()), dim, dim, true);
    this->pMatrixFunc->add(sum, cov, dim, dim, true);

    MFree_2D(cov);
    clusterCount++;
    pHook = pHook->Next;
  }

  //get the determinant and compute the final result
	//jin, kubala & schwartz: res = det * sqrt(clusterCount); also suggested: det + C*sqrt(clusterCount)   or   det + C*log(clusterCount)
	//liu, kubala: res = det * sqrt(clusterCount);
  det = this->pMatrixFunc->det(sum, dim);
	penalty = this->pTweak->distanceMeasure.WCDpenaltyFactor * log((double)(clusterCount));
  res = det + penalty; //for the penalty to work well, it has to be the same order of magnitude than the det!
	//sclib::tupelOut("wcd.txt", det, penalty, this->pTweak);

	if (x != NULL) {
		*x = det;
	}
	if (y != NULL) {
		*y = log((double)(clusterCount));
	}

  MFree_2D(sum);

  return res;
};

//====================================================================================================================
//	The Information Change Rate as defined in 'A Robust Stopping Criterion for Agglomerative Hierarchical Clustering 
//  in a Speaker Diarization System', Han, Narayanan, 2007 (assuming gaussian clusters)
//
//  ICR = 1/(M+N) * ln(GLR)
//  with M, N = number of features in each of the 2 segments
//       GLR = the Generalized Likelihood Ratio (based on single multivariate gaussians) between the two segments
//====================================================================================================================
double SC_DistanceMeasures::ICR(SC_Cluster *pFirst, SC_Cluster *pSecond) {
	double res, **cov1, **cov2, **cov, det1, det2, det;
	int T, M, N, dim;

  //get a matrix-function-object, if it doesn't already exist
  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

	//get the 3 determinats
	cov = pFirst->getCombinedCovar(pSecond, T, dim, false, false); //we get problems if the time model replaced the training data (because of high-dimesional features AND normalization!) but removed code to unwind this because it needs so much prior knowledge of what was done in prior prcessing steps that can't be packed into parameters - solution: just don't replace training data!
	det = this->pMatrixFunc->logDet(cov, dim); //logDet() works better when cov is nearly singular so that det=0 would be returned
	MFree_2D(cov);
	cov1 = pFirst->getCovar(M, dim, false, false);
	det1 = this->pMatrixFunc->logDet(cov1, dim); //logDet() works better when cov is nearly singular so that det=0 would be returned
	MFree_2D(cov1);
	cov2 = pSecond->getCovar(N, dim, false, false);
	det2 = this->pMatrixFunc->logDet(cov2, dim); //logDet() works better when cov is nearly singular so that det=0 would be returned
	MFree_2D(cov2);
	//sclib::matrixOutEx("cov1.txt", cov1, dim, dim, this->pTweak, ",");
	//sclib::matrixOutEx("cov2.txt", cov2, dim, dim, this->pTweak, ",");
	//sclib::matrixOutEx("cov.txt", cov, dim, dim, this->pTweak, ",");

	//GLR is logP(x|l_X)+logP(y|l_Y)-logP(xy|l_XY); logP(x|l_X)=-d/2*log(2pi) - 1/2*log|Sigma_x| - 1/2*d (if the gaussian parameters are estimated using only x); the rest below is algebra
	res = -1.0*(double)(dim)*sclib::log_2pi - (double)(dim); 
	//res -= log(det1) - log(det2) + log(det);
	res -= det1 - det2 + det; //TODO: signs correct?
	res /= 2.0;

	res *= 1.0 / ((double)(M) + (double)(N));

	return res;
}

//====================================================================================================================
//	Function to calcualte the divergence-shape distance between to segments of speech by using only their 
//	covariance-matrices
//	According to: 'speaker recognition: a tutorial', j.campbell, 1997
//	              'real-time unsupervised speaker change detection', l.lu, h.j.zhang, 2002
//
//	the formula in short (the mean-parts are discarded, because they are easily biased by noise etc.):
//	D(1,2) = trace( (cov_1 - cov_2) * (inv(cov_2) - inv(cov_1)) )
//
//	own additions:
//		-	the abs-value of the distance is returned, so that the result is positive in any case
//    - a return-value of -1 indicates an error (matrix couldn't be inverted)
//    - a minimum variance can be specified to avoid variances to be too close to 0 in any direction (+/-)
//====================================================================================================================
double SC_DistanceMeasures::divergenceShape(double** cov_1, double** cov_2, unsigned long int dim, double minVariance) {
	double dist = 0.0, tmp1, tmp2;
	double **diff_1, **diff_2, **inv_1, **inv_2;
	int res;
	unsigned long int i, j;
  
  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

	if (dim == 1) {
		tmp1 = (cov_1[0][0] == 0.0) ? 1.0/minVariance : 1.0/cov_1[0][0];
		tmp2 = (cov_2[0][0] == 0.0) ? 1.0/minVariance : 1.0/cov_2[0][0];
		dist = (cov_1[0][0] - cov_2[0][0]) * (tmp2 - tmp1);
	} else {
		MArray_2D(diff_1, (long int)(dim), (long int)(dim), double, "SC_DistanceMeasures.divergenceShape: diff_1");
		MArray_2D(diff_2, (long int)(dim), (long int)(dim), double, "SC_DistanceMeasures.divergenceShape: diff_2");
		MArray_2D(inv_1, (long int)(dim), (long int)(dim), double, "SC_DistanceMeasures.divergenceShape: inv_1");
		MArray_2D(inv_2, (long int)(dim), (long int)(dim), double, "SC_DistanceMeasures.divergenceShape: inv_2");

		for (i = 0; i < dim; i++) { //subtract cov2 form cov1, create copies of them as input to the inv()-methods, regard minVariance
			for (j = 0; j < dim; j++) {
				tmp1 = sclib::sg(cov_1[i][j]) * sclib::max(minVariance, fabs(cov_1[i][j]));
				tmp2 = sclib::sg(cov_2[i][j]) * sclib::max(minVariance, fabs(cov_2[i][j]));
				diff_1[i][j] = tmp1 - tmp2;
				inv_1[i][j] = tmp1;
				inv_2[i][j] = tmp2;
			}
		}

		res	= this->pMatrixFunc->inv(inv_1, dim); //build first inverse, check for errors
		if (res < 0) {
			dist = -1.0;
		} else {
			res = this->pMatrixFunc->inv(inv_2, dim); //build second inverse, check for errors
			if (res < 0) {
				dist = -1.0;
			} else {
				
				dist = 0.0;
				for (i = 0; i < dim; i++) { //for all elements on the diagonal
					for (j = 0; j < dim; j++) { //for all elements within one row/column
						dist += diff_1[i][j] * (inv_2[j][i] - inv_1[j][i]); //3in1: difference of the 2 inverses, multiplication of the 2 differences (only for the necessary cells on the main diagonal, respectively) and trace of the result
					}
				}

				if (dist < 0.0) {
					//sclib::matrixOut("cov_1.txt", cov_1, dim, dim, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
					//sclib::matrixOut("cov_2.txt", cov_2, dim, dim, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
					//sclib::matrixOut("diff_1.txt", diff_1, dim, dim, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
					//sclib::matrixOut("inv_1.txt", inv_1, dim, dim, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
					//sclib::matrixOut("inv_2.txt", inv_2, dim, dim, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
					REPORT_ERROR(SVLIB_Fail, "negative divergence shape distance occured");
					dist = fabs(dist); //TODO: email to j.p.campbell written on 29.05.2007 regarding this issue... answer: shouldn't become negative, at least he can't remember he had that problem
				}

			}
		}

		MFree_2D(inv_1);
		MFree_2D(inv_2);
		MFree_2D(diff_1);
		MFree_2D(diff_2);
	}

	return dist;
}

//====================================================================================================================
//	The Earth Mover's Distance for two Speaker Models
//====================================================================================================================
double SC_DistanceMeasures::EMD(SC_Model* pModel_1, SC_Model* pModel_2) {
  SC_Signature *model1 = pModel_1->toSignature(), *model2 = pModel_2->toSignature();
  double dist = EMD(model1, model2, this->pTweak);

	pModel_1->killSignature(model1);
  pModel_2->killSignature(model2);
 
  return dist;
}

//====================================================================================================================
//	The Earth Mover's Distance directly for two signatures
//  This is a wrapper around the SC_EMD class; the EMD-implementation is the one from it's author, Yossi Rubner
//  See also: Y.Rubner, C.Tomasi, L.J.Guibas, "The Earth Mover's Distance as a Metric for Image Retrieval", 
//            International Journal of Computer Vision 40(2) 99-121, Kluwer Academic Publishing, 2000
//====================================================================================================================
double SC_DistanceMeasures::EMD(SC_Signature* pSignature1, SC_Signature* pSignature2, SC_TweakableParameters *pTweak) {
	SC_EMD earthMover(pTweak);
  return earthMover.emd(pSignature1, pSignature2, NULL, NULL);
}

//====================================================================================================================
//	Beigi distance measure between 2 gaussian models with diagonal covariances
//  The actual distance measure used internally is chosen by the centroids in the signatures towhich the models are 
//  converted and is controlled by the tweakable parameters
//====================================================================================================================
double SC_DistanceMeasures::beigi(SC_Model* pModel_1, SC_Model* pModel_2) {
  SC_Signature *model1 = pModel_1->toSignature(), *model2 = pModel_2->toSignature();
	double res = beigi(model1, model2);

  pModel_1->killSignature(model1);
  pModel_2->killSignature(model2);
  
  return res;
}

//====================================================================================================================
//	Beigi distance measure directly between to signatures
//  According to: H.S.M.Beigi, S.H.Maes, J.S.Sorensen, "A Distance Measure Between Collections of Distributions and 
//                it's Application to Speaker Recognition"
//====================================================================================================================
double SC_DistanceMeasures::beigi(SC_Signature* pSignature1, SC_Signature* pSignature2) {
  unsigned short int i, j, N = pSignature2->getN(), M = pSignature1->getN();
  double dist, *Wn, *Wm, wm = 0.0, wn = 0.0, sumWeight1 = 0.0, sumWeight2 = 0.0, res;

  MArray_1D(Wm, M, double, "SC_DistanceMeasures.beigi: Wm");
  MArray_1D(Wn, N, double, "SC_DistanceMeasures.beigi: Wn");
  for (j = 0; j < N; j++) {Wn[j] = numeric_limits<double>::max();}

  for (i = 0; i < M; i++) {
		Wm[i] = numeric_limits<double>::max();

    for (j = 0; j < N; j++) {
			dist = (double)(pSignature1->getCentroid(i)->getDistance(pSignature2->getCentroid(j)));
      if (dist < Wm[i]) {Wm[i] = dist;}
      if (dist < Wn[j]) {Wn[j] = dist;}
    }

		wm += Wm[i] * pSignature1->getWeight(i);
    sumWeight1 += pSignature1->getWeight(i);
  }

  for (j = 0; j < N; j++) {
    wn += Wn[j] * pSignature2->getWeight(j);
    sumWeight2 += pSignature2->getWeight(j);
  }
  MFree_1D(Wm);
  MFree_1D(Wn);

  //the denominator is the sum of all the counts of feature-vectors which build up each mixture;
  //normally our "counts" are weights, whichs sum up to 1, this yields 2 (all weights of model 1(=1) + all weight of model 2(=1) = 2
  //but who knows what happens in the future, so we really compute the sum...
  //the maskLevel happened in the future... ;-)
  res = (wm + wn) / (sumWeight1 + sumWeight2); 
  
  return res;
}

//====================================================================================================================
//	Beigi distance measure between a signature and a vector
//  According to: H.S.M.Beigi, S.H.Maes, J.S.Sorensen, "A Distance Measure Between Collections of Distributions and 
//                it's Application to Speaker Recognition"
//====================================================================================================================
double SC_DistanceMeasures::beigi(SC_Signature* signature, float *vector, int dim, double *Ws, double *Wv) {
  unsigned short int i, j, N = signature->getN(), M = dim;
  double dist, wm = 0.0, wn = 0.0, sumWeight1 = 0.0, sumWeight2 = 0.0, res;
	int squareSize = (int)(sqrt((double)dim));
	double coord[2];
	double *Wm = Wv, *Wn = Ws;
	SC_Centroid_Point* pPoint = NULL;

	if (signature->getCentroid(0)->getCentroidType() != sclib::centroidPoint) {
		return -1.0;
	}

  for (j = 0; j < N; j++) {
		Wn[j] = numeric_limits<double>::max();
	}

  for (i = 0; i < M; i++) { //components of the vector
		Wm[i] = numeric_limits<double>::max();
		if (vector[i] > 0.0) {
			coord[0] = (double)(i / squareSize);
			coord[1] = (double)(i % squareSize);

			for (j = 0; j < N; j++) { //centroids in the signature
				pPoint = (SC_Centroid_Point*)(signature->getCentroid(j));
				dist = euclid(pPoint->getCoordinate(), coord, 2);
				if (dist < Wm[i]) {Wm[i] = dist;}
				if (dist < Wn[j]) {Wn[j] = dist;}
			}

			wm += Wm[i] * vector[i];
			sumWeight1 += vector[i];
		}
  }

  for (j = 0; j < N; j++) {
    wn += Wn[j] * signature->getWeight(j);
    sumWeight2 += signature->getWeight(j);
  }

  //the denominator is the sum of all the counts of feature-vectors which build up each mixture;
  //normally our "counts" are weights, whichs sum up to 1, this yields 2 (all weights of model 1(=1) + all weight of model 2(=1) = 2
  //but who knows what happens in the future, so we really compute the sum...
  //the maskLevel happened in the future... ;-)
  res = (wm + wn) / (sumWeight1 + sumWeight2); 
  
  return res;
}

//====================================================================================================================
//	The std. Euclidean Distance measure between two vectors in libSVM format
//====================================================================================================================
double SC_DistanceMeasures::euclid(SC_SVMnode *x, SC_SVMnode *y) {
	double res = 0.0;
	int	d = 0, dd = 0;

	while (x[d].index!=-1 || y[dd].index!=-1) {
		if (x[d].index == y[dd].index) {
			res += pow(x[d++].value-y[dd++].value, 2.0);
		} else if (x[d].index == -1) {
			res += pow(y[dd++].value, 2.0);
		} else if (y[dd].index == -1) {
			res += pow(x[d++].value, 2.0);
		} else if (x[d].index > y[dd].index) {
			res += pow(y[dd++].value, 2.0);
		} else if (y[dd].index > x[d].index) {
			res += pow(x[d++].value, 2.0);
		}
	}
	res = sqrt(res);

	return res;
}

//====================================================================================================================
//	2 methods to deal with a 'flat' matrix as created by the above getPairWiseDistance()-method
//====================================================================================================================
void SC_DistanceMeasures::set1DMatrixValue(double *vectorMatrix, unsigned long int length, unsigned long int row, unsigned long int column, double value) {
  if (row > column) {
    set1DMatrixValue(vectorMatrix, length, column, row, value);
  } else if (row < column) {
    vectorMatrix[row*length + column - (row+1)*(row+2)/2] = value;
  } //ignore row == column

  return;
}

double SC_DistanceMeasures::get1DMatrixValue(double *vectorMatrix, unsigned long int length, unsigned long int row, unsigned long int column) {
  if (row == column) {
    return 0.0;
  } else if (row > column) {
    return get1DMatrixValue(vectorMatrix, length, column, row);
  } else { // for explanation of index formula, see set method
    return vectorMatrix[row*length + column - (row+1)*(row+2)/2]; 
  }
}

//====================================================================================================================
//	Distance between the support of two one-class SVMs (given in libSVMs SVMmodel structure) as proposed in the paper
//  "An Online Kernel Change Detection Algorithm", Desobry, Davy, Doncarli, 2005
//  The pSVM parameter is an SVM object containing the algorithm used for computation, not a trained SVM itself!
//
//  The formulas can be looked up entirely in the paper ((10), (14), (15) and (16) are relevant), but the idea is as 
//  follows: A one-class SVM doesn't solve the pdf-estimation problem, but the simpler one (hence less data necessary)
//  of finding the level set (probable area of the pdf) of the data; this region is represented as a section of the 
//  unit circle in svm-space to which the original data is mapped by the kernel function; this distance now measures
//  the amount of overlap between two such sections represented in two svm-models by taking the arc distance between 
//  the centers of the segments normalized by their width (spread, variance). The below code is how this can be 
//  calculated via the methods and parameters of the libSVM.
//====================================================================================================================
double SC_DistanceMeasures::svmArcDistance(SC_SVM *pSVM, SC_SVMmodel *oneClassSVM_1, SC_SVMmodel *oneClassSVM_2) {
	int i, j;
	double val, crossTerm = 0.0, sqrtTerm_1 = 0.0, sqrtTerm_2 = 0.0, res = 0.0, numerator, denominator;
	double eps1 = 0.000001, eps2 = 0.000001; //do avoid division by zero as suggested in the paper (remark 1)

	if (oneClassSVM_1->param.gamma != oneClassSVM_2->param.gamma) {
		REPORT_ERROR(SVLIB_BadData, "Kernel parameters for arc distance computation have to be equal");
	}

	//TODO: check parameters in libsvm with variables in formulas - why sv_coef[0]? why the SVs instead of all training data?
	if (oneClassSVM_1->param.svm_type == SCLIB_SVM_TYPE_ONECLASS && oneClassSVM_2->param.svm_type == SCLIB_SVM_TYPE_ONECLASS) {
		for (i = 0; i < oneClassSVM_1->l; i++) {
			val = 0.0;
			for (j = 0; j < oneClassSVM_2->l; j++) {
				val += oneClassSVM_2->sv_coef[0][j] * pSVM->evaluateKernel(oneClassSVM_1->SV[i], oneClassSVM_2->SV[j], oneClassSVM_2->param);
			}
			crossTerm += oneClassSVM_1->sv_coef[0][i] * val;
		}

		for (i = 0; i < oneClassSVM_1->l; i++) {
			val = 0.0;
			for (j = 0; j < oneClassSVM_1->l; j++) {
				val += oneClassSVM_1->sv_coef[0][j] * pSVM->evaluateKernel(oneClassSVM_1->SV[i], oneClassSVM_1->SV[j], oneClassSVM_2->param);
			}
			sqrtTerm_1 += oneClassSVM_1->sv_coef[0][i] * val;
		}

		for (i = 0; i < oneClassSVM_2->l; i++) {
			val = 0.0;
			for (j = 0; j < oneClassSVM_2->l; j++) {
				val += oneClassSVM_2->sv_coef[0][j] * pSVM->evaluateKernel(oneClassSVM_2->SV[i], oneClassSVM_2->SV[j], oneClassSVM_2->param);
			}
			sqrtTerm_2 += oneClassSVM_2->sv_coef[0][i] * val;
		}

		numerator = sclib::min(crossTerm / sqrt(sqrtTerm_1 * sqrtTerm_2), 1.0); //min() to avoid error in case of rounding errors
		sqrtTerm_1 = sclib::min(oneClassSVM_1->rho[0] / sqrt(sqrtTerm_1), 1.0);
		sqrtTerm_2 = sclib::min(oneClassSVM_2->rho[0] / sqrt(sqrtTerm_2), 1.0);

		if (numerator > 1.0 || sqrtTerm_1 > 1.0 || sqrtTerm_2 > 1.0) {
			REPORT_ERROR(SVLIB_BadData, "Arccos() cannot take arguments outside the range [-1, 1]");
			res = -1.0;
		} else {
			numerator = acos(numerator); //arc-distance between segment-centers (arc-distance means: dist. between two points not measured along a straight line but a piece of an arc)
			denominator = acos(sqrtTerm_1) * acos(sqrtTerm_2); //spread of both regions to scale the distance accordingly

			if (numerator == 0.0 && denominator == 0.0) {
				res = 0.0; //otherwise distance would be 1.0 which is against intuition if both segments equal 
			} else {
				res = (numerator + eps1) / (denominator + eps2);
			}
		}
	}
	return res;
}

//====================================================================================================================
//	"New" weighted euclidean distance measure between two (consecutive) segments of features to detect a changepoint
//  between both; the weights are higher for those attributes (feature-dims) that have a higher between-segment 
//  variance as compared to the sum of the two inter-segment variances (and are even teared further apart by a 
//  sigmoidal function). The idea is due to (and the formulas can be looked up in):
//  "Speaker Change Detection Using a New Weighted Distance Measure", Kwon, Narayanan, ICASSP 2002
//
//  The implementation here relies on the pure data of the segments and not on their means (they are used for the pure
//  distance) so that the weights can be computed inside, too
//====================================================================================================================
double SC_DistanceMeasures::newWeightedEuclideanDistance(SV_Data *pData_1, SV_Data *pData_2) {
	int d, dim = pData_1->Col;
	double p = 200.0, res = 0.0;
	double *totalVar, *seg1var, *seg2var, *seg1mean, *seg2mean, *weight;
	double weightMean = 1.0 / (double)(dim), norm = 0.0;
	SV_Data *pCompleteData, *pHook = pData_1->Next;
	double eps1 = 0.000001, eps2 = 0.000001; //to avoid division by zero

  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

	if (pData_1->Col == pData_2->Col) {
		MArray_1D(weight, dim, double, "SC_DistanceMeasures.newWeightedEuclideanDistance: weight");

		pData_1->Next = pData_2;
		pCompleteData = pData_1->MergeData(2);
		pData_1->Next = pHook;

		seg1mean = this->pMatrixFunc->mean(pData_1->Mat, pData_1->Row, dim);
		seg2mean = this->pMatrixFunc->mean(pData_2->Mat, pData_2->Row, dim);
		seg1var = this->pMatrixFunc->variance(pData_1->Mat, pData_1->Row, dim, seg1mean);
		seg2var = this->pMatrixFunc->variance(pData_2->Mat, pData_2->Row, dim, seg1mean);
		totalVar = this->pMatrixFunc->variance(pCompleteData->Mat, pCompleteData->Row, dim);

		for (d = 0; d < dim; d++) {
			weight[d] = (totalVar[d] + eps1) / (seg1var[d] + seg2var[d] + eps2);
			norm += weight[d];
		}
		for (d = 0; d < dim; d++) {
			weight[d] = sclib::sigmoid(weight[d]/norm - weightMean, p); //sigmoidal scaling, p was given in the paper with 150 for news and 200 for movie data... TODO: ok?
			res += weight[d] * (seg1mean[d] - seg2mean[d]) * (seg1mean[d] - seg2mean[d]); //weighted squared euclidean distance of the means of the two segments
		}

		MFree_1D(weight);
		MFree_1D(totalVar);
		MFree_1D(seg1var);
		MFree_1D(seg2var);
		MFree_1D(seg1mean);
		MFree_1D(seg2mean);
		MFree_0D(pCompleteData);
	} else {
		res = -1.0;
	}
	
	return res;
}

//====================================================================================================================
//	The T^2 statistic is a 1st order statistical test if the mean of the left side (0..b-1) equals the mean of the 
//  right side (b..T-1), gaussianity assumed. If b is >0, the value of T^2 for each possible b between the two 
//  segments is returned, otherwise only the specified value, each time in a newly allocated array. See also Zhou, 
//  Hansen, "Efficient Audio Stream Segmentation via the Combined T^2 Statistic and Bayesian Information Criterion", 
//  2005
//====================================================================================================================
double* SC_DistanceMeasures::tSquareStatistic(SV_Data *pSegment_1, SV_Data *pSegment_2, int b) {
	double *res = NULL, **data = NULL;
	int i, j, len = pSegment_1->Row+pSegment_2->Row;
	double *m1 = NULL, *m2 = NULL, **cov = NULL;

	assert(pSegment_1->Col == pSegment_2->Col);
	assert(b < pSegment_1->Row+pSegment_2->Row);

	//copy the data of both segments into a single container
	MArray_2D(data, pSegment_1->Row+pSegment_2->Row, pSegment_1->Col, double, "SC_DistanceMeasures.tSquareStatistic: data");
	for (j = 0; j < pSegment_1->Col; j++) {
		for (i = 0; i < pSegment_1->Row; i++) {
			data[i][j] = pSegment_1->Mat[i][j];
		}
		for (i = 0; i < pSegment_2->Row; i++) {
			data[pSegment_1->Row+i][j] = pSegment_2->Mat[i][j];
		}
	}

	if (b >= 0) {
		MArray_1D(res, 1, double, "SC_DistanceMeasures.tSquareStatistic: res");
		res[0] = tSquareStatistic(data, len, pSegment_1->Col, b, m1, m2, cov);
	} else {
		MArray_1D(res, len-1, double, "SC_DistanceMeasures.tSquareStatistic: res");
		for (i = 1; i < len; i++) {
			res[i-1] = tSquareStatistic(data, len, pSegment_1->Col, i, m1, m2, cov);
		}
	}

	MFree_1D(m1);
	MFree_1D(m2);
	MFree_2D(cov);
	MFree_2D(data);

	return res;
}

//====================================================================================================================
//	The T^2 statistic is a 1st order statistical test if the mean of the left side (0..b-1) equals the mean of the 
//  right side (b..T-1), gaussianity assumed. If b is <0, the value of T^2 for each possible b between the two 
//  segments is returned, otherwise only the specified value, each time in a newly allocated array. By start and end 
//  the segment can be specied that is tested; by setting middle1 and middle2 >0 and start<middle1<middle2<end a 
//  segment can be constructed without the middle part, i.e. start->middle1, middle2->end. See also Zhou, Hansen, 
//  "Efficient Audio Stream Segmentation via the Combined T^2 Statistic and Bayesian Information Criterion", 2005
//====================================================================================================================
double* SC_DistanceMeasures::tSquareStatistic(SV_Data *pSegment, int b, int start, int middle1, int middle2, int end) {
	double *res = NULL, **data = NULL;
	int i, j, len;
	int s, m1, m2, e;
	double *ml = NULL, *mr = NULL, **cov = NULL;

	s = (start >= 0) ? start : 0;
	e = (end >= 0 && end < pSegment->Row) ? end : pSegment->Row-1;
	m1 = (middle1 > s) ? middle1 : (e-s)/2;
	m2 = (middle2 > m1 && middle2 < e) ? middle2 : ((e-s)/2)+1;
	len = m1-s+1 + e-m2+1;

	assert(b<0 || sclib::isBetween(s, b, m1) || sclib::isBetween(m2, b, e));

	//copy the data of both segments into a single container
	MArray_2D(data, len, pSegment->Col, double, "SC_DistanceMeasures.tSquareStatistic: data");
	for (j = 0; j < pSegment->Col; j++) {
		for (i = s; i <= m1; i++) {
			data[i-s][j] = pSegment->Mat[i][j];
		}
		for (i = m2; i <= e; i++) {
			data[i-m2+m1+1][j] = pSegment->Mat[i][j];
		}
	}

	if (b >= 0) {
		MArray_1D(res, 1, double, "SC_DistanceMeasures.tSquareStatistic: res");
		res[0] = tSquareStatistic(data, len, pSegment->Col, b, ml, mr, cov);
	} else {
		MArray_1D(res, len-1, double, "SC_DistanceMeasures.tSquareStatistic: res");
		for (i = 1; i < len; i++) {
			res[i-1] = tSquareStatistic(data, len, pSegment->Col, i, ml, mr, cov);
		}
	}

	MFree_1D(ml);
	MFree_1D(mr);
	MFree_2D(cov);
	MFree_2D(data);

	return res;
}

//====================================================================================================================
//	Gives the cross likelihood ratio between two mutlivariate, diagonal-covariance gaussians fitted to both datasets
//  of equal size according to the paper "A Simple but Effective Approach to Speaker Tracking in Broadcast News", 
//  Rodríguez, Penagarikano, Bordel, 2007
//====================================================================================================================
double SC_DistanceMeasures::crossLikelihoodRatio(SV_Data *pData_1, SV_Data *pData_2) {
	int d, t, T = pData_1->Row, D = pData_1->Col;
	double res = 0.0;
	double *seg1var, *seg2var, *seg1mean, *seg2mean, *seg1scaled, *seg2scaled;
	double ll_11 = 0.0, ll_12 = 0.0, ll_22 = 0.0, ll_21 = 0.0;
	const double minVar = 0.001;

  if (this->pMatrixFunc == NULL) {
    this->pMatrixFunc = new SC_MatrixFunctions();
    this->linkedMatrixFunc = false;
  }

	if (pData_1->Col==pData_2->Col && pData_1->Row==pData_2->Row) {
		seg1mean = this->pMatrixFunc->mean(pData_1->Mat, pData_1->Row, D);
		seg2mean = this->pMatrixFunc->mean(pData_2->Mat, pData_2->Row, D);
		seg1var = this->pMatrixFunc->variance(pData_1->Mat, pData_1->Row, D, seg1mean);
		seg2var = this->pMatrixFunc->variance(pData_2->Mat, pData_2->Row, D, seg1mean);

		MArray_1D(seg1scaled, D, double, "SC_DistanceMeasures.crossLikelihoodRatio: seg1scaled");
		MArray_1D(seg2scaled, D, double, "SC_DistanceMeasures.crossLikelihoodRatio: seg2scaled");
		for (d = 0; d < D; d++) { //precompute first term in gaussian function
			if (seg1var[d] < minVar) { //apply a variance floor to avoid numerical difficulties
				seg1var[d] = minVar;
			}
			if (seg2var[d] < minVar) { //apply a variance floor to avoid numerical difficulties
				seg2var[d] = minVar;
			}
			seg1scaled[d] = log(sclib::one_div_sqrt_2pi * sqrt(seg1var[d]));
			seg2scaled[d] = log(sclib::one_div_sqrt_2pi * sqrt(seg2var[d]));
		}

		for (t = 0; t < T; t++) {
			for (d = 0; d < D; d++) { //compute univariate gaussians, because of diag. covar it can be added together to multivariate version
				ll_11 += seg1scaled[d] + (((double)(pData_1->Mat[t][d]) - seg1mean[d]) * ((double)(pData_1->Mat[t][d]) - seg1mean[d]) / (-2.0 * seg1var[d]));
				ll_12 += seg2scaled[d] + (((double)(pData_1->Mat[t][d]) - seg2mean[d]) * ((double)(pData_1->Mat[t][d]) - seg2mean[d]) / (-2.0 * seg2var[d]));
				ll_22 += seg2scaled[d] + (((double)(pData_2->Mat[t][d]) - seg2mean[d]) * ((double)(pData_2->Mat[t][d]) - seg2mean[d]) / (-2.0 * seg2var[d]));
				ll_21 += seg1scaled[d] + (((double)(pData_2->Mat[t][d]) - seg1mean[d]) * ((double)(pData_2->Mat[t][d]) - seg1mean[d]) / (-2.0 * seg1var[d]));
			}
		}
		res = -1.0 * ((ll_12+ll_21) - (ll_11+ll_22)); //ll_* are the loglikelihoods, first number gives the segmentnumber of testdata, second number is the number of the gaussian parameters

		MFree_1D(seg1var);
		MFree_1D(seg2var);
		MFree_1D(seg1mean);
		MFree_1D(seg2mean);
		MFree_1D(seg1scaled);
		MFree_1D(seg2scaled);
	} else {
		res = -1.0;
	}
	
	return res;
}

//====================================================================================================================
//	The information measure as defined in "Data Mining", Witten & Frank, 2005, pp100
//====================================================================================================================
double SC_DistanceMeasures::information(double p1, double p2) {
	double sumP = p1 + p2;
	double p1Ratio, p2Ratio, info = 0.0; //information in bits

	if (sumP > 0.0) {
		p1Ratio = p1/sumP;
		p2Ratio = p2/sumP;
		info -= (p1Ratio > 0.0) ? p1Ratio*sclib::log2(p1Ratio) : 0.0;
		info -= (p2Ratio > 0.0) ? p2Ratio*sclib::log2(p2Ratio) : 0.0;
	}

	return info;
}

//====================================================================================================================
//	The information measure as defined in "Data Mining", Witten & Frank, 2005, pp100
//  just as above, but getting an dim-dimensional array of arguments
//====================================================================================================================
double SC_DistanceMeasures::information(double *p_i, int dim) {
	double sumP = 0.0, info;
	int i;
	
	for (i = 0; i < dim; i++) {
		sumP += p_i[i];
	}

	info = sumP * sclib::log2(sumP);
	for (i = 0; i < dim; i++) {
		info -= p_i[i] * sclib::log2(p_i[i]);
	}
	info /= sumP;

	return info; //information in bits
}

//====================================================================================================================
//	The information measure as defined in "Data Mining", Witten & Frank, 2005, pp100
//  implementation of witten&frank's notation info([p11, p12], [p21, p22])
//====================================================================================================================
double SC_DistanceMeasures::information(double p11, double p12, double p21, double p22) {
	double sumP1 = p11+p12, sumP2 = p21+p22, sumP = sumP1+sumP2;
	double info = 0.0;

	if (sumP > 0.0) {
		info = (sumP1/sumP)*information(p11, p12) + (sumP2/sumP)*information(p21, p22); //information in bits
	}

	return info;
}

//====================================================================================================================
//	The MDL criterion as defined in "Data Mining", Witten & Frank, 2005, pp301, for the task of evaluating a split in 
//  attribute discretization: gain is the information gain gained by the split, N is the number of instances, k the 
//  number of classes, k1 and k2 the number of classes appearing in the first and second half created by the split, E 
//  the entropy of the instances, E1 and E2 the entropy of the instances in the two subintervals; the boolean return 
//  value indicates if this split should be done
//====================================================================================================================
bool SC_DistanceMeasures::splittingMDL(double gain, int N, int k, int k1, int k2, double E, double E1, double E2) {
	double MDL = sclib::log2((double)(N)-1.0)/(double)(N) + (sclib::log2(pow(3.0, (double)(k))) - (double)(k*E + k1*E1 + k2*E2))/(double)(N);

	return (gain > MDL) ? true : false;
}
