/**************************************************************************/
/*	A class to hold the implementations of all used distance-measures			*/
/*	throughout the SC_* project																						*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.04.2005																								*/
/**************************************************************************/

#ifndef __SC_DistanceMeasures_H__
#define __SC_DistanceMeasures_H__

#include "SC_TweakableParameters.h"
#include "SC_ModelHandler.h"
#include "SC_MatrixFunctions.h"
#include "SC_EMD.h"
#include "SC_Cluster.h"
#include "SC_Partition.h"
#include "SC_SVM.h"
#include <SV_Data.h>
#include <SV_Error.h>

class SC_DistanceMeasures {

	private :

	protected:

		SC_TweakableParameters*	pTweak;	//for tweakable parameters
    SC_MatrixFunctions* pMatrixFunc; //for some matrix calculations in BIC/divergenceShape
    SC_ModelHandler* pModelHandler; //for some model-operations
    
    bool verbose; //to control some console output...
    bool linkedMatrixFunc; //shows whether pMatrixFunc was linked in the constructor or created in this class (for destruction)

	public :

		//====================================================================================================================
		//	constructor/destructor
		//====================================================================================================================
		SC_DistanceMeasures(SC_TweakableParameters* pTweak = NULL, SC_MatrixFunctions* pMatrixFunc = NULL, bool verbose = true);
		virtual ~SC_DistanceMeasures();

		//====================================================================================================================
		//	Cross-Likelihood Ratio distance measure between this and the second cluster
    //	According to: 'Unsupervised Speaker Indexing Using Speaker Model Selektion Based On BIC', Nishida, Kawahara, 
    //								 2003 (IEEE ICASSP)
    //								'Blind Clustering of speech utterances based on speaker and language characteristics', Reynolds,
    //								 Singer, Carlson, McLaughlin, O'Leary, Zissman, 1998 (IEEE ICSLP)
    //  Needs (as parameters) for both models to test against each other:
    //    - the model
    //    - (linked list of) background-models
    //    - (linked list of) training feature vectors
    //    - count of list-elements to merge 
		//====================================================================================================================
		double CLR(SC_Model* pModel_1, SC_Model* pBkModels_1, SV_Data* pData_1, unsigned long int segmentCount_1, SC_Model* pModel_2, SC_Model* pBkModels_2, SV_Data* pData_2, unsigned long int segmentCount_2);

		//====================================================================================================================
		//	Generalized Likelihood Ratio distance measure between this and the second cluster
    //	According to: 'Online Speaker Clustering', Liu, Kubala, 2004 (IEEE ICASSP)
    //								'Segregation Of Speakers For Speech Recognition And Speaker Identification', Gish, Schmidt, 1991
    //								 (IEEE ICASSP)
    //  Needs (as parameters) for both models to test against each other:
    //    - the model
    //    - (linked list of) background-models
    //    - (linked list of) training feature vectors
    //    - count of list-elements to merge 
		//====================================================================================================================
		double GLR(SC_Model* pModel_1, SC_Model* pBkModels_1, SV_Data* pData_1, unsigned long int segmentCount_1, SC_Model* pModel_2, SC_Model* pBkModels_2, SV_Data* pData_2, unsigned long int segmentCount_2);

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
		double BIC(SC_Model* pCompleteModel, SC_Model* pCompleteBkModels, SV_Data* pCompleteData, unsigned long int completeSegmentCount, SC_Model* pModel_1, SC_Model* pBkModels_1, SV_Data* pData_1, unsigned long int segmentCount_1, SC_Model* pModel_2, SC_Model* pBkModels_2, SV_Data* pData_2, unsigned long int segmentCount_2, double alpha = 0.0);

    //====================================================================================================================
    //	The standard BIC for one concrete Model and it's test-result:
    //  
    //  BIC(model, data) = -2*log(likelihood(data|model)) + #(model)*sclib::ln(#(data))
    //
    //  where #() is the count-operator (number of parameters/datapoints)
    //  smaller values are indicate a better model here because of the factor of -2.0
    //====================================================================================================================
    double BIC(SC_Model *pModel, double logLikelihood, unsigned long int dataCount);

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
    double BIC(double** completeCov, double** cov_1, unsigned long int trainingDataCount_1, double** cov_2, unsigned long int trainingDataCount_2, unsigned short int dim, double lambda);
		double BIC(SV_Data *pFirst, SV_Data *pSecond);

    //====================================================================================================================
		//	The Bayesian Information Criterion as defined in 'Speaker, Environment And Channel Change Detection And Clustering 
		//  Via The BIC', Chen, Gopalakrishnan, 1998. It assumes gaussian clusters and doesn't use the cluster's models, but 
		//  works on the covars of original features only. If x and y are non-NULL, they get filled with scalars representing 
		//  the two terms in BIC=x-alpha*y, needed for tuning the alpha value outside this class
    //====================================================================================================================
    double BIC(SC_Cluster *pFirst, SC_Cluster *pSecond, double alpha, double *x = NULL, double *y = NULL);

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
    double withinClusterDispersion(SC_Partition *pPartition, double *x = NULL, double *y = NULL);

    //====================================================================================================================
    //	The Information Change Rate as defined in 'A Robust Stopping Criterion for Agglomerative Hierarchical Clustering 
		//  in a Speaker Diarization System', Han, Narayanan, 2007 (assuming gaussian clusters)
    //
    //  ICR = 1/(M+N) * ln(GLR)
    //  with M, N = number of features in each of the 2 segments
    //       GLR = the Generalized Likelihood Ratio (based on single multivariate gaussians) between the two segments
    //====================================================================================================================
    double ICR(SC_Cluster *pFirst, SC_Cluster *pSecond);

    //====================================================================================================================
    //	Function to calcualte the divergence-shape distance between to segments of speech by using only their 
    //	covariance-matrices
    //	According to: 'speaker recognition: a tutorial', j.campbell, 1997
    //	              'real-time unsupervised speaker change detection', l.lu, h.j.zhang, 2002
    //
    //	the formula in short (the mean-parts are discarded, because they are easily biased by noise etc.):
    //	D(i,j) = trace( (cov_1 - cov_2) * (inv(cov_1) - inv(cov_2)) )
    //
    //	own additions:
    //		-	the abs-value of the distance is returned, so that the result is positive
    //    - a return-value of -1 indicates an error (matrix couldn't be inverted)
    //====================================================================================================================
    double divergenceShape(double** cov_1, double** cov_2, unsigned long int dim, double minVariance = 0.0);

    //====================================================================================================================
    //	The Earth Mover's Distance for two Speaker Models (GMM-Like)
    //  This is a wrapper around the SC_EMD class; the EMD-implementation is the one from it's autor, Yossi Rubner
    //  See also: Y.Rubner, C.Tomasi, L.J.Guibas, "The Earth Mover's Distance as a Metric for Image Retrieval", 
    //            International Journal of Computer Vision 40(2) 99-121, Kluwer Academic Publishing, 2000
    //====================================================================================================================
    double EMD(SC_Model* pModel_1, SC_Model* pModel_2);
		static double EMD(SC_Signature* pSignature1, SC_Signature* pSignature2, SC_TweakableParameters *pTweak);

    //====================================================================================================================
    //	Beigi distance measure between 2 gaussian models with diagonal covariances
    //  The actual distance measure used internally is chosen by the centroids in the signatures towhich the models are 
		//  converted and is controlled by the tweakable parameters
    //  According to: H.S.M.Beigi, S.H.Maes, J.S.Sorensen, "A Distance Measure Between Collections of Distributions and 
    //                it's Application to Speaker Recognition"
    //====================================================================================================================
    double beigi(SC_Model* pModel_1, SC_Model* pModel_2);
    static double beigi(SC_Signature* pSignature1, SC_Signature* pSignature2);
    double beigi(SC_Signature* signature, float *vector, int dim, double *Ws, double *Wv);
  
    //====================================================================================================================
    //	The Mahalanobis distance measure for gaussians with diagonal covariances as suggested in "Couvreur, Boite, 
		//  'Speaker Tracking In Broadcast Audio Material In The Framework Of The THISL Project', 1999";
		//  we divide by var1*var2, not just var1 as in the Beigi-paper; and we normalize with 1/d
    //====================================================================================================================
		template<class T> static double mahalanobis(T* mean_1, T* mean_2, T* variance_1, T* variance_2, unsigned short int dim) {
			unsigned short int d;
			double res = 0.0;
		  
			//the mahalanobis-distance is a simple sum for a diagonal covariance matrix
			for (d = 0; d < dim; d++) {
				//assert(variance_1[d]*variance_2[d] > 0.0);
				res += ((mean_2[d] - mean_1[d]) * (mean_2[d] - mean_1[d])) / (sqrt(variance_1[d]) * sqrt(variance_2[d]));
			}
		  
			return res/(double)(d);
		}

    //====================================================================================================================
    //	The Bhattacharyya distance measure for gaussians with diagonal covariances as suggested in 
		//  "http://en.wikipedia.org/wiki/Bhattacharyya_distance";
    //====================================================================================================================
		template<class T> static double bhattacharyya(T* mean_1, T* mean_2, T* variance_1, T* variance_2, unsigned short int dim) {
			unsigned short int d;
			double res1 = 0.0, res2 = 1.0, res3 = 1.0;
		  
			//the bhattacharyya-distance is a simple sum for a diagonal covariance matrix
			for (d = 0; d < dim; d++) {
				//assert(variance_1[d]*variance_2[d] > 0.0);
				res1 += ((mean_1[d]-mean_2[d]) * (mean_1[d]-mean_2[d])) / variance_1[d];
				res2 *= (variance_1[d]+variance_2[d]) / 2.0;
				res3 *= variance_1[d] * variance_2[d];
			}
		  
			return 0.125*res1 + 0.5*sclib::ln(res2/sqrt(res3));
		}

    //====================================================================================================================
    //  The Kullback-Leibler distance as suggested in 
		//	"http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Kullback.E2.80.93Leibler_divergence"
    //====================================================================================================================
		template<class T> static double kullbackLeibler(T* mean_1, T* mean_2, T* variance_1, T* variance_2, unsigned short int dim) {
			unsigned short int d;
			double res1 = 0.0, res2 = 0.0, res3 = 1.0, res4 = 1.0;
		  
			for (d = 0; d < dim; d++) {
				//corrupt formula from "Couvreur, Boite, 'Speaker Tracking In Broadcast Audio Material In The Framework Of The THISL Project', 1999":
				//assert(variance_1[d] > 0.0 && variance_2[d] > 0.0);
				//res1 += ((mean_1[d] - mean_2[d]) * (mean_1[d] - mean_2[d])) / (variance_1[d] + variance_2[d]);
				//res2 += (variance_2[d] / variance_1[d]) + (variance_1[d] / variance_2[d]) - 2.0;
				
				//corrupted (??) formula from "Unsupervised Speaker Change Detection for Broadcast News Segmentation", Jorgensen, Molgaard, Hansen, 2006:
				//res1 += (variance_1[d]-variance_2[d]) / (variance_2[d]-variance_1[d]);
				//res2 += ((mean_1[d]-mean_2[d])*(mean_1[d]-mean_2[d])) / (variance_1[d]+variance_2[d]);

				res3 *= variance_2[d]; //det. var2
				res4 *= variance_1[d]; //det. var1
				res1 += variance_1[d] / variance_2[d]; //trace var2^-1 * var1
				res2 += ((mean_2[d]-mean_1[d])*(mean_2[d]-mean_1[d])) / variance_2[d];
			}

			return 0.5 * (sclib::ln(res3/res4) + res1 + res2 - (double)(dim));  //0.5*res1 + 0.5*res2;
		}

    //====================================================================================================================
    //	The std. Euclidean Distance measure for the Beigi distance
    //====================================================================================================================
    template<class T> static double euclid(T* x, T* y, unsigned short int dim) {
      unsigned short int d;
      double res = 0.0;
      
      for (d = 0; d < dim; d++) {
        res += (double)((x[d] - y[d]) * (x[d] - y[d]));
      }
      res = sqrt(res);
      
      return res;
    }

    //====================================================================================================================
    //	The squared Eculidean distance or "squared error" between x and y
    //====================================================================================================================
		template<class T> static double squaredError(T* x, T* y, unsigned short int dim) {
      unsigned short int d;
      double res = 0.0;
      
      for (d = 0; d < dim; d++) {
        res += (double)((x[d] - y[d]) * (x[d] - y[d]));
      }
      
      return res;
		}

    //====================================================================================================================
    //	The std. Euclidean Distance measure between two integral-type scalars
    //====================================================================================================================
    template<class T> static double euclid(T x, T y) {
      return sqrt((double)((x-y) * (x-y)));
    }

    //====================================================================================================================
    //	The std. Euclidean Distance measure between two vectors in libSVM format
    //====================================================================================================================
    static double euclid(SC_SVMnode *x, SC_SVMnode *y);

    //====================================================================================================================
    //	Compute all pairwise distances of an integral-type matrix and store them in a flat vector of the  following form:
    //    usually, it would be stored in a matrix at position (t, tt), that is t*T + tt
    //    but all positions with x are not present:
    //      x 0 1 2 3
    //      x x 4 5 6
    //      x x x 7 8
    //      x x x x 9
    //      x x x x x
    //   so, (t+1)*(t+2)/2 counts the x's that need to be subtracted
    //====================================================================================================================
    template<class T> double* getPairwiseDistances(T data, unsigned long int length, unsigned short int dim) {
      double *distances = NULL;
      unsigned long int t, tt, T = length, D = dim;

      //make room for n*(n-1)/2 pairwise distances
      MArray_1D(distances, T*(T-1)/2, double, "getPairwiseDistances: distances");

      for (t = 0; t < T-1; t++) {
        for (tt = t+1; tt < T; tt++) {
          distances[t*T + tt - (t+1)*(t+2)/2] = euclid(data[t], data[tt], dim);
        }
      }

      return distances;
    }

    //====================================================================================================================
    //	Compute all pairwise distances of one column of an integral-type matrix and store them in a flat vector of the 
    //  following form:
    //    usually, it would be stored in a matrix at position (t, tt), that is t*T + tt
    //    but all positions with x are not present:
    //      x 0 1 2 3
    //      x x 4 5 6
    //      x x x 7 8
    //      x x x x 9
    //      x x x x x
    //   so, (t+1)*(t+2)/2 counts the x's that need to be subtracted
    //====================================================================================================================
    template<class T> double* getPairwiseDistances(T data, unsigned long int length, unsigned long int dim, unsigned long int whichDimension) {
      double *distances = NULL;
      unsigned long int t, tt, T = length, d = whichDimension, D = dim;

      //make room for n*(n-1)/2 pairwise distances
      MArray_1D(distances, T*(T-1)/2, double, "getPairwiseDistances: distances");

      for (t = 0; t < T-1; t++) {
        for (tt = t+1; tt < T; tt++) {
          distances[t*T + tt - (t+1)*(t+2)/2] = euclid(data[t][d], data[tt][d]);
        }
      }

      return distances;
    }

    //====================================================================================================================
    //	2 methods to deal with a 'flat' matrix as created by the above getPairWiseDistance()-method
    //====================================================================================================================
    void set1DMatrixValue(double *vectorMatrix, unsigned long int length, unsigned long int row, unsigned long int column, double value);
    double get1DMatrixValue(double *vectorMatrix, unsigned long int length, unsigned long int row, unsigned long int column);

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
		double svmArcDistance(SC_SVM *pSVM, SC_SVMmodel *oneClassSVM_1, SC_SVMmodel *oneClassSVM_2);

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
		double newWeightedEuclideanDistance(SV_Data *pSegment_1, SV_Data *pSegment_2);

		//====================================================================================================================
		//	The T^2 statistic is a 1st order statistical test if the mean of the left side (0..b-1) equals the mean of the 
		//  right side (b..T-1), gaussianity assumed. If b is >0, the value of T^2 for each possible b between the two 
		//  segments is returned, otherwise only the specified value, each time in a newly allocated array. See also Zhou, 
		//  Hansen, "Efficient Audio Stream Segmentation via the Combined T^2 Statistic and Bayesian Information Criterion", 
		//  2005
		//====================================================================================================================
		double* tSquareStatistic(SV_Data *pSegment_1, SV_Data *pSegment_2, int b = -1);

		//====================================================================================================================
		//	The T^2 statistic is a 1st order statistical test if the mean of the left side (0..b-1) equals the mean of the 
		//  right side (b..T-1), gaussianity assumed. The T^2 value for the specified position b is returned. 
		//  Because this method is ofton called several times in a row for increasing values of b, a method has been 
		//  implemented for speed up: the internally computed means and the covariance are returned to the caller (who needs
		//  to destruct those 3!); when the mehtod is called for the next b value (next means: b=b+1), those are only updated
		//  and note computed from scratch! See also Zhou, Hansen, "Efficient Audio Stream Segmentation via the Combined T^2 
		//  Statistic and Bayesian Information Criterion", 2005
		//====================================================================================================================
		template<class T> double tSquareStatistic(T **data, int len, int dim, unsigned int b, double* &leftMean, double* &rightMean, double** &completeInvertedCov) {
			int i, j;
			double t2 = 0.0, tmp;
			double *m = NULL, *m1 = NULL, *m2 = NULL, **cov = NULL, *mDiff = NULL;
			
			if (this->pMatrixFunc == NULL) {
				this->pMatrixFunc = new SC_MatrixFunctions();
				this->linkedMatrixFunc = false;
			}

			//compute both side means only for the 1st time completely, after that the cache from last round can be used and only needs update of 1 vector
			if (leftMean != NULL) {
				m1 = leftMean;
				for (i = 0; i < dim; i++) { //update old mean by adding the next vector
					m1[i] = ((double)(b-1)/(double)(b))*m1[i] + (1.0/(double)(b))*data[b-1][i];
				}
			} else {
				m1 = this->pMatrixFunc->mean(data, len, dim, 0, b);
			}
			if (rightMean != NULL) {
				m2 = rightMean;
				for (i = 0; i < dim; i++) { //update old mean by removing the leftmost vector
					m2[i] = ((double)(len-b+1)/(double)(len-b))*m2[i] - (1.0/(double)(len-b))*data[b-1][i];
				}				
			} else {
				m2 = this->pMatrixFunc->mean(data, len, dim, b, len);
			}

			//use complete cov from last run or compute it for the first time
			if (completeInvertedCov != NULL) {
				cov = completeInvertedCov;
			} else {
				MArray_1D(m, dim, double, "SC_DistanceMeasures.tSquareStatistic: m");
				for (i = 0; i < dim; i++) {
					m[i] = ((double)(b)/(double)(len))*m1[i] + ((double)(len-b)/(double)(len))*m2[i];
				}
				cov = this->pMatrixFunc->cov(data, len, dim, m);
				this->pMatrixFunc->inv(cov, dim); //cov := cov^-1
				MFree_1D(m);
			}

			//compute mean difference
			MArray_1D(mDiff, dim, double, "SC_DistanceMeasures.tSquareStatistic: mDiff");
			for (i = 0; i < dim; i++) {
				mDiff[i] = m1[i] - m2[i];
			}

			//this is the t^2 computation
			for (i = 0; i < dim; i++) { //(m1-m2)' * cos^-1 * (m1-m2)
				tmp = 0.0;
				for (j = 0; j < dim; j++) {
					tmp += mDiff[j] * cov[i][j];
				}
				t2 += tmp * mDiff[i];
			}
			t2 *= (double)((b+1)*(len-b+1)) / (double)(len); //+1 because in the original paper, the segments are from 1..b and b+1..T instead of zero-based
			MFree_1D(mDiff);

			//do some caching to speed up next round
			leftMean = m1;
			rightMean = m2;
			completeInvertedCov = cov;

			return t2;
		}

		//====================================================================================================================
		//	The T^2 statistic is a 1st order statistical test if the mean of the left side (0..b-1) equals the mean of the 
		//  right side (b..T-1), gaussianity assumed. If b is <0, the value of T^2 for each possible b between the two 
		//  segments is returned, otherwise only the specified value, each time in a newly allocated array. By start and end 
		//  the segment can be specied that is tested; by setting middle1 and middle2 >0 and start<middle1<middle2<end a 
		//  segment can be constructed without the middle part, i.e. start->middle1, middle2->end. See also Zhou, Hansen, 
		//  "Efficient Audio Stream Segmentation via the Combined T^2 Statistic and Bayesian Information Criterion", 2005
		//====================================================================================================================
		double* tSquareStatistic(SV_Data *pSegment, int b = -1, int start = -1, int middle1 = -1, int middle2 = -1, int end = -1);

		//====================================================================================================================
		//	Gives the cross likelihood ratio between two mutlivariate, diagonal-covariance gaussians fitted to both datasets
		//  of equal size according to the paper "A Simple but Effective Approach to Speaker Tracking in Broadcast News", 
		//  Rodríguez, Penagarikano, Bordel, 2007
		//====================================================================================================================
		double crossLikelihoodRatio(SV_Data *pData_1, SV_Data *pData_2);

    //====================================================================================================================
    //	The information measure as defined in "Data Mining", Witten & Frank, 2005, pp100
    //====================================================================================================================
		static double information(double p1, double p2);
		static double information(double *p_i, int dim);
		static double information(double p11, double p12, double p21, double p22);

    //====================================================================================================================
		//	The MDL criterion as defined in "Data Mining", Witten & Frank, 2005, pp301, for the task of evaluating a split in 
		//  attribute discretization: gain is the information gain gained by the split, N is the number of instances, k the 
		//  number of classes, k1 and k2 the number of classes appearing in the first and second half created by the split, E 
		//  the entropy of the instances, E1 and E2 the entropy of the instances in the two subintervals; the boolean return 
		//  value indicates if this split should be done
    //====================================================================================================================
		static bool splittingMDL(double gain, int N, int k, int k1, int k2, double E, double E1, double E2);

		//====================================================================================================================
		//	access to protected member
		//====================================================================================================================
		SC_MatrixFunctions* getMatrixFunc(void);
};

#endif
