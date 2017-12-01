/**************************************************************************/
/*    Derived from:																												*/
/*      - GN_Matrix																												*/
/*			                                                                  */
/*    Matrices are meant to have the form m[len][dim], where len is       */
/*    refereed to as y, dim is refereed to as x                           */
/*                                                                        */
/*    Responsibility:																											*/
/*      - encapsulates some commonly used mathematical functions					*/
/*			- most of them are matrix-related, so this class is derived from	*/
/*				GN_Matrix to augment it's functionality													*/
/*			- there are also some statistical functions like covariance				*/
/*				and mean calculation																						*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.03.2004																								*/
/**************************************************************************/

#ifndef __SC_MatrixFunctions_H__
#define __SC_MatrixFunctions_H__

#include <assert.h>
#include <math.h>
#include <limits>
#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include <GN_Matrix.h>
#include <SV_Error.h>
#include <SV_DataIO.h>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

class SCLIB_API SC_MatrixFunctions {

	private :

	protected:

		GN_Matrix *pSolver; //solves more sophisticated matrix tasks like det, inv, svd, ...
		
		//====================================================================================================================
		// to hide the sclib namespace in this header
		//====================================================================================================================
		double getRandom(double mean, double variance);

	public :

		SC_MatrixFunctions();
		virtual ~SC_MatrixFunctions();

		//====================================================================================================================
		// Interface to more sophisticated methods provided by the GN_Matrix object
		//====================================================================================================================
		double det(double** squareMatrix, int dim);
		int inv(double** squareMatrix, int dim);
		void eigenDecomposition(double **symmetricMatrix, int dim, double **eigenVectors, double *eigenValues);

		//====================================================================================================================
		// This version of eigenDecomposition overwrites the covariance-matrix to save the eigentvectors (thus saving some 
		// memory); the code herein is merely copied from GN_Matrix.Eigen()
		//====================================================================================================================
		void eigenDecomposition(double **symmetricMatrix, int dim, double *eigenValues);

    ////====================================================================================================================
    ////  Symmetrizes a square matrix by copying the upper triangle upon the lower one
    ////====================================================================================================================
    //void      symmetrize(double** matrix, unsigned long int dim);

		//====================================================================================================================
		// Returns the log(det(squareMatrix)) via the Cholesky decomposition when squareMatrix is near singular and computing
		// it the standard way would produce log(0); squareMatrix needs to be positive definite, like e.g. a covariance matrix
		//====================================================================================================================
		double logDet(double **squareMatrix, int dim, bool checkAssumptions = false);

		//====================================================================================================================
    //  Cholesky decomposition: A=LL^T as described in "numerical recipes in c"
    //  returns the Cholesky factor L, a lower triangular matrix also known as A's square-root
    //  A must be symmetric and positive definite; if checkAssumptions==true, this properties are first checked 
    //  and in case of a negative result NULL is returned.
		//====================================================================================================================
    double** choleskyDecomposition(double** matrix, unsigned long int dim, bool checkAssumptions = false);

    //====================================================================================================================
    //  QR decomposition: A=QR where A is len*dim, Q (len*dim) is orthogonal and R (dim*dim) is upper triangular
    //  taken from (but modifyed):
    //====================================================================================================================
    void qrDecomposition(double** A, unsigned long int len, unsigned long int dim, double** &Q, double** &R);
    
    //====================================================================================================================
    //  This function provides an interface to the SV_Lib SVD-function, which also returns the len*dim-matrix U instead 
    //  of just U*S. If U*S is not invertible, U is returned as NULL; calculating U can be switched off by setting needU
		//  to false.
    //====================================================================================================================
    void SVD(double** M, double** &U, double** &US, double** &S, double** &V, unsigned long int len, unsigned long int dim, bool needU = true);

    //====================================================================================================================
    //	Checks whether the (real) matrix is positive definite; it is iff all eigenvalues of  it's symmetric part are 
    //  positive
    //====================================================================================================================
    bool isPositiveDefinite(double **matrix, unsigned long int dim);

    //====================================================================================================================
    //	returns the symmetric part of the square matrix: Ms = 1/2 * (M + M')
    //====================================================================================================================
    template<class T> T** getSymmetricPart(T** matrix, unsigned long int dim) {
      T **symm, **mt;

      mt = transpose(matrix, dim, dim);
      symm = add(matrix, mt, dim, dim);
      mult(symm, (double)0.5, dim, dim, true);

      MFree_2D(mt);

      return symm;
    }

    //====================================================================================================================
    //	returns the anti-symmetric part of the square matrix: Ms = 1/2 * (M - M')
    //====================================================================================================================
    template<class T> T** getAntiSymmetricPart(T** matrix, unsigned long int dim) {
      T **asymm, **mt;

      mt = transpose(matrix, dim, dim);
      asymm = sub(matrix, mt, dim, dim);
      mult(asymm, (double)0.5, dim, dim, true);

      MFree_2D(mt);

      return asymm;
    }

    //====================================================================================================================
    //	Checks whether the matrix is symmetric (m_ij==m_ji forall i,j)
    //====================================================================================================================
    template<class T> bool isSymmetric(T** matrix, unsigned long int dim) {
      for (unsigned long int y = 0; y < dim; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          if (matrix[y][x] != matrix[x][y]) {
            return false;
          }
        }
      }
      return true;
    }

    //====================================================================================================================
    //	add two matrices; the space for the new matrix is allocated by this function
    //  if replace==true, m_1 then contains the new values and the pointer to it is returned (it hasn't changed)
    //====================================================================================================================
    template<class T> T** add(T** m_1, T** m_2, unsigned long int len, unsigned long int dim, bool replace = false) {
	    T** sum;
    	
      if (replace == false) {
        MArray_2D(sum, (long int)len, dim, T, "SC_MatrixFunctions.add: sum");
      } else {
        sum = m_1;
      }
	
	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    sum[y][x] = m_1[y][x] + m_2[y][x];
		    }
	    }

	    return sum;
    }

    //====================================================================================================================
    //	add a scalar to each cell of a matrix; the space for the new matrix is allocated by this function
    //  if replace==true, m_1 then contains the new values and the pointer to it is returned (it hasn't changed)
    //====================================================================================================================
    template<class T> T** add(T** m_1, T scalar, unsigned long int len, unsigned long int dim, bool replace = false) {
	    T** sum;

      if (replace == false) {
        MArray_2D(sum, (long int)len, dim, T, "SC_MatrixFunctions.add: sum");
      } else {
        sum = m_1;
      }
	    
	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    sum[y][x] = m_1[y][x] + scalar;
		    }
	    }

	    return sum;
    }

    //====================================================================================================================
    //	multiplicate two matrices componentwise; the space for the new matrix is allocated by this function
    //  if replace==true, m_1 then contains the new values and the pointer to it is returned (it hasn't changed)
    //====================================================================================================================
    template<class T> T** multComponentWise(T** m_1, T** m_2, unsigned long int len, unsigned long int dim, bool replace = false) {
	    T** prod;

      if (replace == false) {
        MArray_2D(prod, (long int)len, dim, T, "SC_MatrixFunctions.mult: prod");
      } else {
        prod = m_1;
      }

      for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    prod[y][x] = m_1[y][x] * m_2[y][x];
		    }
	    }

	    return prod;
    }

    //====================================================================================================================
    //	multiplicate two matrices in the mathematical sense; the space for the new matrix is allocated by this function
    //  Falk-Scheme for this:
    //                dim2
    //             _________
    //            |        | len2 = dim1
    //        ___ |________|
    //       |  |
    //  len1 |  |
    //       |__|
    //       dim1
    //
    //  replaces original matrix m_1 and returns NULL if replace is switched on
    //====================================================================================================================
    template<class T> T** mult(T** &m_1, T** m_2, unsigned long int len1, unsigned long int dim1, unsigned long int len2, unsigned long int dim2, bool replace = false) {
	    T** prod;
    	
      assert(dim1 == len2);
	    MArray_2D(prod, (long int)len1, dim2, T, "SC_MatrixFunctions.mult: prod");

	    for (unsigned long int y = 0; y < len1; y++) { //rows of m_1
		    for (unsigned long int x = 0; x < dim2; x++) { //cols of m_2
          prod[y][x] = 0.0;
          for (unsigned long int d = 0; d < dim1; d++) {
            prod[y][x] += m_1[y][d] * m_2[d][x];
          }
		    }
	    }

      if (replace == true) {
        MFree_2D(m_1);
        m_1 = prod;
        prod = NULL;
      }

	    return prod;
    }

    //====================================================================================================================
    //	same as above, but with m_2 being a column vector
    //====================================================================================================================
    template<class T> T* mult(T** m, T* v, unsigned long int len1, unsigned long int dim1, unsigned long int dim2) {
	    T* prod;
    	
      assert(dim1 == dim2);
	    MArray_1D(prod, len1, T, "SC_MatrixFunctions.mult: prod");
      

	    for (unsigned long int y = 0; y < len1; y++) { 		    
        prod[y] = dotProduct(m[y], v, dim1);
	    }

	    return prod;
    }

    //====================================================================================================================
    //	multiplicate a matrix with a scalar; the space for the new matrix is allocated by this function
    //  in case of replace==true, the original matrix-elements are altered and the (untouched) pointer to m_1 is returned
    //====================================================================================================================
    template<class T> T** mult(T** m_1, T scalar, unsigned long int len, unsigned long int dim, bool replace = false) {
      T** prod;

      if (replace == false) {
        MArray_2D(prod, (long int)len, dim, T, "SC_MatrixFunctions.mult: prod");
      } else {
        prod = m_1;
      }

	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    prod[y][x] = m_1[y][x] * scalar;
		    }
	    }

      return prod;
    }

		//====================================================================================================================
		//	cov, mean and std calculation are pretty much the same than in Jialong's SV_Model-class, but instead of a 
		//	SV_Data-object, they just work with arbitrary integral**-matrices; names have been changed to refelct sclib naming 
		// conventions, though
		//====================================================================================================================

    //====================================================================================================================
    //	function to calculate covariance-matrix from a given set of feature-vectors (organized in a float-matrix)
    //	memory for cov is allocated by this function
    //  TODO: normalize by len-1 instead of len for the best unbiased estimator for the covar?!? what effects does this 
    //        have on code already written & existent & "in production"?
    //  => DONE, any negative effects???
    //  with the parameters min/max-len/dim, a lower-dimensional submatrix can be selected to calculate it's covariance
    //====================================================================================================================
    template<class T> double** cov(T** features, unsigned long int len, unsigned long int dim, double* meanV = NULL, unsigned long int minLen = 0, unsigned long int maxLen = 0, unsigned long int minDim = 0, unsigned long int maxDim = 0, bool diagonal = false) {
	    double **cov, *mean, *diff, norm;
	    unsigned long int y;
	    unsigned long int x, x_2;
      bool meanCalculated;
      unsigned long int mxDim = maxDim, mxLen = maxLen, length, dimension;

      if (mxLen == 0) {mxLen = len;}
      if (mxDim == 0) {mxDim = dim;}
      length = mxLen - minLen; 
			dimension = mxDim - minDim;
			norm = (double)(length-1);

	    MArray_1D(diff, dimension, double, "SC_MatrixFunctions.cov: diff");
      cov = zeros(dimension, dimension);
      if (meanV == NULL) {
        mean = this->mean(features, len, dim, minLen, maxLen, minDim, maxDim);
        meanCalculated = true;
      } else {
        mean = meanV;
        meanCalculated = false;
      }

	    for (y = minLen; y < mxLen; y++) {
				if (diagonal == false) {
					for (x = minDim; x < mxDim; x++) {
						diff[x-minDim] = (double)(features[y][x]) - mean[x-minDim];
					}
					for (x = minDim; x < mxDim; x++) {
						for (x_2= minDim ; x_2 < mxDim; x_2++) {
							cov[x-minDim][x_2-minDim] += (diff[x-minDim] * diff[x_2-minDim]) / norm;
						}
					}
				} else {
					for (x = minDim; x < mxDim; x++) {
						cov[x-minDim][x-minDim] += pow((double)(features[y][x]) - mean[x-minDim], 2.0) / norm;
					}
				}
	    }

      if (meanCalculated == true) {
				MFree_1D(mean);
			}
	    MFree_1D(diff);

	    return cov;
    }

    //====================================================================================================================
    //	This is a function for recursively compute the sample covariance matrix and sample mean: We already have estimates
    //  for the first N samples and no want to update them with the next sample
    //  From: Campbell, "Speaker Recognition - A Tutorial", IEEE, 1997
    //  TODO: Not tested yet, cov-formula correct???
    //====================================================================================================================
    template<class T> void updateCovAndMean(double** oldCov, double* oldMean, unsigned long int oldN, T* features, unsigned long int len, unsigned long int dim) {
      //update cov (use the old mean for this!)
      //cov_(n+1) = ((n-1)/n) * cov_n  +  (1/(n+1)) * (x_(n+1) - mean_n) * transpose(x_(n+1) - mean_n)
      mult(oldCov, (double)(oldN-1)/(double)(oldN), dim, dim, true);
      add(oldCov, (1.0/(double)(oldN+1))*dotProduct(features, oldMean, dim), dim, dim, true);

      //update mean
      //mean_(n+1) = mean_n  +  (1/(n+1)) * (x_(n+1) - mean_n) ( = mean_(n+1) = n/(n+1)*mean_n + 1/(n+1)*x_(n+1) )
      double *diff = *(sub(&features, &oldMean, 1, dim));
      mult(&diff, 1.0/(double)(oldN+1), 1, dim, true);
      add(&oldMean, &diff, 1, dim, true);

	    return;
    }

    //====================================================================================================================
    //	function to calculate mean-vector from a given set of feature-vectors (organized in a float-matrix)
    //	memory for mean is allocated by this function
    //  with the parameters min/max-len/dim, a lower-dimensional submatrix can be selected to calculate its mean
    //====================================================================================================================
    template<class T> double* mean(T** features, unsigned long int len, unsigned long int dim, unsigned long int minLen = 0, unsigned long int maxLen = 0, unsigned long int minDim = 0, unsigned long int maxDim = 0) {	
	    double* mean;
      unsigned long int mxDim = maxDim, mxLen = maxLen, length, x;

      if (mxLen == 0) {mxLen = len;}
      if (mxDim == 0) {mxDim = dim;}
      length = mxLen - minLen; 

	    MArray_1D(mean, mxDim - minDim, double, "SC_MatrixFunctions.mean: mean");

	    for (x = 0; x < (mxDim - minDim); x++) {
		    mean[x] = 0.0;
	    }

	    for (unsigned long int y = minLen; y < mxLen; y++) {
		    for (x = minDim; x < mxDim; x++) {
			    mean[x-minDim] += (double)(features[y][x]) / (double)length;
		    }
	    }

			return mean;
    }

    //====================================================================================================================
    //	function to calculate mean-vector from a given set of features (organized in a float-vector)
    //  with the parameters min/max-dim, a lower-dimensional subvector can be selected to calculate it's mean
    //====================================================================================================================
    template<class T> double mean(T* features, unsigned long int dim, unsigned long int minDim = 0, unsigned long int maxDim = 0) {	
	    double mean = 0.0;
      unsigned long int mxDim = maxDim, length;

      if (mxDim == 0) {
				mxDim = dim;
			}
			length = mxDim - minDim + 1;

	    for (unsigned long int y = minDim; y < mxDim; y++) {
		    mean += (double)(features[y]) / (double)(length);
	    }

	    return mean;
    }

    //====================================================================================================================
    //	function to calculate standard-deviation-vector (=sqrt(variance) from a given set of feature-vectors
    //	memory for std is allocated by this function
    //  with the parameters min/max-len/dim, a lower-dimensional submatrix can be selected to calculate it's mean
    //====================================================================================================================
    template<class T> double* std(T** features, unsigned long int len, unsigned long int dim, double* meanV = NULL, unsigned long int minLen = 0, unsigned long int maxLen = 0, unsigned long int minDim = 0, unsigned long int maxDim = 0) {	
	    double *std, *mean, xc, length1;
      bool meanCalculated;
      unsigned long int mxDim = maxDim, mxLen = maxLen, length, x;

      if (mxLen == 0) {mxLen = len;}
      if (mxDim == 0) {mxDim = dim;}
      length = mxLen - minLen;
			length1 = (double)(length - 1);

	    MArray_1D(std, mxDim - minDim, double, "SC_MatrixFunctions.std: std");
      if (meanV == NULL) {
        mean = this->mean(features, len, dim, minLen, maxLen, minDim, maxDim);
        meanCalculated = true;
      } else {
        mean = meanV;
        meanCalculated = false;
      }

      for (x = 0; x < (mxDim - minDim); x++) {
		    std[x] = 0.0;
	    }

	    for (unsigned long int y = minLen; y < mxLen; y++) {
		    for (x = minDim; x < mxDim; x++) {
          xc = (double)(features[y][x]) - mean[x-minDim];
			    std[x-minDim] += xc * xc / length1;
		    }
	    }

      for (x = 0; x < (mxDim - minDim); x++) {
        std[x] = sqrt(std[x]);
      }

      if (meanCalculated == true) {
				MFree_1D(mean);
			}

	    return std;
    }

    //====================================================================================================================
    //	function to calculate variance-vector (=sd^2) from a given set of feature-vectors
    //	memory for the variance is allocated by this function
    //====================================================================================================================
    template<class T> double* variance(T** features, unsigned long int len, unsigned long int dim, double* meanV = NULL) {	
	    double *variance, *mean, xc;
      bool meanCalculated;
			double len1 = (double)(len - 1);

	    MArray_1D(variance, dim, double, "SC_MatrixFunctions.variance: std");
      if (meanV == NULL) {
        mean = this->mean(features, len, dim);
        meanCalculated = true;
      } else {
        mean = meanV;
        meanCalculated = false;
      }

      for (unsigned long int x = 0; x < dim; x++) {
		    variance[x] = 0.0;
	    }

	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
          xc = (double)(features[y][x]) - mean[x];
			    variance[x] += (xc * xc) / len1;
		    }
	    }

      if (meanCalculated == true) {
				MFree_1D(mean);
			}

	    return variance;
    }

    //====================================================================================================================
    //	normally, all this gaussian stuff is handled in SC_Gauss, but only for the diagonal covariance case. this 
		//  multivariate gaussian full-covariance log-likelihood is an exception! (TODO)
		//  if featuresBoreParameters==true (i.e. exactly these feature where used to estimate mean and cov), things can be 
		//  speeded up.
    //====================================================================================================================
		template<class T> double gaussLogLikelihood(T** features, unsigned long int len, unsigned long int dim, double *mean, double **cov, bool featuresBoreParameters = false) {
			double ll = 0.0;
			double detCov = det(cov, dim);
			double constPart = -1.0*(double)(len) * (((double)(dim)/2.0)*1.8378770664093454835606594728112 + 0.5*log(detCov)); //-#features * (d/2*log(2pi) + 1/2*log|cov|)
			double **invCov, *diff, *tmp;
			
			if (featuresBoreParameters == true) {//-#features * (d/2*log(2pi) + 1/2*log|cov|) - 1/2*dim
				ll = constPart - 0.5*(double)(dim)*(double)(len-1); //the exp() part reduces mathematically exact to this, yeah!
			} else { //-#features * (d/2*log(2pi) + 1/2*log|cov|) - 1/2*SUM[(x-mean)'*cov^-1*(x-mean)]
				if (detCov > 0) {
					invCov = copy(cov, dim, dim);
					inv(invCov, dim);
					MArray_1D(diff, dim, double, "SC_MatrixFunctions.gaussLogLikelihood: diff");
					MArray_1D(tmp, dim, double, "SC_MatrixFunctions.gaussLogLikelihood: tmp");

					for (unsigned long int y = 0; y < len; y++) { //ll = SUM[(x-mean)'*cov^-1*(x-mean)]
						sub(features[y], mean, diff, dim);
						for (unsigned long int x = 0; x < dim; x++) {
							tmp[x] = dotProduct(invCov[x], diff, dim);
						}
						ll += dotProduct(diff, tmp, dim);
					}

					MFree_2D(invCov);
					MFree_1D(tmp);
					MFree_1D(diff);

					ll = constPart - 0.5*ll;
				} else {
					ll = constPart;
				}
			}

			return ll;
		}

    //====================================================================================================================
    //	function to calculate variance (=sd^2) from a given set of features
		//  meanV should be a pointer to a scalar; it's just a pointer so that NULL can be passed to indicate there is no 
		//  pre-computed mean
    //====================================================================================================================
    template<class T> double variance(T* features, unsigned long int dim, double* meanV = NULL) {	
	    double variance = 0.0, mean, xc;
			double dim1 = (double)(dim - 1);

      if (meanV == NULL) {
        mean = this->mean(features, dim);
      } else {
        mean = *meanV;
      }

	    for (unsigned long int y = 0; y < dim; y++) {
        xc = (double)(features[y]) - mean;
			  variance += xc * xc / dim1;
	    }

	    return variance;
    }

    //====================================================================================================================
    //	return the MLE of the inverse gaussian's lambda parameter per dimension; the second parameter of this pdf is the 
		//  sample mean as also used by the gaussian distribution
		//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
    //====================================================================================================================
    template<class T> double* lambda(T** features, unsigned long int len, unsigned long int dim, double* meanV = NULL) {	
	    double *lambda, *mean;
      bool meanCalculated;
			double invD = 1.0 / (double)(len);

	    MArray_1D(lambda, dim, double, "SC_MatrixFunctions.lambda: lambda");
      if (meanV == NULL) {
        mean = this->mean(features, len, dim);
        meanCalculated = true;
      } else {
        mean = meanV;
        meanCalculated = false;
      }

      for (unsigned long int x = 0; x < dim; x++) {
		    lambda[x] = 0.0;
	    }

	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
          lambda += invD * (1.0/(double)(features[y][x]) - 1.0/mean[x]);
		    }
	    }

			for (unsigned long int x = 0; x < dim; x++) {
				if (lambda[x] > 0) {
					lambda[x] = 1.0 / lambda[x];
				} else {
					lambda[x] = 0.0;
				}
			}

      if (meanCalculated == true) {
				MFree_1D(mean);
			}

	    return lambda;
    }

    //====================================================================================================================
    //	return the MLE of the inverse gaussian's lambda parameter; the second parameter of this pdf is the sample mean
		//  as also used by the gaussian distribution
		//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
    //====================================================================================================================
		template<class T> double lambda(T* features, unsigned long int dim, double *meanV = NULL) {
	    double lambda = 0.0, mean;
			double invD = 1.0 / (double)(dim);

      if (meanV == NULL) {
        mean = this->mean(features, dim);
      } else {
        mean = *meanV;
      }

	    for (unsigned long int y = 0; y < dim; y++) {
        lambda += invD * (1.0/(double)(features[y]) - 1.0/mean);
	    }

			if (lambda > 0) {
				lambda = 1.0 / lambda;
			} else {
				lambda = 0.0;
			}

	    return lambda;
		}

    //====================================================================================================================
    //	return the standard deviation according to the parameters of the given inverse gaussian distribution
		//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
    //====================================================================================================================
		double igStd(double mean, double lambda);	
		double* igStd(double *mean, double *lambda, unsigned long int dim);

    //====================================================================================================================
    //	return the log-likelihood of the feature given the parameters of an inverse gaussian (ig) pdf
		//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
    //====================================================================================================================
		template<class T> double igLogLikelihood(T feature, double mean, double lambda) {
			double difference = feature-mean;

			return log(sqrt(lambda / (6.283185307179586476925286766559*feature*feature*feature))) + ((-1.0*lambda*difference*difference) / (2.0*mean*mean*feature)); //sclib::two_pi
		}

    //====================================================================================================================
    //	return the log-likelihood of the feature-vector given the parameters of an inverse gaussian (ig) pdf
		//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
    //====================================================================================================================
		template<class T> double igLogLikelihood(T *feature, double *mean, double *lambda, unsigned long int dim) {
			double difference;
			double ll = 0.0;

			for (unsigned long int x = 0; x < dim; x++) {
				difference = feature[x] - mean[x];
				ll += log(sqrt(lambda[x] / (6.283185307179586476925286766559*feature[x]*feature[x]*feature[x]))) + ((-1.0*lambda[x]*difference*difference) / (2.0*mean[x]*mean[x]*feature[x])); //sclib::two_pi
			}

			return ll;
		}

    //====================================================================================================================
    //	matrix-trace-function ("Spur")
    //	returns the sum of all elements on the main diagonal of a square matrix
    //====================================================================================================================
    template<class T> T trace(T** matrix, int dim) {
	    T sum = (T)(0.0);

	    for (int x = 0; x < dim; x++) {
		    sum += matrix[x][x];
	    }

	    return sum; 
    }

    //====================================================================================================================
    //	subtract two matrices; the space for the new matrix is allocated by this function
    //====================================================================================================================
    template<class T> T** sub(T** m_1, T** m_2, unsigned long int len, unsigned long int dim) {
	    T** diff;
    	
	    MArray_2D(diff, (long int)len, dim, T, "SC_MatrixFunctions.sub: diff");

			for (unsigned long int y = 0; y < len; y++) {
				for (unsigned long int x = 0; x < dim; x++) {
					diff[y][x] = m_1[y][x] - m_2[y][x];
				}
			}

			return diff;
    }

    //====================================================================================================================
    //	subtract two vectors; the space for the new vector is allocated by this function
    //====================================================================================================================
    template<class T> T* sub(T* v_1, T* v_2, unsigned long int dim) {
	    T* diff;
    	
	    MArray_1D(diff, dim, T, "SC_MatrixFunctions.sub: diff");

			for (unsigned long int x = 0; x < dim; x++) {
				diff[x] = v_1[x] - v_2[x];
			}

			return diff;
    }

    //====================================================================================================================
    //	subtract two vectors and place the result in the third
    //====================================================================================================================
    template<class C, class T> void sub(C* v_1, T* v_2, T* v_res, unsigned long int dim) {
			for (unsigned long int x = 0; x < dim; x++) {
				v_res[x] = v_1[x] - v_2[x];
			}

			return;
    }

    //====================================================================================================================
    //	compute the dot-product of two vectors
    //====================================================================================================================
    template<class T> T dotProduct(T* v_1, T* v_2, unsigned long int dim) {
      T result = (T)(0.0);

      for (unsigned long int x = 0; x < dim; x++) {
        result += v_1[x] * v_2[x];
      }

      return result;
    }

    //====================================================================================================================
    //  the component-wise absolute values of a matrix
    //  just as an example how to use execComponentWise()
    //====================================================================================================================
    //void SC_MatrixFunctions::mabs(double** matrix, unsigned long int len, unsigned long int dim) {
    //  execComponentWise(matrix, len, dim, &fabs);
    //
    //  return;
    //}

    //====================================================================================================================
    //  Executes the function f() (getting one double argument, returning a double) for each cell of the matrix
    //====================================================================================================================
    template<class T> void execComponentWise(T** matrix, unsigned long int len, unsigned long int dim, T (*f)(T m)) {
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          matrix[y][x] = f(matrix[y][x]);
        }
      }

      return;  
    }

    //====================================================================================================================
    //  Executes the function f(arg1, arg2) (getting arg1/2 as it's argument, returning the new value for each cell
    //  if arg1==NULL, the old cell-entry is taken as the first argument
    //====================================================================================================================
    template<class T, class S> void execComponentWise(T** matrix, unsigned long int len, unsigned long int dim, T (*f)(S p1, S p2), S *arg1, S *arg2) {
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          matrix[y][x] = f((arg1 == NULL) ? matrix[x][y] : *arg1, *arg2);
        }
      }

      return;  
    }

    //====================================================================================================================
    //  return an allocated matrix with elements set according to the funxtion parameter "value"
    //====================================================================================================================
    template<class T> T** initMatrix(unsigned long int len, unsigned long int dim, T value) {
      T **matrix;

      MArray_2D(matrix, (long)(len), dim, T, "SC_MatrixFunctions.initMatrix: matrix");
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          matrix[y][x] = value;
	      }
      }

      return matrix;
    }
    //template<class T> T** zeros(unsigned long int len, unsigned long int dim, T destType) {return initMatrix(len, dim, (T)(0.0));};
    //template<class T> T** ones(unsigned long int len, unsigned long int dim, T destType) {return initMatrix(len, dim, (T)(1.0));};
    double** zeros(unsigned long int len, unsigned long int dim) {return initMatrix(len, dim, (double)(0.0));};
    double** ones(unsigned long int len, unsigned long int dim) {return initMatrix(len, dim, (double)(1.0));};
    
    //====================================================================================================================
    //sets the elemants of the matrix to zero (whatever this is in the context of the given type...)
    //====================================================================================================================
		template<class T> void clear(T** matrix, unsigned long int len, unsigned long int dim, T clearValue) {
			unsigned long int y, x;

			for (y = 0; y < len; y++) {
				for (x = 0; x < dim; x++) {
					matrix[y][x] = clearValue;
				}
			}

			return;
		}

    //====================================================================================================================
    //  return an allocated vector with elements set according to the funxtion parameter "value"
    //====================================================================================================================
    template<class T> T* initVector(unsigned long int dim, T value) {
      T *vector;

      MArray_1D(vector, dim, T, "SC_MatrixFunctions.initVector: vector");
      for (unsigned long int x = 0; x < dim; x++) {
        vector[x] = value;
	    }

      return vector;
    }
    //template<class T> T* zeros(unsigned long int dim, T destType) {return initVector(dim, (T)(0.0));};
    //template<class T> T* ones(unsigned long int dim, T destType) {return initVector(dim, (T)(1.0));};
    double* zeros(unsigned long int dim) {return initVector(dim, (double)(0.0));};
    double* ones(unsigned long int dim) {return initVector(dim, (double)(1.0));};

    //====================================================================================================================
    //make an exact copy of a matrix; the space for the result matrix is allocated by this function
    //destination is needed by the compiler to determine the type of S
    //====================================================================================================================
    template<class T, class S> S** scopy(T** matrix, unsigned long int len, unsigned long int dim, S destType) {
	    S** newMatrix;
      
	    MArray_2D(newMatrix, (long int)len, dim, S, "SC_MatrixFunctions.copy: newMatrix");

	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    newMatrix[y][x] = (S)(matrix[y][x]);
		    }
	    }

	    return newMatrix;
    }

    //====================================================================================================================
    //make an exact copy of a vector; the space for the result vector is allocated by this function
    //destination is needed by the compiler to determine the type of S
    //====================================================================================================================
    template<class T, class S> S* scopy(T* vector, unsigned long int dim, S destType) {
	    S* newVector = NULL;
      
			if (vector != NULL && dim > 0) {
				MArray_1D(newVector, dim, S, "SC_MatrixFunctions.copy: newVector");

				for (unsigned long int x = 0; x < dim; x++) {
					newVector[x] = (S)(vector[x]);
				}
			}

	    return newVector;
    }
    
    //====================================================================================================================
    //make an exact copy of a matrix; the space for the result matrix is allocated by this function
    //====================================================================================================================
    template<class T> T** copy(T** matrix, unsigned long int len, unsigned long int dim) {
	    T** newMatrix = NULL;
      
			if (matrix != NULL && len > 0 && dim > 0) {
				MArray_2D(newMatrix, (long int)len, dim, T, "SC_MatrixFunctions.copy: newMatrix");

				for (unsigned long int y = 0; y < len; y++) {
					for (unsigned long int x = 0; x < dim; x++) {
						newMatrix[y][x] = (T)(matrix[y][x]);
					}
				}
			}

	    return newMatrix;
    }

    //====================================================================================================================
    //make an exact copy of a vector; the space for the result vector is allocated by this function
    //====================================================================================================================
    template<class T> T* copy(T* vector, unsigned long int dim) {
	    T* newVector = NULL;
      
			if (vector != NULL && dim > 0) {
				MArray_1D(newVector, dim, T, "SC_MatrixFunctions.copy: newVector");

				for (unsigned long int x = 0; x < dim; x++) {
					newVector[x] = (T)(vector[x]);
				}
			}

	    return newVector;
    }

    //====================================================================================================================
    //make an exact copy of a matrix; it is placed in dest, which's pointers aren't altered!
    //====================================================================================================================
    template<class T> T** insert(T** src, unsigned long int len, unsigned long int dim, T** dest) {
	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    dest[y][x] = (T)(src[y][x]);
		    }
	    }

	    return dest;
    }

    //====================================================================================================================
    //make an exact copy of a matrix; it is placed in dest, which's pointers aren't altered!
    //====================================================================================================================
    template<class T> T* insert(T* src, unsigned long int dim, T* dest) {
		  for (unsigned long int x = 0; x < dim; x++) {
			  dest[x] = (T)(src[x]);
		  }

	    return dest;
    }

    //====================================================================================================================
    //fills each cell of the matrix with value
    //====================================================================================================================
    template<class T> void fillWithValue(T** matrix, unsigned long int len, unsigned long int dim, T value) {
	    for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
			    matrix[y][x] = value;
		    }
	    }

	    return;
    }

    //====================================================================================================================
    //fills each cell of the vector with value
    //====================================================================================================================
    template<class T> void fillWithValue(T* vector, unsigned long int dim, T value) {
		  for (unsigned long int x = 0; x < dim; x++) {
			  vector[x] = value;
		  }

	    return;
    }

    //====================================================================================================================
    //returns a vector containing the minima of the matrix' cols
    //====================================================================================================================
    template<class T> double* min(T** matrix, unsigned long int len, unsigned long int dim) {
      double *minimum;

      MArray_1D(minimum, dim, double, "SC_MatrixFunctions.min: minimum");
      fillWithValue(&minimum, 1, dim, std::numeric_limits<double>::max());
      for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
          if (matrix[y][x] < minimum[x]) {minimum[x] = matrix[y][x];}
		    }
	    }

	    return minimum;
    }

    //====================================================================================================================
    //returns a scalar containing the minimum of the given col of the matrix
    //====================================================================================================================
    template<class T> double minOfCol(T** matrix, unsigned long int len, unsigned long int whichDim) {
      double minimum = std::numeric_limits<double>::max();

      for (unsigned long int y = 0; y < len; y++) {
        if (matrix[y][whichDim] < minimum) {
					minimum = matrix[y][whichDim];
				}
	    }

	    return minimum;
    }

    //====================================================================================================================
    //returns a scalar containing the minimum of the vectors elements
    //if the pointer to a scalar called index is !=NULL, it will be filled with the index of the minimum value
		//if the pointer to a scalar called thresholdValue !=NULL, the minimum value of vector greater than that threshold 
		//will be returned instead of the overall minimal value
    //====================================================================================================================
    template<class T> double min(T* vector, unsigned long int dim, unsigned long int* index = NULL, T* thresholdValue = NULL) {
      double minimum;

      minimum = std::numeric_limits<double>::max();
	    for (unsigned long int x = 0; x < dim; x++) {
        if (vector[x] < minimum && (thresholdValue == NULL || vector[x] > *thresholdValue)) {
          minimum = vector[x];
          if (index != NULL) {
						*index = x;
					}
        }
	    }

	    return minimum;
    }

    //====================================================================================================================
    //returns a matrix containing the minima of the matrices elements and the value 
    //====================================================================================================================
    template<class T> T** min(T** matrix, unsigned long int len, unsigned long int dim, T value) {
      double **minimum;

      MArray_2D(minimum, len, dim, double, "SC_MatrixFunctions.min: minimum");
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          minimum[y][x] = (matrix[y][x] < value) ? matrix[y][x] : value;
		    }
      }

	    return minimum;
    }

    //====================================================================================================================
    //returns a matrix containing the elementwise minimum of the two source-matrices
    //replaces original matrix m_1 and returns NULL if replace is switched on
    //====================================================================================================================
    template<class T> T** min(T** m_1, T** m_2, unsigned long int len, unsigned long int dim, bool replace = false) {
      T **minimum;

      MArray_2D(minimum, (long)len, dim, T, "SC_MatrixFunctions.max: minimum");
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          minimum[y][x] = (m_1[y][x] < m_2[y][x]) ? m_1[y][x] : m_2[y][x];
		    }
      }

      if (replace == true) {
        for (unsigned long int y = 0; y < len; y++) {
          for (unsigned long int x = 0; x < dim; x++) {
            m_1[y][x] = minimum[y][x];
		      }
        }
        MFree_2D(minimum);
      }

	    return minimum;
    }

    //====================================================================================================================
    //returns a vector containing the maxima of the matrix' cols
    //====================================================================================================================
    template<class T> double* max(T** matrix, unsigned long int len, unsigned long int dim) {
      double *maximum;

      MArray_1D(maximum, dim, double, "SC_MatrixFunctions.max: maximum");
      fillWithValue(&maximum, 1, dim, -1.0 * std::numeric_limits<double>::max());
      for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
          if (matrix[y][x] > maximum[x]) {
						maximum[x] = matrix[y][x];
					}
		    }
	    }

	    return maximum;
    }

    //====================================================================================================================
    //returns a scalar containing the maximum of the given col of the matrix
    //====================================================================================================================
    template<class T> double maxOfCol(T** matrix, unsigned long int len, unsigned long int whichDim) {
      double maximum = -1.0 * std::numeric_limits<double>::max();

      for (unsigned long int y = 0; y < len; y++) {
        if (matrix[y][whichDim] > maximum) {
					maximum = matrix[y][whichDim];
				}
	    }

	    return maximum;
    }

    //====================================================================================================================
    //returns a scalar containing the maximum of the vectors elements
    //if the pointer to a scalar called index is !=NULL, it will be filled with the index of the minimum value
    //====================================================================================================================
    template<class T> double max(T* vector, unsigned long int dim, unsigned long int* index = NULL) {
      double maximum;

      maximum = std::numeric_limits<double>::min();
	    for (unsigned long int x = 0; x < dim; x++) {
        if (vector[x] > maximum) {
					maximum = vector[x];
          if (index != NULL) {
						*index = x;
					}
				}
	    }

	    return maximum;
    }

    //====================================================================================================================
    //returns a matrix containing the maximum of the matrices elements and the value 
    //====================================================================================================================
    template<class T> T** max(T** matrix, unsigned long int len, unsigned long int dim, T value) {
      T **maximum;

      MArray_2D(maximum, (long)len, dim, T, "SC_MatrixFunctions.max: maximum");
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          maximum[y][x] = (matrix[y][x] > value) ? matrix[y][x] : value;
		    }
      }

	    return maximum;
    }

    //====================================================================================================================
    //returns a matrix containing the elementwise maximum of the two source-matrices
    //replaces original matrix m_1 and returns NULL if replaces is switched on
    //====================================================================================================================
    template<class T> T** max(T** m_1, T** m_2, unsigned long int len, unsigned long int dim, bool replace = false) {
      T **maximum;

      MArray_2D(maximum, (long)len, dim, T, "SC_MatrixFunctions.max: maximum");
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          maximum[y][x] = (m_1[y][x] > m_2[y][x]) ? m_1[y][x] : m_2[y][x];
		    }
      }

      if (replace == true) {
        for (unsigned long int y = 0; y < len; y++) {
          for (unsigned long int x = 0; x < dim; x++) {
            m_1[y][x] = maximum[y][x];
		      }
        }
        MFree_2D(maximum);
      }

	    return maximum;
    }

    //====================================================================================================================
    //returns the index of the maximal element in the vector
    //====================================================================================================================
    template<class T> unsigned long int maxIdx(T* vector, unsigned long int dim) {
      unsigned long int i = 0, d;
      T maximum = std::numeric_limits<T>::max() * -1;

      for (d = 0; d < dim; d++) {
        if (vector[d] > maximum) {
          i = d;
					maximum = vector[d];
        }
      }

      return i;
    }

    //====================================================================================================================
    //returns the index of the minimal element in the vector
    //====================================================================================================================
    template<class T> unsigned long int minIdx(T* vector, unsigned long int dim) {
      unsigned long int i = 0, d;
      T minimum = std::numeric_limits<T>::max();

      for (d = 0; d < dim; d++) {
        if (vector[d] < minimum) {
          i = d;
					minimum = vector[d];
        }
      }

      return i;
    }

    //====================================================================================================================
    //divides each element of the matrix by value
    //====================================================================================================================
    template<class T> void divide(T** matrix, unsigned long int len, unsigned long int dim, T value) {
      assert(value != 0.0);
      
      for (unsigned long int y = 0; y < len; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          matrix[y][x] /= value;
		    }
      }

	    return;
    }

    //====================================================================================================================
    //	returns a dim*dim matrix with zero elements and the vector elements on the main diagonal
    //====================================================================================================================
    template<class T> T** diag(T* vector, unsigned long int dim) {
      double **matrix;

      MArray_2D(matrix, (long)dim, (long)dim, double, "SC_MatrixFunctions.diag: matrix");
      for (unsigned long int y = 0; y < dim; y++) {
        for (unsigned long int x = 0; x < dim; x++) {
          matrix[y][x] = (x == y) ? vector[x] : (T)(0.0);
	      }
      }

      return matrix;
    }

    //====================================================================================================================
    //	returns a vector of size dim with the elements on the main diagonal of the dimxdim-matrix
    //====================================================================================================================
    template<class T> T* diag(T** matrix, unsigned long int dim) {
      double *vector;

      MArray_1D(vector, dim, double, "SC_MatrixFunctions.diag: vector");
      for (unsigned long int y = 0; y < dim; y++) {
        vector[y] = matrix[y][y];
      }

      return vector;
    }

    //====================================================================================================================
    //	inserts the vector into the main diagonal of the matrix
    //====================================================================================================================
    template<class T> void diag(T** matrix, T* vector, unsigned long int dim) {
      for (unsigned long int y = 0; y < dim; y++) {
        matrix[y][y] = vector[y];
      }

      return;
    }

    //====================================================================================================================
    //returns a vector containing the sum of the matrix' cols (if cols==true) or the matrix' rows
    //====================================================================================================================
    template<class T> T* sum(T** matrix, unsigned long int len, unsigned long int dim, bool cols = true) {
      T *sumVec;
      long int vecDim = (cols == true) ? dim : len;

      MArray_1D(sumVec, vecDim, T, "SC_MatrixFunctions.sum: sumVec");
      fillWithValue(&sumVec, 1, vecDim, (T)(0.0));
      for (unsigned long int y = 0; y < len; y++) {
		    for (unsigned long int x = 0; x < dim; x++) {
          sumVec[(cols == true) ? x : y] += matrix[y][x];
		    }
	    }

	    return sumVec;
    }

    //====================================================================================================================
    //returns a scalar containing the sum of the vectors elements
    //====================================================================================================================
    template<class T> T sum(T* vector, unsigned long int dim) {
      T res = (T)(0.0);

		  for (unsigned long int x = 0; x < dim; x++) {
        res += vector[x];
		  }

	    return res;
    }

    //====================================================================================================================
    //returns a scalar containing the sum of the matrix' col (if col==true) or row indicated by idx
    //====================================================================================================================
    template<class T> T sum(T** matrix, unsigned long int len, unsigned long int dim, unsigned long int idx, bool col = true) {
      T res = (T)(0.0);

			if ((idx < dim && col == true) || (idx < len && col != true)) {
				for (unsigned long int z = 0; z < ((col == true) ? len : dim); z++) {
					res += (col == true) ? matrix[z][idx] : matrix[idx][z];
				}
			}

	    return res;
    }


    //====================================================================================================================
    //returns a transposed (cols and rows switched) matrix
    //replaces original matrix and returns NULL if replace is switched on
    //====================================================================================================================
    template<class T> T** transpose(T** matrix, unsigned long int len, unsigned long int dim, bool replace = false) {
      T **transposed = NULL, tmp;
      
			if (replace == false) {
				MArray_2D(transposed, (long)dim, (long)len, T, "SC_MatrixFunctions.transpose: transposed");
				for (unsigned long int y = 0; y < len; y++) {
					for (unsigned long int x = 0; x < dim; x++) {
						transposed[x][y] = matrix[y][x];
					}
				}
			} else {
        for (unsigned long int y = 0; y < len; y++) {
		      for (unsigned long int x = y+1; x < dim; x++) {
						tmp = matrix[x][y];
            matrix[x][y] = matrix[y][x];
						matrix[y][x] = tmp;
		      }
	      }
      }

	    return transposed;
    }

    //====================================================================================================================
    //returns the upper triangular part (incl main diagonal) of the matrix; all other elements are set to zero
    //the original matrix is replaced
    //====================================================================================================================
    template<class T> void triu(T** matrix, unsigned long int dim) {
      for (unsigned long int y = 1; y < dim; y++) {
		    for (unsigned long int x = 0; x < y; x++) {
          matrix[y][x] = (T)(0.0);
		    }
	    }

	    return;
    }

    //====================================================================================================================
    //returns the lower triangular part (incl main diagonal) of the matrix; all other elements are set to zero
    //the original matrix is replaced
    //====================================================================================================================
    template<class T> void tril(T** matrix, unsigned long int dim) {
      for (unsigned long int y = 0; y < dim-1; y++) {
		    for (unsigned long int x = y+1; x < dim; x++) {
          matrix[y][x] = (T)(0.0);
		    }
	    }

	    return;
    }

    //====================================================================================================================
    //returns a vector conating a copy of the specified col of the source-matrix
    //====================================================================================================================
    template<class T> T* getCol(T** matrix, unsigned long int len, unsigned long int dim, unsigned long int col) {
      T* vector;

      MArray_1D(vector, len, T, "SC_MatrixFunctions.getCol: vector");
      for (unsigned long int y = 0; y < len; y++) {
        vector[y] = matrix[y][col];
      }

      return vector;
    }

    //====================================================================================================================
    //returns a vector containg a copy of the specified row of the source-matrix
    //====================================================================================================================
    template<class T> void setCol(T** matrix, T* vector, unsigned long int len, unsigned long int dim, unsigned long int col) {
      for (unsigned long int y = 0; y < len; y++) {
        matrix[y][col] = vector[y];
      }

      return;
    }

    //====================================================================================================================
    //returns a vector containg a copy of the specified row of the source-matrix
    //====================================================================================================================
    template<class T> T* getRow(T** matrix, unsigned long int len, unsigned long int dim, unsigned long int row) {
      T* vector;

      MArray_1D(vector, dim, T, "SC_MatrixFunctions.getRow: vector");
      for (unsigned long int x = 0; x < dim; x++) {
        vector[x] = matrix[row][x];
      }

      return vector;
    }

    //====================================================================================================================
    //sets the elemants of the matrix row'th row to the values of the corresponding elements of the vector
    //====================================================================================================================
    template<class T> void setRow(T** matrix, T* vector, unsigned long int len, unsigned long int dim, unsigned long int row) {
      for (unsigned long int x = 0; x < dim; x++) {
        matrix[row][x] = vector[x];
      }

      return;
    }

    //====================================================================================================================
    //exchanges the values of row1 and row 2 of the given matrix 
    //====================================================================================================================
    template<class T> void exchangeRows(T** matrix, unsigned long int len, unsigned long int dim, unsigned long int row1, unsigned long int row2) {
			T tmp;

			if (row1 != row2) {
				for (unsigned long int x = 0; x < dim; x++) {
					tmp = matrix[row1][x];
					matrix[row1][x] = matrix[row2][x];
					matrix[row2][x] = tmp;
				}
			}

      return;
    }

    //====================================================================================================================
    //returns true if all values of the vector equal the given compareValue
    //====================================================================================================================
    template<class T> bool equals(T* vector, unsigned long int dim, T compareValue) {
			bool res = true;

			for (unsigned long int x = 0; x < dim; x++) {
				if (vector[x] != compareValue) {
					res = false;
					break;
				}
      }

      return res;
    }

		//====================================================================================================================
		//check, if the two given vectors are equal. Added from Bing Shi
		//====================================================================================================================
		template<class T> bool equals(T* vector, unsigned long int dim, T* compareVector) {
			bool res = true;

			for (unsigned long int x = 0; x < dim; x++) {
				if (vector[x] != compareVector[x]) {
					res = false;
					break;
				}
			}

			return res;
		}

    //====================================================================================================================
    //returns a new matrix m which is the concatenation of m1 and m2
    //if vertical==true, m is (len1+len2)*dim [where dim==dim1==dim2], else m is len*(dim1+dim2) [where len==len1==len2]
    //====================================================================================================================
    template<class T> T** concat(T** m1, T** m2, unsigned long int len1, unsigned long int dim1, unsigned long int len2, unsigned long int dim2, bool vertical = true) {
      unsigned long int x, y, X, Y;
      T** m;
  
      if (vertical == true) {
        if (!(dim1 == dim2)) {return NULL;}
        X = dim1;
        Y = len1 + len2;
      } else {
        if (!(len1 == len2)) {return NULL;}
        X = dim1 + dim2;
        Y = len1;
      }
      MArray_2D(m, (long)Y, (long)X, T, "SC_MatrixFunctions.concat");

      for (y = 0; y < Y; y++) {
        for (x = 0; x < X; x++) {
          if (vertical == true) {
            m[y][x] = (y < len1) ? m1[y][x] : m2[y-len1][x];
          } else {
            m[y][x] = (x < dim1) ? m1[y][x] : m2[y][x-dim1];
          }
        }
      }

      return m;
    }

    //====================================================================================================================
    // Gets as input a scatter matrix which describes the result of an classification process. In it's most general form,
    // it looks like this:
    //
    //               GroundTruth (classes: 1..ca)
    //  Hypothesized [ x_1_1  | x_1_2 | ... | x_1_ca  ] (the correct clustered CEs appear on the main diagonal of the 
    //  (clusters:   [ x_2_1  | x_2_2 |     | ...     ]  matrix, i.e. x_i_i for i=1..ca; the rest in this upper part shows
    //   1..cu)      [ ...    | ...   |     |         ]  the spread of errors)
    //               [ x_cu_1 |       |     | x_cu_ca ] (cu >= ca)
    //               [Pop_1   | ...   |     | Pop_ca  ] (overall number of CEs belonging to this class, maybe not all where clustered)
    //               [Pos_1   | ...   |     | Pos_ca  ] (overall number of CEs belonging to this class in the possible set, maybe more than where clustered)
    //
    // The last two rows are optional and can be specified as special by the parameter lastSpecialRows (should be set to 
    // the number of last rows that are to be treated in a special way). With the parameter redundantClustersAreErrors it 
    // can be specified if more than one cluster (row) for the same class (col) should be treated as an error or not. See
    // SC_Score.h for more details on how theses scatterMatrixes are to be build up.
    //
    // What is done by the method is the follwing: The classification scatter matrix is transformed to an detection 
    // scatter matrix: The interestingDim col/row-cell is regarded as TP, everything else in this row/col as FP/FN, all
    // other as cells as TN (as above, see SC_Score.h for further explanations). The lastSpecialRows are transformed to
    // rows with only cols: 1 is the inteestingDim, the other is the sum of the other cells in the original row. The 
    // result may look like this:
    //
    //               GroundTruth
    //  Hypothesized [ TP | FP ] (number of true positives, false positives of the current CE)
    //               [ FN | TN ] (number of false negatives, true negatives of the current CE)
    //               [PPop|NPop] (number of positives/negatives according to ground-truth in the population)
    //               [PPos|NPos] (number of positives/negatives according to ground-truth in the possible set)
    //
    //====================================================================================================================
    template<class T> T** sumScatterMatrix(T** scatterMatrix, unsigned long int dim, unsigned long int length, unsigned long int interestingDim, unsigned long int lastSpecialRows = 0, bool redundantClustersAreErrors = true) {
      T **m = NULL;
      unsigned long int x, y;

      //initialize detection scatter matrix with maybe additional rows
      MArray_2D(m, (long int)(2+lastSpecialRows), 2, T, "SC_MatrixFunctions.sumScatterMatrix: m");
      for (y = 0; y < 2+lastSpecialRows; y++) {
        for (x = 0; x < 2; x++) {
          m[y][x] = 0;
        }
      }

      //TP are just the correctly classified CEs
      m[0][0] = scatterMatrix[interestingDim][interestingDim];
      if (length-lastSpecialRows > dim && redundantClustersAreErrors == false) { //care for the redundant clusters, if any & wished
        for (y = dim; y < length-lastSpecialRows; y++) { //search only in the redundant clusters for rows containing additional material of the class corresponding to the interstingDim
          if (maxIdx(scatterMatrix[y], dim) == interestingDim) { //this additional material is found in rows among the redundant clusters where the cell with the greatest content has the index interestingDim
            m[0][0] += scatterMatrix[y][interestingDim];
          }
        }
      }
      
      //FP are all CEs in the same row as the interestingDim except in that special cell
      for (x = 0; x < dim ; x++) {
        if (x != interestingDim) {
          m[0][1] += scatterMatrix[interestingDim][x];
        }
      }
      if (length-lastSpecialRows > dim && redundantClustersAreErrors == false) { //care for the redundant clusters, if any & wished
        for (y = dim; y < length-lastSpecialRows; y++) { //search only in the redundant clusters for rows containing additional material of the class corresponding to the interstingDim
          if (maxIdx(scatterMatrix[y], dim) == interestingDim) { //this additional material is found in rows among the redundant clusters where the cell with the greatest content has the index interestingDim
            for (x = 0; x < dim ; x++) {
              if (x != interestingDim) {
                m[0][1] += scatterMatrix[y][x];
              }
            }
          }
        }
      }

      //FN are all CEs in the same col as the interestingDim except in that special cell, and, if wished, except the redundant clusters
      if (length-lastSpecialRows > dim && redundantClustersAreErrors == false) { //there are redundant clusters and they should not be traeted as errors
        for (y = 0; y < length-lastSpecialRows; y++) {
          if (y != interestingDim && maxIdx(scatterMatrix[y], dim) != interestingDim) { //the first condition is also included in the second one; it is stated for clearness: the correct cluster & the redundant clusters should not be counted
            m[1][0] += scatterMatrix[y][interestingDim];
          }
        }
      } else {
        for (y = 0; y < length-lastSpecialRows; y++) {
          if (y != interestingDim) {
            m[1][0] += scatterMatrix[y][interestingDim];
          }
        }
      }

      //TN are all other, till now not regarded cells without the optional lastSpecialRows
      if (length-lastSpecialRows > dim && redundantClustersAreErrors == false) { //there are redundant clusters and they should not be traeted as errors
        for (y = 0; y < length-lastSpecialRows; y++) {
          for (x = 0; x < dim; x++) {
            if (y != interestingDim && x != maxIdx(scatterMatrix[y], dim)) {
              m[1][1] += scatterMatrix[y][x];
            }
          }
        }
      } else {
        for (y = 0; y < length-lastSpecialRows; y++) {
          for (x = 0; x < dim; x++) {
            if (y != interestingDim && x != interestingDim) {
              m[1][1] += scatterMatrix[y][x];
            }
          }
        }
      }
  
      //care for special rows (if any):
      if (lastSpecialRows > 0) {
        for (y = length-lastSpecialRows; y < length; y++) {
          m[2+y-length+lastSpecialRows][0] = scatterMatrix[y][interestingDim];
          for (x = 0; x < dim; x++) {
            if (x != interestingDim) {
              m[2+y-length+lastSpecialRows][1] += scatterMatrix[y][x];
            }
          }
        }
      }

      return m;
    }

    //====================================================================================================================
    // returns a submatrix containing only of those columns marked with the "selected" character in the col- and row-
		// selection array (which corresponds with the corresponding matrix dimensions). finalCols and finalRows will be 
		// filled with those values on correct execution.
    //====================================================================================================================
		template<class T> T** getSubMatrix(T **matrix, unsigned long int dim, unsigned long int len, char *colSelection, char *rowSelection, unsigned long int &finalCols, unsigned long int &finalRows, char selected = 'x') {
			T **subMatrix = NULL;
			unsigned long int x, y, i, j = 0;

			//get dimensionality of sub matrix
			finalRows = 0;
			finalCols = 0;
			for (x = 0; x < dim; x++) {
				if (colSelection[x] == selected) {
					finalCols++;
				}
			}
			for (y = 0; y < len; y++) {
				if (rowSelection[y] == selected) {
					finalRows++;
				}
			}

			//create & fill submatrix
			if (finalCols > 0 && finalRows > 0) {
				MArray_2D(subMatrix, (long int)(finalRows), (long int)(finalCols), double, "SC_MatrixFunsions.getSubMatrix: subMatrix");
				for (y = 0; y < len; y++) {
					if (rowSelection[y] == selected) {
						i = 0;
						for (x = 0; x < dim; x++) {
							if (colSelection[x] == selected) {
								subMatrix[j][i++] = matrix[y][x];
							}
						}
						j++;
					}
				}
			} else {
				finalRows = 0;
				finalCols = 0;
			}

			return subMatrix;
		}

    //====================================================================================================================
		// replace zeros in the matrix with the replacement-value (whatever zero and the replacement is in the datatype of
		// this matrix); if the random==true, a random number with 0 mean and "replacement" variance is taken to replace zeros
		// the number of replacements is returned
    //====================================================================================================================
		template<class T> unsigned long int replaceZeros(T **matrix, unsigned long int len, unsigned long int dim, bool random = false, double replacement = 0.00001) {
			unsigned long int zeroCounter = 0;
			
			for (unsigned long int l = 0; l < len; l++) {
				for (unsigned long int d = 0; d < dim; d++) {
					if (matrix[l][d] == (T)(0.0)) {
						if (random == true) {
							matrix[l][d] = (T)(getRandom(0.0, replacement)); //N(0,replacement) distributed random numbers
						} else {
							matrix[l][d] = (T)(replacement);
						}
						zeroCounter++;
					}
				}
			}

			return zeroCounter;
		}

    //====================================================================================================================
		// adds random N(0, variance)-distributed "noise" to all cells of the given matrix
    //====================================================================================================================
		template<class T> void addNoise(T **matrix, unsigned long int len, unsigned long int dim, double variance = 1.0) {
			for (unsigned long int l = 0; l < len; l++) {
				for (unsigned long int d = 0; d < dim; d++) {
					matrix[l][d] = (T)(matrix[l][d] + getRandom(0.0, variance)); //N(0,variance) distributed random numbers
				}
			}

			return;
		}

		//-------------------------------------------------------------------------------------------------------------------
		//	differentiate the given array: array[x] = array[x] - weight*array[x-1]
		//-------------------------------------------------------------------------------------------------------------------
		template<class T> void differentiate(T* array, unsigned long int dim, double weight = 1.0) {
			for (unsigned long int d = dim-1; d > 0; d--) {
				array[d] -= (T)(weight*array[d-1]);
			}

			return;
		}

		//-------------------------------------------------------------------------------------------------------------------
		//	differentiate the given Matrix: matrix[t][d] = matrix[t][d] - weight*matrix[t-1][d]
		//	if onlyThisCol>=0, only the selected column of the matrix is differentiated 
		//-------------------------------------------------------------------------------------------------------------------
		template<class T> void differentiate(T** matrix, unsigned long int len, unsigned long int dim, long int onlyThisCol = -1, double weight = 1.0) {
			unsigned long int dMin = ((onlyThisCol > 0) ? onlyThisCol : 0), dMax = ((onlyThisCol >= 0) ? onlyThisCol+1 : dim);

			for (unsigned long int t = len-1; t > 0; t--) {
				for (unsigned long int d = dMin; d < dMax; d++) {
					matrix[t][d] -= (T)(weight*matrix[t-1][d]);
				}
			}

			return;
		}
};

#endif
