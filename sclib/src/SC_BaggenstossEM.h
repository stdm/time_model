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

#ifndef __SC_BaggenstossEM_H__
#define __SC_BaggenstossEM_H__

#include "SC_TweakableParameters.h"
#include "SC_MatrixFunctions.h"
#include <SV_Data.h>
#include <SV_DataIO.h>

class SC_BaggenstossEM {

	private :

	protected:

    bool verbose;                //controls whether full cholsky-covars are plotted via '<<' or just their determinants

    unsigned short int  maxMixtureCount;        //upper bound for mixture-splitting
	  unsigned short int	mixtureCount;						//count of mixture-components
		unsigned short int	dim;										//dimension of feature-vectors
    double*             minStd;                 //minimum standard-deviation
	  double*							weight;		              //mixture weight vector
	  double**						mean;			              //array of mean-vectors
	  double***						choleskyCovar;	        //array of cholesky-covariance-matrices; variance = std^2; std = sqrt(variance)
		unsigned long int		trainingDataCount;			//count of feature-vectors used for training
		double              ***U, ***V, ***S;       //cache for the results of Singular Value Decomposition of the cholesky-covariances; is filled on demand and set to zero if a choleskyCovar is altered
		double              ***transposedInvertedCholeskyCovar, *logDet; //cache for the transposed-inverted covar and the covar's log-determinant to speed up lqrEval(); is computed on demand and set to NULL whenever the cholesky-covar is changed

    SC_TweakableParameters* pTweak;
    SC_MatrixFunctions* pMatrixFunc;
    
    double lqrEval(double *x, double *mean, double **choleskyCovar, unsigned short dim);
    double lqrEval(float *x, double *mean, double **choleskyCovar, unsigned short dim);
		double lqrEval(float *x, double *mean, double **choleskyCovar, double** &transposedInvertedCholeskyCovar, unsigned short dim, double &logDet);
		double lqrEval(double *x, double *mean, double **choleskyCovar, double** &transposedInvertedCholeskyCovar, unsigned short dim, double &logDet);
    virtual double mixtureCloseness(double* mean1, double* mean2, double** choleskyCovar1, double** choleskyCovar2, unsigned short int dim, double** &S1, double** &S2, double** &V1, double** &V2, double** &U1, double** &U2, double** transposedInvertedCholesky1 = NULL, double** transposedInvertedCholesky2 = NULL, double logDet1 = 0.0, double logDet2 = 0.0);
    virtual double emStep(SV_Data *pData, bool bias = false);
    virtual unsigned short int deflate(double minWeightOne, double minWeightAll);
    virtual unsigned short int merge(double maxCloseness);
    virtual unsigned short int split(SV_Data *pData, double threshold = 1.0);

	public :

		SC_BaggenstossEM(SC_TweakableParameters *pTweak = NULL, bool verbose = false);
		SC_BaggenstossEM(const SC_BaggenstossEM& pParent);
		virtual ~SC_BaggenstossEM();

    virtual void init(SV_Data *pData, unsigned short int mixtureCount, double minStd, bool randInit = true);
    virtual void init(SV_Data *pData, unsigned short int mixtureCount, double* minStd, bool randInit = true);
    virtual double train(SV_Data *pData, unsigned int &actualIterations, unsigned long int maxIterations = 150, long int samplesPerMode = -1, bool bias = false, double maxCloseness = 1.0, bool addModes = true, double splitThresh = -1.0);
    virtual double test(SV_Data *pData);
    virtual double** getCovar(unsigned short int mixture); //computes and returnes the covariance-matrix out of its stored cholesky-decomposition
    virtual void setMaxMixtureCount(unsigned short int maxMixtureCount) {this->maxMixtureCount = maxMixtureCount;}
    virtual void setVerboseMode(bool verbose) {this->verbose = verbose;}
	  friend ostream& operator<<(ostream& os, SC_BaggenstossEM& Data);

		virtual unsigned short int getMaxMixtureCount(void) {return this->maxMixtureCount;}
		virtual unsigned short int getMixtureCount(void) {return this->mixtureCount;}
		virtual unsigned short int getDim(void) {return this->dim;}
		virtual double* getMinStd(void) {return this->minStd;}
		virtual double* getWeight(void) {return this->weight;}
		virtual double** getMean(void) {return this->mean;}
		virtual double*** getCholeskyCovar(void) {return this->choleskyCovar;}
		virtual SC_MatrixFunctions* getMatrixFunctions(void) {return this->pMatrixFunc;}
		virtual void setMean(double **newMean) {MFree_2D(this->mean); this->mean = newMean;}
		virtual void setWeight(double *newWeight) {MFree_1D(this->weight); this->weight = newWeight;}

		virtual unsigned int write(fstream *file);
		virtual unsigned int read(fstream *file, SV_DataIO::SV_DatatypeSizes *fileSizes);
};

#endif

