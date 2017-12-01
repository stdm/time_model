/**************************************************************************/
/*    Inspired by:																												*/
/*      - SC_BaggenstossEM 																								*/
/*																																				*/
/*    Responsibility:																											*/
/*      - this is a variant of SC_BaggenstossEM, which uses only diagonal	*/
/*        covariances and is therefore somewhat optimized for not using   */
/*        so much matrix algebra                                  				*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 23.02.2005																								*/
/**************************************************************************/

#ifndef __SC_BaggenstossEMex_H__
#define __SC_BaggenstossEMex_H__

#include "SC_BaggenstossEM.h"

class SC_BaggenstossEMex : public SC_BaggenstossEM {

	private :

	protected:

    double* minVariance; //minimum variance in contrast to min. standard-deviation in basic class.
	  double** variance; //array of variance-vectors (diagonal covariance matrices)

    double fullLogGaussian(double *x, double *mean, double *variance, unsigned short dim);
    double fullLogGaussian(float *x, double *mean, double *variance, unsigned short dim);
    virtual double mixtureCloseness(double* mean1, double* mean2, double* variance1, double* variance2, unsigned short int dim);
    virtual double emStep(SV_Data *pData, bool bias = false);
    virtual unsigned short int deflate(double minWeightOne, double minWeightAll);
    virtual unsigned short int merge(double maxCloseness);
    virtual unsigned short int split(SV_Data *pData, double threshold = 1.0);

	public :

		SC_BaggenstossEMex(SC_TweakableParameters *pTweak = NULL, bool verbose = false);
		SC_BaggenstossEMex(const SC_BaggenstossEMex& pParent);
		virtual ~SC_BaggenstossEMex();

    virtual void init(SV_Data *pData, unsigned short int mixtureCount, double minVariance, bool randInit = true);
    virtual void init(SV_Data *pData, unsigned short int mixtureCount, double* minVariance, bool randInit = true);
    virtual double train(SV_Data *pData, unsigned int &actualIterations, unsigned long int maxIterations = 150, long int samplesPerMode = -1, bool bias = false, double maxCloseness = 1.0, bool addModes = true, double splitThresh = -1.0);
    virtual double test(SV_Data *pData);
    virtual double** getCovar(unsigned short int mixture); //computes and returnes the covariance-matrix out of its stored variance-vector
	  friend ostream& operator<<(ostream& os, SC_BaggenstossEMex& Data);

		virtual double* getMinVariance(void) {return this->minVariance;}
		virtual double** getVariance(void) {return this->variance;}
		virtual void setVariance(double **newVariance) {MFree_2D(this->variance); this->variance = newVariance;}

		virtual unsigned int write(fstream *file);
		virtual unsigned int read(fstream *file, SV_DataIO::SV_DatatypeSizes *fileSizes);
};

#endif

