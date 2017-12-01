/**************************************************************************/
/*	This class implements gaussian functions: normal distribution,        */
/*  errror-function, etc.                                                 */ 
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.11.2005																								*/
/**************************************************************************/

#ifndef __SC_Gauss_H__
#define __SC_Gauss_H__

#include "SC_Aux.h"
#include <GN_Func.h>

class SC_Gauss {

	private :

    double *pdf, *log_pdf;
    double *cdf, *log_cdf;
    
    double minX, maxX, step;
    unsigned long int uBound; //highest index into the arrays
		GN_Func F;

	protected:

    void initTables(void);

	public :

	  SC_Gauss(double tableSpan = 80.0, double tableStep = 0.001);
		virtual ~SC_Gauss();

    inline unsigned long int getIdx(double x) {double idx = (x - this->minX) / this->step; return (idx > this->uBound) ? this->uBound : (idx > 0) ? (unsigned long int)(idx) : 0;}

    inline double tabledGaussian(double x, double mean, double sd) {double scaledX = sclib::zTransform(x, mean, sd); return this->pdf[(this->getIdx)(scaledX)] / sd;}
    inline double tabledGaussian(double scaledX, double sd) {return this->pdf[(this->getIdx)(scaledX)] / sd;}
    inline double tabledGaussian(unsigned long int idx, double sd) {return this->pdf[idx] / sd;}
        
    inline double tabledLogGaussian(double x, double mean, double sd) {double scaledX = sclib::zTransform(x, mean, sd); return this->log_pdf[(this->getIdx)(scaledX)] - log(sd);}
    inline double tabledLogGaussian(double x, double mean, double sd, double logSd) {double scaledX = sclib::zTransform(x, mean, sd); return this->log_pdf[(this->getIdx)(scaledX)] - logSd;}
    inline double tabledLogGaussian(double scaledX, double logSd) {return this->log_pdf[(this->getIdx)(scaledX)] - logSd;}
    inline double tabledLogGaussian(unsigned long int idx, double logSd) {return this->log_pdf[idx] - logSd;}
   
    inline double tabledErf(double x, double mean, double sd) {double scaledX = sclib::zTransform(x, mean, sd); return this->cdf[(this->getIdx)(scaledX)];}
    inline double tabledErf(double scaledX) {return this->cdf[(this->getIdx)(scaledX)];}
    inline double tabledErf(unsigned long int idx) {return this->cdf[idx];}

    inline double tabledLogErf(double x, double mean, double sd) {double scaledX = sclib::zTransform(x, mean, sd); return this->log_cdf[(this->getIdx)(scaledX)];}
    inline double tabledLogErf(double scaledX) {return this->log_cdf[(this->getIdx)(scaledX)];}
    inline double tabledLogErf(unsigned long int idx) {return this->log_cdf[idx];}

    inline void tabledGaussianAndErf(unsigned long int idx, double sd, double &pdf, double &cdf) {pdf = this->pdf[idx] / sd; cdf = this->cdf[idx]; return;}
    inline void tabledGaussianAndLogErf(unsigned long int idx, double logSd, double &pdf, double &cdf) {pdf = this->pdf[idx] - logSd; cdf = this->log_cdf[idx]; return;}
    inline void tabledLogGaussianAndErf(unsigned long int idx, double logSd, double &pdf, double &cdf) {pdf = this->log_pdf[idx] - logSd; cdf = this->cdf[idx]; return;}
    inline void tabledLogGaussianAndLogErf(unsigned long int idx, double logSd, double &pdf, double &cdf) {pdf = this->log_pdf[idx] - logSd; cdf = this->log_cdf[idx]; return;}
    
    inline double gaussian(double x, double mean, double variance) {double t = x - mean, exponent = -0.5 * t * t / variance; return (1.0 / (sclib::sqrt_2pi * sqrt(variance))) * sclib::sExp(exponent);}
    inline double gaussian(double x, double mean, double variance, double sd) {double t = x - mean, exponent = -0.5 * t * t / variance; return (1.0 / (sclib::sqrt_2pi * sd)) * sclib::sExp(exponent);}

    inline double logGaussian(double x, double mean, double variance) {double t = x - mean; return log(1.0 / (sclib::sqrt_2pi * sqrt(variance))) + (-0.5 * t * t / variance);}
    inline double logGaussian(double x, double mean, double variance, double sd) {double t = x - mean; return log(1.0 / (sclib::sqrt_2pi * sd)) + (-0.5 * t * t / variance);}

    inline double fastGaussian(double x, double mean, double variance, double scaledVariance) {double t = x - mean; return scaledVariance * sclib::sExp(-0.5 * t * t / variance);}		
    inline double fastLogGaussian(double x, double mean, double variance, double scaledVariance) {double t = x - mean; return scaledVariance + (-0.5 * t * t / variance);}

    double erf(double x, double mean, double sd);
		double logErf(double x, double mean, double sd);
    double approxErf(double x, double xDensity, double mean, double variance, double step);
    double fastApproxErf(double x, double xDensity, double mean, double variance, double scaledVariance, double step);

    double multivariateGaussian(float *x, double *mean, double *variance, double *sd, unsigned short int dim);
		double multivariateLogGaussian(float *x, double *mean, double *variance, double *scaledVariance, unsigned short int dim);
    double multivariateErf(float *x, double *mean, double *sd, unsigned short int dim);
		double multivariateAverageErf(float *x, double *mean, double *sd, unsigned short int dim);

		double* scaleVariance(double *variance, int dim, bool logarithmize = false);
		double* variance2sd(double *variance, int dim);
};

#endif
