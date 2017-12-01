/**************************************************************************/
/*	This class implements gaussian functions: normal distribution,        */
/*  errror-function, etc.                                                 */ 
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.11.2005																								*/
/**************************************************************************/

#include <math.h>
#include "SC_Gauss.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Gauss::SC_Gauss(double tableSpan, double tableStep) {
  this->minX = 0.0 - tableSpan/2.0;
  this->maxX = 0.0 + tableSpan/2.0;
  this->uBound = (unsigned long int)(tableSpan / tableStep);
  this->step = tableStep;

  initTables();
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Gauss::~SC_Gauss() {
  MFree_1D(this->pdf);
  MFree_1D(this->log_pdf);
  MFree_1D(this->cdf);
  MFree_1D(this->log_cdf);
}

//====================================================================================================================
//	initialize the lookup-tables
//====================================================================================================================
void SC_Gauss::initTables(void) {
  const double unitVarianceFactor = 0.39894228040143267793994605993438;
  const double logUnitVarianceFactor = -0.91893853320467274178032973640562;
  unsigned long int i = 0;
  double t;

  MArray_1D(this->pdf, this->uBound + 1, double, "SC_Gauss.initTables: pdf");
  MArray_1D(this->log_pdf, this->uBound + 1, double, "SC_Gauss.initTables: log_pdf");
  MArray_1D(this->cdf, this->uBound + 1, double, "SC_Gauss.initTables: cdf");
  MArray_1D(this->log_cdf, this->uBound + 1, double, "SC_Gauss.initTables: log_cdf");

  for (double x = this->minX; x <= this->maxX; x += this->step) {    
    t = (x * x) / -2.0;
    this->pdf[i] = unitVarianceFactor * sclib::sExp(t);
    this->log_pdf[i] = logUnitVarianceFactor + t;

    if (i == 0) {
      this->cdf[i] = this->pdf[i];
    } else {
      this->cdf[i] = this->cdf[i-1] + ((this->step * this->pdf[i-1]) + (this->step * fabs(this->pdf[i] - this->pdf[i-1]) / 2.0)); //previous_sum + rectangle + triangle
    }
    if (this->cdf[i] > 1.0) {
      this->cdf[i] = 1.0;
    }

    this->log_cdf[i] = sclib::sLog(this->cdf[i]);

    i++;
  }

  //set start- and end-values
  this->pdf[0] = this->pdf[this->uBound] = 0.0;
  this->log_pdf[0] = this->log_pdf[this->uBound] = -705.0;
  this->cdf[0] = 0.0;
  this->log_cdf[0] = -705.0;
  this->cdf[this->uBound] = 1.0;
  this->log_cdf[this->uBound] = 0.0;

  return;
}

//====================================================================================================================
// Function to evaluate cumulative gaussian density B(z), where:
// B(z) = Integral(b(x), dx, -INF, z) = 1/2*sclib::sg(variance)*erf(sqrt(2)*(z-mean)/2*variance)+1/2
//====================================================================================================================
double SC_Gauss::erf(double x, double mean, double sd) {
	double erf = this->F.erf((sclib::sqrt_2 * (x - mean)) / (2.0 * sd));
	
  return 0.5 * sclib::sg(sd) * erf + 0.5;
}

//====================================================================================================================
// Function to evaluate the cumulative distribution function approximately if the analytical function gives 0 due to
// numerical problems: Here, it is tried to represent the integral as the sum of the areas of the trapezoids between 
// the density-value of an x1 and the density-value of an x2 = x1-step:
//
// y1 ------------------>_________|_/
// y2 ---->_|___________/         | <------- triangle
// _______/ |                     |	<------- rectangle
// ---------|---------------------|---------
//         x2                    x1
//
//    /\ x														___
//   |																\			((area of the rectangle: (x1-x2)*y2)) + 
//   |  gauss(x, mean, variance) dx ~	/__		 (area of the triangle: (x1-x2)*(y1-y2)/2))
//   |																		 
// \/  -inf
//====================================================================================================================
double SC_Gauss::approxErf(double x, double xDensity, double mean, double variance, double step) {
	double x1, x2, x1Density, x2Density;
	double rectangle, triangle;
	double sum = 0.0, oldSum;
  unsigned long int count = 0;
	
	x2 = x;
	x2Density = xDensity;

	do {
    if (++count >= 100) { //TODO: pTweak?!?!
      break;
    }

		oldSum = sum;
		x1 = x2;
		x1Density	= x2Density;

		x2 = x1 - step;
		x2Density = gaussian(x2, mean , variance);

    if (sclib::isFinite(x2Density) == true) {
		  rectangle	= (x1 - x2) * x2Density;
		  triangle = ((x1 - x2) * (x1Density - x2Density)) / 2.0;

		  sum += rectangle + triangle;
    }
	} while (sum != oldSum);

	return sum;
}

//====================================================================================================================
// Approximate cumulative gaussian as above, but using the scaledVariance to compute the normal density with the 
// fastGaussian()-method
//====================================================================================================================
double SC_Gauss::fastApproxErf(double x, double xDensity, double mean, double variance, double scaledVariance, double step) {
	double x1, x2, x1Density, x2Density;
	double rectangle, triangle;
	double sum = 0.0, oldSum;
  unsigned long int count = 0;
	
	x2 = x;
	x2Density = xDensity;

	do {
    if (++count >= 100) { //TODO: pTweak!?!?
      break;
    }

		oldSum = sum;
		x1 = x2;
		x1Density	= x2Density;

		x2 = x1 - step;
		x2Density = fastGaussian(x2, mean , variance, scaledVariance);

    if (sclib::isFinite(x2Density) == true) {
		  rectangle	= (x1 - x2) * x2Density;
		  triangle = ((x1 - x2) * (x1Density - x2Density)) / 2.0;

		  sum += rectangle + triangle;
    }
	} while (sum != oldSum);

	return sum;
}

//====================================================================================================================
// Function to evaluate the log of a cumulative gaussian density B(z)
// The Formula for log(a+b)=log(1+exp(log(a)-log(b)))+log(b) comes originaly from kingsbury&rayner and is cited by 
// nadas, nahamoo & picheny in "speech recognition using noise-adaptive prototypes"
//====================================================================================================================
double SC_Gauss::logErf(double x, double mean, double sd) {
	double erf = this->F.erf((sclib::sqrt_2 * (x - mean)) / (2.0 * sd));
	
  return sclib::sLog( 0.5 * sclib::sg(sd) * erf + 0.5);
}

//====================================================================================================================
// Multivariate gaussian pdf (originally Wesley-style, but checked for mathematical correctness now)
//====================================================================================================================
double SC_Gauss::multivariateGaussian(float *x, double *mean, double *variance, double *sd, unsigned short int dim) {
  double p = 0.0;
  int i;

	for (i = 0; i < dim; i++) {
    //p += log(1.0 / (sd[i] * sclib::sqrt_2pi)) + (((double)(x[i]) - mean[i]) * ((double)(x[i]) - mean[i]) / (-2.0 * variance[i]));
		p += log(sclib::one_div_sqrt_2pi * sd[i]) + (((double)(x[i]) - mean[i]) * ((double)(x[i]) - mean[i]) / (-2.0 * variance[i]));
	}

  return exp(p);
}

//====================================================================================================================
// Multivariate gaussian pdf, logarithmized
//====================================================================================================================
double SC_Gauss::multivariateLogGaussian(float *x, double *mean, double *variance, double *scaledVariance, unsigned short int dim) {
  double p = 0.0;
  int i;

	for (i = 0; i < dim; i++) {
		p += scaledVariance[i] + (((double)(x[i]) - mean[i]) * ((double)(x[i]) - mean[i]) / (-2.0 * variance[i]));
	}

  return p;
}

//====================================================================================================================
// Multivariate gaussian cdf (originally Wesley-style, but checked for mathematical correctness now)
//====================================================================================================================
double SC_Gauss::multivariateErf(float *x, double *mean, double *sd, unsigned short int dim) {
	int i, ChangeSign;
	double Prod = 1.0, z, Rz;
	
	for (i = 0; i < dim; i ++) {
		ChangeSign = 0;
		z = ((double)x[i] - mean[i]) / sd[i];
		if (z < 0.0) {
			z *=-1.0;
			ChangeSign=1;
		}
	
		if( (z>=0.0)&&(z<0.05) )	Rz=0.500;
		else if( (z>=0.05)&&(z<0.10) )	Rz=0.520;
		else if( (z>=0.10)&&(z<0.15) )	Rz=0.540;
		else if( (z>=0.15)&&(z<0.20) )	Rz=0.560;
		else if( (z>=0.20)&&(z<0.25) )	Rz=0.579;
		else if( (z>=0.25)&&(z<0.30) )	Rz=0.599;
		else if( (z>=0.30)&&(z<0.35) )	Rz=0.618;
		else if( (z>=0.35)&&(z<0.40) )	Rz=0.637;
		else if( (z>=0.40)&&(z<0.45) )	Rz=0.655;
		else if( (z>=0.45)&&(z<0.50) )	Rz=0.674;
		else if( (z>=0.50)&&(z<0.55) )	Rz=0.691;
		else if( (z>=0.55)&&(z<0.60) )	Rz=0.709;
		else if( (z>=0.60)&&(z<0.65) )	Rz=0.726;
		else if( (z>=0.65)&&(z<0.70) )	Rz=0.742;
		else if( (z>=0.70)&&(z<0.75) )	Rz=0.758;
		else if( (z>=0.75)&&(z<0.80) )	Rz=0.773;
		else if( (z>=0.80)&&(z<0.85) )	Rz=0.788;
		else if( (z>=0.85)&&(z<0.90) )	Rz=0.802;
		else if( (z>=0.90)&&(z<0.95) )	Rz=0.816;
		else if( (z>=0.95)&&(z<1.00) )	Rz=0.829;
		else if( (z>=1.00)&&(z<1.05) )	Rz=0.841;
		else if( (z>=1.05)&&(z<1.10) )	Rz=0.853;
		else if( (z>=1.10)&&(z<1.15) )	Rz=0.864;
		else if( (z>=1.15)&&(z<1.20) )	Rz=0.875;
		else if( (z>=1.20)&&(z<1.25) )	Rz=0.885;
		else if( (z>=1.25)&&(z<1.282) )	Rz=0.894;
		else if( (z>=1.282)&&(z<1.30) )	Rz=0.900;
		else if( (z>=1.30)&&(z<1.35) )	Rz=0.903;
		else if( (z>=1.35)&&(z<1.40) )	Rz=0.911;
		else if( (z>=1.40)&&(z<1.45) )	Rz=0.919;
		else if( (z>=1.45)&&(z<1.50) )	Rz=0.926;
		else if( (z>=1.50)&&(z<1.55) )	Rz=0.933;
		else if( (z>=1.55)&&(z<1.60) )	Rz=0.939;
		else if( (z>=1.60)&&(z<1.645) )	Rz=0.945;
		else if( (z>=1.645)&&(z<1.65) )	Rz=0.950;
		else if( (z>=1.65)&&(z<1.70) )	Rz=0.951;
		else if( (z>=1.70)&&(z<1.75) )	Rz=0.955;
		else if( (z>=1.75)&&(z<1.80) )	Rz=0.960;
		else if( (z>=1.80)&&(z<1.85) )	Rz=0.964;
		else if( (z>=1.85)&&(z<1.90) )	Rz=0.968;
		else if( (z>=1.90)&&(z<1.95) )	Rz=0.971;
		else if( (z>=1.95)&&(z<2.00) )	Rz=0.974;
		else if( (z>=2.00)&&(z<2.05) )	Rz=0.977;
		else if( (z>=2.05)&&(z<2.10) )	Rz=0.980;
		else if( (z>=2.10)&&(z<2.15) )	Rz=0.982;
		else if( (z>=2.15)&&(z<2.20) )	Rz=0.984;
		else if( (z>=2.20)&&(z<2.25) )	Rz=0.986;
		else if( (z>=2.25)&&(z<2.30) )	Rz=0.988;
		else if( (z>=2.30)&&(z<2.35) )	Rz=0.989;
		else if( (z>=2.35)&&(z<2.40) )	Rz=0.991;
		else if( (z>=2.40)&&(z<2.45) )	Rz=0.992;
		else if( (z>=2.45)&&(z<2.50) )	Rz=0.993;
		else if( (z>=2.50)&&(z<2.55) )	Rz=0.994;
		else if( (z>=2.55)&&(z<2.60) )	Rz=0.995;
		else if( (z>=2.60)&&(z<2.65) )	Rz=0.995;
		else if( (z>=2.65)&&(z<2.70) )	Rz=0.995;
		else if( (z>=2.70)&&(z<2.75) )	Rz=0.997;
		else if( (z>=2.75)&&(z<2.80) )	Rz=0.997;
		else if( (z>=2.80)&&(z<2.85) )	Rz=0.997;
		else if( (z>=2.85)&&(z<2.90) )	Rz=0.998;
		else if( (z>=2.90)&&(z<2.95) )	Rz=0.998;
		else if( (z>=2.95)&&(z<3.00) )	Rz=0.998;
		else if( (z>=3.00)&&(z<4.00) )	Rz=0.999;
		else	Rz=1.0;

		if (ChangeSign == 1)	Rz=1.0-Rz;
		Prod *=Rz;
  }

	return Prod;
}

//====================================================================================================================
//  Wesley-style "Multivariate gaussian cdf": its the average of the individual dimension's cdfs!!!
//====================================================================================================================
double SC_Gauss::multivariateAverageErf(float *x, double *mean, double *sd, unsigned short int dim) {
	int i, ChangeSign;
	double Prod = 0.0, z, Rz;
	
	for (i = 0; i < dim; i ++) {
		ChangeSign = 0;
		z = ((double)x[i] - mean[i]) / sd[i];
		if (z < 0.0) {
			z *=-1.0;
			ChangeSign=1;
		}
	
		if( (z>=0.0)&&(z<0.05) )	Rz=0.500;
		else if( (z>=0.05)&&(z<0.10) )	Rz=0.520;
		else if( (z>=0.10)&&(z<0.15) )	Rz=0.540;
		else if( (z>=0.15)&&(z<0.20) )	Rz=0.560;
		else if( (z>=0.20)&&(z<0.25) )	Rz=0.579;
		else if( (z>=0.25)&&(z<0.30) )	Rz=0.599;
		else if( (z>=0.30)&&(z<0.35) )	Rz=0.618;
		else if( (z>=0.35)&&(z<0.40) )	Rz=0.637;
		else if( (z>=0.40)&&(z<0.45) )	Rz=0.655;
		else if( (z>=0.45)&&(z<0.50) )	Rz=0.674;
		else if( (z>=0.50)&&(z<0.55) )	Rz=0.691;
		else if( (z>=0.55)&&(z<0.60) )	Rz=0.709;
		else if( (z>=0.60)&&(z<0.65) )	Rz=0.726;
		else if( (z>=0.65)&&(z<0.70) )	Rz=0.742;
		else if( (z>=0.70)&&(z<0.75) )	Rz=0.758;
		else if( (z>=0.75)&&(z<0.80) )	Rz=0.773;
		else if( (z>=0.80)&&(z<0.85) )	Rz=0.788;
		else if( (z>=0.85)&&(z<0.90) )	Rz=0.802;
		else if( (z>=0.90)&&(z<0.95) )	Rz=0.816;
		else if( (z>=0.95)&&(z<1.00) )	Rz=0.829;
		else if( (z>=1.00)&&(z<1.05) )	Rz=0.841;
		else if( (z>=1.05)&&(z<1.10) )	Rz=0.853;
		else if( (z>=1.10)&&(z<1.15) )	Rz=0.864;
		else if( (z>=1.15)&&(z<1.20) )	Rz=0.875;
		else if( (z>=1.20)&&(z<1.25) )	Rz=0.885;
		else if( (z>=1.25)&&(z<1.282) )	Rz=0.894;
		else if( (z>=1.282)&&(z<1.30) )	Rz=0.900;
		else if( (z>=1.30)&&(z<1.35) )	Rz=0.903;
		else if( (z>=1.35)&&(z<1.40) )	Rz=0.911;
		else if( (z>=1.40)&&(z<1.45) )	Rz=0.919;
		else if( (z>=1.45)&&(z<1.50) )	Rz=0.926;
		else if( (z>=1.50)&&(z<1.55) )	Rz=0.933;
		else if( (z>=1.55)&&(z<1.60) )	Rz=0.939;
		else if( (z>=1.60)&&(z<1.645) )	Rz=0.945;
		else if( (z>=1.645)&&(z<1.65) )	Rz=0.950;
		else if( (z>=1.65)&&(z<1.70) )	Rz=0.951;
		else if( (z>=1.70)&&(z<1.75) )	Rz=0.955;
		else if( (z>=1.75)&&(z<1.80) )	Rz=0.960;
		else if( (z>=1.80)&&(z<1.85) )	Rz=0.964;
		else if( (z>=1.85)&&(z<1.90) )	Rz=0.968;
		else if( (z>=1.90)&&(z<1.95) )	Rz=0.971;
		else if( (z>=1.95)&&(z<2.00) )	Rz=0.974;
		else if( (z>=2.00)&&(z<2.05) )	Rz=0.977;
		else if( (z>=2.05)&&(z<2.10) )	Rz=0.980;
		else if( (z>=2.10)&&(z<2.15) )	Rz=0.982;
		else if( (z>=2.15)&&(z<2.20) )	Rz=0.984;
		else if( (z>=2.20)&&(z<2.25) )	Rz=0.986;
		else if( (z>=2.25)&&(z<2.30) )	Rz=0.988;
		else if( (z>=2.30)&&(z<2.35) )	Rz=0.989;
		else if( (z>=2.35)&&(z<2.40) )	Rz=0.991;
		else if( (z>=2.40)&&(z<2.45) )	Rz=0.992;
		else if( (z>=2.45)&&(z<2.50) )	Rz=0.993;
		else if( (z>=2.50)&&(z<2.55) )	Rz=0.994;
		else if( (z>=2.55)&&(z<2.60) )	Rz=0.995;
		else if( (z>=2.60)&&(z<2.65) )	Rz=0.995;
		else if( (z>=2.65)&&(z<2.70) )	Rz=0.995;
		else if( (z>=2.70)&&(z<2.75) )	Rz=0.997;
		else if( (z>=2.75)&&(z<2.80) )	Rz=0.997;
		else if( (z>=2.80)&&(z<2.85) )	Rz=0.997;
		else if( (z>=2.85)&&(z<2.90) )	Rz=0.998;
		else if( (z>=2.90)&&(z<2.95) )	Rz=0.998;
		else if( (z>=2.95)&&(z<3.00) )	Rz=0.998;
		else if( (z>=3.00)&&(z<4.00) )	Rz=0.999;
		else	Rz=1.0;

		if (ChangeSign == 1)	Rz=1.0-Rz;
		Prod +=Rz; //http://www.fenews.com/fen46/risk-reward/risk-reward.htm: normally and, not or => should be a product, not a sum!
  }
	Prod /=(double)dim; //strange, too!

	return Prod;
}

//====================================================================================================================
//  computes what we call the "scalec variance": the first term in the gaussian function; depending on where the 
//  scaled variance should be used later on, the logarithmized version should be choosen (for fastLogGaussian(), e.g.)
//====================================================================================================================
double* SC_Gauss::scaleVariance(double *variance, int dim, bool logarithmize) {
	double *scaledVariance = NULL;
	
	MArray_1D(scaledVariance, dim, double, "SC_Gauss.scaleVariance: scaledVariance");

	for (int i = 0; i < dim; i++) {
		if (logarithmize == true) {
			scaledVariance[i] = log(sclib::one_div_sqrt_2pi * sqrt(variance[i]));
		} else {
			scaledVariance[i] = sclib::one_div_sqrt_2pi * sqrt(variance[i]);
		}
	}

	return scaledVariance;
}

//====================================================================================================================
//  computes what we call the "scalec variance": the first term in the gaussian function
//====================================================================================================================
double* SC_Gauss::variance2sd(double *variance, int dim) {
	double *sd = NULL;
	
	MArray_1D(sd, dim, double, "SC_Gauss.scaleVariance: sd");

	for (int i = 0; i < dim; i++) {
		sd[i] = sqrt(variance[i]);
	}

	return sd;
}
