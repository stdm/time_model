/**************************************************************************/
/*	This class is a centroid representing a diagonal covariance gaussian  */
/*  as used in the SC_MixtureModel classes																*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.11.2006																								*/
/**************************************************************************/

#ifndef __SC_Centroid_Gaussian_H__
#define __SC_Centroid_Gaussian_H__

#include "SC_Centroid.h"

class SC_Centroid_Gaussian : public SC_Centroid {

	private :

	protected:

		int dim;
		double *mean;
		double *variance;

		void killMean(void);
		void killVariance(void);

		//====================================================================================================================
	  // To make output in Christian Beecks' format kind of virtual
		//====================================================================================================================
		virtual std::ostream& centroidOut(std::ostream& os);

	public :

		SC_Centroid_Gaussian(SC_TweakableParameters *pTweak = NULL, int dim = 0, double *mean = NULL, double *variance = NULL, bool justLink = true);
		virtual ~SC_Centroid_Gaussian();

		//====================================================================================================================
	  //  Getters
		//====================================================================================================================
		virtual double* getMean(void) {return this->mean;}
		virtual double getMean(int d) {return (d >= 0 && d < this->dim) ? this->mean[d] : 0.0;}
		virtual double* getVariance(void) {return this->variance;}
		virtual double getVariance(int d) {return (d >= 0 && d < this->dim) ? this->variance[d] : 0.0;}
		virtual int getDim(void) {return this->dim;}

		//====================================================================================================================
	  //  Setters
		//====================================================================================================================
		virtual void setMean(double *newMean);
		virtual void setMean(int d, double newMean);
		virtual void setVariance(double *variance);
		virtual void setVariance(int d, double newVariance);
		virtual void setDim(int newDim);

		//====================================================================================================================
	  //  Computes the ground distance between this and the secondCentroid
		//====================================================================================================================
		virtual double getDistance(SC_Centroid *secondCentroid);
		virtual double getDistance(SC_Centroid_Gaussian *secondCentroid);
};

#endif
