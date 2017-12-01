/**************************************************************************/
/*	This class is a centroid representing a point d-dimensional space     */
/*  as used in SC_Model_Pareto or for images															*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.11.2006																								*/
/**************************************************************************/

#ifndef __SC_Centroid_Point_H__
#define __SC_Centroid_Point_H__

#include "SC_Centroid.h"

class SC_Centroid_Point : public SC_Centroid {

	private :

	protected:

		int dim;
		double *coordinate;

		void killCoordinate(void);

		//====================================================================================================================
	  // To make output in Christian Beecks' format kind of virtual
		//====================================================================================================================
		virtual std::ostream& centroidOut(std::ostream& os);

	public :

		SC_Centroid_Point(SC_TweakableParameters *pTweak = NULL, int dim = 0, double *coordinate = NULL, bool justLink = true);
		virtual ~SC_Centroid_Point();

		//====================================================================================================================
	  //  Getters
		//====================================================================================================================
		virtual double* getCoordinate(void) {return this->coordinate;}
		virtual double getCoordinate(int d) {return (d >= 0 && d < this->dim) ? this->coordinate[d] : 0.0;}// renamed geCoordiate() to getCoordinate() By Bing 23.03.07
		virtual int getDim(void) {return this->dim;}

		//====================================================================================================================
	  //  Setters
		//====================================================================================================================
		virtual void setCoordinate(double *newCoordinate);
		virtual void setCoordinate(int d, double newCoordinate);
		virtual void setDim(int newDim);

		//====================================================================================================================
	  //  Computes the ground distance between this and the secondCentroid
		//====================================================================================================================
		virtual double getDistance(SC_Centroid *secondCentroid);
		virtual double getDistance(SC_Centroid_Point *secondCentroid);
};

#endif
