/**************************************************************************/
/*	This class is the base class for centroids: Sets of cluster-centroid  */
/*  that serve as input for e.g. the Eearth Mover's Distance; more on			*/
/*  signatures can be found in "The Earth Mover's Distance as a Metric    */
/*  for Image Retrieval", Rubner, Tomasi, Guibas (2000)										*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.11.2006																								*/
/**************************************************************************/

#ifndef __SC_Centroid_H__
#define __SC_Centroid_H__

#include <iostream>
#include <stddef.h>
#include "SC_TweakableParameters.h"

class SC_Centroid {

	private :

	protected:

		SC_TweakableParameters *pTweak;
		int centroidType;
		bool justLink; 

		//====================================================================================================================
	  // To make output in Christian Beecks' format kind of virtual
		//====================================================================================================================
		virtual std::ostream& centroidOut(std::ostream& os) = 0;

	public:

		SC_Centroid(SC_TweakableParameters *pTweak = NULL, bool justLink = true);
		virtual ~SC_Centroid();

		//====================================================================================================================
	  //  Getters
		//====================================================================================================================
		virtual int getCentroidType(void) {return this->centroidType;}
		virtual int getDim(void) = 0;

		//====================================================================================================================
	  //  Setters
		//====================================================================================================================
		virtual void setTweak(SC_TweakableParameters *pTweak) {this->pTweak = pTweak; return;}

		//====================================================================================================================
	  //  If the just linked memory should be freed by the destructor anyway, just flip the status...
		//====================================================================================================================
		virtual void setJustLink(bool newStatus) {this->justLink = newStatus; return;}

		//====================================================================================================================
	  //  Computes the ground distance between this and the secondCentroid
		//====================================================================================================================
		virtual double getDistance(SC_Centroid *secondCentroid) = 0;

		//====================================================================================================================
	  //  Output in Christian Beecks' format
		//====================================================================================================================
		friend std::ostream& operator<<(std::ostream& os, SC_Centroid& centroid);
};

#endif
