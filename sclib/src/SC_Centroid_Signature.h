/**************************************************************************/
/*	This class is a centroid represented by whole signature (yes, it's		*/
/*  recursive: the centroid is build up by smaller centroids, e.g. the 		*/
/*  result of a kMeans clustering of a set of images, each represented by */
/*  points)																																*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.11.2006																								*/
/**************************************************************************/

#ifndef __SC_Centroid_Signature_H__
#define __SC_Centroid_Signature_H__

#include "SC_Centroid.h"
#include "SC_Signature.h"

class SC_Centroid_Signature : public SC_Centroid {

	private :

	protected:

		SC_Signature *signature;

		//====================================================================================================================
	  // To make output in Christian Beecks' format kind of virtual
		//====================================================================================================================
		virtual std::ostream& centroidOut(std::ostream& os);

	public :

		//====================================================================================================================
	  //  justLink has no effect here, the signature is always just linked; it just controls the way the destructor works
		//====================================================================================================================
		SC_Centroid_Signature(SC_TweakableParameters *pTweak = NULL, SC_Signature *pSignature = NULL, bool justLink = true);
		virtual ~SC_Centroid_Signature();

		//====================================================================================================================
	  //  Getters
		//====================================================================================================================
		virtual SC_Signature* getSignature(void) {return this->signature;}
		virtual int getDim(void) {return (this->signature!=NULL && this->signature->getCentroids()!=NULL) ? this->signature->getN()*this->signature->getCentroids()[0]->getDim() : 0;}

		//====================================================================================================================
	  //  Setters
		//====================================================================================================================
		virtual void setSignature(SC_Signature *newSignature) {this->signature = newSignature; return;};

		//====================================================================================================================
	  //  Computes the ground distance between this and the secondCentroid
		//====================================================================================================================
		virtual double getDistance(SC_Centroid *secondCentroid);
		virtual double getDistance(SC_Centroid_Signature *secondCentroid);
};

#endif
