/**************************************************************************/
/*	This class implements signatures: Sets of cluster-centroids that			*/
/*  serve as input for e.g. the Eearth Mover's Distance; more on          */
/*  signatures can be found in "The Earth Mover's Distance as a Metric    */
/*  for Image Retrieval", Rubner, Tomasi, Guibas (2000)										*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.11.2006																								*/
/**************************************************************************/

#ifndef __SC_Signature_H__
#define __SC_Signature_H__

#include <ios>
#include "SC_Centroid.h"

class SC_Signature {

	private :

	protected:

		int n;
		SC_Centroid **centroids; //an array of pointers to centoids
	  double *weights;
		double smallestUnnormalizedWeight; //if it is necessary for further processing steps to know not only how the centroids related to each other (via the weights), but also their absolute scaling, the smalles unnormalized weight can be given here.
		bool justLinkWeights;

		void killWeights(void);

	public :

		//====================================================================================================================
	  //  The controids get just linked, but destructed by this class; for weights, it can be decided
		//====================================================================================================================
	  SC_Signature(SC_Centroid **centroids = NULL, double *weights = NULL, int n = 0, bool justLinkWeights = true, double smallestUnnormalizedWeight = 1.0);
		SC_Signature(SC_Centroid **centroids, int n);
		virtual ~SC_Signature();

		//====================================================================================================================
	  //  Getters
		//====================================================================================================================
		virtual double* getWeight(void) {return this->weights;}
		virtual double getWeight(int d) {return (d < this->n) ? this->weights[d] : 0;}
		virtual int getN(void) {return this->n;}
		virtual SC_Centroid** getCentroids(void) {return this->centroids;}
		virtual SC_Centroid* getCentroid(int d) {return (d >= 0 && d < this->n) ? this->centroids[d] : NULL;}
		virtual double getSmallestUnnormalizedWeight(void) {return this->smallestUnnormalizedWeight;}

		//====================================================================================================================
	  //  Setters
		//====================================================================================================================
		virtual void setWeight(double *newWeight);
		virtual void setWeight(int d, double newWeight);
		virtual void setN(int newN);
		virtual void setCentroids(SC_Centroid **newCentroids);
		virtual void setSmallestUnnormalizedWeight(double newSmallstUnnormalizedWeight) {this->smallestUnnormalizedWeight = newSmallstUnnormalizedWeight; return;};

		//====================================================================================================================
	  //  If the just linked memory should be freed by the destructor anyway, just flip the status...
		//====================================================================================================================
		virtual void setJustLink(bool newStatus, bool alsoForCentroids = false);

		//====================================================================================================================
	  //  Output in Christian Beecks' format
		//====================================================================================================================
		friend std::ostream& operator<<(std::ostream& os, SC_Signature& sig);
};

#endif
