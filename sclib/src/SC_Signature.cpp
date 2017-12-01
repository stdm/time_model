/**************************************************************************/
/*	This class implements signatures: Sets of cluster-centroids that			*/
/*  serve as input for e.g. the Eearth Mover's Distance; more on          */
/*  signatures can be found in "The Earth Mover's Distance as a Metric    */
/*  for Image Retrieval", Rubner, Tomasi, Guibas (2000)										*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.11.2006																								*/
/**************************************************************************/

#include "SC_Signature.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	Constructor
//  The centroids get just linked, but destructed by this class; for weights, it can be decided
//====================================================================================================================
SC_Signature::SC_Signature(SC_Centroid **centroids, double *weights, int n, bool justLinkWeights, double smallestUnnormalizedWeight) {
	this->n = n;
	this->justLinkWeights = justLinkWeights;
	this->centroids = centroids;
	this->smallestUnnormalizedWeight = smallestUnnormalizedWeight;

	if (this->justLinkWeights == false) {
		if (this->n > 0) {
			MArray_1D(this->weights, this->n, double, "SC_Signature: weigths");
			for (int i = 0; i < this->n; i++) {
				this->weights[i] = weights[i];
			}
		}
	} else {
		this->weights = weights;
	}
}

//====================================================================================================================
//	Constructor; constructor for equal weights
//  The centroids get just linked, but destructed by this class
//====================================================================================================================
SC_Signature::SC_Signature(SC_Centroid **centroids, int n) {
	this->n = n;
	this->justLinkWeights = false;
	this->centroids = centroids;
	this->smallestUnnormalizedWeight = 1.0;

	if (this->n > 0) {
		MArray_1D(this->weights, this->n, double, "SC_Signature: weigths");
		for (int i = 0; i < this->n; i++) {
			this->weights[i] = 1.0 / (double)(this->n);
		}
	} else {
		this->weights = NULL;
	}
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Signature::~SC_Signature() {
	for (int i = 0; i < this->n; i++) {
		MFree_0D(this->centroids[i]);
	}
	MFree_1D(this->centroids);
	if (this->justLinkWeights == false) {
		MFree_1D(this->weights);
	}
}

//====================================================================================================================
//	Set a complete new weights vector (copy or link according to justLinkWeights); it is assumed that it's dimension 
//  has been already stored in this->n
//====================================================================================================================
void SC_Signature::setWeight(double *newWeight) {
	if (this->justLinkWeights == true) {
		this->weights = newWeight;
	} else {
		MFree_1D(this->weights);
		MArray_1D(this->weights, this->n, double, "SC_Signature.setWeights: weights");

		for (int i = 0; i < this->n; i++) {
			this->weights[i] = newWeight[i];
		}
	}

	return;
}

//====================================================================================================================
//	Set one dim of the weights-vector, if it is in the correct range
//====================================================================================================================
void SC_Signature::setWeight(int d, double newWeight) {
	if (d > 0 && d < this->n) {
		this->weights[d] = newWeight; 
	}

	return;
}

//====================================================================================================================
//	Change the number of centroids; if it really changed, release the old members
//====================================================================================================================
void SC_Signature::setN(int newN) {
	if (newN > 0) {
		if (newN != this->n) {
			killWeights();
			this->centroids = NULL;
		}
		this->n = newN;
	}

	return;
}

//====================================================================================================================
//	Link (in any case, just link) to the new centroids
//====================================================================================================================
void SC_Signature::setCentroids(SC_Centroid **newCentroids) {
	this->centroids = centroids;

	return;
}

//====================================================================================================================
//	Free the weights or just remove the reference to it according to justLinkWeights
//====================================================================================================================
void SC_Signature::killWeights(void) {
	if (this->justLinkWeights == false) {
		MFree_1D(this->weights);
	} else {
		this->weights = NULL;
	}
}

//====================================================================================================================
//  If the just linked memory should be freed by the destructor anyway, just flip the status...
//====================================================================================================================
void SC_Signature::setJustLink(bool newStatus, bool alsoForCentroids) {
	this->justLinkWeights = newStatus;
	if (alsoForCentroids == true) {
		for (int i = 0; i < this->n; i++) {
			this->centroids[i]->setJustLink(newStatus);
		}
	}

	return;
}

//====================================================================================================================
//  Output in Christian Beecks' format
//====================================================================================================================
std::ostream& operator<<(std::ostream& os, SC_Signature& sig) {
	//format:
	//class-string : #centroids #dim (weight: dim1 dim2 ... dimN ) ... #

	if (sig.getN() > 0) {
		os << ": " << sig.getN() << " " << sig.getCentroids()[0]->getDim();
		for (int i = 0; i < sig.getN(); i++) {
			os << " ( " << sig.getWeight(i) << ": " << *(sig.getCentroid(i)) << ")";
		}
		os << " #" << endl;
	}

	return os;
}
