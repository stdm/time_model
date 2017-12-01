/**************************************************************************/
/*	This class is a centroid represented by whole signature (yes, it's		*/
/*  recursive: the centroid is build up by smaller centroids, e.g. the 		*/
/*  result of a kMeans clustering of a set of images, each represented by */
/*  points)																																*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.11.2006																								*/
/**************************************************************************/

#include <limits.h>
#include "SC_Centroid_Signature.h"
#include "SC_DistanceMeasures.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Centroid_Signature::SC_Centroid_Signature(SC_TweakableParameters *pTweak, SC_Signature *pSignature, bool justLink) : SC_Centroid(pTweak, justLink) {
	this->centroidType = sclib::centroidSignature;
	this->signature = pSignature;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Centroid_Signature::~SC_Centroid_Signature() {
	if (justLink == false) {
		MFree_0D(this->signature);
	}
}

//====================================================================================================================
//  Computes the ground distance between this and the secondCentroid, thereby checks for compatibility of the 2nd
//  centroid
//====================================================================================================================
double SC_Centroid_Signature::getDistance(SC_Centroid *secondCentroid) {
	if (secondCentroid->getCentroidType() == this->centroidType) {
		return getDistance((SC_Centroid_Signature*)(secondCentroid));
	} else {
		REPORT_ERROR(SVLIB_BadArg, "SC_Centroid_Signature.getDistance: Incompatibel centroid types!");
    return numeric_limits<double>::max();
	}
}

//====================================================================================================================
//  Computes the ground distance between this and the secondCentroid
//====================================================================================================================
double SC_Centroid_Signature::getDistance(SC_Centroid_Signature *secondCentroid) {
	double res = 0.0;

	if (sclib::bitTest(this->pTweak->distanceMeasure.groundDistance, sclib::dmEMD) == true) {
		res = SC_DistanceMeasures::EMD(this->signature, secondCentroid->getSignature(), this->pTweak);
	} else if (sclib::bitTest(this->pTweak->distanceMeasure.groundDistance, sclib::dmBeigi) == true) {
		res = SC_DistanceMeasures::beigi(this->signature, secondCentroid->getSignature());
	} else {
		REPORT_ERROR(SVLIB_BadArg, "SC_Centroid_Signature.getDistance: Specified ground-distance-measure is unkonwn, at least for this type of centroid!");
    res = numeric_limits<double>::max();
  }

	return res;
}

//====================================================================================================================
// To make output in Christian Beecks' format kind of virtual
//====================================================================================================================
std::ostream& SC_Centroid_Signature::centroidOut(std::ostream& os) {
	//format:
	//class-string : #centroids #dim (weight: dim1 dim2 ... dimN ) ... #

	if (this->getSignature() != NULL) {
		os << "[ ";
		//os << *(this->getSignature());

		//do signature output here different than in SC_Signature: no class-string at the beginning, no comment and newline at the end
		if (this->getSignature()->getN() > 0) {
			os << this->getSignature()->getN() << " " << this->getSignature()->getCentroids()[0]->getDim();
			for (int i = 0; i < this->getSignature()->getN(); i++) {
				os << " ( " << this->getSignature()->getWeight(i) << ": " << *(this->getSignature()->getCentroid(i)) << ")";
			}
			os << " ";
		}

		os << " ] ";
	}

	return os;
}