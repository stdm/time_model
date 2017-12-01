/**************************************************************************/
/*	This class is a centroid representing a point d-dimensional space     */
/*  as used in SC_Model_Pareto or for images															*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.11.2006																								*/
/**************************************************************************/

#include <limits.h>
#include "SC_Centroid_Point.h"
#include "SC_DistanceMeasures.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Centroid_Point::SC_Centroid_Point(SC_TweakableParameters *pTweak, int dim, double *coordinate, bool justLink) : SC_Centroid(pTweak, justLink) {
	this->centroidType = sclib::centroidPoint;
	this->dim = dim;

	if (this->justLink == false) {
		if (this->dim > 0) {
			MArray_1D(this->coordinate, this->dim, double, "SC_Centroid_Point: coordinate");
			for (int i = 0; i < this->dim; i++) {
				this->coordinate[i] = coordinate[i];
			}
		} else {
			this->coordinate = NULL;
		}
	} else {
		this->coordinate = coordinate;
	}
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Centroid_Point::~SC_Centroid_Point() {
	killCoordinate();
}

//====================================================================================================================
//	Set a complete new mean vector (copy or link according to justLink); it is assumed that it's dimension has been 
//  already stored in this->dim
//====================================================================================================================
void SC_Centroid_Point::setCoordinate(double *newCoordinate) {
	if (this->justLink == true) {
		this->coordinate = newCoordinate;
	} else {
		MFree_1D(this->coordinate);
		MArray_1D(this->coordinate, this->dim, double, "SC_Centroid_Point.setMean: coordinate");

		for (int i = 0; i < this->dim; i++) {
			this->coordinate[i] = newCoordinate[i];
		}
	}

	return;
}

//====================================================================================================================
//	Set one dim of the mean-vector, if it is in the correct range
//====================================================================================================================
void SC_Centroid_Point::setCoordinate(int d, double newCoordinate) {
	if (d > 0 && d < this->dim) {
		this->coordinate[d] = newCoordinate; 
	}

	return;
}

//====================================================================================================================
//	Change the dimension of the mean/variance; if it really changed, kill old mean/variance
//====================================================================================================================
void SC_Centroid_Point::setDim(int newDim) {
	if (newDim > 0) {
		if (newDim != this->dim) {
			killCoordinate();
		}
		this->dim = newDim;
	}

	return;
}

//====================================================================================================================
//	Free the mean or just remove the reference to it according to justLink
//====================================================================================================================
void SC_Centroid_Point::killCoordinate(void) {
	if (this->justLink == false) {
		MFree_1D(this->coordinate);
	} else {
		this->coordinate = NULL;
	}
}

//====================================================================================================================
//  Computes the ground distance between this and the secondCentroid, thereby checks for compatibility of the 2nd
//  centroid
//====================================================================================================================
double SC_Centroid_Point::getDistance(SC_Centroid *secondCentroid) {
	if (secondCentroid->getCentroidType() == this->centroidType) {
		return getDistance((SC_Centroid_Point*)(secondCentroid));
	} else {
		REPORT_ERROR(SVLIB_BadArg, "SC_Centroid_Point.getDistance: Incompatibel centroid types!");
    return numeric_limits<double>::max();
	}
}

//====================================================================================================================
//  Computes the ground distance between this and the secondCentroid
//====================================================================================================================
double SC_Centroid_Point::getDistance(SC_Centroid_Point *secondCentroid) {
	double res = SC_DistanceMeasures::euclid(this->coordinate, secondCentroid->getCoordinate(), this->dim);

	return res;
}

//====================================================================================================================
// To make output in Christian Beecks' format kind of virtual
//====================================================================================================================
std::ostream& SC_Centroid_Point::centroidOut(std::ostream& os) {
	//format:
	//class-string : #centroids #dim (weight: dim1 dim2 ... dimN ) ... #

	if (this->getDim() > 0) {
		for (int i = 0; i < this->getDim(); i++) {
			os << this->getCoordinate(i) << " ";
		}
	}

	return os;
}