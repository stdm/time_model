/**************************************************************************/
/*	This class is a centroid representing a diagonal covariance gaussian  */
/*  as used in the SC_MixtureModel classes																*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.11.2006																								*/
/**************************************************************************/

#include <limits.h>
#include "SC_Centroid_Gaussian.h"
#include "SC_DistanceMeasures.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Centroid_Gaussian::SC_Centroid_Gaussian(SC_TweakableParameters *pTweak, int dim, double *mean, double *variance, bool justLink) : SC_Centroid(pTweak, justLink) {
	this->centroidType = sclib::centroidGaussian;
	this->dim = dim;

	if (this->justLink == false) {
		if (this->dim > 0) {
			MArray_1D(this->mean, this->dim, double, "SC_Centroid_Gaussian: mean");
			if (variance != NULL) {
				MArray_1D(this->variance, this->dim, double, "SC_Centroid_Gaussian: variance");
			}
			for (int i = 0; i < this->dim; i++) {
				this->mean[i] = mean[i];
				if (variance != NULL) {
					this->variance[i] = variance[i];
				}
			}
		} else {
			this->mean = NULL;
			this->variance = NULL;
		}
	} else {
		this->mean = mean;
		this->variance = variance;
	}
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Centroid_Gaussian::~SC_Centroid_Gaussian() {
	killMean();
	killVariance();
}

//====================================================================================================================
//	Set a complete new mean vector (copy or link according to justLink); it is assumed that it's dimension has been 
//  already stored in this->dim
//====================================================================================================================
void SC_Centroid_Gaussian::setMean(double *newMean) {
	if (this->justLink == true) {
		this->mean = newMean;
	} else {
		MFree_1D(this->mean);
		MArray_1D(this->mean, this->dim, double, "SC_Centroid_Gaussian.setMean: mean");

		for (int i = 0; i < this->dim; i++) {
			this->mean[i] = newMean[i];
		}
	}

	return;
}

//====================================================================================================================
//	Set one dim of the mean-vector, if it is in the correct range
//====================================================================================================================
void SC_Centroid_Gaussian::setMean(int d, double newMean) {
	if (d > 0 && d < this->dim) {
		this->mean[d] = newMean; 
	}

	return;
}

//====================================================================================================================
//	Set a complete new variance vector (copy or link according to justLink); it is assumed that it's dimension has 
//  been already stored in this->dim
//====================================================================================================================
void SC_Centroid_Gaussian::setVariance(double *newVariance) {
	if (this->justLink == true) {
		this->variance = newVariance;
	} else {
		MFree_1D(this->variance);
		MArray_1D(this->variance, this->dim, double, "SC_Centroid_Gaussian.setVariance: variance");

		for (int i = 0; i < this->dim; i++) {
			this->variance[i] = newVariance[i];
		}
	}

	return;
}

//====================================================================================================================
//	Set one dim of the variance-vector, if it is in the correct range
//====================================================================================================================
void SC_Centroid_Gaussian::setVariance(int d, double newVariance) {
	if (d > 0 && d < this->dim) {
		this->variance[d] = newVariance; 
	}

	return;
}

//====================================================================================================================
//	Change the dimension of the mean/variance; if it really changed, kill old mean/variance
//====================================================================================================================
void SC_Centroid_Gaussian::setDim(int newDim) {
	if (newDim > 0) {
		if (newDim != this->dim) {
			killMean();
			killVariance();
		}
		this->dim = newDim;
	}

	return;
}

//====================================================================================================================
//	Free the mean or just remove the reference to it according to justLink
//====================================================================================================================
void SC_Centroid_Gaussian::killMean(void) {
	if (this->justLink == false) {
		MFree_1D(this->mean);
	} else {
		this->mean = NULL;
	}
}

//====================================================================================================================
//	Free the variance or just remove the reference to it according to justLink
//====================================================================================================================
void SC_Centroid_Gaussian::killVariance(void) {
	if (this->justLink == false) {
		MFree_1D(this->variance);
	} else {
		this->variance = NULL;
	}
}

//====================================================================================================================
//  Computes the ground distance between this and the secondCentroid, thereby checks for compatibility of the 2nd
//  centroid
//====================================================================================================================
double SC_Centroid_Gaussian::getDistance(SC_Centroid *secondCentroid) {
	if (secondCentroid->getCentroidType() == this->centroidType) {
		return getDistance((SC_Centroid_Gaussian*)(secondCentroid));
	} else {
		REPORT_ERROR(SVLIB_BadArg, "SC_Centroid_Gaussian.getDistance: Incompatibel centroid types!");
    return numeric_limits<double>::max();
	}
}

//====================================================================================================================
//  Computes the ground distance between this and the secondCentroid
//====================================================================================================================
double SC_Centroid_Gaussian::getDistance(SC_Centroid_Gaussian *secondCentroid) {
	double res = 0.0;

	if (sclib::bitTest(this->pTweak->distanceMeasure.groundDistance, sclib::dmEuclid) == true) {
		res = SC_DistanceMeasures::euclid(this->mean, secondCentroid->getMean(), this->dim);
  } else  if (sclib::bitTest(this->pTweak->distanceMeasure.groundDistance, sclib::dmMahalanobis) == true) {
    res = SC_DistanceMeasures::mahalanobis(this->mean, secondCentroid->getMean(), this->variance, secondCentroid->getVariance(), this->dim);
	} else if (sclib::bitTest(this->pTweak->distanceMeasure.groundDistance, sclib::dmKullbackLeibler) == true) {
    res = SC_DistanceMeasures::kullbackLeibler(this->mean, secondCentroid->getMean(), this->variance, secondCentroid->getVariance(), this->dim);
	} else if (sclib::bitTest(this->pTweak->distanceMeasure.groundDistance, sclib::dmBhattacharyya) == true) {
    res = SC_DistanceMeasures::bhattacharyya(this->mean, secondCentroid->getMean(), this->variance, secondCentroid->getVariance(), this->dim);
	} else {
		REPORT_ERROR(SVLIB_BadArg, "SC_Centroid_Gaussian.getDistance: Specified ground-distance-measure is unkonwn, at least for this type of centroid!");
    res = numeric_limits<double>::max();
  }

	return res;
}

//====================================================================================================================
// To make output in Christian Beecks' format kind of virtual
//====================================================================================================================
std::ostream& SC_Centroid_Gaussian::centroidOut(std::ostream& os) {
	//format:
	//class-string : #centroids #dim (weight: dim1 dim2 ... dimN ) ... #

	if (this->getDim() > 0) {
		for (int i = 0; i < this->getDim(); i++) {
			os << this->getMean(i) << " ";
		}
		for (int i = 0; i < this->getDim(); i++) {
			os << this->getVariance(i) << " ";
		}
	}

	return os;
}