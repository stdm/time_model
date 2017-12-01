/**************************************************************************/
/*	This class is the base class for centroids: Sets of cluster-centroid */
/*  that serve as input for e.g. the Eearth Mover's Distance; more on			*/
/*  signatures can be found in "The Earth Mover's Distance as a Metric    */
/*  for Image Retrieval", Rubner, Tomasi, Guibas (2000)										*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.11.2006																								*/
/**************************************************************************/

#include "SC_Centroid.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Centroid::SC_Centroid(SC_TweakableParameters *pTweak, bool justLink) {
	this->pTweak = pTweak;
	this->justLink = justLink;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Centroid::~SC_Centroid() {

}

//====================================================================================================================
// To make output in Christian Beecks' format kind of virtual
//====================================================================================================================
std::ostream& operator<<(std::ostream& os, SC_Centroid& centroid) {
	return centroid.centroidOut(os);
}
