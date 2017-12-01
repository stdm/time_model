/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates the simple threshold-based silence detector        */
/*        suggested in "Content-based Audio Classification and            */
/*        Segmentation by Using Support Vector Machines", Lu/Zhang/Li     */
/*        2003                                                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.10.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Silence_LZL_H__
#define __SC_Segmentation_Silence_LZL_H__

#include "SC_Segmentation_Silence.h"

class SC_Segmentation_Silence_LZL : public SC_Segmentation_Silence {
	
  private:

  protected:
	
	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Segmentation_Silence_LZL(SC_TweakableParameters* pTweak);
		virtual ~SC_Segmentation_Silence_LZL();

		//====================================================================================================================
		// analyze the features of the given audio-segment and mark them as silence/non-silence in the ground truth
		// use energy and zcr per frame as the features
		// pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
		// the log of the feature-set constants sclib::feature* as indices into the array)
		//====================================================================================================================
    virtual int markSilence(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
		// analyzes the given features of pure silence and prints statistics about the features; suitable threshold values are
		// the maximae...
		//====================================================================================================================
    virtual int trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
		// are handled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(void);

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureSTE|sclib::featureZCR;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void);

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_Silence_LZL";};
};

#endif
