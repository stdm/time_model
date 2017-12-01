/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for silence detectors                        				*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Silence_H__
#define __SC_Segmentation_Silence_H__

#include "SC_GroundTruth.h"
#include "SC_TweakableParameters.h"

class SC_Segmentation_Silence {
	
  private:

  protected:
	
		//====================================================================================================================
		// There are some tweakable parameters in the SC_Lib library; they can be centraly managed in this class.
		//====================================================================================================================
    SC_TweakableParameters* pTweak;

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Segmentation_Silence(SC_TweakableParameters* pTweak);
		virtual ~SC_Segmentation_Silence();

    //====================================================================================================================
		// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
    // implementation
    // The features needed for training are meant to reside in an array of feature-containers as returned by 
    // SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
		//====================================================================================================================
    virtual int trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {return SVLIB_Ok;};

		//====================================================================================================================
		// analyze the features of the given audio-segment and mark them as silence/non-silence in the ground truth
    // pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    // the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
		//====================================================================================================================
    virtual int markSilence(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) = 0;

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
		// are handled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(void) = 0;

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) = 0;

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void) = 0;

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) = 0;
};

#endif
