/**************************************************************************/
/*    Responsibility:																											*/
/*      - This class implements the standard chzange detector that sets a */
/*        changepoint just at the beginning of each new segment of the    */
/*        given type (speech/noise)                                       */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.02.2009																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Changes_Std_H__
#define __SC_Segmentation_Changes_Std_H__

#include "SC_Segmentation_Changes.h"

class SCLIB_API SC_Segmentation_Changes_Std : public SC_Segmentation_Changes {
  
  private :

  protected :

    //====================================================================================================================
    //	sets a change point at each segment start
    //====================================================================================================================
		int markSegmentStarts(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int audioType, unsigned long int boundaryType);

  public :

		SC_Segmentation_Changes_Std(SC_TweakableParameters* pTweak, int mode = sclib::modeSpeakerChange);
    virtual ~SC_Segmentation_Changes_Std();

    //====================================================================================================================
    //	detect changes in the characteristics of the audio-stream
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants sclib::feature* as indices into the array)
    //====================================================================================================================
    virtual int detectChanges(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
		// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
    // implementation
    // The features needed for training are meant to reside in an array of feature-containers as returned by 
    // SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
		//====================================================================================================================
    virtual int trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {return SVLIB_Ok;}

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
		// are hadnled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(void) {return 0;}

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureNoFeature;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void) {return NULL;}

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_Changes_Std";};
};

#endif
