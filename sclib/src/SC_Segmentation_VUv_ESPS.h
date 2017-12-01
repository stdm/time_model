/**************************************************************************/
/*    Responsibility:																											*/
/*      - implements a voiced/unvoiced detector based on the voicing      */
/*        decision of the ESPS pitch tracker decision. so, in fact,       */
/*        because the voicing decision corresponds with a non-zero pitch, */
/*        any pitch detector can be used that stores the pitch per frame  */
/*        in the first column of the feature matrix; this "algorithm"     */
/*        then reduces to storing this information in form of v/uv labels */
/*        in the ground truth object.                                     */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.04.2008																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_VUv_ESPS_H__
#define __SC_Segmentation_VUv_ESPS_H__

#include "SC_Segmentation_VUv.h"

class SC_Segmentation_VUv_ESPS : public SC_Segmentation_VUv {
	
  private:

  protected:
	
		//====================================================================================================================
    // Idea: Mark those frames with a non-zero pitch (in the first, i.e. 0th column of the feature-set) as voiced, the 
		//       others as unvoiced
		//====================================================================================================================
    int silence2unvoiced(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pPitch);

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Segmentation_VUv_ESPS(SC_TweakableParameters* pTweak);
		virtual ~SC_Segmentation_VUv_ESPS();
		
		//====================================================================================================================
		// analyze the features of the speech-frames in a given audio-segment and mark them as voiced or unvoiced in the 
    // ground truth
    // pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    // the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
    // the return value is an SV_Data container with pitch-frequency per frame or NULL if the concrete algorithm isn't 
    // able to extract pitch
		//====================================================================================================================
    virtual SV_Data* markVoicedUnvoicedSpeech(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
		// are handled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(void) {return 2*this->pTweak->segmentationVUvEsps.pitchParameters.frameStep;};

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featurePitch;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void);

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_VUv_ESPS";};
};

#endif
