/**************************************************************************/
/*    Responsibility:																											*/
/*      - implements a voiced/unvoiced detector based on the adaptive     */
/*        silence detector (ASD) described in 'Content-Based Movie        */
/*        Analysis And Indexing Based On Audio-Visual Cues', Li,          */
/*        Narayanan, Kuo, 2002 (Draft)                                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_VUv_LNK_H__
#define __SC_Segmentation_VUv_LNK_H__

#include "SC_Segmentation_VUv.h"

class SC_Segmentation_VUv_LNK : public SC_Segmentation_VUv {
	
  private:

  protected:
	
    //====================================================================================================================
		// Threshold finding method for ASD
		//====================================================================================================================
    unsigned int otsuThreshold(unsigned int classCount, unsigned int* itemsPerClass, unsigned long int itemCount);

		//====================================================================================================================
    // Idea: A simple method: The LNK-ASD-Agorithm reliablie removes both silence, pauses, and unvoiced frames. We now 
    //       only have to examin the energy/zcr-distribution in the prevously labeled sclib::atSilence-regions: The high-
		//			 energy/zcr-parts within them is probably unvoiced speech. So This unvoiced-"detector" is  pretty much the 
		//       same as the ASD itself, just operating on it's results. And it doesn't remove the unvoiced speech (that has 
		//       the ASD already done), it just labels it to be able later to distinguish between speech-pauses and unvoiced 
		//       speech!
		//====================================================================================================================
    int silence2unvoiced(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pFeatures);

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Segmentation_VUv_LNK(SC_TweakableParameters* pTweak);
		virtual ~SC_Segmentation_VUv_LNK();
		
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
		virtual unsigned long int getUncertaintyRegionWidth(void);

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureSTE|sclib::featureZCR;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void) {return NULL;}

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_VUv_LNK";};
};

#endif
