/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates the Adaptive Silence Detector algorithm explained  */
/*        in 'Content-Based Movie Analysis And Indexing Based On Audio-   */
/*        Visual Cues', Li, Narayanan, Kuo, 2002 (Draft)                  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Silence_LNK_H__
#define __SC_Segmentation_Silence_LNK_H__

#include "SC_Segmentation_Silence.h"

class SC_Segmentation_Silence_LNK : public SC_Segmentation_Silence {
	
  private:

  protected:
	
    //====================================================================================================================
		// Threshold finding method for ASD
		//====================================================================================================================
    unsigned int otsuThreshold(unsigned int classCount, unsigned int* itemsPerClass, unsigned long int itemCount);

		//====================================================================================================================
		// Adaptive silence detection by Li/Narayanan/Kuo (LNK)
    // uses energy (column 0) and zcr (column 1) per frame as the features
		//====================================================================================================================
		int ASD(SC_GroundTruth*	pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data* pEnergy);

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Segmentation_Silence_LNK(SC_TweakableParameters* pTweak);
		virtual ~SC_Segmentation_Silence_LNK();

		//====================================================================================================================
		// analyze the features of the given audio-segment and mark them as silence/non-silence in the ground truth
    // pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    // the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
		//====================================================================================================================
    virtual int markSilence(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

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
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureSTE;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void) {return NULL;}

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_Silence_LNK";};
};

#endif

