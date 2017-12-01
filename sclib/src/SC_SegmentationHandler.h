/**************************************************************************/
/*    Responsibility:																											*/
/*      - Provides simple acces to various audio-segmentation algorithms  */
/*        (silence detection, V/Uv-classification, change detection, ...) */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_SegmentationHandler_H__
#define __SC_SegmentationHandler_H__

#include "SC_Api.h"
#include "SC_GroundTruth.h"
#include "SC_TweakableParameters.h"

class SCLIB_API SC_SegmentationHandler {
	
  private:

  protected:
	
    //====================================================================================================================
		// There are some tweakable parameters in the SC_Lib library; they can be centraly managed in this class.
		//====================================================================================================================
    SC_TweakableParameters* pTweak;

		bool verbose;

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_SegmentationHandler(SC_TweakableParameters* pTweak, bool verbose = true);
		virtual ~SC_SegmentationHandler();

    //====================================================================================================================
    //	Change the parameter-container to use extracted features with different parameter-settings
    //====================================================================================================================
    void setTweak(SC_TweakableParameters *pTweak) {this->pTweak = pTweak; return;}

    //====================================================================================================================
    //	detect and mark silence frames
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants sclib::feature* as indices into the array)
    //====================================================================================================================
    int silenceDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm = sclib::algorithm_sd_LNK);

    //====================================================================================================================
    //	distinguish between different audio classes (non-silences, such as speech, music, noise, ...)
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants sclib::feature* as indices into the array)
    //====================================================================================================================
    int audioClassification(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm = sclib::algorithm_ac_LZL);

    //====================================================================================================================
    //	mark speech as voiced or unvoiced
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants sclib::feature* as indices into the array)
    //  the return value is an SV_Data container with pitch-frequency per frame or NULL if the concrete algorithm isn't 
    //  able to extract pitch
    //====================================================================================================================
    SV_Data* unvoicedDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm = sclib::algorithm_vud_LNK);

    //====================================================================================================================
    //	detect and makr changes in the identity of the current speaker
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants sclib::feature* as indices into the array)
    //====================================================================================================================
    int speakerChangeDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm = sclib::algorithm_cd_LZW);

    //====================================================================================================================
    //	detect and mark changes in the characteristics of the non-speech frames
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
    //====================================================================================================================
    int acousticChangeDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm = sclib::algorithm_cd_LZW);

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to the choosen algorithm are taken into account (factors regarding the ground-truth class
		// are handled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(unsigned short int algorithm);

		//====================================================================================================================
		// Returns an or-concatenated list of sclib::feature* constants for all feature-types used by the given algorithm or
		// by all algorithms selected in the tweakable parameters, if algorithm = sclib::algorithm_nothing; independently from
		// choosing (or not choosing) a specific algorithm, via upToStage and startFromStage parameters features for 
		// algorithms can be selected according to the following stage succession model of audio segmentation:
		//   0 - silence detection
		//   1 - audio type classification
		//   2 - acoustic change detection
		//   3 - v/uv speech classification
		//   4 - speaker change detection
		// e.g., upToStage==2 means that the features-types needed by the choosen silence-detection-, atc- and acd-algorithms 
		// are returned.
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(unsigned short int algorithm = sclib::algorithm_nothing, short int upToStage = -1, short int startFromStage = -1);

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(unsigned short int algorithm = sclib::algorithm_nothing);
};

#endif
