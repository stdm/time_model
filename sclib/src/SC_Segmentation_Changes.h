/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for change detectors (speaker/acoustic)       				*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Changes_H__
#define __SC_Segmentation_Changes_H__

#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Model_Pareto.h"
#include "SC_Api.h"

class SCLIB_API SC_Segmentation_Changes {
  
  private :

  protected :

    SC_TweakableParameters* pTweak; //for user-defined parameters
		int mode; //tells whether to detect speaker- or acoustic changes: sclib::modeSpeakerChange || sclib::modeAcousticChange

    //====================================================================================================================
    //	If any algorithm detected a change at this window according to the previous ones, this function looks where in 
    //  this window the real changepoint may lie (e.g. the speech in this window is connected, then the changepoint is at 
    //  the specified refPoint of the window; if this window consists of many speech-segments, the changepoint is 
		//  somewhere at the segement-borders... we set it at the beginning of the segment that contributes most frames to the 
		//  window
		//  input- and output-values are, as always, sample-based
    //====================================================================================================================
    unsigned long int refineChangePoint(SC_GroundTruth *pGT, unsigned long int windowStart, unsigned long int windowEnd, unsigned long int oldWindowEnd, unsigned long int type, int refPoint = sclib::refpointStart);

  public :

    SC_Segmentation_Changes(SC_TweakableParameters* pTweak, int mode);
    virtual ~SC_Segmentation_Changes();

    //====================================================================================================================
		// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
    // implementation
    // The features needed for training are meant to reside in an array of feature-containers as returned by 
    // SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
		//====================================================================================================================
    virtual int trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {return SVLIB_Ok;};

    //====================================================================================================================
    //	detect changes in the characteristics of the audio-stream
    //====================================================================================================================
    virtual int detectChanges(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) = 0;

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
		// Toggles speaker- or acoustic change detection mode: sclib::modeSpeakerChange or sclib::modeAcousticChange
		//====================================================================================================================
		virtual void setMode(int newMode) {this->mode = newMode; return;}

		//====================================================================================================================
		// Returns a non-parametric density-estimation of the length of speech/noise segments between speaker/acoustic changes
		// type has to be sclib::atSpeech or sclib::atNoise
		//====================================================================================================================
		SC_Model_Pareto* getSegmentDurationDistribution(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int type = sclib::atSpeech);

		//====================================================================================================================
		// Returns a non-parametric density-estimation of the detectability of the changepoints for the given type;
		// type has to be sclib::atSpeech or sclib::atNoise
		//====================================================================================================================
		SC_Model_Pareto* getDetectability(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int type = sclib::atSpeech);

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) = 0;
};

#endif
