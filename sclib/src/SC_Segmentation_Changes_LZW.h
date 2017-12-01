/**************************************************************************/
/*    Responsibility:																											*/
/*      - This class implements the (speaker) change detector introduced  */
/*        in 'Real-Time Unsupervised Speaker Change Detection', L.Lu,     */
/*        H.J.Zhang, 2002 (IEEE)                                          */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Changes_LZW_H__
#define __SC_Segmentation_Changes_LZW_H__

#include "SC_Segmentation_Changes.h"
#include "SC_MatrixFunctions.h"
#include "SC_DistanceMeasures.h"

class SCLIB_API SC_Segmentation_Changes_LZW : public SC_Segmentation_Changes {
  
  private :

  protected :

    SC_DistanceMeasures* pDist; //for divergence-shape calculations 
    SC_MatrixFunctions* pMatrixFunc; //for cov estimations

    //====================================================================================================================
    //	calculates the adaptive treshold according to the paper
    //	'Real-Time Unsupervised Speaker Change Detection', L.Lu, H.J.Zhang, 2002 (IEEE)
    //
    //	it takes into account the last min(maxN, currSeg) distances to compute the threshold, alpha is a amplifying factor
    //
    //	formula: Th = alpha * 0.5 * SUM(distances[n]),
    //					 SUM from n=0 to min(maxN, currSeg)
    //====================================================================================================================
    double adaptiveThreshold(double** distances, unsigned int column, unsigned int currSeg, unsigned int maxN);
    
    //====================================================================================================================
    //	(Speaker) Change Detection as described in the following paper:
    //	'Real-Time Unsupervised Speaker Change Detection', L.Lu, H.J.Zhang, 2002 (IEEE)
    //	Copies all speechFrames between start&end together and moves a window over them
    //  this is more as the original intention of lu/zhang and helps that the adaptiveThreshold works because of more 
    //  last distances
    //
    //	return-value: # of changes
    //====================================================================================================================
    unsigned int LuZhangWindow(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, unsigned long int type = sclib::atSpeech, unsigned long int typesNot = sclib::atPause|sclib::atUnvoiced, unsigned long int boundaryType = sclib::atSpeakerBoundary);

    //====================================================================================================================
    //	estimates the needed prior probabilities for the bayesian fusion engine from the labeled training data and prints
		//  it to a file
    //====================================================================================================================
		int trainLuZhangWindow(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, unsigned long int type = sclib::atSpeech, unsigned long int typesNot = sclib::atPause|sclib::atUnvoiced, unsigned long int boundaryType = sclib::atSpeakerBoundary);

  public :

    SC_Segmentation_Changes_LZW(SC_TweakableParameters* pTweak, int mode);
    virtual ~SC_Segmentation_Changes_LZW();

    //====================================================================================================================
    //	detect changes in the characteristics of the audio-stream
    //  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    //  the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
    //====================================================================================================================
    virtual int detectChanges(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
		// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
    // implementation
    // The features needed for training are meant to reside in an array of feature-containers as returned by 
    // SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
		//====================================================================================================================
    virtual int trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
		// are hadnled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(void);

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureMFCC|sclib::featureLSP|sclib::featurePitch;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void) {return NULL;}

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_Changes_LZW";};
};

#endif
