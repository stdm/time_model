/**************************************************************************/
/*    Responsibility:																											*/
/*      - This class implements the (speaker) change detector introduced  */
/*        in Kotti, Benetos, Kotropoulos, "Computationally Efficient and  */
/*        Robust BIC-Based Speaker Segmentation", 2008                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 29.08.2008																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_Changes_KBK_H__
#define __SC_Segmentation_Changes_KBK_H__

#include "SC_Segmentation_Changes.h"
#include "SC_MatrixFunctions.h"
#include "SC_DistanceMeasures.h"

class SCLIB_API SC_Segmentation_Changes_KBK : public SC_Segmentation_Changes {
  
  private :

		double tolerance; //range (in [s]) around a groundtruth change point in which a hypothesized change point must lie (+- tolerance/2) in order to be considered a correct detection

  protected :

		//====================================================================================================================
		//  The segmentation algorithm from Kotti et al.: Kotti, Benetos, Kotropoulos, "Computationally Efficient and Robust 
		//  BIC-Based Speaker Segmentation", IEEE Transactions on Audio, Speech, and Language Processing, 2008, 16, 920-933
		//  return-value: # of changes
		//  lambda is the data-dependant BIC penalty weighting factor
		//  r is the mean utterance duration in seconds
		//  resultionFactor means: every r/resolutionFacor seconds, a BIC test will be done
		//====================================================================================================================
    unsigned int bicSegmentationKotti(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, double lambda, double r, unsigned long int type = sclib::atSpeech, unsigned long int typesNot = sclib::atPause|sclib::atUnvoiced, unsigned long int boundaryType = sclib::atSpeakerBoundary);

		//====================================================================================================================
		//  The segmentation algorithm from Cettolo et al.: Cettolo, Vescovi, Rizzi, "Evaluation of BIC-based Algorithms for 
		//  Audio Segmentation", Computer Speech and Language, 2005, 19, 147-170
		//  return-value: # of changes
		//  lambda is the data-dependant BIC penalty weighting factor
		//  r is the mean utterance duration in seconds
		//  resultionFactor means: every r/resolutionFacor seconds, a BIC test will be done
		//====================================================================================================================
    unsigned int bicSegmentationCettolo(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, double lambda, double r, unsigned long int type = sclib::atSpeech, unsigned long int typesNot = sclib::atPause|sclib::atUnvoiced, unsigned long int boundaryType = sclib::atSpeakerBoundary);

		//====================================================================================================================
		// Computes all delta-BIC values between start (inclusive) and end (exclusive) and returns the instant of and the 
		// biggest BIC value; calls isTrueChangePoint() for debug output
		//====================================================================================================================
		double getDeltaBIC(SV_Data *pFeatures, int startInclusive, int endExclusive, int resolutionInFrames, int Nmargin, double lambda, int &bestPos, SC_GroundTruth *pGT = NULL, int toleranceInSamples = 0, SV_Data *pSampleNumbers = NULL, bool verbose = false, bool lockToTrueCp = false);

		//====================================================================================================================
		// Returns true if the candidate change point at sampleNr is a true cp; may also write debug output
		//====================================================================================================================
		bool isTrueChangePoint(SC_GroundTruth *pGT, unsigned long int sampleNr, unsigned long int toleranceInSamples, bool verbose = true, double bic = 0.0, double lX = 0.0, double lY = 0.0, double lZ = 0.0, double penalty = 0.0);

		//====================================================================================================================
		// Simple but correct evaluation script to test SC_Score_ChangeDetection
		//====================================================================================================================
		void simpleEvaluation(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pStrippedSampleNumbers, double &recall, double &precision, double &F, double &far, double &mdr);

  public :

		SC_Segmentation_Changes_KBK(SC_TweakableParameters* pTweak, double tolerance, int mode = sclib::modeSpeakerChange);
    virtual ~SC_Segmentation_Changes_KBK();

		//====================================================================================================================
		//  performs the algorithm with an oracle in the background that correctly decides on which occasion a change-point
		//  lies; then, the resulting ineuqalities for lambda are solved in a way that the F1 measure on the change-point
		//  class is optimized; the found lambda is returned; the tolerance is the range in [s] in which (+-tolerance/2) a 
		//  true changepoint must lie in order for the oracle to decide that a given BIC test instance should return a value 
		//  >0
		//====================================================================================================================
		double findLambdaKotti(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, double r, double tolerance = 1.0, unsigned long int type = sclib::atSpeech, unsigned long int typesNot = sclib::noType, unsigned long int boundaryType = sclib::atSpeakerBoundary);

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
		virtual unsigned long int getUncertaintyRegionWidth(void);

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureMFCC;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void);

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_Changes_KBK";};
};

#endif
