/**************************************************************************/
/*    Responsibility:																											*/
/*      - This class implements the standard change detector that sets a  */
/*        changepoint just at the beginning of each new segment of the    */
/*        given type (speech/noise)                                       */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.02.2009																								*/
/**************************************************************************/

#include "SC_Segmentation_Changes_Std.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_Changes_Std::SC_Segmentation_Changes_Std(SC_TweakableParameters* pTweak, int mode) : SC_Segmentation_Changes(pTweak, mode) {

}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_Changes_Std::~SC_Segmentation_Changes_Std() {

}

//====================================================================================================================
//	detect changes in the characteristics of the audio-stream
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
//====================================================================================================================
int SC_Segmentation_Changes_Std::detectChanges(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	if (this->mode == sclib::modeSpeakerChange) {
		return markSegmentStarts(pGT, segmentStart, segmentEnd, sclib::atSpeech, sclib::atSpeakerBoundary);
	} else {
		return markSegmentStarts(pGT, segmentStart, segmentEnd, sclib::atNoise, sclib::atNoiseBoundary);
	}
}

//====================================================================================================================
//	sets a change point at each segment start
//====================================================================================================================
int SC_Segmentation_Changes_Std::markSegmentStarts(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int audioType, unsigned long int boundaryType) {
	long int speechSegStart, speechSegEnd;
	int changes = 0;
	
	pGT->setSegment(segmentStart, segmentEnd, boundaryType, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized, false, true);
	for (unsigned long int y = segmentStart; y < segmentEnd; y++) {
		pGT->getNextSegment(y, speechSegStart, speechSegEnd, audioType, sclib::searchForward);
		if (speechSegStart != sclib::noSegment && speechSegEnd != sclib::noSegment && speechSegStart <= (long int)(segmentEnd)) {
			pGT->setSegment(speechSegStart, speechSegStart, boundaryType, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
			changes++;
			y = speechSegEnd;
		} else {
			break;
		}
	}

	return changes;
}
