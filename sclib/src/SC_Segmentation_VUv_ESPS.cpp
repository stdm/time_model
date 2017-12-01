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

#include "SC_Segmentation_VUv_ESPS.h"
#include "SC_FeatureHandler.h"
#include <SV_Data.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_VUv_ESPS::SC_Segmentation_VUv_ESPS(SC_TweakableParameters* pTweak) : SC_Segmentation_VUv(pTweak) {

}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_VUv_ESPS::~SC_Segmentation_VUv_ESPS(){

}

//====================================================================================================================
// analyze the features of the speech-frames in a given audio-segment and mark them as voiced or unvoiced in the 
// ground truth
// pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
// the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
// the return value is an SV_Data container with pitch-frequency per frame or NULL if the concrete algorithm isn't 
// able to extract pitch
//====================================================================================================================
SV_Data* SC_Segmentation_VUv_ESPS::markVoicedUnvoicedSpeech(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	if (pFeatures[sclib::bitPosition(sclib::featurePitch)] == NULL) {
    REPORT_ERROR(SVLIB_BadArg, "The ESPS voiced/unvoiced speech detector needs a pitch contour to operate on!");
  } else {
    silence2unvoiced(pGT, segmentStart, segmentEnd, pFeatures[sclib::bitPosition(sclib::featurePitch)]);
  }

  return NULL;
}

//====================================================================================================================
// Idea: Mark those frames with a non-zero pitch (in the first, i.e. 0th column of the feature-set) as voiced, the 
//       others as unvoiced
//====================================================================================================================
int SC_Segmentation_VUv_ESPS::silence2unvoiced(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pPitch) {
	int x;
	unsigned long int start, end;

	//remove all previous v/uv labels
	pGT->setSegment(segmentStart, segmentEnd, sclib::atVoiced|sclib::atUnvoiced, false, sclib::noSpeaker, sclib::modeLabelRemove);
	
	//copy the information from the feature-set to the groun-truth object
	for (x = 0; x < pPitch->Row; x++) {
    start = segmentStart + x*pPitch->Hdr.frameStep;
    end = start + pPitch->Hdr.frameSize;
		if (pPitch->Mat[x][0] > 0.0f) {
			pGT->setSegmentIf(start, end, sclib::atSpeech, false, sclib::noType, false, sclib::atVoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
		} else {
			pGT->setSegmentIf(start, end, sclib::atSpeech, false, sclib::noType, false, sclib::atUnvoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
		}
	}

  //mark all speech-segments, which are not unvoiced or pause, as voiced
  pGT->setSegmentIf(segmentStart, segmentEnd, sclib::atSpeech, false, sclib::atUnvoiced|sclib::atPause, false, sclib::atVoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);

	return SVLIB_Ok;
}

//====================================================================================================================
// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
// parameters
//====================================================================================================================
SC_TweakableParameters::SC_FeaturePar* SC_Segmentation_VUv_ESPS::getSpecialFeatureParameters(void) {
	//properly link the needed parameter-sets
	this->pTweak->segmentationVUvEsps.pitchParameters.Next = NULL;

	return &(this->pTweak->segmentationVUvEsps.pitchParameters);
}
