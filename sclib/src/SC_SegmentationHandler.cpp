/**************************************************************************/
/*    Responsibility:																											*/
/*      - Provides simple acces to various audio-segmentation algorithms  */
/*        (silence detection, V/Uv-classification, change detection, ...) */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#include "SC_SegmentationHandler.h"
#include "SC_Segmentation_Silence_LNK.h"
#include "SC_Segmentation_Silence_LZL.h"
#include "SC_Segmentation_VUv_LNK.h"
#include "SC_Segmentation_Changes_LZW.h"
#include "SC_Segmentation_AudioType_LZL.h"
#include "SC_Segmentation_VUv_ESPS.h"
#include "SC_Segmentation_Changes_KBK.h"
#include "SC_Segmentation_Changes_Std.h"
#include "SC_FeatureHandler.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_SegmentationHandler::SC_SegmentationHandler(SC_TweakableParameters* pTweak, bool verbose) {
  this->pTweak = pTweak;
	this->verbose = verbose;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_SegmentationHandler::~SC_SegmentationHandler() {

}

//====================================================================================================================
//	detect and mark silence frames
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants sclib::feature* as indices into the array)
//====================================================================================================================
int SC_SegmentationHandler::silenceDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm) {
  int res = SVLIB_Fail;
  SC_Segmentation_Silence *pDetector = NULL;
	SC_FeatureHandler handler(this->pTweak, false);
  
  switch (algorithm) {
    case sclib::algorithm_sd_LNK: {
      pDetector = new SC_Segmentation_Silence_LNK(this->pTweak);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->markSilence(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_sd_LZL: {
      pDetector = new SC_Segmentation_Silence_LZL(this->pTweak);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->markSilence(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_nothing : {
      res = SVLIB_Ok;
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Silence detection algorithm unknown!");
      break;
    }
  }

  return res;
}

//====================================================================================================================
//	distinguish between different audio classes (non-silences, such as speech, music, noise, ...)
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants sclib::feature* as indices into the array)
//====================================================================================================================
int SC_SegmentationHandler::audioClassification(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm) {
  int res = SVLIB_Fail;
  SC_Segmentation_AudioType *pDetector = NULL;
	SC_FeatureHandler handler(this->pTweak, false);
  
  switch (algorithm) {
    case sclib::algorithm_ac_LZL: {
      pDetector = new SC_Segmentation_AudioType_LZL(this->pTweak);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->classifyAudioType(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_nothing : {
      res = SVLIB_Ok;
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Audio type classification algorithm unknown!");
      break;
    }
  }

  return res;
}

//====================================================================================================================
//	mark speech as voiced or unvoiced
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants sclib::feature* as indices into the array)
//  the return value is an SV_Data container with pitch-frequency per frame or NULL if the concrete algorithm isn't 
//  able to extract pitch
//====================================================================================================================
SV_Data* SC_SegmentationHandler::unvoicedDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm) {
  SV_Data* pRes = NULL;
  SC_Segmentation_VUv *pDetector = NULL;
	SC_FeatureHandler handler(this->pTweak, false);

  switch (algorithm) {
    case sclib::algorithm_vud_LNK: {
      if (this->pTweak->segmentationHandler.silenceDetectorMode != sclib::algorithm_sd_LNK) {
        REPORT_ERROR(SVLIB_BadArg, "The LNK-Voiced/Unvoiced detection algorithm makes only sense in combination with the LNK-Silence-Detector!");
      } else {
        pDetector = new SC_Segmentation_VUv_LNK(this->pTweak);
				handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
        pRes = pDetector->markVoicedUnvoicedSpeech(pGT, segmentStart, segmentEnd, pFeatures);
        MFree_0D(pDetector);
      }
      break;
    }
    case sclib::algorithm_vud_ESPS: {
      pDetector = new SC_Segmentation_VUv_ESPS(this->pTweak);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      pRes = pDetector->markVoicedUnvoicedSpeech(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_nothing : {
			//make all speech voiced if nothing is selected so that subsequent algorithms use maximum data available
			pGT->setSegmentIf(segmentStart, segmentEnd, sclib::atSpeech, false, sclib::noType, false, sclib::atVoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Voiced/Unvoiced detection algorithm unknown!");
      break;
    }
  }

  return pRes;
}

//====================================================================================================================
//	detect and makr changes in the identity of the current speaker
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants sclib::feature* as indices into the array)
//====================================================================================================================
int SC_SegmentationHandler::speakerChangeDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm) {
  int res = SVLIB_Fail;
  SC_Segmentation_Changes *pDetector = NULL;
	SC_FeatureHandler handler(this->pTweak, false);

  switch (algorithm) {
    case sclib::algorithm_cd_LZW: {
      pDetector = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeSpeakerChange);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_cd_KBK: {
      pDetector = new SC_Segmentation_Changes_KBK(this->pTweak, this->pTweak->segmentationChangesKbk.tolerance, sclib::modeSpeakerChange);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
		case sclib::algorithm_cd_Std: {
      pDetector = new SC_Segmentation_Changes_Std(this->pTweak, sclib::modeSpeakerChange);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
		}
    case sclib::algorithm_nothing : {
      //mark too short segments as too short (get's also done within the scd/asd-algorithms)
      pGT->markShortSpeechSegments(segmentStart, segmentEnd, sclib::max(1, pGT->getConverter()->ms2sample(pTweak->general.shortSpeechThreshold, sclib::alignmentEnd)));
			
			//set speaker boundaries at each speech segment's beginning if no groundtruth was loaded
			if (pGT->getGTtype() == sclib::gtMPEG7) {
				pDetector = new SC_Segmentation_Changes_Std(this->pTweak, sclib::modeSpeakerChange);
				handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
				res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
				MFree_0D(pDetector);
			} else {
				res = SVLIB_Ok;
			}

      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Speaker change detection algorithm unknown!");
      break;
    }
  }

  return res;
}

//====================================================================================================================
//	detect and mark changes in the characteristics of the non-speech frames
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants sclib::feature* as indices into the array)
//====================================================================================================================
int SC_SegmentationHandler::acousticChangeDetection(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned short int algorithm) {
  int res = SVLIB_Fail;
  SC_Segmentation_Changes *pDetector = NULL;
	SC_FeatureHandler handler(this->pTweak, false);

  switch (algorithm) {
    case sclib::algorithm_cd_LZW: {
      pDetector = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeAcousticChange);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_cd_KBK: {
      pDetector = new SC_Segmentation_Changes_KBK(this->pTweak, this->pTweak->segmentationChangesKbk.tolerance, sclib::modeAcousticChange);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
    }
		case sclib::algorithm_cd_Std: {
      pDetector = new SC_Segmentation_Changes_Std(this->pTweak, sclib::modeAcousticChange);
			handler.prepareFeatureSet(pFeatures, pDetector->getSpecialFeatureParameters());
      res = pDetector->detectChanges(pGT, segmentStart, segmentEnd, pFeatures);
      MFree_0D(pDetector);
      break;
		}
    case sclib::algorithm_nothing : {
      res = SVLIB_Ok;
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Acoustic change detection algorithm unknown!");
      break;
    }
  }

  return res;
}

//====================================================================================================================
// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
// Only the factors due to the choosen algorithm are taken into account (factors regarding thr ground-truth class
// are handled therein)
//====================================================================================================================
unsigned long int SC_SegmentationHandler::getUncertaintyRegionWidth(unsigned short int algorithm) {
  unsigned long int res;

	switch (algorithm) {
    case sclib::algorithm_nothing : {
      res = 0;
      break;
    }
    case sclib::algorithm_cd_LZW : {
      SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_cd_KBK : {
      SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_KBK(this->pTweak, this->pTweak->segmentationChangesKbk.tolerance, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_cd_Std : {
			SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_Std(this->pTweak, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
		case sclib::algorithm_sd_LNK : {
      SC_Segmentation_Silence *pDetector = new SC_Segmentation_Silence_LNK(this->pTweak);
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_sd_LZL : {
      SC_Segmentation_Silence *pDetector = new SC_Segmentation_Silence_LZL(this->pTweak);
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_ac_LZL : {
      SC_Segmentation_AudioType *pDetector = new SC_Segmentation_AudioType_LZL(this->pTweak);
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_vud_LNK : {
      SC_Segmentation_VUv *pDetector = new SC_Segmentation_VUv_LNK(this->pTweak);
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    case sclib::algorithm_vud_ESPS : {
      SC_Segmentation_VUv *pDetector = new SC_Segmentation_VUv_ESPS(this->pTweak);
			res = pDetector->getUncertaintyRegionWidth();
			MFree_0D(pDetector);
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Choosen algorithm unknown!");
      break;
    }
  }

  return res;
}

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
unsigned long int SC_SegmentationHandler::getUsedFeatures(unsigned short int algorithm, short int upToStage, short int startFromStage) {
  unsigned long int res = 0, individualRes;
	char *featureNames;
	SC_FeatureHandler extractor(this->pTweak, this->verbose);

	if (algorithm == sclib::algorithm_cd_LZW || 
		  (this->pTweak->segmentationHandler.changeDetectorMode == sclib::algorithm_cd_LZW && (algorithm==sclib::algorithm_nothing || upToStage>=2 || (startFromStage<=2 && startFromStage>=0)))) {
    SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_cd_KBK || 
		  (this->pTweak->segmentationHandler.changeDetectorMode == sclib::algorithm_cd_KBK && (algorithm==sclib::algorithm_nothing || upToStage>=2 || (startFromStage<=2 && startFromStage>=0)))) {
    SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_KBK(this->pTweak, this->pTweak->segmentationChangesKbk.tolerance, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_cd_Std || 
		 (this->pTweak->segmentationHandler.changeDetectorMode == sclib::algorithm_cd_Std && (algorithm==sclib::algorithm_nothing || upToStage>=2 || (startFromStage<=2 && startFromStage>=0)))) {
    SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_Std(this->pTweak, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_sd_LNK || 
		  (this->pTweak->segmentationHandler.silenceDetectorMode == sclib::algorithm_sd_LNK && (algorithm==sclib::algorithm_nothing || upToStage>=0 || (startFromStage<=0 && startFromStage>=0)))) {
    SC_Segmentation_Silence *pDetector = new SC_Segmentation_Silence_LNK(this->pTweak);
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_sd_LZL || 
		  (this->pTweak->segmentationHandler.silenceDetectorMode == sclib::algorithm_sd_LZL && (algorithm==sclib::algorithm_nothing || upToStage>=0 || (startFromStage<=0 && startFromStage>=0)))) {
    SC_Segmentation_Silence *pDetector = new SC_Segmentation_Silence_LZL(this->pTweak);
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_ac_LZL || 
		  (this->pTweak->segmentationHandler.audioTypeMode == sclib::algorithm_ac_LZL && (algorithm==sclib::algorithm_nothing || upToStage>=1 || (startFromStage<=1 && startFromStage>=0)))) {
    SC_Segmentation_AudioType *pDetector = new SC_Segmentation_AudioType_LZL(this->pTweak);
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_vud_LNK || 
		  (this->pTweak->segmentationHandler.vUvDetectorMode == sclib::algorithm_vud_LNK && (algorithm==sclib::algorithm_nothing || upToStage>=3 || (startFromStage<=3 && startFromStage>=0)))) {
    SC_Segmentation_VUv *pDetector = new SC_Segmentation_VUv_LNK(this->pTweak);
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_vud_ESPS || 
		  (this->pTweak->segmentationHandler.vUvDetectorMode == sclib::algorithm_vud_ESPS && (algorithm==sclib::algorithm_nothing || upToStage>=3 || (startFromStage<=3 && startFromStage>=0)))) {
    SC_Segmentation_VUv *pDetector = new SC_Segmentation_VUv_ESPS(this->pTweak);
		individualRes = pDetector->getUsedFeatures();
		res |= individualRes;
		if (this->verbose == true) {
			featureNames = extractor.getFeatureName(individualRes);
			printf("%s:%s;", pDetector->getName(), featureNames);	
			MFree_1D(featureNames);
		}
		MFree_0D(pDetector);
  }

  return res;
}

//====================================================================================================================
// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
// parameters
//====================================================================================================================
SC_TweakableParameters::SC_FeaturePar* SC_SegmentationHandler::getSpecialFeatureParameters(unsigned short int algorithm) {
  SC_TweakableParameters::SC_FeaturePar *pParameterList = NULL;

	if (algorithm == sclib::algorithm_cd_LZW || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.changeDetectorMode == sclib::algorithm_cd_LZW)) {
    SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_LZW(this->pTweak, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_cd_KBK || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.changeDetectorMode == sclib::algorithm_cd_KBK)) {
    SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_KBK(this->pTweak, this->pTweak->segmentationChangesKbk.tolerance, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_cd_Std || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.changeDetectorMode == sclib::algorithm_cd_Std)) {
    SC_Segmentation_Changes *pDetector = new SC_Segmentation_Changes_Std(this->pTweak, sclib::modeSpeakerChange); //the mode has no effect here, so pick sc-mode...
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_sd_LNK || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.silenceDetectorMode == sclib::algorithm_sd_LNK)) {
    SC_Segmentation_Silence *pDetector = new SC_Segmentation_Silence_LNK(this->pTweak);
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_sd_LZL || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.silenceDetectorMode == sclib::algorithm_sd_LZL)) {
    SC_Segmentation_Silence *pDetector = new SC_Segmentation_Silence_LZL(this->pTweak);
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_ac_LZL || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.audioTypeMode == sclib::algorithm_ac_LZL)) {
    SC_Segmentation_AudioType *pDetector = new SC_Segmentation_AudioType_LZL(this->pTweak);
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_vud_LNK || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.vUvDetectorMode == sclib::algorithm_vud_LNK)) {
    SC_Segmentation_VUv *pDetector = new SC_Segmentation_VUv_LNK(this->pTweak);
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }
	if (algorithm == sclib::algorithm_vud_ESPS || (algorithm == sclib::algorithm_nothing && this->pTweak->segmentationHandler.vUvDetectorMode == sclib::algorithm_vud_ESPS)) {
    SC_Segmentation_VUv *pDetector = new SC_Segmentation_VUv_ESPS(this->pTweak);
		sclib::addToList(pParameterList, pDetector->getSpecialFeatureParameters());
		MFree_0D(pDetector);
  }

  return pParameterList;
}
