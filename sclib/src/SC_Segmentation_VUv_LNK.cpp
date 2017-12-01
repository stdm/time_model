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

#include "SC_Segmentation_VUv_LNK.h"
#include "SC_FeatureHandler.h"
#include <SV_Data.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_VUv_LNK::SC_Segmentation_VUv_LNK(SC_TweakableParameters* pTweak) : SC_Segmentation_VUv(pTweak) {

}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_VUv_LNK::~SC_Segmentation_VUv_LNK(){

}

//====================================================================================================================
// analyze the features of the speech-frames in a given audio-segment and mark them as voiced or unvoiced in the 
// ground truth
// pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
// the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
// the return value is an SV_Data container with pitch-frequency per frame or NULL if the concrete algorithm isn't 
// able to extract pitch
//====================================================================================================================
SV_Data* SC_Segmentation_VUv_LNK::markVoicedUnvoicedSpeech(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	SV_Data *pEnergyZCR;
	SC_FeatureHandler handler(this->pTweak, false);

	if (pFeatures[sclib::bitPosition(sclib::featureSTE)] == NULL || pFeatures[sclib::bitPosition(sclib::featureSTE)]->Hdr.ID != sclib::featureSTE ||
  		pFeatures[sclib::bitPosition(sclib::featureZCR)] == NULL || pFeatures[sclib::bitPosition(sclib::featureZCR)]->Hdr.ID != sclib::featureZCR ||
		  pFeatures[sclib::bitPosition(sclib::featureSTE)]->Hdr.frameSize != pFeatures[sclib::bitPosition(sclib::featureZCR)]->Hdr.frameSize ||
			pFeatures[sclib::bitPosition(sclib::featureSTE)]->Hdr.frameStep != pFeatures[sclib::bitPosition(sclib::featureZCR)]->Hdr.frameStep ||
			pFeatures[sclib::bitPosition(sclib::featureSTE)]->Row != pFeatures[sclib::bitPosition(sclib::featureZCR)]->Row) {
    REPORT_ERROR(SVLIB_BadArg, "The LNK voiced/unvoiced speech detector needs synchronized STE and ZCR features to operate on!");
  } else {
		pEnergyZCR = handler.combineFeatureVectors(pFeatures, handler.getFeatureCount(), sclib::featureSTE|sclib::featureZCR);
    silence2unvoiced(pGT, segmentStart, segmentEnd, pEnergyZCR);
		MFree_0D(pEnergyZCR);
  }

  return NULL;
}

//====================================================================================================================
// Idea: A simple method: The LNK-ASD-Agorithm reliablie removes both silence, pauses, and unvoiced frames. We now 
//       only have to examin the energy/zcr-distribution in the prevously labeled sclib::atSilence-regions: The high-
//			 energy/zcr-parts within them is probably unvoiced speech. So This unvoiced-"detector" is  pretty much the 
//       same as the ASD itself, just operating on it's results. And it doesn't remove the unvoiced speech (that has 
//       the ASD already done), it just labels it to be able later to distinguish between speech-pauses and unvoiced 
//       speech!
//====================================================================================================================
int SC_Segmentation_VUv_LNK::silence2unvoiced(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pEnergyZCR) {
  unsigned long int y;
  long int x, start, end;
	float unvoicedThreshold = 0, maximum, minimum;
	float* averageFeature = new float[pTweak->segmentationVUvLnk.energyQuantizationLevel]; //average energy per class
	unsigned int* quantiCount = new unsigned int[pTweak->segmentationVUvLnk.energyQuantizationLevel]; //nr. of frames per class	
  unsigned int thresholdClass;
	unsigned short int feature;
  double prob;
	SV_Data* pFrames;	//the pFrames-matrix has the following semantics:
                    //	pFrames->Mat[][0] => frame-nr (casted to float)
										//	pFrames->Mat[][1]	=> energy of this frame, as computed earlier
										//	pFrames->Mat[][2]	=> zcr of this frame, as computed earlier
										//	pFrames->Mat[][3] => nr. of class this frame is quantized to (casted to float), as computed in the progress of this algorithm

	//remove all previous v/uv labels
	pGT->setSegment(segmentStart, segmentEnd, sclib::atVoiced|sclib::atUnvoiced, false, sclib::noSpeaker, sclib::modeLabelRemove);
	
	//make all speech segments voiced because the unvoiced parts were prevously labeled as silence
	pGT->setSegmentIf(segmentStart, segmentEnd, sclib::atSpeech, false, sclib::noType, false, sclib::atVoiced);

	//maybe we want do run this algorihtm more than one time to remove some more unvoiced-like frames...
	for (int c = 0; c < 1; c++) {

		//copy alle frames of this scene in an array, sort it by energy...
		pFrames = pGT->copyFramesTogether(pEnergyZCR, segmentStart, segmentStart, segmentEnd, sclib::atPause|sclib::atSilence, sclib::atNoise|sclib::atUnvoiced, 2, 1, true);
		
		if (pFrames != NULL) { // maybe nothing was labeled as sclib::atSilence/sclib::atPause
			for (feature = 1; feature <= 2; feature++) { //do the search for high-ernergy- and high-zcr-regions
				maximum	= 0.0;
				minimum = std::numeric_limits<float>::max();

				//...copy alle frames of this scene in an array, sort it by current feature
				for (x = 0; x < pFrames->Row; x++) { //'<' because the sceneEnd belong to the scene, but due to windowed analysis we loose one frame
					if (pFrames->Mat[x][feature] > maximum) {
						maximum = pFrames->Mat[x][feature];
					}
					if (pFrames->Mat[x][feature] < minimum) {
						minimum = pFrames->Mat[x][feature];
					}
				}
        sclib::quickSort(pFrames->Mat, 0, pFrames->Row-1, 4, feature); //parameter 3 doesn't mean the length of pFrames, but it's last element!!!

				//quantize the frames into energyQuantizationLevel classes, compute average of each class
				for (x = 0; x < pTweak->segmentationVUvLnk.energyQuantizationLevel; x++) {
					averageFeature[x] = 0.0;
					quantiCount[x] = 0;
				}
				for (x = 0; x < pFrames->Row; x++) {
					for (y = 1; y <= pTweak->segmentationVUvLnk.energyQuantizationLevel; y++) {
						if ((pFrames->Mat[x][feature] >= minimum + (maximum-minimum)*(float)(y - 1)/(float)pTweak->segmentationVUvLnk.energyQuantizationLevel) &&
								(pFrames->Mat[x][feature] <	 minimum + (maximum-minimum)*(float)y/(float)pTweak->segmentationVUvLnk.energyQuantizationLevel)) {
							pFrames->Mat[x][3] = (float)y;
							averageFeature[y-1] += pFrames->Mat[x][feature];
							quantiCount[y-1]++;
						}
					}
				}
				for (x = 0; x < pTweak->segmentationVUvLnk.energyQuantizationLevel; x++) {
					averageFeature[x] /= (quantiCount[x] > 0) ? (float)(quantiCount[x]) : (float)(1.0);
				}
				
				//calculate the Threshold
				thresholdClass = otsuThreshold(pTweak->segmentationVUvLnk.energyQuantizationLevel, quantiCount, pFrames->Row);
				unvoicedThreshold = averageFeature[thresholdClass];
				if (unvoicedThreshold == 0.0) {
					for (x = thresholdClass-1; x >= 0; x--) {
						unvoicedThreshold = averageFeature[x];
						if (unvoicedThreshold > 0.0) {
							break;
						}
					}
					if (unvoicedThreshold == 0.0) {
						for (x = thresholdClass+1; x <= this->pTweak->segmentationVUvLnk.energyQuantizationLevel; x++) {
							unvoicedThreshold = averageFeature[x];
							if (unvoicedThreshold > 0.0) {
								break;
							}
						}
					}
				}

				//mark unvoiced frames
				for (x = 0; x < pFrames->Row; x++) {
          start = (unsigned long)floor(pFrames->Mat[x][0]);
          end = (unsigned long)floor(pFrames->Mat[x][0]) + pEnergyZCR->Hdr.frameSize - 1;
					if (pFrames->Mat[x][feature] > unvoicedThreshold) {
						pGT->setSegment(start, end, sclib::atUnvoiced|sclib::atSpeech|sclib::atNoisySpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
            prob = 1.0 - (unvoicedThreshold / pFrames->Mat[x][feature]); //pseudo-"probability" of this frame being unvoiced
          } else {
            prob = (pFrames->Mat[x][feature] / unvoicedThreshold); //pseudo-"probability" of this frame being unvoiced
          }
          pGT->setProbability(start, end, sclib::atUnvoiced, prob);
				}

			} //for feature...
		} //pFrames != NULL

		MFree_0D(pFrames);
	} //for c
	
	MFree_1D(averageFeature);
	MFree_1D(quantiCount);

  //mark all silence-segments as silence alone by re-setting the silence-tag with rules switched on (this is postprocessing for the lnk-silence-detector where it was omitted to make this algorithm possible)
  pGT->setSegmentIf(segmentStart, segmentEnd, sclib::atSilence, false, sclib::noType, false, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);

  //mark all speech-segments, which are not unvoiced or pause, as voiced
  pGT->setSegmentIf(segmentStart, segmentEnd, sclib::atSpeech, false, sclib::atUnvoiced|sclib::atPause, false, sclib::atVoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);

	return SVLIB_Ok;	
}

//==================================================================================================================== 
//	Method to find optimal threshold for separating classed data into 2 distinct classes by minimization of 
//  inter-class variance (originally for binarising grey-scale-images), according to the following paper:
//  'A threshold selection method from gray level', N.Otsu, 1979, IEEE Trans. Systems, Man and Cybernetics
//
//	Meaning of parameters:
//		- classCount: total count of classes
//    - itemsPerClass: vector of length "classCount", containing count of items per class
//    - itemCount: total count of items
//    - return value: class providing the optimal threshold
//==================================================================================================================== 
unsigned int SC_Segmentation_VUv_LNK::otsuThreshold(unsigned int classCount, unsigned int* itemsPerClass, unsigned long int itemCount) {
  unsigned int	k, threshold = classCount / 2; //if everything goes wrong, still a meaningfull value will be returned
  double totalMean	= 0.0, variance = 0.0, maxVariance = 0.0;
  double zerothCumuMoment = 0.0, firstCumuMoment = 0.0;

  for(k = 1; k <= classCount; k++) {
	  totalMean += k * (itemsPerClass[k-1] / (double)itemCount);
  }

  for(k = 1; k <= classCount; k++) {
	  zerothCumuMoment += itemsPerClass[k-1] / (double)itemCount;
	  firstCumuMoment += k * itemsPerClass[k-1] / (double)itemCount;

	  variance = (totalMean * zerothCumuMoment) - firstCumuMoment;
	  variance *= variance;
    variance /= zerothCumuMoment * (1 - zerothCumuMoment) + std::numeric_limits<double>::epsilon();

	  if(variance > maxVariance) {
		  maxVariance = variance;
		  threshold = k;
	  }
  }

  return threshold;
}

//====================================================================================================================
// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
// are handled therein)
//====================================================================================================================
unsigned long int SC_Segmentation_VUv_LNK::getUncertaintyRegionWidth(void) {
	unsigned long int width = 0;

	//factors due to frame-based analysis
	width += 2 * this->pTweak->featureSte.frameStep;

	return width;
}
