/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates the speech-signal-processing algorithms like				*/
/*				- silence detection																							*/
/*				- speaker chage detection																				*/
/*				- ...																														*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.03.2004																								*/
/**************************************************************************/

#include "SC_Aux.h"
#include "SC_Segmentation_Silence_LNK.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_Silence_LNK::SC_Segmentation_Silence_LNK(SC_TweakableParameters* pTweak) : SC_Segmentation_Silence(pTweak) {
  
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_Silence_LNK::~SC_Segmentation_Silence_LNK() {

}

//====================================================================================================================
// analyze the features of the given audio-segment and mark them as silence/non-silence in the ground truth
// use energy (column 0) and zcr (column 1) per frame as the features
// pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
// the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
//====================================================================================================================
int SC_Segmentation_Silence_LNK::markSilence(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
  if (pFeatures[sclib::bitPosition(sclib::featureSTE)] == NULL || pFeatures[sclib::bitPosition(sclib::featureSTE)]->Hdr.ID != sclib::featureSTE) {
    REPORT_ERROR(SVLIB_BadArg, "The LNK silence detector needs ShortTimeEnergy features to operate on!");
    return SVLIB_Fail;
  } else {
    return ASD(pGT, segmentStart, segmentEnd, pFeatures[sclib::bitPosition(sclib::featureSTE)]);
  }
}

//==================================================================================================================== 
//	Adaptive Silence Detector as described in the following paper:
//	'Content-Based Movie Analysis And Indexing Based On AudioVisual Cues', Li, Narayanan, Kuo, 2002 (Draft)
//
//	It works like this:
//		- assume that speech has a higher ernergy-level than background (noise, hopefully silence)
//		-	take all audioframes of a segment, sort them in an array by their frame-energy (desc.)
//		- the lowest value indicates the lower bound of noise, the highest one the upper bound of speech energy
//		- the threshold to seperate both classes from each other lies somewhere between them
//		- quantize the values in energyQuantizationLevel classes
//		- take the sum of the average energys of the first and last 3 classes and the signal-to-noise-ratio to 
//			compute threshold (the paper doesn't tell how exactly use theses components)
//      ATTENTION: here's a change to the proposed method: we take the threshold-function by Otsu (see below for
//      description)
//==================================================================================================================== 
int	SC_Segmentation_Silence_LNK::ASD(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data* pEnergy) {
  unsigned short int y;
	unsigned long int x, start, end, frameCount = pEnergy->Row; //frameCount = pGT->getAudioFrameCountInScene(segmentStart, segmentEnd, pEnergyZCR->Hdr.frameSize, pEnergyZCR->Hdr.frameStep); //TODO: test if this is equal to pEnergyZCR->Row (should be...)
  float maximum	= 0.0, minimum = std::numeric_limits<float>::max();
	float* averageEnergy = new float[pTweak->segmentationSilenceLnk.energyQuantizationLevel]; //average energy per class
	unsigned int thresholdClass, *quantiCount = new unsigned int[pTweak->segmentationSilenceLnk.energyQuantizationLevel]; //nr. of frames per class	
	float silenceThreshold = 0;
  double prob; 
	SV_Data* pFrames;	//the pFrames-matrix has the following semantics:
					  				//	pFrames->Mat[][0] => sample-nr of first sample in frame (casted to float)
						  			//	pFrames->Mat[][1]	=> energy of this frame, as computed earlier
					  	  		//	pFrames->Mat[][2] => nr. of class this frame is quantized to (casted to float)

	//remove all previous silence tags
	pGT->setSegment(segmentStart, segmentEnd, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd);

	//copy alle frames of this scene in an array, sort it by energy
	pFrames = new SV_Data(frameCount, 3);
	for (x = 0; x < frameCount; x++) { //'<' because the sceneEnd belong to the scene, but due to windowed analysis we loose one frame ?!?
		pFrames->Mat[x][0] = (float)(segmentStart + pGT->getConverter()->audioFrame2sample(x, pEnergy->Hdr.frameSize, pEnergy->Hdr.frameStep));
		pFrames->Mat[x][1] = pEnergy->Mat[x][0];
		if (pFrames->Mat[x][1] > maximum) {maximum = pFrames->Mat[x][1];}
		if (pFrames->Mat[x][1] < minimum) {minimum = pFrames->Mat[x][1];}
	}
  sclib::quickSort(pFrames->Mat, 0, frameCount-1, 3, 1); //parameter 3 doesn't mean the length of pFrames, but it's last element!!!

	//quantize the frames into pTweak->energyQuantizationLevel energy-classes, compute average of each class
	for (x = 0; x < pTweak->segmentationSilenceLnk.energyQuantizationLevel; x++) {
		averageEnergy[x] = 0.0;
		quantiCount[x] = 0;
	}
	for (x = 0; x < frameCount; x++) {
		for (y = 1; y <= pTweak->segmentationSilenceLnk.energyQuantizationLevel; y++) {
			if ((pFrames->Mat[x][1] >= minimum + (maximum-minimum)*(float)(y - 1)/(float)pTweak->segmentationSilenceLnk.energyQuantizationLevel) &&
					(pFrames->Mat[x][1] <	 minimum + (maximum-minimum)*(float)y/(float)pTweak->segmentationSilenceLnk.energyQuantizationLevel)) {
				pFrames->Mat[x][2] = (float)y;
				averageEnergy[y-1] += pFrames->Mat[x][1];
				quantiCount[y-1]++;
			}
		}
	}
	for (x = 0; x < pTweak->segmentationSilenceLnk.energyQuantizationLevel; x++) {
		averageEnergy[x] /= (quantiCount > 0) ? (float)(quantiCount[x]) : (float)(1.0);
	}
	
	//calculate the Threshold
  thresholdClass = otsuThreshold(pTweak->segmentationSilenceLnk.energyQuantizationLevel, quantiCount, frameCount);
  silenceThreshold = averageEnergy[thresholdClass];
	if (silenceThreshold == 0.0) {
		for (x = thresholdClass-1; x >= 0; x--) {
			silenceThreshold = averageEnergy[x];
			if (silenceThreshold > 0.0) {
				break;
			}
		}
		if (silenceThreshold == 0.0) {
			for (x = thresholdClass+1; x <= this->pTweak->segmentationSilenceLnk.energyQuantizationLevel; x++) {
				silenceThreshold = averageEnergy[x];
				if (silenceThreshold > 0.0) {
					break;
				}
			}
		}
	}

	//mark silence
	for (x = 0; x < frameCount; x++) {
    start = (unsigned long)floor(pFrames->Mat[x][0]);
    end = (unsigned long)floor(pFrames->Mat[x][0]) + pEnergy->Hdr.frameSize - 1;

		if (pFrames->Mat[x][1] < silenceThreshold) {
      pGT->setSegment(start, end, sclib::atSilence);

      //co-occurence of speech/noise/silence must be prohibited... but the lnk-v/uv-detector is quite cloesely related to this 
      //silence detector (it labels the more energetic parts of what is here called "silence" more correctly as unvoiced speech)
      //thus, it still needs the speech/non-speech labels maybe present, whereas in all other situations, they must be removed
      //the removal of co-occuring speech and silence is cared for in the silence2pause method
      if (this->pTweak->segmentationHandler.vUvDetectorMode != sclib::algorithm_vud_LNK) {
				pGT->setSegment(start, end, SC_GroundTruth::getNoiseRelatedTypes(), false, sclib::noSpeaker, sclib::modeLabelRemove);
				pGT->setSegment(start, end, SC_GroundTruth::getSpeechRelatedTypes(), false, sclib::noSpeaker, sclib::modeLabelRemove);
      }
      
      prob = 1.0 - (pFrames->Mat[x][1] / silenceThreshold); //pseudo-"probability" of this frame being silence
    } else {
      prob = (silenceThreshold / pFrames->Mat[x][1]); //pseudo-"probability" of this frame being silence (yes, NOT non-silence; probs shall always tell be for the label being positive...)
    }

    pGT->setProbability(start, end, sclib::atSilence, prob);
	}

  MFree_1D(averageEnergy);
  MFree_1D(quantiCount);
	MFree_0D(pFrames);
	  
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
unsigned int SC_Segmentation_Silence_LNK::otsuThreshold(unsigned int classCount, unsigned int* itemsPerClass, unsigned long int itemCount) {
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
unsigned long int SC_Segmentation_Silence_LNK::getUncertaintyRegionWidth(void) {
	unsigned long int width = 0;

	//factors due to frame-based analysis
	width += 2 * this->pTweak->featureSte.frameStep;

	return width;
}
