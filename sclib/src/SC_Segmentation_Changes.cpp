/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for change detectors (speaker/acoustic)       				*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#include "SC_Segmentation_Changes.h"
#include "SC_Aux.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_Changes::SC_Segmentation_Changes(SC_TweakableParameters* pTweak, int mode) {
  this->pTweak = pTweak;
	this->mode = mode;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_Changes::~SC_Segmentation_Changes() {

}

//====================================================================================================================
//	If any algorithm detected a change at this window according to the previous ones, this function looks where in 
//  this window the real changepoint may lie (e.g. the speech in this window is connected, then the changepoint is at 
//  the specified refPoint of the window; if this window consists of many speech-segments, the changepoint is 
//  somewhere at the segement-borders... we set it at the beginning of the segment that contributes most frames to the 
//  window
//  input- and output-values are, as always, sample-based
//====================================================================================================================
unsigned long int SC_Segmentation_Changes::refineChangePoint(SC_GroundTruth *pGT, unsigned long int windowStart, unsigned long int windowEnd, unsigned long int oldWindowEnd, unsigned long int type, int refPoint) {
	unsigned long int changePoint = (refPoint == sclib::refpointMiddle) ? windowStart + (windowEnd-windowStart+1)/2 : windowStart; //default value
  long int segmentStart, segmentEnd, contribution, maxContribution = 0;
  //bool justPause = true;
	bool firstTime = true;
	
	pGT->getNextSegment(windowStart, segmentStart, segmentEnd, type, sclib::searchMiddle);

	//Search for the segment contributing most samples to the current window and set the changepoint at it's beginning (but inside the window)
	//if none is found (should never hit due to the way the windows are constructed from features) the changepoint is set to be at the window's refPoint
	while (segmentStart <= (long int)(windowEnd) && segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
		contribution = sclib::min(segmentEnd, windowEnd) - sclib::max(segmentStart, windowStart) + 1;
		if (contribution > maxContribution) {
			maxContribution = contribution;
			if (firstTime == true && refPoint == sclib::refpointMiddle) {
				changePoint = sclib::max(windowStart, segmentStart) + (sclib::min(segmentEnd, windowEnd)-sclib::max(segmentStart, windowStart)+1)/2;
			} else {
				changePoint = sclib::max(windowStart, segmentStart);
			}
		}
		pGT->getNextSegment(segmentEnd+1, segmentStart, segmentEnd, type, sclib::searchForward);
		firstTime = false;
	}

	/*
  if (segmentEnd < (long int)windowEnd) {
    if (segmentStart > (long int)oldWindowEnd) { //don't move the changepoint into the scope of the last window because the change was detected in this one!!!
      changePoint = segmentStart;
      maxContribution = sclib::min(segmentEnd, windowEnd) - sclib::max(segmentStart, windowStart) + 1;
    }

    //loop over all segments which may start within this window
    while (segmentEnd < (long int)windowEnd && segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
      pGT->getNextSegment(segmentEnd+1, segmentStart, segmentEnd, type, sclib::searchForward);    

      if (segmentStart < (long int)windowEnd && segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
        if (sclib::min(segmentEnd, windowEnd) - sclib::max(segmentStart, windowStart) > maxContribution) {
          maxContribution = sclib::min(segmentEnd, windowEnd) - sclib::max(segmentStart, windowStart) + 1;
          changePoint = segmentStart;
        }
      }
    }
  }

  //if considering a speaker-change, move the change-point forward if it is near the beginning of a speech-segment, 
  //but preceeded by some pause-frames: then, set it at the segment's beginning
  //this problem arises because for building up the windows only non-pause frames are used
  if (type == sclib::atSpeech) {
    pGT->getNextSegment(changePoint, segmentStart, segmentEnd, type, sclib::searchMiddle);
    if (segmentStart < (long int)changePoint && segmentEnd > (long int)changePoint && segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
      if (pGT->testSegment(segmentStart, segmentEnd, true, sclib::atSpeech|sclib::atPause, true) == 0) {
        justPause = false;
      }
      if (justPause == true && segmentStart > (long int)oldWindowEnd) { //don't move the changepoint into the scope of the last window because the change was detected in this one!!!
        changePoint = segmentStart;
      }
    }
  }
	*/

  return changePoint;

	//return windowStart + ((windowEnd - windowStart + 1) / 2);
}

//====================================================================================================================
// Returns a non-parametric density-estimation of the length of speech/noise segments between speaker/acoustic changes
// type has to be sclib::atSpeech or sclib::atNoise
//====================================================================================================================
SC_Model_Pareto* SC_Segmentation_Changes::getSegmentDurationDistribution(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int type) {
	long int boundaryType = (type == sclib::atSpeech) ? sclib::atSpeakerBoundary : sclib::atNoiseBoundary;
	long int types = (boundaryType == sclib::atSpeakerBoundary) ? sclib::atSpeech : sclib::atNoise;
	long int typesNot = (types == sclib::atSpeech) ? sclib::atPause|sclib::atSilence : sclib::noType;
	long int y, boundaryStart, boundaryEnd, boundaryCount = pGT->testSegment(segmentStart, segmentEnd, false, boundaryType, false, sclib::noType, false, sclib::modeGroundtruth) / pGT->getInternalFrameSize(); 
	SC_Model_Pareto *pModel = NULL;
	SV_Data *pDurations = NULL;

	if (boundaryCount > 0) {
		pDurations = new SV_Data(boundaryCount, 1);

		boundaryCount = 0;
		for (y = (long int)(segmentStart); y <= (long int)(segmentEnd); y++) {
			pGT->getNextBoundary(y, boundaryStart, boundaryEnd, boundaryType, sclib::searchForward, sclib::modeGroundtruth);
			if (boundaryStart != sclib::noSegment && boundaryEnd != sclib::noSegment && boundaryStart <= (long int)(segmentEnd)) {
				pDurations->Mat[boundaryCount][0] = (float)(pGT->testSegment(boundaryStart, boundaryEnd, false, types, false, typesNot, false, sclib::modeGroundtruth));
				pDurations->Mat[boundaryCount][0] = (float)(pGT->getConverter()->sample2ms((unsigned long int)(pDurations->Mat[boundaryCount][0])) / 1000.0); //in seconds

				boundaryCount++;
				y = boundaryEnd;
			} else {
				break;
			}
		}

		pModel = new SC_Model_Pareto(this->pTweak, NULL);
		pModel->TrainModel(pDurations);
	}

	if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true) {
		sclib::matrixOut("durations.txt", pDurations->Mat, pDurations->Row, 1, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
	}

	MFree_0D(pDurations);

	return pModel;
}

//====================================================================================================================
// Returns a non-parametric density-estimation of the detectability of the changepoints for the given type;
// type has to be sclib::atSpeech or sclib::atNoise
//====================================================================================================================
SC_Model_Pareto* SC_Segmentation_Changes::getDetectability(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int type) {
	long int boundaryType = (type == sclib::atSpeech) ? sclib::atSpeakerBoundary : sclib::atNoiseBoundary;
	long int types = (boundaryType == sclib::atSpeakerBoundary) ? sclib::atSpeech : sclib::atNoise;
	long int typesNot = (types == sclib::atSpeech) ? sclib::atPause|sclib::atSilence : sclib::noType;
	long int y, boundaryStart, boundaryEnd, boundaryCount = pGT->testSegment(segmentStart, segmentEnd, false, boundaryType, false, sclib::noType, false, sclib::modeGroundtruth) / pGT->getInternalFrameSize(); 
	SC_Model_Pareto *pModel = NULL;
	SV_Data *pDetectability = NULL;

	if (boundaryCount > 0) {
		pDetectability = new SV_Data(boundaryCount, 1);

		//get the length of all segments of given type between 2 boundaries
		boundaryCount = 0;
		for (y = (long int)(segmentStart); y <= (long int)(segmentEnd); y++) {
			pGT->getNextBoundary(y, boundaryStart, boundaryEnd, boundaryType, sclib::searchForward, sclib::modeGroundtruth);
			if (boundaryStart != sclib::noSegment && boundaryEnd != sclib::noSegment && boundaryStart <= (long int)(segmentEnd)) {
				pDetectability->Mat[boundaryCount][0] = (float)(pGT->testSegment(boundaryStart, boundaryEnd, false, types, false, typesNot, false, sclib::modeGroundtruth));
				pDetectability->Mat[boundaryCount][0] = (float)(pGT->getConverter()->sample2ms((unsigned long int)(pDetectability->Mat[boundaryCount][0])) / 1000.0); //in seconds

				boundaryCount++;
				y = boundaryEnd;
			} else {
				break;
			}
		}

		//calculate the detectability as suggested in S.S. Chen, P.S. Gopalakrishnanm "Speaker, Environment and Channel Change Detection and Clustering via the Bayesian Information Criterion", Proc. DARPA Broadcast News Transcription and Understanding Workshop, Landsdowne, VA, 1998, pp. 127-132
		for (y = 0; y < pDetectability->Row-1; y++) {
			pDetectability->Mat[y][0] = sclib::min(pDetectability->Mat[y][0], pDetectability->Mat[y+1][0]);
		}
		pDetectability->Mat[pDetectability->Row-1][0] = 0.0;

		pModel = new SC_Model_Pareto(this->pTweak, NULL);
		pModel->TrainModel(pDetectability);
	}

	if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbSCD) == true) {
		sclib::matrixOut("detectability.txt", pDetectability->Mat, pDetectability->Row, 1, this->pTweak, 0, 0, 0, 0, sclib::matlabSyntax);
	}

	MFree_0D(pDetectability);

	return pModel;
}
