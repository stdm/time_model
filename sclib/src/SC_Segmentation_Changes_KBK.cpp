/**************************************************************************/
/*    Responsibility:																											*/
/*      - This class implements the (speaker) change detector introduced  */
/*        in Kotti, Benetos, Kotropoulos, "Computationally Efficient and  */
/*        Robust BIC-Based Speaker Segmentation", 2008                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 29.08.2008																								*/
/**************************************************************************/

#include "SC_Segmentation_Changes_KBK.h"
#include <list>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_Changes_KBK::SC_Segmentation_Changes_KBK(SC_TweakableParameters* pTweak, double tolerance, int mode) : SC_Segmentation_Changes(pTweak, mode) {
	this->tolerance = tolerance;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_Changes_KBK::~SC_Segmentation_Changes_KBK() {

}

//====================================================================================================================
// the bigger lambda (the BIC threshold) is, the less changepoints are found - this method is meant to be called 
// several times in a row to find an optimal lambda on training data: starting with a very big lambda (e.g.
// std::numeric_limits<double>::max()) that will report no changepoint at all, this method returns a new lambda that 
// will report one changepoint more (the next one reachable). then, F1 (or another global performance measure) can be 
// evaluated on the result and the method can be called once more with that new lambda, and so on. finally, pick the 
// lambda that maximized F1. lambda can't be tweaked analytically because of the dynamic way in which a newly found
// changepoint affects the window-sizes so that everything after it changes.
//====================================================================================================================
//TODO: this method only works correct if in sync with all editings in bicSegmentationKotti()!
double SC_Segmentation_Changes_KBK::findLambdaKotti(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, double r, double tolerance, unsigned long int type, unsigned long int typesNot, unsigned long int boundaryType) {
	SV_Data *pStrippedFeatures, *pStrippedSampleNumbers;
	unsigned long int a, b, r_frames, toleranceInSamples, t, v;
	double lX, lY, lZ, parameterCount, bic, logTerm;
	SC_MatrixFunctions matFunc;
	double *xMean, *yMean, *zMean, **xCov, **yCov, **zCov;
	double d, lambda = std::numeric_limits<double>::max(),  nextLambda = std::numeric_limits<double>::max()*-1.0;
	unsigned long int testedSample, start, end, changePoint;
	long int lastTrueChange = 0, tmp;
	
	//because we want do detect boundaries here, delete all previously set boundary-labels.
	pGT->setSegment(segmentStart, segmentEnd, boundaryType, false, sclib::noSpeaker, sclib::modeLabelRemove);

  //collect all desired frames in this segment
  segmentEnd = sclib::min(segmentEnd, pGT->getConverter()->audioFrame2sample(pGT->getLastAudioFrameNrInSegment(segmentStart, segmentEnd, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep), pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameSize, pFeatures[sclib::bitPosition(sclib::featureMFCC)]->Hdr.frameStep, sclib::alignmentEnd));
  pStrippedFeatures = pGT->copyFramesTogether(pFeatures[this->getUsedFeatures()], segmentStart, segmentStart, segmentEnd, type, typesNot, 0, 0, true);
	pStrippedSampleNumbers = pStrippedFeatures->Next;
	pStrippedFeatures->Next = NULL;
	
	if (pStrippedFeatures != NULL) { //no frames of desired type here

		r_frames = pGT->getConverter()->ms2audioFrame(sclib::round(r*1000.0), pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		toleranceInSamples = pGT->getConverter()->ms2sample(sclib::round(tolerance*1000.0), sclib::alignmentStart);
		a = 0;
		b = 2 * r_frames;
		d = (double)(pStrippedFeatures->Col);
		parameterCount = d + (d*(d+1.0))/2.0;

		while (b <= (unsigned long int)(pStrippedFeatures->Row)) { //the main loop

			//prepare the complete interval's data and model
			zMean = matFunc.mean(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, a, b);
			zCov = matFunc.cov(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, zMean, a, b);
			lZ = ((double)(b-a)/2.0)*sclib::sLog(matFunc.det(zCov, pStrippedFeatures->Col));

			for (t = a+r_frames; t < b-r_frames; t+=r_frames) { //do BIC tests within [a, b] at a resolution (step-size) of "resolution" frames
				//prepare the datasets and models of the individual segments
				xMean = matFunc.mean(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, a, t);
				xCov = matFunc.cov(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, xMean, a, t);
				lX = ((double)(t-a)/2.0)*sclib::sLog(matFunc.det(xCov, pStrippedFeatures->Col));
				yMean = matFunc.mean(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, t, b);
				yCov = matFunc.cov(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, yMean, t, b);
				lY = ((double)(b-t)/2.0)*sclib::sLog(matFunc.det(yCov, pStrippedFeatures->Col));

				//do BIC test
				logTerm = log((double)(b-a));
				bic = lZ-lX-lY - (lambda/2.0)*(parameterCount)*logTerm; //equation 4 from Dellacourt&Wellekens, 2000 (taken *-1 because there a changepoint is found for bic<0!)
				MFree_1D(xMean);
				MFree_2D(xCov);
				MFree_1D(yMean);
				MFree_2D(yCov);

				testedSample = sclib::round(pStrippedSampleNumbers->Mat[t][0]); //pGT->getConverter()->audioFrame2sample(t, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep);
				start = sclib::max(lastTrueChange, (long int)(testedSample-toleranceInSamples/2));
				end = testedSample + toleranceInSamples/2 + toleranceInSamples%2;
				//sclib::tupelOut("se.txt", start, end, this->pTweak);
				if (pGT->testSegment(start, end, false, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) { //a current MD
					//sclib::tupelOut("cp.txt", start, end, this->pTweak);
					pGT->getNextBoundary(start, lastTrueChange, tmp, sclib::atSpeakerBoundary, sclib::searchForward, sclib::modeGroundtruth); //we only need the start of that segment in lastTrueChange, tmp is discarded
					if (bic <= 0.0) {
						lambda = sclib::decrementDouble(2.0 * (lZ-lX-lY) / (parameterCount*logTerm)); //this equality (without the decremt) yields a bic==0; lambda needs to be a little bit smaller to give a bic>0.0
						while (lZ-lX-lY - (lambda/2.0)*(parameterCount)*logTerm <= 0.0) { //make sure that the smallest possible decrement in lambda also gives an increment in bic; if not, iterate decrementing lambda stepwise
							lambda = sclib::decrementDouble(lambda);
						}
						nextLambda = sclib::min(sclib::max(nextLambda, lambda), this->pTweak->segmentationChangesKbk.lambda); //we want the biggest lambda smaller than the current one
					}
				}

				if (bic > 0) { //change point found
					changePoint = testedSample;
					pGT->setSegment(changePoint, changePoint, boundaryType);
					break;
				}
			} //test at all resolution instances t in a, b
			MFree_1D(zMean);
			MFree_2D(zCov);
			
			if (bic > 0) { //change point found
				v = t; //Kotti et al. used (a+b)/2, that misses change points!
				a = v;
				b = sclib::min(v+2*r_frames, pStrippedFeatures->Row-1);
			} else { //no change point found
				b = sclib::min(b+r_frames, pStrippedFeatures->Row-1);
			}
		} //while b not at the end

	} //if there are stripped featrues
	MFree_0D(pStrippedFeatures);
	MFree_0D(pStrippedSampleNumbers);

	if (nextLambda == std::numeric_limits<double>::max()*-1.0) { //no better lambda found => return old one, which can be used as a termination criterion for the loop this method is called in
		nextLambda = this->pTweak->segmentationChangesKbk.lambda;
	}

	return nextLambda;
}

//====================================================================================================================
//	detect changes in the characteristics of the audio-stream
//  pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
//  the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
//====================================================================================================================
int SC_Segmentation_Changes_KBK::detectChanges(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	if (pFeatures[sclib::bitPosition(getUsedFeatures())] == NULL) {
    return SVLIB_Fail;
  } else {
		if (this->mode == sclib::modeSpeakerChange) {
			//return bicSegmentationKotti(pGT, segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationChangesKbk.lambda, this->pTweak->segmentationChangesKbk.r, sclib::noType, sclib::noType, sclib::atSpeakerBoundary); //sclib::atSpeech, sclib::atPause|sclib::atSilence
			return bicSegmentationCettolo(pGT, segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationChangesKbk.lambda, this->pTweak->segmentationChangesKbk.r, sclib::noType, sclib::noType, sclib::atSpeakerBoundary); //sclib::atSpeech, sclib::atPause|sclib::atSilence
		} else {
			return bicSegmentationKotti(pGT, segmentStart, segmentEnd, pFeatures, this->pTweak->segmentationChangesKbk.lambda, this->pTweak->segmentationChangesKbk.r, sclib::atNoise, sclib::noType, sclib::atNoiseBoundary);
		}
  }
}

//====================================================================================================================
//  The segmentation algorithm from Kotti et al.: Kotti, Benetos, Kotropoulos, "Computationally Efficient and Robust 
//  BIC-Based Speaker Segmentation", IEEE Transactions on Audio, Speech, and Language Processing, 2008, 16, 920-933
//  return-value: # of changes
//  lambda is the data-dependant BIC penalty weighting factor
//  r is the mean utterance duration in seconds
//  resultionFactor means: every r/resolutionFacor seconds, a BIC test will be done
//====================================================================================================================
unsigned int SC_Segmentation_Changes_KBK::bicSegmentationKotti(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, double lambda, double r, unsigned long int type, unsigned long int typesNot, unsigned long int boundaryType) {
	SV_Data *pStrippedFeatures, *pStrippedSampleNumbers;
	unsigned long int a, b, r_frames, t, v, toleranceInSamples = pGT->getConverter()->ms2sample(sclib::round(this->tolerance*1000.0));
	double parameterCount, d;
	unsigned int changes = 0;
	SC_MatrixFunctions matFunc;

	//because we want do detect boundarys here, delete all previously set boundary-labels.
	pGT->setSegment(segmentStart, segmentEnd, boundaryType, false, sclib::noSpeaker, sclib::modeLabelRemove);
  
  //mark too short segments as too short (doesn't matter if it has been done before)
	pGT->markShortSpeechSegments(segmentStart, segmentEnd, sclib::max(1, pGT->getConverter()->ms2sample(this->pTweak->general.shortSpeechThreshold, sclib::alignmentEnd)));

  //collect all desired frames in this segment
  segmentEnd = sclib::min(segmentEnd, pGT->getConverter()->audioFrame2sample(pGT->getLastAudioFrameNrInSegment(segmentStart, segmentEnd, pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameSize, pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameStep), pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameSize, pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameStep, sclib::alignmentEnd));
  pStrippedFeatures = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(getUsedFeatures())], segmentStart, segmentStart, segmentEnd, type, typesNot, 0, 0, true);
	pStrippedSampleNumbers = pStrippedFeatures->Next;
	pStrippedFeatures->Next = NULL;

	//construct simple ground truth
	SV_Data *pList = new SV_Data(pStrippedSampleNumbers->Row, 4);
	pList->Hdr = pStrippedSampleNumbers->Hdr;
	for (t = 0; t < (unsigned long int)(pList->Row); t++) {
		pList->Mat[t][0] = pStrippedSampleNumbers->Mat[t][0]; //sample-number
		pList->Mat[t][1] = 0.0f; //1 if a hypothesized change point
		if (pGT->testSegment(sclib::round(pList->Mat[t][0]), sclib::round(pList->Mat[t][0])+pList->Hdr.frameSize, false, sclib::atSpeakerBoundary, false, sclib::atArtificialBoundary, false, sclib::modeGroundtruth) > 0) {
			if (t > 0) {
				pList->Mat[t-1][2] = 0.0f; //delete multiple marked changepoints (due to frame overlap -> take the last one)
			}
			pList->Mat[t][2] = 1.0f; //1 if a groundtruth change point
		} else {
			pList->Mat[t][2] = 0.0f;
		}
		pList->Mat[t][3] = pGT->getConverter()->sample2ms(sclib::round(pList->Mat[t][0])) / 1000.0f; //sample number converted to seconds
	}
	
	if (pStrippedFeatures != NULL) { //no frames of desired type here


		//the resolution "r" in terms of audio-frames instead of seconds
		unsigned long int Nmargin = pGT->getConverter()->ms2audioFrame(150, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		r_frames = pGT->getConverter()->ms2audioFrame(sclib::round(r*1000.0), pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		a = 0; //inclusive
		b = 2 * r_frames; //exclusive
		d = (double)(pStrippedFeatures->Col); //dimensionality of feature vectors
		parameterCount = d + (d*(d+1.0))/2.0; //number of parameters of model in BIC's penalty factor

		while (b < (unsigned long int)(pStrippedFeatures->Row)) { //the main loop
			int Imax = -1;
			double BICimax = 0.0;
			BICimax = getDeltaBIC(pStrippedFeatures, a, b, Nmargin, Nmargin, lambda, Imax, pGT, toleranceInSamples, pStrippedSampleNumbers, true, false);

			/*
			//prepare the complete interval's data and model
			double *zMean = matFunc.mean(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, a, b);
			double **zCov = matFunc.cov(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, zMean, a, b);
			double lZ = ((double)(b-a)/2.0) * sclib::sLog(matFunc.det(zCov, pStrippedFeatures->Col));
			MFree_1D(zMean);
			MFree_2D(zCov);

			double penalty = (lambda/2.0) * parameterCount * log((double)(b-a)); //penalty term for the BIC
			double bic = 0.0;

			for (t = a+Nmargin; t < b-Nmargin; t+=r_frames) { //do BIC tests within [a, b] at a resolution (step-size) of r_frames frames without areas around the borders
				//prepare the datasets and models of the individual segments
				double *xMean = matFunc.mean(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, a, t); //t+1
				double **xCov = matFunc.cov(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, xMean, a, t); //t+1
				double lX = ((double)(t-a)/2.0) * sclib::sLog(matFunc.det(xCov, pStrippedFeatures->Col));
				double *yMean = matFunc.mean(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, t, b); //t+1, b+1
				double **yCov = matFunc.cov(pStrippedFeatures->Mat, pStrippedFeatures->Row, pStrippedFeatures->Col, yMean, t, b); //t+1, b+1
				double lY = ((double)(b-t)/2.0) * sclib::sLog(matFunc.det(yCov, pStrippedFeatures->Col));

				//do BIC test
				bic = lZ-lX-lY - penalty; //equation 4 from Dellacourt&Wellekens, 2000 (taken *-1 because there a changepoint is found for bic<0!); equation (8) from Cettolo et al., 2005
				MFree_1D(xMean);
				MFree_2D(xCov);
				MFree_1D(yMean);
				MFree_2D(yCov);

				//debug output involving groundtruth, acutal algorithmic decision and the two BIC terms
				unsigned long int sample = sclib::round(pStrippedSampleNumbers->Mat[t][0]); // pGT->getConverter()->audioFrame2sample(t, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep);
				unsigned long int st = sclib::max(0, (long int)(sample-toleranceInSamples/2));
				unsigned long int en = sample + toleranceInSamples/2 + toleranceInSamples%2;
				if (pGT->testSegment(st, en, false, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) {
					//                             ground-truth, actual decision, likelihood-ratio, penalty-term
					sclib::quadrupelOut("bic.txt", 1,            bic>0,           lZ-lX-lY,         -penalty, this->pTweak);
				} else {
					sclib::quadrupelOut("bic.txt", 0, bic>0, lZ-lX-lY, -penalty, this->pTweak);
				}

				if (bic > 0) { //change point found
					//BICimax = bic;
					//Imax = t;
					pList->Mat[t][1] = 1.0;
					pGT->setSegment(sample, sample, boundaryType);
					changes++;
					break;
				}
			} //test at all resolution instances t in a, b
			*/
			
			//SC_SignalHandler sig(this->pTweak, sclib::stWave);
			//sig.storeSignal("segment.wav", segmentStart+pGT->getConverter()->audioFrame2sample(a, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep), segmentStart+pGT->getConverter()->audioFrame2sample(b, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd), pGT->getSignalPrototype(), pGT, sclib::noType, sclib::noType, sclib::modeGroundtruth, sclib::atSpeakerBoundary);

			if (Imax>0 && BICimax>0.0) { //change point found
				unsigned long int sample = sclib::round(pStrippedSampleNumbers->Mat[Imax][0]); // pGT->getConverter()->audioFrame2sample(t, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep);
				pGT->setSegment(sample, sample, boundaryType);
				changes++;
				pList->Mat[Imax][1] = 1.0f;
				v = Imax; //Kotti et al. used (a+b)/2, that misses change points! => use Imax as Dellacourt&Wellekens (2000) do
				a = v;
				b = sclib::min(v+2*r_frames, pStrippedFeatures->Row); //Kotti et al. use v+r, but this just yields an empty loop until next increase´´
			} else { //no change point found
				b = sclib::min(b+r_frames, pStrippedFeatures->Row); //2*r_frames
			}
		} //while b not at the end

	} //if there are stripped featrues

	//simple evaluation
	double recall, precision, F, far, mdr;
	simpleEvaluation(pGT, segmentStart, segmentEnd, pList, recall, precision, F, far, mdr);
	printf("\nREC=%f\tPRC=%f\tF1=%f\tMDR=%f\tFAR=%f\n", recall, precision, F, mdr, far);
	char buf[sclib::bufferSize];
	sprintf(buf, "%f\t%f\t%f\t%f\t%f", recall, precision, F, mdr, far);
	sclib::stringOut("res_simple.txt", buf, this->pTweak);

	MFree_0D(pList);
	MFree_0D(pStrippedFeatures);
	MFree_0D(pStrippedSampleNumbers);

	return changes;
}

//====================================================================================================================
//  The segmentation algorithm from Cettolo et al.: Cettolo, Vescovi, Rizzi, "Evaluation of BIC-based Algorithms for 
//  Audio Segmentation", Computer Speech and Language, 2005, 19, 147-170
//  return-value: # of changes
//  lambda is the data-dependant BIC penalty weighting factor
//  r is the mean utterance duration in seconds
//  resultionFactor means: every r/resolutionFacor seconds, a BIC test will be done
//====================================================================================================================
unsigned int SC_Segmentation_Changes_KBK::bicSegmentationCettolo(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data** pFeatures, double lambda, double r, unsigned long int type, unsigned long int typesNot, unsigned long int boundaryType) {
	SV_Data *pStrippedFeatures, *pStrippedSampleNumbers;
	int a, b, changePoint, toleranceInSamples = pGT->getConverter()->ms2sample(sclib::round(this->tolerance*1000.0));
	unsigned int changes = 0;

	//because we want do detect boundarys here, delete all previously set boundary-labels.
	pGT->setSegment(segmentStart, segmentEnd, boundaryType, false, sclib::noSpeaker, sclib::modeLabelRemove);
  
  //mark too short segments as too short (doesn't matter if it has been done before)
	pGT->markShortSpeechSegments(segmentStart, segmentEnd, sclib::max(1, pGT->getConverter()->ms2sample(this->pTweak->general.shortSpeechThreshold, sclib::alignmentEnd)));

  //collect all desired frames in this segment
  segmentEnd = sclib::min(segmentEnd, pGT->getConverter()->audioFrame2sample(pGT->getLastAudioFrameNrInSegment(segmentStart, segmentEnd, pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameSize, pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameStep), pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameSize, pFeatures[sclib::bitPosition(getUsedFeatures())]->Hdr.frameStep, sclib::alignmentEnd));
  pStrippedFeatures = pGT->copyFramesTogether(pFeatures[sclib::bitPosition(getUsedFeatures())], segmentStart, segmentStart, segmentEnd, type, typesNot, 0, 0, true);
	pStrippedSampleNumbers = pStrippedFeatures->Next;
	pStrippedFeatures->Next = NULL;

	//construct simple ground truth
	SV_Data *pList = new SV_Data(pStrippedSampleNumbers->Row, 4);
	pList->Hdr = pStrippedSampleNumbers->Hdr;
	for (int t = 0; t < pList->Row; t++) {
		pList->Mat[t][0] = pStrippedSampleNumbers->Mat[t][0]; //sample-number
		pList->Mat[t][1] = 0.0f; //1 if a hypothesized change point
		if (pGT->testSegment(sclib::round(pList->Mat[t][0]), sclib::round(pList->Mat[t][0])+pList->Hdr.frameSize, false, sclib::atSpeakerBoundary, false, sclib::atArtificialBoundary, false, sclib::modeGroundtruth) > 0) {
			if (t > 0) {
				pList->Mat[t-1][2] = 0.0f; //delete multiply marked changepoints (due to frame overlap -> take the last one)
			}
			pList->Mat[t][2] = 1.0f; //1 if a groundtruth change point
		} else {
			pList->Mat[t][2] = 0.0f;
		}
		pList->Mat[t][3] = pGT->getConverter()->sample2ms(sclib::round(pList->Mat[t][0])) / 1000.0f; //sample number converted to seconds
	}
	
	if (pStrippedFeatures != NULL) { //no frames of desired type here

		//parameters according to Cettolo et al.
		int Nmin = 200; //pGT->getConverter()->ms2audioFrame(1000, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		int Nmax = 500; //pGT->getConverter()->ms2audioFrame(20000, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		int Nmargin = 50; //pGT->getConverter()->ms2audioFrame(50, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		int Nshift = 100; //pGT->getConverter()->ms2audioFrame(100, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		int Ngrow = 50; //pGT->getConverter()->ms2audioFrame(75, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		int Nsecond = 400; //pGT->getConverter()->ms2audioFrame(10000, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep, sclib::alignmentEnd);
		int delta1 = 25;
		int delta2 = delta1 / 5;

		a = 0; //inclusive
		b = Nmin; //exclusive
		double d = (double)(pStrippedFeatures->Col); //dimensionality of feature vectors
		int T = pStrippedFeatures->Row;

		//sclib::classOut("features.txt", pStrippedFeatures, this->pTweak);

		while (b <= T) { //the main loop
			int iMax;
			double bic = getDeltaBIC(pStrippedFeatures, a, b, delta1, Nmargin, lambda, iMax, pGT, toleranceInSamples, pStrippedSampleNumbers);
			
			while (bic<=0.0 && b-a<Nmax && b<T) {
				b = sclib::min(b+Ngrow, pStrippedFeatures->Row);
				bic = getDeltaBIC(pStrippedFeatures, a, b, delta1, Nmargin, lambda, iMax, pGT, toleranceInSamples, pStrippedSampleNumbers);
			}

			while (bic<=0.0 && b<T) {
				b = sclib::min(b+Nshift, T);
				a = sclib::min(a+Nshift, b-Nmin);
				bic = getDeltaBIC(pStrippedFeatures, a, b, delta1, Nmargin, lambda, iMax, pGT, toleranceInSamples, pStrippedSampleNumbers);
			}

			if (bic > 0) { //change point found
				//printf("\n%d", iMax);

				//centering of window
				if (b-a > Nsecond) {
					b = a + Nsecond;
				}
				a = sclib::max(a, iMax - ((b-a)/2));
				b = sclib::min(iMax+(iMax-a), T);

				bic = getDeltaBIC(pStrippedFeatures, a, b, delta2, Nmargin, lambda, iMax, pGT, toleranceInSamples, pStrippedSampleNumbers);
				
				if (bic > 0) {
					//output changepoint
					pList->Mat[iMax][1] = 1.0;
					changePoint = sclib::round(pStrippedSampleNumbers->Mat[iMax][0]); // pGT->getConverter()->audioFrame2sample(t, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep);
					pGT->setSegment(changePoint, changePoint, boundaryType);
					changes++;

					a = iMax + 1; 
					b = sclib::min(a+Nmin, T);
				} else {
					a = sclib::max(a+1, iMax-Nmargin+1);
					b = sclib::min(a+Nmin, T);
				}				
			}

			if (b == T) {
				break;
			}
		} //while b not at the end

	} //if there are stripped featrues

	double recall, precision, F, far, mdr;
	simpleEvaluation(pGT, segmentStart, segmentEnd, pList, recall, precision, F, far, mdr);
	printf("\nREC=%f\tPRC=%f\tF1=%f\tMDR=%f\tFAR=%f\n", recall, precision, F, mdr, far);
	char buf[sclib::bufferSize];
	sprintf(buf, "%f\t%f\t%f\t%f\t%f", recall, precision, F, mdr, far);
	sclib::stringOut("res_simple.txt", buf, this->pTweak);

	MFree_0D(pList);
	MFree_0D(pStrippedFeatures);
	MFree_0D(pStrippedSampleNumbers);

	return changes;
}

//====================================================================================================================
// Computes all delta-BIC values between start (inclusive) and end (exclusive) and returns the instant of and the 
// biggest BIC value; calls isTrueChangePoint() for debug output
//====================================================================================================================
double SC_Segmentation_Changes_KBK::getDeltaBIC(SV_Data *pFeatures, int startInclusive, int endExclusive, int resolutionInFrames, int Nmargin, double lambda, int &bestPos, SC_GroundTruth *pGT, int toleranceInSamples, SV_Data *pSampleNumbers, bool verbose, bool lockToTrueCp) {
	int t = 0, a = startInclusive, b = endExclusive; 
	unsigned int changes = 0;
	SC_MatrixFunctions matFunc;
	double d = pFeatures->Col;

	//prepare the complete interval's data and model
	double *zMean = matFunc.mean(pFeatures->Mat, pFeatures->Row, pFeatures->Col, a, b);
	double **zCov = matFunc.cov(pFeatures->Mat, pFeatures->Row, pFeatures->Col, zMean, a, b);
	double lZ = ((double)(b-a)/2.0) * sclib::sLog(matFunc.det(zCov, pFeatures->Col));
	MFree_2D(zCov);
	MFree_1D(zMean);

	double parameterCount = d + (d*(d+1.0))/2.0; //number of parameters of model in BIC's penalty factor
	double penalty = (lambda/2.0) * parameterCount * log((double)(b-a)); //penalty term for the BIC
	double BICimax = std::numeric_limits<double>::max()*-1.0;
	bestPos = -1;
	
	for (t = a+Nmargin; t <= b-Nmargin; t+=resolutionInFrames) { //do BIC tests within [a, b] at a resolution (step-size) of r_frames frames without areas around the borders
		//prepare the datasets and models of the individual segments
		double *xMean = matFunc.mean(pFeatures->Mat, pFeatures->Row, pFeatures->Col, a, t); //t+1
		double **xCov = matFunc.cov(pFeatures->Mat, pFeatures->Row, pFeatures->Col, xMean, a, t); //t+1
		double lX = ((double)(t-a)/2.0) * sclib::sLog(matFunc.det(xCov, pFeatures->Col));
		double *yMean = matFunc.mean(pFeatures->Mat, pFeatures->Row, pFeatures->Col, t, b); //t+1, b+1
		double **yCov = matFunc.cov(pFeatures->Mat, pFeatures->Row, pFeatures->Col, yMean, t, b); //t+1, b+1
		double lY = ((double)(b-t)/2.0) * sclib::sLog(matFunc.det(yCov, pFeatures->Col));

		//do BIC test
		double bic = lZ-lX-lY - penalty; //equation 4 from Dellacourt&Wellekens, 2000 (taken *-1 because there a changepoint is found for bic<0!); equation (8) from Cettolo et al., 2005
		MFree_1D(xMean);
		MFree_2D(xCov);
		MFree_1D(yMean);
		MFree_2D(yCov);

		if (bic > BICimax) { //change point found
			BICimax = bic;
			bestPos = t;
		}
		if (isTrueChangePoint(pGT, sclib::round(pSampleNumbers->Mat[t][0]), toleranceInSamples, verbose, bic, lX, lY, lZ, penalty) == true) { //just for debug output if lockToTrueCp==false
			if (lockToTrueCp == true) {
				BICimax = bic;
				bestPos = t;
				break;
			}
		}
	} //test at all resolution instances t in a, b

	return BICimax;
}

//====================================================================================================================
// Returns true if the candidate change point at sampleNr is a true cp; may also write debug output
//====================================================================================================================
bool SC_Segmentation_Changes_KBK::isTrueChangePoint(SC_GroundTruth *pGT, unsigned long int sampleNr, unsigned long int toleranceInSamples, bool verbose, double bic, double lX, double lY, double lZ, double penalty) {
	bool res = false;

	//debug output involving groundtruth, acutal algorithmic decision and the two BIC terms
	//unsigned long int sample = sclib::round(pStrippedSampleNumbers->Mat[t][0]); // pGT->getConverter()->audioFrame2sample(t, pStrippedFeatures->Hdr.frameSize, pStrippedFeatures->Hdr.frameStep);
	unsigned long int st = sclib::max(0, (long int)(sampleNr-toleranceInSamples/2));
	unsigned long int en = sampleNr + toleranceInSamples/2 + toleranceInSamples%2;
	if (pGT->testSegment(st, en, false, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) {
		if (verbose == true) {
			//                             ground-truth, actual decision, likelihood-ratio, penalty-term
			sclib::quadrupelOut("bic.txt", 1,            bic>0,           lZ-lX-lY,         -penalty, this->pTweak);
		}
		res = true;
	} else {
		if (verbose == true) {
			sclib::quadrupelOut("bic.txt", 0, bic>0, lZ-lX-lY, -penalty, this->pTweak);
		}
		res = false;
	}

	return res;
}

//====================================================================================================================
// Simple but correct evaluation script to test SC_Score_ChangeDetection
//====================================================================================================================
void SC_Segmentation_Changes_KBK::simpleEvaluation(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pStrippedSampleNumbers, double &recall, double &precision, double &F, double &far, double &mdr) {
	unsigned long int start, end, changePoints = 0;
	unsigned int frameLength = pStrippedSampleNumbers->Hdr.frameSize;
	int scatterMatrix[2][2] = {{0, 0}, {0, 0}}; //[TP FP]
	                                            //[FN TN]

	//sclib::classOut("list.txt", pStrippedSampleNumbers, this->pTweak);

	//this experimental part expects a strippedSampleNumber Matrix not as returned by copyFramesTogether(), but with additional
	//3 columns: [1] 1<=>hypo-cp, [2] 1<=>gt-cp, [3] time-index of this frame in seconds
	//additional assumption: the frames must!!! be concesutive!
	int toleranceInFrames = pGT->getConverter()->ms2audioFrame(sclib::round(this->tolerance*1000.0), pStrippedSampleNumbers->Hdr.frameSize, pStrippedSampleNumbers->Hdr.frameStep);
	int negativeTolerance = toleranceInFrames/2, positiveTolerance = toleranceInFrames/2 + toleranceInFrames%2;
	long int lastCheck = -1;
	for (int t = 0; t < pStrippedSampleNumbers->Row; t++) {
		if (sclib::round(pStrippedSampleNumbers->Mat[t][2]) > 0) {  //a gt cp
			//find next gt cp
			int t_next;
			for (t_next = t+1; t_next < pStrippedSampleNumbers->Row; t_next++) {
				if (sclib::round(pStrippedSampleNumbers->Mat[t_next][2]) > 0) {
					break;
				}
			}
		
			//the tolerance area around the current true positive
			//+-half the tolerance, but:
			//  maximally halfway betwen this and the next true change point
			//  minimally starting at the last tolerance range end
			long int toleranceStart = sclib::max(lastCheck+1, t-negativeTolerance);
			long int toleranceEnd = sclib::min(pStrippedSampleNumbers->Row-1, t+positiveTolerance);
			if (toleranceEnd>=t_next && t_next<pStrippedSampleNumbers->Row) {
				toleranceEnd = (t+t_next) / 2;
			}

			//the false alarms that occured between the last boundary (+tolerance) until the current boundary (-tolerance)
			start = lastCheck + 1;
			end = sclib::max(lastCheck+1, t-negativeTolerance-1);
			for (unsigned long int tt = start; tt <= end; tt++) {
				if (sclib::round(pStrippedSampleNumbers->Mat[tt][1]) > 0) {
					scatterMatrix[0][1]++;
					//sclib::scalarOut("fa.txt", pStrippedSampleNumbers->Mat[tt][0], this->pTweak);
				}
			}

			//true positives (and more false alarms) in the tolerance range
			//starting from the true cp, hop look left and right (in the bounds) for the nearest hypo-cp and mark it as true
			//all others are false alarms
			start = toleranceStart;
			end = toleranceEnd;
			int count = 0;
			int current = t;
			int sign = -1;
			int offset = 0;
			bool rightOut = false;
			bool leftOut = false;
			while (!rightOut && !leftOut) {
				if (current > toleranceEnd) {
					rightOut = true;
				} else if (current < toleranceStart) {
					leftOut = true;
				} else if (sclib::round(pStrippedSampleNumbers->Mat[current][1]) > 0) {
					if (count == 0) {
						scatterMatrix[0][0]++;
						//sclib::scalarOut("tp.txt", pStrippedSampleNumbers->Mat[current][0], this->pTweak);
					} else {
						scatterMatrix[0][1]++;
						//sclib::scalarOut("fa.txt", pStrippedSampleNumbers->Mat[current][0], this->pTweak);
					}
					count++;					
				}
				sign *= -1;
				offset++;
				current += sign*offset;
			}

			//false negatives (=missed detections)
			if (count == 0) {
				scatterMatrix[1][0]++;
				//sclib::scalarOut("fn.txt", pStrippedSampleNumbers->Mat[t][0], this->pTweak);
			}

			changePoints++;
			lastCheck = toleranceEnd;
		}
	}

	scatterMatrix[1][1] = pStrippedSampleNumbers->Row - changePoints;
	
	//hits = pGT->testSegment(segmentStart, segmentEnd, false, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeHypothesized) / pGT->getInternalFrameSize();
	//int hits2 = pGT->testSegment(segmentStart, segmentEnd, false, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeGroundtruth) / pGT->getInternalFrameSize();
				
	//compute the scores;
	recall = (double)(scatterMatrix[0][0]) / (double)(scatterMatrix[0][0]+scatterMatrix[1][0]); //REC = CFC/(CFC+MD), equation (14), CFC=true positive
	precision = (double)(scatterMatrix[0][0]) / (double)(scatterMatrix[0][0]+scatterMatrix[0][1]); //PRC = CFC/(CFC+FA), equation (14), CFC=true positive;
	F = ((recall+precision)>0) ? ((2.0*precision*recall) / (precision+recall)) : 0.0; //equation (14)
	far = (double)(scatterMatrix[0][1]) / (double)(changePoints + scatterMatrix[0][1]); //FAR = FA/(GT+FA), equation (13) from Kotti's paper
	mdr = (double)(scatterMatrix[1][0]) / (double)(changePoints); //MDR = MD/GT, equation (13)

	return;
}

//====================================================================================================================
// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
// Only the factors due to this specific algorithm are taken into account (factors regarding the ground-truth class
// are handled therein)
//====================================================================================================================
unsigned long int SC_Segmentation_Changes_KBK::getUncertaintyRegionWidth(void) {
	unsigned long int width = 0;
	SC_Conversion converter;

	//factors due to frame-based analysis
	width += 2 * this->pTweak->segmentationChangesKbk.mfccParameters.frameStep;

	//factors due to windowed analysis
	width += 2 * sclib::round(1000.0 * this->pTweak->segmentationChangesKbk.r);

	return sclib::round(1000 * this->tolerance); //sclib::max(width, 1000.0*this->tolerance);
}

//====================================================================================================================
// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
// parameters
//====================================================================================================================
SC_TweakableParameters::SC_FeaturePar* SC_Segmentation_Changes_KBK::getSpecialFeatureParameters(void) {
	this->pTweak->segmentationChangesKbk.mfccParameters.Next = NULL;

	return &(this->pTweak->segmentationChangesKbk.mfccParameters);
}
