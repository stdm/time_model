/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes overall speaker-classification score, so to speak the  */
/*        final results if the SCiVo application.													*/
/*        To accomplish this, this results are no more dependant on the   */
/*        results of previous algorithms; rather, a direct comparison of  */
/*        ground-truth speech segments and algorthmic results is done     */
/*        (of course also for different CEs).                             */
/*																																				*/
/*      - ATTENTION: It is assumed that proper speaker-mappings already   */
/*                   exist in the ground-truth class (created via         */
/*                   speaker-clustering-/-id-scoring)!                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.12.2006																								*/
/**************************************************************************/

#include "SC_SpeakerScore_Classification.h"
#include "SC_Aux.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_SpeakerScore_Classification::SC_SpeakerScore_Classification(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT) : SC_SpeakerScore(pTweak, pGT) {
	this->ca = this->pGT->getRealSpeakerCount(true) + 1; //+1 for the NO_SPEAKER (subsumes all non-speech- and falsely classified speech-segments)
	
	for (unsigned int ce = 0; ce < 4; ce++) {
		this->scatterMatrix[ce] = NULL;
		this->averageOmission[ce] = 0.0;
		this->averageCommission[ce] = 0.0;
		this->missclassificationRate[ce]= 0.0;
		this->kappaStatistic[ce] = 0.0;
		this->recall[ce] = NULL;
		this->precision[ce] = NULL;
		this->missRate[ce] = NULL;
		this->falseAlarmRate[ce] = NULL;
		this->errorRate[ce] = NULL;
		this->specificity[ce] = NULL;
		this->accuracy[ce] = NULL;
	}
	this->DER = 0.0;

	sprintf(this->purpose, "%s\0", "Final (speaker-classification) Assessment");
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_SpeakerScore_Classification::~SC_SpeakerScore_Classification() {
	for (unsigned int ce = 0; ce < 4; ce++) {
		MFree_2D(this->scatterMatrix[ce]);
		MFree_1D(this->recall[ce]);
		MFree_1D(this->precision[ce]);
		MFree_1D(this->missRate[ce]);
		MFree_1D(this->falseAlarmRate[ce]);
		MFree_1D(this->errorRate[ce]);
		MFree_1D(this->specificity[ce]);
		MFree_1D(this->accuracy[ce]);
	}
}

//====================================================================================================================
//	Fills the internal score-variables by computing their values according to the frameList in pGT, so that the
//  get*()-Functions return reasonable values (before calling calcScores(), they return all 0)
//  "start" and "end" refer to sample-numbers so that the area of the frameList for which scores shall be computed can
//  be specified; this way, scores can be calculated only for parts of the video/corpus, e.g. for a scene.
//  In algorithmicUncertaintyDiameter a value [in samples!] can be given which describes the precision with which the
//  specific algorithm responsible for the results (and only the algorithm, not the gt...) can predict the place of 
//  event-on- and -offsets
//====================================================================================================================
void SC_SpeakerScore_Classification::calcScores(unsigned long int start, unsigned long int end, unsigned long int algorithmicUncertaintyDiameter) {
	long int gtSegStart, gtSegEnd, hypoSegStart, hypoSegEnd, oldGtSegStart = sclib::noSegment, gtID, hypoID;
	unsigned long int softBoundaryDiameter = this->pGT->getUncertaintyRegionWidth(true) + algorithmicUncertaintyDiameter;
	long int idxGt, idxHypo;
	unsigned long int x, y, scene = 0, lastScene = 0, shot = 0, lastShot = 0;
  unsigned long int msPerGaussian, modelThreshold, clusterThreshold, segLength, surroundingSegLength;
	int ce;
	bool **hitsInScene, **hitsInShot, gtEndIsShotBoundary, firstHypoSegment, correct;
	SC_MatrixFunctions *pFunc = NULL;

	//initializations
	this->start = start;
	this->end = (end > 0 && end < this->pGT->getAudioSampleCount()) ? end : this->pGT->getAudioSampleCount()-1;
	this->ca = this->pGT->getRealSpeakerCount(true, this->start, this->end) + 1; //+1 for the NO_SPEAKER (subsumes all non-speech- and falsely classified speech-segments)
	pFunc = new SC_MatrixFunctions();
	for (ce = 0; ce < 4; ce++) {
		MFree_2D(this->scatterMatrix[ce]);
		this->scatterMatrix[ce] = pFunc->initMatrix(this->ca+2, this->ca, (long int)(0));
	}
	MArray_2D(hitsInScene, (long int)(this->ca)+2, this->ca, bool, "SC_SpeakerScore_Classification.calcScores: hitsInScene");
	MArray_2D(hitsInShot, (long int)(this->ca)+2, this->ca, bool, "SC_SpeakerScore_Classification.calcScores: hitsInShot");
	for (y = 0; y < this->ca+2; y++) {
		for (x = 0; x < this->ca; x++) {
			hitsInScene[y][x] = false;
			hitsInShot[y][x] = false;
		}
	}
	msPerGaussian = this->pGT->getConverter()->ms2sample(this->pTweak->modelHandler.msPerGaussian); //this corresponds with complete speaker-segments, i.e. all single speech-segments between consecutive speaker-boundaries
	modelThreshold = this->pGT->getConverter()->ms2sample(this->pTweak->general.shortSpeechThreshold); //this corresponds with single speech-segments
	clusterThreshold = this->pGT->getConverter()->ms2sample(this->pTweak->speakerClusterer.speechSegLengthThreshold); //this also corresponds with complete speaker-segments

	//we compare gt- and hypo-based speaker-id's ce-wise
	for (y = this->start; y < this->end; y++) {
		this->pGT->getNextSegment(y, gtSegStart, gtSegEnd, sclib::noType, sclib::searchForward, sclib::modeGroundtruth, true, true, sclib::atShotBoundary); //get a segment in which the gt-speaker-id is homogeniuos (can be sclib::noSpeaker)
		if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment && gtSegStart <= (long)(this->end)){
			gtID = this->pGT->getSamplesSpeakerID(gtSegStart, sclib::modeGroundtruth);

			//is it long enough so that it got evaluated?
			segLength = (unsigned long int)(gtSegEnd - gtSegStart + 1);
			surroundingSegLength = this->pGT->getSurroundingSpeakersSegmentSize(gtSegStart, gtSegEnd);
			if ((segLength > modelThreshold || modelThreshold == 0) && 
					(surroundingSegLength > clusterThreshold || clusterThreshold == 0) &&
					(surroundingSegLength > msPerGaussian)) {

				//check for scene- and shot-changes
				scene = this->pGT->sample2scene(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastScene);
				if (scene > this->pTweak->general.lastScene) {
					break;
				}

				//has this scene been processed by any algorithm?
				if (scene >= this->pTweak->general.firstScene && (this->pTweak->general.sceneSelection == 0 || sclib::bitTest(sclib::bit(scene), pTweak->general.sceneSelection) == true) && scene <= this->pTweak->general.lastScene) {

					idxGt = class2idx(gtID);
					gtEndIsShotBoundary = (gtSegEnd == this->pGT->getAudioSampleCount()-1 || this->pGT->testSegment(gtSegEnd+1, gtSegEnd+1, true, sclib::atShotBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) ? true : false;
					firstHypoSegment = true;

					if (scene != lastScene) {
						pFunc->clear(hitsInScene, this->ca+2, this->ca, false); //clear rememberance of hits
					}
					shot = this->pGT->sample2shot(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastShot);
					if (shot != lastShot) {
						pFunc->clear(hitsInShot, this->ca+2, this->ca, false); //clear rememberance of hits
					}

					//handle population/possible sets
					//population and possible set are the same for the u/uv task 
					//(if we take all samples between start and end as population, how shal we split them to u/uc columns?)
					//(and actually we decided to just score segments whoch are SPEECH in both gt and hypo, so really pop==pos within this decision)
					this->scatterMatrix[sclib::ceSample][this->ca][idxGt] += gtSegEnd - gtSegStart + 1;
					this->scatterMatrix[sclib::ceSample][this->ca+1][idxGt] += gtSegEnd - gtSegStart + 1;
					this->scatterMatrix[sclib::ceSegment][this->ca][idxGt]++;
					this->scatterMatrix[sclib::ceSegment][this->ca+1][idxGt]++;
					if (hitsInShot[this->ca][idxGt] == false) {
						this->scatterMatrix[sclib::ceShot][this->ca][idxGt]++;
						hitsInShot[this->ca][idxGt] = true;
					}
					if (hitsInShot[this->ca+1][idxGt] == false) {
						this->scatterMatrix[sclib::ceShot][this->ca+1][idxGt]++;
						hitsInShot[this->ca+1][idxGt] = true;
					}
					if (hitsInScene[this->ca][idxGt] == false) {
						this->scatterMatrix[sclib::ceScene][this->ca][idxGt]++;
						hitsInScene[this->ca][idxGt] = true;
					}
					if (hitsInScene[this->ca+1][idxGt] == false) {
						this->scatterMatrix[sclib::ceScene][this->ca+1][idxGt]++;
						hitsInScene[this->ca+1][idxGt] = true;
					}

					for (x = (unsigned long int)(gtSegStart); x <= (unsigned long int)(gtSegEnd); x++) {
						this->pGT->getNextSegment(x, hypoSegStart, hypoSegEnd, sclib::noType, sclib::searchWithin, sclib::modeHypothesized, true, true, sclib::atShotBoundary);
						if (hypoSegStart != sclib::noSegment && hypoSegEnd != sclib::noSegment && hypoSegStart <= (long)(gtSegEnd)) {
							hypoID = this->pGT->getSamplesSpeakerID(hypoSegStart, sclib::modeHypothesized);
							gtID = this->pGT->getSpeakerGIDFromHID(hypoID, correct); //get corresponding gt-id
							idxHypo = (correct == true) ? class2idx(gtID) : class2idx(sclib::noSpeaker); //count the segment as a failure if the mapping is not "correct" (see SC_SpeakerScore_Clutsering for explanation of this term)

							//care for soft boundaries
							if (abs(gtSegStart - hypoSegStart) <= (long int)(softBoundaryDiameter) && firstHypoSegment == true) {
								hypoSegStart = gtSegStart;
							}
							if (abs(gtSegEnd - hypoSegEnd) <= (long int)(softBoundaryDiameter) && gtEndIsShotBoundary == false) {
								hypoSegEnd = gtSegEnd;
							}

							//fill the scatter-matrix for all 4 CEs
							if (gtSegStart == hypoSegStart && gtSegEnd == hypoSegEnd) { //here we still need the knowledge if the segments fitted correctly
								this->scatterMatrix[sclib::ceSegment][idxHypo][idxGt]++;
							}
							hypoSegEnd = sclib::min(hypoSegEnd, gtSegEnd); //from here, the hypoSegment must not exceed the gtSegment!
							this->scatterMatrix[sclib::ceSample][idxHypo][idxGt] += (hypoSegEnd - hypoSegStart + 1);
							if (hitsInShot[idxHypo][idxGt] == false) {
								this->scatterMatrix[sclib::ceShot][idxHypo][idxGt]++;
								hitsInShot[idxHypo][idxGt] = true;
							}
							if (hitsInScene[idxHypo][idxGt] == false) {
								this->scatterMatrix[sclib::ceScene][idxHypo][idxGt]++;
								hitsInScene[idxHypo][idxGt] = true;
							}

							firstHypoSegment = false;
							x = hypoSegEnd;
						} else {
							break;
						}					
					}

				} //scene was to be analyzed
		
				lastScene = scene;
				lastShot = shot;
				oldGtSegStart = gtSegStart;

			} //speech-segments was long enough to get evaluated

			y = gtSegEnd;
		} else {
			break;
		}
	}

	for (ce = 0; ce < 4; ce++) {
		//calc scores for each class independantly
		calcSupervisedClassificationScores(this->scatterMatrix[ce], this->ca, this->averageOmission[ce], this->averageCommission[ce], this->missclassificationRate[ce], this->kappaStatistic[ce], this->recall[ce], this->precision[ce], this->missRate[ce], this->falseAlarmRate[ce], this->errorRate[ce], this->specificity[ce], this->accuracy[ce], 2);
	}
	this->DER = calcDER(softBoundaryDiameter); //uses the groundtruth object instead of the scatter matrix to calculate the score because non-speech segments are not included in the scatter matrixes but are scored in DER

	MFree_2D(hitsInScene);
	MFree_2D(hitsInShot);
	MFree_0D(pFunc);

  return;
}

//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_SpeakerScore_Classification::output(ostream& OutS) {
	unsigned long int x, y, ce;
	char buf[sclib::bufferSize];
	
	//overall classification
	OutS << "Classification result:\n";
	OutS << "-----------------------\n";
	OutS << "Scatter-matrix:\n";
	for (ce = 0; ce < 4; ce++) { //loop over all CEs
		OutS << "\t";
		for (x = 0; x < this->ca; x++) { //headline
			idx2class(x, buf, true); //get the class' short name
			OutS << "\t" << buf;
		}
		ce2text(ce, buf);
		OutS << " (per " << buf << ")\n";
		for (y = 0; y < this->ca; y++) { //table
			idx2class(y, buf, true); //get the class' short name
			OutS << "\t" << buf << "\t";
			for (x = 0; x < this->ca; x++) {
				OutS << this->scatterMatrix[ce][y][x] << "\t";
			}
			OutS << "\n";
		}
		OutS << "\n";
	}
	OutS << "                          Sample-\tSeg.-\tShot-\tScene-based\n";
	OutS << "Average omission:         " << setprecision(4) << this->averageOmission[sclib::ceSample] << "\t" << this->averageOmission[sclib::ceSegment] << "\t" << this->averageOmission[sclib::ceShot] << "\t" << this->averageOmission[sclib::ceScene] << "\n";
	OutS << "Average commission:       " << setprecision(4) << this->averageCommission[sclib::ceSample] << "\t" << this->averageCommission[sclib::ceSegment] << "\t" << this->averageCommission[sclib::ceShot] << "\t" << this->averageCommission[sclib::ceScene] << "\n";
	OutS << "Missclassification rate:  " << setprecision(4) << this->missclassificationRate[sclib::ceSample] << "\t" << this->missclassificationRate[sclib::ceSegment] << "\t" << this->missclassificationRate[sclib::ceShot] << "\t" << this->missclassificationRate[sclib::ceScene] << "\n";
	OutS << "Kappa statistic:          " << setprecision(4) << this->kappaStatistic[sclib::ceSample] << "\t" << this->kappaStatistic[sclib::ceSegment] << "\t" << this->kappaStatistic[sclib::ceShot] << "\t" << this->kappaStatistic[sclib::ceScene] << "\n";
	OutS << "K-hat index:              see kappa statistic\n";
	OutS << "Diarization Error Rate:   " << setprecision(4) << this->DER << "\n";

	//detection singleclass vs. rest
	OutS << "\n";
	OutS << "Detection results for each class vs. rest:\n";
	OutS << "-------------------------------------------\n";
	OutS << "                          Sample-\tSeg.-\tShot-\tScene-based\n";
	for (x = 0; x < this->ca; x++) { //loop over all classes
		this->idx2class(x, buf, false); //get the class' long name
		OutS << "-> " << buf << " vs. rest:\n";
		OutS << "    Recall:              " << setprecision(4) << this->recall[sclib::ceSample][x] << "\t" << this->recall[sclib::ceSegment][x] << "\t" << this->recall[sclib::ceShot][x] << "\t" << this->recall[sclib::ceScene][x] << "\n";
		OutS << "    Precision:           " << setprecision(4) << this->precision[sclib::ceSample][x] << "\t" << this->precision[sclib::ceSegment][x] << "\t" << this->precision[sclib::ceShot][x] << "\t" << this->precision[sclib::ceScene][x] << "\n";
		OutS << "    Specificity:         " << setprecision(4) << this->specificity[sclib::ceSample][x] << "\t" << this->specificity[sclib::ceSegment][x] << "\t" << this->specificity[sclib::ceShot][x] << "\t" << this->specificity[sclib::ceScene][x] << "\n";
		OutS << "    Omission:            " << setprecision(4) << 1.0-this->recall[sclib::ceSample][x] << "\t" << 1.0-this->recall[sclib::ceSegment][x] << "\t" << 1.0-this->recall[sclib::ceShot][x] << "\t" << 1.0-this->recall[sclib::ceScene][x] << "\n";
		OutS << "    Commission:          " << setprecision(4) << 1.0-this->precision[sclib::ceSample][x] << "\t" << 1.0-this->precision[sclib::ceSegment][x] << "\t" << 1.0-this->precision[sclib::ceShot][x] << "\t" << 1.0-this->precision[sclib::ceScene][x] << "\n";
		OutS << "    Miss rate:           " << setprecision(4) << this->missRate[sclib::ceSample][x] << "\t" << this->missRate[sclib::ceSegment][x] << "\t" << this->missRate[sclib::ceShot][x] << "\t" << this->missRate[sclib::ceScene][x] << "\n";
		OutS << "    False alarm rate:    " << setprecision(4) << this->falseAlarmRate[sclib::ceSample][x] << "\t" << this->falseAlarmRate[sclib::ceSegment][x] << "\t" << this->falseAlarmRate[sclib::ceShot][x] << "\t" << this->falseAlarmRate[sclib::ceScene][x] << "\n";
		OutS << "    Error rate:          " << setprecision(4) << this->errorRate[sclib::ceSample][x] << "\t" << this->errorRate[sclib::ceSegment][x] << "\t" << this->errorRate[sclib::ceShot][x] << "\t" << this->errorRate[sclib::ceScene][x] << "\n";
		OutS << "    Accuracy:            " << setprecision(4) << this->accuracy[sclib::ceSample][x] << "\t" << this->accuracy[sclib::ceSegment][x] << "\t" << this->accuracy[sclib::ceShot][x] << "\t" << this->accuracy[sclib::ceScene][x] << "\n";
		OutS << "    Fidelity:            " << setprecision(4) << 1.0-this->errorRate[sclib::ceSample][x] << "\t" << 1.0-this->errorRate[sclib::ceSegment][x] << "\t" << 1.0-this->errorRate[sclib::ceShot][x] << "\t" << 1.0-this->errorRate[sclib::ceScene][x] << "\n";
		OutS << "    Sensitivity:         see recall\n";
		OutS << "    Producer's accuracy: see omission\n";
		OutS << "    User's accuracy:     see commission\n";
		OutS << "    True positive rate:  see recall\n";
		OutS << "    False negative rate: see miss rate\n";
		OutS << "    False positive rate: see false alarm rate\n";
	}

  return OutS;
}
