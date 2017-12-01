/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the audio-type    */
/*        classification process																					*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 27.09.2006																								*/
/**************************************************************************/

#include "SC_Score_AudioTypeClassification.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Score_AudioTypeClassification::SC_Score_AudioTypeClassification(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT) : SC_Score(pTweak, pGT) {
	this->ca = 8;
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

		this->sn_scatterMatrix[ce] = NULL;
		this->sn_averageCommission[ce] = 0.0;
		this->sn_averageOmission[ce] = 0.0;
		this->sn_kappaStatistic[ce] = 0.0;
		this->sn_missclassificationRate[ce] = 0.0;
	}
	sprintf(this->purpose, "%s\0", "Audio-Type Classification Assessment");
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_Score_AudioTypeClassification::~SC_Score_AudioTypeClassification() {
	for (unsigned int ce = 0; ce < 4; ce++) {
		MFree_2D(this->scatterMatrix[ce]);
		MFree_1D(this->recall[ce]);
		MFree_1D(this->precision[ce]);
		MFree_1D(this->missRate[ce]);
		MFree_1D(this->falseAlarmRate[ce]);
		MFree_1D(this->errorRate[ce]);
		MFree_1D(this->specificity[ce]);
		MFree_1D(this->accuracy[ce]);

		MFree_2D(this->sn_scatterMatrix[ce]);
	}
}

//====================================================================================================================
//	Convert class-tags as used in the ground-truth (sclib::atSpeech etc.) to indices into the scatter-matrix/result-
//  vectors; SVLIB_Fail is returned if the mapping can't be established
//====================================================================================================================
long int SC_Score_AudioTypeClassification::class2idx(unsigned long int classTag) {
	long int res = SVLIB_Fail;

	switch (classTag) {
		case sclib::atSilence:
			res = 0;
			break;
		case sclib::atPureSpeech:
			res = 1;
			break;
		case sclib::atNoisySpeech:
			res = 2;
			break;
		case sclib::atBackground:
			res = 3;
			break;
		case sclib::atMusic:
			res = 4;
			break;
		case sclib::atAction:
			res = 5;
			break;
		case sclib::atBreath:
			res = 6;
			break;
		case sclib::atUndefined:
			res = 7;
			break;
		case sclib::atSpeech:
			res = 8;
			break;
		case sclib::atNoise:
			res = 9;
			break;
	}

	return res;
}

//====================================================================================================================
//	Convert class indices into the scatter-matrix/result-vectors to class-tags (sclib::atSpeech etc.) as used in the 
//  ground-truth; SVLIB_Fail is returned if the mapping can't be established; if className is != NULL, a string 
//  giving the name of the class (verbose or short according to shortName) is also returned in that variable
//====================================================================================================================
long int SC_Score_AudioTypeClassification::idx2class(unsigned long int classIdx, char *className, bool shortName) {
	long int res = SVLIB_Fail;

	switch (classIdx) {
		case 0:
			res = sclib::atSilence;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "Silnc");
				} else {
					sprintf(className, "%s\0", "Silence");
				}
			}
			break;
		case 1:
			res = sclib::atPureSpeech;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "PSpch");
				} else {
					sprintf(className, "%s\0", "Pure Speech");
				}
			}
			break;
		case 2:
			res = sclib::atNoisySpeech;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "NSpch");
				} else {
					sprintf(className, "%s\0", "Noisy Speech");
				}
			}
			break;
		case 3:
			res = sclib::atBackground;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "BG");
				} else {
					sprintf(className, "%s\0", "Background");
				}
			}
			break;
		case 4:
			res = sclib::atMusic;
			if (className != NULL) {
				sprintf(className, "%s\0", "Music");
			}
			break;
		case 5:
			res = sclib::atAction;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "Actn");
				} else {
					sprintf(className, "%s\0", "Action");
				}
			}
			break;
		case 6:
			res = sclib::atBreath;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "Brth");
				} else {
					sprintf(className, "%s\0", "Breathing");
				}
			}
			break;
		case 7:
			res = sclib::atUndefined;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "Undef");
				} else {
					sprintf(className, "%s\0", "Undefined");
				}
			}
			break;
		case 8:
			res = sclib::atSpeech;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "Spch");
				} else {
					sprintf(className, "%s\0", "Speech");
				}
			}
			break;
		case 9:
			res = sclib::atNoise;
			if (className != NULL) {
				sprintf(className, "%s\0", "Noise");
			}
			break;
	}

	return res;
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
void SC_Score_AudioTypeClassification::calcScores(unsigned long int start, unsigned long int end, unsigned long int algorithmicUncertaintyDiameter) {
	unsigned int ce, x, y, types = 0;
	unsigned long int gtType, hypoType, softBoundaryDiameter = this->pGT->getUncertaintyRegionWidth(true) + algorithmicUncertaintyDiameter;
	unsigned long int idxX, idxY, scene = 0, lastScene = 0, shot = 0, lastShot = 0;
	unsigned long int idxPureSpeech = class2idx(sclib::atPureSpeech), idxNoisySpeech = class2idx(sclib::atNoisySpeech);
	bool **hitsInScene, **hitsInShot;
	long int **detectionScatterMatrix = NULL, gtSegStart, gtSegEnd, hypoSegStart, hypoSegEnd, oldGtSegEnd = sclib::noSegment, oldGtSegStart = sclib::noSegment;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

	this->start = start;
	this->end = (end > 0 && end < this->pGT->getAudioSampleCount()) ? end : this->pGT->getAudioSampleCount()-1;
	
	//get all types for which this class is responsible (to minimize number of hardcoded lists)
	for (x = 0; x < this->ca; x++) {
		types |= idx2class(x);
	}

	for (ce = 0; ce < 4; ce++) {
		MFree_2D(this->scatterMatrix[ce]);
		this->scatterMatrix[ce] = pFunc->initMatrix(this->ca+2, this->ca, (long int)(0));
	}

	MArray_2D(hitsInScene, (long int)(this->ca)+2, this->ca, bool, "SC_Score_AudioTypeClassification.calcScores: hitsInScene");
	MArray_2D(hitsInShot, (long int)(this->ca)+2, this->ca, bool, "SC_Score_AudioTypeClassification.calcScores: hitsInShot");
	for (y = 0; y < this->ca+2; y++) {
		for (x = 0; x < this->ca; x++) {
			hitsInScene[y][x] = false;
			hitsInShot[y][x] = false;
		}
	}

	//basically, partition all available samples into gt-segments and match (parts of) hypo-segments to them
	for (y = this->start; y <= this->end; y++) {
		//get continuous segments of only one of the desired types out of the groundtruth
		//ATTENTION: it is assumed that there is no gap between those segments, so that all false positive values are found 
		//           as false negatives in other gt-based segments (and we don't have to look hypo-based into into the gaps...)
		gtType = this->pGT->getClosestSegment(y, gtSegStart, gtSegEnd, types, false, sclib::searchForward, sclib::modeGroundtruth, false, true, sclib::atShotBoundary);
		if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment && gtSegStart < (long)(this->end)) {
			idxX = class2idx(gtType);

			//check for scene- and shot-changes
			scene = this->pGT->sample2scene(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastScene);
			if (scene > this->pTweak->general.lastScene) {
				break;
			}
			if (scene != lastScene) {
				pFunc->clear(hitsInScene, this->ca+2, this->ca, false); //clear rememberance of hits
			}
			shot = this->pGT->sample2shot(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastShot);
			if (shot != lastShot) {
				pFunc->clear(hitsInShot, this->ca+2, this->ca, false); //clear rememberance of hits
			}

			//has this scene been processed by any algorithm?
			if (scene >= this->pTweak->general.firstScene && (this->pTweak->general.sceneSelection == 0 || sclib::bitTest(sclib::bit(scene), pTweak->general.sceneSelection) == true) && scene <= this->pTweak->general.lastScene) {

				//fill population- and possible set (they are the same for the ATC task)
				this->scatterMatrix[sclib::ceSample][this->ca][idxX] += gtSegEnd - gtSegStart + 1; //population 
				this->scatterMatrix[sclib::ceSegment][this->ca][idxX] += 1;
				if (hitsInShot[this->ca][idxX] == false) {
					this->scatterMatrix[sclib::ceShot][this->ca][idxX] += 1;
					hitsInShot[this->ca][idxX] = true;
				}
				if (hitsInScene[this->ca][idxX] == false) {
					this->scatterMatrix[sclib::ceScene][this->ca][idxX] += 1;
					hitsInScene[this->ca][idxX] = true;
				}
				this->scatterMatrix[sclib::ceSample][this->ca+1][idxX] += gtSegEnd - gtSegStart + 1; //possible set
				this->scatterMatrix[sclib::ceSegment][this->ca+1][idxX] += 1;
				if (hitsInShot[this->ca+1][idxX] == false) {
					this->scatterMatrix[sclib::ceShot][this->ca+1][idxX] += 1;
					hitsInShot[this->ca+1][idxX] = true;
				}
				if (hitsInScene[this->ca+1][idxX] == false) {
					this->scatterMatrix[sclib::ceScene][this->ca+1][idxX] += 1;
					hitsInScene[this->ca+1][idxX] = true;
				}

				if (oldGtSegEnd == sclib::noSegment || gtSegStart == oldGtSegEnd+1) { //check no-gap assumption
					for (x = (unsigned long int)(gtSegStart); x <= (unsigned long int)(gtSegEnd); x++) {
						//SEARCH_WITHIN is used so that segments are found even when they don't match the exact gt-boundaries
						//it is NOT tried to find the best fitting segment if boudaries don't match, because this is not needed 
						//(for samples, it is automatically cared for, and for segments, all is lost if the boundaries don't match anyway)
						hypoType = this->pGT->getClosestSegment(x, hypoSegStart, hypoSegEnd, types, false, sclib::searchWithin, sclib::modeHypothesized, false, true, sclib::atShotBoundary);
						if (hypoSegStart != sclib::noSegment && hypoSegEnd != sclib::noSegment && hypoSegStart < gtSegEnd) {
							idxY = class2idx(hypoType);
		
							//care for soft boundaries near to the groundtruth on-/offset
							if (abs(hypoSegStart - gtSegStart) <= (long)(softBoundaryDiameter)) {
								hypoSegStart = gtSegStart;
							}
							if (abs(hypoSegEnd - gtSegEnd) <= (long)(softBoundaryDiameter)) {
								hypoSegEnd = gtSegEnd;
							}

							//CE==sample
							this->scatterMatrix[sclib::ceSample][idxY][idxX] += min(hypoSegEnd, gtSegEnd) - hypoSegStart + 1;
							
							//CE==segment
							if (hypoSegStart == gtSegStart && hypoSegEnd == gtSegEnd) {
								this->scatterMatrix[sclib::ceSegment][idxY][idxX] += 1;
							}

							//CE==shot
							if (hitsInShot[idxY][idxX] == false) {
								this->scatterMatrix[sclib::ceShot][idxY][idxX] += 1;
								hitsInShot[idxY][idxX] = true; //remember that we already counted a hit in this shot for these two types
							}

							//CE==scene
							if (hitsInScene[idxY][idxX] == false) {
								this->scatterMatrix[sclib::ceScene][idxY][idxX] += 1;
								hitsInScene[idxY][idxX] = true; //remember that we already counted a hit in this scene for these two types
							}

							//truncate too big hypo-segments at the end of the gt-segment
							if (hypoSegEnd > gtSegEnd) {
								hypoSegEnd = gtSegEnd;
							}

							x = hypoSegEnd;
						} else { //still space left in the gt-segment, but no further hypo-segments
							if (gtSegEnd-x > 0) { 
								//no more hypo-segments, so possibly remaining unmatched samples are counted as being labeled UNDEF in hypo
								idxY = class2idx(sclib::atUndefined);
								this->scatterMatrix[sclib::ceSample][idxY][idxX] += gtSegEnd - x; 
								this->scatterMatrix[sclib::ceSegment][idxY][idxX] += 1;
								if (hitsInShot[idxY][idxX] == false) {
									this->scatterMatrix[sclib::ceShot][idxY][idxX] += 1;
									hitsInShot[idxY][idxX] = true;
								}
								if (hitsInScene[idxY][idxX] == false) {
									this->scatterMatrix[sclib::ceScene][idxY][idxX] += 1;
									hitsInScene[idxY][idxX] = true;
								}
							}
							break;
						} //valid hypo-segment found (or not)
					} //find hypo-segments for one given gt-segment
				} else { //there was a gap between the previous and this segment 
					REPORT_ERROR(SVLIB_BadData, "Gap found between consecutive gt-segments during ATC-scoring - this must not be!");
				}
			} //scene was not to evaluate

			lastScene = scene;
			lastShot = shot;
      oldGtSegEnd = gtSegEnd;
			oldGtSegStart = gtSegStart;
			y = gtSegEnd;
		} else {
			break;
		}
	}

	MFree_2D(hitsInShot);
	MFree_2D(hitsInScene);

	for (ce = 0; ce < 4; ce++) {
		//calc scores for each class independantly
		calcSupervisedClassificationScores(this->scatterMatrix[ce], this->ca, this->averageOmission[ce], this->averageCommission[ce], this->missclassificationRate[ce], this->kappaStatistic[ce], this->recall[ce], this->precision[ce], this->missRate[ce], this->falseAlarmRate[ce], this->errorRate[ce], this->specificity[ce], this->accuracy[ce], 2);

		//create additional entrys for speech/non-speech discrimination
		detectionScatterMatrix = pFunc->initMatrix(4, 2, (long int)(0));
		for (y = 0; y < this->ca+2; y++) { //scatter matrix with speech=P, non-speech=N
			for (x = 0; x < this->ca; x++) {
				if (y < this->ca) { //TP/TN/FP/FN
					if (y == idxPureSpeech || y == idxNoisySpeech) { //speech rows
						if (x == idxPureSpeech || x == idxNoisySpeech) { //speech cols
							detectionScatterMatrix[0][0] += this->scatterMatrix[ce][y][x]; //TPs
						} else { //noise cols
							detectionScatterMatrix[0][1] += this->scatterMatrix[ce][y][x]; //FPs
						}
					} else { //noise rows
						if (x == idxPureSpeech || x == idxNoisySpeech) { //speech cols
							detectionScatterMatrix[1][0] += this->scatterMatrix[ce][y][x]; //FNs
						} else { //noise cols
							detectionScatterMatrix[1][1] += this->scatterMatrix[ce][y][x]; //TNs
						}
					}
				} else { //population/possible
					if (x == idxPureSpeech || x == idxNoisySpeech) { //speech cols
						detectionScatterMatrix[2+y-this->ca][0] += this->scatterMatrix[ce][y][x]; //TPs
					} else { //noise cols
						detectionScatterMatrix[2+y-this->ca][1] += this->scatterMatrix[ce][y][x]; //FPs
					}
				}
			}
		}
		idxX = class2idx(sclib::atSpeech);
		calcDetectionScores(detectionScatterMatrix, this->recall[ce][idxX], this->precision[ce][idxX], this->missRate[ce][idxX], this->falseAlarmRate[ce][idxX], this->errorRate[ce][idxX], this->specificity[ce][idxX], this->accuracy[ce][idxX]);
		
		MFree_2D(this->sn_scatterMatrix[ce]);
		this->sn_scatterMatrix[ce] = detectionScatterMatrix;
		calcSupervisedClassificationScores(this->sn_scatterMatrix[ce], 2, this->sn_averageOmission[ce], this->sn_averageCommission[ce], this->sn_missclassificationRate[ce], this->sn_kappaStatistic[ce]);
		
		//scatter matrix with speech=N, non-speech=P
		detectionScatterMatrix = pFunc->copy(this->sn_scatterMatrix[ce], 4, 2);
		flipDetectionScatterMatrix(detectionScatterMatrix);
		idxX = class2idx(sclib::atNoise);
		calcDetectionScores(detectionScatterMatrix, this->recall[ce][idxX], this->precision[ce][idxX], this->missRate[ce][idxX], this->falseAlarmRate[ce][idxX], this->errorRate[ce][idxX], this->specificity[ce][idxX], this->accuracy[ce][idxX]);

		MFree_2D(detectionScatterMatrix);
	}

	MFree_0D(pFunc);

  return;
}

//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_Score_AudioTypeClassification::output(ostream& OutS) {
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
	OutS << "                  \t\tSample-\tSeg.-\tShot-\tScene-based\n";
	OutS << "Average omission: \t\t" << setprecision(4) << this->averageOmission[sclib::ceSample] << "\t" << this->averageOmission[sclib::ceSegment] << "\t" << this->averageOmission[sclib::ceShot] << "\t" << this->averageOmission[sclib::ceScene] << "\n";
	OutS << "Average commission: \t\t" << setprecision(4) << this->averageCommission[sclib::ceSample] << "\t" << this->averageCommission[sclib::ceSegment] << "\t" << this->averageCommission[sclib::ceShot] << "\t" << this->averageCommission[sclib::ceScene] << "\n";
	OutS << "Missclassification rate: \t" << setprecision(4) << this->missclassificationRate[sclib::ceSample] << "\t" << this->missclassificationRate[sclib::ceSegment] << "\t" << this->missclassificationRate[sclib::ceShot] << "\t" << this->missclassificationRate[sclib::ceScene] << "\n";
	OutS << "Kappa statistic:  \t\t" << setprecision(4) << this->kappaStatistic[sclib::ceSample] << "\t" << this->kappaStatistic[sclib::ceSegment] << "\t" << this->kappaStatistic[sclib::ceShot] << "\t" << this->kappaStatistic[sclib::ceScene] << "\n";
	OutS << "K-hat index:      \t\tsee kappa statistic\n";

	//detection singleclass vs. rest
	OutS << "\n";
	OutS << "Detection results for each class vs. rest:\n";
	OutS << "-------------------------------------------\n";
	OutS << "             \t\t\tSample-\tSeg.-\tShot-\tScene-based\n";
	for (x = 0; x < this->ca; x++) { //loop over all classes
		this->idx2class(x, buf, false); //get the class' long name
		OutS << "-> " << buf << " vs. rest:\n";
		OutS << "    Recall: \t\t\t" << setprecision(4) << this->recall[sclib::ceSample][x] << "\t" << this->recall[sclib::ceSegment][x] << "\t" << this->recall[sclib::ceShot][x] << "\t" << this->recall[sclib::ceScene][x] << "\n";
		OutS << "    Precision: \t\t\t" << setprecision(4) << this->precision[sclib::ceSample][x] << "\t" << this->precision[sclib::ceSegment][x] << "\t" << this->precision[sclib::ceShot][x] << "\t" << this->precision[sclib::ceScene][x] << "\n";
		OutS << "    Specificity: \t\t" << setprecision(4) << this->specificity[sclib::ceSample][x] << "\t" << this->specificity[sclib::ceSegment][x] << "\t" << this->specificity[sclib::ceShot][x] << "\t" << this->specificity[sclib::ceScene][x] << "\n";
		OutS << "    Omission: \t\t\t" << setprecision(4) << 1.0-this->recall[sclib::ceSample][x] << "\t" << 1.0-this->recall[sclib::ceSegment][x] << "\t" << 1.0-this->recall[sclib::ceShot][x] << "\t" << 1.0-this->recall[sclib::ceScene][x] << "\n";
		OutS << "    Commission: \t\t" << setprecision(4) << 1.0-this->precision[sclib::ceSample][x] << "\t" << 1.0-this->precision[sclib::ceSegment][x] << "\t" << 1.0-this->precision[sclib::ceShot][x] << "\t" << 1.0-this->precision[sclib::ceScene][x] << "\n";
		OutS << "    Miss rate: \t\t\t" << setprecision(4) << this->missRate[sclib::ceSample][x] << "\t" << this->missRate[sclib::ceSegment][x] << "\t" << this->missRate[sclib::ceShot][x] << "\t" << this->missRate[sclib::ceScene][x] << "\n";
		OutS << "    False alarm rate: \t\t" << setprecision(4) << this->falseAlarmRate[sclib::ceSample][x] << "\t" << this->falseAlarmRate[sclib::ceSegment][x] << "\t" << this->falseAlarmRate[sclib::ceShot][x] << "\t" << this->falseAlarmRate[sclib::ceScene][x] << "\n";
		OutS << "    Error rate: \t\t" << setprecision(4) << this->errorRate[sclib::ceSample][x] << "\t" << this->errorRate[sclib::ceSegment][x] << "\t" << this->errorRate[sclib::ceShot][x] << "\t" << this->errorRate[sclib::ceScene][x] << "\n";
		OutS << "    Accuracy: \t\t\t" << setprecision(4) << this->accuracy[sclib::ceSample][x] << "\t" << this->accuracy[sclib::ceSegment][x] << "\t" << this->accuracy[sclib::ceShot][x] << "\t" << this->accuracy[sclib::ceScene][x] << "\n";
		OutS << "    Fidelity: \t\t\t" << setprecision(4) << 1.0-this->errorRate[sclib::ceSample][x] << "\t" << 1.0-this->errorRate[sclib::ceSegment][x] << "\t" << 1.0-this->errorRate[sclib::ceShot][x] << "\t" << 1.0-this->errorRate[sclib::ceScene][x] << "\n";
		OutS << "    Sensitivity: \t\tsee recall\n";
		OutS << "    Producer's accuracy: \tsee omission\n";
		OutS << "    User's accuracy: \t\tsee commission\n";
		OutS << "    True positive rate: \tsee recall\n";
		OutS << "    False negative rate: \tsee miss rate\n";
		OutS << "    False positive rate: \tsee false alarm rate\n";
	}

	//speech/non-speech
	OutS << "\n";
	OutS << "Speech vs. noise detection result:\n";
	OutS << "-----------------------------------\n";
	OutS << "Scatter-matrix:\n";
	for (ce = 0; ce < 4; ce++) { //loop over all CEs
		OutS << "\t";
		for (x = this->ca; x < this->ca+2; x++) { //headline
			idx2class(x, buf, true); //get the class' short name
			OutS << "\t" << buf;
		}
		ce2text(ce, buf);
		OutS << " (per " << buf << ")\n";
		for (y = 0; y < 2; y++) { //table
			//if (y < 2) {
				idx2class(y+this->ca, buf, true); //get the class' short name
				OutS << "\t" << buf << "\t";
			//} else {
			//	OutS << "\t\t";
			//}
			for (x = 0; x < 2; x++) {
				OutS << this->sn_scatterMatrix[ce][y][x] << "\t";
			}
			OutS << "\n";
		}
		OutS << "\n";
	}
	OutS << "                  \t\tSample-\tSeg.-\tShot-\tScene-based\n";
	OutS << "Average omission: \t\t" << setprecision(4) << this->averageOmission[sclib::ceSample] << "\t" << this->averageOmission[sclib::ceSegment] << "\t" << this->averageOmission[sclib::ceShot] << "\t" << this->averageOmission[sclib::ceScene] << "\n";
	OutS << "Average commission: \t\t" << setprecision(4) << this->averageCommission[sclib::ceSample] << "\t" << this->averageCommission[sclib::ceSegment] << "\t" << this->averageCommission[sclib::ceShot] << "\t" << this->averageCommission[sclib::ceScene] << "\n";
	OutS << "Missclassification rate: \t" << setprecision(4) << this->missclassificationRate[sclib::ceSample] << "\t" << this->missclassificationRate[sclib::ceSegment] << "\t" << this->missclassificationRate[sclib::ceShot] << "\t" << this->missclassificationRate[sclib::ceScene] << "\n";
	OutS << "Kappa statistic:  \t\t" << setprecision(4) << this->kappaStatistic[sclib::ceSample] << "\t" << this->kappaStatistic[sclib::ceSegment] << "\t" << this->kappaStatistic[sclib::ceShot] << "\t" << this->kappaStatistic[sclib::ceScene] << "\n";
	OutS << "K-hat index:      \t\tsee kappa statistic\n\n";
	
	for (x = this->ca; x < this->ca+2; x++) { //loop over all classes
		this->idx2class(x, buf, false); //get the class' long name
		OutS << "-> " << buf << "-detection:\n";
		OutS << "    Recall: \t\t\t" << setprecision(4) << this->recall[sclib::ceSample][x] << "\t" << this->recall[sclib::ceSegment][x] << "\t" << this->recall[sclib::ceShot][x] << "\t" << this->recall[sclib::ceScene][x] << "\n";
		OutS << "    Precision: \t\t\t" << setprecision(4) << this->precision[sclib::ceSample][x] << "\t" << this->precision[sclib::ceSegment][x] << "\t" << this->precision[sclib::ceShot][x] << "\t" << this->precision[sclib::ceScene][x] << "\n";
		OutS << "    Specificity: \t\t" << setprecision(4) << this->specificity[sclib::ceSample][x] << "\t" << this->specificity[sclib::ceSegment][x] << "\t" << this->specificity[sclib::ceShot][x] << "\t" << this->specificity[sclib::ceScene][x] << "\n";
		OutS << "    Omission: \t\t\t" << setprecision(4) << 1.0-this->recall[sclib::ceSample][x] << "\t" << 1.0-this->recall[sclib::ceSegment][x] << "\t" << 1.0-this->recall[sclib::ceShot][x] << "\t" << 1.0-this->recall[sclib::ceScene][x] << "\n";
		OutS << "    Commission: \t\t" << setprecision(4) << 1.0-this->precision[sclib::ceSample][x] << "\t" << 1.0-this->precision[sclib::ceSegment][x] << "\t" << 1.0-this->precision[sclib::ceShot][x] << "\t" << 1.0-this->precision[sclib::ceScene][x] << "\n";
		OutS << "    Miss rate: \t\t\t" << setprecision(4) << this->missRate[sclib::ceSample][x] << "\t" << this->missRate[sclib::ceSegment][x] << "\t" << this->missRate[sclib::ceShot][x] << "\t" << this->missRate[sclib::ceScene][x] << "\n";
		OutS << "    False alarm rate: \t\t" << setprecision(4) << this->falseAlarmRate[sclib::ceSample][x] << "\t" << this->falseAlarmRate[sclib::ceSegment][x] << "\t" << this->falseAlarmRate[sclib::ceShot][x] << "\t" << this->falseAlarmRate[sclib::ceScene][x] << "\n";
		OutS << "    Error rate: \t\t" << setprecision(4) << this->errorRate[sclib::ceSample][x] << "\t" << this->errorRate[sclib::ceSegment][x] << "\t" << this->errorRate[sclib::ceShot][x] << "\t" << this->errorRate[sclib::ceScene][x] << "\n";
		OutS << "    Accuracy: \t\t\t" << setprecision(4) << this->accuracy[sclib::ceSample][x] << "\t" << this->accuracy[sclib::ceSegment][x] << "\t" << this->accuracy[sclib::ceShot][x] << "\t" << this->accuracy[sclib::ceScene][x] << "\n";
		OutS << "    Fidelity: \t\t\t" << setprecision(4) << 1.0-this->errorRate[sclib::ceSample][x] << "\t" << 1.0-this->errorRate[sclib::ceSegment][x] << "\t" << 1.0-this->errorRate[sclib::ceShot][x] << "\t" << 1.0-this->errorRate[sclib::ceScene][x] << "\n";
		OutS << "    Sensitivity: \t\tsee recall\n";
		OutS << "    Producer's accuracy: \tsee omission\n";
		OutS << "    User's accuracy: \t\tsee commission\n";
		OutS << "    True positive rate: \tsee recall\n";
		OutS << "    False negative rate: \tsee miss rate\n";
		OutS << "    False positive rate: \tsee false alarm rate\n";
	}

  return OutS;
}
