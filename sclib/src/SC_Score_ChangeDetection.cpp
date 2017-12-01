/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the change        */
/*        detction (speaker and acoustic) process													*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.11.2006																								*/
/**************************************************************************/

#include "SC_Score_ChangeDetection.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Score_ChangeDetection::SC_Score_ChangeDetection(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT) : SC_Score(pTweak, pGT) {
	unsigned int x, y;
	
	this->ca = 2;
	for (x = 0; x < 4; x++) {
		MArray_1D(this->scatterMatrix[x], this->ca, long int**, "SC_Score_ChangeDetection: scatterMatrix[x]");
		for (y = 0; y < this->ca; y++) {
			this->scatterMatrix[x][y] = NULL;
		}
		MArray_1D(this->recall[x], this->ca, double, "SC_Score_ChangeDetection: recall[x]");
		MArray_1D(this->precision[x], this->ca, double, "SC_Score_ChangeDetection: precision[x]");
		MArray_1D(this->missRate[x], this->ca, double, "SC_Score_ChangeDetection: missRate[x]");
		MArray_1D(this->falseAlarmRate[x], this->ca, double, "SC_Score_ChangeDetection: falseAlarmRate[x]");
		MArray_1D(this->errorRate[x], this->ca, double, "SC_Score_ChangeDetection: errorRate[x]");
		MArray_1D(this->specificity[x], this->ca, double, "SC_Score_ChangeDetection: specificity[x]");
		MArray_1D(this->accuracy[x], this->ca, double, "SC_Score_ChangeDetection: accuracy[x]");
		for (y = 0; y < this->ca; y++) {
			this->recall[x][y] = 0.0;
			this->precision[x][y] = 0.0;
			this->missRate[x][y] = 0.0;
			this->falseAlarmRate[x][y] = 0.0;
			this->errorRate[x][y] = 0.0;
			this->specificity[x][y] = 0.0;
			this->accuracy[x][y] = 0.0;
		}
	}
	for (x = 0; x < this->ca; x++) {
		this->gtBoundaryList[x] = NULL;
		this->hypoFNlist[x] = NULL;
		this->hypoFPlist[x] = NULL;
		this->hypoTPlist[x] = NULL;
	}
	sprintf(this->purpose, "%s\0", "Change Detection Assessment");
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_Score_ChangeDetection::~SC_Score_ChangeDetection() {
	unsigned int x, y;

	for (x = 0; x < 4; x++) {
		for (y = 0; y < this->ca; y++) {
			MFree_2D(this->scatterMatrix[x][y]);
		}
		MFree_1D(this->scatterMatrix[x]);
		MFree_1D(this->recall[x]);
		MFree_1D(this->precision[x]);
		MFree_1D(this->missRate[x]);
		MFree_1D(this->falseAlarmRate[x]);
		MFree_1D(this->errorRate[x]);
		MFree_1D(this->specificity[x]);
		MFree_1D(this->accuracy[x]);
	}
	for (x = 0; x < this->ca; x++) {
		sclib::destructLinkedList(this->gtBoundaryList[x]);
		sclib::destructLinkedList(this->hypoFNlist[x]);
		sclib::destructLinkedList(this->hypoFPlist[x]);
		sclib::destructLinkedList(this->hypoTPlist[x]);
	}
}

//====================================================================================================================
//	Convert class-tags as used in the ground-truth (sclib::atSpeech etc.) to indices into the scatter-matrix/result-
//  vectors; SVLIB_Fail is returned if the mapping can't be established
//====================================================================================================================
long int SC_Score_ChangeDetection::class2idx(unsigned long int classTag) {
	long int res = SVLIB_Fail;

	switch (classTag) {
		case sclib::atSpeakerBoundary:
			res = 0;
			break;
		case sclib::atNoiseBoundary:
			res = 1;
			break;
	}

	return res;
}

//====================================================================================================================
//	Convert class indices into the scatter-matrix/result-vectors to class-tags (sclib::atSpeech etc.) as used in the 
//  ground-truth; SVLIB_Fail is returned if the mapping can't be established; if className is != NULL, a string 
//  giving the name of the class (verbose or short according to shortName) is also returned in that variable
//====================================================================================================================
long int SC_Score_ChangeDetection::idx2class(unsigned long int classIdx, char *className, bool shortName) {
	long int res = SVLIB_Fail;

	switch (classIdx) {
		case 0:
			res = sclib::atSpeakerBoundary;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "SpkBd");
				} else {
					sprintf(className, "%s\0", "Speaker-Boundary");
				}
			}
			break;
		case 1:
			res = sclib::atNoiseBoundary;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "AcuBd");
				} else {
					sprintf(className, "%s\0", "Acoustic-Boundary");
				}
			}
			break;
	}

	return res;
}

//====================================================================================================================
//  This method returns for a given class index the types (or typesNot, according to "typesNot"-parameter) that 
//  correspond with the analyzed parts of the signal for this class
//  TODO: this should probably be placed into the algorithm's class itself and given to this class as a parameter for 
//  better decoupling...
//====================================================================================================================
long int SC_Score_ChangeDetection::getRelatedTypes(long int idx, bool typesNot) {
	long int res = sclib::noType;

	switch (idx) {
		case 0: //sclib::atSpeakerBoundary
			res = (typesNot == false) ? sclib::atSpeech : sclib::atPause|sclib::atUnvoiced;
			break;
		case 1: //sclib::atNoiseBoundary
			res = (typesNot == false) ? sclib::atNoise : sclib::noType;
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
void SC_Score_ChangeDetection::calcScores(unsigned long int start, unsigned long int end, unsigned long int algorithmicUncertaintyDiameter) {
	unsigned int c, ce, y, hits, b;
	unsigned long int softBoundaryDiameter = this->pGT->getUncertaintyRegionWidth(true) + algorithmicUncertaintyDiameter;
	unsigned long int scene, lastScene, shot, lastShot, lastShotSample;
	bool **hitsInScene, **hitsInShot, lastBoundaryValid, nextBoundaryValid, isArtificialBoundary;
	long int lastBoundary, nextBoundary, gtSegStart, gtSegEnd, sceneStart, sceneEnd, regionStart, regionEnd;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

	this->start = start;
	this->end = (end > 0 && end < this->pGT->getAudioSampleCount()) ? end : this->pGT->getAudioSampleCount()-1;
	
	//initialize scatter matrixes
	for (ce = 0; ce < 4; ce++) {
		for (c = 0; c < this->ca; c++) {
			MFree_2D(this->scatterMatrix[ce][c]);
			this->scatterMatrix[ce][c] = pFunc->initMatrix(4, 2, (long int)(0)); //standard detection scatter matrix: 2x2 + population + possible
		}
	}

	MArray_2D(hitsInScene, 4, 2, bool, "SC_Score_ChangeDetection.calcScores: hitsInScene");
	MArray_2D(hitsInShot, 4, 2, bool, "SC_Score_ChangeDetection.calcScores: hitsInShot");
	
	//basically, do the following:
	// for each boundary-type:
	// 1. init: lastBoundary = start
	// 2. nextBoundary = getNextBoundary()
	// 3. mark all hypo-boundaries betweeen lastBoundary and nextBoundary (respecting uncertaintyDiameter) as FP, all non-boundaries als TN
	// 4. mark max. one hypo-boundary within the uncertaintyDiameter around nextBoundary as a TP, all others as FP, all non-boundaries as TN; if there is none, count one FN
	// 5. lastBoundary = nextBoundary, GOTO 2 UNTIL end reached
	for (c = 0; c < this->ca; c++) {
		//1.
		nextBoundary = sclib::noSegment;
		scene = 0;
		shot = 0;
		lastShotSample = 0;
	
		for (y = this->start; y <= this->end; y++) {
			//analyse scene-wise
			this->pGT->getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward, sclib::modeGroundtruth); 
			if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment && sceneStart < (long int)(this->end)) {

				//1.
				sceneEnd = sclib::min(sceneEnd, this->end);
				lastScene = scene;
				lastBoundaryValid = false; //the lastBoundary from the last scene is valid no more in the next scene
				nextBoundaryValid = false; //because this is fed into the lastBoundaryValid variable later...
				pFunc->clear(hitsInScene, 4, 2, false); //clear rememberance of hits

				//has this scene been processed by any algorithm?
				scene = this->pGT->sample2scene(sceneStart, (y == this->start) ? 0 : y-1, lastScene);
				if (scene > this->pTweak->general.lastScene) {
					break;
				}	
				if (scene >= this->pTweak->general.firstScene && (this->pTweak->general.sceneSelection == 0 || sclib::bitTest(sclib::bit(scene), pTweak->general.sceneSelection) == true) && scene <= this->pTweak->general.lastScene) {

					//find all gt-boundaries in this scene
					for (b = (unsigned long int)sceneStart; b <= (unsigned long int)sceneEnd; b++) {
						
						//handle shot-changes //TODO: what if a shot changes in between one boundary?!?
						lastShot = shot;
						shot = this->pGT->sample2shot(b, lastShotSample, lastShot);
						lastShotSample = b;
						if (shot != lastShot) {
							pFunc->clear(hitsInShot, 4, 2, false);
						}

						//2.
						this->pGT->getNextBoundary(b, gtSegStart, gtSegEnd, idx2class(c), sclib::searchForward, sclib::modeGroundtruth);
						gtSegEnd = sclib::min(gtSegEnd, sceneEnd); //because boundaries are always set at the beginning of the new type's segment, they may lie in a new scene
						lastBoundary = nextBoundary; //5.: the nextBoundary from the last "round"
						lastBoundaryValid = nextBoundaryValid;
						if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment) {
							isArtificialBoundary = (this->pGT->testSegment(gtSegStart, gtSegStart, true, sclib::atArtificialBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) ? true : false;
						} else {
							isArtificialBoundary = false;
						}
						if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment && gtSegStart <= sceneEnd && isArtificialBoundary == false) { //also ignore artificial changes
							nextBoundaryValid = true;
							nextBoundary = gtSegStart;
						} else {
							nextBoundaryValid = false;
						}

						if (nextBoundaryValid == true || isArtificialBoundary == false) { //only count something if this is (a) a valid boundary or (b) there is no boundary in this scene; in other words: just skip artificial boundaries!
							addBoundary(this->gtBoundaryList[c], this->pGT->getConverter()->sample2videoFrame(nextBoundary)); //rember true boundary-positions for output

							//handle population rows
							this->scatterMatrix[sclib::ceSample][c][this->ca][0] += (nextBoundaryValid == true) ? 1 : 0; //PPop
							if (hitsInShot[2][0] == false && nextBoundaryValid == true) {
								hitsInShot[2][0] = true;
								this->scatterMatrix[sclib::ceShot][c][this->ca][0] += 1;
							}
							if (hitsInScene[2][0] == false && nextBoundaryValid == true) {
								hitsInScene[2][0] = true;
								this->scatterMatrix[sclib::ceScene][c][this->ca][0] += 1;
							}
							regionStart = (lastBoundaryValid == true) ? lastBoundary + 1 : sceneStart; //NPop
							regionEnd = (nextBoundaryValid == true) ? sclib::max(0, nextBoundary-1) : sceneEnd;
							this->scatterMatrix[sclib::ceSample][c][this->ca][1] += this->pGT->getConverter()->sample2audioFrame(regionEnd-regionStart+1-((nextBoundaryValid==true)?1:0), this->pGT->getInternalFrameSize(), this->pGT->getInternalFrameSize()); 
							if (hitsInShot[2][1] == false) { //there will always be negative population once a shot
								hitsInShot[2][1] = true;
								this->scatterMatrix[sclib::ceShot][c][this->ca][1] += 1;
							}
							if (hitsInScene[2][1] == false) { //there will always be negative population once a scene
								hitsInScene[2][1] = true;
								this->scatterMatrix[sclib::ceScene][c][this->ca][1] += 1;
							}

							//handle possible-set rows
							this->scatterMatrix[sclib::ceSample][c][this->ca+1][0] += (nextBoundaryValid == true) ? 1 : 0; //PPos
							if (hitsInShot[3][0] == false && nextBoundaryValid == true) {
								hitsInShot[3][0] = true;
								this->scatterMatrix[sclib::ceShot][c][this->ca+1][0] += 1;
							}
							if (hitsInScene[3][0] == false && nextBoundaryValid == true) {
								hitsInScene[3][0] = true;
								this->scatterMatrix[sclib::ceScene][c][this->ca+1][0] += 1;
							}
							regionStart = (lastBoundaryValid == true) ? lastBoundary + 1 : sceneStart; //NPos
							regionEnd = (nextBoundaryValid == true) ? sclib::max(0, nextBoundary-1) : sceneEnd; //TODO: hits is now really sample-based and no more a FLI-count, any problems?
							hits = this->pGT->testSegment(regionStart, regionEnd, false, getRelatedTypes(c, false), true, getRelatedTypes(c, true), false, sclib::modeGroundtruth);
							this->scatterMatrix[sclib::ceSample][c][this->ca+1][1] += hits / this->pGT->getInternalFrameSize(); //(hits-((hits>0&&nextBoundaryValid==true)?1:0))
							if (hitsInShot[3][1] == false && hits > 0) {
								hitsInShot[3][1] = true;
								this->scatterMatrix[sclib::ceShot][c][this->ca+1][1] += 1;
							}
							if (hitsInScene[3][1] == false && hits > 0) {
								hitsInScene[3][1] = true;
								this->scatterMatrix[sclib::ceScene][c][this->ca+1][1] += 1;
							}

							//3.
							regionStart = (lastBoundaryValid == true) ? lastBoundary + 1 + softBoundaryDiameter/2 : sceneStart;
							regionEnd = (nextBoundaryValid == true) ? nextBoundary - 1 - softBoundaryDiameter/2 : sceneEnd;
							hits = this->pGT->testSegment(regionStart, regionEnd, false, idx2class(c), false, sclib::atArtificialBoundary, false, sclib::modeHypothesized); //don't count artifical boundaries set by the algorithms
							if (hits > 0) { //remember the actual positions of the boundaries together with their classification (correct, wrong, ...) for output
								for (long int k = regionStart; k <= regionEnd; k++) {
									if (this->pGT->testSegment(k, k, false, idx2class(c), false, sclib::atArtificialBoundary, false, sclib::modeHypothesized) > 0) {
										addBoundary(this->hypoFPlist[c], this->pGT->getConverter()->sample2videoFrame(k));
										//if (c==0) {
										//	sclib::scalarOut("fa_s.txt", k, this->pTweak);
										//}
									}
									k += this->pGT->getInternalFrameSize() - 1;
								}
							}
							this->scatterMatrix[sclib::ceSample][c][0][1] += hits / this->pGT->getInternalFrameSize(); //FPs
							if (hitsInShot[0][1] == false && hits > 0) {
								hitsInShot[0][1] = true;
								this->scatterMatrix[sclib::ceShot][c][0][1] += 1;
							}
							if (hitsInScene[0][1] == false && hits > 0) {
								hitsInScene[0][1] = true;
								this->scatterMatrix[sclib::ceScene][c][0][1] += 1;
							}
							this->scatterMatrix[sclib::ceSample][c][1][1] += this->pGT->getConverter()->sample2audioFrame(sclib::max(regionEnd-regionStart+1, 0), this->pGT->getInternalFrameSize(), this->pGT->getInternalFrameSize()); //TNs
							if (hitsInShot[1][1] == false && (regionEnd-regionStart+1) > 0) {
								hitsInShot[1][1] = true;
								this->scatterMatrix[sclib::ceShot][c][1][1] += 1;
							}
							if (hitsInScene[1][1] == false && (regionEnd-regionStart+1) > 0) {
								hitsInScene[1][1] = true;
								this->scatterMatrix[sclib::ceScene][c][1][1] += 1;
							}

							//4.
							if (nextBoundaryValid == true) {
								regionStart = nextBoundary - softBoundaryDiameter/2;
								regionEnd = nextBoundary + softBoundaryDiameter/2;
								hits = this->pGT->testSegment(regionStart, regionEnd, false, idx2class(c), false, sclib::atArtificialBoundary, false, sclib::modeHypothesized); //don't count artifical boundaries set by the algorithms
								if (hits > 0) { //remember the actual positions of the boundaries together with their classification (correct, wrong, ...) for output
									bool found = false;
									for (long int k = regionStart; k <= regionEnd; k++) {
										if (this->pGT->testSegment(k, k, false, idx2class(c), false, sclib::atArtificialBoundary, false, sclib::modeHypothesized) > 0) {
											addBoundary((found == false)?this->hypoTPlist[c]: this->hypoFPlist[c], this->pGT->getConverter()->sample2videoFrame(k));
											//if (c==0) {
											//	sclib::scalarOut((found == false)?"tp_s.txt":"fa_s.txt", k, this->pTweak);
											//}
											found = true;
										}
										k += this->pGT->getInternalFrameSize() - 1;
									}
								} else {
									addBoundary(this->hypoFNlist[c], this->pGT->getConverter()->sample2videoFrame(nextBoundary));
									//if (c==0) {
									//	sclib::scalarOut("fn_s.txt", nextBoundary, this->pTweak);
									//}
								}
								this->scatterMatrix[sclib::ceSample][c][0][0] += (hits > 0) ? 1 : 0; //TPs
								if (hitsInShot[0][0] == false && hits > 0) {
									hitsInShot[0][0] = true;
									this->scatterMatrix[sclib::ceShot][c][0][0] += 1;
								}
								if (hitsInScene[0][0] == false && hits > 0) {
									hitsInScene[0][0] = true;
									this->scatterMatrix[sclib::ceScene][c][0][0] += 1;
								}
								this->scatterMatrix[sclib::ceSample][c][0][1] += (hits > this->pGT->getInternalFrameSize()) ? (hits/this->pGT->getInternalFrameSize())-1 : 0; //hits-1 : 0; //FPs
								if (hitsInShot[0][1] == false && hits > 1) {
									hitsInShot[0][1] = true;
									this->scatterMatrix[sclib::ceShot][c][0][1] += 1;
								}
								if (hitsInScene[0][1] == false && hits > 1) {
									hitsInScene[0][1] = true;
									this->scatterMatrix[sclib::ceScene][c][0][1] += 1;
								}
								this->scatterMatrix[sclib::ceSample][c][1][0] += (hits == 0) ? 1 : 0; //FNs
								if (hitsInShot[1][0] == false && hits == 0) {
									hitsInShot[1][0] = true;
									this->scatterMatrix[sclib::ceShot][c][1][0] += 1;
								}
								if (hitsInScene[1][0] == false && hits == 0) {
									hitsInScene[1][0] = true;
									this->scatterMatrix[sclib::ceScene][c][1][0] += 1;
								}
								this->scatterMatrix[sclib::ceSample][c][1][1] += sclib::max(this->pGT->getConverter()->sample2audioFrame(regionEnd-regionStart+1, this->pGT->getInternalFrameSize(), this->pGT->getInternalFrameSize())-(hits/this->pGT->getInternalFrameSize()), 0); //sclib::max(regionEnd-regionStart+1 - ((hits > this->pGT->getInternalFrameSize()) ? hits-1 : 0), 0); //TNs
								if (hitsInShot[1][1] == false && (regionEnd-regionStart+1 - ((hits > this->pGT->getInternalFrameSize()) ? hits-1 : 0)) > 0) {
									hitsInShot[1][1] = true;
									this->scatterMatrix[sclib::ceShot][c][1][1] += 1;
								}
								if (hitsInScene[1][1] == false && (regionEnd-regionStart+1 - ((hits > this->pGT->getInternalFrameSize()) ? hits-1 : 0)) > 0) {
									hitsInScene[1][1] = true;
									this->scatterMatrix[sclib::ceScene][c][1][1] += 1;
								}
							}
						} //nextBoundaryValid == true || isArtificialBoundary == false

						if (nextBoundaryValid == true || isArtificialBoundary == true) { //continue loop over this scene if the boundary was dubbed invalid just due to artificiality
							b = gtSegEnd;
						} else {
							break;
						}
					} //for b...

				} //scene is to be processed
						
				y = sceneEnd;
			} else { //no more scenes
				break;
			}
		} //for y...
	} //for c...

	MFree_2D(hitsInShot);
	MFree_2D(hitsInScene);

	//now the scatter-matrixes are filled, so calculate the scores
	for (ce = 0; ce < 4; ce++) {
		for (c = 0; c < this->ca; c++) {
			calcDetectionScores(this->scatterMatrix[ce][c], this->recall[ce][c], this->precision[ce][c], this->missRate[ce][c], this->falseAlarmRate[ce][c], this->errorRate[ce][c], this->specificity[ce][c], this->accuracy[ce][c]);
		}
	}

	MFree_0D(pFunc);

  return;
}

//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_Score_ChangeDetection::output(ostream& OutS) {
	unsigned long int c, x, y, ce;
	char buf[sclib::bufferSize];
	SC_Score_ChangeDetection::SC_Boundary *hook;
	double F1;
			
	for (c = 0; c < this->ca; c++) { //loop over all classes
		this->idx2class(c, buf, false); //get the class' long name
		if (c > 0) {
			OutS << "\n";
		}
		OutS << buf << " detection result:\n";
		OutS << "------------------------------------\n";
		OutS << "Scatter-matrix:\n";
		for (ce = 0; ce < 4; ce++) { //loop over all CEs
			if (ce != sclib::ceSegment) { //segment-bases scores are the same as sample-based for chage detection tasks, so don't put them out to avoid redundancy
				ce2text(ce, buf);
				OutS << "\t\tP\tN (per " << buf << ")\n";
				for (y = 0; y < 2; y++) { //table
					OutS << "\t" << ((y == 0) ? "P" : "N") << "\t";
					for (x = 0; x < 2; x++) {
						OutS << this->scatterMatrix[ce][c][y][x] << "\t";
					}
					OutS << "\n";
				}
				OutS << "\n";
			}
		}
		OutS << "                        \tSample-\t\tShot-\t\tScene-based\n";
		OutS << "Recall:                 \t" << setprecision(4) << this->recall[sclib::ceSample][c] << "\t\t" << this->recall[sclib::ceShot][c] << "\t\t" << this->recall[sclib::ceScene][c] << "\n";
		OutS << "Precision:              \t" << setprecision(4) << this->precision[sclib::ceSample][c] << "\t\t" << this->precision[sclib::ceShot][c] << "\t\t" << this->precision[sclib::ceScene][c] << "\n";
		F1 = (this->recall[sclib::ceSample][c]+this->precision[sclib::ceSample][c] > 0.0) ? ((2.0*this->recall[sclib::ceSample][c]*this->precision[sclib::ceSample][c])/(this->recall[sclib::ceSample][c]+this->precision[sclib::ceSample][c])) : 0.0;
		OutS << "F1:                     \t" << setprecision(4) << F1 << "\t\t";
		F1 = (this->recall[sclib::ceShot][c]+this->precision[sclib::ceShot][c] > 0.0) ? ((2.0*this->recall[sclib::ceShot][c]*this->precision[sclib::ceShot][c])/(this->recall[sclib::ceShot][c]+this->precision[sclib::ceShot][c])) : 0.0;
		OutS << F1 << "\t\t";
		F1 = (this->recall[sclib::ceScene][c]+this->precision[sclib::ceScene][c] > 0.0) ? ((2.0*this->recall[sclib::ceScene][c]*this->precision[sclib::ceScene][c])/(this->recall[sclib::ceScene][c]+this->precision[sclib::ceScene][c])) : 0.0;
		OutS << F1 << "\n";
		OutS << "Miss rate (MDR):        \t" << setprecision(4) << this->missRate[sclib::ceSample][c] << "\t\t" << this->missRate[sclib::ceShot][c] << "\t\t" << this->missRate[sclib::ceScene][c] << "\n";
		OutS << "False alarm rate (FAR): \t" << setprecision(4) << this->falseAlarmRate[sclib::ceSample][c] << "\t\t" << this->falseAlarmRate[sclib::ceShot][c] << "\t\t" << this->falseAlarmRate[sclib::ceScene][c] << "\n";
		OutS << "Specificity:            \t" << setprecision(4) << this->specificity[sclib::ceSample][c] << "\t\t" << this->specificity[sclib::ceShot][c] << "\t\t" << this->specificity[sclib::ceScene][c] << "\n";
		OutS << "Omission:               \t" << setprecision(4) << 1.0-this->recall[sclib::ceSample][c] << "\t\t" << 1.0-this->recall[sclib::ceShot][c] << "\t\t" << 1.0-this->recall[sclib::ceScene][c] << "\n";
		OutS << "Commission:             \t" << setprecision(4) << 1.0-this->precision[sclib::ceSample][c] << "\t\t" << 1.0-this->precision[sclib::ceShot][c] << "\t\t" << 1.0-this->precision[sclib::ceScene][c] << "\n";
		OutS << "Error rate:             \t" << setprecision(4) << this->errorRate[sclib::ceSample][c] << "\t\t" << this->errorRate[sclib::ceShot][c] << "\t\t" << this->errorRate[sclib::ceScene][c] << "\n";
		OutS << "Accuracy:               \t" << setprecision(4) << this->accuracy[sclib::ceSample][c] << "\t\t" << this->accuracy[sclib::ceShot][c] << "\t\t" << this->accuracy[sclib::ceScene][c] << "\n";
		OutS << "Fidelity:               \t" << setprecision(4) << 1.0-this->errorRate[sclib::ceSample][c] << "\t\t" << 1.0-this->errorRate[sclib::ceShot][c] << "\t\t" << 1.0-this->errorRate[sclib::ceScene][c] << "\n";
		OutS << "Sensitivity:            \tsee recall\n";
		OutS << "Producer's accuracy:    \tsee omission\n";
		OutS << "User's accuracy:        \tsee commission\n";
		OutS << "True positive rate:     \tsee recall\n";
		OutS << "False negative rate:    \tsee miss rate\n";
		OutS << "False positive rate:    \tsee false alarm rate\n";
		
		//Output raw boundary data according to it's classification to allow better error assessment
		OutS << "\nRaw boundary data output for error assessment (what's reached, what's missed?):\n";
		OutS << "gt:";
		hook = this->gtBoundaryList[c];
		while (hook != NULL) {
			OutS << "\t" << hook->boundary << "\n";
			hook = hook->Next;
		}
		OutS << "TP:";
		hook = this->hypoTPlist[c];
		while (hook != NULL) {
			OutS << "\t" << hook->boundary << "\n";
			hook = hook->Next;
		}
		OutS << "FN:";
		hook = this->hypoFNlist[c];
		while (hook != NULL) {
			OutS << "\t" << hook->boundary << "\n";
			hook = hook->Next;
		}
		OutS << "FP:";
		hook = this->hypoFPlist[c];
		while (hook != NULL) {
			OutS << "\t" << hook->boundary << "\n";
			hook = hook->Next;
		}
	}

  return OutS;
}

//====================================================================================================================
//  Add a boundary to the given list (see above), so that it can be outputted when printing the report
//====================================================================================================================
void SC_Score_ChangeDetection::addBoundary(SC_Score_ChangeDetection::SC_Boundary* &boundaryList, unsigned long int boundary) {
	SC_Score_ChangeDetection::SC_Boundary *hook, *newEntry = new SC_Score_ChangeDetection::SC_Boundary(boundary);
	
	if (boundaryList == NULL) {
		boundaryList = newEntry;
	} else {
		hook = sclib::getLastInList(boundaryList);
		hook->Next = newEntry;
	}

	return;
}
