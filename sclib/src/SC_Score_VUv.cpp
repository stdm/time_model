/**************************************************************************/
/*    Responsibility:																											*/
/*		  - scores for voiced/unvoiced speech decision                      */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.03.2006																								*/
/**************************************************************************/

#include <stdio.h>
#include <iomanip>
#include "SC_Aux.h"
#include "SC_Score_VUv.h"
#include "SC_GroundTruth_TIMIT.h"
#include "SC_MatrixFunctions.h"

#define SCLIB_MODE_V2UV_ERROR 1
#define SCLIB_MODE_UV2V_ERROR 2
#define SCLIB_MODE_OCCURENCE  3
#define SCLIB_MODE_DISCRETE_OCCURENCE 4

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Score_VUv::SC_Score_VUv(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT) : SC_Score(pTweak, pGT) {
	unsigned int x;
	
	//for the standard detection measures for v and uv and the classification measures between this two classes
	this->ca = 2; //voiced & unvoiced speech
	for (x = 0; x < 4; x++) {
		this->scatterMatrix[x] = NULL;
		this->recall[x] = NULL;
		this->precision[x] = NULL;
		this->missRate[x] = NULL;
		this->falseAlarmRate[x] = NULL;
		this->errorRate[x] = NULL;
		this->specificity[x] = NULL;
		this->accuracy[x] = NULL;
		this->averageOmission[x] = 0.0;
		this->averageCommission[x] = 0.0;
		this->missclassificationRate[x] = 0.0;
		this->kappaStatistic[x] = 0.0;
	}
	sprintf(this->purpose, "%s\0", "Voiced/Unvoiced Speech Classification Assessment");

	//for phone-specific details
  this->phoneListSize = 59;
	MArray_2D(this->phoneResultList, this->phoneListSize, 5, long int, "SC_Score_VUv: phoneResultList"); 
	this->phoneResultList[0][0] = sclib::phone_b;
	this->phoneResultList[1][0] = sclib::phone_d;    
	this->phoneResultList[2][0] = sclib::phone_g;       
	this->phoneResultList[3][0] = sclib::phone_p;
	this->phoneResultList[4][0] = sclib::phone_t ;      
	this->phoneResultList[5][0] = sclib::phone_k;     
 	this->phoneResultList[6][0] = sclib::phone_dx;      
	this->phoneResultList[7][0] = sclib::phone_q;    
	this->phoneResultList[8][0] = sclib::phone_bcl;    
	this->phoneResultList[9][0] = sclib::phone_dcl;     
	this->phoneResultList[10][0] = sclib::phone_gcl;    
 	this->phoneResultList[11][0] = sclib::phone_pcl;    
 	this->phoneResultList[12][0] = sclib::phone_tck;     
 	this->phoneResultList[13][0] = sclib::phone_kcl;     
 	this->phoneResultList[14][0] = sclib::phone_tcl;     
 	this->phoneResultList[15][0] = sclib::phone_jh;    
 	this->phoneResultList[16][0] = sclib::phone_ch;      
 	this->phoneResultList[17][0] = sclib::phone_s;     
 	this->phoneResultList[18][0] = sclib::phone_sh;      
 	this->phoneResultList[19][0] = sclib::phone_z;     
 	this->phoneResultList[20][0] = sclib::phone_zh;      
 	this->phoneResultList[21][0] = sclib::phone_f;     
 	this->phoneResultList[22][0] = sclib::phone_th;      
 	this->phoneResultList[23][0] = sclib::phone_v;     
 	this->phoneResultList[24][0] = sclib::phone_dh;      
 	this->phoneResultList[25][0] = sclib::phone_m;      
 	this->phoneResultList[26][0] = sclib::phone_n;     
 	this->phoneResultList[27][0] = sclib::phone_ng;    
 	this->phoneResultList[28][0] = sclib::phone_em;      
 	this->phoneResultList[29][0] = sclib::phone_en;    
 	this->phoneResultList[30][0] = sclib::phone_eng;     
 	this->phoneResultList[31][0] = sclib::phone_nx;      
 	this->phoneResultList[32][0] = sclib::phone_l;      
 	this->phoneResultList[33][0] = sclib::phone_r;     
 	this->phoneResultList[34][0] = sclib::phone_w;       
 	this->phoneResultList[35][0] = sclib::phone_y;       
 	this->phoneResultList[36][0] = sclib::phone_hh;     
 	this->phoneResultList[37][0] = sclib::phone_hv;      
 	this->phoneResultList[38][0] = sclib::phone_el;     
 	this->phoneResultList[39][0] = sclib::phone_iy;     
 	this->phoneResultList[40][0] = sclib::phone_ih;    
 	this->phoneResultList[41][0] = sclib::phone_eh;     
 	this->phoneResultList[42][0] = sclib::phone_ey;     
 	this->phoneResultList[43][0] = sclib::phone_ae;     
 	this->phoneResultList[44][0] = sclib::phone_aa;     
 	this->phoneResultList[45][0] = sclib::phone_aw;      
 	this->phoneResultList[46][0] = sclib::phone_ay;     
 	this->phoneResultList[47][0] = sclib::phone_ah;     
 	this->phoneResultList[48][0] = sclib::phone_ao;      
 	this->phoneResultList[49][0] = sclib::phone_oy;     
 	this->phoneResultList[50][0] = sclib::phone_ow;      
 	this->phoneResultList[51][0] = sclib::phone_uh;  
 	this->phoneResultList[52][0] = sclib::phone_uw;      
 	this->phoneResultList[53][0] = sclib::phone_ux;      
 	this->phoneResultList[54][0] = sclib::phone_er;   
 	this->phoneResultList[55][0] = sclib::phone_ax;      
 	this->phoneResultList[56][0] = sclib::phone_ix;   
 	this->phoneResultList[57][0] = sclib::phone_axr;     
	this->phoneResultList[58][0] = sclib::phone_ax_h;
	//exlcude ohter, non-phonemic symbols
	for(x = 0;  x < (unsigned int)(this->phoneListSize); x++){
		this->phoneResultList[x][1] = 0;
		this->phoneResultList[x][2] = 0;
		this->phoneResultList[x][3] = 0;
		this->phoneResultList[x][4] = 0;
	}
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_Score_VUv::~SC_Score_VUv() {
	for (unsigned int x = 0; x < 4; x++) {
		MFree_2D(this->scatterMatrix[x]);
		MFree_1D(this->scatterMatrix[x]);
		MFree_1D(this->recall[x]);
		MFree_1D(this->precision[x]);
		MFree_1D(this->missRate[x]);
		MFree_1D(this->falseAlarmRate[x]);
		MFree_1D(this->errorRate[x]);
		MFree_1D(this->specificity[x]);
		MFree_1D(this->accuracy[x]);
	}
	MFree_2D(this->phoneResultList);
}

//====================================================================================================================
//	Convert class-tags as used in the ground-truth (sclib::atVoiced etc.) to indices into the scatter-matrix/result-
//  vectors; SVLIB_Fail is returned if the mapping can't be established
//====================================================================================================================
long int SC_Score_VUv::class2idx(unsigned long int classTag) {
	long int res = SVLIB_Fail;

	switch (classTag) {
		case sclib::atVoiced:
			res = 0;
			break;
		case sclib::atUnvoiced:
			res = 1;
			break;
	}

	return res;
}

//====================================================================================================================
//	Convert class indices into the scatter-matrix/result-vectors to class-tags (sclib::atVoiced etc.) as used in the 
//  ground-truth; SVLIB_Fail is returned if the mapping can't be established; if className is != NULL, a string 
//  giving the name of the class (verbose or short according to shortName) is also returned in that variable
//====================================================================================================================
long int SC_Score_VUv::idx2class(unsigned long int classIdx, char *className, bool shortName) {
	long int res = SVLIB_Fail;

	switch (classIdx) {
		case 0:
			res = sclib::atVoiced;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "V");
				} else {
					sprintf(className, "%s\0", "Voiced speech");
				}
			}
			break;
		case 1:
			res = sclib::atUnvoiced;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "Uv");
				} else {
					sprintf(className, "%s\0", "Unvoiced speech");
				}
			}
			break;
		case 2: //this is not really another class, but just helps to give the Pop-Row of the scatter-matrix a label during outputting...
			res = sclib::noType;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%s\0", "#");
				} else {
					sprintf(className, "%s\0", "Overall available");
				}
			}
			break;
	}

	return res;
}

//====================================================================================================================
// Insert values according to mode into phoneResultList according to the phoneCode
// mode may be:
//  - SCLIB_MODE_V2UV_ERROR
//  - SCLIB_MODE_UV2V_ERROR
//  - SCLIB_MODE_OCCURENCE
//  - SCLIB_MODE_DISCRETE_OCCURENCE
//====================================================================================================================
void SC_Score_VUv::insertIntoPhoneResultList(unsigned long int phoneCode, int mode, long int count) {
	int idx = 3;
	
	switch (mode) {
		case SCLIB_MODE_OCCURENCE: 
			idx = 3;
			break;
		case SCLIB_MODE_V2UV_ERROR:
			idx = 1;
			break;
		case SCLIB_MODE_UV2V_ERROR:
			idx = 2;
			break;
		case SCLIB_MODE_DISCRETE_OCCURENCE:
			idx = 4;
			break;
	}

	for (int i = 0; i < this->phoneListSize; i++) {
    if (this->phoneResultList[i][0] == phoneCode) {
      this->phoneResultList[i][idx] += count;
      break;
    }
	}

	return;
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
void SC_Score_VUv::calcScores(unsigned long int start, unsigned long int end, unsigned long int algorithmicUncertaintyDiameter){
	long int gtSegStart, gtSegEnd, hypoSegStart, hypoSegEnd, oldGtSegStart = sclib::noSegment;
	unsigned long int gtType, hypoType, softBoundaryDiameter = this->pGT->getUncertaintyRegionWidth(true) + algorithmicUncertaintyDiameter;
	long int idxGt, idxHypo, idxV = class2idx(sclib::atVoiced), idxUv = class2idx(sclib::atUnvoiced);
	unsigned long int x, y, z, internalFrameSize = this->pGT->getInternalFrameSize(), scene = 0, lastScene = 0, shot = 0, lastShot = 0;
	int *phoneList = NULL, listLength, ce;
	bool **hitsInScene, **hitsInShot, gtEndIsShotBoundary, firstHypoSegment;
	SC_MatrixFunctions *pFunc = NULL;
	long int lastPhone = sclib::noPhone;

	//initializations
	this->start = start;
	this->end = (end > 0 && end < this->pGT->getAudioSampleCount()) ? end : this->pGT->getAudioSampleCount()-1;
	pFunc = new SC_MatrixFunctions();
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

	//here, only segments are considered which are labeled as SPEECH in both groundtruth and algorithmic results (hypo)
	//this way, the errors already counted by SC_Score_AudioTypeClassification are not considered once more
	//for each gt-v/uv-segment and all the hypo-v/uv-segments therin, the correspondences of the givel labels are checked
	for (y = this->start; y < this->end; y++) {
		gtType = this->pGT->getClosestSegment(y, gtSegStart, gtSegEnd, sclib::atVoiced|sclib::atUnvoiced, false, sclib::searchForward, sclib::modeGroundtruth, false, true, sclib::atShotBoundary);
		if (gtType != sclib::noType && gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment && gtSegStart <= (long)(this->end)){
			idxGt = class2idx(gtType);
			gtEndIsShotBoundary = (this->pGT->testSegment(gtSegEnd+1, gtSegEnd+1, true, sclib::atShotBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) ? true : false;
			firstHypoSegment = true;

			//check for scene- and shot-changes
			scene = this->pGT->sample2scene(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastScene);
			if (scene > this->pTweak->general.lastScene) {
				break;
			}

			//has this scene been processed by any algorithm?
			if (scene >= this->pTweak->general.firstScene && (this->pTweak->general.sceneSelection == 0 || sclib::bitTest(sclib::bit(scene), pTweak->general.sceneSelection) == true) && scene <= this->pTweak->general.lastScene) {

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
					hypoType = this->pGT->getClosestSegment(x, hypoSegStart, hypoSegEnd, sclib::atVoiced|sclib::atUnvoiced, false, sclib::searchWithin, sclib::modeHypothesized, false, true, sclib::atShotBoundary);
					if (hypoType != sclib::noType && hypoSegStart != sclib::noSegment && hypoSegEnd != sclib::noSegment && hypoSegStart <= (long)(gtSegEnd)) {
						idxHypo = class2idx(hypoType);

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

						//fill the phone-based error statistics (also only if there is SPEECH detected both in hypo & gt)
						//but we can only calculate phone-based scores with phone-groundtruth available, and this is only true (at the moment) for the TIMIT data...
						if (this->pGT->getGTtype() == sclib::gtTIMIT) {
							phoneList = ((SC_GroundTruth_TIMIT*)this->pGT)->getPhones(hypoSegStart, hypoSegEnd, listLength, sclib::modeGroundtruth);
							for (z = 0; z < (unsigned int)(listLength); z++) {
								insertIntoPhoneResultList(phoneList[z], SCLIB_MODE_OCCURENCE, internalFrameSize); //in [samples]
								if (gtType != hypoType) {
									insertIntoPhoneResultList(phoneList[z], (gtType == sclib::atVoiced) ? SCLIB_MODE_V2UV_ERROR : SCLIB_MODE_UV2V_ERROR, internalFrameSize); //in [samples]
								}
								if (phoneList[z] != lastPhone) { //if the phone in the current FLI differs from the previous one (also from previous segment via lastPhone), a new discrete phone occurence for the current phone is found
									insertIntoPhoneResultList(phoneList[z], SCLIB_MODE_DISCRETE_OCCURENCE, 1);
								}
								lastPhone = phoneList[z]; //for next loop as well as next segment
							}
							MFree_1D(phoneList);
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
			y = gtSegEnd;
		} else {
			break;
		}
	}

	//calculate the scores from the previously fetched raw data
	for (ce = 0; ce < 4; ce++) {
		calcSupervisedClassificationScores(this->scatterMatrix[ce], this->ca, this->averageOmission[ce], this->averageCommission[ce], this->missclassificationRate[ce], this->kappaStatistic[ce], this->recall[ce], this->precision[ce], this->missRate[ce], this->falseAlarmRate[ce], this->errorRate[ce], this->specificity[ce], this->accuracy[ce]);
	}
	
	//sort the phoneResultList
	sclib::quickSort(this->phoneResultList, 0, this->phoneListSize-1, 5, 1);

	MFree_2D(hitsInScene);
	MFree_2D(hitsInShot);
	MFree_0D(pFunc);

	return;
}

//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_Score_VUv::output(ostream& OutS) {
	unsigned long int x, y, ce;
	char buf[sclib::bufferSize];
	double occurrencesErrrorRate, errorsErrorRate, errSum = 0.0, occurenceSum = 0.0, avgLength;
	
	//voiced/unvoiced classification scores
	OutS << "\n";
	OutS << "Voiced vs. unvoiced speech classification result:\n";
	OutS << "--------------------------------------------------\n";
	
	if (this->scatterMatrix[sclib::ceSample] != NULL) {	
		OutS << "Scatter-matrix:\n";
		for (ce = 0; ce < 4; ce++) { //loop over all CEs
			OutS << "\t";
			for (x = 0; x < this->ca; x++) { //headline
				idx2class(x, buf, true); //get the class' short name
				OutS << "\t" << buf;
			}
			ce2text(ce, buf);
			OutS << " (per " << buf << ")\n";
			for (y = 0; y < this->ca+1; y++) { //table+pop+pos
				idx2class(y, buf, true); //get the class' short name
				OutS << "\t" << buf << "\t";
				for (x = 0; x < 2; x++) {
					OutS << this->scatterMatrix[ce][y][x] << "\t";
				}
				OutS << "\n";
			}
			OutS << "\n";
		}
		OutS << "                         Sample- Seg.-   Shot-   Scene-based\n";
		OutS << "Average omission:        " << left << setprecision(4) << setw(8) << this->averageOmission[sclib::ceSample] << left << setprecision(4) << setw(8) << this->averageOmission[sclib::ceSegment] << left << setprecision(4) << setw(8) << this->averageOmission[sclib::ceShot] << left << setprecision(4) << setw(8) << this->averageOmission[sclib::ceScene] << "\n";
		OutS << "Average commission:      " << left << setprecision(4) << setw(8) << this->averageCommission[sclib::ceSample]<< left  << setprecision(4) << setw(8) << this->averageCommission[sclib::ceSegment] << left << setprecision(4) << setw(8) << this->averageCommission[sclib::ceShot] << left << setprecision(4) << setw(8) << this->averageCommission[sclib::ceScene] << "\n";
		OutS << "Missclassification rate: " << left << setprecision(4) << setw(8) << this->missclassificationRate[sclib::ceSample] << left << setprecision(4) << setw(8) << this->missclassificationRate[sclib::ceSegment] << left << setprecision(4) << setw(8) << this->missclassificationRate[sclib::ceShot] << left << setprecision(4) << setw(8) << this->missclassificationRate[sclib::ceScene] << "\n";
		OutS << "Kappa statistic:         " << left << setprecision(4) << setw(8) << this->kappaStatistic[sclib::ceSample] << left << setprecision(4) << setw(8) << this->kappaStatistic[sclib::ceSegment] << left << setprecision(4) << setw(8) << this->kappaStatistic[sclib::ceShot] << left << setprecision(4) << setw(8) << this->kappaStatistic[sclib::ceScene] << "\n";
		OutS << "K-hat index:             see kappa statistic\n";
		
		//detection scores for each class separately
		OutS << "\n";
		OutS << "Detection results for each class vs. rest:\n";
		OutS << "-------------------------------------------\n";
		OutS << "             \t\t\tSample-\tSeg.-\tShot-\tScene-based\n";
		for (x = 0; x < this->ca; x++) { //loop over all classes
			this->idx2class(x, buf, false); //get the class' long name
			OutS << "-> " << buf << " detection:\n";
			OutS << "    Recall: \t\t\t" << setprecision(4) << setw(6) << this->recall[sclib::ceSample][x] << "\t" << this->recall[sclib::ceSegment][x] << "\t" << this->recall[sclib::ceShot][x] << "\t" << this->recall[sclib::ceScene][x] << "\n";
			OutS << "    Precision: \t\t\t" << setprecision(4) << setw(6) << this->precision[sclib::ceSample][x] << "\t" << this->precision[sclib::ceSegment][x] << "\t" << this->precision[sclib::ceShot][x] << "\t" << this->precision[sclib::ceScene][x] << "\n";
			OutS << "    Specificity: \t\t" << setprecision(4) << setw(6) << this->specificity[sclib::ceSample][x] << "\t" << this->specificity[sclib::ceSegment][x] << "\t" << this->specificity[sclib::ceShot][x] << "\t" << this->specificity[sclib::ceScene][x] << "\n";
			OutS << "    Omission: \t\t\t" << setprecision(4) << setw(6) << 1.0-this->recall[sclib::ceSample][x] << "\t" << 1.0-this->recall[sclib::ceSegment][x] << "\t" << 1.0-this->recall[sclib::ceShot][x] << "\t" << 1.0-this->recall[sclib::ceScene][x] << "\n";
			OutS << "    Commission: \t\t" << setprecision(4) << setw(6) << 1.0-this->precision[sclib::ceSample][x] << "\t" << 1.0-this->precision[sclib::ceSegment][x] << "\t" << 1.0-this->precision[sclib::ceShot][x] << "\t" << 1.0-this->precision[sclib::ceScene][x] << "\n";
			OutS << "    Miss rate: \t\t\t" << setprecision(4) << setw(6) << this->missRate[sclib::ceSample][x] << "\t" << this->missRate[sclib::ceSegment][x] << "\t" << this->missRate[sclib::ceShot][x] << "\t" << this->missRate[sclib::ceScene][x] << "\n";
			OutS << "    False alarm rate: \t\t" << setprecision(4) << setw(6) << this->falseAlarmRate[sclib::ceSample][x] << "\t" << this->falseAlarmRate[sclib::ceSegment][x] << "\t" << this->falseAlarmRate[sclib::ceShot][x] << "\t" << this->falseAlarmRate[sclib::ceScene][x] << "\n";
			OutS << "    Error rate: \t\t" << setprecision(4) << setw(6) << this->errorRate[sclib::ceSample][x] << "\t" << this->errorRate[sclib::ceSegment][x] << "\t" << this->errorRate[sclib::ceShot][x] << "\t" << this->errorRate[sclib::ceScene][x] << "\n";
			OutS << "    Accuracy: \t\t\t" << setprecision(4) << setw(6) << this->accuracy[sclib::ceSample][x] << "\t" << this->accuracy[sclib::ceSegment][x] << "\t" << this->accuracy[sclib::ceShot][x] << "\t" << this->accuracy[sclib::ceScene][x] << "\n";
			OutS << "    Fidelity: \t\t\t" << setprecision(4) << setw(6) << 1.0-this->errorRate[sclib::ceSample][x] << "\t" << 1.0-this->errorRate[sclib::ceSegment][x] << "\t" << 1.0-this->errorRate[sclib::ceShot][x] << "\t" << 1.0-this->errorRate[sclib::ceScene][x] << "\n";
			OutS << "    Sensitivity: \t\tsee recall\n";
			OutS << "    Producer's accuracy: \tsee omission\n";
			OutS << "    User's accuracy: \t\tsee commission\n";
			OutS << "    True positive rate: \tsee recall\n";
			OutS << "    False negative rate: \tsee miss rate\n";
			OutS << "    False positive rate: \tsee false alarm rate\n";
		}

		//errors per phone, if we have phone-groundtruth
		if (pGT->getGTtype() == sclib::gtTIMIT) {
			OutS << "\n";
			OutS << "Errors per phone:\n";
			OutS << "------------------\n";
			OutS << "P-Error: Percentage of Errors [samples] as measured by occurences [samples] of this Phone\n";
			OutS << "O-Error: Precentage of Errors [samples] of this phone as measured by Overall number of errors [samples]\n\n";
			OutS << "  Phone    V->UV err.  UV->V err.  Occurences [samples]  Occurences [discrete]  avg. length [ms]  P-Error [%]  O-Error [%]\n";
			
			for (y = 0; y < (unsigned int)(this->phoneListSize); y++) { //calculate overall error based on the phoneResultList
				errSum += this->phoneResultList[y][1] + this->phoneResultList[y][2];
				occurenceSum += this->phoneResultList[y][3];
			}

			for (y = 0; y < (unsigned int)(this->phoneListSize); y++) {
				occurrencesErrrorRate = (this->phoneResultList[y][3] > 0.0) ? ((this->phoneResultList[y][1] + this->phoneResultList[y][2]) / (double)(this->phoneResultList[y][3])) * 100.0 : 0.0;
				errorsErrorRate = (errSum > 0.0) ? ((this->phoneResultList[y][1] + this->phoneResultList[y][2]) / errSum) * 100 : 0.0;
				avgLength = (this->phoneResultList[y][4] > 0) ? this->pGT->getConverter()->sample2ms(sclib::round((double)(this->phoneResultList[y][3]) / (double)(this->phoneResultList[y][4]))) : 0.0;

				OutS << "  "	<< left											<< setw(9)	<< ((SC_GroundTruth_TIMIT*)pGT)->phoneType2string(this->phoneResultList[y][0])
					 						<< left << setprecision(4)	<< setw(12) << this->phoneResultList[y][1]
 											<< left << setprecision(4)	<< setw(12) << this->phoneResultList[y][2]
											<< left											<< setw(22) << this->phoneResultList[y][3]
											<< left                     << setw(23) << this->phoneResultList[y][4]
											<< left << setprecision(4)  << setw(18) << avgLength
											<< left << setprecision(4)	<< setw(13) << occurrencesErrrorRate
											<< left << setprecision(4)							<< errorsErrorRate << "\n";
			}

			OutS << "\n  Overall erros: " <<  errSum << " / " << occurenceSum << " = " << errSum/occurenceSum * 100.0 << "%\n";
		}
	} else { //no scores
		OutS << "Not available.\n";
	}

  return OutS;
}

//====================================================================================================================
//	returns the overall error rate (the one specified by errorType: v2uv or vice versa) for the given phone; -1 is 
//  returned if the given phoneType is unknown
//====================================================================================================================
double SC_Score_VUv::getPhoneError(unsigned long int phoneCode, int errorType) {
	int errIdx = (errorType == sclib::etV2Uv) ? 1 : 2;
	double res = -1.0;
	
	for (int i = 0; i < this->phoneListSize; i++) {
    if (this->phoneResultList[i][0] == phoneCode) {
      res = (double)(this->phoneResultList[i][errIdx]) / (double)(this->phoneResultList[i][3]);
      break;
    }
	}

	return res;
}

//====================================================================================================================
//	returns the number of occurences for the given phones (0 for unknown phoness)
//====================================================================================================================
unsigned long int SC_Score_VUv::getPhoneOccurences(unsigned long int phoneCode) {
	unsigned long int res = 0;
	
	for (int i = 0; i < this->phoneListSize; i++) {
    if (this->phoneResultList[i][0] == phoneCode) {
      res = this->phoneResultList[i][3];
      break;
    }
	}

	return res;
}
