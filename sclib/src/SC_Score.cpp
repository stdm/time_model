/**************************************************************************/
/*    Responsibility:																											*/
/*		  - base class for computing scores (missmatch between groundtruth  */
/*        and algorithmic results) according to different measures and    */
/*        different tasks                                                 */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.03.2006																								*/
/**************************************************************************/

#include <stdio.h>
#include <iomanip>
#include <time.h>
#include "SC_Aux.h"
#include "SC_Score.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Score::SC_Score(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT) {
  this->pTweak = pTweak;
  if (this->pTweak == NULL) {
    REPORT_ERROR(SVLIB_BadArg, "Tweakable parameters are needed by SC_Score*!");
  }
	this->pGT = pGT;
  if (this->pGT == NULL) {
    REPORT_ERROR(SVLIB_BadArg, "Ground-truth is needed by SC_Score*!");
  }
	this->Next = NULL;
	this->start = 0;
	this->end = this->pGT->getAudioSampleCount() - 1;
	sprintf(this->purpose, "\0");
}

//====================================================================================================================
//	The destructor; doesn't destruct the linked list in which this object might be a member
//====================================================================================================================
SC_Score::~SC_Score() {

}

//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_Score::output(ostream& OutS) {
	//this method (instead of putting it into the operator<<()-friend) makes this operator kind of virtual
	return OutS; //no specialized output in the base class because here nothing is measured
}
ostream& operator<< (ostream& OutS, SC_Score& pScore) {
  return pScore.output(OutS);
}

//====================================================================================================================
//	Presents the results, also care for the output of the complete framelist when wished in the tweakable parameters
//  this method works bestthe following way: a linked list of score-objects in the order of processing should be 
//  created (e.g. ATC->V/Uv->CD->Clustering->ID for speaker diarization) and then the printReport() metod of the first
//  list element should be invoked with parameters such sthat the introduxctory header, the complete framelist and the 
//  other scores are printed. This will pwoduce a complete estimony of the overall analyzation run and all scoreable
//  algorithms in a suitable form. This procedure can on wish be repeated for different parts of the analyzed data: 
//  For the whole thing at once (start=0, end=pGT->getAudioSampleCount()-1) and for each scene separately...
//====================================================================================================================
void SC_Score::printReport(const char* fileName, bool introductoryHeader, unsigned long int elapsedSeconds, bool allowFrameListOut, bool ignoreNextScore) {
  unsigned long sceneNr = 1;
  bool printedSth = false;
	char* fName = NULL;
	fstream fileOut;
	time_t epoches;
	tm *now;

	fName = new char[strlen(this->pTweak->debug.debugDir) + strlen(fileName) + 1];
	sprintf(fName, "%s%s\0", this->pTweak->debug.debugDir, fileName);
	fileOut.open(fName, ios_base::out|ios_base::app);
  MFree_1D(fName);

	//print an introductory header including the filename of the analyzed data, from where to where it was measured, how long that took, when the run did take place, and the used parameter-settings
	if (introductoryHeader == true) {
		time(&epoches);
		now = localtime(&epoches);
		
		//print out header-information
		fileOut << "\n================================================================================\n";
		fileOut << "File: " << this->pGT->getAudioFileName() << "\n";
		fileOut << "Run: " << now->tm_year + 1900 << "-" << setw(2) << now->tm_mon + 1 << "-" << setw(2) << now->tm_mday << " " << setw(2) << now->tm_hour << ":" << setw(2) << now->tm_min << ":" << setw(2) << now->tm_sec << " UTC\n";
		fileOut << "Elapsed time: " << setw(5) << (elapsedSeconds / 60) << ":" << setw(2) << (elapsedSeconds % 60) << "\n";

		//print out the current parameter-set in verbose form
		fileOut << "Parameter-settings: " << "\n\n";
		if (this->pTweak != NULL) {
			this->pTweak->toggleVerboseMode(true);
			fileOut << *(this->pTweak);
		}

		//print out explanations of different measures:
		fileOut << "\n--------------------------------------------------------------------------------\n";
		fileOut << "Meaning of the implemented measures: \n";
		fileOut << "  Detection scores: \n";
		fileOut << "    recall                     := likelihood of a positive item to be detected (relevant item to be retrieved)\n";
		fileOut << "    precision                  := likelihood of a detected item to be positive\n";
		fileOut << "    miss rate                  := likelihood of a positive item to be not detected\n";
		fileOut << "    FA rate                    := likelihood of a detected item to be positive\n";
		fileOut << "    specificity                := likelihood of a not-detected item to be negative\n";
		fileOut << "    error rate                 := likelihood of an incorrect result based on all processed items\n";
		fileOut << "    accuracy                   := likelihood of a correct result result based to all items in the analyzed part\n";
		fileOut << "    fidelity                   := reciprocal of error rate; like accuracy, but maybe based on fewer items\n";
		fileOut << "    other measures are common synomyms to the abovementioned scores\n";
		fileOut << "  Supervised classification scores: \n";
		fileOut << "    average omission           := likelihood of an item of class X to be recognized as belonging to any other class (exclusion error)\n";
		fileOut << "    average commission         := likelihood of an item not belonging to the class it was classified as (erroneous inclusion error) \n";
		fileOut << "    missclassification rate    := likelihood of an item to be missclassified based on all possible items\n";
		fileOut << "    kappa statistic            := likelihood of the agreement between the results and the groundtruth to be only due to chance\n";
		fileOut << "  Unsupervised classification scores: \n";
		fileOut << "    overall recall             := likelihood of an item being in the correct cluster (the one with most items from it's class\n";
		fileOut << "    overall precision          := likelihood of an item being in a cluster where items of it's own class are in the majority\n";
		fileOut << "    missclassification rate    := likelihood of an item to be missclassified based on all possible items\n";
		fileOut << "    average cluster purity     := likelihood of the items of one cluster belonging together\n";
		fileOut << "    average class purity       := likelihood of the items classified as belonging to class X really doing so\n";
		fileOut << "    average average precision  := if the agglomerative clustering process is seen as a special case of retrieval, this evaluates the partition/result list in terms of recall, precision & relevance ranking\n";
		fileOut << "    purity                     := average of class- and cluster purity\n";
		fileOut << "    rand index                 := decreases with the number of correctly clustered items (=> unnormalized, but smaller is better)\n";
		fileOut << "    BBN metric                 := increases with the number of big, pure clusters (=> unnormalized, but greater is better)\n";
		fileOut << "    diarization error rate     := ratio of samples assigned to the wrong speaker, including speaker error time, missed speaker time, and false alarm speaker time\n";
		fileOut << "  Meaning of different bases:\n";
		fileOut << "    possible items             := all items that could have been analyzed (can be greater as the next one due to algorithmic reasons)\n";
		fileOut << "    processed items            := all items that really have been analyzed (the smallest set)\n";
		fileOut << "    items in the analyzed part := pure count of items between analyzation-start and -end (the biggest set)\n";
	}

	//print out the specialized measures of the derived classes
	fileOut << "\n--------------------------------------------------------------------------------\n";
	if (this->purpose != NULL) {
		fileOut << this->purpose << "\n";
	} else {
		fileOut << "Score sub-report; no purpose given\n";
	}
	fileOut << "Analyzed section [samples]: " << this->start << "-" << this->end << "\n\n";
	
  //print out scores
  fileOut << *(this);

	fileOut.close();

	//shall the next score out of a linked list be printed, too?
	if (ignoreNextScore == false && this->Next != NULL) {
		this->Next->printReport(fileName, false, 0, false, false);
	}
	
  //print out complete framelist on whish
	if (allowFrameListOut == true) {
		if (this->pTweak->debug.debugMode & sclib::dbCompleteResults) {
			fName = new char[strlen(this->pTweak->debug.debugDir) + strlen("frameList.txt") + 1];
			sprintf(fName, "%s%s\0", this->pTweak->debug.debugDir, "frameList.txt");
  		this->pGT->frameListOut(fName, 0, pGT->getAudioSampleCount() - 1, 0, 0);
			MFree_1D(fName);
		}
	}
  
	return;
}

//====================================================================================================================
//	This is rather simple: Based on the number of overall indetifyable items and the number of correct 
//  identifications, the id-accuracy (or recognition rate) is returned as the percentage of correctItems; the return
//  value may indicvate error by returning SVLIB_Fail, otherwise SVLIB_Ok
//====================================================================================================================
int SC_Score::calcIdentificationScores(long int overallItems, long int correctItems, double &accuracy) {
	int res = SVLIB_Ok;

	accuracy = (overallItems > 0) ? ((double)(correctItems) / (double)(overallItems)) : 1.0;

	return res;
}

//====================================================================================================================
//	Receives a scatter matrix of the following format:
//
//               GroundTruth
//  Hypothesized [ TP | FP ] (number of true positives, false positives of the current CE)
//               [ FN | TN ] (number of false negatives, true negatives of the current CE)
//               [PPop|NPop] (number of positives/negatives according to ground-truth in the population)
//               [PPos|NPos] (number of positives/negatives according to ground-truth in the possible set)
//
//  the "population" is the count of overall available CEs in the area under consideration
//  the "possible set" is the count of CEs actually subject to the detection task (e.g. for change detection, only 
//  segment boundaries are actually searched for changes, where the overall population would also include the parts in 
//  the segments)
//
//  Returns recall (=sensitivity =true positive rate), precision, specificity (=fallout), miss rate (=false negative 
//  rate=MDR), false alarm rate (=FAR=false positive rate), error rate (=1-fidelity) and accuracy as parameters and an 
//  integer indicating an error (if <> SVLIB_OK)
//====================================================================================================================
int SC_Score::calcDetectionScores(long int** scatterMatrix, double& recall, double& precision, double& missRate, double& falseAlarmRate, double& errorRate, double& specificity, double& accuracy) {
  int res = SVLIB_Ok;
  long int TP = scatterMatrix[0][0], FP = scatterMatrix[0][1], FN = scatterMatrix[1][0], TN = scatterMatrix[1][1];
  long int population = scatterMatrix[2][0]+scatterMatrix[2][1], possible = scatterMatrix[3][0]+scatterMatrix[3][1];
  long int processed = TP + FP + FN + TN;

	recall = (TP+FN > 0) ? (double)(TP) / ((double)(TP)+(double)(FN)) : 1.0; //=sensitivity =true positive rate
	precision = (TP+FP > 0) ? (double)(TP) / ((double)(TP)+(double)(FP)) : 1.0;
	missRate = (TP+FN > 0) ? (double)(FN) / ((double)(TP)+(double)(FN)) : 0.0; //=false negative rate =MDR
	//falseAlarmRate = (FP+TN > 0) ?(double)(FP) / (double)(FP+TN) : 0.0; //=false positive rate
	falseAlarmRate = (TP+FN+FP > 0) ?(double)(FP) / ((double)(TP)+(double)(FN)+(double)(FP)) : 0.0; //=FAR as defined by Kotti et al.
	errorRate = (processed > 0) ?((double)(FP)+(double)(FN)) / (double)(processed) : 0.0; //=1-fidelity
	specificity = (TN+FP > 0) ? (double)(TN) / ((double)(TN)+(double)(FP)) : 1.0; //=fallout
	accuracy = (population > 0) ? (((double)(TP)+(double)(TN)) / (double)(population)) : (((double)(FP)+(double)(FN) > 0.0) ? 0.0 : 1.0);
	//F1 = (recall+precision > 0.0) ? ((2.0*recall*precision) / (recall + precision)) : 0.0;

  return res;
}

//====================================================================================================================
//	Receives a scatter matrix ((cu+2) x ca) of the following format:
//
//               GroundTruth (classes: 1..ca)
//  Hypothesized [ x_1_1  | x_1_2 | ... | x_1_ca  ] (the correct classified CEs appear on the main diagonal of the 
//  (classes:    [ x_2_1  | x_2_2 |     | ...     ]  matrix, i.e. x_i_i for i=1..ca; the rest in this upper part shows
//   1..ca)      [ ...    | ...   |     |         ]  the spread of errors)
//               [ x_ca_1 |       |     | x_ca_ca ] 
//               [Pop_1   | ...   |     | Pop_ca  ] (overall number of CEs belonging to this class, maybe not all where clustered)
//               [Pos_1   | ...   |     | Pos_ca  ] (overall number of CEs belonging to this class in the possible set, maybe more than where clustered)
//
//  the "population" is the count of overall available CEs in the area under consideration
//  the "possible" is the count of CEs actually subject to the classification task (e.g. for change detection, only 
//  segment boundaries are actually searched for changes, where the overall population would also include the parts in
//  the segments)
//
//  The correspondence betweeen correct clusters and classes is provided by the arrangement of the rows and cols in
//  the scatter matrix: It is such that cluster/row i belongs to class/col i; 
//
//  Returned are average omission & commission, the missclassification rate and the kappa statistic (=k-hat index) to 
//  characterize the overall result; for ech class separately, recall (=1-omission =1-producer's accuracy) & precision 
//  (=1-commission =1-user's accuracy), miss & FA & error rate, specificity and accuracy (for further details, see 
//  calcDetectionScores() comment) are provided
//
//  If additionalVectorCols>0, the per-class result vectors are made such an amount bigger (may be needed by the 
//  calling function)
//
//  An integer indicating an error (if <> SVLIB_OK) is returned
//====================================================================================================================
int SC_Score::calcSupervisedClassificationScores(long int** scatterMatrix, unsigned long int ca, double& averageOmission, double& averageCommission, double& missclassificationRate, double& kappaStatistic, double*& recall, double*& precision, double*& missRate, double*& falseAlarmRate, double*& errorRate, double*& specificity, double*& accuracy, unsigned int additionalVectorFields) {
  int res = SVLIB_Ok;
	unsigned long int x;
	long int **reducedScatterMatrix = NULL;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

	//calculate values per class
	MFree_1D(recall);
	MArray_1D(recall, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: recall");
	MFree_1D(precision);
	MArray_1D(precision, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: precision");
	MFree_1D(missRate);
	MArray_1D(missRate, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: missRate");
	MFree_1D(falseAlarmRate);
	MArray_1D(falseAlarmRate, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: falseAlarmRate");
	MFree_1D(errorRate);
	MArray_1D(errorRate, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: errorRate");
	MFree_1D(specificity);
	MArray_1D(specificity, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: specificity");
	MFree_1D(accuracy);
	MArray_1D(accuracy, ca+additionalVectorFields, double, "SC_Score.calcSupervisedClassificationScores: accuracy");

	for (x = 0; x < ca; x++) {
		reducedScatterMatrix = pFunc->sumScatterMatrix(scatterMatrix, ca, ca+2, x, 2);
		if (reducedScatterMatrix != NULL) {
			res = calcDetectionScores(reducedScatterMatrix, recall[x], precision[x], missRate[x], falseAlarmRate[x], errorRate[x], specificity[x], accuracy[x]);
			MFree_2D(reducedScatterMatrix);
		}
		if (res != SVLIB_Ok) {
			break;
		}
	}

	//calculate values describing the whole result
	if (res == SVLIB_Ok) {
		res = calcSupervisedClassificationScores(scatterMatrix, ca, averageOmission, averageCommission, missclassificationRate, kappaStatistic);
	}

	MFree_0D(pFunc);

  return res;
}

//====================================================================================================================
//	Calculate just the scores describing the complete supervised classification result (see comment above...)
//====================================================================================================================
int SC_Score::calcSupervisedClassificationScores(long int** scatterMatrix, unsigned long int ca, double& averageOmission, double& averageCommission, double& missclassificationRate, double& kappaStatistic) {
	unsigned long int x, y;
	long int *row, *col;
	double pe, po, possible, correct, overall; //double to overcome biggest number limitiation of ints
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	
	row = pFunc->sum(scatterMatrix, ca, ca, false);
	col = pFunc->sum(scatterMatrix, ca, ca, true);

	possible = 0.0;
	correct = 0.0;
	overall = 0.0;
	averageCommission = 0.0;
	averageOmission = 0.0;
	for (y = 0; y < ca; y++) { //here, x==y (because matrix is square when not considering pop and pos), so the character itself has no semantic meaning regarding rows or cols...
		averageCommission += (row[y] != 0) ? (row[y] - scatterMatrix[y][y]) / (double)(row[y]) : 0.0;
		averageOmission += (col[y] != 0) ? (col[y] - scatterMatrix[y][y]) / (double)(col[y]) : 0.0;
		possible += (double)(scatterMatrix[ca+1][y]); //Pos_* row
		correct += (double)(scatterMatrix[y][y]);
		for ( x = 0; x < ca; x++) {
			overall += (double)(scatterMatrix[y][x]);
		}
	}
	averageCommission /= (ca > 0) ? (double)(ca) : 0.0;
	averageOmission /= (ca > 0) ?(double)(ca) : 0.0; 
	missclassificationRate = (possible > 0) ? (possible - correct) / possible : 0.0;

	po = correct / overall; //observed agreement
	pe = 0.0; //expected agreement
	for (y = 0; y < ca; y++) {
		pe += ((double)(row[y]) / overall) * ((double)(col[y]) / overall);
	}
	kappaStatistic = (po - pe) / (1.0 - pe);
	
	MFree_1D(row);
	MFree_1D(col);
	MFree_0D(pFunc);

  return SVLIB_Ok;
}

//====================================================================================================================
//	Receives a scatter matrix ((cu+2) x ca) of the following format:
//
//               GroundTruth (classes: 1..ca)
//  Hypothesized [ x_1_1  | x_1_2 | ... | x_1_ca  ] (the correct clustered CEs appear on the main diagonal of the 
//  (clusters:   [ x_2_1  | x_2_2 |     | ...     ]  matrix, i.e. x_i_i for i=1..ca; the rest in this upper part shows
//   1..cu)      [ ...    | ...   |     |         ]  the spread of errors)
//               [ x_cu_1 |       |     | x_cu_ca ] (cu >= ca)
//               [Pop_1   | ...   |     | Pop_ca  ] (overall number of CEs belonging to this class, maybe not all where clustered)
//               [Pos_1   | ...   |     | Pos_ca  ] (overall number of CEs belonging to this class in the possible set, maybe more than where clustered)
//
//  the "population" is the count of overall available CEs in the area under consideration
//  the "possible" is the count of CEs actually subject to the classification task (e.g. for change detection, only 
//  segment boundaries are actually searched for changes, where the overall population would also include the parts in
//  the segments)
//
//  The correspondence betweeen correct clusters and classes is provided by the arrangement of the rows and cols in
//  the scatter matrix: It is such that cluster/row i belongs to class/col i; the problem for the just fitting 
//  clusters (redundant clusters, where CEs of one class are in the majority and therefore this class assigns the 
//  cluster's label, but there is a bigger one which is labeled as correct) is solved by ordering them this way:
//  after row ca comes the the biggest cluster fitting (but not being correct) to class 1, followed by the 2nd biggest
//  etc. after all fitting clusters of class 1 comes the biggest cluster fitting to class 2 and so on.
//
//  In addition, some information concerning the order of the classification process, if it was hierarchical, is 
//  expected to compute average precision as a measure of quality of the cluster merging process in an information 
//  retrieval sense. This information is included in the array mergeHistory, where each cell i holds an array j=1..J 
//  of 1/0 content telling if j'th merge was good with respect to cluster i; because each merge-history can have a 
//  different length J, each array is terminated after an entry with the value of "-1". a merge is considered as 
//  "good/1", if the  merged clusters both where fitting to the final id. After the "-1" entry there follows one 
//  last element in the array containing the (gt-based) total number of correct segments for this cluster, to include 
//  recall into the average-precision measure
//  
//                       1.. j ..J
//  1..    cluster_1  -> 1|0|0|1|1|-1|3
//    i    cluster_2  -> 1|1|-1|2
//    ..cu cluster_cu -> 1|0|1|1|1|1|1|0|-1|7
//
//  Returned are overall recall & precision, the missclassification rate (=error rate), average class & cluster 
//  purity, average average precision, purity, the rand indexa and the BBN metric as measures of the overall goodness.
//  To assess goodness in every single class, the following measures are provided per class (sometimes also per 
//  cluster), in the order the classes/clusters appear in the scatter matrix, respectively: average precision, class &
//  cluster purity, all values returned by calcDetectionScores() (as each class on its own can be viewed as a
//  detection problem).
//
//  An integer indicating an error (if <> SVLIB_OK) is returned
//====================================================================================================================
int SC_Score::calcUnsupervisedClassificationScores(long int** scatterMatrix, unsigned long int cu, unsigned long int ca, long int** mergeHistory, long int maxMergeHistoryLength, double& overallRecall, double& overallPrecision, double& missclassificationRate, double& averageClusterPurity, double& averageClassPurity, double& average2precision, double& purity, double& randIndex, double& BBNmetric, double*& averagePrecision, double*& classPurity, double*& clusterPurity, double*& recall, double*& precision, double*& missRate, double*& falseAlarmRate, double*& errorRate, double*& specificity, double*& accuracy) {
  int res = SVLIB_Ok;
	unsigned long int ownerIdx, x, y, i, possible = 0, processed = 0, correct = 0, fitting;
	long int **reducedScatterMatrix = NULL, *row = NULL, *col = NULL;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	double precisionAtC;

	//calculate values per class
	MFree_1D(recall);
	MArray_1D(recall, ca, double, "SC_Score.calcUnsupervisedClassificationScores: recall");
	MFree_1D(precision);
	MArray_1D(precision, ca, double, "SC_Score.calcUnsupervisedClassificationScores: precision");
	MFree_1D(missRate);
	MArray_1D(missRate, ca, double, "SC_Score.calcUnsupervisedClassificationScores: missRate");
	MFree_1D(falseAlarmRate);
	MArray_1D(falseAlarmRate, ca, double, "SC_Score.calcUnsupervisedClassificationScores: falseAlarmRate");
	MFree_1D(errorRate);
	MArray_1D(errorRate, ca, double, "SC_Score.calcUnsupervisedClassificationScores: errorRate");
	MFree_1D(specificity);
	MArray_1D(specificity, ca, double, "SC_Score.calcUnsupervisedClassificationScores: specificity");
	MFree_1D(accuracy);
	MArray_1D(accuracy, ca, double, "SC_Score.calcUnsupervisedClassificationScores: accuracy");
	MFree_1D(clusterPurity);
	MArray_1D(clusterPurity, cu, double, "SC_Score.calcUnsupervisedClassificationScores: clusterPurity");
	MFree_1D(classPurity);
	MArray_1D(classPurity, ca, double, "SC_Score.calcUnsupervisedClassificationScores: classPurity");
	MFree_1D(averagePrecision);
	MArray_1D(averagePrecision, cu, double, "SC_Score.calcUnsupervisedClassificationScores: averagePrecision");

	col = pFunc->sum(scatterMatrix, cu, ca, true);
	for (x = 0; x < ca; x++) {
		reducedScatterMatrix = pFunc->sumScatterMatrix(scatterMatrix, ca, cu+2, x, 2, true);
		if (reducedScatterMatrix != NULL) {
      res = calcDetectionScores(reducedScatterMatrix, recall[x], precision[x], missRate[x], falseAlarmRate[x], errorRate[x], specificity[x], accuracy[x]);
			MFree_2D(reducedScatterMatrix);
		}
		if (res != SVLIB_Ok) {
			break;
		}
		classPurity[x] = 0.0;
		for (y = 0; y < cu; y++) {
			classPurity[x] += (col[x] > 0) ? ((double)(scatterMatrix[y][x]) * (double)(scatterMatrix[y][x])) / ((double)(col[x]) * (double)(col[x])) : 0.0;
		}
		possible += scatterMatrix[cu+1][x]; //Pos_* row
		correct += scatterMatrix[x][x];
	}

	//calculate fitting ces
	fitting = correct; //all correct ces are also fitting
	for (x = ca; x < cu; x++) {
		ownerIdx = pFunc->maxIdx(scatterMatrix[x], ca); //column-index of the owner (the class contributing most ces to it) of this cluster
		fitting += scatterMatrix[x][ownerIdx];
	}

	if (res == SVLIB_Ok) {
		row = pFunc->sum(scatterMatrix, cu, ca, false);
		for (y = 0; y < cu; y++) {
			clusterPurity[y] = 0.0;
			for (x = 0; x < ca; x++) {
				clusterPurity[y] += (row[y] > 0) ? ((double)(scatterMatrix[y][x]) * (double)(scatterMatrix[y][x])) / ((double)(row[y]) * (double)(row[y])) : 0.0;
				processed += scatterMatrix[y][x];
			}
			
			x = 0;
			averagePrecision[y] = 0.0;
			if (mergeHistory != NULL && maxMergeHistoryLength>0) { //could have been switched off
				while (mergeHistory[y][x]!=-1 && x<(unsigned long int)(maxMergeHistoryLength)) {
					i = 0;
					precisionAtC = 0.0;
					while (i <= x) {
						precisionAtC += mergeHistory[y][i];
						i++;
					}
					precisionAtC /= (double)(x+1);
					averagePrecision[y] += precisionAtC * ((mergeHistory[y][x] >= 0) ? mergeHistory[y][x] : 0);
					x++;
				}
				if (mergeHistory[y][0] == -1) { //for a cluster with no mergences, the average-precision is just determined by the recall
					averagePrecision[y] = (x+1 < (unsigned long int)(maxMergeHistoryLength)) ?  1.0/mergeHistory[y][x+1] : 0.0;
				} else {
					if (x+1 < (unsigned long int)(maxMergeHistoryLength)) {
						averagePrecision[y] /= (double)(sclib::max(1.0, mergeHistory[y][x+1]-1)); //the (x+1)'st entry of the merge-history-rows contains the total number of relevant/correct segments that could have been merged if one big pure cluster would have been created 
					}
				}
			}
		}
	}

	//calculate values describing the whole result
	if (res == SVLIB_Ok) {
		overallRecall = (possible > 0) ? (double)(correct) / (double)(possible) : 0.0;
		overallPrecision = (processed > 0) ? (double)(fitting) / (double)(processed) : 0.0;
		missclassificationRate = (possible > 0) ? (double)(possible - correct) / (double)(possible) : 0.0;

		BBNmetric = -1.0 * this->pTweak->score.BBNmetricLambda * (double)(cu); //penalty for more unpure but bigger clusters, typically = 0.5
		randIndex = 0.0;
		averageClassPurity = 0.0;
		averageClusterPurity = 0.0;
		average2precision = 0.0;
		for (y = 0; y < cu; y++) {
			BBNmetric += row[y] * clusterPurity[y];
			randIndex += (0.5 - clusterPurity[y]) * row[y] * row[y];
			averageClusterPurity += clusterPurity[y] * row[y];
			average2precision += averagePrecision[y];
			if (y < ca) {
				randIndex += 0.5 * col[y] *col[y];
				averageClassPurity += classPurity[y] * col[y];
			}
		}
		averageClassPurity = (possible > 0) ? averageClassPurity / (double)(possible) : 0.0;
		averageClusterPurity = (possible > 0) ? averageClusterPurity / (double)(possible) : 0.0;
		average2precision /= (double)(cu);
		purity = sqrt(averageClassPurity * averageClusterPurity);
	}

	MFree_1D(row);
	MFree_1D(col);
	MFree_0D(pFunc);

  return res;
}

//====================================================================================================================
//	Copies a string containing a description of the given CE or "" if the CE is invalid into the elsewhere alocated 
//  buffer
//====================================================================================================================
void SC_Score::ce2text(unsigned int countedEntity, char *buffer) {
	if (buffer != NULL) {
		switch (countedEntity) {
			case sclib::ceSample:
				sprintf(buffer, "%s\0", "Sample");
				break;
			case sclib::ceSegment:
				sprintf(buffer, "%s\0", "Segment");
				break;
			case sclib::ceShot:
				sprintf(buffer, "%s\0", "Shot");
				break;
			case sclib::ceScene:
				sprintf(buffer, "%s\0", "Scene");
				break;
			default:
				sprintf(buffer, "%s\0", "");
				break;
		}
	}

	return;
}

//====================================================================================================================
//	Flips the elements of a 4x2 detection scatter matrix such that the other class' values are now the positives
//====================================================================================================================
void SC_Score::flipDetectionScatterMatrix(long int** scatterMatrix) {
	long int tmp;

	//               GroundTruth			
	//  Hypothesized [ TP | FP ]			[ TN | FN ]
	//               [ FN | TN ]  =>	[ FP | TP ]
	//               [PPop|NPop]			[NPop|PPop]
	//               [PPos|NPos]			[NPos|PPos]

	tmp = scatterMatrix[0][0];
	scatterMatrix[0][0] = scatterMatrix[1][1];
	scatterMatrix[1][1] = tmp;

	tmp = scatterMatrix[0][1];
	scatterMatrix[0][1] = scatterMatrix[1][0];
	scatterMatrix[1][0] = scatterMatrix[0][1];

	tmp = scatterMatrix[2][0];
	scatterMatrix[2][0] = scatterMatrix[2][1];
	scatterMatrix[2][1] = tmp;

	tmp = scatterMatrix[3][0];
	scatterMatrix[3][0] = scatterMatrix[3][1];
	scatterMatrix[3][1] = tmp;

	return;
}

//====================================================================================================================
//	ranking1 and ranking2 represent one classification experiment per line, with the score for different clases in 
//  the columns; both experiments are meant to are identical, except that the models/classifiers used to create the
//  scores shall differ. The concordance score gives now a value in the interval [0..1] of how much both 
//  classifiers/models agree in their ranking (classification) on all possible pairs of classes (columns) for each
//  row. 1 means perfect concordance (match), 0 means total discordance
//====================================================================================================================
double SC_Score::getConcordance(double** ranking1, double** ranking2, int len, int dim) {
	int y, i, j, c;
	double concordance = 0.0;
	bool **pairs1, **pairs2; //just store the raw values if one likes to plot them...
	SC_MatrixFunctions mFunc;

	pairs1 = mFunc.initMatrix(len, dim*(dim-1)/2, false);
	pairs2 = mFunc.initMatrix(len, dim*(dim-1)/2, false);

	for (y = 0; y < len; y++) {
		c = 0;
		for (i = 0; i < dim; i++) {
			for (j = i+1; j < dim; j++) {
				pairs1[y][c] = (ranking1[y][i] < ranking1[y][j]) ? true : false; //choose ith rank over jth rank if ith value is smaller
				pairs2[y][c] = (ranking2[y][i] < ranking2[y][j]) ? true : false;
				concordance += (pairs1[y][c] == pairs2[y][c]) ? 1.0 : 0.0; //the ranking is concordant if both models give the same rank, else its discordant
				c++; //number/index of the current pair (i,j) between 0 and dim*(dim-1)/2 (overall numer of distinct pairs)
			}
		}
	}

	MFree_2D(pairs1);
	MFree_2D(pairs2);

	return concordance / (double)(len * dim*(dim-1.0)/2.0); //normaliz to the interval [0..1]
}
