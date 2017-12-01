/**************************************************************************/
/*    Responsibility:																											*/
/*      - Implements a decision-stump classifier that spearates two       */
/*        classes by one single (the best) thresholded feature            */
/*      - In the Viola and Jones AdaBoost face-detection framework this   */
/*        is called the "weak classifier"                                 */
/*      - "Robust Real-Time Face Detection", Viola, Jones, Int. J. Comp.  */
/*        Vision 57(2), 137-154, 2004                                     */
/*      - the returned probability-estimates are not proper normalized,   */
/*        so they are better viewed as confidence-scores where greater is */
/*        more confident                                                  */
/*      - TODO: Enhance multiway-splitting to a complete multiclass       */
/*              scheme; that would include managing a dynamic mapping of  */
/*              classes to indexes and vice versa                         */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.04.2007																								*/
/**************************************************************************/

#include "SC_Classifier_DecisionStump.h"
#include "SC_Aux.h"
#include "SC_DistanceMeasures.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Classifier_DecisionStump::SC_Classifier_DecisionStump(SC_TweakableParameters* pTweak, bool verbose, bool optimizeF1pos) : SC_ClassifierWithWeights(pTweak, false, verbose) { //no scaling here; we rely on just one feature!
	this->classifierType = sclib::ctDecisionStump;
	this->intervalCount = 0;
	this->feature = -1;
	this->threshold = NULL;
	this->label = NULL;
	this->confidence = NULL;
	this->margin = NULL;
	this->trainingError = -1.0;
	this->optimizeF1pos = optimizeF1pos;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Classifier_DecisionStump::~SC_Classifier_DecisionStump() {
	MFree_1D(this->threshold);
	MFree_1D(this->label);
	MFree_1D(this->confidence);
	MFree_1D(this->margin);
}

//====================================================================================================================
//  auxiliary method for indexMeregeSort(); the basic idea is from http://en.wikipedia.org/wiki/Merge_sort
//====================================================================================================================
/*void SC_Classifier_DecisionStump::indexMergeCombine(float** index, double **associatedMatrix, int rowLBound, int rowUBound, int pivot, int dataSortBy, bool sortByIndex) {
	int k, i = rowLBound, j = pivot+1;
	float *indexTemp;
	double *associatedTemp;
	
	while (j != rowUBound+1 && i != j) { //continue until either list runs out
		if ((sortByIndex==true && index[j][dataSortBy]<=index[i][dataSortBy]) ||
			  (sortByIndex==false && associatedMatrix[j][dataSortBy]<=associatedMatrix[i][dataSortBy])) { //Move the jth element in front of the ith element
			indexTemp = index[j]; //save for later insertion
			associatedTemp = associatedMatrix[j];

			for (k = j; k > i; k--) {
				index[k] = index[k-1]; //Shifts elements from i to j one step forward
				associatedMatrix[k] = associatedMatrix[k-1];
			}

			index[i] = indexTemp; //Puts final element in the right place
			associatedMatrix[i] = associatedTemp;

			i++;
			j++;
		} else {
			i++; // Skip to the next element
		}
	}

	return;
}*/

//====================================================================================================================
//  a special variant of mergesort that uses knowledge of the datasructures in trainTwoClass(): given is an index-
//  array that holds per row a pointer to the rows of the actual data-matrix (not necessarily one single matrix, rows 
//  can be distributed over many matrixes), and a matrix associated with the data and index (rows correspond). now, 
//  the index and the associated matrix get sorted while the original data is left as is. 
//  rowLBound is lower bound, rowUBound is upper bound of the index-array and data-matrix' 1st dimension (rows, from 
//  i=rowLBound;i<=rowUBound;i++). dataSortBy is the data-column which gives the reference for sorting; if 
//  sortByIndex==false, all is sorted by the dataSortBy column of the associatedMatrix instead of the index.
//  we use mergesort now because it is stable, i.e. rows with equal values remain in their previous order (needed to 
//  sort by different columns subsequently)
//====================================================================================================================
void SC_Classifier_DecisionStump::indexMergeSort(float **index, double **associatedMatrix, float **tmpIdx, double **tmpMat, int rowLBound, int rowUBound, int dataSortBy, bool sortByIndex) {
	int i, j, k, m;

	if (rowUBound > rowLBound) {
		m = rowLBound + ((rowUBound-rowLBound)/2); //middle of left and right, complicated to avoid overflow (see http://googleresearch.blogspot.com/2006/06/extra-extra-read-all-about-it-nearly.html for details)
		
		indexMergeSort(index, associatedMatrix, tmpIdx, tmpMat, rowLBound, m, dataSortBy, sortByIndex);
		indexMergeSort(index, associatedMatrix, tmpIdx, tmpMat, m+1, rowUBound, dataSortBy, sortByIndex);
		for (i = m+1; i > rowLBound; i--) {
			tmpIdx[i-1] = index[i-1];
			tmpMat[i-1] = associatedMatrix[i-1];
		}
		for (j = m; j < rowUBound; j++) {
			tmpIdx[rowUBound+m-j] = index[j+1];
			tmpMat[rowUBound+m-j] = associatedMatrix[j+1];
		}
		if (sortByIndex == true) {
			for (k = rowLBound; k <= rowUBound; k++) {
				if (tmpIdx[i][dataSortBy] < tmpIdx[j][dataSortBy]) {
					index[k] = tmpIdx[i];
					associatedMatrix[k] = tmpMat[i++];
				} else {
					index[k] = tmpIdx[j];
					associatedMatrix[k] = tmpMat[j--];
				}
			}
		} else {
			for (k = rowLBound; k <= rowUBound; k++) {
				if (tmpMat[i][dataSortBy] < tmpMat[j][dataSortBy]) {
					index[k] = tmpIdx[i];
					associatedMatrix[k] = tmpMat[i++];
				} else {
					index[k] = tmpIdx[j];
					associatedMatrix[k] = tmpMat[j--];
				}
			}
		}
	}

	return;
	
	/*	int pivot = rowLBound + ((rowUBound-rowLBound)/2); //middle of left and right, complicated to avoid overflow (see http://googleresearch.blogspot.com/2006/06/extra-extra-read-all-about-it-nearly.html for details)

	if (rowLBound != rowUBound) {
		indexMergeSort(index, associatedMatrix, rowLBound, pivot, dataSortBy); //First half
		indexMergeSort(index, associatedMatrix, pivot+1, rowUBound, dataSortBy); //Second half
		//indexMergeCombine(index, associatedMatrix, rowLBound, rowUBound, pivot, dataSortBy, sortByIndex); //Combines two sorted arrays (left -> pivot) and (pivot+1 -> right)

		int i = rowLBound, j = pivot+1;
		float *indexTemp;
		double *associatedTemp;

		if (sortByIndex == true) { //write the algorithm twice, onece for each version of the comparison statement, so that the inner loop is as fast as possible
			while (j != rowUBound+1 && i != j) { //continue until either list runs out
				if (index[j][dataSortBy] <= index[i][dataSortBy]) { //Move the jth element in front of the ith element
					indexTemp = index[j]; //save for later insertion
					associatedTemp = associatedMatrix[j];

					for (int k = j; k > i; k--) {
						index[k] = index[k-1]; //Shifts elements from i to j one step forward
						associatedMatrix[k] = associatedMatrix[k-1];
					}

					index[i] = indexTemp; //Puts final element in the right place
					associatedMatrix[i] = associatedTemp;

					i++;
					j++;
				} else {
					i++; // Skip to the next element
				}
			}
		} else {
			while (j != rowUBound+1 && i != j) { //continue until either list runs out
				if (associatedMatrix[j][dataSortBy] <= associatedMatrix[i][dataSortBy]) { //Move the jth element in front of the ith element
					indexTemp = index[j]; //save for later insertion
					associatedTemp = associatedMatrix[j];

					for (int k = j; k > i; k--) {
						index[k] = index[k-1]; //Shifts elements from i to j one step forward
						associatedMatrix[k] = associatedMatrix[k-1];
					}

					index[i] = indexTemp; //Puts final element in the right place
					associatedMatrix[i] = associatedTemp;

					i++;
					j++;
				} else {
					i++; // Skip to the next element
				}
			}
		}
	}

	return;*/
}

//====================================================================================================================
//  a special variant of quicksort that uses knowledge of the datasructures in trainTwoClass(): given is an index-
//  array that holds per row a pointer to the rows of the actual data-matrix (not necessarily one single matrix, rows 
//  can be distributed over many mazrixes), and a matrix associated with the data and index (rows correspond). now, 
//  the index and the associated matrix get sorted while the original data is left as is. 
//  rowLBound is lower bound, rowUBound is upper bound of the index-array and data-matrix' 1st dimension (rows, from 
//  i=rowLBound;i<=rowUBound;i++).dataSortBy is the data-column which gives the reference for sorting
//====================================================================================================================
void SC_Classifier_DecisionStump::indexQuickSort(float** index, double **associatedMatrix, int rowLBound, int rowUBound, int dataSortBy) {
	int i, j;
	float x; //storage for comparison-element
	double *associatedTemp; //storage for temp-element for triangle-exchange
	float *indexTemp;
	
	if (rowLBound >= rowUBound) { //nothing to do here
    return;
  } 

	i=rowLBound;
	j=rowUBound;
	x = index[(rowLBound+rowUBound)/2][dataSortBy];

	//start of division-step
	while (i <= j) {
		while (index[i][dataSortBy] < x) {i++;}
		while (index[j][dataSortBy] > x) {j--;}
		if (i<=j)	{ //do the triangle-exchange
			indexTemp = index[i]; //exchange the indixes
			index[i] = index[j];
			index[j] = indexTemp;

			associatedTemp = associatedMatrix[i]; //exchange the associated matrix
			associatedMatrix[i] = associatedMatrix[j];
			associatedMatrix[j] = associatedTemp;

			i++; 
			j--;
		}
	}

	indexQuickSort(index, associatedMatrix, rowLBound, j, dataSortBy);
	indexQuickSort(index, associatedMatrix, i, rowUBound, dataSortBy);

	return;
}

//====================================================================================================================
//	called by trainTwoClass() to find a single (2-way) split that minimizes the training error
//  the method receives meta-data (with T rows of dimensionality dim, describing both the positive and negative 
//  examples in their original matrixes) of the form [label|weight|posWeightSum|negWeightSum] and a corresponding 
//  index-array with pointers (in it's single-column rows) to the original untouched two datasets.
//  the result is directly inserted into the class-members label, threshold, feature and confidence, the return-value 
//  gives the training error (-1.0 to indicat error)
//====================================================================================================================
double SC_Classifier_DecisionStump::findSingleBestSplit(double **metaData, float **index, int T, int dim, double totalPositives, double totalNegatives) {
	int t, d, t1, dim4 = dim+4, T1 = T-1, dim1 = dim-1;
	double thresh[3], error, lowPosError, lowNegError, minError = 1.0, lastPercentage = 0.0, currentMargin = 0.0, maxMargin = std::numeric_limits<double>::min();
	bool lowPos; //label the lower part positive or not?
	double positives[2], negatives[2], eps = 1.0 / (double)(T); //epsilon to smooth confidence scores as proposed in Schapire, Singer, "Improved Boosting Algorithms Using Confidence-rated Predictions", Machine Learning, 37(3): 297-336, 1999
	//double score, minScore = std::numeric_limits<double>::max();
	double **tmpMeta;
	float **tmpIdx;

	//new parameters
  this->intervalCount = 2; //lower bound -> thresh, thresh -> upper bound
	MArray_1D(this->label, this->intervalCount, int, "SC_Classifier_DecisionStump.findSingleBestSplit: label");
	MArray_1D(this->threshold, this->intervalCount+1, double, "SC_Classifier_DecisionStump.findSingleBestSplit: threshold"); //threshs + upper bound (max) + lower bound (min)
	MArray_1D(this->confidence, this->intervalCount, double, "SC_Classifier_DecisionStump.findSingleBestSplit: confidence");
	MArray_1D(this->margin, this->intervalCount-1, double, "SC_Classifier_DecisionStump.findSingleBestSplit: margin"); //one for each inner thresh, i.e. 1

	MArray_1D(tmpMeta, T, double*, "SC_Classifier_DecisionStump.findSingleBestSplit: tmpMeta");
	MArray_1D(tmpIdx, T, float*, "SC_Classifier_DecisionStump.findSingleBestSplit: tmpIdx");
	indexMergeSort(index, metaData, tmpIdx, tmpMeta, 0, T1, 0, false); //pre-sort the two matrixes by the label so that examples with equal label and equal value are in a block later on

	//search for the best feature, i.e. the one that provides lowest class prediction error when classifiying with a single threshold
	for (d = 0; d < dim; d++) {
		//this is a version of quicksort that operates on the two original datasets via the united index, and also sorts the data-matrix accordingly; the original data is untouched while this saves lots of memory
		indexMergeSort(index, metaData, tmpIdx, tmpMeta, 0, T1, d, true); //sort in ascending order by the single feature currently to evaluate

		error = 1.0; //initialize minimal error
		thresh[0] = index[0][d];//min. or lower bound for current feature
		thresh[2] = index[T1][d]; //max. or upper bound for current feature

		//do first loop (t==0) extra to save if()-clauses inside the inner loop... this saves a lot of statements here, too
		//semantics: put the split *before* the first entry
		metaData[0][2] = 0.0;
		metaData[0][3] = 0.0;
		//lowPos = totalPositives < totalNegatives; //to which class does the lower interval belong (true => positives)
		//error = (lowPos == true) ? totalPositives : totalNegatives; //error on training data when threshold is put here
		//thresh[1] = index[0][d];

		//just in case the best thresh (only happens if data is class-homogenious) is before the first example:
		thresh[1] = thresh[0];
		currentMargin = thresh[2] - thresh[0];
		lowPos = false;
		positives[0] = 0.0; //positive weights in lower half
		positives[1] = totalPositives; //positive weights in upper half
		negatives[0] = 0.0; //negative weights in lower half
		negatives[1] = totalNegatives; //negative weights in upper half
		
		//now the complete loop
		t1 = 0;
		for (t = 1; t < T; t++) { //optimized threshold position: candidates are halfway between last und current example
			if (metaData[t1][0] == (double)(sclib::labelPositive)) {
				metaData[t][2] = metaData[t1][2] + metaData[t1][1]; //sum of so far seen positive weights (without the current example)
				metaData[t][3] = metaData[t1][3];										// -- " --           negative
			} else {
				metaData[t][2] = metaData[t1][2];
				metaData[t][3] = metaData[t1][3] + metaData[t1][1];
			}

			if (metaData[t1][0] != metaData[t][0] && index[t1][d] != index[t][d]) { //only look for a potential split point between examples of different class as suggested by witten&frank (pp. 300) for the multiway-split case, and also don't split between examples with the same value
				lowPosError = metaData[t][3] + (totalPositives - metaData[t][2]); //error of labeling all examples before this one (excluding this one) positive and the following (including the current) ones negative = nr. of negatives in lower part + nr. of positives in upper part
				lowNegError = metaData[t][2] + (totalNegatives - metaData[t][3]); //labeling vice versa

				if (sclib::min(lowPosError, lowNegError) < error) { //if this is the best position for the threshold (between current and previous example) seen so far...
					lowPos = lowPosError < lowNegError; //to which class does the lower interval belong (true => positives)
					thresh[1] = (index[t][d] + index[t1][d]) / 2.0; //the thresh is halfway between this and the last example
					currentMargin = thresh[1] - index[t1][d]; //the margin of the threshold
					error = (lowPos == true) ? lowPosError : lowNegError; //error on training data when threshold is put here
					positives[0] = metaData[t][2]; //positive weights in lower half
					positives[1] = totalPositives - metaData[t][2]; //positive weights in upper half
					negatives[0] = metaData[t][3]; //negative weights in lower half
					negatives[1] = totalNegatives - metaData[t][3]; //negative weights in upper half
				}
			}
			
			t1++; //t1==t-1 to speed up things a little bit...
		}

		//score = sclib::scaleToInterval(error, 0.0, 0.5, 0.0001, 1.0) * 
		//				sclib::scaleToInterval((thresh[2]-thresh[0])/currentMargin, 2.0, 5000.0, 0.01, 1.0); //invert and normalize margin to make it comparabel across different possibly unscaled features
		//if (score < minScore) {
			//minScore = score;
		if ((error < minError) ||
			  (error==minError && currentMargin>maxMargin)) { //the current feature performed better than all previous single features
			minError = error;
			maxMargin = currentMargin;
			this->feature = d;
			this->threshold[0] = thresh[0];
			this->threshold[1] = thresh[1];
			this->threshold[2] = thresh[2];
			if (lowPos == true) {
				this->label[0] = sclib::labelPositive;
				this->label[1] = sclib::labelNegative;
			} else {
				this->label[0] = sclib::labelNegative;
				this->label[1] = sclib::labelPositive;
			}
			this->confidence[0] = fabs(0.5 * sclib::ln((positives[0]+eps) / (negatives[0]+eps))); //smooth the confidence estimates as proposed by Schapire&Singer
			this->confidence[1] = fabs(0.5 * sclib::ln((positives[1]+eps) / (negatives[1]+eps))); //take the absolute value here 'cause the label (which equates to the sign of this expression) is already encoded in the label-variable
			this->margin[0] = currentMargin;
		}

		if (this->verbose == true) {
			lastPercentage = sclib::printPercentage(dim1, d, lastPercentage, 1.0, d==0 ? true : false);
		}
	}

	MFree_1D(tmpMeta);
	MFree_1D(tmpIdx);

	return minError;
}

//====================================================================================================================
//	as above, but with the F1 measure on the combined recall&precision of the positive class as the optimization
//  criterion, in contrast to mere classification error. an additional constarint here: all samples above the
//  threshold will be regarded as negative (opposed to: the predicted class of an interval is assumed to be the
//  majority class within => this method is best used to compute the lambda penalty factor for BIC computations). 
//  the best F1 is returned.
//====================================================================================================================
double SC_Classifier_DecisionStump::findSingleBestSplitF1pos(double **metaData, float **index, int T, int dim, double totalPositives, double totalNegatives) {
	int t, d, t1, dim4 = dim+4, T1 = T-1, dim1 = dim-1;
	double thresh[3], lastPercentage = 0.0, currentMargin = 0.0, maxMargin = std::numeric_limits<double>::min();
	double positives[2], negatives[2], eps = 1.0 / (double)(T); //epsilon to smooth confidence scores as proposed in Schapire, Singer, "Improved Boosting Algorithms Using Confidence-rated Predictions", Machine Learning, 37(3): 297-336, 1999
	double **tmpMeta;
	float **tmpIdx;
	double recall, precision, bestF1 = 0.0, maxF1, f1; //, correct, found, gt = totalPositives;

	//new parameters
  this->intervalCount = 2; //lower bound -> thresh, thresh -> upper bound
	MArray_1D(this->label, this->intervalCount, int, "SC_Classifier_DecisionStump.findSingleBestSplitF1pos: label");
	MArray_1D(this->threshold, this->intervalCount+1, double, "SC_Classifier_DecisionStump.findSingleBestSplitF1pos: threshold"); //threshs + upper bound (max) + lower bound (min)
	MArray_1D(this->confidence, this->intervalCount, double, "SC_Classifier_DecisionStump.findSingleBestSplitF1pos: confidence");
	MArray_1D(this->margin, this->intervalCount-1, double, "SC_Classifier_DecisionStump.findSingleBestSplitF1pos: margin"); //one for each inner thresh, i.e. 1

	MArray_1D(tmpMeta, T, double*, "SC_Classifier_DecisionStump.findSingleBestSplitF1pos: tmpMeta");
	MArray_1D(tmpIdx, T, float*, "SC_Classifier_DecisionStump.findSingleBestSplitF1pos: tmpIdx");
	indexMergeSort(index, metaData, tmpIdx, tmpMeta, 0, T1, 0, false); //pre-sort the two matrixes by the label so that examples with equal label and equal value are in a block later on

	//search for the best feature, i.e. the one that provides lowest class prediction error when classifiying with a single threshold
	for (d = 0; d < dim; d++) {
		//this is a version of quicksort that operates on the two original datasets via the united index, and also sorts the data-matrix accordingly; the original data is untouched while this saves lots of memory
		indexMergeSort(index, metaData, tmpIdx, tmpMeta, 0, T1, d, true); //sort in ascending order by the single feature currently to evaluate

		maxF1 = 0.0; //initialize maximal f1
		thresh[0] = index[0][d];//min. or lower bound for current feature
		thresh[2] = index[T1][d]; //max. or upper bound for current feature

		//do first loop (t==0) extra to save if()-clauses inside the inner loop... this saves a lot of statements here, too
		//semantics: put the split *before* the first entry
		metaData[0][2] = 0.0;
		metaData[0][3] = 0.0;

		//just in case the best thresh (only happens if data is class-homogenious) is before the first example:
		thresh[1] = thresh[0];
		currentMargin = thresh[2] - thresh[0];
		positives[0] = 0.0; //positive weights in lower half
		positives[1] = totalPositives; //positive weights in upper half
		negatives[0] = 0.0; //negative weights in lower half
		negatives[1] = totalNegatives; //negative weights in upper half
		
		//now the complete loop
		t1 = 0;
		for (t = 1; t < T; t++) { //optimized threshold position: candidates are halfway between last und current example
			if (metaData[t1][0] == (double)(sclib::labelPositive)) {
				metaData[t][2] = metaData[t1][2] + metaData[t1][1]; //sum of so far seen positive weights (without the current example)
				metaData[t][3] = metaData[t1][3];										// -- " --           negative
			} else {
				metaData[t][2] = metaData[t1][2];
				metaData[t][3] = metaData[t1][3] + metaData[t1][1];
			}

			if (metaData[t1][0] != metaData[t][0] && index[t1][d] != index[t][d]) { //only look for a potential split point between examples of different class as suggested by witten&frank (pp. 300) for the multiway-split case, and also don't split between examples with the same value
				//lowPosError = metaData[t][3] + (totalPositives - metaData[t][2]); //error of labeling all examples before this one (excluding this one) positive and the following (including the current) ones negative = nr. of negatives in lower part + nr. of positives in upper part
				//lowNegError = metaData[t][2] + (totalNegatives - metaData[t][3]); //labeling vice versa
				//found = metaData[t][2] + metaData[t][3]; //sum of all weights that will be ascribed the positive label with this threshold: everything before the current sample
				//correct = metaData[t][2]; //summ of all weights that will correctly be ascribed the positive label: all positive already seen
				recall = metaData[t][2] / totalPositives; //correct / gt;
				precision = metaData[t][2] / (metaData[t][2] + metaData[t][3]); //correct / found;
				f1 = (recall+precision > 0.0) ? ((2.0*recall*precision) / (recall+precision)) : 0.0;

				if (f1 > maxF1) { //if this is the best position for the threshold (between current and previous example) seen so far...
					thresh[1] = (index[t][d] + index[t1][d]) / 2.0; //the thresh is halfway between this and the last example
					currentMargin = thresh[1] - index[t1][d]; //the margin of the threshold
					maxF1 = f1; //F1 on positive class on training data when threshold is put here
					positives[0] = metaData[t][2]; //positive weights in lower half
					positives[1] = totalPositives - metaData[t][2]; //positive weights in upper half
					negatives[0] = metaData[t][3]; //negative weights in lower half
					negatives[1] = totalNegatives - metaData[t][3]; //negative weights in upper half
				}
			}
			
			t1++; //t1==t-1 to speed up things a little bit...
		}

		if ((maxF1 > bestF1) ||
			  (maxF1==bestF1 && currentMargin>maxMargin)) { //the current feature performed better than all previous single features
			bestF1 = maxF1; //bestF1 is globally for all feature-columns, maxF1 is maximum of current column
			maxMargin = currentMargin;
			this->feature = d;
			this->threshold[0] = thresh[0];
			this->threshold[1] = thresh[1];
			this->threshold[2] = thresh[2];
			this->label[0] = sclib::labelPositive;
			this->label[1] = sclib::labelNegative;
			this->confidence[0] = fabs(0.5 * sclib::ln((positives[0]+eps) / (negatives[0]+eps))); //smooth the confidence estimates as proposed by Schapire&Singer
			this->confidence[1] = fabs(0.5 * sclib::ln((positives[1]+eps) / (negatives[1]+eps))); //take the absolute value here 'cause the label (which equates to the sign of this expression) is already encoded in the label-variable
			this->margin[0] = currentMargin;
		}

		if (this->verbose == true) {
			lastPercentage = sclib::printPercentage(dim1, d, lastPercentage, 1.0, d==0 ? true : false);
		}
	}

	MFree_1D(tmpMeta);
	MFree_1D(tmpIdx);

	return bestF1;
}

//====================================================================================================================
//  called by findMultiwaySplit() to recursively split the given data-set (see calling method for parameter info) as 
//  long as the MDL criterion promotes the found splits; returns the training error on found splitting, or -1.0 on 
//  error/failure to find split point; details for the found split points are returned in the linked list of 
//  thresholds.
//====================================================================================================================
double SC_Classifier_DecisionStump::splitInterval(double **metaData, float **index, int startRow, int stopRow, int column, double totalPositives, double totalNegatives, SC_Classifier_DecisionStump::SC_ThresholdList* &thresholdList) {
	double thresh, noSplitInfo = SC_DistanceMeasures::information(totalPositives, totalNegatives), lowerInfo, upperInfo, splitInfo, infoGain, maxGain = 0.0;
	double remainingPos, remainingNeg, allWeights = totalPositives+totalNegatives;
	double error = -1.0, lowerError, upperError;
	int t, t1, splitIdx, N = stopRow-startRow+1, k1, k2;
	bool lowPos, upPos;
	SC_Classifier_DecisionStump::SC_ThresholdList *result = NULL;
	
	if (stopRow <= startRow || totalPositives == 0.0 || totalNegatives == 0.0) { //data-set is just 1 row in size or class-homogenious
		return -1.0;
	}

	//do first loop (t==startRow) extra to save if()-clauses inside the inner loop... this saves a lot of statements here, too
	//semantics: put the split *before* the first entry
	//what if this is the best split (happens only if the data is class homogenious...)?
	// => can't happen 'cause then totalPositives or totalNegatives would be zero and the methods would have returned already with a value of -1.0
	metaData[startRow][2] = 0.0;
	metaData[startRow][3] = 0.0;
	//lowPos = totalPositives < totalNegatives; //to which class does the lower interval (i.e. everything below the startRow) belong (true => positives)
	//splitInfo = noSplitInfo; //the information when setting a split just before the first example is the same as the information with no split at all
	//thresh = index[startRow][column];

	//now the complete loop
	t1 = startRow;
	for (t = startRow+1; t <= stopRow; t++) { //optimized threshold position: candidates are halfway between last und current example
		if (metaData[t1][0] == (double)(sclib::labelPositive)) { //last example was positive
			metaData[t][2] = metaData[t1][2] + metaData[t1][1]; //sum of so far seen positive weights (without the current example)
			metaData[t][3] = metaData[t1][3];										// -- " --           negative
		} else {
			metaData[t][2] = metaData[t1][2];
			metaData[t][3] = metaData[t1][3] + metaData[t1][1];
		}

		if (metaData[t1][0] != metaData[t][0] && index[t1][column] != index[t][column]) { //only look for a possible split between examples with different labels to speed up as suggested in witten&frank, pp. 300; the second test is to avoid putting a threshold between examples with the same value
			remainingPos = sclib::max(totalPositives - metaData[t][2], 0.0); //positives in upper part (including the current example); max() to avoid negative numbers due to rounding errors
			remainingNeg = sclib::max(totalNegatives - metaData[t][3], 0.0); //negatives -- " --
			lowerInfo = SC_DistanceMeasures::information(metaData[t][2], metaData[t][3]); //information (or entropy) in the lower part
			upperInfo = SC_DistanceMeasures::information(remainingPos, remainingNeg); //information (or entropy) in the upper part
			splitInfo = ((metaData[t][2]+metaData[t][3])/allWeights)*lowerInfo + ((remainingPos+remainingNeg)/allWeights)*upperInfo; //information of both parts
			infoGain = noSplitInfo - splitInfo; //information gain := difference between information without split and information with split here (witten&frank, pp. 299)
			k1 = (metaData[t][2]==0.0 || metaData[t][3]==0.0) ? 1 : 2; //number of distinct classes in lower part
			k2 = (remainingPos==0.0 || remainingNeg==0.0) ? 1 : 2; //number of distinct classes in upper part

			if (infoGain > maxGain && SC_DistanceMeasures::splittingMDL(infoGain, N, 2, k1, k2, noSplitInfo, lowerInfo, upperInfo) == true) { //we tested before that there are both classes available in the complete data-set
				maxGain = infoGain;
				thresh = (index[t][column] + index[t1][column]) / 2.0; //the thresh is halfway between this and the last example
				lowPos = metaData[t][2] > metaData[t][3]; //shall the lower part get the positive label due to a majority of the positive class therein?
				upPos = remainingPos > remainingNeg; //both parts can get the same label but may be splitted due to distribution changes (witten&frank, pp. 300)
				splitIdx = t; //split is *before* this example

				if (result == NULL) {
					result = new SC_Classifier_DecisionStump::SC_ThresholdList(thresh, infoGain, lowPos?sclib::labelPositive:sclib::labelNegative, metaData[t][2], metaData[t][3], upPos?sclib::labelPositive:sclib::labelNegative, remainingPos, remainingNeg, index[t][column]-thresh);
				} else {
					result->setValues(thresh, infoGain, lowPos?sclib::labelPositive:sclib::labelNegative, metaData[t][2], metaData[t][3], upPos?sclib::labelPositive:sclib::labelNegative, remainingPos, remainingNeg, index[t][column]-thresh);
				}
			}
		}

		t1++;
	}

	//do recursion on upper and lower part to see if further spilts are needed if a valid split was found previously
	if (result != NULL) {
		if (thresholdList != NULL) { //connect the current result with the result-list or start one
			sclib::getLastInList(thresholdList)->Next = result; 
		} else {
			thresholdList = result;
		}

		//recursion on lower part
		lowerError = splitInterval(metaData, index, startRow, splitIdx-1, column, metaData[splitIdx][2], metaData[splitIdx][3], thresholdList);
		if (lowerError >= 0.0) {
			result->lowerPartLabel = 0; //"delete" label if it is further specified by a child-split
		} else {
			lowerError = (result->lowerPartLabel == sclib::labelPositive) ? result->lowerPartClassCount[label2idx(sclib::labelNegative)] : result->lowerPartClassCount[label2idx(sclib::labelPositive)]; //the error on the upper half of the data set when just the current split is considered
		}

		//recursion on upper part
		upperError = splitInterval(metaData, index, splitIdx, stopRow, column, totalPositives-metaData[splitIdx][2], totalNegatives-metaData[splitIdx][3], thresholdList);
		if (upperError >= 0.0) {
			result->upperPartLabel = 0; //"delete" label if it is further specified by a child-split
		} else {
			upperError = (result->upperPartLabel == sclib::labelPositive) ? result->upperPartClassCount[label2idx(sclib::labelNegative)] : result->upperPartClassCount[label2idx(sclib::labelPositive)]; //the error on the upper half of the data set when just the current split is considered
		}

		error = lowerError + upperError; //the missclassified weights when considerung this- and all child splits
	}

	return error; //true if a valid splitting point was found
}

//====================================================================================================================
//  takes a list of threshold-facts and returns the complete margin of the classifier composed of all the thresholds 
//  in it as a weighted average (weighted by number of instances in the corresponding intervals) of the individual 
//  margins (scaled so that margins on different features are comparable despite the feature where scaled previously 
//  or not)
//====================================================================================================================
double SC_Classifier_DecisionStump::getCompleteMargin(SC_Classifier_DecisionStump::SC_ThresholdList *pThreshList, int listCount, double minValue, double maxValue) {
	int i;
	double completeMargin = 0.0, nrOfExamples = 0.0;
	double greaterThan, smallerThan, *margins, *examples, weight;
	SC_Classifier_DecisionStump::SC_ThresholdList *pHook = pThreshList, *pCurrent = NULL;

	MArray_1D(margins, listCount, double, "SC_ClassifierDecisionStump.getCompleteMargin: margins");
	MArray_1D(examples, listCount+1, double, "SC_ClassifierDecisionStump.getCompleteMargin: examples");

	for (i = 0; i < listCount; i++) {
		//search for the ith smallest threshold in the list by finding the minimal threshold that is greater than the last one
		greaterThan = (i == 0) ? minValue-std::numeric_limits<double>::epsilon() : pCurrent->threshold;
		smallerThan = maxValue;
		pHook = pThreshList;
		while (pHook != NULL) {
			if (pHook->threshold > greaterThan && pHook->threshold < smallerThan) {
				smallerThan = pHook->threshold; //to find the minimal threshold that exceeds greaterThan
				pCurrent = pHook;
			}
			pHook = pHook->Next;
		}

		//to make completeMargins on different features (that have not been scaled to a common interval beforehand) comparable, we can just
		//scale the margins directly by dividing them by the original dynamic range (this is the same as scaoling all features by 
		//(x-min)/(max-min) before threshold-learning).
		margins[i] = pCurrent->margin / (maxValue - minValue); //i.e. no standardization (-min) necessary because it's a distance, hence only normalization (/(max-min))

		//fill number of examples variable
		if (pCurrent->lowerPartLabel != 0) { //a 0 here indicates that there is a more detailed (child) subinterval that further splits this lower part and therefore holds the example-counts
			//examples[i] holds the number of examples in the interval just below the threshold
			examples[i] = pCurrent->lowerPartClassCount[label2idx(sclib::labelPositive)] + pCurrent->lowerPartClassCount[label2idx(sclib::labelNegative)];
			nrOfExamples += examples[i];
		}
		if (pCurrent->upperPartLabel != 0) { //as above
			examples[i+1] = pCurrent->upperPartClassCount[label2idx(sclib::labelPositive)] + pCurrent->upperPartClassCount[label2idx(sclib::labelNegative)];
			nrOfExamples += examples[i+1];
		}
	} //for i

	//compute the weighted average
	for (i = 0; i < listCount; i++) {
		if (listCount == 1) {
			weight = examples[0] + examples[1];
		} else if (i == 0) {
			weight = examples[0] + examples[1]/2.0;
		} else if (i == listCount-1) {
			weight = examples[i]/2.0 + examples[listCount];
		} else {
			weight = (examples[i] + examples[i+1]) / 2.0; //all but the first an last interval are counted twice (as an upper and lower part, respectively), thus the division by 2
		}
		completeMargin += margins[i] * weight;
	}
	completeMargin /= nrOfExamples;
		
	MFree_1D(margins);
	MFree_1D(examples);

	return completeMargin;
}

//====================================================================================================================
//	called by train*Class() to find the multiway-split that minimizes the entropy
//  the method receives meta-data (with T rows of dimensionality dim, discribing both the positive and negative 
//  examples in their original matrixes) of the form [label|weight|posWeightSum|negWeightSum] and a corresponding 
//  index-array with pointers (int it's single-column rows) to the original untouched two datasets.
//  the result is directly inserted into the class-members label, threshold, feature and confidence, the return-value 
//  gives the training error (-1.0 to indicat error)
//====================================================================================================================
double SC_Classifier_DecisionStump::findMultiwaySplit(double **metaData, float **index, int T, int dim, double totalPositives, double totalNegatives) {
	SC_Classifier_DecisionStump::SC_ThresholdList *thresholdList, *bestThresholdList = NULL, *currentThreshold;
	double error, minError = 1.0, bestMinThresh, bestMaxThresh, lastPercentage = 0.0, greaterThan, smallerThan, maxMargin = std::numeric_limits<double>::min();
	int bestColumn, T1 = T-1, dim1 = dim-1, thresholdCount, minThresholdCount = T;
	double eps = 1.0 / (double)(T); //epsilon to smooth confidence scores as proposed in Schapire, Singer, "Improved Boosting Algorithms Using Confidence-rated Predictions", Machine Learning, 37(3): 297-336, 1999
	char debugOutput[sclib::bufferSize]; //for debug output only
	double *debugPos = NULL, *debugNeg = NULL; //for debug output only
	double completeMargin; //, score, minScore = std::numeric_limits<double>::max();
	double **tmpMeta;
	float **tmpIdx;

	MArray_1D(tmpMeta, T, double*, "SC_Classifier_DecisionStump.findMultiwaySplit: tmpMeta");
	MArray_1D(tmpIdx, T, float*, "SC_Classifier_DecisionStump.findMultiwaySplit: tmpIdx");
	indexMergeSort(index, metaData, tmpIdx, tmpMeta, 0, T1, 0, false); //pre-sort the two matrixes by the label so that examples with equal label and equal value are in a block later on

	//search for the best feature, i.e. the one that provides lowest class prediction error when classifiying with a single threshold
	for (int d = 0; d < dim; d++) {
		//this is a version of quicksort that operates on the two original datasets via the united index, and also sorts the data-matrix accordingly; the original data is untouched while this saves lots of memory
		indexMergeSort(index, metaData, tmpIdx, tmpMeta, 0, T1, d, true); //sort in ascending order by the single feature currently to evaluate

		//sclib::matrixOut("index.txt", index, T, dim, this->pTweak, 0, 0, 0, 1);
		//sclib::matrixOut("metaData.txt", metaData, T, 4, this->pTweak, 0, 0, 0, 2);

		thresholdList = NULL;
		error = splitInterval(metaData, index, 0, T1, d, totalPositives, totalNegatives, thresholdList);
		if (error >= 0.0) {
			thresholdCount = sclib::getListCount(thresholdList);
			completeMargin = getCompleteMargin(thresholdList, thresholdCount, index[0][d], index[T1][d]);

			//score = sclib::scaleToInterval(error, 0.0, 0.5, 0.01, 1.0) * 
			//	      sclib::scaleToInterval((double)(thresholdCount), 1.0, 10.0, 0.001, 1.0) * 
			//				sclib::scaleToInterval(1.0/completeMargin, 2.0, 5000.0, 0.0001, 1.0);
			//if (score < minScore) {
				//minScore = score;
			if (error >= 0.0 && ((error < minError) || 
													(error==minError && thresholdCount<minThresholdCount) || 
													(error==minError && thresholdCount==minThresholdCount && completeMargin>maxMargin))) { //error==-1.0 indicates failure in making any split on this feature due to MDL criterion
				minError = error;
				minThresholdCount = thresholdCount; //we want minimum error achieved with a minimal count of thresholds
				maxMargin = completeMargin;
				sclib::destructLinkedList(bestThresholdList);
				bestThresholdList = thresholdList;
				bestMinThresh = index[0][d]; //min. or lower bound for current feature
				bestMaxThresh = index[T1][d]; //max. or upper bound for current feature
				bestColumn = d;
				
				//if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
				//	sclib::matrixOut("bestFeature.txt", index, T, dim, this->pTweak, 0, 0, d, d+1);
				//	sclib::matrixOut("bestLabels.txt", metaData, T, 4, this->pTweak, 0, 0, 0, 1);
				//}
			} else {
				sclib::destructLinkedList(thresholdList);
			}
		} else {
			sclib::destructLinkedList(thresholdList);
		}

		if (this->verbose == true) {
			lastPercentage = sclib::printPercentage(dim1, d, lastPercentage, 1.0, d==0 ? true : false);
		}
	}

	MFree_1D(tmpMeta);
	MFree_1D(tmpIdx);

	if (bestThresholdList != NULL) { //splitting is found, decode the tree in the linked result-list to set class parameters
		this->intervalCount = sclib::getListCount(bestThresholdList) + 1; //lower bound -> thresh1, thresh1 -> thresh2, ...,  threshN -> upper bound
		
		MArray_1D(this->label, this->intervalCount, int, "SC_Classifier_DecisionStump.findMultiwaySplit: label");
		MArray_1D(this->threshold, this->intervalCount+1, double, "SC_Classifier_DecisionStump.findMultiwaySplit: threshold"); //threshs + upper bound (max) + lower bound (min)
		MArray_1D(this->confidence, this->intervalCount, double, "SC_Classifier_DecisionStump.findMultiwaySplit: confidence");
		MArray_1D(this->margin, this->intervalCount-1, double, "SC_Classifier_DecisionStump.findMultiwaySplit: margin"); //one for each inner threshold

		this->threshold[0] = bestMinThresh;
		this->threshold[this->intervalCount] = bestMaxThresh;
		this->feature = bestColumn;

		if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
			MArray_1D(debugPos, this->intervalCount, double, "SC_Classifier_DecisionStump.findMultiwaySplit: debugPos");
			MArray_1D(debugNeg, this->intervalCount, double, "SC_Classifier_DecisionStump.findMultiwaySplit: debugNeg");
			sprintf(debugOutput, "\n%s \t %s \t %s \t %s", "thresh", "lower label", "lower positives", "lower negatives");
			sclib::stringOut("multiwaySplit.txt", debugOutput, this->pTweak, "\n");
			sclib::scalarOut("multiwaySplit.txt", this->threshold[0], this->pTweak, true, "\n");
		}

		//decode the thresholdList
		for (int i = 1; i < this->intervalCount; i++) {
			//search for the ith smallest threshold in the list by finding the minimal threshold that is greater than the last one
			greaterThan = this->threshold[i-1];
			smallerThan = this->threshold[this->intervalCount];
			thresholdList = bestThresholdList;
			while (thresholdList != NULL) {
				if (thresholdList->threshold > greaterThan && thresholdList->threshold < smallerThan) {
					smallerThan = thresholdList->threshold; //to find the minimal threshold that exceeds greaterThan
					currentThreshold = thresholdList;
				}
				thresholdList = thresholdList->Next;
			}

			//fill ith class members
			this->threshold[i] = currentThreshold->threshold;
			this->margin[i-1] = currentThreshold->margin; //ith threshold corresponds with i-1th margin
			if (currentThreshold->lowerPartLabel != 0) { //a 0 here indicates that there is a more detailed (child) subinterval that further splits this lower part and therefore holds the label and confidence measures
				this->label[i-1] = currentThreshold->lowerPartLabel; //label for the interval below the threshold
				this->confidence[i-1] = fabs(0.5 * sclib::ln((currentThreshold->lowerPartClassCount[label2idx(sclib::labelPositive)]+eps) / (currentThreshold->lowerPartClassCount[label2idx(sclib::labelNegative)]+eps))); //smooth the confidence estimates as proposed by Schapire&Singer
				if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
					debugPos[i-1] = currentThreshold->lowerPartClassCount[label2idx(sclib::labelPositive)];
					debugNeg[i-1] = currentThreshold->lowerPartClassCount[label2idx(sclib::labelNegative)];
				}
			}
			if (currentThreshold->upperPartLabel != 0) { //as above
				this->label[i] = currentThreshold->upperPartLabel; //label for the interval above (and including) the threshold
				this->confidence[i] = fabs(0.5 * sclib::ln((currentThreshold->upperPartClassCount[label2idx(sclib::labelPositive)]+eps) / (currentThreshold->upperPartClassCount[label2idx(sclib::labelNegative)]+eps))); //smooth the confidence estimates as proposed by Schapire&Singer
				if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
					debugPos[i] = currentThreshold->upperPartClassCount[label2idx(sclib::labelPositive)];
					debugNeg[i] = currentThreshold->upperPartClassCount[label2idx(sclib::labelNegative)];
				}
			}

			if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
				sprintf(debugOutput, "%f \t %d \t %f \t %f", this->threshold[i], this->label[i-1], debugPos[i-1], debugNeg[i-1]);
				sclib::stringOut("multiwaySplit.txt", debugOutput, this->pTweak, "\n");
			}
		} //for i

		if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
			sprintf(debugOutput, "%f \t %d \t %f \t %f", this->threshold[this->intervalCount], this->label[this->intervalCount-1], debugPos[this->intervalCount-1], debugNeg[this->intervalCount-1]);
			sclib::stringOut("multiwaySplit.txt", debugOutput, this->pTweak, "\n");
			MFree_1D(debugPos);
			MFree_1D(debugNeg);
		}
	} else {
		printf("\nNo valid split found in multiway splitting (features are crap!), using findSingleBestSplit() instead!\n");
		minError = findSingleBestSplit(metaData, index, T, dim, totalPositives, totalNegatives);
	}

	return minError;
}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//  in the two-class-case, the constants SCIB_LABEL_POSITIVE/sclib::labelNegative should be used to indicate the 
//  classes, and they normaly should evaluate to +1/-1, respectively
//  the weights in the two arrays together (corresponding to rows in the SV_Data objects) are meant to sum up to 1
//====================================================================================================================
int SC_Classifier_DecisionStump::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double *positiveWeights, double *negativeWeights) {
	int t, T = pPositive->Row + pNegative->Row, dim = pPositive->Col;
	double **metaData, **metaIndex; // label|weight|posWeightSum|negWeightSum for each row in both datasets
	float **dataIndex; //rows correspond to data; includes in it's single column a pointer to a row of the original dataset
	double totalPositives = 0.0, totalNegatives = 0.0;

  //destroy previously trained classifier (and corresponding data), if any
	MFree_1D(this->threshold);
	MFree_1D(this->label);  
	MFree_1D(this->confidence);
	MFree_1D(this->margin);
	this->isTrained = false;
	this->classCount = 0;
	this->trainingError = -1.0;
  
	//form one list of labels, weights and patterns
	MArray_2D(metaData, T, 4, double, "SC_Classifier_DecisionStump.trainTwoClass: metaData");
	MArray_1D(metaIndex, T, double*, "SC_Classifier_DecisionStump.trainTwoClass: metaIndex");
	MArray_1D(dataIndex, T, float*, "SC_Classifier_DecisionStump.trainTwoClass: dataIndex");
	for (t = 0; t < T; t++) {
		metaIndex[t] = metaData[t]; //needed for later fast sorting
		if (t < pPositive->Row) {
			metaData[t][0] = (double)(sclib::labelPositive);
			metaData[t][1] = positiveWeights[t];
			totalPositives += positiveWeights[t];
			dataIndex[t] = pPositive->Mat[t];
		} else {
			metaData[t][0] = (double)(sclib::labelNegative);
			metaData[t][1] = negativeWeights[t-pPositive->Row];
			totalNegatives += negativeWeights[t-pPositive->Row];
			dataIndex[t] = pNegative->Mat[t-pPositive->Row];
		}//ignore to fill the positive and negative sum fields (idx 2 and 3) here because later sorting would require to do it again
	}

	//call suitable learning method...
	if (this->pTweak->classifierDecisionStump.splitMultiway == true) {
		this->trainingError = findMultiwaySplit(metaIndex, dataIndex, T, dim, totalPositives, totalNegatives);
	} else {
		if (this->optimizeF1pos == true) {
			this->trainingError = findSingleBestSplitF1pos(metaIndex, dataIndex, T, dim, totalPositives, totalNegatives); // a special case suitable to learn the best lambda in BIC speaker segmentation... shame on me that it is not more general here
		} else {
			this->trainingError = findSingleBestSplit(metaIndex, dataIndex, T, dim, totalPositives, totalNegatives);
		}
	}
	
	//...and eat the result
	if (this->trainingError < 0.0) {
		MFree_1D(this->threshold);
		MFree_1D(this->label);  
		MFree_1D(this->confidence);
		MFree_1D(this->margin);
		this->isTrained = false;
		this->classCount = 0;
	} else {
		this->isTrained = true;
		this->classCount = 2;
	}

	//cleanup & finish
	MFree_2D(metaData);
	MFree_1D(metaIndex);
	MFree_1D(dataIndex);

	return (this->trainingError >= 0.0) ? SVLIB_Ok : SVLIB_Fail;
}

//====================================================================================================================
//	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
//  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
//
//  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
//====================================================================================================================
int* SC_Classifier_DecisionStump::classify(SV_Data *pData, SV_Data* &pProbabilities) {
	SV_Data *pScaledData;
	int *classes = NULL, bestIdx = -1, i, i1;
	double thresh0 = -1.0 * std::numeric_limits<double>::max(), threshN = std::numeric_limits<double>::max();
	double thresh0Save, threshNsave, confidenceDueToCloseness, middle;

  if (this->isTrained == true) {
    MFree_0D(pProbabilities); 
		pProbabilities = new SV_Data(pData->Row, this->classCount);
    MArray_1D(classes, pData->Row, int, "SC_Classifier_DecisionStump.classify: classes");
		
		thresh0Save = this->threshold[0]; //replace the (informative) min/max values of the training set with the ranges of the datatype to span the whole range of possible values with the intervals
		this->threshold[0] = thresh0;
		threshNsave = this->threshold[this->intervalCount];
		this->threshold[this->intervalCount] = threshN;

    pScaledData = scaleFeatures(pData, this->pScale, -1, true);

    for (long int y = 0; y < pData->Row; y++) { //classify each feature vector separately
			i1 = 0; //i-1 to speed up things a little
			for (i = 1; i <= this->intervalCount; i++) { //determine in which interval the given example lies (TODO: the last interval isn't open due to "<" instead of "<="...)
				if ((double)(pScaledData->Mat[y][this->feature]) >= this->threshold[i1] && (double)(pScaledData->Mat[y][this->feature]) < this->threshold[i]) { //if the feature-value falls into this interval, its confidence and label is returned
					classes[y] = this->label[i1];
					for (int p = 0; p < this->classCount; p++) { //go trough all possible classes and assign confidence values accordingly
						if (idx2label(p) == classes[y]) {
							//confidences can be computed by the means of cleanness of the segment the example falls into (this is advocated by Schapire&Singer in "Improved Boosting Algorithms Using Confidence-Rated Predictions");
							//but another indicator for the confidence of a decision is the closeness of the example to the threshold itself; 
							//we multiply both measures together here in order to make a logical 'AND' between them: a segment should be clean and the example far away from it's borders to give a high confidence

							//first, we compute this new closeness-to-border-measure (smaller with decreasing distance to the threshold)
							if (pScaledData->Mat[y][this->feature]-this->threshold[i1] < this->threshold[i]-pScaledData->Mat[y][this->feature]) { //the example is closer to the lower threshold
								if (i == 1) {
									confidenceDueToCloseness = 1.0; //if this is the first (lowest) segment and the example is nearer to -inf than to the threshold, the confidence due to closeness can be assumed to be 1.0
								} else if (i == this->intervalCount) {
									middle = (threshNsave + this->threshold[i1]) / 2.0; //middle of current interval
									confidenceDueToCloseness = sclib::sigmoid(pScaledData->Mat[y][this->feature], this->threshold[i1], this->threshold[i1]+this->margin[i1-1], middle);
								} else {
									middle = (this->threshold[i] + this->threshold[i1]) / 2.0; //middle of current interval
									confidenceDueToCloseness = sclib::sigmoid(pScaledData->Mat[y][this->feature], this->threshold[i1], this->threshold[i1]+this->margin[i1-1], middle);
								}
							} else { //the example is closer to the upper border
								if (i == this->intervalCount) {
									confidenceDueToCloseness = 1.0; //same as above: if this is the last (highest) segment and the example is close to the upper border (+inf), confidence is 1.0
								} else if (i == 1) {
									middle = (this->threshold[i] + thresh0Save) / 2.0; //middle of current interval
									confidenceDueToCloseness = sclib::sigmoid(pScaledData->Mat[y][this->feature], this->threshold[i], this->threshold[i]-this->margin[i-1], middle);
								} else {
									middle = (this->threshold[i] + this->threshold[i1]) / 2.0; //middle of current interval
									confidenceDueToCloseness = sclib::sigmoid(pScaledData->Mat[y][this->feature], this->threshold[i], this->threshold[i]-this->margin[i-1], middle);
								}
							}

							//combine it with a logical 'AND'
							pProbabilities->Mat[y][p] = (float)(this->confidence[i1] * confidenceDueToCloseness); //this->confidence[i1] could also be named 'confidenceDueToPurity'
						} else {
							//because the confidence-scores are no real probabilities, we can't just give the converse probability in the other class' 
							//column: the predicted label has the associated confidence, and the non-predicted label has confidence 0 (otherwise, 
							//because its frequency-based and not merely closeness-to-threshold-based, it could be higher than the confidence of the 
							//actually predicted label)
							pProbabilities->Mat[y][p] = (float)(0.0); //in the multiclass-case, also all non-predicted labels have to have confidence==0.0
						}
					} //for p
					break; //fitting interval found, break
				} //if this is the correct interval
				i1++;
			} //for i
    } //for y

		if (pScaledData != pData) {
			MFree_0D(pScaledData);
		}

		this->threshold[0] = thresh0Save; //change back to informative min/max values
		this->threshold[this->intervalCount] = threshNsave;
  }

  return classes;
}

//====================================================================================================================
//	save a trained classifier to a file
//====================================================================================================================
int SC_Classifier_DecisionStump::saveClassifier(const char *fileName) {
  int i, res = SVLIB_Fail;

  if (strlen(fileName) > 0) {
    if (this->isTrained == true) {
			sclib::scalarOut(fileName, this->feature, NULL, false, "\n");
			sclib::scalarOut(fileName, this->intervalCount, NULL, false, "\n");
			for (i = 0; i < this->intervalCount+1; i++) {
				sclib::scalarOut(fileName, this->threshold[i], NULL, false, "\n");
			}
			for (i = 0; i < this->intervalCount; i++) {
				sclib::scalarOut(fileName, this->label[i], NULL, false, "\n");
			}
			for (i = 0; i < this->intervalCount; i++) {
				sclib::scalarOut(fileName, this->confidence[i], NULL, true, "\n");
			}
			for (i = 0; i < this->intervalCount-1; i++) {
				sclib::scalarOut(fileName, this->margin[i], NULL, true, "\n");
			}
			sclib::scalarOut(fileName, this->trainingError, NULL, false, "\n");
			res = SVLIB_Ok;
    }
  }

  return res;
}

//====================================================================================================================
//	load a trained classifier from a file
//====================================================================================================================
int SC_Classifier_DecisionStump::loadClassifier(const char *fileName) {
  int bytes, res = SVLIB_Ok;
  char buffer[sclib::bufferSize];
  FILE* inFile = NULL; 

  //destroy previously trained classifier, if any
	MFree_1D(this->threshold);
	MFree_1D(this->label);  
	this->isTrained = false;

  //load classifier
  if ((inFile = fopen(fileName, "r")) != NULL) {
    bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read feature-index
		if (bytes > 0) {
			this->feature = atoi(buffer);
		} else {
			res = SVLIB_Fail;
			return res;
		}

    bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read interval count
		if (bytes > 0) {
			this->intervalCount = atoi(buffer);
		} else {
			res = SVLIB_Fail;
			return res;
		}
		
		MArray_1D(this->threshold, this->intervalCount+1, double, "SC_Classifier_DecisionStump.loadClassifier: threshold");
    for (int i = 0; i < this->intervalCount+1; i++) {
			bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read thresholds
			if (bytes > 0) {
			this->threshold[i] = atof(buffer);
			} else {
				res = SVLIB_Fail;
				break;
			}
		}
		
		if (res != SVLIB_Fail) {
			MArray_1D(this->label, this->intervalCount, int, "SC_Classifier_DecisionStump.loadClassifier: label");
			for (int i = 0; i < this->intervalCount; i++) {
				bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read labels
				if (bytes > 0) {
					this->label[i] = atoi(buffer);
				} else {
					res = SVLIB_Fail;
					break;
				}
			}
		}

		if (res != SVLIB_Fail) {
			MArray_1D(this->confidence, this->intervalCount, double, "SC_Classifier_DecisionStump.loadClassifier: confidence");
			for (int i = 0; i < this->intervalCount; i++) {
				bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read confidence
				if (bytes > 0) {
					this->confidence[i] = atof(buffer);
				} else {
					res = SVLIB_Fail;
					break;
				}
			}
		}

		if (res != SVLIB_Fail) {
			MArray_1D(this->margin, this->intervalCount-1, double, "SC_Classifier_DecisionStump.loadClassifier: margin");
			for (int i = 0; i < this->intervalCount-1; i++) {
				bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read margin
				if (bytes > 0) {
					this->margin[i] = atof(buffer);
				} else {
					res = SVLIB_Fail;
					break;
				}
			}
		}

		if (res != SVLIB_Fail) {
			bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read training error
			if (bytes > 0) {
				this->trainingError = atof(buffer);
			} else {
				res = SVLIB_Fail;
			}
		}

		fclose(inFile);
  }

  if (res == SVLIB_Ok) {
    this->isTrained = true;
		this->classCount = getDistinctClassCount(this->label, this->intervalCount);
  } else {
		this->intervalCount = 0;
		MFree_1D(this->threshold);
		MFree_1D(this->label);  
		MFree_1D(this->confidence);
		MFree_1D(this->margin);
		this->isTrained = false;
		this->classCount = 0;
		this->trainingError = -1.0;
  }

  return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns -1 if something goes wrong
//====================================================================================================================
long int SC_Classifier_DecisionStump::label2idx(long int label) {
	long int res = -1; 

	switch (label) {
		case sclib::labelPositive: {
			res = 0;
			break;
		}
		case sclib::labelNegative: {
			res = 1;
			break;
		}
	}

	return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns sclib::noType if something goes wrong
//====================================================================================================================
long int SC_Classifier_DecisionStump::idx2label(long int idx) {
	long int res = sclib::noType;

	switch (idx) {
		case 0: {
			res = sclib::labelPositive;
			break;
		}
		case 1: {
			res = sclib::labelNegative;
			break;
		}
	}

	return res;
}

//====================================================================================================================
//	easy access to the class-separating threshold in case of a two-class problem
//====================================================================================================================
double SC_Classifier_DecisionStump::getTwoClassThreshold(void) {
	double thresh = -1.0;

	if (this->classCount == 2) {
		thresh = this->threshold[1];
	}

	return thresh;
}
