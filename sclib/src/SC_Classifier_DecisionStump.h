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

#ifndef __SC_Classifier_DecisionStump_H__
#define __SC_Classifier_DecisionStump_H__

#include "SC_Aux.h"
#include "SC_ClassifierWithWeights.h"
#include <SV_Error.h>

class SCLIB_API SC_Classifier_DecisionStump : public SC_ClassifierWithWeights {
  
  private :

		class SC_ThresholdList {
			public:
				SC_ThresholdList(double thresh = 0.0, double info = 0.0, int lowerLabel = sclib::labelNegative, double lowerPos = 0, double lowerNeg = 0, int upperLabel = sclib::labelPositive, double upperPos = 0, double upperNeg = 0, double margin = 0.0) {
					setValues(thresh, info, lowerLabel, lowerPos, lowerNeg, upperLabel, upperPos, upperNeg, margin);
					this->Next = NULL;
				}
				void setValues(double thresh, double info, int lowerLabel, double lowerPos, double lowerNeg, int upperLabel, double upperPos, double upperNeg, double margin) {
					this->threshold = thresh;
					this->informationGain = info;
					this->lowerPartLabel = lowerLabel;
					this->lowerPartClassCount[0] = lowerPos;
					this->lowerPartClassCount[1] = lowerNeg;
					this->upperPartLabel = upperLabel;
					this->upperPartClassCount[0] = upperPos;
					this->upperPartClassCount[1] = upperNeg;
					this->margin = margin;
				}
				double margin;
				double threshold;
				double informationGain;
				int lowerPartLabel;
				double lowerPartClassCount[2]; //idx==0 => positive, idx==1 => negative, as stated in idx2label()
				int upperPartLabel;
				double upperPartClassCount[2]; //idx==0 => positive, idx==1 => negative, as stated in idx2label()
				SC_ThresholdList *Next;
        int Valid(void) {return 1;}
		};

  protected :

		int intervalCount; //nr. of intervals created by the thresholds
		int feature; //index of the selected column (feature, attribute) of the training-dataset for classification
		double *threshold; //an array of thresholds in preparation for the multi-class case; meant to include upper and lower boundaries, so that in the two-class case there are 3 threshs: lower and upper bound as well as the learned threshold somewhere in between
		int *label; //the thresholds create #threshold-1 intervals; here are the class-labels associated with each interval
		double *confidence; //an array of confidence-scores (not normalized probabilities, but greater is more confident) for the prediction in each interval
		double *margin; //an array of margins, one for each inner threshold (i.e. not the min/max threshs): margin[0] corresponds with threshold[1] and so forth
		double trainingError; //if known, the training error can be saved here (automatically set by trainXClass())
		bool optimizeF1pos; //if true, F1 on positive class is optimzed instead of training error in the two-class single-split case

		//====================================================================================================================
		//  a special variant of quicksort that uses knowledge of the datasructures in trainTwoClass(): given is an index-
		//  array that holds per row a pointer to the rows of the actual data-matrix (not necessarily one single matrix, rows 
		//  can be distributed over many mazrixes), and a matrix associated with the data and index (rows correspond). now, 
		//  the index and the associated matrix get sorted while the original data is left as is. 
		//  rowLBound is lower bound, rowUBound is upper bound of the index-array and data-matrix' 1st dimension (rows, from 
		//  i=rowLBound;i<=rowUBound;i++).dataSortBy is the data-column which gives the reference for sorting
		//====================================================================================================================
		void indexQuickSort(float **index, double **associatedMatrix, int rowLBound, int rowUBound, int dataSortBy);

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
		void indexMergeSort(float **index, double **associatedMatrix, float **tmpIdx, double **tmpMat, int rowLBound, int rowUBound, int dataSortBy, bool sortByIndex = true);
		//void indexMergeCombine(float **index, double **associatedMatrix, int rowLBound, int rowUBound, int pivot, int dataSortBy, bool sortByIndex = true);

    //====================================================================================================================
    //	called by trainTwoClass() to find a single (2-way) split that minimizes the training error
		//  the method receives meta-data (with T rows of dimensionality dim, describing both the positive and negative 
		//  examples in their original matrixes) of the form [label|weight|posWeightSum|negWeightSum] and a corresponding 
		//  index-array with pointers (in it's single-column rows) to the original untouched two datasets.
		//  the result is directly inserted into the class-members label, threshold, feature and confidence, the return-value 
		//  gives the training error (-1.0 to indicat error)
    //====================================================================================================================
		double findSingleBestSplit(double **metaData, float **index, int T, int dim, double totalPositives, double totalNegatives);

		//====================================================================================================================
		//	as above, but with the F1 measure on the combined recall&precision of the positive class as the optimization
		//  criterion, in contrast to mere classification error. an additional constarint here: all samples above the
		//  threshold will be regarded as negative (opposed to: the predicted class of an interval is assumed to be the
		//  majority class within => this method is best used to compute the lambda penalty factor for BIC computations). 
		//  the best F1 is returned.
		//====================================================================================================================
		double findSingleBestSplitF1pos(double **metaData, float **index, int T, int dim, double totalPositives, double totalNegatives);

    //====================================================================================================================
    //	called by train*Class() to find the multiway-split that minimizes the entropy
		//  the method receives meta-data (with T rows of dimensionality dim, discribing both the positive and negative 
		//  examples in their original matrixes) of the form [label|weight|posWeightSum|negWeightSum] and a corresponding 
		//  index-array with pointers (int it's single-column rows) to the original untouched two datasets.
		//  the result is directly inserted into the class-members label, threshold, feature and confidence, the return-value 
		//  gives the training error (-1.0 to indicat error)
    //====================================================================================================================
		double findMultiwaySplit(double **metaData, float **index, int T, int dim, double totalPositives, double totalNegatives);

		//====================================================================================================================
		//  called by findMultiwaySplit() to recursively split the given data-set (see calling method for parameter info) as 
		//  long as the MDL criterion promotes the found splits; returns the training error on found splitting, or -1.0 on 
		//  error/failure to find split point; details for the found split points are returned in the linked list of 
		//  thresholds.
		//====================================================================================================================
		double splitInterval(double **metaData, float **index, int startRow, int stopRow, int column, double totalPositives, double totalNegatives, SC_Classifier_DecisionStump::SC_ThresholdList* &thresholdList);

		//====================================================================================================================
		//  takes a list of threshold-facts and returns the complete margin of the classifier composed of all the thresholds 
		//  in it as a weighted average (weighted by number of instances in the corresponding intervals) of the individual 
		//  margins (scaled so that margins on different featrues are comparable despite the feature where scaled previously 
		//  or not)
		//====================================================================================================================
		double getCompleteMargin(SC_Classifier_DecisionStump::SC_ThresholdList *pThreshList, int listCount, double minValue, double maxValue);

  public :

    SC_Classifier_DecisionStump(SC_TweakableParameters* pTweak, bool verbose = true, bool optimizeF1pos = false);
    virtual ~SC_Classifier_DecisionStump();

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //====================================================================================================================
    //virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative);
		//uses the base-class method!

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //  in the two-class-case, the constants sclib::labelPositive/sclib::labelNegative should be used to indicate the 
    //  classes, and they normaly should evaluate to +1/-1, respectively
		//  the weights in the two arrays together (corresponding to rows in the SV_Data objects) are meant to sum up to 1
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double *positiveWeights, double *negativeWeights);
		
    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes) {return SVLIB_Fail;};

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
		//  the weights (corresponding to rows in the SV_Data object) are meant to sum up to 1
    //====================================================================================================================
		virtual int trainMultiClass(SV_Data *pData, int *classes, double *weights) {return SVLIB_Fail;};

    //====================================================================================================================
    //	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
    //  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
		//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
    //
    //  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
    //====================================================================================================================
    virtual int* classify(SV_Data *pData, SV_Data* &pProbabilities);

    //====================================================================================================================
    //	save and load a trained classifier to/from a file
    //====================================================================================================================
    virtual int saveClassifier(const char *fileName);
    virtual int loadClassifier(const char *fileName);

    //====================================================================================================================
    //	to convert between given labels and indices into the probability-parameter of the classifiy()-method
    //====================================================================================================================
		virtual long int label2idx(long int label);
		virtual long int idx2label(long int idx);

    //====================================================================================================================
    //	getter/setter
    //====================================================================================================================
		int getFeature(void) {return this->feature;};
		void setFeature(int newFeature) {this->feature = (newFeature >= 0) ? newFeature : this->feature; return;}
		int getIntervalCount(void) {return this->intervalCount;};
		double getTwoClassThreshold(void); //easy access to the class-separating threshold in case of a two-class problem
		double getF1(void) {return (this->optimizeF1pos==true && this->isTrained==true) ? this->trainingError : -1.0;}
};

#endif
