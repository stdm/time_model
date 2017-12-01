/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for classifiers                                      */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Classifier_H__
#define __SC_Classifier_H__

#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include <SV_Data.h>

//now in pTweak
//#define sclib::labelNegative (-1)
//#define sclib::labelPositive (+1)

class SCLIB_API SC_Classifier {
  
  private :

  protected :

    SC_TweakableParameters* pTweak; //for user-defined parameters
    SV_Data *pScale; //scaling parameters for the test data, conditioned on the training data
    bool isTrained; //is true only if the classifier has been trained or loaded
    bool doScaling; //controls whether the traing- and testdata should be scaled to [0..1] prior to be subject to the SVM algorithms
		int classifierType; //a label giving the type of classifier
		int classCount; //number of classes that can be distinguished by the current trained classifier

  public :

		SC_Classifier *Next; //to form a linked list of classifiers

    SC_Classifier(SC_TweakableParameters* pTweak, bool doScaling);
		SC_Classifier(const SC_Classifier& pParent);
    virtual ~SC_Classifier();

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //  in the two-class-case, the constants SCIB_LABEL_POSITIVE/sclib::labelNegative should be used to indicate the 
    //  classes, and they normaly should evaluate to +1/-1, respectively
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative) = 0;

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes) = 0;

    //====================================================================================================================
    //	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
    //  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
		//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
    //
    //  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
    //====================================================================================================================
    virtual int* classify(SV_Data *pData, SV_Data* &pProbabilities) = 0;
    
    //====================================================================================================================
    //	save and load a trained classifier to/from a file
    //====================================================================================================================
    virtual int saveClassifier(const char *fileName) = 0;
    virtual int loadClassifier(const char *fileName) = 0;

    //====================================================================================================================
    //	to convert between given labels and indices into the probability-parameter of the classifiy()-method
    //====================================================================================================================
		virtual long int label2idx(long int label) = 0;
		virtual long int idx2label(long int idx) = 0;

    //====================================================================================================================
    //	methods to find the training error for a just build classifier given the training data in the appropriate form;
		//  the actual classifications of each example are also returned in the "labels" parameter; weights for the examples
		//  can be provided, otherwise uniform weighting is assumed
    //====================================================================================================================
		double getTrainingError(SV_Data *pPositive, SV_Data *pNegative, int* &labels, double* &confidence, double *positiveWeights = NULL, double *negativeWeights = NULL);
		double getTrainingError(SV_Data *pData, int *classes, int* &labels, double* &confidence, double *weights = NULL);

    //====================================================================================================================
    //	given the set of labels, return number of disitinct classes to set classCount
    //====================================================================================================================
		int getDistinctClassCount(int *classes, long int length);
		int getDistinctClassCount(int *classes, long int length, int* &statistics);

    //====================================================================================================================
    //	getter
    //====================================================================================================================
		int getClassifierType(void) {return this->classifierType;};
		int getClassCount(void) {return this->classCount;};
		bool doesScaling(void) {return this->doScaling;};
		SV_Data *getScalingParameters(void) {return this->pScale;};

    //====================================================================================================================
    //	setter
    //====================================================================================================================
		void setDoScaling(bool doScaling) {this->doScaling = doScaling; return;}

    //====================================================================================================================
    //	find min and max values of each column of the feature vectors, so they can later be scaled to the range [0..1] 
    //  using this values
		//  TODO: could be made more robust by using 10th percentile instead of min and 50th percentile instead of max to 
		//        scale between 0 and 0.5, levelling off below 0 and above 1 (would need to retrain already built classifiers, 
		//        though)
    //====================================================================================================================
    SV_Data* findScalingParameters(SV_Data *pData);
		SV_Data* findScalingParameters(SV_Data *pPositive, SV_Data *pNegative);

    //====================================================================================================================
    //	scale the given features with the given scaling parameters, return a new result set and don't touch the original
    //  features; if the scaling parameters are NULL, return a copy of the original data
    //  if row!=-1, only the specified row get's copied
		//  exception: if forceNoCopy==true and row isn't specified, the dataset gets not copied but just the pointer to pData 
		//  is returned in case of no scaling parameters given; if scaling-parameters are given and all rows shall be treated 
		//  and anyhow no copy should be created, evenOverwrite can be set to true, and the the original data is modified 
		//  (scaled) and a pointer to it returned.
		//  if scaling actually applies, each example is multiplied by factor *after* scaling (this way, pos. and neg. ex-
		//  amples can be both scaled between 0..1 with the same parameters and the stretched to -1..1 afterwards)
    //====================================================================================================================
    SV_Data* scaleFeatures(SV_Data *pData, SV_Data *pScalingParameters, long int row = -1, bool forceNoCopy = false, bool evenOverwrite = false, double factor = 1.0);

		/*
    //====================================================================================================================
    //	some tests towards learning to sort things by using binary classification
    //====================================================================================================================
		virtual int trainToSort(SV_Data *pOrderedFeatures);
		virtual SV_Data* sort(SV_Data *pFeatures);
		*/
};

#endif
