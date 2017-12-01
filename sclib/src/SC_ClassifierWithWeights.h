/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for classifiers that are able to handle weighted     */
/*        instances                                                       */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.04.2007																								*/
/**************************************************************************/

#ifndef __SC_ClassifierWithWeights_H__
#define __SC_ClassifierWithWeights_H__

#include "SC_Classifier.h"

class SCLIB_API SC_ClassifierWithWeights : public SC_Classifier {
  
  private :

  protected :

		bool verbose;

  public :

    SC_ClassifierWithWeights(SC_TweakableParameters* pTweak, bool doScaling, bool verbose);
		SC_ClassifierWithWeights(const SC_ClassifierWithWeights& pParent);
    virtual ~SC_ClassifierWithWeights();

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //  in the two-class-case, the constants SCIB_LABEL_POSITIVE/sclib::labelNegative should be used to indicate the 
    //  classes, and they normaly should evaluate to +1/-1, respectively
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative);

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //  in the two-class-case, the constants SCIB_LABEL_POSITIVE/sclib::labelNegative should be used to indicate the 
    //  classes, and they normaly should evaluate to +1/-1, respectively
		//  the weights in the two arrays together (corresponding to rows in the SV_Data objects) are meant to sum up to 1
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double *positiveWeights, double *negativeWeights) = 0;

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes);

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
		//  the weights (corresponding to rows in the SV_Data object) are meant to sum up to 1
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes, double *weights) = 0;

    //====================================================================================================================
    //	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
    //  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
		//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
    //
    //  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
    //====================================================================================================================
    virtual int* classify(SV_Data *pData, SV_Data *&pProbabilities) = 0;
    
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
};

#endif
