/**************************************************************************/
/*    Responsibility:																											*/
/*      - implements adaboost (AdaC2) as described in "Cost-Sensitive     */
/*        Boosting for Classification of Imbalanced Data", Y. Sun, M.S.   */
/*        Kamel, A.K.C. Wong, Y. Wang, 2007 for the 2-class case and      */
/*        standard adaboost for the multiclass case.                      */
/*      - other ideas are from "AdaBoost.M1" as described in "Data        */
/*        Mining - Practical Machine Learning Tools and Techniques",      */
/*        Frank, Witten, 2nd Edition 2005, Elsevier, p. 322; "Robust      */
/*        Real-Time Face Detection", Viola, Jones, Int. J. Comp. Vision   */
/*        57(2), 137-154, 2004; "Real AdaBoost.MH" in Schapire, Singer,   */
/*        "Improved Boosting Algorithms Using Confidence-rated            */
/*        Predictions",                                                   */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.04.2007																								*/
/**************************************************************************/

#ifndef __SC_Classifier_AdaBoost_H__
#define __SC_Classifier_AdaBoost_H__

#include "SC_Aux.h"
#include "SC_Classifier.h"
#include "SC_ClassifierWithWeights.h"
#include <SV_Error.h>

class SCLIB_API SC_Classifier_AdaBoost : public SC_Classifier {
  
  private :

  protected :

		bool verbose;
		int weakClassifierType;
		int weakClassifierCount;
		double *alpha;
		SC_ClassifierWithWeights **pWeakClassifier;
		double trainingError; //if known, the training error can be saved here (automatically set by trainXClass())
		int *trainingExamplesCount; //if known, the nr of examples per class can be saved here (done by trainXClass()); there are 2*classCount+1 dimensions, always label and corresponding count in successive columns, plus the nr. of featrues in the last one

    //====================================================================================================================
    //	auxiliary method for freeing the weak classifiers
    //====================================================================================================================
		void freeWeakClassifiers(void);

  public :

    SC_Classifier_AdaBoost(SC_TweakableParameters* pTweak, int classifierType = sclib::ctDecisionStump, bool doScaling = false, bool verbose = true);
    virtual ~SC_Classifier_AdaBoost();

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative);
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double positiveMissclassificationCost, double negativeMissclassificationCost);
		
    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
		virtual int trainMultiClass(SV_Data *pData, int *classes);
		virtual int trainMultiClass(SV_Data *pData, int *classes, double *missclassificationCosts);

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
    //	to convert between given labels and indices into the probability-parameter of the classifiy()-method; at the 
		//  moment, probability-output is only possible in the case of two-class classification with probability-estimates on 
		//  the weak-learner side.
    //====================================================================================================================
		virtual long int label2idx(long int label);
		virtual long int idx2label(long int idx);

    //====================================================================================================================
    //	getter
    //====================================================================================================================
		int getWeakClassifieryCount(void) {return this->weakClassifierCount;}
		int getWeakClassifieryType(void) {return this->weakClassifierType;}
		SC_ClassifierWithWeights* getWeakClassifier(int idx) {return (idx >= 0 && idx < this->weakClassifierCount) ? this->pWeakClassifier[idx] : NULL;}
};

#endif
