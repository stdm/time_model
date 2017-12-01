/**************************************************************************/
/*    Responsibility:																											*/
/*		  - provides the possbility to build classifiers of all implemented */
/*        types                                                           */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 23.04.2007																								*/
/**************************************************************************/

#ifndef __SC_ClassifierHandler_H__
#define __SC_ClassifierHandler_H__

#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_Classifier.h"

class SCLIB_API SC_ClassifierHandler {
	private :

  protected :

    SC_TweakableParameters *pTweak;

  public :
		
    SC_ClassifierHandler(SC_TweakableParameters *pTweak);
    virtual ~SC_ClassifierHandler();

    //====================================================================================================================
    //	build a classifier of the desired type
    //====================================================================================================================
		SC_Classifier* buildClassifier(int classifierType, SC_TweakableParameters *pTweak, bool doScaling, bool verbose);

    //====================================================================================================================
    //	returns true if the classifier can make use of a set of *weighted* examples for training
    //====================================================================================================================
		bool canHandleWeights(SC_Classifier *pClassifier);
		bool canHandleWeights(int classifierType);

    //====================================================================================================================
    //	returns true if the classifier returns also probability-estimates on testing
    //====================================================================================================================
		bool canHandleProbabilities(SC_Classifier *pClassifier);
		bool canHandleProbabilities(int classifierType);

    //====================================================================================================================
    //	returns true if the classifier is able to distinguish between more than 2 classes
		//  (every classifier must provide a corresponding method to fullfill the interface, but this doesn't necessaryly mean
		//  that it actually does something inside...)
    //====================================================================================================================
		bool canHandleMulticlassData(SC_Classifier *pClassifier);
		bool canHandleMulticlassData(int classifierType);
};

#endif
