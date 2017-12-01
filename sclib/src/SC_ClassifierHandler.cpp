/**************************************************************************/
/*    Responsibility:																											*/
/*		  - provides the possbility to build classifiers of all implemented */
/*        types                                                           */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 23.04.2007																								*/
/**************************************************************************/

#include "SC_ClassifierHandler.h"
#include "SC_Classifier_AdaBoost.h"
#include "SC_Classifier_DecisionStump.h"
#include "SC_Classifier_ML.h"
#include "SC_Classifier_SVM.h"
#include "SC_ClassifierTree.h"
#include "SC_ClassifierWithWeights.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_ClassifierHandler::SC_ClassifierHandler(SC_TweakableParameters *pTweak) {
  this->pTweak = pTweak;
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_ClassifierHandler::~SC_ClassifierHandler() {

}

//====================================================================================================================
//	build a classifier of the desired type
//====================================================================================================================
SC_Classifier* SC_ClassifierHandler::buildClassifier(int classifierType, SC_TweakableParameters *pTweak, bool doScaling, bool verbose) {
	SC_Classifier *pClassifier = NULL;

	switch (classifierType) {
		case sclib::ctAdaBoost: {
			pClassifier = new SC_Classifier_AdaBoost(pTweak, pTweak->classifierAdaBoost.weakClassifierType, doScaling);
			break;
		}
		case sclib::ctDecisionStump: {
			pClassifier = new SC_Classifier_DecisionStump(pTweak, verbose); //no caling here because it relies on just one feature and scaling is only needed to give individual features equal weight
			break;
		}
		case sclib::ctML: {
			pClassifier = new SC_Classifier_ML(pTweak, doScaling);
			break;
		}
		case sclib::ctSVM: {
			pClassifier = new SC_Classifier_SVM(pTweak, pTweak->classifierSvm.doCV, doScaling, verbose);
			break;
		}
		default: {
			break;
		}
	}

	return pClassifier;
}

//====================================================================================================================
//	returns true if the classifier can make use of a set of *weighted* examples for training
//====================================================================================================================
bool SC_ClassifierHandler::canHandleWeights(SC_Classifier *pClassifier) {
	bool res = false;

	if (typeid(pClassifier) == typeid(SC_Classifier_DecisionStump*) ||
			typeid(pClassifier) == typeid(SC_ClassifierWithWeights*) ||
			typeid(pClassifier) == typeid(SC_Classifier_SVM*)) {
		res = true;
	}

	return res;
}

//====================================================================================================================
//	returns true if the classifier can make use of a set of *weighted* examples for training
//====================================================================================================================
bool SC_ClassifierHandler::canHandleWeights(int classifierType) {
	bool res = false;

	if (classifierType == sclib::ctDecisionStump ||
		  classifierType == sclib::ctSVM) {
		res = true;
	}

	return res;
}

//====================================================================================================================
//	returns true if the classifier returns also probability-estimates on testing
//====================================================================================================================
bool SC_ClassifierHandler::canHandleProbabilities(SC_Classifier *pClassifier) {
	bool res = false;

	if (typeid(pClassifier) == typeid(SC_Classifier_SVM*) ||
		  typeid(pClassifier) == typeid(SC_Classifier_ML*) ||
			typeid(pClassifier) == typeid(SC_Classifier_DecisionStump*)) {
		res = true;
	} else if (typeid(pClassifier) == typeid(SC_ClassifierTree*)) {
		res = canHandleProbabilities(((SC_ClassifierTree*)pClassifier)->getClassifierType());
	}

	return res;
}

//====================================================================================================================
//	returns true if the classifier returns also probability-estimates on testing
//====================================================================================================================
bool SC_ClassifierHandler::canHandleProbabilities(int classifierType) {
	bool res = false;

	if (classifierType == sclib::ctSVM ||
		  classifierType == sclib::ctML ||
			classifierType == sclib::ctDecisionStump) {
		res = true;
	}

	return res;
}

//====================================================================================================================
//	returns true if the classifier is able to distinguish between more than 2 classes
//  (every classifier must provide a corresponding method to fullfill the interface, but this doesn't necessaryly mean
//  that it actually does something inside...)
//====================================================================================================================
bool SC_ClassifierHandler::canHandleMulticlassData(SC_Classifier *pClassifier) {
	bool res = false;

	if (typeid(pClassifier) == typeid(SC_Classifier_SVM*) ||
		  typeid(pClassifier) == typeid(SC_Classifier_ML*) ||
			typeid(pClassifier) == typeid(SC_Classifier_AdaBoost*)) {
		res = true;
	} else if (typeid(pClassifier) == typeid(SC_ClassifierTree*)) {
		res = canHandleMulticlassData(((SC_ClassifierTree*)pClassifier)->getClassifierType());
	}


	return res;
}

//====================================================================================================================
//	returns true if the classifier is able to distinguish between more than 2 classes
//  (every classifier must provide a corresponding method to fullfill the interface, but this doesn't necessaryly mean
//  that it actually does something inside...)
//====================================================================================================================
bool SC_ClassifierHandler::canHandleMulticlassData(int classifierType) {
	bool res = false;

	if (classifierType == sclib::ctML ||
		  classifierType == sclib::ctSVM ||
			classifierType == sclib::ctAdaBoost) {
		res = true;
	}

	return res;
}
