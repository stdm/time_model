/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for classifiers that are able to handle weighted     */
/*        instances                                                       */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.04.2007																								*/
/**************************************************************************/

#include "SC_ClassifierWithWeights.h"
#include "SC_MatrixFunctions.h"
#include "SC_Aux.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_ClassifierWithWeights::SC_ClassifierWithWeights(SC_TweakableParameters* pTweak, bool doScaling, bool verbose) : SC_Classifier(pTweak, doScaling) {
	this->verbose = verbose;
}

//====================================================================================================================
//	copy-constructor
//====================================================================================================================
SC_ClassifierWithWeights::SC_ClassifierWithWeights(const SC_ClassifierWithWeights& pParent) : SC_Classifier(pParent) {
	this->verbose = pParent.verbose;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_ClassifierWithWeights::~SC_ClassifierWithWeights() {

}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//====================================================================================================================
int SC_ClassifierWithWeights::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative) {
	int res, T = pPositive->Row + pNegative->Row;
	double *posWeights, *negWeights;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

	posWeights = pFunc->initVector(pPositive->Row, 1.0 / (double)(T)); //create uniform weights if none are given
	negWeights = pFunc->initVector(pNegative->Row, 1.0 / (double)(T));

	res = trainTwoClass(pPositive, pNegative, posWeights, negWeights);

	MFree_1D(posWeights);
	MFree_1D(negWeights);
	MFree_0D(pFunc);

	return res;
}

//====================================================================================================================
//	train a classifier for distinguishing between the classes, for which examples are given in the SV_Data object
//====================================================================================================================
int SC_ClassifierWithWeights::trainMultiClass(SV_Data *pData, int *classes) {
	int res;
	double *weights;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

	weights = pFunc->initVector(pData->Row, 1.0 / (double)(pData->Row)); //create uniform weights if none are given

	res = trainMultiClass(pData, classes, weights);

	MFree_1D(weights);
	MFree_0D(pFunc);

	return res;
}
