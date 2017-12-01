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

#include <list>
#include <map>
#include "SC_Classifier_AdaBoost.h"
#include "SC_Classifier_DecisionStump.h"
#include "SC_Aux.h"
#include "SC_ClassifierHandler.h"
#include <SV_DataIO.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Classifier_AdaBoost::SC_Classifier_AdaBoost(SC_TweakableParameters* pTweak, int classifierType, bool doScaling, bool verbose) : SC_Classifier(pTweak, doScaling) {
	this->classifierType = sclib::ctAdaBoost;
	this->weakClassifierType = classifierType;
	this->verbose = verbose;
	this->weakClassifierCount = 0;
	this->alpha = NULL;
	this->pWeakClassifier = NULL;
	this->trainingError = -1.0; //some invalid value
	this->trainingExamplesCount = NULL;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Classifier_AdaBoost::~SC_Classifier_AdaBoost() {
	freeWeakClassifiers();
	MFree_1D(this->alpha);
	MFree_1D(this->trainingExamplesCount);
}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//====================================================================================================================
int SC_Classifier_AdaBoost::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative) {
	return trainTwoClass(pPositive, pNegative, 1.0, 1.0);
}

//====================================================================================================================
//	AdaC2 from Yanmin Sun (2007) et al.
//====================================================================================================================
int SC_Classifier_AdaBoost::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double positiveMissclassificationCost, double negativeMissclassificationCost) {
	int res = SVLIB_Ok, *labels = NULL, classifierCount = 0, t, T = pPositive->Row + pNegative->Row;
	double *posWeights, *negWeights, error, weightSum, *confidence = NULL;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	SC_ClassifierWithWeights *pFirst = NULL, *pHook = NULL;
	SC_ClassifierHandler handler(this->pTweak);
	SV_Data *pScaledPositive = NULL, *pScaledNegative = NULL;
	std::list<double> errors;
	std::list<double> alphas;
	double alpha, correct, wrong, initialWeightSum, Cp, Cn; //different costs of errors for positive and negative examples

	if (positiveMissclassificationCost > negativeMissclassificationCost) {
		Cp = 1.0; //fixed to 1.0
		Cn = negativeMissclassificationCost / positiveMissclassificationCost; //[0..1]
	} else {
		Cn = 1.0; //fixed to 1.0
		Cp = positiveMissclassificationCost / negativeMissclassificationCost; //[0..1]
	}

	//free previously trained classifier
	MFree_1D(this->alpha);
	freeWeakClassifiers();
	this->isTrained = false;
	this->classCount = 0;
	MFree_1D(this->trainingExamplesCount);

	//test the assumption
	if (handler.canHandleWeights(this->weakClassifierType) != true) {
		return SVLIB_Fail;
	}

	//find scaling parameters
	MFree_0D(this->pScale);
  if (this->doScaling == true) { 
    this->pScale = findScalingParameters(pPositive, pNegative);
	}
	pScaledPositive = scaleFeatures(pPositive, this->pScale, -1, true, true); //save some memory by just linking also in case of scaling; TODO
	pScaledNegative = scaleFeatures(pNegative, this->pScale, -1, true, true);

	//initialze the weights
	initialWeightSum = Cp*(double)(pScaledPositive->Row)/(double)(T) + Cn*(double)(pScaledNegative->Row)/(double)(T);
	posWeights = pFunc->initVector(pScaledPositive->Row, (double)((Cp * 1.0/(double)(T)) / initialWeightSum));
	negWeights = pFunc->initVector(pScaledNegative->Row, (double)((Cn * 1.0/(double)(T)) / initialWeightSum));

	//add weak classifiers till the error raises or falls too much
	while (classifierCount < (int)(this->pTweak->classifierAdaBoost.maxWeakClassifiers)) {

		//if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
		//	sclib::vectorOut("posWeights.txt", posWeights, pScaledPositive->Row, true, this->pTweak);
		//	sclib::vectorOut("negWeights.txt", negWeights, pScaledNegative->Row, true, this->pTweak);
		//}

		if (pFirst == NULL) {
			pFirst = (SC_ClassifierWithWeights*)(handler.buildClassifier(this->weakClassifierType, this->pTweak, false, this->verbose)); //scaling is handled here, not therein!
			pHook = pFirst;
		} else {
			pHook->Next = (SC_ClassifierWithWeights*)(handler.buildClassifier(this->weakClassifierType, this->pTweak, false, this->verbose)); //scaling is handled here, not therein!
			pHook = (SC_ClassifierWithWeights*)(pHook->Next);
		}
		
    res = pHook->trainTwoClass(pScaledPositive, pScaledNegative, posWeights, negWeights);
		if (res == SVLIB_Fail) {
			REPORT_ERROR(SVLIB_Fail, "Error during build of weak classifier");
			break;
		} else {
			error = pHook->getTrainingError(pScaledPositive, pScaledNegative, labels, confidence, posWeights, negWeights);
			if (this->verbose == true) {
				printf(" (error: %f) ", error);
			}
		}
		errors.push_back(error);

		if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClassifierTraining) == true) {
			sclib::scalarOut("error.txt", error, this->pTweak);
			if (this->weakClassifierType == sclib::ctDecisionStump) {
				sclib::scalarOut("feature.txt", ((SC_Classifier_DecisionStump*)pHook)->getFeature(), this->pTweak);
			}
		}

		//"choose" alpha_t so that the costs are regarded and the new set is maximally difficult (i.e. after reweighting, the sum of all weights of previously wrong classified examples is 0.5)
		correct = 0.0;
		wrong = 0.0;
		for (t = 0; t < pScaledPositive->Row; t++) {
			if (labels[t] == sclib::labelPositive) {
				correct += Cp * posWeights[t];
			} else {
				wrong += Cp * posWeights[t];
			}
		}
		for (t = 0; t < pScaledNegative->Row; t++) {
			if (labels[t+pScaledPositive->Row] == sclib::labelNegative) {
				correct += Cn * negWeights[t];
			} else {
				wrong += Cn * negWeights[t];
			}
		}
		alpha = 0.5 * log(correct / wrong); //this together with the reweighting formula assures that the next round gets maximaly difficult, i.e. wrong-weights=correct-weights=0.5
		alphas.push_back(alpha);

		//stopp adding further weak classifiers if error is too high or zero (cause then always the same weak learner would be added again and again due to no change in waits at zero error)
		//if (error >= 0.5 || error <= 0.0) { 
		if (error >= 0.5 || error <= 0.0 || wrong > correct) { 
			if (this->verbose == true) {
				printf(" stopped adding weak learner due to error=%f", error);
			}
			if (classifierCount > 0) { //discard a worse-than-chance classifier if there are already some others
				MFree_0D(pHook);
			} else { //special treatment if this would be the only classifier
				classifierCount++;
			}
			MFree_1D(confidence);
			MFree_1D(labels);
			break; //leaf the enclosing while-loop that adds more weak classifiers
		}

		//change weights of false negatives
		weightSum = 0.0;
		for (t = 0; t < pScaledPositive->Row; t++) { 
			if (labels[t] == sclib::labelPositive) { //this uses the knowledge that getTrainingError() for the two-class-case puts the labels of the positive examples first
				posWeights[t] *= Cp * exp(-1.0 * alpha); //exp(-1.0 * confidence[t]); //alpha is already folded into the hypothesis as in Schapire&Singer, 1999
			} else {
				posWeights[t] *= Cp * exp(+1.0 * alpha); //exp(+1.0 * confidence[t]);
			}
			weightSum += posWeights[t];
		}

		//change weights of false positives
		for (t = 0; t < pScaledNegative->Row; t++) { 
			if (labels[t+pScaledPositive->Row] == sclib::labelNegative) {
				negWeights[t] *= Cn * exp(-1.0 * alpha); //exp(-1.0 * confidence[t+pScaledPositive->Row]);
			} else {
				negWeights[t] *= Cn * exp(+1.0 * alpha); //exp(+1.0 * confidence[t+pScaledPositive->Row]);
			}
			weightSum += negWeights[t];
		}
		MFree_1D(confidence);
		MFree_1D(labels);

		//re-normalize the weights
		for (t = 0; t < T; t++) { 
			if (t < pScaledPositive->Row) {
				posWeights[t] /= weightSum;
			} else {
				negWeights[t-pPositive->Row] /= weightSum;
			}
		}

		//finish this round
		classifierCount++;

		if (this->verbose == true) {
			printf(".");
		}
	}

	//form final strong classifier
	if (res != SVLIB_Fail) {
		this->weakClassifierCount = classifierCount;
		MArray_1D(this->alpha, this->weakClassifierCount, double, "SC_Classifier_AdaBoost.trainTwoClass: alpha");
		MArray_1D(this->pWeakClassifier, this->weakClassifierCount, SC_ClassifierWithWeights*, "SC_Classifier_AdaBoost.trainTwoClass: pWeakClassifier");
		pHook = pFirst;
		for (t = 0; t < this->weakClassifierCount; t++) {
			error = errors.front();
			errors.pop_front();
			alpha = alphas.front();
			alphas.pop_front();

			if (error > 0.0)  { //standard requirement
				this->alpha[t] = alpha;
			} else { //special treatment for one single perfect classifier
				this->alpha[t] = 1.0;
			}

			this->pWeakClassifier[t] = pHook;
			pHook = (SC_ClassifierWithWeights*)(pHook->Next);
		}
		this->isTrained = true;
		this->classCount = 2;
		MArray_1D(this->trainingExamplesCount, 2*this->classCount+1, int, "SC_Classifier_AdaBoost.trainTwoClass: trainingExamplesCount");
		this->trainingExamplesCount[0] = sclib::labelPositive;
		this->trainingExamplesCount[1] = pScaledPositive->Row;
		this->trainingExamplesCount[2] = sclib::labelNegative;
		this->trainingExamplesCount[3] = pScaledNegative->Row;
		this->trainingExamplesCount[4] = pScaledPositive->Col;
		
		//fill the training error, therefore reinitialize the weights
		MFree_1D(posWeights);
		MFree_1D(negWeights);
		posWeights = pFunc->initVector(pScaledPositive->Row, (double)((Cp * 1.0/(double)(T)) / initialWeightSum));
		negWeights = pFunc->initVector(pScaledNegative->Row, (double)((Cn * 1.0/(double)(T)) / initialWeightSum));
		this->trainingError = getTrainingError(pScaledPositive, pScaledNegative, labels, confidence, posWeights, negWeights);
		MFree_1D(labels);
		MFree_1D(confidence);
	}

	//clean up
	MFree_1D(posWeights);
	MFree_1D(negWeights);
	MFree_0D(pFunc);
	if (pScaledPositive != pPositive) {
		MFree_0D(pScaledPositive);
	}
	if (pScaledNegative != pNegative) {
		MFree_0D(pScaledNegative);
	}

	return res;
}

//====================================================================================================================
//	train a classifier for distinguishing between several classes
//  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
//  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
//  respective row of pData.
//====================================================================================================================
int SC_Classifier_AdaBoost::trainMultiClass(SV_Data *pData, int *classes) {
	return trainMultiClass(pData, classes, NULL);;
}

//====================================================================================================================
//	AdaC2 from Yanmin Sun (2007) et al.
//  The missclassification costs are per sample (as the labels), not per class
//====================================================================================================================
int SC_Classifier_AdaBoost::trainMultiClass(SV_Data *pData, int *classes, double *missclassificationCosts) {
	int res = SVLIB_Ok, *labels, classifierCount = 0, t;
	double *weights, error, weightSum, *confidence;
	SC_ClassifierWithWeights *pFirst = NULL, *pHook = NULL;
	SC_ClassifierHandler handler(this->pTweak);
	SV_Data *pScaledData = NULL;
	SC_MatrixFunctions mFunc;
	std::list<double> errors;
	std::map<int, int> classesCount;
	std::map<int, int>::const_iterator i;
	double *costs, maxCost, correct, wrong, alpha, initialWeightSum = 0.0;
	std::list<double> alphas;

	//free previously trained classifier
	MFree_0D(this->pScale);
	MFree_1D(this->alpha);
	freeWeakClassifiers();
	this->isTrained = false;
	this->classCount = 0;
	MFree_1D(this->trainingExamplesCount);
	this->trainingError = -1.0;

	//test assumptions
	if (handler.canHandleWeights(this->weakClassifierType) != true ||
		  handler.canHandleMulticlassData(this->weakClassifierType) != true) {
		return SVLIB_Fail;
	}

	//find scaling parameters
  if (this->doScaling == true) { 
    this->pScale = findScalingParameters(pData);
  }
	pScaledData = scaleFeatures(pData, this->pScale, -1, true); //save some memory by just linking in case of no scaling

	//initialze the weights and the training examples statistics
	MArray_1D(weights, pScaledData->Row, double, "SC_Classifier_AdaBoost.trainMultiClass: weights");
	MArray_1D(costs, pScaledData->Row, double, "SC_Classifier_AdaBoost.trainMultiClass: costs");
	for (t = 0; t < pScaledData->Row; t++) {
		classesCount[classes[t]] += 1; //count number of instances in each class
	}
	this->classCount = (int)(classesCount.size()); //number of distinct class labels
	MArray_1D(this->trainingExamplesCount, this->classCount+1, int, "SC_Classifier_AdaBoost.trainMultiClass: trainingExamplesCount");
	t = 0;
	for (i = classesCount.begin(); i != classesCount.end(); ++i) {
		this->trainingExamplesCount[t++] = i->first;
		this->trainingExamplesCount[t++] = i->second;
	}
	maxCost = (missclassificationCosts != NULL) ? mFunc.max(missclassificationCosts, pScaledData->Row) : 0.0;
	for (t = 0; t < pScaledData->Row; t++) {
		costs[t] = (maxCost > 0.0) ? missclassificationCosts[t]/maxCost : 1.0; //normalize costs so that the maximum cost is 1.0
		weights[t] = costs[t] * 1.0/(double)(pScaledData->Row); //incorporate costs into first weights
		initialWeightSum += weights[t]; //used also later to get complete training error
	}
	for (t = 0; t < pScaledData->Row; t++) { //renormalize the weights to sum up to 1.0
		weights[t] /= initialWeightSum;
	}

	//add weak classifiers till the error raises or falls too much
	while (true) {
		if (pFirst == NULL) {
			pFirst = (SC_ClassifierWithWeights*)(handler.buildClassifier(this->weakClassifierType, this->pTweak, false, this->verbose)); //scaling is handled here, not therein!
			pHook = pFirst;
		} else {
			pHook->Next = (SC_ClassifierWithWeights*)(handler.buildClassifier(this->weakClassifierType, this->pTweak, false, this->verbose)); //scaling is handled here, not therein!
			pHook = (SC_ClassifierWithWeights*)(pHook->Next);
		}
		
    res = pHook->trainMultiClass(pScaledData, classes, weights);
		if (res == SVLIB_Fail) {
			break;
		} else {
			error = pHook->getTrainingError(pScaledData, classes, labels, confidence, weights);
			MFree_1D(confidence);
		}
		errors.push_back(error);

		//"choose" alpha_t so that the costs are regarded and the new set is maximally difficult (i.e. after reweighting, the sum of all weights of previously wrong classified examples is 0.5)
		correct = 0.0;
		wrong = 0.0;
		for (t = 0; t < pScaledData->Row; t++) {
			if (labels[t] == classes[t]) {
				correct += costs[t] * weights[t];
			} else {
				wrong += costs[t] * weights[t];
			}
		}
		alpha = 0.5 * log(correct / wrong); //0.5 has nothing to do with the class count - i prooved it!
		alphas.push_back(alpha);

		//stopp adding further weak classifiers if error is too low or too high
		if (error <= 0.0 || error >= 0.5 || wrong > correct) { 
			if (classifierCount > 0) { //discard a perfect or worse-than-chance classifier if there are already some others
				MFree_0D(pHook);
			} else { //special treatment for first (perfect) classifier
				classifierCount++;
			}
			MFree_0D(labels);
			break;
		}

		//reweight the dataset: increase weights of falsely classified examples, decrease weights of correct ones, regard costs
		weightSum = 0.0;
		for (t = 0; t < pData->Row; t++) {
			if (labels[t] == classes[t]) {
        weights[t] *= costs[t] * exp(-1.0 * alpha);
			} else {
				weights[t] *= costs[t] * exp(+1.0 * alpha);
			}
			weightSum += weights[t];
		}
		MFree_1D(labels);

		//re-normalize the weights
		for (t = 0; t < pData->Row; t++) { 
			weights[t] /= weightSum;
		}

		classifierCount++;
	}

	//form final strong classifier
	if (res != SVLIB_Fail) {
		this->weakClassifierCount = classifierCount;
		MArray_1D(this->alpha, this->weakClassifierCount, double, "SC_Classifier_AdaBoost.trainTwoClass: alpha");
		MArray_1D(this->pWeakClassifier, this->weakClassifierCount, SC_ClassifierWithWeights*, "SC_Classifier_AdaBoost.trainTwoClass: pWeakClassifier");
		pHook = pFirst;
		for (t = 0; t < this->weakClassifierCount; t++) {
			error = errors.front();
			errors.pop_front();
			alpha = alphas.front();
			alphas.pop_front();

			if (error > 0.0)  { //standard requirement
				this->alpha[t] = alpha; //-1.0 * log(error / (1.0 - error)); // == log(1.0 / (error / (1.0 - error))) as viola&jones state it
			} else { //special treatment for one single perfect classifier
				this->alpha[t] = 1.0;
			}
			this->pWeakClassifier[t] = pHook;
			pHook = (SC_ClassifierWithWeights*)(pHook->Next);
		}
		this->isTrained = true;
		this->classCount = classCount;
		this->trainingExamplesCount[this->classCount] = pScaledData->Col;

		//fill training error variable, therefore reinitialize the weights
		for (t = 0; t < pScaledData->Row; t++) {
			weights[t] = costs[t] * (1.0/(double)(pScaledData->Row)) / initialWeightSum;
		}
		this->trainingError = getTrainingError(pScaledData, classes, labels, confidence, weights);
		MFree_1D(labels);
		MFree_1D(classes);
	} else {
		this->classCount = 0;
		MFree_1D(this->trainingExamplesCount);
	}

	//clean up
	MFree_1D(costs);
	MFree_1D(weights);
  if (pScaledData != pData) {
    MFree_0D(pScaledData);
  }

	return SVLIB_Ok;
}

//====================================================================================================================
//	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
//  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
//
//  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
//====================================================================================================================
int* SC_Classifier_AdaBoost::classify(SV_Data *pData, SV_Data* &pProbabilities) {
	double maxWeight;
	int *labels = NULL, *label = NULL;
	int nrOfWeakLearners = sclib::min(this->weakClassifierCount, this->pTweak->classifierAdaBoost.maxWeakClassifiers);
	SV_Data *pScaledExample, *pProbability = NULL;
	std::map<int, double> classesCount;
	std::map<int, double> confidence;
	std::map<int, double>::const_iterator i;
	SC_ClassifierHandler handler(this->pTweak);
	bool giveConfidence = handler.canHandleProbabilities(this->weakClassifierType);
	
  if (this->doScaling == true && this->pScale == NULL) {
    REPORT_ERROR(SVLIB_BadData, "Can't do scaling if no scaling parameters where found");
  }

	if (this->isTrained == true) {
		MFree_0D(pProbabilities);
    MArray_1D(labels, pData->Row, int, "SC_Classifier_AdaBoost.classify: classes");
		if (giveConfidence == true) {
			pProbabilities = new SV_Data(pData->Row, this->classCount);
		}

    for (long int y = 0; y < pData->Row; y++) { //classify each feature vector separately
      pScaledExample = scaleFeatures(pData, this->pScale, ((pData->Row > 1) ? y : -1), true); //save time for copying in case of just one row

			for (int c = 0; c < nrOfWeakLearners; c++) { //get all votes from all weak classifiers
				label = this->pWeakClassifier[c]->classify(pScaledExample, pProbability);
				classesCount[label[0]] += this->alpha[c];
				if (giveConfidence == true) {
					confidence[label[0]] += this->alpha[c] * pProbability->Mat[0][label2idx(label[0])];
				}
				MFree_0D(pProbability);
				MFree_0D(label);
			}
			maxWeight = -1.0;
			for (i = classesCount.begin(); i != classesCount.end(); ++i) { //take the majority vote
				if (i->second > maxWeight) {
					maxWeight = i->second; //the value (weight)
					labels[y] = i->first; //the key (class-label)
				}
			}
			classesCount.clear();

			if (giveConfidence == true) {
				for (i = confidence.begin(); i != confidence.end(); ++i) {
					pProbabilities->Mat[y][label2idx(i->first)] = (float)(i->second);
				}
				confidence.clear();
			}

			if (pScaledExample != pData) {
				MFree_0D(pScaledExample);
			}
    }
  }

  return labels;
}

//====================================================================================================================
//	save a trained classifier to a file
//  write one summaryfile conatining the filenames of the single models
//====================================================================================================================
int SC_Classifier_AdaBoost::saveClassifier(const char *fileName) {
  int i, res = SVLIB_Ok;
  char *scaleFileName = NULL, temp[sclib::bufferSize];
  SV_DataIO io;

  if (strlen(fileName) > 0) {
    if (this->isTrained == true) {
			sclib::scalarOut(fileName, (this->doScaling == true) ? 1 : 0, NULL, false, "\n");
			sclib::scalarOut(fileName, (this->verbose == true) ? 1 : 0, NULL, false, "\n");
			sclib::scalarOut(fileName, this->weakClassifierType, NULL, false, "\n");
			sclib::scalarOut(fileName, this->weakClassifierCount, NULL, false, "\n");
			for (i = 0; i < this->weakClassifierCount; i++) {
				sclib::scalarOut(fileName, this->alpha[i], NULL, false, "\n");
				
				sprintf(temp, "%s_wc%d", fileName, i+1); //write the weak classifiers
				res = this->pWeakClassifier[i]->saveClassifier(temp);
				if (res == SVLIB_Fail) {
					break;
				}
			}
			sclib::scalarOut(fileName, this->trainingError, NULL, false, "\n"); //write trainign error
			for (i = 0; i < 2*this->classCount+1; i++) { //write labels and correspinding nr. of training examples as well as nr. of training-columns
				sclib::scalarOut(fileName, this->trainingExamplesCount[i], NULL, false, "\n");
			}

      if (res != SVLIB_Fail && this->doScaling == true) { //save also the scaling parameters of this training set to use it with the test-data to come
        scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
        io.OpenFile(scaleFileName, WRITE_REC);
		    res = io.PutDataRec(*this->pScale);
        res = (res == 0) ? SVLIB_Fail : SVLIB_Ok; //translate between different traditions to report errrors or success...
        io.CloseFile();
        MFree_1D(scaleFileName);
      }
    }
  }

  return res;
}

//====================================================================================================================
//	load a trained classifier from a file
//====================================================================================================================
int SC_Classifier_AdaBoost::loadClassifier(const char *fileName) {
  int bytes, i, res = SVLIB_Ok;
  char *scaleFileName = NULL, buffer[sclib::bufferSize];
  SV_DataIO IO;
	FILE* inFile = NULL; 
	SC_ClassifierHandler *pHandler = new SC_ClassifierHandler(this->pTweak);

  //destroy previously trained classifier, if any
	MFree_1D(this->alpha);
	freeWeakClassifiers();
	this->isTrained = false;
	this->classCount = 0;
	this->trainingError = -1.0;
	MFree_1D(this->trainingExamplesCount);

  //load classifier
  if ((inFile = fopen(fileName, "r")) != NULL) { 

    bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read doScaling
		if (bytes > 0) {
			i = atoi(buffer);
			this->doScaling = (i != 0) ? true : false;
		} else {
			res = SVLIB_Fail;
		}

		if (res != SVLIB_Fail) {
			bytes = sclib::readline(inFile, buffer, sclib::bufferSize); //read verbose
			if (bytes > 0) {
				i = atoi(buffer);
				this->verbose = (i != 0) ? true : false;
			} else {
				res = SVLIB_Fail;
			}
		}    

		if (res != SVLIB_Fail) { //read weakClassifierType
			bytes = sclib::readline(inFile, buffer, sclib::bufferSize); 
			if (bytes > 0) {
				this->weakClassifierType = atoi(buffer);
				if (pHandler->canHandleWeights(this->weakClassifierType) != true) {
					res = SVLIB_Fail;
				}
			} else {
				res = SVLIB_Fail;
			}
		}

		if (res != SVLIB_Fail) { //read weakClassifierCount
			bytes = sclib::readline(inFile, buffer, sclib::bufferSize); 
			if (bytes > 0) {
				this->weakClassifierCount = atoi(buffer);
				if (this->weakClassifierCount <= 0) {
					res = SVLIB_Fail;
				}
			} else {
				res = SVLIB_Fail;
			}
		}

		if (res != SVLIB_Fail) { //read alphas and weak classifiers
			MArray_1D(this->alpha, this->weakClassifierCount, double, "SC_Classifier_AdaBoost.loadClassifier: alpha");
			MArray_1D(this->pWeakClassifier, this->weakClassifierCount, SC_ClassifierWithWeights*, "SC_Classifier_AdaBoost.loadClassifier: pWeakClassifier");
			for (i = 0; i < this->weakClassifierCount; i++) {
				bytes = sclib::readline(inFile, buffer, sclib::bufferSize); 
				if (bytes > 0) {
					this->alpha[i] = atof(buffer);
				} else {
					res = SVLIB_Fail;
					break;
				}
				this->pWeakClassifier[i] = (SC_ClassifierWithWeights*)(pHandler->buildClassifier(this->weakClassifierType, this->pTweak, false, this->verbose)); //create weak classifier
				if (this->pWeakClassifier[i] == NULL) {
					res = SVLIB_Fail;
					break;
				} else {
					sprintf(buffer, "%s_wc%d", fileName, i+1);
					res = this->pWeakClassifier[i]->loadClassifier(buffer); //load weak classifier
					if (res == SVLIB_Fail) {
						break;
					}
				}
			}
			this->classCount = (this->weakClassifierCount > 0) ? this->pWeakClassifier[0]->getClassCount() : 0;
		}

		if (res != SVLIB_Fail) { //read training error (may be missing in file due to late introduction of parameter)
			bytes = sclib::readline(inFile, buffer, sclib::bufferSize); 
			if (bytes > 0) {
				this->trainingError = atof(buffer);
				for (i = 0; i < 2*this->classCount+1; i++) { //write labels and correspinding nr. of training examples as well as nr. of training-columns
					bytes = sclib::readline(inFile, buffer, sclib::bufferSize);
					if (bytes > 0) {
						if (i == 0) {
							MArray_1D(this->trainingExamplesCount, 2*this->classCount+1, int, "SC_Classifier_AdaBoost.loadClassifier: trainingExamplesCount");
						}
						this->trainingExamplesCount[i] = atoi(buffer);
					} else {
						MFree_1D(this->trainingExamplesCount);
						break;
					}
				}
			}
		}

		fclose(inFile);
	}

  if (res == SVLIB_Ok) {
    this->isTrained = true;
  } else {
		this->weakClassifierCount = 0;
		MFree_1D(this->alpha);
		MFree_1D(this->pWeakClassifier);
		this->isTrained = false;
		this->classCount = 0;
  }

  //also try to load (if exists) the scaling parameters:
  MFree_0D(this->pScale);
  if (res == SVLIB_Ok && this->doScaling == true) {
    scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
    if (sclib::fileExists(scaleFileName) == true) {
      IO.OpenFile(scaleFileName, READ_REC);
      this->pScale = IO.GetAllRec();
      IO.CloseFile();
    }
    MFree_1D(scaleFileName);
  }

  MFree_0D(pHandler);

  return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns -1 if something goes wrong; at the moment, probability-output is only possible in the case of two-class 
//  classification with probability-estimates on the weak-learner side.
//====================================================================================================================
long int SC_Classifier_AdaBoost::label2idx(long int label) {
	long int res = -1;

	if (this->isTrained == true && this->classCount == 2) {
		if (label == sclib::labelPositive) {
			res = 0;
		} else if (label == sclib::labelNegative) {
			res = 1;
		}
	}

	return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns sclib::noType if something goes wrong; at the moment, probability-output is only possible in the case of 
//  two-class classification with probability-estimates on the weak-learner side.
//====================================================================================================================
long int SC_Classifier_AdaBoost::idx2label(long int idx) {
	long int res = sclib::noType;
	
	if (this->isTrained == true && this->classCount == 2) {
		if (idx == 0) {
			res = sclib::labelPositive;
		} else if (idx == 1) {
			res = sclib::labelNegative;
		}
	}

	return res;
}

//====================================================================================================================
//	auxiliary method for freeing the weak classifiers
//====================================================================================================================
void SC_Classifier_AdaBoost::freeWeakClassifiers(void) {
	for (int i = 0; i < this->weakClassifierCount; i++) {
		MFree_0D(this->pWeakClassifier[i]);
	}
	MFree_1D(this->pWeakClassifier);
	this->weakClassifierCount = 0;

	return;
}
