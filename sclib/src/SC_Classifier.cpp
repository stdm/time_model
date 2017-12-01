/**************************************************************************/
/*    Responsibility:																											*/
/*      - Base class for classifiers                                      */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.03.2006																								*/
/**************************************************************************/

#include <limits.h>
#include <map>
#include <vector>
#include "SC_Classifier.h"
#include "SC_Aux.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Classifier::SC_Classifier(SC_TweakableParameters* pTweak, bool doScaling) {
  this->classCount = 0;
	this->classifierType = sclib::ctUndefined;
	this->pTweak = pTweak;
  this->isTrained = false;
  this->doScaling = doScaling;
  this->pScale = NULL;
	this->Next = NULL;
}

//====================================================================================================================
//	copy-constructor
//====================================================================================================================
SC_Classifier::SC_Classifier(const SC_Classifier& pParent) {
	this->classCount = pParent.classCount;
	this->classifierType = pParent.classifierType;
	this->pTweak = pParent.pTweak;
	this->isTrained = pParent.isTrained;
	this->doScaling = pParent.doScaling;
	this->pScale = (pParent.pScale != NULL) ? new SV_Data((SV_Data&)(*(pParent.pScale)), false) : NULL;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Classifier::~SC_Classifier() {
  MFree_0D(this->pScale);
}

//====================================================================================================================
//	find min and max values of each column of the feature vectors, so they can later be scaled to the range [0..1] 
//  using this values
//  TODO: could be made more robust by using 10th percentile instead of min and 50th percentile instead of max to 
//        scale between 0 and 0.5, levelling off below 0 and above 1 (would need to retrain already built classifiers, 
//        though)
//====================================================================================================================
SV_Data* SC_Classifier::findScalingParameters(SV_Data *pData) {
  long int x, y;
  SV_Data *pScale = new SV_Data(2, pData->Col); //first row takes minima of each column of pData, second row takes maxima
 
  //initialize the min/max values
  for (x = 0; x < pScale->Col; x++) {
    pScale->Mat[0][x] = std::numeric_limits<float>::max();
    pScale->Mat[1][x] = std::numeric_limits<float>::min();
  }

  //compute minima/maxima
  for (y = 0; y < pData->Row; y++) {
    for (x = 0; x < pData->Col; x++) {
      if (pData->Mat[y][x] < pScale->Mat[0][x]) { //find minimum per column
        pScale->Mat[0][x] = pData->Mat[y][x];
      }
      if (pData->Mat[y][x] > pScale->Mat[1][x]) { //find maximum per column
        pScale->Mat[1][x] = pData->Mat[y][x];
      }
    }
  }

  return pScale;
}

//====================================================================================================================
//	find min and max values of each column of the feature vectors in the two datasets, so they can later be scaled to 
//  the range [0..1] using this values
//  this method makes it unnecessary to copy both datasets together to find scaling parameters and therefor saves
//  memory and computational time
//====================================================================================================================
SV_Data* SC_Classifier::findScalingParameters(SV_Data *pPositive, SV_Data *pNegative) {
  long int x, y;
  SV_Data *pScale = NULL;
	
	if (pPositive->Col = pNegative->Col) {
		pScale = new SV_Data(2, pPositive->Col); //first row takes minima of each column of pData, second row takes maxima
	 
		//initialize the min/max values
		for (x = 0; x < pScale->Col; x++) {
			pScale->Mat[0][x] = std::numeric_limits<float>::max();
			pScale->Mat[1][x] = std::numeric_limits<float>::min();
		}

		//compute minima/maxima over both datasets
		for (y = 0; y < pPositive->Row; y++) {
			for (x = 0; x < pPositive->Col; x++) {
				if (pPositive->Mat[y][x] < pScale->Mat[0][x]) { //find minimum per column
					pScale->Mat[0][x] = pPositive->Mat[y][x];
				}
				if (pPositive->Mat[y][x] > pScale->Mat[1][x]) { //find maximum per column
					pScale->Mat[1][x] = pPositive->Mat[y][x];
				}
			}
		}
		for (y = 0; y < pNegative->Row; y++) {
			for (x = 0; x < pNegative->Col; x++) {
				if (pNegative->Mat[y][x] < pScale->Mat[0][x]) { //find minimum per column
					pScale->Mat[0][x] = pNegative->Mat[y][x];
				}
				if (pNegative->Mat[y][x] > pScale->Mat[1][x]) { //find maximum per column
					pScale->Mat[1][x] = pNegative->Mat[y][x];
				}
			}
		}
	}

  return pScale;
}

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
SV_Data* SC_Classifier::scaleFeatures(SV_Data *pData, SV_Data *pScalingParameters, long int row, bool forceNoCopy, bool evenOverwrite, double factor) {
  SV_Data *pScaledData = pData;
  int start, end, count = 0;

  if (pScalingParameters != NULL) { //scaling parameters given?
		if (row >= 0 && row < pData->Row) {
			start = row;
			end = row+1;
		} else {
			start = 0;
			end = pData->Row;
		}

		if (forceNoCopy == false || evenOverwrite == false || start < 0 || end > pData->Row) { //only link and modify original data if both parameters are true and all rows shall be treated
			pScaledData = new SV_Data(end-start, pData->Col);
			pScaledData->Hdr.frameSize = pData->Hdr.frameSize; //copy important header info
			pScaledData->Hdr.frameStep = pData->Hdr.frameStep;
			pScaledData->Hdr.sampleRate = pData->Hdr.sampleRate;
			pScaledData->Hdr.ID = pData->Hdr.ID;
		} else {
			pScaledData = pData;
		}

    for (int y = start; y < end; y++) {
      for (int x = 0; x < pScaledData->Col; x++) {
        pScaledData->Mat[count][x] = (pData->Mat[y][x] - pScalingParameters->Mat[0][x]) / (pScalingParameters->Mat[1][x] - pScalingParameters->Mat[0][x]); //(x-min)/(max-min)
				pScaledData->Mat[count][x] *= (float)(factor);
      }
      count++;
    }
  } else { //no scaling parameters available
    if (row == -1) {
      if (forceNoCopy == true) {
        pScaledData = pData;
      } else {
        pScaledData = new SV_Data((SV_Data&)*pData);
      }
    } else {
      pScaledData = new SV_Data(1, pData->Col);
      for (int x = 0; x < pScaledData->Col; x++) {
        pScaledData->Mat[0][x] = pData->Mat[row][x];
      }
    }
  }

  return pScaledData;  
}

//====================================================================================================================
//	methods to find the training error for a just build classifier given the training data in the two-class form
//  the actual classifications of each example are also returned in the "labels" parameter; the confidence 
//  (probability) of the predicted class is also returned in the "confidence" parameter
//====================================================================================================================
double SC_Classifier::getTrainingError(SV_Data *pPositive, SV_Data *pNegative, int* &labels, double* &confidence, double *positiveWeights, double *negativeWeights) {
	double error = 0.0, increment = 1.0 / (double)(pPositive->Row+pNegative->Row); //assuming uniform weights...
	int *predictions, i;
	SV_Data *pProbabilities = NULL;

	if (this->isTrained == true) {
		MFree_1D(labels);
		MArray_1D(labels, pPositive->Row+pNegative->Row, int, "SC_Classifier.getTrainingError: labels");
		MFree_1D(confidence);
		MArray_1D(confidence, pPositive->Row+pNegative->Row, double, "SC_Classifier.getTrainingError: confidence");

		predictions = classify(pPositive, pProbabilities);
		for (i = 0; i < pPositive->Row; i++) { //attention: there are some derived classes that internally use the knowledge that this method puts the labels of the positive examples first in the array and then the negatives!
			labels[i] = predictions[i];
			confidence[i] = pProbabilities->Mat[i][label2idx(labels[i])];
			if (predictions[i] != sclib::labelPositive) {
				error += (positiveWeights == NULL) ? increment : positiveWeights[i];
			}
		}
		MFree_0D(pProbabilities);
		MFree_1D(predictions);

		predictions = classify(pNegative, pProbabilities);
		for (i = 0; i < pNegative->Row; i++) {
			labels[i+pPositive->Row] = predictions[i];
			confidence[i+pPositive->Row] = pProbabilities->Mat[i][label2idx(labels[i+pPositive->Row])];
			if (predictions[i] != sclib::labelNegative) {
				error += (negativeWeights == NULL) ? increment : negativeWeights[i];
			}
		}
		MFree_0D(pProbabilities);
		MFree_1D(predictions);
	}	

	return error;
}

//====================================================================================================================
//	methods to find the training error for a just build classifier given the training data in the multi-class form
//  the actual classifications of each example are also returned in the "labels" parameter; weights for the examples
//  can be provided, otherwise uniform weighting is assumed; the confidence (probability) of the predicted class is 
//  also returned in the "confidence" parameter
//====================================================================================================================
double SC_Classifier::getTrainingError(SV_Data *pData, int *classes, int* &labels, double* &confidence, double *weights) {
	double error = 0.0, increment = 1.0 / (double)(pData->Row); //assuming uniform weights...
	int i;
	SV_Data *pProbabilities = NULL;

	if (this->isTrained == true) {
		MFree_1D(labels);
		MFree_1D(confidence);
		MArray_1D(confidence, pData->Row, double, "SC_Classifier.getTrainingError: confidence");
		labels = classify(pData, pProbabilities);
		for (i = 0; i < pData->Row; i++) {
			confidence[i] = pProbabilities->Mat[i][label2idx(labels[i])];
			if (labels[i] != classes[i]) {
				error += (weights == NULL) ? increment : weights[i];
			}
		}
		MFree_0D(pProbabilities);
	}	

	return error;
}

//====================================================================================================================
//	given the set of labels, return number of disitinct classes to set classCount
//====================================================================================================================
int SC_Classifier::getDistinctClassCount(int *classes, long int length) {
	std::map<int, int> counter;
	
	for (int t = 0; t < length; t++) {
		counter[classes[t]]++;
	}

	return (int)(counter.size());
}

//====================================================================================================================
//	given the set of labels, return number of disitinct classes to set classCount; also, in statistics, an array of
//  integers is created and returned having 2*classCount+1 columns: always the label and corresponding nr. of examples
//  in successive columns, plus one last column for extra data to be filled from outside (e.g. nr. of features)
//====================================================================================================================
int SC_Classifier::getDistinctClassCount(int *classes, long int length, int* &statistics) {
	int t, T;
	std::map<int, int> counter;
	std::map<int, int>::const_iterator i;
	
	for (t = 0; t < length; t++) {
		counter[classes[t]]++;
	}
	T = (int)(counter.size());

	MFree_1D(statistics);
	MArray_1D(statistics, 2*T+1, int, "SC_Classifier.getDistinctClassCount: statistics");
	t = 0;
	for (i = counter.begin(); i != counter.end(); ++i) {
		statistics[t++] = i->first;
		statistics[t++] = i->second;
	}

	return T;
}

/*
//====================================================================================================================
//	some tests towards learning to sort things by using binary classification
//====================================================================================================================
int SC_Classifier::trainToSort(SV_Data *pOrderedFeatures) {
	int i, j, res = SVLIB_Fail;
	SV_Data *pPositive, *pNegative;

	pPositive = new SV_Data(pOrderedFeatures->Row-1, pOrderedFeatures->Col); //*2);
	pNegative = new SV_Data(pOrderedFeatures->Row-1, pOrderedFeatures->Col); //*2);

	for (i = 0; i < pOrderedFeatures->Row-1; i++) {
		for (j = 0; j < pOrderedFeatures->Col; j++) {
			pPositive->Mat[i][j] = pOrderedFeatures->Mat[i][j] - pOrderedFeatures->Mat[i+1][j];
			//pPositive->Mat[i][pOrderedFeatures->Col+j] = pOrderedFeatures->Mat[i+1][j];
			pNegative->Mat[i][j] = pOrderedFeatures->Mat[i+1][j] - pOrderedFeatures->Mat[i][j];
			//pNegative->Mat[i][pOrderedFeatures->Col+j] = pOrderedFeatures->Mat[i][j];
		}
	}

	res = trainTwoClass(pPositive, pNegative);

	MFree_0D(pPositive);
	MFree_0D(pNegative);

	return res;
}

SV_Data* SC_Classifier::sort(SV_Data *pFeatures) {
	SV_Data *pSorted = NULL, *pProbe, *pProbability = NULL;;
	double **mayFollow, rowSum, bestRowSum;
	int d, t, *result, bestIdx, lastIdx = -1;
	unsigned int i, j;
	vector<int> availableIdx;
	bool allZero;

	MArray_2D(mayFollow, pFeatures->Row, pFeatures->Row, double, "SC_Classifier.sort: mayFollow");
	pSorted = new SV_Data(pFeatures->Row, pFeatures->Col);
	pProbe = new SV_Data(1, pFeatures->Col); //*2);
	pSorted->Hdr = pFeatures->Hdr;

	//create a matrix: may frame j follow frame i?
	for (i = 0; i < (unsigned int)(pFeatures->Row); i++) {
		for (j = 0; j < (unsigned int)(pFeatures->Row); j++) {
			if (i == j) {
				mayFollow[i][j] = 0.0;
			} else {
				for (d = 0; d < pFeatures->Col; d++) {
					pProbe->Mat[0][d] = pFeatures->Mat[i][d] - pFeatures->Mat[j][d];
					//pProbe->Mat[0][pFeatures->Col+d] = pFeatures->Mat[j][d];
				}
				result = classify(pProbe, pProbability);

				//if (pProbability != NULL) {
				//	mayFollow[i][j] = pProbability->Mat[0][label2idx(sclib::labelPositive)];
				//} else {
					mayFollow[i][j] = (result[0] == sclib::labelPositive) ? 1.0 : 0.0;
				//}

				MFree_1D(result);
				MFree_0D(pProbability);					
			}
		}
		availableIdx.push_back(i);
	}

	for (i = 0; i < pFeatures->Row; i++) {
		sclib::scalarOut("mayFollow.txt", pFeatures->Mat[i][0], this->pTweak, false, ";");
	}
	sclib::matrixOut("mayFollow.txt", mayFollow, pFeatures->Row, pFeatures->Row, this->pTweak);

	//find a suitable last frame
	for (i = 0; i < availableIdx.size(); i++) {
		for (j = 0; j < availableIdx.size(); j++) {
			allZero = true;
			if (mayFollow[availableIdx[i]][availableIdx[j]] > 0.0) {
				allZero = false;
			}
		} //for j
		if (allZero == true) {
			bestIdx = i;
			break;
		}
	} //for i

	//use the above created matrix to decide, from back to front, which is the next frame in the sorted list
	for (t = pFeatures->Row-1; t >= 0; t--) {
		if (t < pFeatures->Row-1) {
			bestIdx = 0;
			bestRowSum = std::numeric_limits<double>::max();
			for (i = 0; i < availableIdx.size(); i++) {
				rowSum = 0.0;
				for (j = 0; j < availableIdx.size(); j++) {
					rowSum += mayFollow[availableIdx[i]][availableIdx[j]];
				}
				if (mayFollow[availableIdx[i]][lastIdx] > mayFollow[availableIdx[bestIdx]][lastIdx] || rowSum < bestRowSum) { //find the best matching previous frame for the lastIdx' frame; it is best if the prob. is highest or the chance to couple with others is lowest
					bestIdx = i;
					bestRowSum = rowSum;
				}
			}
			if (mayFollow[availableIdx[bestIdx]][lastIdx] <= 0.0) { //the lastIdx' frame has no favoured predecessor, so find a frame also has no connections
				for (i = 0; i < availableIdx.size(); i++) {
					for (j = 0; j < availableIdx.size(); j++) {
						allZero = true;
						if (mayFollow[availableIdx[i]][availableIdx[j]] > 0.0) {
							allZero = false;
						}
						if (allZero == true) {
							bestIdx = i;
							break;
						}
					} //for j
				} //for i
			}
		}
		lastIdx = availableIdx[bestIdx];

		for (d = 0; d < pFeatures->Col; d++) {
			pSorted->Mat[t][d] = pFeatures->Mat[availableIdx[bestIdx]][d];
		}
		availableIdx.erase(availableIdx.begin()+bestIdx);
	}

	MFree_2D(mayFollow);	
	MFree_0D(pProbe);

	return pSorted;
}
*/
