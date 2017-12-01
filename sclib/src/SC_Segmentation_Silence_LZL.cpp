/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates the simple threshold-based silence detector        */
/*        suggested in "Content-based Audio Classification and            */
/*        Segmentation by Using Support Vector Machines", Lu/Zhang/Li     */
/*        2003                                                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.10.2006																								*/
/**************************************************************************/

#include "SC_Aux.h"
#include "SC_Segmentation_Silence_LZL.h"
#include "SC_FeatureHandler.h"
#include "SC_Model_Pareto.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Segmentation_Silence_LZL::SC_Segmentation_Silence_LZL(SC_TweakableParameters* pTweak) : SC_Segmentation_Silence(pTweak) {
  
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Segmentation_Silence_LZL::~SC_Segmentation_Silence_LZL(){

}

//====================================================================================================================
// analyze the features of the given audio-segment and mark them as silence/non-silence in the ground truth
// use energy and zcr per frame as the features
// pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
// the log of the feature-set constants sclib::feature* as indices into the array)
//====================================================================================================================
int SC_Segmentation_Silence_LZL::markSilence(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
  unsigned long int y, idxSTE = sclib::bitPosition(sclib::featureSTE), idxZCR = sclib::bitPosition(sclib::featureZCR), start, end;
  int res = SVLIB_Fail;
	double prob, nrgProb, zcrProb;
	double nrgThresh, zcrThresh;

	//include specificity as provided by the user to alter trained thresholds
	nrgThresh = this->pTweak->segmentationSilenceLzl.energySilenceThreshold * sclib::invSigmoid(this->pTweak->segmentationSilenceLzl.specificity, 0.0, 1.0, 10.0, 0.000001, 0.0); //TODO: use p=0.3 for reasonable sigmoid instead of linear scaling
	zcrThresh = this->pTweak->segmentationSilenceLzl.zcrSilenceThreshold * sclib::invSigmoid(this->pTweak->segmentationSilenceLzl.specificity, 0.0, 1.0, 10.0, 0.000001, 0.0); 

  if (pFeatures[idxSTE] == NULL || pFeatures[idxSTE]->Hdr.ID != sclib::featureSTE ||
		  pFeatures[idxZCR] == NULL || pFeatures[idxZCR]->Hdr.ID != sclib::featureZCR ||
			pFeatures[idxSTE]->Hdr.frameSize != pFeatures[idxZCR]->Hdr.frameSize || pFeatures[idxSTE]->Hdr.frameStep != pFeatures[idxZCR]->Hdr.frameStep || pFeatures[idxSTE]->Row != pFeatures[idxZCR]->Row) {
    REPORT_ERROR(SVLIB_BadArg, "The LNK silence detector needs synchronized STE and ZCR features to operate on!");
	} else {
		//remove all previous silence tags
		pGT->setSegment(segmentStart, segmentEnd, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelRemove);

		/*
		sclib::matrixOutEx("e_zcr.txt", pFeatures[idx]->Mat, pFeatures[idx]->Row, 2, this->pTweak);
		sclib::scalarOut("thresh.txt", this->pTweak->segmentationSilenceLzl.energySilenceThreshold, this->pTweak);
		sclib::scalarOut("thresh.txt", this->pTweak->segmentationSilenceLzl.zcrSilenceThreshold, this->pTweak);
		*/

		for (y = 0; y < (unsigned long int)(pFeatures[idxSTE]->Row); y++) {
			start = segmentStart + (y * pFeatures[idxSTE]->Hdr.frameStep);
      end = start + pFeatures[idxSTE]->Hdr.frameSize - 1;

			//according to Lu, Zhang, Lu, "Content-based Audio Classification and Segmentation by Using Support Vector Machines" (2003),
			//both ZCR and STE must be below a certain threshold
			if (pFeatures[idxSTE]->Mat[y][0] < nrgThresh && //energy
				pFeatures[idxZCR]->Mat[y][0] < zcrThresh) { //zcr
				//use the following as a pseudo-probability for this frame being silence: the normalized distance of the feature 
				//closest to the threshold, therefore responsible for the strength of the decision
				nrgProb = 1.0 - (pFeatures[idxSTE]->Mat[y][0] / nrgThresh);
				zcrProb = 1.0 - (pFeatures[idxZCR]->Mat[y][0] / zcrThresh);
				prob = (min(nrgProb, zcrProb) / 2.0) + 0.5; //scale data between 0.5 and 1.0; 0-0.4999 is reserved for the frames being non-silence

				pGT->setSegment(start, end, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true);
			} else {
				if (pFeatures[idxSTE]->Mat[y][0] < nrgThresh) {
					nrgProb = 1.0 - (pFeatures[idxSTE]->Mat[y][0] / nrgThresh);
				} else {
					nrgProb = (nrgThresh / pFeatures[idxSTE]->Mat[y][0]);
				}
				if (pFeatures[idxZCR]->Mat[y][0] < zcrThresh) {
					zcrProb = 1.0 - (pFeatures[idxZCR]->Mat[y][0] / zcrThresh);
				} else {
					zcrProb = (zcrThresh / pFeatures[idxZCR]->Mat[y][0]);
				}
				prob = min(nrgProb, zcrProb) / 2.0; //scale data between 0 and 0.499999 because we don't have a silence-frame here but want to express the probability for it being silence
			}

			pGT->setProbability(start, end, sclib::atSilence, prob);
    }

    res = SVLIB_Ok;
  }

  return res;
}

//====================================================================================================================
// analyzes the given features of pure silence and prints statistics about the features; suitable threshold values are
// the maximae...
//====================================================================================================================
int SC_Segmentation_Silence_LZL::trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures) {
	unsigned long int idx[2] = {sclib::bitPosition(sclib::featureSTE), sclib::bitPosition(sclib::featureZCR)};
	char *name, *fileName;
	int res = SVLIB_Ok;
	double *tmp;
	SC_MatrixFunctions matFunc;
	SC_FeatureHandler handler(this->pTweak, false);
	SC_Model_Pareto *pModel;

	for (int i = 0; i < 2; i++) {
		name = handler.getFeatureName(sclib::bitValue(idx[i]));
		printf("\nStatistics for feature %s:", name);

		tmp = matFunc.mean(pFeatures[idx[i]]->Mat, pFeatures[idx[i]]->Row, pFeatures[idx[i]]->Col);
		printf("\n\tMean: %f", tmp[0]);
		MFree_0D(tmp);

		tmp = matFunc.std(pFeatures[idx[i]]->Mat, pFeatures[idx[i]]->Row, pFeatures[idx[i]]->Col);
		printf("\n\tSd: %f", tmp[0]);
		MFree_0D(tmp);

		tmp = matFunc.min(pFeatures[idx[i]]->Mat, pFeatures[idx[i]]->Row, pFeatures[idx[i]]->Col);
		printf("\n\tMinimum: %f", tmp[0]);
		MFree_0D(tmp);

		tmp = matFunc.max(pFeatures[idx[i]]->Mat, pFeatures[idx[i]]->Row, pFeatures[idx[i]]->Col);
		printf("\n\tMaximum: %f\n", tmp[0]);
		MFree_0D(tmp);

		printf("\n\tDistribution: see '%s_distribution.txt'\n", name);
		pModel = new SC_Model_Pareto(this->pTweak, NULL);
		pModel->TrainModel(pFeatures[idx[i]]);
		fileName = sclib::exchangeFileExtension(name, "_distribution.txt");
		sclib::classOut(fileName, pModel, this->pTweak);
		MFree_0D(fileName);
		MFree_0D(pModel);
		MFree_0D(name);
	}

	return res;
}

//====================================================================================================================
// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
// are handled therein)
//====================================================================================================================
unsigned long int SC_Segmentation_Silence_LZL::getUncertaintyRegionWidth(void) {
	unsigned long int width = 0;

	//factors due to frame-based analysis
	width += 2 * this->pTweak->featureSte.frameStep;

	return width;
}

//====================================================================================================================
// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
// parameters
//====================================================================================================================
SC_TweakableParameters::SC_FeaturePar* SC_Segmentation_Silence_LZL::getSpecialFeatureParameters(void) {
	//properly link the needed parameter-sets
	this->pTweak->segmentationSilenceLzl.zcrParameters.Next = &(this->pTweak->segmentationSilenceLzl.steParameters);
	this->pTweak->segmentationSilenceLzl.steParameters.Next = NULL;

	return &(this->pTweak->segmentationSilenceLzl.zcrParameters);
}
