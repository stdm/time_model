/**************************************************************************/
/*    Responsibility:																											*/
/*      - Performs the speaker idetification of a linked list of features */
/*        against a set of speaker-clusters                               */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 26.04.2006																								*/
/**************************************************************************/

#include "SC_SpeakerIdentificator.h"
#include "SC_DistanceMeasures.h"
#include "SC_ModelHandler.h"
#include "SC_Partition.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_SpeakerIdentificator::SC_SpeakerIdentificator(SC_TweakableParameters *pTweak, SC_MixtureModel *pBackgroundSpeakerModel, bool verbose){
  this->pTweak = pTweak;
  this->pBackgroundSpeakerModel = pBackgroundSpeakerModel;
  this->verbose = verbose;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_SpeakerIdentificator::~SC_SpeakerIdentificator(){

}

//====================================================================================================================
// Performs score normalization (accounting for unknown impostors and imperfect features/models) under the "distance 
// normalization" paradigm published in "A Monte-Carlo Method for Score Normalization in Automatic Speaker 
// Verification using Kullback-Leibler Distances" by Ben, Blouet, Bimbot on ICASSP 2002
// DNorm is comparable (in terms of recognition performance) to ZNorm and TNorm but dosn't need extra training data
// or speaker models to learn the normalization parameters.
//====================================================================================================================
double SC_SpeakerIdentificator::dNorm(double clientScore, SC_Model *pClientModel, SC_MixtureModel *pWorldModel) {
  SV_Data *pClientData, *pImpostorData;
  SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak, false);
  double normalizedScore, KLc, KLw, KL2, pcc, pcw, pww, pwc;

	//TODO: doesn't work well... problem with 2nd normalization (after bg-normalization) or log-likelihoods and division?!?

  //draw samples from the two distributions in a monte-carlo manner
  pClientData = pClientModel->drawSamplesFromDistribution(this->pTweak->speakerIdentification.dNormSampleCount);
  pImpostorData = pWorldModel->drawSamplesFromDistribution(this->pTweak->speakerIdentification.dNormSampleCount);

  //calc different scores
  pcc = pModelHandler->testModel(pClientModel, pClientData, 1, true, NULL, false);
  pwc = pModelHandler->testModel(pWorldModel, pClientData, 1, true, NULL, false);
  pww = pModelHandler->testModel(pWorldModel, pImpostorData, 1, true, NULL, false);
  pcw = pModelHandler->testModel(pClientModel, pImpostorData, 1, true, NULL, false);

  MFree_0D(pClientData);
  MFree_0D(pImpostorData);
  MFree_0D(pModelHandler);

  //compute final (symmetrized) KL-distance (this is essentially the CLR...)
  KLc = pcc - pwc;
  KLw = pww - pcw;
  KL2 = KLc + KLw;

  //do normalization
  normalizedScore = clientScore / KL2;

  return normalizedScore;
}

//====================================================================================================================
// Identify to which cluster (if any) each feature-set in the linked list of candidates belongs; the candidates are 
// also organized in "clusters" to enhance the features with segment-border information (no speaker-models are needed
// therein!); the speakerList of clusters is assumed to contain *real* clusters as obtained by the speaker-clusterer 
// object
// updateClusters should be set ==true so that correct speaker-information is stroed in the groundtruth
// TODO: - threshold-finding to exclude "impostors"; 0.0 is just the ideal threshold when the world-model scores 
//         higher than any speaker-model...
//       - at the moment, identification is performed based on the merged speaker model; integrate id based on other
//         linkage methods, such as average, complete or single linkage based on single segment models
//====================================================================================================================
SC_SpeakerScore_Identification* SC_SpeakerIdentificator::identifySpeakers(SC_GroundTruth *pGT, SC_Cluster *pCandidateList, SC_Cluster *pSpeakerList, double impostorThreshold, bool updateClusters, bool retrainClusters, long int progressCount, long int progressMax, bool (*doProceed)(long int &curr, long int max, unsigned int cnt)) {
  SC_Cluster *pCandidateHook = pCandidateList, *pSpeakerHook, *pBestSpeaker, *pListEnd = NULL;
  SC_Partition *pPartition = NULL;
  SC_ModelHandler *pModelHandler = new SC_ModelHandler(this->pTweak, false);
  SC_SpeakerScore_Identification *pScore = NULL;
  double backgroundScore, bestScore, score;
  bool needBackgroundScore;
	bool immediateReturn = false;
  unsigned long int x, uncertaintyDiameter;

	//make sure merged speaker models exist
	unsigned int oldMm = this->pTweak->cluster.mergeMode;
	this->pTweak->cluster.mergeMode = sclib::mergeRetrain;
	pSpeakerHook = pSpeakerList;
	while (pSpeakerHook != NULL) { //build GMM-UBM models for all speakers found during clustering because they are best for id with short segments as in pallSHORTclusters
		if (this->pTweak->speakerIdentification.useUBMs == true && (pSpeakerList->getMergedModel()==NULL || pSpeakerList->getMergedModel()->Hdr.ModelType!=sclib::mtGMM_UBM)) {
			pSpeakerHook->retrainMergedModel(sclib::mtGMM_UBM); //convert all speaker-models to UBMs if wished
		} else if (pSpeakerHook->getMergedModel() == NULL) {
			pSpeakerHook->reBuildMergedModel();
		}
		pSpeakerHook = pSpeakerHook->Next;
	}
	this->pTweak->cluster.mergeMode = oldMm;
	needBackgroundScore = (this->pBackgroundSpeakerModel == NULL || pSpeakerList->getMergedModel()->Hdr.ModelType == sclib::mtGMM_UBM) ? false : true; //we don't need extra background-speaker normalization if this is done within the model itself or if we don't have any background-speaker model...

  while (pCandidateHook != NULL) { //loop over all speech-segments for which speaker-identity is to be identified
    //calculate background-score, if needed
    if (needBackgroundScore == true) {
      backgroundScore = pModelHandler->testModel(this->pBackgroundSpeakerModel, pCandidateHook->getSpeechFrames(), pCandidateHook->getSegmentCount(), false);
    }

    //init some variables for this loop
    pSpeakerHook = pSpeakerList;
    bestScore = -1.0 * std::numeric_limits<float>::max();
    pBestSpeaker = NULL;

    //find to which speaker the current frames are most similar
    while (pSpeakerHook != NULL) { //loop over all possible speakers
      if (pSpeakerHook->getMergedModel() != NULL) { //there can be clusters without a speaker model due to the updateClusters parameter
        score = pModelHandler->testModel(pSpeakerHook->getMergedModel(), pCandidateHook->getSpeechFrames(), pCandidateHook->getSegmentCount(), false);
        if (needBackgroundScore == true) {
          score -= backgroundScore; //the scores are log-likeliohoods, so '-' instead of '/' here for correct arithmetic
        }

				//normalize the score to cope with varying conditions/unoptimal models/impostors, if wished
				switch (this->pTweak->speakerIdentification.normalizationMode) {
					case sclib::normalizationDNorm: {
						if (this->pBackgroundSpeakerModel != NULL) {
							score = dNorm(score, pSpeakerHook->getMergedModel(), this->pBackgroundSpeakerModel);
						} else {
							REPORT_ERROR(SVLIB_BadArg, "Background-Model for score-normalization not available");
						}
						break;
					}
					case sclib::normalizationNone: {
						break;
					}
				}				
				
				if (score > bestScore) { //found a new best-scoring speaker
          bestScore = score;
          pBestSpeaker = pSpeakerHook;
        }
        if (this->verbose == true) {printf(".");}
      }

      pSpeakerHook = pSpeakerHook->Next;
    } //loop over all speakers

    pCandidateHook->setIDconfidenceScore(bestScore);
    if (bestScore > impostorThreshold) { //no "impostor", meaning the segment was spoken by the best-scoring known speaker
      pCandidateHook->setBestFitID(pBestSpeaker->getSpeakerID()); //remember best-scoring speaker in the cluster-object of the candidate
      if (updateClusters == true) {
        for (x = 0; x < pCandidateHook->getSegmentCount(); x++) {
          pBestSpeaker->addSegment(pCandidateHook->getSpeechFrames(x), pCandidateHook->getSegmentList(x, 1), pCandidateHook->getSegmentList(x, 2), true, false);
        }
				if (retrainClusters == true) {
					pBestSpeaker->retrainMergedModel(pBestSpeaker->getMergedModel()->Hdr.ModelType);
				}
      }
    } else {//otherwise, an "impostor" is found, meaning the segment was spoken by a previously unknown speaker; because of that, the previous speaker-id/speaker-name of that cluster just remains being active
			pCandidateHook->setBestFitID(sclib::noSpeaker); //remember that this cluster was found to correspond with none from the speaker-list...
      if (updateClusters == true) {
        pListEnd = sclib::getLastInList(pSpeakerList);
        pListEnd->Next = new SC_Cluster((SC_Cluster&)*pCandidateHook, true);
        pListEnd->Next->Next = NULL;
        if (retrainClusters == true) {
          pListEnd->Next->retrainMergedModel(pSpeakerList->getMergedModel()->Hdr.ModelType); //TODO: this forces model-building even for too short segments!
        }
      }  
    }

    pCandidateHook = pCandidateHook->Next;
		if (doProceed != NULL && progressMax > 0) { //care for proper progress display in SC_Lib::audioSegmentation(), which is called by Videana
			immediateReturn = (*doProceed)(progressCount, progressMax, 1);
			if (immediateReturn == true) {
				break;
			}
		}
    if (this->verbose == true) {
			printf(":");
		}
  } //loop over all candidates

	if (immediateReturn == false) {
		//store the results in the GT
		if (updateClusters == true) {
			pPartition = new SC_Partition(pSpeakerList, 0.0, 0, true); //all information is contained in the continuously updated speaker-list
		} else {
			pCandidateHook = sclib::getLastInList(pSpeakerList);
			pCandidateHook->Next = pCandidateList;
			pPartition = new SC_Partition(pSpeakerList, 0.0, 0, true); //all information is present if the speaker- and candidate-list are concatented
		}
		pPartition->storeIDs(pGT, true, true, true);
		if (updateClusters == false) {
			pCandidateHook->Next = NULL; //release the concatenation of the speaker- and candidate-list
		}

		//calculate scores
		pScore = new SC_SpeakerScore_Identification(this->pTweak, pGT, pCandidateList, pSpeakerList);
		uncertaintyDiameter = pPartition->getClusters()->getSpeechFrames()->Hdr.frameStep * 2; 
		pScore->calcScores(0, 0, uncertaintyDiameter);

		pPartition->setClusters(NULL); //so that the destruction afterwards doesn't kill the clusters in pCandidateList
		MFree_0D(pPartition);
	}

	MFree_0D(pModelHandler);

  return pScore;
}
