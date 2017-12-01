/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the speaker       */
/*        identification according to the ground truth                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 01.05.2006																								*/
/**************************************************************************/

#include <assert.h>
#include "SC_SpeakerScore_Identification.h"
#include "SC_Aux.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_SpeakerScore_Identification::SC_SpeakerScore_Identification(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT, SC_Cluster *pCandidateList, SC_Cluster *pSpeakerList) : SC_SpeakerScore(pTweak, pGT) {
	this->pCandidatePartition = new SC_Partition(pCandidateList, 0.0, 0, true); //just to provide random (id-based) access on the list of speakers
	this->pSpeakerPartition = new SC_Partition(pSpeakerList, 0.0, 0, true);

	this->accuracy = 0.0;
	this->pSpeakerMapping = NULL;

	sprintf(this->purpose, "%s\0", "Speaker-Identification Assessment");
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_SpeakerScore_Identification::~SC_SpeakerScore_Identification() {
	sclib::destructLinkedList(this->pSpeakerMapping);

	this->pSpeakerPartition->setClusters(NULL); //so the linked clutsers don't get destructed together with the partition-container
	MFree_0D(this->pSpeakerPartition);

	this->pCandidatePartition->setClusters(NULL);
	MFree_0D(this->pCandidatePartition);
}

//====================================================================================================================
//	Set a new candidate-list
//====================================================================================================================
void SC_SpeakerScore_Identification::setCandidateList(SC_Cluster *pCandidates) {
	this->pCandidatePartition->setClusters(NULL);
	MFree_0D(this->pCandidatePartition);

	this->pCandidatePartition = new SC_Partition(pCandidates, 0.0, 0, true);
	
	return;
}

//====================================================================================================================
//	Set a new speaker-list
//====================================================================================================================
void SC_SpeakerScore_Identification::setSpeakerList(SC_Cluster *pSpeakers) {
	this->pSpeakerPartition->setClusters(NULL);
	MFree_0D(this->pSpeakerPartition);

	this->pSpeakerPartition = new SC_Partition(pSpeakers, 0.0, 0, true);

	return;
}

//====================================================================================================================
//  Adds a mapping between candidates and speakers
//====================================================================================================================
void SC_SpeakerScore_Identification::addSpeakerMapping(long int speakerId, long int candidateId) {
  SC_SpeakerScore_Identification::SC_IndexMapping *pMapping, *pHook = this->pSpeakerMapping, *pLast = NULL;
	bool found = false;

  //is the mapping already established?
  while (pHook != NULL) {
    if (pHook->idx1 == speakerId && pHook->idx2 == candidateId) {
      found = true;
			break;
    }
    pLast = pHook;
    pHook = pHook->Next;
  }

  //if not found, make a new entry
  if (found == false) {
    pMapping = new SC_SpeakerScore_Identification::SC_IndexMapping(speakerId, candidateId);
    if (pLast == NULL) {
      this->pSpeakerMapping = pMapping;
    } else {
      pLast->Next = pMapping;
    }
  }	

	return;
}

//====================================================================================================================
//	Fills the internal score-variables by computing their values according to the frameList in pGT, so that the
//  get*()-Functions return reasonable values (before calling calcScores(), they return all 0)
//  "start" and "end" refer to sample-numbers so that the area of the frameList for which scores shall be computed can
//  be specified; this way, scores can be calculated only for parts of the video/corpus, e.g. for a scene.
//  In algorithmicUncertaintyDiameter a value [in samples!] can be given which describes the precision with which the
//  specific algorithm responsible for the results (and only the algorithm, not the gt...) can predict the place of 
//  event-on- and -offsets
//====================================================================================================================
void SC_SpeakerScore_Identification::calcScores(unsigned long int start, unsigned long int end, unsigned long int algorithmicUncertaintyDiameter) {
  unsigned long int available, correct;
	long int candidateOwnerID, speakerOwnerID;
	bool idFound; 
  SC_Cluster *pCandidateHook = this->pCandidatePartition->getClusters(), *pSpeakerHook, *pBestFittingSpeaker;

  //init some variables
	this->start = start;
	this->end = (end > 0 && end < this->pGT->getAudioSampleCount()) ? end : this->pGT->getAudioSampleCount()-1; //not so very relevant in the id context...
  available = sclib::getListCount(this->pCandidatePartition->getClusters());
  correct = 0;

  //gather information for each candidate separately
  while (pCandidateHook != NULL) {
		candidateOwnerID = determineClusterOwner(pCandidateHook, algorithmicUncertaintyDiameter);
		pCandidateHook->setOwnerID(candidateOwnerID);
		addSpeakerMapping(pCandidateHook->getBestFitID(), pCandidateHook->getSpeakerID());

		if (pCandidateHook->getBestFitID() == sclib::noSpeaker) { //this candidate was rated an impostor... is that correct (i.e. is there no cluster of this owner in the speaker-list)?

			pSpeakerHook = this->pSpeakerPartition->getClusters();
			idFound = false;

			while (pSpeakerHook != NULL) {
				//we want to know if the identification was correct at the time it took place; 
				//therefor, exclude all segments from influencing the ownerID that didn't belong to that cluster at the time of identification;
				//this segments are those that belong to candidates later in the linked list because thay may only have been included on the cluster later on
				speakerOwnerID = determineClusterOwner(pSpeakerHook, algorithmicUncertaintyDiameter, sclib::ceSample, pCandidateHook);
				if (speakerOwnerID == candidateOwnerID) {
					idFound = true;
					break;
				}
				pSpeakerHook = pSpeakerHook->Next;
			}

			if (idFound == false) {
				correct++;
				this->pGT->addSpeakerMapping(candidateOwnerID, pCandidateHook->getSpeakerID(), true); //add the mapping; because the candidate was corectly identified as an impostor, the correct-flag is set
			} else {
				this->pGT->addSpeakerMapping(candidateOwnerID, pCandidateHook->getSpeakerID(), false); //the cluster should have been identified as one of the speakers, so the correct-flag is ==false
			}

		} else { //the candidate was identified to belong to the best-fitting speaker... do they have the same owner?
			pBestFittingSpeaker = this->pSpeakerPartition->getClusterByID(pCandidateHook->getBestFitID()); //this cluster must exist!
			speakerOwnerID = determineClusterOwner(pBestFittingSpeaker, algorithmicUncertaintyDiameter, sclib::ceSample, pCandidateHook); //see above for explanations on the exlude-list
			if (speakerOwnerID == candidateOwnerID) {
				correct++;
				this->pGT->addSpeakerMapping(speakerOwnerID, pCandidateHook->getSpeakerID(), pBestFittingSpeaker->getOwnerIsCorrect()); //add the mapping with the correct-status of the identified best speaker (assumes that clustering-scoring was performed beforehand)
			} else {
				this->pGT->addSpeakerMapping(speakerOwnerID, pCandidateHook->getSpeakerID(), false); //the identification failed, so correct==false
			}
		}

    pCandidateHook = pCandidateHook->Next;
  }

  //calc score
	calcIdentificationScores(available, correct, this->accuracy);

  return;
}


//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_SpeakerScore_Identification::output(ostream& OutS) {
	char *speakerName, buffer[sclib::bufferSize];
	bool isCorrect;
	SC_SpeakerScore::SC_IndexMapping *pHook = this->pSpeakerMapping;
	SC_Cluster *pCandidateHook; 
	
	OutS << "\n";
	OutS << "Identification scores:\n";
	OutS << "-----------------------\n";
	OutS << "Accuracy:         " << this->accuracy << "\n";
	OutS << "Recognition rate: " << "see accuracy\n";

	OutS << "\n";
	OutS << "Speaker mapping:\n";
	OutS << "-----------------\n";
	OutS << "An asterisk (*) after a speaker-name indicates that the candidate was not identified as\n"
	     << "belonging to the biggest cluster for that speaker; but, at least, it was not identified as\n"
			 << "a totally different speaker.\n";

	OutS << "\nCandidate ID    mapped to ID    =>hypothesized speaker    =>real speaker\n";
	while (pHook != NULL) {
		OutS << left << setw(16) << pHook->idx2;

		sprintf(buffer, "%d\0", pHook->idx1);
		OutS << left << setw(16) << ((pHook->idx1 != sclib::noSpeaker) ? buffer : "impostor");

		speakerName = const_cast<char*>(this->pGT->getSpeakerName(this->pGT->getSpeakerGIDFromHID((pHook->idx1 != sclib::noSpeaker) ? pHook->idx1 : pHook->idx2, isCorrect))); //if the speaker was identified as an impostor, look for a mapping with his own id to a gt-id; otherwise, look for a mapping between the bestFitID and a gt-id
		sprintf(buffer, "%s%s\0", speakerName, ((isCorrect == true) ? "" : "*"));
		OutS << left << setw(26) << buffer;

		pCandidateHook = this->pCandidatePartition->getClusterByID(pHook->idx2);
		speakerName = const_cast<char*>(this->pGT->getSpeakerName(pCandidateHook->getOwnerID()));
		OutS <<	speakerName;

		OutS << "\n";

		pHook = pHook->Next;
	}

  return OutS;
}
