/**************************************************************************/
/*	A Partition stores the result of a clustering algorithm 							*/
/*																																				*/
/*  Attention: It really holds copys of the clusters, so it's very        */
/*             memory-consuming if the clusters where not created with    */
/*             the 'justLink'-option                                      */
/*																																				*/
/*    Responsibility:																											*/
/*      - encapsulates a complete partition:                            	*/
/*				- the linked list of clusters     															*/
/*				- methods to organize them  																		*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.09.2005																								*/
/**************************************************************************/

#include "SC_Partition.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor to maintain a copy of the cluster-list
//  if justLink=true, the new objects gets just a pointer to the given list of clusters; otherwise, it makes a copy
//  and stores it. if noWorkData==true (and justLink==false), the speech-frames and models in the clusters aren't 
//  disposed at all (like calling the constructor and realeaseWorkData() directly after each other)
//====================================================================================================================
SC_Partition::SC_Partition(SC_Cluster* pClusters, double globalCriterion, unsigned long int partitionNr, bool justLink, bool noWorkData) : partitionNr(partitionNr) {
  if (justLink == true) {
    this->pClusters = pClusters;
  } else {
		this->pClusters = sclib::copyLinkedList(pClusters, true, noWorkData);
	}

  this->globalCriterion = globalCriterion;
  this->classSig = 196274;
	this->Next = NULL;
	this->pFather = NULL;
	this->pMother = NULL;
	this->pChild = NULL;
}

//====================================================================================================================
//	copy-constructor, see above for comments
//====================================================================================================================
SC_Partition::SC_Partition(const SC_Partition& pParent, bool justLink) {
  if (justLink == true) {
    this->pClusters = pParent.pClusters;
  } else {
    this->pClusters = sclib::copyLinkedList(pParent.pClusters);
  }

  this->globalCriterion = pParent.globalCriterion;
  this->classSig = 196274;
	this->Next = pParent.Next;
	this->pFather = pParent.pFather;
	this->pMother = pParent.pMother;
	this->pChild = pParent.pChild;
	this->partitionNr = pParent.partitionNr;
}

//====================================================================================================================
//	destructor
//  the clusters get destructed, regardless if they where copys or just linked
//====================================================================================================================
SC_Partition::~SC_Partition() {
  sclib::destructLinkedList(this->pClusters);
  this->pClusters = NULL;
  this->Next = NULL;
  this->classSig = 0;
}

//====================================================================================================================
//	Needed by sclib::destructLinkedList()
//====================================================================================================================
bool SC_Partition::Valid(void) {
  return (this->classSig == 196274) ? true: false;
}

//====================================================================================================================
//	Return the pointer to the (first) cluster in the linked list of clusters, which has the given speaker name
//====================================================================================================================
SC_Cluster* SC_Partition::getClusterByName(const char* speakerName) {
  SC_Cluster *pRes = NULL, *pHook = this->pClusters;

  while (pHook != NULL) {
    if (strncmp(speakerName, pHook->getSpeakerName(), strlen(speakerName)) == 0) {
      pRes = pHook;
      break;
    }

    pHook = pHook->Next;
  }

  return pRes;
}
    
//====================================================================================================================
//	Return the pointer to the (first) cluster in the linked list of clusters, which has the given (hypothesized) 
//  speaker ID
//====================================================================================================================
SC_Cluster* SC_Partition::getClusterByID(unsigned long int speakerID) {
  SC_Cluster *pRes = NULL, *pHook = this->pClusters;

  while (pHook != NULL) {
    if (speakerID == pHook->getSpeakerID()) {
      pRes = pHook;
      break;
    }

    pHook = pHook->Next;
  }

  return pRes;
}

//====================================================================================================================
//	Return the index into the linked list of the cluster with given ID; SVLIB_Fail, if not found
//====================================================================================================================
int SC_Partition::getClusterIdxByID(unsigned long int speakerID) {
	int c = 0, idx = SVLIB_Fail;
	SC_Cluster *pHook = this->getClusters();

	while (pHook != NULL) {
		if (pHook->getSpeakerID() == speakerID) {
			idx = c;
			break;
		}
		c++;
		pHook = pHook->Next;
	}

	return idx;
}

//====================================================================================================================
// to store the ID-information found in the cluster-objects in the frameList of the groundtruth-object
// that means: remove all previous speaker-boundary/speech-seg-start/speech-seg-end/speaker-id labels and set them 
// according to the information stored in the clusters of this partition
//====================================================================================================================
void SC_Partition::storeIDs(SC_GroundTruth *pGT, bool removeOldSpeakerBoundarys, bool removeOldSpeakerIDs, bool removeOldSpeechSegments) {
  unsigned long int x;
  long int segmentStart, segmentEnd;
  long int currentSpeakerID, oldSpeakerID = sclib::noSpeaker;
  SC_Cluster* pHook;

  //remove all speaker-boundary-labels, speaker-id-tags and speech-segment-start/end-tags, if wished
  if (removeOldSpeakerBoundarys == true) {
    pGT->setSegment(0, pGT->getAudioSampleCount(), sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized);
  }
  if (removeOldSpeakerIDs == true) {
    pGT->setSegment(0, pGT->getAudioSampleCount(), sclib::noType, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized, true);
  }
  if (removeOldSpeechSegments == true) {
    pGT->setSegment(0, pGT->getAudioSampleCount(), sclib::atSpeechSegmentStart, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized);
    pGT->setSegment(0, pGT->getAudioSampleCount(), sclib::atSpeechSegmentEnd, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized);
  }

  //1. store the segment-boundaries and speaker-ids
  //the segments here are maybe composed of smaller subsegments due to the speaker change detection; this gives rise to the 3rd block
	pHook = this->getClusters();
	while (pHook != NULL) {
		for (x = 0; x < pHook->getSegmentCount(); x++) {
      if (pHook->getSpeakerID() != sclib::noSpeaker) {
        pGT->setSegment(pHook->getSegmentList(x, 1), pHook->getSegmentList(x, 2), sclib::noType, false, pHook->getSpeakerID()); //speaker-id for the whole segment
				pGT->setSegment(pHook->getSegmentList(x, 1), pHook->getSegmentList(x, 1), sclib::atSpeechSegmentStart); //speech-seg-start
				pGT->setSegment(pHook->getSegmentList(x, 2), pHook->getSegmentList(x, 2), sclib::atSpeechSegmentEnd); //speech-seg-end (all speaker-id's have been set in the 1st block)
      } else {
        pGT->setSegment(pHook->getSegmentList(x, 1), pHook->getSegmentList(x, 2), sclib::noType, false, pHook->getSpeakerID(), sclib::modeLabelRemove, sclib::modeHypothesized, true); //remove the speaker-id for the whole segment
      }
		}
		pHook = pHook->Next;
	}

  //2. store the speaker-boundaries
  for (x = 0; x <= pGT->getAudioSampleCount(); x++) {
    pGT->getNextSegment(x, segmentStart, segmentEnd, sclib::atSpeechSegmentStart);
    if ((segmentStart != sclib::noSegment) && (segmentEnd != sclib::noSegment)) {
      pGT->setSegment(segmentStart, segmentEnd, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized); //remove speaker-boundaries within this segment
      currentSpeakerID = pGT->getSamplesSpeakerID(segmentStart, sclib::modeHypothesized);
      if (currentSpeakerID != oldSpeakerID && currentSpeakerID != sclib::noSpeaker) {
        pGT->setSegment(segmentStart, segmentStart, sclib::atSpeakerBoundary); //speaker-boundary (speaker-change)
        oldSpeakerID = currentSpeakerID;
      }
      x = segmentEnd;
    } else {
      break;
    }
  }

  //3. give each separated speech-segment it's own speech-segment-start/end
  //the effect of this: in gardener (the vb-frontend), you see each single short speech segment 
  //and not just one big segment from speaker-boundary to speaker-boundary
  for (x = 0; x <= pGT->getAudioSampleCount(); x++) {
		pGT->getNextSegment(x, segmentStart, segmentEnd, sclib::atSpeech, sclib::searchForward, sclib::modeHypothesized, true);
    if ((segmentStart != sclib::noSegment) && (segmentEnd != sclib::noSegment)) {
      pGT->setSegment(segmentStart, segmentStart, sclib::atSpeechSegmentStart); //speech-seg-start
      pGT->setSegment(segmentEnd, segmentEnd, sclib::atSpeechSegmentEnd); //speech-seg-end (all speaker-id's have been set in the 1st block)
      x = segmentEnd;
    } else {
      break;
    }
  }

  return;
}

//====================================================================================================================
//	return the uncertainty diameter due to clustering, i.e. two framesteps of the used features
//====================================================================================================================
unsigned long int SC_Partition::getUncertaintyDiameter(void) {
	unsigned long int diameter = 0;
	SC_Cluster *pHook = this->pClusters;
	SV_Data *pFrames;

	while (pHook != NULL) {
		pFrames = pHook->getSpeechFrames();
		if (pFrames != NULL) {
			diameter = pFrames->Hdr.frameStep * 2;
			break;
		}

		pHook = pHook->Next;
	}

	return diameter;
}

//====================================================================================================================
//	destructs all the basic data-structures ((background-)models and speechFrames) of all clusters in this partition
//  this saves much of memory when all computations are done and the object is just holded to store the results in a 
//  partition
//====================================================================================================================
void SC_Partition::releaseWorkData(void) {
	SC_Cluster *pHook = this->pClusters;

	while (pHook != NULL) {
		pHook->releaseWorkData();
		pHook = pHook->Next;
	}

	return;
}

//====================================================================================================================
//	if this is not the first partition in a succession of partitions created during clustering, it doesn't hold any
//  work data, e.g. speech frames; but the first partition of the process does; this method returns, for a given id
//  of one of its clusters, a linked list of *new* objetcs (the objects are new and need destruction afterwards, the 
//  data itself is just linked to save space) representing all frames for this cluster, which can be found in this
//  first partition via the segment-boundaries, which don't get changed during clustering
//====================================================================================================================
SV_Data* SC_Partition::getSpeechFramesOfCluster(unsigned long int speakerID, SC_Partition *pFirstPartition) {
	SV_Data *pSpeechFrames = NULL, *pSpeechHook;
	SC_Cluster *pCluster = getClusterByID(speakerID), *pClusterHook;
	unsigned long int start, end;
	bool found;

	for (unsigned long int i = 0; i < pCluster->getSegmentCount(); i++) {
		//we need this segment from the given cluster
		start = pCluster->getSegmentList(i, 1);
		end = pCluster->getSegmentList(i, 2);

		pClusterHook = pFirstPartition->getClusters();
		found = false;
		while (pClusterHook != NULL) {
			for (unsigned long int j = 0; j < pClusterHook->getSegmentCount(); j++) {
				if (pClusterHook->getSegmentList(j, 1)==start && pClusterHook->getSegmentList(j, 2)==end) {
					//we found it when we find a segment in the first partition with corresponding boundaries
					found = true;
					if (pSpeechFrames == NULL) {
						pSpeechFrames = new SV_Data(*pClusterHook->getSpeechFrames(j), true);
						pSpeechHook = pSpeechFrames;
					} else {
						pSpeechHook->Next = new SV_Data(*pClusterHook->getSpeechFrames(j), true);
						pSpeechHook = pSpeechHook->Next;
					}
					break;
				}
			}
			if (found == true) {
				break;
			}
			pClusterHook = pClusterHook->Next;
		}

		if (found == false) { //lament on failure -> shouldn't happen at all, but for a bug in prvious stages (clustering)
			REPORT_ERROR(SVLIB_Fail, "can't find segment in first partition");
		}
	}

	return pSpeechFrames;
}

//====================================================================================================================
//	as above, for all clusters in the current partition; the speech frames are directly stored in the corresponding
//  cluster objects
//====================================================================================================================
void SC_Partition::getSpeechFrames(SC_Partition *pFirstPartition) {
	SC_Cluster *pHook = this->pClusters;
	SV_Data *pFrames;

	if (this != pFirstPartition) { //can't copy from oneself!
		while (pHook != NULL) {
			if (pHook->getSpeechFrames() == NULL) { //only create new speech frames if old ones aren't present
				pFrames = getSpeechFramesOfCluster(pHook->getSpeakerID(), pFirstPartition);
				pHook->setSpeechFrames(pFrames);
			}
			pHook = pHook->Next;
		}
	}

	return;
}
