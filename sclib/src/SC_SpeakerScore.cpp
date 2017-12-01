/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Abstarct base class for speaker-related scores (because there   */
/*        is some common functionality among them													*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.12.2006																								*/
/**************************************************************************/

#include <assert.h>
#include "SC_SpeakerScore.h"
#include "SC_Aux.h"
#include "SC_MatrixFunctions.h"

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_SpeakerScore::SC_SpeakerScore(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT) : SC_Score(pTweak, pGT) {
  this->pClassMapping = NULL;
	for (unsigned int ce = 0; ce < 4; ce++) {
		this->pClusterMapping[ce] = NULL;
	}

	this->ca = pGT->getRealSpeakerCount(true);
}

//====================================================================================================================
//	The destructor
//  destructs all the linked partitions, too!!!
//====================================================================================================================
SC_SpeakerScore::~SC_SpeakerScore() {
  sclib::destructLinkedList(this->pClassMapping);

	for (unsigned int ce = 0; ce < 4; ce++) {
		sclib::destructLinkedList(this->pClusterMapping[ce]);
	}
}

//====================================================================================================================
//	Convert a class-label (i.e. gt-speaker-id's) to a valid col-index into the scatter-matrix; SVLIB_Fail is returned 
//  if the mapping can't be established
//====================================================================================================================
long int SC_SpeakerScore::class2idx(long int classTag) {
	long int res = SVLIB_Fail;
	unsigned long int count = 0;
  SC_SpeakerScore::SC_IndexMapping *pMapping, *pHook = this->pClassMapping, *pLast = NULL;

  //find the index
  while (pHook != NULL) {
    if (pHook->idx1 == classTag) {
      res = pHook->idx2;
			break;
    }
    count++;
    pLast = pHook;
    pHook = pHook->Next;
  }

  //if not found, make a new entry and return the new index
  if (res == SVLIB_Fail) { //this maybe creates more mappings than there are real classes (this->ca); it is ok because it may be called from a method called by createMergeHistory(), and in the history there will most likely exist more clusters than classes
    pMapping = new SC_SpeakerScore::SC_IndexMapping(classTag, count);
    if (pLast == NULL) {
      this->pClassMapping = pMapping;
    } else {
      pLast->Next = pMapping;
    }
    res = pMapping->idx2;
  }
  
	return res;
}

//====================================================================================================================
//	Convert a valid col-index of the scatter-matrix into a gt-speaker-id; SVLIB_Fail is returned on failure (unknown
//  speaker-id); return (on wish) the speaker-id as the shortName or the speaker-realname as the className
//====================================================================================================================
long int SC_SpeakerScore::idx2class(unsigned long int classIdx, char* className, bool shortName) {
	long int res = SVLIB_Fail;
	SC_SpeakerScore::SC_IndexMapping *pHook = this->pClassMapping;

	while (pHook != NULL) {
		if (pHook->idx2 == classIdx) {
			res = pHook->idx1;
			if (className != NULL) {
				if (shortName == true) {
					sprintf(className, "%d\0", pHook->idx1); //take the speaker-id as a short-name
				} else {
					sprintf(className, "%s\0", ((pHook->idx1 == sclib::noSpeaker) ? "NO_SPEAKER" : this->pGT->getSpeakerName(pHook->idx1))); //care for a proper name of the "no-speaker"
				}
			}
			break;
		}
		pHook = pHook->Next;
	}

	return res;
}

//====================================================================================================================
//  Returns true if the given segment-borders do not appear in any of the clusters in the exlude-list
//====================================================================================================================
bool SC_SpeakerScore::segmentIsUnknwon(long int segmentStart, long int segmentEnd, SC_Cluster *pExludeList) {
	bool res = true;
	SC_Cluster *pHook = pExludeList;

	while (pHook != NULL) {
		for (unsigned long int x = 0; x < pHook->getSegmentCount(); x++) {
			if (segmentStart == pHook->getSegmentList(x, 1) && segmentEnd == pHook->getSegmentList(x, 2)) {
				res = false;
				break;
			}
		}

		if (res == false) {
			break;
		}

		pHook = pHook->Next;
	}

	return res;
}

//====================================================================================================================
//	Determines and returns the owner (gt-speaker-)id of the given cluster; if the linked list of cluster-objects in
//  pSegmentExludeList is not empty, only segments from pClutser are considered for determining the owner that do not
//  occur in one of the clusters of that list (needed to do id-scoring); in differentCa, a class-count different to
//  the member this->ca can be given to assure correct array sizes if called by a method that temporarily manages
//  more classes (as createMergeHistory() in SC_SpeakerScore_Clustering does)
//====================================================================================================================
long int SC_SpeakerScore::determineClusterOwner(SC_Cluster *pCluster, unsigned long int algorithmicUncertaintyDiameter, int countedEntity, SC_Cluster *pSegmentExcludeList, unsigned long int differentCa) {
	unsigned long int x, y, scene = 0, lastScene = 0, shot = 0, lastShot = 0, **orderedSegments = NULL, softBoundaryDiameter = this->pGT->getUncertaintyRegionWidth(true) + algorithmicUncertaintyDiameter, newCa = sclib::max(this->ca, differentCa);
	long int hypoSegStart, hypoSegEnd, gtSegStart, gtSegEnd, oldGtSegStart = sclib::noSegment, gtID, hypoID, idxX, ownerID = sclib::noSpeaker, *scatterVector = NULL;
	bool *hits = NULL;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

	MArray_1D(scatterVector, newCa, long int, "SC_SpeakerScore.determineClusterOwner: scatterVector");
	if (countedEntity == sclib::ceShot || countedEntity == sclib::ceScene) {
		hits = pFunc->initVector(newCa, false);
	}
	
	pCluster->uniteNeighboringSegments(); //to allow correct segment-based scoring
	orderedSegments = pCluster->getOrderedSegments(); //we need a ordeed list of segments for correct scene/shot-handling
	hypoID = pCluster->getSpeakerID();

	for (y = 0; y < pCluster->getSegmentCount(); y++) {
		hypoSegStart = orderedSegments[y][0];
		hypoSegEnd = orderedSegments[y][1];

		if (segmentIsUnknwon(hypoSegStart, hypoSegEnd, pSegmentExcludeList) == true) {
			for (x = hypoSegStart; x <= (unsigned long int)(hypoSegEnd); x++) {
				this->pGT->getNextSegment(x, gtSegStart, gtSegEnd, sclib::atSpeech, sclib::searchWithin, sclib::modeGroundtruth, true, true, sclib::atShotBoundary);
				
				if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment && gtSegStart <= hypoSegEnd && gtSegStart <= gtSegEnd) {
					gtID = this->pGT->getSamplesSpeakerID(gtSegStart, sclib::modeGroundtruth);
					if (gtID != sclib::noSpeaker) {
						idxX = class2idx(gtID);

						//care for soft boundaries
						if (abs(gtSegStart - hypoSegStart) < (long int)(softBoundaryDiameter)) {
							gtSegStart = hypoSegStart;
						}
						if (abs(gtSegEnd - hypoSegEnd) < (long int)(softBoundaryDiameter)) {
							gtSegEnd = hypoSegEnd;
						}
						
						//truncate overlapping segments
						if (countedEntity != sclib::ceSegment && gtSegEnd > hypoSegEnd) { 
							gtSegEnd = hypoSegEnd;
						}

						//fill the scatter-vector of the desired CE
						if (countedEntity == sclib::ceSample) {
							scatterVector[idxX] += gtSegEnd - gtSegStart + 1;
						} else if (countedEntity == sclib::ceSegment) {
							scatterVector[idxX] += (gtSegStart == hypoSegStart && gtSegEnd == hypoSegEnd) ? 1 : 0;
						} else if (countedEntity == sclib::ceShot) {
							shot = this->pGT->sample2shot(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastShot);
							if (shot != lastShot) {
								pFunc->clear(&hits, 1, newCa, false); //clear rememberance of hits
							}
							if (hits[idxX] == false) {
								scatterVector[idxX] += 1;
								hits[idxX] = true;
							}
						} else if (countedEntity == sclib::ceScene) {
							scene = this->pGT->sample2scene(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastScene);
							if (scene != lastScene) {
								pFunc->clear(&hits, 1, newCa, false); //clear rememberance of hits
							}
							if (hits[idxX] == false) {
								scatterVector[idxX] += 1;
								hits[idxX] = true;
							}
						}

					} //no gt-speaker
				} else { //no more gt-segments within this cluster-segment
					break;
				}

				oldGtSegStart = gtSegStart;
				lastScene = scene;
				lastShot = shot;
				x = gtSegEnd;
			} //for all gt-segments within a cluster-segment
		} //current hypo-segment is unknown in the exclude-segment-list
	} //for all cluster-segments

	idxX = pFunc->maxIdx(scatterVector, newCa);
	ownerID = class2idx(idxX);
	
	MFree_0D(pFunc);
	MFree_1D(scatterVector);
	MFree_1D(hits);
	MFree_2D(orderedSegments);

	return ownerID;
}

//====================================================================================================================
//	Returns the number of segments that belong to clutsers with the given owner-id
//====================================================================================================================
unsigned long int SC_SpeakerScore::getNrOfSegmentsWithOwner(SC_Cluster *pClusters, long int ownerID) {
	unsigned long int res = 0;
	SC_Cluster * pHook = pClusters;

	while (pHook != NULL) {
		if (pHook->getOwnerID() == sclib::noSpeaker) {
			pHook->setOwnerID(determineClusterOwner(pHook));
		}
		if (pHook->getOwnerID() == ownerID) {
			res += pHook->getSegmentCount();
		}
		pHook = pHook->Next;
	}
	
	return res;
}

//====================================================================================================================
//  Calculate the Diarization error Rate (DER) according to Huang, Marchertet, Visweswariah, Potamianos, "The IBM RT07 
//  Evaluation System for Speaker Diarization on Lecutre Meetings", 2007:
//    "In accordance to NIST scoring, results are reported in terms of diarization error rate (DER). DER is calculated 
//     by first finding the optimal one-to-one mapping between reference speakers and the hypothesized ones, and then 
//     computing the percentage of time that is wrongly assigned according to the optimal mapping. DER includes 
//     speaker error time, missed speaker time, and false alarm speaker time, thus also taking SAD errors into 
//     account" (SAD=speech activity detection)
//  ATTENTION: Computation moved to SC_GroundTruth for speed reasons
//====================================================================================================================
double SC_SpeakerScore::calcDER(unsigned long int softBoundaryDiameter) {
	return this->pGT->calcDER(softBoundaryDiameter);
	/*unsigned long int wrong = 0;
	long int gtId, hypoId, hypoGtId, start, end, speechStart, speechEnd; //, tmp;
	bool correct;

	speechStart = this->pGT->getNextOccasion(0, sclib::atSpeechSegmentStart, sclib::searchForward, sclib::modeHypothesized, sclib::alignmentStart); //get the boundaries for the first speech segment: segment-start...
	speechEnd = this->pGT->getNextOccasion(speechStart, sclib::atSpeechSegmentEnd, sclib::searchForward, sclib::modeHypothesized, sclib::alignmentEnd); //...and the respective segment-end
	//this->pGT->getNextSegment(0, speechStart, speechEnd, sclib::atSpeechSegmentStart, sclib::searchForward, sclib::modeHypothesized); //get the boundaries for the first speech segment: segment-start...
	//this->pGT->getNextSegment(speechStart+1, tmp, speechEnd, sclib::atSpeechSegmentEnd, sclib::searchForward, sclib::modeHypothesized); //...and the respective segment-end
	start = 0; //this are the boundaries of the first nonspeech-segment, then: it just preceeds the first speech segment
	end = speechStart - 1;

	for (unsigned long int s = 0; s < this->pGT->getAudioSampleCount(); s++) {
		if ((long int)(s) > end) { //if the current speech- or nonspeech-segment is exhausted, look for the next one:
			if (end < speechEnd) { //we where in a nonspeech-segment, so we already now the boundaries of the next speech segment
				start = speechStart;
				end = speechEnd;
			} else { //the current speech segment is over, so get the boundaries of the next one but first account for the noise-segment in between
				start = speechEnd + 1;
				speechStart = this->pGT->getNextOccasion(s, sclib::atSpeechSegmentStart, sclib::searchForward, sclib::modeHypothesized, sclib::alignmentStart);
				//this->pGT->getNextSegment(s, speechStart, speechEnd, sclib::atSpeechSegmentStart, sclib::searchForward, sclib::modeHypothesized);
				if (speechStart == sclib::noSegment) {
					end = this->pGT->getAudioSampleCount() - 1;
				} else {
					speechEnd = this->pGT->getNextOccasion(speechStart+1, sclib::atSpeechSegmentEnd, sclib::searchForward, sclib::modeHypothesized, sclib::alignmentEnd);
					//this->pGT->getNextSegment(speechStart+1, tmp, speechEnd, sclib::atSpeechSegmentEnd, sclib::searchForward, sclib::modeHypothesized);
					end = speechStart - 1;
				}
			}
		}
		
		if (sclib::isBetween(start+softBoundaryDiameter, s, end-softBoundaryDiameter) == true) { //account for soft boundaries by not counting any wrong decisions made in the diameter around noise/speech segment boundaries
			gtId = this->pGT->getSamplesSpeakerID(s, sclib::modeGroundtruth);
			hypoId = this->pGT->getSamplesSpeakerID(s, sclib::modeHypothesized);
			hypoGtId = this->pGT->getSpeakerGIDFromHID(hypoId, correct);

			if ((hypoId!=sclib::noSpeaker && gtId==sclib::noSpeaker) ||
					(hypoId==sclib::noSpeaker && gtId!=sclib::noSpeaker) ||
					(hypoGtId!=gtId) ||
					(hypoGtId!=sclib::noSpeaker && hypoGtId==gtId && correct == false)) {
				wrong++;
			}
		}
	}

	return (double)(wrong) / (double)(this->pGT->getAudioSampleCount());*/
}
