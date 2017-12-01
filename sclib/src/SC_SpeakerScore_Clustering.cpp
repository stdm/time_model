/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the created       */
/*        clustering according to the ground truth                  			*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 13.09.2005																								*/
/**************************************************************************/

#include <assert.h>
#include "SC_SpeakerScore_Clustering.h"
#include "SC_Aux.h"
#include "SC_MatrixFunctions.h"
#include "SC_Timer.h"

//====================================================================================================================
//	The constructor
//  the final partition must be a pointer to an onject within the partition-list in order to get destructed!
//====================================================================================================================
SC_SpeakerScore_Clustering::SC_SpeakerScore_Clustering(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT, SC_Partition *pPartitionList, SC_Partition *pFinalPartition, bool computeMergeHistory) : SC_SpeakerScore(pTweak, pGT) {
  this->pPartitionList = pPartitionList; 
  this->pFinalPartition = pFinalPartition;
  this->pClassMapping = NULL;
	this->computeMergeHistory = computeMergeHistory;

	this->ca = pGT->getRealSpeakerCount(true);
	this->cu = (this->pFinalPartition != NULL) ? sclib::getListCount(this->pFinalPartition->getClusters()) : 0;

	for (unsigned int ce = 0; ce < 4; ce++) {
		this->scatterMatrix[ce] = NULL;
		this->recall[ce] = NULL;
		this->precision[ce] = NULL;
		this->missRate[ce] = NULL;
		this->falseAlarmRate[ce] = NULL;
		this->errorRate[ce] = NULL;
		this->specificity[ce] = NULL;
		this->accuracy[ce] = NULL;
		this->overallRecall[ce] = 0.0;
		this->overallPrecision[ce] = 0.0;
		this->missclassificationRate[ce] = 0.0;
		this->averageClusterPurity[ce] = 0.0;
		this->averageClassPurity[ce] = 0.0;
		this->purity[ce] = 0.0;
		this->randIndex[ce] = 0.0;
		this->BBNmetric[ce] = 0.0; 
		this->classPurity[ce] = NULL;
		this->clusterPurity[ce] = NULL;
	}
	this->mergeHistory = NULL;
	this->averagePrecision = NULL;
	this->average2precision = 0.0;
	this->mergeHistoryIsValid = true;
	this->DER = 0.0;

	sprintf(this->purpose, "%s\0", "Speaker-Clustering Assessment");
}

//====================================================================================================================
//	The destructor
//  destructs all the linked partitions, too!!!
//====================================================================================================================
SC_SpeakerScore_Clustering::~SC_SpeakerScore_Clustering() {
  sclib::destructLinkedList(this->pPartitionList);

	for (unsigned int ce = 0; ce < 4; ce++) {
		MFree_2D(this->scatterMatrix[ce]);
		MFree_1D(this->recall[ce]);
		MFree_1D(this->precision[ce]);
		MFree_1D(this->missRate[ce]);
		MFree_1D(this->falseAlarmRate[ce]);
		MFree_1D(this->errorRate[ce]);
		MFree_1D(this->specificity[ce]);
		MFree_1D(this->accuracy[ce]);
		MFree_1D(this->classPurity[ce]);
		MFree_1D(this->clusterPurity[ce]);
	}

	MFree_2D(this->mergeHistory);
	MFree_1D(this->averagePrecision);
}

//====================================================================================================================
//	Convert a cluster-label (i.e. hypo-speaker-id's) to a valid row-index into the scatter-matrix; SVLIB_Fail is 
//  returned if the mapping can't be established
//====================================================================================================================
long int SC_SpeakerScore_Clustering::cluster2idx(unsigned long int clusterTag, int countedEntity) {
	long int res = SVLIB_Fail;
	unsigned long int count = 0;
  SC_SpeakerScore::SC_IndexMapping *pMapping, *pHook = this->pClusterMapping[countedEntity], *pLast = NULL;

  //find the index
  while (pHook != NULL) {
    if (pHook->idx1 == clusterTag) {
      res = pHook->idx2;
			break;
    }
    count++;
    pLast = pHook;
    pHook = pHook->Next;
  }

  //if not found, make a new entry and return the new index
  if (res == SVLIB_Fail && count < this->cu) { //don't create more mappings than there are clusters
    pMapping = new SC_SpeakerScore::SC_IndexMapping(clusterTag, count);
    if (pLast == NULL) {
      this->pClusterMapping[countedEntity] = pMapping;
    } else {
      pLast->Next = pMapping;
    }
    res = pMapping->idx2;
  }
  
	return res;
}

//====================================================================================================================
//	Convert a valid row-index of the scatter-matrix into a hypo-speaker-id; SVLIB_Fail is returned on failure (unknown
//  speaker-id); return the speaker-id as either shortName or clusterName
//====================================================================================================================
long int SC_SpeakerScore_Clustering::idx2cluster(unsigned long int clusterIdx, int countedEntity, char* clusterName, bool shortName) {
	long int res = SVLIB_Fail;
	SC_SpeakerScore::SC_IndexMapping *pHook = this->pClusterMapping[countedEntity];

	while (pHook != NULL) {
		if (pHook->idx2 == clusterIdx) {
			res = pHook->idx1;
			if (clusterName != NULL) {
				if (shortName == true) {
					sprintf(clusterName, "%d\0", pHook->idx1); //take the speaker-id as a short-name
				} else {
					sprintf(clusterName, "%d\0", pHook->idx1); //take the speaker-id also as a cluster-name
				}
			}
			break;
		}
		pHook = pHook->Next;
	}

	return res;
}

//====================================================================================================================
//	Exchanges for two arbitrary cluster-tags the indexes to which they are mapped
//====================================================================================================================
void SC_SpeakerScore_Clustering::exchangeClusterMappings(unsigned long int clusterIdx1, unsigned long int clusterIdx2, int countedEntity) {
	SC_SpeakerScore::SC_IndexMapping *pHook = this->pClusterMapping[countedEntity];
	unsigned long int exchangedClusterID;

	if (clusterIdx1 != clusterIdx2) {
		while (pHook != NULL) {
			if (pHook->idx2 == clusterIdx1) {
				pHook->idx2 = clusterIdx2;
				exchangedClusterID = pHook->idx1;
				break;
			}
			pHook = pHook->Next;
		}

		pHook = this->pClusterMapping[countedEntity];
		while (pHook != NULL) {
			if (pHook->idx2==clusterIdx2 && pHook->idx1!=exchangedClusterID) { //test both indices because clusterIdx2 now exists 2 times in the list
				pHook->idx2 = clusterIdx1;
				break;
			}
			pHook = pHook->Next;
		}
	}

	return;
}

//====================================================================================================================
//	if a global criterion has been computed during clustering, but hasn't been used as the termination criterion (i.e. 
//  sclib::tcTrue has been used), we have a final partition with only one cluster but due to the global criterion this
//  class thinks that "cu" distinct speakers have been found; this method picks and returns the partition from the
//  list that has "cu" clusters and thus is the wanted "virtual" final partition.
//====================================================================================================================
SC_Partition* SC_SpeakerScore_Clustering::findVirtualFinalPartition(SC_Partition *pPartitionList, SC_Partition *pFinalPartition, unsigned int virtualSpeakerCount) {
	SC_Partition *pVirtual = pFinalPartition;
	unsigned int partitionCount, finalSpeakerCount = sclib::getListCount(pFinalPartition->getClusters());

	if (finalSpeakerCount < virtualSpeakerCount) {
		partitionCount = sclib::getListCount(pPartitionList);
		pVirtual = sclib::getListWithIndex(pPartitionList, partitionCount-(virtualSpeakerCount-finalSpeakerCount)-1);
	}

	return pVirtual;
}

//====================================================================================================================
//	Creates (recursively) the merge-history for the given clusterID and stores the findings in the given resultVector
//====================================================================================================================
void SC_SpeakerScore_Clustering::createMergeHistory(long int clusterID, SC_Partition *pPartitionList, SC_Partition *pCurrentPartition, long int finalOwnerID, long int* resultVector) {
	SC_Cluster *pMotherHook, *pFatherHook, *pChildHook = pCurrentPartition->getClusterByID(clusterID);
	long int childID = clusterID, motherID = pChildHook->getParentID(0), fatherID = pChildHook->getParentID(1);
	long int motherOwnerID, fatherOwnerID;
	long int searchPartitionIdx = sclib::getListItemIndex(pPartitionList, pCurrentPartition) - 1; //index of the predecessor of the currentPartition in the list of successive partitions (starting with 0 for the first partition)
	unsigned long int resultPosition;
	SC_Partition *pSearchPartition;

	//search for the parents; therefore, step back in the partition-history to the partition just before the mergence happened and both parents still existed separately
	if (searchPartitionIdx >= 0 && motherID != sclib::noSpeaker && fatherID != sclib::noSpeaker) {
		do {
			pSearchPartition = sclib::getListWithIndex(pPartitionList, searchPartitionIdx--);
			pMotherHook = pSearchPartition->getClusterByID(motherID);
			pFatherHook = pSearchPartition->getClusterByID(fatherID);
		} while ((pMotherHook == NULL || pFatherHook == NULL) && searchPartitionIdx >= 0);
	} else {
		//the child is already an atomic cluster and the currentPartition is the first one, so no more mergences are to be evaluated
		return;
	}

	//now we can get all info about the parents (ownerIDs we need...) from that searchPartition; we also save it if this information is not already stored in the cluster object
	if (pMotherHook->getOwnerID() == sclib::noSpeaker) {
		motherOwnerID = determineClusterOwner(pMotherHook, 0, sclib::ceSample, NULL, sclib::maxSpeakers);
		pMotherHook->setOwnerID(motherOwnerID);
	} else {
		motherOwnerID = pMotherHook->getOwnerID();
	}
	if (pFatherHook->getOwnerID() == sclib::noSpeaker) {
		fatherOwnerID = determineClusterOwner(pFatherHook, 0, sclib::ceSample, NULL, sclib::maxSpeakers);
		pFatherHook->setOwnerID(fatherOwnerID);
	} else {
		fatherOwnerID = pFatherHook->getOwnerID();
	}

	//a mergence is considered correct iif both parents have the same owner and also the same as the final cluster
	resultPosition = pChildHook->getMergenceCount() - 1; //segmentcount can't be used to guess nr. of mergences because of uniteNeighboringSegments(), so we introduced this direct measure; -1 to make it a 0-based index
	resultVector[resultPosition] = (motherOwnerID == fatherOwnerID && motherOwnerID == finalOwnerID) ? 1 : 0;

	//recursively create merge-history-entries for the parent-clusters if we didn't already reach the bottom of the list
	if (pSearchPartition != pPartitionList) {
		createMergeHistory(motherID, pPartitionList, pSearchPartition, finalOwnerID, resultVector);
		createMergeHistory(fatherID, pPartitionList, pSearchPartition, finalOwnerID, resultVector);
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
void SC_SpeakerScore_Clustering::calcScores(unsigned long int start, unsigned long int end, unsigned long int algorithmicUncertaintyDiameter) {
	unsigned long int **CEsFromOwner[4], maxCEs, minIdx, mergeHistoryVectorLength = 0;
	unsigned long int x, y, z, scene = 0, lastScene = 0, shot = 0, lastShot = 0, **orderedSegments = NULL, softBoundaryDiameter = this->pGT->getUncertaintyRegionWidth(true) + algorithmicUncertaintyDiameter;
	long int hypoSegStart, hypoSegEnd, gtSegStart = 0, gtSegEnd, oldGtSegStart = sclib::noSegment, gtID = sclib::noSpeaker, hypoID, idxX, idxY, lastGtID = sclib::noSpeaker;
	bool **hitsInScene, **hitsInShot, correct, hasSpeakerBoundary;
	int ce;
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	SC_Cluster *pHook;
	SC_Partition *pVirtualFinalPartition;
	long int subSegStart, subSegEnd;
  
	if (this->pFinalPartition == NULL) { //we urgently need this for most of the calculations, so... aborted clustering => no scoring
		return;
	}

  //initializations...
	this->start = start;
	this->end = (end > 0 && end < this->pGT->getAudioSampleCount()) ? end : this->pGT->getAudioSampleCount()-1; //not so very relevant in the clustering context...
	sclib::destructLinkedList(this->pClassMapping);
	this->ca = this->pGT->getRealSpeakerCount(true, this->start, this->end);
	this->cu = sclib::max(this->ca, sclib::getListCount(this->pFinalPartition->getClusters())); //we need minimum as much clusters as classes to fill the scatter-matrix that way that the upper square part contains on it's diagonal the non-errors
	if (this->cu <= 0 || this->ca <= 0) { //don't calc anything if there is no gt
		this->cu = 0;
		this->ca = 0;
		MFree_0D(pFunc);
		return;
	}
	for (ce = 0; ce < 4; ce++) {
		sclib::destructLinkedList(this->pClusterMapping[ce]);
		MFree_2D(this->scatterMatrix[ce]);
		this->scatterMatrix[ce] = pFunc->initMatrix(this->cu+2, this->ca, (long int)(0));
		MArray_2D(CEsFromOwner[ce], (long int)(this->cu), 2, unsigned long int, "SC_SpeakerScore_Clustering.calcScores: CEsFromOwner"); //format per row: |ownerID=gtSpeakerId|#CEs from this speaker| (idx of each row is class2idx(clusterID))
	}
	MArray_2D(hitsInScene, (long int)(this->cu)+2, this->ca, bool, "SC_SpeakerScore_Clustering.calcScores: hitsInScene");
	MArray_2D(hitsInShot, (long int)(this->cu)+2, this->ca, bool, "SC_SpeakerScore_Clustering.calcScores: hitsInShot");
	for (y = 0; y < this->cu+2; y++) {
		for (x = 0; x < this->ca; x++) {
			hitsInScene[y][x] = false;
			hitsInShot[y][x] = false;
		}
	}

	//remove all previous speaker-mappings in the ground-truth
  this->pGT->removeAllSpeakerMappings();

  //calc available gt segment-count/-size (pos/pop-set)
  for (y = 0; y < this->pGT->getAudioSampleCount(); y++) {
		this->pGT->getNextSegment(y, gtSegStart, gtSegEnd, sclib::atSpeech, sclib::searchWithin, sclib::modeGroundtruth, true, true, sclib::atShotBoundary); //SEARCH_WITHIN instead of SEARCH_FORWARD, so that if a boundary cuts a speech-segment in two pieces, they are both counted
		if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment) {

			//check for scene- and shot-changes
			scene = this->pGT->sample2scene(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastScene);
			if (scene > this->pTweak->general.lastScene) {
				break;
			}

			//has this scene been processed by any algorithm?
			if (scene >= this->pTweak->general.firstScene && (this->pTweak->general.sceneSelection == 0 || sclib::bitTest(sclib::bit(scene), pTweak->general.sceneSelection) == true) && scene <= this->pTweak->general.lastScene) {

				if (scene != lastScene) {
					pFunc->clear(hitsInScene, this->cu+2, this->ca, false); //clear rememberance of hits
				}
				shot = this->pGT->sample2shot(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastShot);
				if (shot != lastShot) {
					pFunc->clear(hitsInShot, this->cu+2, this->ca, false); //clear rememberance of hits
				}

				//do the real counting work, all around is just framework...
				if (pGT->getSurroundingSpeakersSegmentSize(gtSegStart, gtSegEnd) >= this->pGT->getConverter()->ms2sample(this->pTweak->speakerClusterer.speechSegLengthThreshold)) { //not the segment itself has to exceed the theshld, but the complete content inside the speaker-boundaries it lies in
					lastGtID = gtID;
					gtID = this->pGT->getSamplesSpeakerID(gtSegStart, sclib::modeGroundtruth); //speaker-id according to ground-truth
					hasSpeakerBoundary = (this->pGT->testSegment(gtSegStart, gtSegStart, true, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0);
					idxX = class2idx(gtID);

					this->scatterMatrix[sclib::ceSample][this->cu][idxX] += gtSegEnd - gtSegStart + 1; //Population
					this->scatterMatrix[sclib::ceSample][this->cu+1][idxX] += gtSegEnd - gtSegStart + 1; //Possible set
					if (gtID!=lastGtID || hasSpeakerBoundary==true) { //regard speaker-homogenious segments as segments here, not just speech separated by silence
						this->scatterMatrix[sclib::ceSegment][this->cu][idxX] += 1;
						this->scatterMatrix[sclib::ceSegment][this->cu+1][idxX] += 1;
					}
					if (hitsInShot[this->cu][idxX] == false) {
						this->scatterMatrix[sclib::ceShot][this->cu][idxX] += 1;
						hitsInShot[this->cu][idxX] = true;
					}
					if (hitsInShot[this->cu+1][idxX] == false) {
						this->scatterMatrix[sclib::ceShot][this->cu+1][idxX] += 1;
						hitsInShot[this->cu+1][idxX] = true;
					}
					if (hitsInScene[this->cu][idxX] == false) {
						this->scatterMatrix[sclib::ceScene][this->cu][idxX] += 1;
						hitsInScene[this->cu][idxX] = true;
					}
					if (hitsInScene[this->cu+1][idxX] == false) {
						this->scatterMatrix[sclib::ceScene][this->cu+1][idxX] += 1;
						hitsInScene[this->cu+1][idxX] = true;
					}
				}
			} //scene was to evaluate

			lastScene = scene;
			lastShot = shot;
      oldGtSegStart = gtSegStart;
			y = gtSegEnd;
		} else { //no valid segment
			break;
		}
	}

	//process each cluster, fill it's scatter-matrix-row (showing the spread of gt-speaker-id's in it's segments) and determine it's owner
	//this->pGT->frameListOut("frames.txt", 0, this->pGT->getAudioSampleCount()-1);
	pHook = this->pFinalPartition->getClusters();
	while (pHook != NULL) {
		pHook->uniteNeighboringSegments(); //to allow correct segment-based scoring
		orderedSegments = pHook->getOrderedSegments(); //we need an ordered list of segments for correct scene/shot-handling
		hypoID = pHook->getSpeakerID();

		oldGtSegStart = sclib::noSegment; //care for correct initialization of scene/shot-handling
		lastScene = 0;
		lastShot = 0;
		gtID = sclib::noSpeaker;

		for (y = 0; y < pHook->getSegmentCount(); y++) {
			hypoSegStart = orderedSegments[y][0];
			hypoSegEnd = orderedSegments[y][1];

			for (x = hypoSegStart; x <= (unsigned long int)(hypoSegEnd); x++) {
				//a problem with segment-based scoring is that initial speaker models are not build per speech-segment, but per speaker-segment 
				//(i.e. alle segments between [hypothesized] speaker boundaries); so, because a speaker-segment is ofton composed of several 
				//speech portions separated by silence, the old check below of equal borders between cluster segments (speaker-based) and
				//class segments (speech-based) most ofton return false and thus segment-based scores are awful. what we do now is dividing
				//a speaker-based (hypo) segment into its (hypothesized) speech-segment-parts and try to match those with the gt-segments
				//TODO: here we search for gt-segments within one hypo-sub-segment; what if it is the other way round and we need to search 
				//      for hypo-segments within gt-segments?				
				this->pGT->getNextSegment(x, subSegStart, subSegEnd, sclib::atSpeech, sclib::searchWithin,sclib::modeHypothesized, true, true, sclib::atShotBoundary);
				if (subSegStart != sclib::noSegment && subSegEnd != sclib::noSegment && subSegStart <= hypoSegEnd && subSegStart <= subSegEnd) {
					for (z = subSegStart; z <= (unsigned long int)(subSegEnd); z++) {
						this->pGT->getNextSegment(z, gtSegStart, gtSegEnd, sclib::atSpeech, sclib::searchWithin, sclib::modeGroundtruth, true, true, sclib::atShotBoundary);
						if (gtSegStart != sclib::noSegment && gtSegEnd != sclib::noSegment && gtSegStart <= hypoSegEnd && gtSegStart <= gtSegEnd) {
							lastGtID = gtID;
							gtID = this->pGT->getSamplesSpeakerID(gtSegStart, sclib::modeGroundtruth);
							hasSpeakerBoundary = (this->pGT->testSegment(gtSegStart, gtSegStart, true, sclib::atSpeakerBoundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0);
							if (gtID != sclib::noSpeaker) {
								idxX = class2idx(gtID);

								//care for soft boundaries
								if (abs(gtSegStart - subSegStart) < (long int)(softBoundaryDiameter)) {
									gtSegStart = subSegStart;
								}
								if (abs(gtSegEnd - subSegEnd) < (long int)(softBoundaryDiameter)) {
									gtSegEnd = subSegEnd;
								}

								//care for shot/scene-changes
								scene = this->pGT->sample2scene(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastScene);
								if (scene != lastScene) {
									pFunc->clear(hitsInScene, this->cu+2, this->ca, false); //clear rememberance of hits
								}
								shot = this->pGT->sample2shot(gtSegStart, (oldGtSegStart == sclib::noSegment) ? 0 : oldGtSegStart, lastShot);
								if (shot != lastShot) {
									pFunc->clear(hitsInShot, this->cu+2, this->ca, false); //clear rememberance of hits
								}

								//fill scatter-matrix
								idxY = cluster2idx(hypoID, sclib::ceSegment);
								if (gtID!=lastGtID || hasSpeakerBoundary==true) { //regard speaker-homogenious segments as segments here, not just speech separated by silence
									this->scatterMatrix[sclib::ceSegment][idxY][idxX] += (gtSegStart == subSegStart && gtSegEnd == subSegEnd) ? 1 : 0;
								}
								if (gtSegEnd > subSegEnd) { //truncate overlapping segments
									gtSegEnd = subSegEnd;
								}
								idxY = cluster2idx(hypoID, sclib::ceSample);
								this->scatterMatrix[sclib::ceSample][idxY][idxX] += gtSegEnd - gtSegStart + 1;
								idxY = cluster2idx(hypoID, sclib::ceShot);
								if (hitsInShot[idxY][idxX] == false) {
									this->scatterMatrix[sclib::ceShot][idxY][idxX] += 1;
									hitsInShot[idxY][idxX] = true;
								}
								idxY = cluster2idx(hypoID, sclib::ceScene);
								if (hitsInScene[idxY][idxX] == false) {
									this->scatterMatrix[sclib::ceScene][idxY][idxX] += 1;
									hitsInScene[idxY][idxX] = true;
								}
							} //no gt-speaker
						} else { //no more gt-segments within this cluster-sub-segment
							break;
						}

						oldGtSegStart = gtSegStart;
						lastScene = scene;
						lastShot = shot;
						z = gtSegEnd;
					} //for all gt-segments within a cluster-sub-segment (z)
				} else { //no more sub-segments within this cluster-segment
					break;
				} 

				x = subSegEnd;
			} //for all subsgements within a cluster-segment (x)
		} //for all cluster-segments (y)

		for (ce = 0; ce < 4; ce++) {
			idxY = cluster2idx(hypoID, ce); //row index (in scatter matrix) of current cluster
			idxX = pFunc->maxIdx(this->scatterMatrix[ce][idxY], this->ca); //column index with the most ce's
			gtID = idx2class(idxX); //gt-id corresponding with the found max. column
			CEsFromOwner[ce][idxY][0] = gtID;
			CEsFromOwner[ce][idxY][1] = this->scatterMatrix[ce][idxY][idxX]; //to facilitate finding the correct owner later
			if (ce == sclib::ceSample) { //the real owner is computed only sample-based!
				pHook->setOwnerID(gtID); //determine the owner of this cluster (=speaker that utterred most of the included samples according to gt)
				this->pGT->addSpeakerMapping(gtID, hypoID, false); //add a speaker-mapping in the groundtruth, but make no statement about correctness of the ownership yet
			}
		}

		MFree_2D(orderedSegments); //clean up & go on
		pHook = pHook->Next;
	}

	//determine the correctness of the ownerships of clusters and reorder rows in the upper square part of the scatter matrix
	for (ce = 0; ce < 4; ce++) {
		//sclib::matrixOut("CEsFromOwner.txt", CEsFromOwner[ce], this->cu, 2, this->pTweak);
		//sclib::matrixOut("scatterMatrix.txt", this->scatterMatrix[ce], this->cu+2, this->ca, this->pTweak);

		for (x = 0; x < this->ca; x++) { //for each possible speaker (not cluster!)
			maxCEs = 0;
			idxY = sclib::noSpeaker;
			gtID = idx2class(x);

			for (y = 0; y < this->cu; y++) { //step through the samplesFromOwner-list, find cluster with most samples from the current speaker => it is the correct owner
				if (CEsFromOwner[ce][y][0] == gtID && maxCEs < CEsFromOwner[ce][y][1]) {
					maxCEs = CEsFromOwner[ce][y][1];
					hypoID = this->pGT->getSpeakerHIDFromGID(gtID, correct);
					idxY = cluster2idx(hypoID, ce);
				}
			}
			if (idxY == sclib::noSpeaker || hypoID == sclib::noSpeaker) {
				maxCEs = 0;
			}	
			//sclib::quadrupelOut("idx.txt", x, gtID, this->pGT->getSpeakerHIDFromGID(gtID, correct), idxY, this->pTweak);
			//clusterMappingOut("mapping.txt", ce);

			if (maxCEs > 0 && ce == sclib::ceSample) {
				hypoID = idx2cluster(idxY, ce);
				if (hypoID != sclib::noSpeaker) {
					pHook = this->pFinalPartition->getClusterByID(hypoID);
					pHook->setOwnerIsCorrect(true);
					this->pGT->addSpeakerMapping(gtID, hypoID, true); 
				} else {
					maxCEs = 0;
				}
			}

			if (maxCEs > 0) { //we found a correct cluster for speaker 'x', so remember it in the cluster-object and the ground-truth' speaker-mapping if we are in sclib::ceSample, otherwise only exchange rows
				//exchange rows in the upper (correct) square of the matrix
				pFunc->exchangeRows(this->scatterMatrix[ce], this->cu+2, this->ca, idxY, x); //so that the max value of this row occurs at colIdx==rowIdx
				//sclib::matrixOut("scatterMatrix.txt", this->scatterMatrix[ce], this->cu+2, this->ca, this->pTweak);
				exchangeClusterMappings(idxY, x, ce); //also alter the mapping accordingly
			} else {
				if (pFunc->equals(this->scatterMatrix[sclib::ceSample][x], this->ca, (long int)(0)) == false) { //make the 'correct' row for this class ==0, if not already so
					for (z = x+1; z < this->cu; z++) { //find the next all-0-row to exchange with the current row
						if (pFunc->equals(this->scatterMatrix[sclib::ceSample][x], this->ca, (long int)(0)) == true) {
							pFunc->exchangeRows(this->scatterMatrix[ce], this->cu+2, this->ca, z, x);
							exchangeClusterMappings(z, x, ce);
						}
					}
				}
			}

		}//for each possible speaker

		MFree_2D(CEsFromOwner[ce]);
	}

	//reorder the below-square part of the scatter-matrix according to the desires of the calculation-methods of the base-class:
	//*first* rows after the square part contain the rows of just fitting clusters with the same owner as the *first* rows of the square part and so on...
	for (ce = 0; ce < 4; ce++) {
		for (x = this->ca; x < this->cu; x++) { //after each loop, the x'th row is properly filled
			//sclib::matrixOut("scatterMatrix.txt", this->scatterMatrix[ce], this->cu+2, this->ca, this->pTweak);
			minIdx = this->cu + 1;
			maxCEs = 0;
			idxY = x;

			//find the minimum col-idx of all per-row-maxima; take as the next storable row the one with max. CEs in this min-idx col
			for (y = x; y < this->cu; y++) {
				idxX = pFunc->minIdx(this->scatterMatrix[ce][y], this->ca);
				if (idxX < (long int)(minIdx)) {
					minIdx = idxX;
					maxCEs = 0; //a smaller minIdx is found => reset the maxCount and directly refill it in the next block
				}
				if (idxX == (long int)(minIdx) && this->scatterMatrix[ce][y][idxX] > (long int)(maxCEs)) { //find the row with the max. nr of CEs at the minIdx col-position
					maxCEs = this->scatterMatrix[ce][y][idxX];
					idxY = y; //row-idx of the row>=x with the biggest nr auf CEs located in a column at leftmost position
				}
			}

			pFunc->exchangeRows(this->scatterMatrix[ce], this->cu+2, this->ca, x, idxY);
			exchangeClusterMappings(x, idxY, ce);
		}
	}

	//calculate the merge-history
	//those mergences are considered correct which merge two clusters both belonging to the final (according to this->pFinalPartition) ownerID of a cluster
	//ATTENTION: This doesn't capture "mergences" done by the SC_SpeakerIdentificator class because they are not noted in a succession-of-partitions-way!
	MFree_2D(this->mergeHistory);
	if (this->computeMergeHistory == true) {
		mergeHistoryVectorLength = sclib::getListCount(this->pPartitionList) + 2;
		this->mergeHistory = pFunc->initMatrix(this->cu, mergeHistoryVectorLength, (long int)(-2));
		pVirtualFinalPartition = findVirtualFinalPartition(this->pPartitionList, this->pFinalPartition, this->cu);
		pHook = pVirtualFinalPartition->getClusters();
		//pHook = this->pFinalPartition->getClusters(); //can't always be used, see comment on findVirtualFinalPartition()!
		this->mergeHistoryIsValid = true;
		while (pHook != NULL) { //analyze the merge-history of each cluster separately
			hypoID = pHook->getSpeakerID();
			idxY = cluster2idx(hypoID);

			if (idxY != sclib::noSpeaker) {
				createMergeHistory(hypoID, this->pPartitionList, pVirtualFinalPartition, pHook->getOwnerID(), this->mergeHistory[idxY]);
				
				idxX = pHook->getSegmentCount() - 2; //there is a mergence less than segments in the final cluster; minus one more to convert the number to an index;
				this->mergeHistory[idxY][idxX+1] = -1; //terminate "list" after the last mergence-history entry
				this->mergeHistory[idxY][idxX+2] = getNrOfSegmentsWithOwner(this->pFinalPartition->getClusters(), pHook->getOwnerID()); //total number of available segments of this owner

				//check if merge-hisory is valid
				for (idxX = 0; idxX < (long int)(mergeHistoryVectorLength); idxX++) {
					if (this->mergeHistory[idxY][idxX] == -1) {
						break;
					}
					if (this->mergeHistory[idxY][idxX] == -2) { //an uninitialized value before the list-termination means that there where more segments in the cluster than ancestors in the previous partitions => history not complete => not valid
						this->mergeHistoryIsValid = false;
						break;
					}
				}
			} else {
				this->mergeHistoryIsValid = false;
			}

			if (this->mergeHistoryIsValid != true) {
				MFree_2D(this->mergeHistory);
				break;
			}

			pHook = pHook->Next;
		}
	} else { 
		this->mergeHistoryIsValid = false;
	} //computeMergeHistory == false


	//really calculate scores now
	for (ce = 0; ce < 4; ce++) {
		calcUnsupervisedClassificationScores(this->scatterMatrix[ce], this->cu, this->ca, this->mergeHistory, mergeHistoryVectorLength, this->overallRecall[ce], this->overallPrecision[ce], this->missclassificationRate[ce], this->averageClusterPurity[ce], this->averageClassPurity[ce], this->average2precision, this->purity[ce], this->randIndex[ce], this->BBNmetric[ce], this->averagePrecision, this->classPurity[ce], this->clusterPurity[ce], this->recall[ce], this->precision[ce], this->missRate[ce], this->falseAlarmRate[ce], this->errorRate[ce], this->specificity[ce], this->accuracy[ce]);
	}
	this->DER = calcDER(softBoundaryDiameter); //uses the groundtruth object instead of the scatter matrix to calculate the score because non-speech segments are not included in the scatter matrixes but are scored in DER

	MFree_0D(pFunc);
	MFree_2D(hitsInScene);
	MFree_2D(hitsInShot);

  return;
}

//====================================================================================================================
//	Output
//====================================================================================================================
ostream& SC_SpeakerScore_Clustering::output(ostream& OutS) {
	unsigned long int x, y, ce;
	char buf[sclib::bufferSize], buf2[sclib::bufferSize];
	bool printedSth = false;
	SC_Cluster* pHook;
	long int sceneStart, sceneEnd;
	unsigned int sceneNr = 1;

	if (this->pFinalPartition == NULL || this->scatterMatrix[sclib::ceSample] == NULL) {
		OutS << "\nNo scores available due to a missing final partition... sorry!\n";
	} else {
		//print out the speaker-names and id's
		OutS << "\n";
		OutS << "IDs and names of hypothesized speakers:\n";
		OutS << "----------------------------------------\n";
		pHook = this->pFinalPartition->getClusters();
		while (pHook != NULL) {
			OutS << pHook->getSpeakerID() << ": " << pHook->getSpeakerName() << "\n"; //" [owner correct: " << pHook->getOwnerIsCorrect() << "]\n";
			pHook = pHook->Next;
		}

		//overall classification
		OutS << "\n";
		OutS << "Clustering result:\n";
		OutS << "-------------------\n";
		OutS << "Scatter-matrix:\n";
		for (ce = 0; ce < 4; ce++) { //loop over all CEs
			OutS << "\t";
			for (x = 0; x < this->ca; x++) { //headline
				idx2class(x, buf, true); //get the class' short name (gt-speaker-id)
				OutS << "\t" << buf;
			}
			ce2text(ce, buf);
			OutS << " (per " << buf << ")\n";
			for (y = 0; y < this->cu; y++) { //table
				idx2cluster(y, ce, buf, true); //get the cluster's short name (hypo-speaker-id)
				OutS << "\t" << buf << "\t";
				for (x = 0; x < this->ca; x++) {
					OutS << this->scatterMatrix[ce][y][x] << "\t";
				}
				OutS << "\n";
			}
			OutS << "\n";
		}
		OutS << "                           Sample-\tSeg.-\tShot-\tScene-based\n";
		OutS << "Overall recall:            " << setprecision(4) << this->overallRecall[sclib::ceSample] << "\t" << this->overallRecall[sclib::ceSegment] << "\t" << this->overallRecall[sclib::ceShot] << "\t" << this->overallRecall[sclib::ceScene] << "\n";
		OutS << "Overall precision:         " << setprecision(4) << this->overallPrecision[sclib::ceSample] << "\t" << this->overallPrecision[sclib::ceSegment] << "\t" << this->overallPrecision[sclib::ceShot] << "\t" << this->overallPrecision[sclib::ceScene] << "\n";
		OutS << "Missclassification rate:   " << setprecision(4) << this->missclassificationRate[sclib::ceSample] << "\t" << this->missclassificationRate[sclib::ceSegment] << "\t" << this->missclassificationRate[sclib::ceShot] << "\t" << this->missclassificationRate[sclib::ceScene] << "\n";
		OutS << "Average cluster purity:    " << setprecision(4) << this->averageClusterPurity[sclib::ceSample] << "\t" << this->averageClusterPurity[sclib::ceSegment] << "\t" << this->averageClusterPurity[sclib::ceShot] << "\t" << this->averageClusterPurity[sclib::ceScene] << "\n";
		OutS << "Average class purity:      " << setprecision(4) << this->averageClassPurity[sclib::ceSample] << "\t" << this->averageClassPurity[sclib::ceSegment] << "\t" << this->averageClassPurity[sclib::ceShot] << "\t" << this->averageClassPurity[sclib::ceScene] << "\n";
		OutS << "Purity:                    " << setprecision(4) << this->purity[sclib::ceSample] << "\t" << this->purity[sclib::ceSegment] << "\t" << this->purity[sclib::ceShot] << "\t" << this->purity[sclib::ceScene] << "\n";
		OutS << "Rand index:                " << setprecision(4) << this->randIndex[sclib::ceSample] << "\t" << this->randIndex[sclib::ceSegment] << "\t" << this->randIndex[sclib::ceShot] << "\t" << this->randIndex[sclib::ceScene] << "\n";
		OutS << "BBN metric:                " << setprecision(4) << this->BBNmetric[sclib::ceSample] << "\t" << this->BBNmetric[sclib::ceSegment] << "\t" << this->BBNmetric[sclib::ceShot] << "\t" << this->BBNmetric[sclib::ceScene] << "\n";
		if (this->mergeHistoryIsValid == true) {
			OutS << "Average average precision: " << setprecision(4) << this->average2precision << "\n";
		} else {
			OutS << "Average average precision: " << "not available" << "\n";
		}
		OutS << "Diarization Error Rate:    " << setprecision(4) << this->DER << "\n";

		//detection of gt-speakers (vs. rest)
		OutS << "\n";
		OutS << "Class (GT-based speaker) detection results:\n";
		OutS << "--------------------------------------------\n";
		OutS << "                         Sample-\tSeg.-\tShot-\tScene-based\n";
		for (x = 0; x < this->ca; x++) { //loop over all classes
			this->idx2class(x, buf, false); //get the class' long name
			this->idx2class(x, buf2, true); //get the class' short name (gt-speaker-id)
			OutS << "-> " << buf << " (" << buf2 << ") vs. rest:\n";
			OutS << "    Recall:              " << setprecision(4) << this->recall[sclib::ceSample][x] << "\t" << this->recall[sclib::ceSegment][x] << "\t" << this->recall[sclib::ceShot][x] << "\t" << this->recall[sclib::ceScene][x] << "\n";
			OutS << "    Precision:           " << setprecision(4) << this->precision[sclib::ceSample][x] << "\t" << this->precision[sclib::ceSegment][x] << "\t" << this->precision[sclib::ceShot][x] << "\t" << this->precision[sclib::ceScene][x] << "\n";
			OutS << "    Specificity:         " << setprecision(4) << this->specificity[sclib::ceSample][x] << "\t" << this->specificity[sclib::ceSegment][x] << "\t" << this->specificity[sclib::ceShot][x] << "\t" << this->specificity[sclib::ceScene][x] << "\n";
			OutS << "    Omission:            " << setprecision(4) << 1.0-this->recall[sclib::ceSample][x] << "\t" << 1.0-this->recall[sclib::ceSegment][x] << "\t" << 1.0-this->recall[sclib::ceShot][x] << "\t" << 1.0-this->recall[sclib::ceScene][x] << "\n";
			OutS << "    Commission:          " << setprecision(4) << 1.0-this->precision[sclib::ceSample][x] << "\t" << 1.0-this->precision[sclib::ceSegment][x] << "\t" << 1.0-this->precision[sclib::ceShot][x] << "\t" << 1.0-this->precision[sclib::ceScene][x] << "\n";
			OutS << "    Miss rate:           " << setprecision(4) << this->missRate[sclib::ceSample][x] << "\t" << this->missRate[sclib::ceSegment][x] << "\t" << this->missRate[sclib::ceShot][x] << "\t" << this->missRate[sclib::ceScene][x] << "\n";
			OutS << "    False alarm rate:    " << setprecision(4) << this->falseAlarmRate[sclib::ceSample][x] << "\t" << this->falseAlarmRate[sclib::ceSegment][x] << "\t" << this->falseAlarmRate[sclib::ceShot][x] << "\t" << this->falseAlarmRate[sclib::ceScene][x] << "\n";
			OutS << "    Error rate:          " << setprecision(4) << this->errorRate[sclib::ceSample][x] << "\t" << this->errorRate[sclib::ceSegment][x] << "\t" << this->errorRate[sclib::ceShot][x] << "\t" << this->errorRate[sclib::ceScene][x] << "\n";
			OutS << "    Accuracy:            " << setprecision(4) << this->accuracy[sclib::ceSample][x] << "\t" << this->accuracy[sclib::ceSegment][x] << "\t" << this->accuracy[sclib::ceShot][x] << "\t" << this->accuracy[sclib::ceScene][x] << "\n";
			OutS << "    Fidelity:            " << setprecision(4) << 1.0-this->errorRate[sclib::ceSample][x] << "\t" << 1.0-this->errorRate[sclib::ceSegment][x] << "\t" << 1.0-this->errorRate[sclib::ceShot][x] << "\t" << 1.0-this->errorRate[sclib::ceScene][x] << "\n";
			OutS << "    Class purity:        " << setprecision(4) << this->classPurity[sclib::ceSample][x] << "\t" << this->classPurity[sclib::ceSegment][x] << "\t" << this->classPurity[sclib::ceShot][x] << "\t" << this->classPurity[sclib::ceScene][x] << "\n";
			OutS << "    Sensitivity:         see recall\n";
			OutS << "    Producer's accuracy: see omission\n";
			OutS << "    User's accuracy:     see commission\n";
			OutS << "    True positive rate:  see recall\n";
			OutS << "    False negative rate: see miss rate\n";
			OutS << "    False positive rate: see false alarm rate\n";
		}

		//measures per cluster
		OutS << "\n";
		OutS << "Results per cluster (hypothesized speaker):\n";
		OutS << "--------------------------------------------\n";
		OutS << "                           Sample-\tSeg.-\tShot-\tScene-based\n";
		for (x = 0; x < this->cu; x++) { //loop over all classes
			this->idx2class(x, buf, false); //get the cluster's long name
			this->idx2class(x, buf2, true); //get the cluster's short name (hypo-speaker-id)
			OutS << "-> " << buf << " (" << buf2 << "):\n";
			OutS << "    Cluster purity:      " << setprecision(4) << this->clusterPurity[sclib::ceSample][x] << "\t" << this->clusterPurity[sclib::ceSegment][x] << "\t" << this->clusterPurity[sclib::ceShot][x] << "\t" << this->clusterPurity[sclib::ceScene][x] << "\n";
			if (this->mergeHistoryIsValid == true) {
				OutS << "    Average precision:   " << setprecision(4) << this->averagePrecision[x] << "\n";
				OutS << "    Merge-history:       ";
				y = 0;
				while (this->mergeHistory[x][y] != -1) {
					OutS << ((this->mergeHistory[x][y] == 1) ? "ok" : "." ) << ((this->mergeHistory[x][y+1] != -1) ? " -> " : "");
					y++;
				}
				if (y == 0) { //if there's only one segment in the cluster, there were no mergences... it's ok, though
					OutS << "(ok)";
				}
				OutS << "\n";
			} else {
				OutS << "    Average precision:   " << "not available" << "\n";
				OutS << "    Merge-history:       " << "not available" << "\n";
			}
		}
	} //scores available

	//print out the segments with speaker-id's
  if (this->pTweak != NULL) {
		OutS << "\n";
		OutS << "Start- and endpoints of resulting speaker-segments:\n";
		OutS << "----------------------------------------------------\n";
    for (unsigned long int x = 0; x < this->pGT->getAudioSampleCount(); x++) {
		  this->pGT->getNextBoundary(x, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward, sclib::modeGroundtruth); 
      if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
        //is this scene to evaluate?
        if ((sceneNr <= this->pTweak->general.lastScene) && (sceneNr >= this->pTweak->general.firstScene) && ((sclib::bit(sceneNr) & this->pTweak->general.sceneSelection) || (this->pTweak->general.sceneSelection == 0))) {
					this->pGT->output(OutS, sceneStart, sceneEnd, sclib::atSpeechSegmentStart|sclib::atSpeechSegmentEnd, sclib::noType, sclib::modeHypothesized, printedSth);
          printedSth = true;
        }
        sceneNr++;
        x = sceneEnd;
        if (sceneNr > this->pTweak->general.lastScene) {
          break;
        }
      } else {
        break;
      }
    }
  }

  return OutS;
}

//====================================================================================================================
//	debug output of the internal mapping of hypo-ids to scatter-matrix rows
//====================================================================================================================
void SC_SpeakerScore_Clustering::clusterMappingOut(const char *fileName, unsigned int ce) {
	SC_SpeakerScore::SC_IndexMapping *pHook = this->pClusterMapping[ce];

	while (pHook != NULL) {
		sclib::tupelOut(fileName, pHook->idx1, pHook->idx2, this->pTweak);
		pHook = pHook->Next;
	}
	sclib::stringOut(fileName, "", this->pTweak);

	return;
}
