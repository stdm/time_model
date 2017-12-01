/**************************************************************************/
/*    Responsibility:																											*/
/*      - Performs the speaker clustering process                         */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#include <list>
#include <vector>
#include "SC_SpeakerClusterer.h"
#include "SC_ModelHandler.h"
#include "SC_Model_SVM.h"
#include "SC_Model_Time.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_SpeakerClusterer::SC_SpeakerClusterer(SC_TweakableParameters* pTweak, bool verbose){
  this->pTweak = pTweak;
  this->verbose = verbose;
	this->pDist = new SC_DistanceMeasures(this->pTweak, NULL, this->verbose);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_SpeakerClusterer::~SC_SpeakerClusterer(){
	MFree_0D(this->pDist);
}

//====================================================================================================================
//	Caculates the global criterion evaluating the partition according to the termination criterion given in the
//  tweakable parameters
//  Parameters: - pPartition is the pointer to the first element of the linked list of clusters forming the partition
//              - pParent1 and pParent2 are the last merged clusters
//              - pChild is the result of this mergence
//              - the numerical result of the algorithm is returned
//====================================================================================================================
double SC_SpeakerClusterer::getClusterGlobalCriterion(SC_Partition *pPartition, SC_Cluster *pParent1, SC_Cluster *pParent2, SC_Cluster *pChild) {
	double res;

	switch (this->pTweak->speakerClusterer.globalCriterion) {
		case sclib::gcBIC:
			res = this->pDist->BIC(pParent1, pParent2, this->pTweak->distanceMeasure.BICpenaltyFactor);
			break;
    case sclib::gcWCD:
      res = this->pDist->withinClusterDispersion(pPartition);
      break;
    case sclib::gcNone:
			res = 0.0;
			break;
		case sclib::gcICR:
			res = this->pDist->ICR(pParent1, pParent2);
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getClusterGlobalCriterion: Specified global criterion is unkonwn!");
      res = 0.0;
			break;
	}

  return res;
}

//====================================================================================================================
//	Decide if it is good to proceed cluster-mergence according to the criterion given in the tweakable parameters
//====================================================================================================================
bool SC_SpeakerClusterer::getClusterContinuation(double globalCriterion, unsigned int clusterCount, unsigned int speakerCount) {
  bool res;
  
  switch (this->pTweak->speakerClusterer.terminationCriterion) {
		case sclib::tcGc:
      switch (this->pTweak->speakerClusterer.globalCriterion) {
        case sclib::gcBIC:
          res = (globalCriterion >= 0.0) ? true : false;
          break;
        case sclib::gcWCD:
          res = true; //all partitions must be computed, then a final one will be picked
          break;
        case sclib::gcNone:
          res = true;
          break;
				case sclib::gcICR:
					res = true; //all partitions must be computed, then a final one will be picked
					break;
        default:
          REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getClusterContinuation: Specified global criterion is unkonwn!");
			    res = false;
          break;
      }
			break;
		case sclib::tcTrue:
			res = true;
      break;
		case sclib::tcFalse:
			res = false;
      break;
    case sclib::tcKnowledge:
      res = (clusterCount > speakerCount) ? true : false;
      break;
		case sclib::tcOptimal:
			res = true; //all partitions must be compiuted, then a final one will be picked
			break;
		case sclib::tcGcOptimal:
			res = true; //all partitions must be compiuted, then a final one will be picked
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getClusterContinuation: Specified termination-criterion is unkonwn!");
			res = false;
      break;
	}

  return res;
}

//====================================================================================================================
//	Pick and return a pointer to the partition form the linked list of all created partitions, which can be regarded 
//  as the final result according to the choosen global criterion
//====================================================================================================================
SC_Partition* SC_SpeakerClusterer::getClusterFinalPartition(SC_GroundTruth *pGT, SC_Partition *pFirstPartition) {
  SC_Partition *pRes = NULL, *pHook, *pPrevious;
  double current = 0.0, min = std::numeric_limits<double>::max(), *der;
	unsigned long int count, c;
	SC_SpeakerScore_Clustering *pScore = NULL;
	int finalNumber; //index of the final partition to pick

  switch (this->pTweak->speakerClusterer.terminationCriterion) {
		case sclib::tcGc:
      switch (this->pTweak->speakerClusterer.globalCriterion) {
        case sclib::gcBIC:
          pRes = sclib::getLastInList(pFirstPartition);
					finalNumber = sclib::getListCount(pFirstPartition)-1;
          break;

        case sclib::gcWCD: //find the partition with the minimum gc
          pPrevious = pFirstPartition;
					c = 0;
					if (pFirstPartition->Next != NULL) {
            pHook = pFirstPartition->Next;
          } else {
            pHook = pFirstPartition;
          }
          while (pHook != NULL) {
            current = pHook->getGlobalCriterion();
            if (current < min) {
              min = current;
              pRes = pPrevious; //because in clusterSpeakers() the gc is always attached to the newly-formed partition, we actually have to take the predecessor of that partition with minimum gc
							finalNumber = c;
            }
						pPrevious = pPrevious->Next;
						c++;
            pHook = pHook->Next;
          }
          break;

        case sclib::gcNone:
          pRes = sclib::getLastInList(pFirstPartition);
					finalNumber = sclib::getListCount(pFirstPartition)-1;
          break;

				case sclib::gcICR: //find (in reverse order, i.e. starting from the last partition with just 1 big cluster) the last partition where the gc doesn't exceed the ICR threshold
					count = sclib::getListCount(pFirstPartition);
					c = count-1;
					pHook = sclib::getLastInList(pFirstPartition);
					while (pHook->getGlobalCriterion() >= this->pTweak->distanceMeasure.ICRthreshold && c > 0) {
						pHook = sclib::getListWithIndex(pFirstPartition, --c); //get the previous partition
					}
					pRes = pHook;
					finalNumber = c;
					break;

        default:
          REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getClusterFinalPartition: Specified global criterion is unkonwn!");
			    pRes = sclib::getLastInList(pFirstPartition);
					finalNumber = sclib::getListCount(pFirstPartition)-1;
          break;
      }
			break;

		case sclib::tcTrue:
      pRes = sclib::getLastInList(pFirstPartition);
			finalNumber = sclib::getListCount(pFirstPartition)-1;
			break;

		case sclib::tcFalse:
      pRes = sclib::getLastInList(pFirstPartition);
			finalNumber = sclib::getListCount(pFirstPartition)-1;
      break;

    case sclib::tcKnowledge:
      pRes = sclib::getLastInList(pFirstPartition);
			finalNumber = sclib::getListCount(pFirstPartition)-1;
      break;

		case sclib::tcOptimal: //look for the partition (among those in the list) that yields minimum diarization error rate
			pScore = new SC_SpeakerScore_Clustering(this->pTweak, pGT, pFirstPartition, NULL, false);
			pHook = pFirstPartition;
			c = 0;
			while (pHook != NULL) {
				pHook->storeIDs(pGT, true, true, true);
				pScore->setFinalPartition(pHook);
				pScore->calcScores(0, 0, pHook->getUncertaintyDiameter());
				current = pScore->getDER();
				if (current < min) {
					min = current;
					pRes = pHook;
					finalNumber = c;
				}
				pHook = pHook->Next;
				c++;
			}
			pScore->setPartitionList(NULL); //so that no partition gets destructed
			MFree_0D(pScore);
			break;

		case sclib::tcGcOptimal:
			if (this->verbose ==true) { //this is a time-consuming task, so give the user some feedback...
				printf("\nScoring all partitions to optimize global criterion:       ");
			}
			pScore = new SC_SpeakerScore_Clustering(this->pTweak, pGT, pFirstPartition, NULL, false);
			count = sclib::getListCount(pFirstPartition);
			MArray_1D(der, count, double, "SC_SpeakerClusterer.getClusterFinalPartition: der");
			pHook = pFirstPartition;
			c = 0;
			while (pHook != NULL) { //first, cache the DER for all partitions
				pHook->storeIDs(pGT, true, true, true);
				pScore->setFinalPartition(pHook);
				pScore->calcScores(0, 0, pHook->getUncertaintyDiameter());
				der[c++] = pScore->getDER(); //remmber in the array, at idx position c, the DER for list element c
				if (this->verbose == true) {
					current = sclib::printPercentage(count, c, current, 0.0, c==1);
				}
				pHook = pHook->Next;
			}
			pScore->setPartitionList(NULL); //so that no partition gets destructed
			MFree_0D(pScore);
			if (this->verbose == true) {
				sclib::printPercentage(1, 1, 0.0, 0.0, false);
				printf("done, starting optimization\n");
			}
			//sclib::vectorOut("der.txt", der, c, false, this->pTweak);
			switch (this->pTweak->speakerClusterer.globalCriterion) { //then, optimize the penalities for the selected gc in order to get the partition with the minimum gc that is reachable
	      case sclib::gcBIC:
					//return predecessor of the first partition (starting from the beginning of the list) with gc < 0, due to this->pTweak->distanceMeasure.BICpenaltyFactor
					this->pTweak->distanceMeasure.BICpenaltyFactor = tuneBIC(pFirstPartition, der, c);
					pRes = sclib::getListWithIndex(pFirstPartition, c);
					finalNumber = c;
					break;

				case sclib::gcICR:
					//return the first partition (starting from the end of the list) with gc < this->pTweak->distanceMeasure.ICRthreshold, tune this threshold
					this->pTweak->distanceMeasure.ICRthreshold = tuneICR(pFirstPartition, der, c);
					pRes = sclib::getListWithIndex(pFirstPartition, c);
					finalNumber = c;
					break;

				case sclib::gcWCD:
					//return partition with minimum gc, due to this->pTweak->distanceMeasure.WCDpenaltyFactor
					this->pTweak->distanceMeasure.WCDpenaltyFactor = tuneWCD(pFirstPartition, der, c);
					pRes = sclib::getListWithIndex(pFirstPartition, c);
					finalNumber = c;
					break;

				case sclib::gcNone:
          pRes = sclib::getLastInList(pFirstPartition);
					finalNumber = sclib::getListCount(pFirstPartition)-1;
					break;

				default:
          REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getClusterFinalPartition: Specified global criterion is unkonwn!");
			    pRes = sclib::getLastInList(pFirstPartition);
					finalNumber = sclib::getListCount(pFirstPartition)-1;
					break;
			}
			MFree_1D(der);
			break;

		default:
			REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getClusterContinuation: Specified termination-criterion is unkonwn!");
      pRes = sclib::getLastInList(pFirstPartition);
			finalNumber = sclib::getListCount(pFirstPartition)-1;
      break;
  }

	//include the user's whish to get more or less clusters
	if (this->pTweak->speakerClusterer.specificity != 0.5) {
		if (this->pTweak->speakerClusterer.specificity > 1.0) {
			finalNumber = 0;
			pHook = pFirstPartition;
			while (pHook!=NULL) {
				if (sclib::getListCount(pHook->getClusters()) < (unsigned long int)(this->pTweak->speakerClusterer.specificity)) {
					finalNumber--;
					break;
				}
				finalNumber++;
				pHook = pHook->Next;
			}
		} else {
			int maxNumber = sclib::getListCount(pFirstPartition)-1;
			double score = sclib::invSigmoid(this->pTweak->speakerClusterer.specificity, 1.0, (double)(finalNumber), (double)(maxNumber), 0.000001, 0.0);
			finalNumber = sclib::getBetween(0, (this->pTweak->speakerClusterer.specificity>0.5) ? sclib::round(score+0.5) : sclib::round(score-0.5), maxNumber); //round up or down depending on the specificity value
		}
		pRes = sclib::getListWithIndex(pFirstPartition, finalNumber);
	}

  return pRes;
}

//====================================================================================================================
//	Returns the distance between the two clusters, based on the precomputed distances between the single segments in 
//  the clusters, stored in the distance matrix according to the choosen linkage criterion (doesn't apply to 
//  linkageMerge though, this is computed directly in getClusterDistanceMatrix(); the segmentMap maps segments 
//  (segmentStart) to rows/cols of the distance matrix
//====================================================================================================================
double SC_SpeakerClusterer::getClusterDistance(SC_Cluster *pFirst, SC_Cluster *pSecond, double **segmentDistanceMatrix, std::map<unsigned long int, unsigned long int> segmentMap) {
	double res = 0.0;
	double min = std::numeric_limits<double>::max(), max = -1.0*std::numeric_limits<double>::max();
	unsigned long int x, y, t, count = 0;

	switch (this->pTweak->speakerClusterer.linkageMode) {
		case sclib::linkageSingle:
			for (unsigned long int i = 0; i < pFirst->getSegmentCount(); i++) {
				for (unsigned long int j = 0; j < pSecond->getSegmentCount(); j++) {
					x = segmentMap[pFirst->getSegmentList(i, 1)];
					y = segmentMap[pSecond->getSegmentList(j, 1)];
					if (x > y) { //exchange indexes if pFirst has greater index than pSecond (because we only compute upper part below)
						t = x;
						x = y;
						y = t;
					}
					if (segmentDistanceMatrix[x][y] < min) {
						min = segmentDistanceMatrix[x][y]; //find the minimum pairwise distance between all the pairs of one segment of the first and one segment of the second cluster
					}
				}
			}
			res = min;
			break;

		case sclib::linkageComplete:
			for (unsigned long int i = 0; i < pFirst->getSegmentCount(); i++) {
				for (unsigned long int j = 0; j < pSecond->getSegmentCount(); j++) {
					x = segmentMap[pFirst->getSegmentList(i, 1)];
					y = segmentMap[pSecond->getSegmentList(j, 1)];
					if (x > y) { //exchange indexes if pFirst has greater index than pSecond (because we only compute upper part below)
						t = x;
						x = y;
						y = t;
					}
					if (segmentDistanceMatrix[x][y] > max) {
						max = segmentDistanceMatrix[x][y]; //find the maximum pairwise distance between all the pairs of one segment of the first and one segment of the second cluster
					}
				}
			}
			res = max;
			break;

		case sclib::linkageAverage:
			for (unsigned long int i = 0; i < pFirst->getSegmentCount(); i++) {
				for (unsigned long int j = 0; j < pSecond->getSegmentCount(); j++) {
					x = segmentMap[pFirst->getSegmentList(i, 1)];
					y = segmentMap[pSecond->getSegmentList(j, 1)];
					if (x > y) { //exchange indexes if pFirst has greater index than pSecond (because we only compute upper part below)
						t = x;
						x = y;
						y = t;
					}
					res += segmentDistanceMatrix[x][y]; //aggregate to find the average distance between all pairs of segments of pFirst and pSecond
					count++;
				}
			}
			res /= (double)(count);
			break;

		default:
			REPORT_ERROR(SVLIB_BadArg, "given linkage mode is unknown");
			res = std::numeric_limits<double>::max();
			break;
	}

	return res;
}

//====================================================================================================================
//	Compute the distance between to single segment models of 2 clusters
//====================================================================================================================
double SC_SpeakerClusterer::getModelDistance(SC_Model *pModel1, SC_Model *pBackground1, SV_Data *pFeatures1, unsigned long int segmentCount1, SC_Model *pModel2, SC_Model *pBackground2, SV_Data *pFeatures2, unsigned long int segmentCount2) {
	double res;
	unsigned short int oldGd = this->pTweak->distanceMeasure.groundDistance;

	if (this->pTweak->speakerClusterer.distanceMeasure == sclib::dmGLR) {
	  res = this->pDist->GLR(pModel1, pBackground1, pFeatures1, segmentCount1, pModel2, pBackground2, pFeatures2, segmentCount2);
  } else if (this->pTweak->speakerClusterer.distanceMeasure == sclib::dmCLR) {
  	res = this->pDist->CLR(pModel1, pBackground1, pFeatures1, segmentCount1, pModel2, pBackground2, pFeatures2, segmentCount2);
  } else if (sclib::bitTest(this->pTweak->speakerClusterer.distanceMeasure, sclib::dmBeigi) == true) {
		this->pTweak->distanceMeasure.groundDistance = this->pTweak->speakerClusterer.distanceMeasure;
    res = this->pDist->beigi(pModel1, pModel2);
		this->pTweak->distanceMeasure.groundDistance = oldGd;
  } else if (sclib::bitTest(this->pTweak->speakerClusterer.distanceMeasure, sclib::dmEMD) == true) {
		this->pTweak->distanceMeasure.groundDistance = this->pTweak->speakerClusterer.distanceMeasure;
		res = this->pDist->EMD(pModel1, pModel2);
		this->pTweak->distanceMeasure.groundDistance = oldGd;
	} else if (this->pTweak->speakerClusterer.distanceMeasure==sclib::mtSVM &&
		         (pModel1->Hdr.ModelType==sclib::mtSVM && pModel2->Hdr.ModelType==sclib::mtSVM) ||
						 (pModel1->Hdr.ModelType==sclib::mtTime && pModel2->Hdr.ModelType==sclib::mtTime && ((SC_Model_Time*)pModel1)->getSubModel()->Hdr.ModelType==sclib::mtSVM && ((SC_Model_Time*)pModel2)->getSubModel()->Hdr.ModelType==sclib::mtSVM)) {
		SV_Data *pScore = NULL;
		if (pModel1->Hdr.ModelType == sclib::mtSVM) {
			if (pModel2->Hdr.ModelType == sclib::mtSVM) {
				pScore = ((SC_Model_SVM*)pModel1)->TestModel((SC_Model_SVM*)(pModel2));
			} else {
				pScore = ((SC_Model_SVM*)pModel1)->TestModel(((SC_Model_SVM*)((SC_Model_Time*)pModel2)->getSubModel()));
			}
		} else {
			if (pModel2->Hdr.ModelType == sclib::mtSVM) {
				pScore = ((SC_Model_SVM*)((SC_Model_Time*)pModel1)->getSubModel())->TestModel((SC_Model_SVM*)(pModel2));
			} else {
				pScore = ((SC_Model_SVM*)((SC_Model_Time*)pModel1)->getSubModel())->TestModel(((SC_Model_SVM*)((SC_Model_Time*)pModel2)->getSubModel()));
			}
		}
		res = pScore->Mat[0][0];
		MFree_0D(pScore);
  } else {
		REPORT_ERROR(SVLIB_BadArg, "SC_SpeakerClusterer.getModelDistance: Specified distance-measure is unkonwn!");
    res = numeric_limits<double>::max();
	}

	return res;
}

//====================================================================================================================
//	Compute all pairwise distances between two single segment models, return a distance matrix (upper triangle) and
//  a map that maps segmentStart samples to rows/cols of the distance-matrix
//====================================================================================================================
double** SC_SpeakerClusterer::getSegmentDistanceMatrix(SC_Partition *pPartition, std::map<unsigned long int, unsigned long int> &segmentMap) {
	unsigned long int dim = 0, idx1 = 0, idx2 = 0;
	unsigned long int x = 0, y = 0, count = 0, overall;
	double **distanceMatrix, lastPercentage = 0.0;
  double minDistance = std::numeric_limits<double>::max();
	SC_Cluster *pFirstCluster, *pSecondCluster;
	SC_Model *pFirstModel, *pSecondModel;
	SC_MatrixFunctions matFunc;
	
  if (this->verbose == true) {printf("      ");}
	segmentMap.clear();

	//compute model-count in the complete partition
	pFirstCluster = pPartition->getClusters();
	while (pFirstCluster != NULL) {
		dim += pFirstCluster->getSegmentCount();
		pFirstCluster = pFirstCluster->Next;
	}

	if (dim > 0) {
		distanceMatrix = matFunc.zeros(dim, dim);
    overall = (dim * (dim-1)) / 2;

	  pFirstCluster = pPartition->getClusters();
		while (pFirstCluster != NULL) {
			idx1 = 0; //index of the current first segment in its cluster
			pFirstModel = pFirstCluster->getSpeakerModels();
			while (pFirstModel != NULL) {
				pSecondCluster = pFirstCluster;
				pSecondModel = (SC_Model*)(pFirstModel->Next);
				x = y + 1;
				while (pSecondCluster != NULL) {
					while (pSecondModel != NULL) {
						distanceMatrix[y][x] = getModelDistance(pFirstModel, pFirstCluster->getBackgroundModels(idx1), pFirstCluster->getSpeechFrames(idx1), pFirstCluster->getSegmentCount(), pSecondModel, pSecondCluster->getBackgroundModels(idx2), pSecondCluster->getSpeechFrames(idx2), pSecondCluster->getSegmentCount());
						x++;
						count++;
						if (this->verbose == true) {
							lastPercentage = sclib::printPercentage(overall, count, lastPercentage, 1.0, !count);
						}
						pSecondModel = (SC_Model*)(pSecondModel->Next);
					} //for all second models, beginning with current first's successor
					pSecondCluster = pSecondCluster->Next;
					idx2 = 0; //index of the current second segment in its cluster
					if (pSecondCluster != NULL) {
						pSecondModel = pSecondCluster->getSpeakerModels();
					}
				} //for all second clusters, beginning with the cirrent first one
				segmentMap[pFirstCluster->getSegmentList(idx1, 1)] = y;
				y++;
				idx1++;
				pFirstModel = (SC_Model*)(pFirstModel->Next);
			} //for all first models in first cluster
			pFirstCluster = pFirstCluster->Next;
		} //for all first clusters
	} //dim>0

  if (this->verbose == true) {
		sclib::printPercentage(1, 1, 0.0, 0.0, false);
  }

	return distanceMatrix;
}

//====================================================================================================================
// (re-)compute the distance-matrix needed by the clustering algorithm based on the merged models in the clusters;
// the new distance matrix is returned together with the x/y row/col-indexes of the closest pair of clusters and a 
// mapping
//====================================================================================================================
double** SC_SpeakerClusterer::getClusterDistanceMatrix(SC_Partition *pPartition, unsigned long int minDistanceIndex[2], double** &segmentBasedDistances, std::map<unsigned long int, unsigned long int> &segmentMap, double **oldDistanceMatrix) {
	unsigned long int x, y, count = 0, overall, dim = sclib::getListCount(pPartition->getClusters());
	unsigned long int oldIndexX, oldIndexY;
	double **distanceMatrix, lastPercentage = 0.0;
  double minDistance = std::numeric_limits<double>::max();
	SC_Cluster *pFirst, *pSecond;
	
	//compute initial distance matrix
	if (oldDistanceMatrix == NULL) { 
		
		//if average/single/complete linkage is choosen, we need all pairwise segment distances across all clusters
		if (this->pTweak->speakerClusterer.linkageMode != sclib::linkageMerge) {
			MFree_2D(segmentBasedDistances);
			segmentBasedDistances = getSegmentDistanceMatrix(pPartition, segmentMap);
		}

		if (this->verbose == true) {printf("      ");}

		if (dim >= 1) {
      MArray_2D(distanceMatrix, (long)(dim), (long)(dim), double, "SC_SpeakerClusterer.getClusterDistanceMatrix: distanceMatrix");
      overall = (dim * (dim-1)) / 2;

		  pFirst = pPartition->getClusters();
		  for (x = 0; x < dim; x++) {
			  pSecond = pFirst->Next;
        for (y = 0; y < dim; y++) { //if we have a symmteric distance measure (d(x,y)=d(y,x), as we assume/have), we only need to compute a triangular matrix without the main diagonal!
				  if (y > x) {
            assert(pFirst != NULL); 
            assert(pSecond != NULL);
						if (this->pTweak->speakerClusterer.linkageMode == sclib::linkageMerge) {
							distanceMatrix[x][y] = getModelDistance(pFirst->getMergedModel(), pFirst->getBackgroundModels(), pFirst->getSpeechFrames(), pFirst->getSegmentCount(), pSecond->getMergedModel(), pSecond->getBackgroundModels(), pSecond->getSpeechFrames(), pSecond->getSegmentCount());
						} else {
							distanceMatrix[x][y] = getClusterDistance(pFirst, pSecond, segmentBasedDistances, segmentMap);
						}

            if (this->verbose == true) {
							lastPercentage = sclib::printPercentage(overall, count, lastPercentage, 1.0, !count);
							count++;
            }

					  if (((x == 0) && (y == 1)) || (distanceMatrix[x][y] < minDistance)) { //initialize or update the minimum-distance
						  minDistance = distanceMatrix[x][y];
						  minDistanceIndex[0] = x;
						  minDistanceIndex[1] = y;
					  }
					  pSecond = pSecond->Next; 
          } else {
					  distanceMatrix[x][y] = 0.0;
          }
			  }
			  pFirst = pFirst->Next;
		  }
    }

	//recompute distance matrix under the assumption, that 2 rows/cols specified by minDistanceIndex[0]/minDistanceIndex[1] 
	//have been removed and a new one has been introduced at first position
	} else { 
		
	  if (this->verbose == true) {printf("      ");}
		MArray_2D(distanceMatrix, (long)(dim), (long)(dim), double, "SC_SpeakerClusterer.getClusterDistanceMatrix: distanceMatrix");
    overall = dim - 1;

		//only the distances involving the new first element must really be computed...
		pFirst = pPartition->getClusters();
		pSecond = pFirst->Next;
    for (x = 0; x < dim; x++) {
      distanceMatrix[x][0] = 0.0;
    }
		for (y = 1; y < dim; y++) { //if we have a symmteric distance measure (d(x,y)=d(y,x), as we assume/have), we only need to compute a triangular matrix without the main diagonal
			if (this->pTweak->speakerClusterer.linkageMode == sclib::linkageMerge) {
				distanceMatrix[0][y] = getModelDistance(pFirst->getMergedModel(), pFirst->getBackgroundModels(), pFirst->getSpeechFrames(), pFirst->getSegmentCount(), pSecond->getMergedModel(), pSecond->getBackgroundModels(), pSecond->getSpeechFrames(), pSecond->getSegmentCount());
			} else {
				distanceMatrix[0][y] = getClusterDistance(pFirst, pSecond, segmentBasedDistances, segmentMap);
			}
			pSecond = pSecond->Next; 

      if (this->verbose == true) {
				lastPercentage = sclib::printPercentage(overall, y, lastPercentage, 1.0, y==1);
      }
		}

		//...all the other distances are already there in the oldDistanceMatrix, they only may need to change places
		for (oldIndexX = 0; oldIndexX < dim+1; oldIndexX++) {
			for (oldIndexY = 0; oldIndexY < dim+1; oldIndexY++) {
        if (oldIndexX != minDistanceIndex[0] && oldIndexX != minDistanceIndex[1] && oldIndexY != minDistanceIndex[0] && oldIndexY != minDistanceIndex[1]) {
          x = oldIndexX + 1; //because there is a new zeroth element
          y = oldIndexY + 1;

          x = x - (oldIndexX > minDistanceIndex[0] ? 1 : 0) - (oldIndexX > minDistanceIndex[1] ? 1 : 0); //because the old minDistancePair is deleted
          y = y - (oldIndexY > minDistanceIndex[0] ? 1 : 0) - (oldIndexY > minDistanceIndex[1] ? 1 : 0); 

				  distanceMatrix[x][y]	= oldDistanceMatrix[oldIndexX][oldIndexY];
        }
			}
		}

		//compute the minimum-distance
		for (x = 0; x < dim; x++) {
			for (y = x+1; y < dim; y++) {
				if (((x == 0) && (y == 1)) || (distanceMatrix[x][y] < minDistance)) { //initialize or update the minimum-distance
					minDistance = distanceMatrix[x][y];
					minDistanceIndex[0] = x;
					minDistanceIndex[1] = y;
				}
			}
		}

		MFree_2D(oldDistanceMatrix);
	} //else

  if (this->verbose == true) {
		sclib::printPercentage(1, 1, 0.0, 0.0, false);
  }

	return distanceMatrix;
}

//====================================================================================================================
// Do the clustering with the initial segmentation found in pPartition; 
// return a score-object including the final partition and the list of all created partitions
// This function also alters the framelist in the groundtruth-object
// Invalid clusters (those who's segments are too short for clustering) are removed from pAllClusters in the 
// beginning and returned as a separate list in pInvalidClusters; they can be combined with the invalidClusters from
// the modeling process and used for identification...
//====================================================================================================================
SC_SpeakerScore_Clustering* SC_SpeakerClusterer::clusterSpeakers(SC_GroundTruth *pGT, SC_Cluster *&pAllClusters, SC_Cluster *&pInvalidClusters, long int progressCount, long int progressMax, bool (*doProceed)(long int &curr, long int max, unsigned int cnt)) {
	unsigned long int	uncertaintyDiameter, actualClusterCount, realClusterCount = pGT->getRealSpeakerCount(true);
	unsigned long int	x, minDistanceIndex[2]; //the x/y-indices of the smallest value of the distanceMatrix
	unsigned int mergeCounter = 1;
	double gc, **distanceMatrix = NULL, **segmentBasedDistances = NULL;
	bool proceed = true, immediateReturn = false;
	SC_Cluster *pFirst, *pSecond, *pNewCluster;
  SC_Partition *pFirstPartition, *pFinalPartition = NULL, *pPartitionHook, *pLastPartition;
  SC_SpeakerScore_Clustering *pScore = NULL;
	std::map<unsigned long int, unsigned long int> segmentMap;

	if (this->pTweak->speakerClusterer.linkageMode == sclib::linkageInteractive) {
		pFirstPartition = new SC_Partition(pAllClusters, 0.0, 0, true);
		pScore = new SC_SpeakerScore_Clustering(this->pTweak, pGT, pFirstPartition, NULL);
		pFinalPartition = interactiveClustering(pFirstPartition);
		pFinalPartition->storeIDs(pGT, true, true, true);
	} else {
		if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClustering) == true) {
			sclib::listOut("original_clusters.txt", pAllClusters, this->pTweak);
		}
		if (this->pTweak->speakerClusterer.constructNonOverlappingClusters == true) {
			constructNonOverlappingClusters(pAllClusters);
		}
		unsigned int oldMm = this->pTweak->cluster.mergeMode;
		if (this->pTweak->speakerClusterer.linkageMode != sclib::linkageMerge) {
			this->pTweak->cluster.mergeMode = sclib::mergeNone; //no merging needed if scoring is done based on the individual models!
		}

		pInvalidClusters = splitClusterList(pGT, pAllClusters);
		actualClusterCount = sclib::getListCount(pAllClusters);
		
		if (actualClusterCount >= 1) { //maybe all initial clusters where invalid...
			pFirstPartition = new SC_Partition(pAllClusters, 0.0, 0, true);
			pPartitionHook = pFirstPartition;
			pScore = new SC_SpeakerScore_Clustering(this->pTweak, pGT, pFirstPartition, NULL);
			
			//initialize the distance-matrix, remember min. distance, maybe do loading/saving
			char fileName[sclib::bufferSize];
			if (this->pTweak->speakerClusterer.firstDistanceMatrixPrefix!=NULL && strncmp(this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "", 1)!=0) {
				sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "distanceMat.dat");
			} else {
				sprintf(fileName, "");
			}
			if (sclib::fileExists(fileName) == true) {
				unsigned long int len, dim;
				distanceMatrix = sclib::loadMatrix(fileName, (double)(1.0), len, len);
				sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "segBasedDistanceMat.dat");
				if (sclib::fileExists(fileName) == true) { //exists only for complete/single/average linkage!
					segmentBasedDistances = sclib::loadMatrix(fileName, (double)(1.0), dim, dim);
				}
				sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "minDistanceIdx.dat");
				int *minDist = sclib::loadVector(fileName, (int)(1), dim);
				sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "segmentMap.dat");
				bool success = true;
				if (sclib::fileExists(fileName) == true) { //exists only for complete/single/average linkage!
					success = sclib::loadScalarMap(fileName, segmentMap);
				}
				if (distanceMatrix==NULL || len!=actualClusterCount || dim!=2 || minDist==NULL || success==false) {
					MFree_2D(distanceMatrix);
					MFree_2D(segmentBasedDistances);
					MFree_1D(minDist);
				} else {
					minDistanceIndex[0] = (unsigned long int)(minDist[0]);
					minDistanceIndex[1] = (unsigned long int)(minDist[1]);
					MFree_1D(minDist);
				}
			} 
			if (distanceMatrix == NULL) {
				distanceMatrix = getClusterDistanceMatrix(pFirstPartition, minDistanceIndex, segmentBasedDistances, segmentMap);
				if (fileName!=NULL && strncmp(fileName, "", 1)!=0) { //save matrix if wished
					sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "distanceMat.dat");
					sclib::saveMatrix(fileName, distanceMatrix, actualClusterCount, actualClusterCount);
					if (segmentBasedDistances != NULL) {
						sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "segBasedDistanceMat.dat");
						sclib::saveMatrix(fileName, segmentBasedDistances, actualClusterCount, actualClusterCount);
					}
					sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "minDistanceIdx.dat");
					int *minDist;
					MArray_1D(minDist, 2, int, "SC_SpeakerClusterer.clusterSpeakers: minDist");
					minDist[0] = minDistanceIndex[0];
					minDist[1] = minDistanceIndex[1];
					sclib::saveVector(fileName, minDist, 2);
					MFree_1D(minDist);
					if (segmentMap.size() > 0) {
						sprintf(fileName, "%s_%s\0", this->pTweak->speakerClusterer.firstDistanceMatrixPrefix, "segmentMap.dat");
						sclib::saveScalarMap(fileName, segmentMap);
					}
				}
			}
			if (distanceMatrix == NULL) {proceed = false;} //there was nothing to cluster!
			if ((this->pTweak->debug.debugMode & sclib::dbClustering) && (actualClusterCount > 1)) {sclib::vectorOut("globalCriterion.txt", "-----------", 10, false, this->pTweak);}

			//do the hierarchical clustering
			while ((proceed == true) && (actualClusterCount > 1))	{
		    
				//create a new partition as a copy of the last one, change it during the loop
				pPartitionHook->Next = new SC_Partition(pPartitionHook->getClusters(), 0.0, mergeCounter, false);
				if (pPartitionHook != pFirstPartition) {
					pPartitionHook->releaseWorkData(); //at least the speech-frames can, if needed, be claimed back later via SC_Partition.getSpeechFrames*()
				}
				pLastPartition = pPartitionHook;
				pPartitionHook = pPartitionHook->Next;
				
				//merge the two clusters whose distance is minimum, check if this is of any good
				pFirst = sclib::getListWithIndex(pPartitionHook->getClusters(), minDistanceIndex[0]);
				pSecond = sclib::getListWithIndex(pPartitionHook->getClusters(), minDistanceIndex[1]);
				pLastPartition->setFather(sclib::getListWithIndex(pLastPartition->getClusters(), minDistanceIndex[0])); //remember the history of mergences in the partitions
				pLastPartition->setMother(sclib::getListWithIndex(pLastPartition->getClusters(), minDistanceIndex[1]));
				pNewCluster = pFirst->mergeClusters(pSecond, mergeCounter++);
				pNewCluster->setIDconfidenceScore(distanceMatrix[minDistanceIndex[0]][minDistanceIndex[1]]); //use the distance between the parents as the confidence...
				gc = getClusterGlobalCriterion(pPartitionHook, pFirst, pSecond, pNewCluster);
				pPartitionHook->setGlobalCriterion(gc); //alway attache the value of the global criterion (gc) to the newly formed partition, even if it is in fact (because at the monent the new partition is just an exact copy of the last one) evaluated on the last one in case the whole partition is used and not just the 3 clusters (as for the WCD)
				proceed = getClusterContinuation(gc, actualClusterCount, realClusterCount);
		    			
				//some debug output
				if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClustering) == true) {sclib::matrixOut("dist.txt", distanceMatrix, actualClusterCount, actualClusterCount, this->pTweak);}
				if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbClustering) == true) {sclib::scalarOut("globalCriterion.txt", gc, this->pTweak);}

				//remove the two parent-clusters from the linked list, proceed with one fewer clusters
				if (proceed	== true) { 
					//remove the parents from this world
					assert(minDistanceIndex[0] < minDistanceIndex[1]);
					for (x = 0; x < 2; x++) {
						pPartitionHook->setClusters(sclib::removeFromList(pPartitionHook->getClusters(), minDistanceIndex[x]-x, actualClusterCount-x));
					}
					MFree_0D(pFirst);
					MFree_0D(pSecond);

					//put the new cluster in front of the linked list
					pNewCluster->Next = pPartitionHook->getClusters();
					pPartitionHook->setClusters(pNewCluster);
					pPartitionHook->setChild(pNewCluster);
					actualClusterCount--;

					//rearrange the distance-matrix
					distanceMatrix = getClusterDistanceMatrix(pPartitionHook, minDistanceIndex, segmentBasedDistances, segmentMap, distanceMatrix);
				} else {
					MFree_0D(pNewCluster);
					MFree_0D(pPartitionHook); //last mergence was bad according to the gc => remove this current partition and stay with the last one
					pLastPartition->Next = NULL;
				}

				if (doProceed != NULL && progressMax > 0) { //care for proper progress display in SC_Lib.audioSegmentation(), which is called by Videana
					immediateReturn = (*doProceed)(progressCount, progressMax, 1);
					if (immediateReturn == true) {
						break;
					}
				}
				if(this->verbose == true) {
					printf(".");
				}
			}

			MFree_2D(distanceMatrix);
			MFree_2D(segmentBasedDistances);
			segmentMap.clear();

			if (immediateReturn == false) {
				//pick the final partition according to the evaluated global criterion
				pFinalPartition = getClusterFinalPartition(pGT, pFirstPartition);
				pFinalPartition->getSpeechFrames(pFirstPartition); //reclaim speech frames so that a later started identification run can rebuild the models
				pFinalPartition->storeIDs(pGT, true, true, true);
			} else {
				pFinalPartition = pPartitionHook; //computation aborted, so just store the last partition as the final one (doesn't matter, we aborted anyway...)
			}
		} else {
			pFinalPartition = NULL;
  		pScore = new SC_SpeakerScore_Clustering(this->pTweak, pGT, NULL, NULL);
		}

		if (this->pTweak->speakerClusterer.linkageMode != sclib::linkageMerge) {
			this->pTweak->cluster.mergeMode = oldMm;
		}
	} //not interactive

  pScore->setFinalPartition(pFinalPartition);
	if (immediateReturn == false) { //this is time-consuming, so skip it if aborted (nobody wants the result anyway...)
		uncertaintyDiameter = pFinalPartition->getUncertaintyDiameter();
		pScore->calcScores(0, 0, uncertaintyDiameter);
	}

  MFree_0D(pDist);

  return pScore;
}

//====================================================================================================================
// Splits the given list of cluster-objetcs: all clusters containing altogether a too short speach portion (as 
// compared to the tweakable parameter speakerClusterer.speechSegLengthThreshold) are removed from the list and 
// returned as a newly linked list of "invalid" clusters
//====================================================================================================================
SC_Cluster* SC_SpeakerClusterer::splitClusterList(SC_GroundTruth *pGT, SC_Cluster *&pAllClusters) {
	unsigned long int threshold = pGT->getConverter()->ms2sample(this->pTweak->speakerClusterer.speechSegLengthThreshold);
	SC_Cluster *pInvalidClusters = NULL, *pValidClusters = pAllClusters, *pInvHook = NULL, *pVHook = NULL, *pHook = pAllClusters, *pTemp;

	if (threshold > 0) { //0 means no threshold here
		while (pHook != NULL) {
			pTemp = pHook->Next; //save
			pHook->Next = NULL; //destroy list link

			//depending on the overall samplecount in the cluster, attach it to the new valid or invalid list
			if (pHook->getSampleCount() < threshold) {
				if (pInvHook == NULL) {
					pInvalidClusters = pHook;
					pInvHook = pInvalidClusters;
				} else {
					pInvHook->Next = pHook; 
					pInvHook = pInvHook->Next;
				}

			} else {
				if (pVHook == NULL) {
					pValidClusters = pHook;
					pVHook = pValidClusters;
				} else {
					pVHook->Next = pHook; 
					pVHook = pVHook->Next;
				}
			}

			pHook = pTemp;
		}
	}

	pAllClusters = pValidClusters;
	return pInvalidClusters;
}

//====================================================================================================================
// Implements the idea from Kwon, Narayanan, "Robust Speaker Identification based on Selective use of Feature 
// Vectors", 2006: Having all speaker-models and recspective training-sets at hand, rebuild the models based solely on
// those vectors that have the maximum likelihood (among all given models) with the model they where used to 
// construct. Thus, models are hoped to be non-overlapping and short utterances might be scored more reliably.
//====================================================================================================================
void SC_SpeakerClusterer::constructNonOverlappingClusters(SC_Cluster *&pAllClusters) {
	SC_Cluster *pHook, *pCurrentCluster = pAllClusters;
	SC_Model *pModelHook, *pNonOverlapModelsHook, *pBackgroundHook;
	SV_Data *pFeatures, **pNonOverlapFrames, *pNonOverlapHook, *pScore, *pFeature;
	int t, d, maxCluster, clusterIdx, currentIdx = 0, cnt, clusterCount = sclib::getListCount(pAllClusters);
	double maxScore;
	std::list<int> overlapIndices, nonOverlapIndices;
	SC_ModelHandler handler(this->pTweak, true);
	SC_Model **pNonOverlapModels;

	pFeature = new SV_Data(1, pAllClusters->getSpeechFrames()->Col);
	pFeature->Hdr = pAllClusters->getSpeechFrames()->Hdr;

	//save the new findings till the end of the process so all vectors are evaluated on the same basis; then, change the original clusters
	MArray_1D(pNonOverlapFrames, clusterCount, SV_Data*, "SC_SpeakerClusterer.constructNonOverlappingClusters: pNonOverlapFrames");
	MArray_1D(pNonOverlapModels, clusterCount, SC_Model*, "SC_SpeakerClusterer.constructNonOverlappingClusters: pNonOverlapModels");

	while (pCurrentCluster != NULL) {
		//find, for the current speaker (cluster), all feature vectors that genuinely belong to it and not to others; preserve segment-structure!
		pFeatures = pCurrentCluster->getSpeechFrames();
		pBackgroundHook = pCurrentCluster->getBackgroundModels();
		pNonOverlapFrames[currentIdx] = NULL;
		pNonOverlapModels[currentIdx] = NULL;
		while (pFeatures != NULL) { //we have a linked list of segments here
			nonOverlapIndices.clear();
			
			for (t = 0; t < pFeatures->Row; t++) {
				for (d = 0; d < pFeatures->Col; d++) {
					pFeature->Mat[0][d] = pFeatures->Mat[t][d];
				}

				//find maximum likelihood model for this vector
				maxScore = std::numeric_limits<double>::max()*-1.0;
				pHook = pAllClusters;
				clusterIdx = 0;
				while (pHook != NULL) {
					pModelHook = pHook->getSpeakerModels();
					while (pModelHook != NULL) { //test against all single models in the current cluster
						pScore = pModelHook->TestModel(pFeature);
						if (pScore->Mat[0][0] > maxScore) {
							maxScore = pScore->Mat[0][0];
							maxCluster = clusterIdx;
						}
						MFree_0D(pScore);
						pModelHook = (SC_Model*)(pModelHook);
					}
					
					pHook = pHook->Next;
					clusterIdx++;
				}

				//remember which rows did or did not overlap with other speaker's models
				if (maxCluster == currentIdx) {
					nonOverlapIndices.push_back(t);
				}
			} //for t

			//build a featureset of non-overlapping frames and include it in a linked list that corresponds with segment-boundaries stored in the cluster
			if (pNonOverlapFrames[currentIdx] == NULL) {
				pNonOverlapFrames[currentIdx] = (nonOverlapIndices.empty()==false) ? new SV_Data((int)(nonOverlapIndices.size()), pFeatures->Col) : new SV_Data(); //create a "dummy" if no frames didn't overlap
				pNonOverlapHook = pNonOverlapFrames[currentIdx];
			} else {
				pNonOverlapHook->Next = (nonOverlapIndices.empty()==false) ? new SV_Data((int)(nonOverlapIndices.size()), pFeatures->Col) : new SV_Data();
				pNonOverlapHook = pNonOverlapHook->Next;
			}
			pNonOverlapHook->Hdr = pFeatures->Hdr;
			cnt = 0;
			while (nonOverlapIndices.empty() == false) {
				for (d = 0; d < pNonOverlapHook->Col; d++) {
					pNonOverlapHook->Mat[cnt][d] = pFeatures->Mat[nonOverlapIndices.front()][d];
				}
				nonOverlapIndices.pop_front();
				cnt++;
			}

			//re-build corresponding list of segment-models
			if (pNonOverlapModels[currentIdx] == NULL) {
				pNonOverlapModels[currentIdx] = handler.buildModel(pNonOverlapFrames[currentIdx], pBackgroundHook, sclib::modeForeground, 1);
				pNonOverlapModelsHook = pNonOverlapModels[currentIdx];
			} else {
				pNonOverlapModelsHook->Next = handler.buildModel(pNonOverlapFrames[currentIdx], pBackgroundHook, sclib::modeForeground, 1);
				pNonOverlapModelsHook = (SC_Model*)(pNonOverlapModelsHook->Next);
			}

			pFeatures = pFeatures->Next;
			pBackgroundHook = (SC_Model*)(pBackgroundHook->Next);
			if ((pFeatures!=NULL && pBackgroundHook==NULL) || (pFeatures==NULL && pBackgroundHook!=NULL)) { //segment-based data, i.e. feature-sets and background-models, have to correspond
				REPORT_ERROR(SVLIB_BadData, "the list of features and the list of background-models in the cluster is out of sync");
				break;
			}
		} //while pFeatures!=NULL

		pCurrentCluster = pCurrentCluster->Next;
		currentIdx++;
	}

	//exchange speech frames and models in the original clusters, re-build merged models in the clusters
	pCurrentCluster = pAllClusters;
	currentIdx = 0;
	while (pCurrentCluster != NULL) {
		pCurrentCluster->setSpeakerModels(pNonOverlapModels[currentIdx]);
		sclib::destructLinkedList(pNonOverlapModels[currentIdx]); //is really copyied above
		pCurrentCluster->setSpeechFrames(pNonOverlapFrames[currentIdx]); //cleaning up what is of no more need is done inside here
		pCurrentCluster->reBuildMergedModel(); //build the overall cluster model according to the choosen criterion

		pCurrentCluster = pCurrentCluster->Next;
		currentIdx++;
	}

	MFree_0D(pFeature);
	MFree_1D(pNonOverlapFrames);
	MFree_1D(pNonOverlapModels);

	return;
}

//====================================================================================================================
//  This method returns, for a given partition-list (that depicts the succession of a clustering process), the alpha
//  (BIC penalty factor) that, for the given values of an objective function per list element that must be minimized, 
//  yields the best result when BIC is used as a termination criterion for clustering. In finalPartitionIdx, the index 
//  of the partition that would be choosen as final by BIC with the computed alpha is returned.
//====================================================================================================================
double SC_SpeakerClusterer::tuneBIC(SC_Partition *pPartitionList, double *objectiveFunction, unsigned long int &finalPartitionIdx) {
	SC_MatrixFunctions matFunc;
	SC_Partition *pHook = pPartitionList;
	SC_Cluster *pFirst, *pSecond;
	SV_Data *pFrames;
	double alpha = 1.0, *x, *y, bic, minObj, lastMinObj;
	unsigned long int i, j, idx = 0, listCount = sclib::getListCount(pPartitionList);
	bool success;

	//get the two parts/terms, x and y, of the bic formula, bic=x-alpha*y, for each partition except the first (has no predecessor, thus no bic)
	MArray_1D(x, listCount, double, "SC_DistanceMeasures.tuneBIC: x");
	MArray_1D(y, listCount, double, "SC_DistanceMeasures.tuneBIC: y");
	x[0] = 0.0;
	y[0] = 0.0;
	while (pHook != NULL) {
		pFirst = pHook->getFather();
		pSecond = pHook->getMother();
		if (pFirst!=NULL && pSecond!=NULL) {
			if (pHook != pPartitionList) { //if this is not the first partition, the clusters have no workdata, so copy it from the first partition and destroy it afterwards
				pFrames = pHook->getSpeechFramesOfCluster(pFirst->getSpeakerID(), pPartitionList);
				pFirst->setSpeechFrames(pFrames);
				pFrames = pHook->getSpeechFramesOfCluster(pSecond->getSpeakerID(), pPartitionList);
				pSecond->setSpeechFrames(pFrames);
			}
			bic = this->pDist->BIC(pFirst, pSecond, alpha, &x[idx+1], &y[idx+1]); //the result of merging these two clusters is the bic value (global criterion) of the next partition; the last one has no parents, so we don't get in trouble with the indexes here
			if (pHook != pPartitionList) {
				pFirst->releaseWorkData();
				pSecond->releaseWorkData();
			}
		}
		idx++;
		pHook = pHook->Next;
	}	
	//sclib::vectorOut("x.txt", x, listCount, false, this->pTweak);
	//sclib::vectorOut("y.txt", y, listCount, false, this->pTweak);

	//for each partition in the order of their objective function value (minimum first), see if it is "reachable"
	for (i = 0; i < listCount; i++) {
		//find the i'th minimum objective function value and the partition (idx) in which it occurs		
		minObj = matFunc.min(objectiveFunction, listCount, &idx, ((i==0)?NULL:&lastMinObj)); 

		if (idx == listCount-1) { //the last partition (i.e. one big cluster) is best, so return an alpha that will always produce such a result
			alpha = -1.0 * std::numeric_limits<double>::max();
			finalPartitionIdx = listCount - 1;
			break;
		} else {
			//compute alpha needed in bic computation to reach a bic-value <0.0 in the successor of the idx'th partition
			alpha = sclib::incrementDouble(x[idx+1]/y[idx+1]); // + std::numeric_limits<double>::epsilon(); //solved the inequality "x-alpha*y < 0" => "alpha > x/y"; the epsilon accounts for the greater-sign

			//see if this is really the first partition in the list where bic falls below 0 with this alpha;
			//if so, we finished optimization;
			//otherwise, we sadly cannot "reach" this minimum objective function partition with the bic stopping criterion and have to look for the 2ns (3rd, 4th, ...) best solution
			success = true;
			for (j = 1; j <= idx; j++) { //start with the 2nd partition because the first one doesn't have a predecessor and thus no bic
				bic = x[j] - alpha*y[j];
				if (bic < 0.0) {
					success = false;
					break;
				}
			}
			if (success == true) {
				finalPartitionIdx = idx;
				break;
			}
		}

		lastMinObj = minObj;
	}

	MFree_1D(x);
	MFree_1D(y);

	return alpha;
}

//====================================================================================================================
//  This method returns, for a given partition-list (that depicts the succession of a clustering process), the ICR 
//  threshold that, for the given values of an objective function per list element that must be minimized, 
//  yields the best result when ICR is used as a termination criterion for clustering. In finalPartitionIdx, the index 
//  of the partition that would be choosen as final by ICR with the computed threshold is returned.
//====================================================================================================================
double SC_SpeakerClusterer::tuneICR(SC_Partition *pPartitionList, double *objectiveFunction, unsigned long int &finalPartitionIdx) {
	SC_MatrixFunctions matFunc;
	SC_Partition *pHook = pPartitionList;
	SC_Cluster *pFirst, *pSecond;
	SV_Data *pFrames;
	double threshold = 1.0, *icr, minObj, lastMinObj;
	unsigned long int i, j, idx = 0, listCount = sclib::getListCount(pPartitionList);
	bool success;

	//get the icr for each partition
	MArray_1D(icr, listCount, double, "SC_SpeakerClusterer.tuneICR: icr");
	icr[0] = 0.0;
	while (pHook != NULL) {
		pFirst = pHook->getFather();
		pSecond = pHook->getMother();
		if (pFirst!=NULL && pSecond!=NULL) {
			if (pHook != pPartitionList) { //if this is not the first partition, the clusters have no workdata, so copy it from the first partition and destroy it afterwards
				pFrames = pHook->getSpeechFramesOfCluster(pFirst->getSpeakerID(), pPartitionList);
				pFirst->setSpeechFrames(pFrames);
				pFrames = pHook->getSpeechFramesOfCluster(pSecond->getSpeakerID(), pPartitionList);
				pSecond->setSpeechFrames(pFrames);
			}
			icr[idx+1] = this->pDist->ICR(pFirst, pSecond); //the result of merging these two clusters is the icr value (global criterion) of the next partition; the last one has no parents, so we don't get in trouble with the indexes here
			if (pHook != pPartitionList) {
				pFirst->releaseWorkData();
				pSecond->releaseWorkData();
			}
		}
		idx++;
		pHook = pHook->Next;
	}	
	//sclib::vectorOut("gc.txt", icr, idx, false, this->pTweak);

	//going from last to first partition (i.e. in reverse order), ICR returns the partition whichs gc value is below the threshold
	//for each partition in the order of their objective function value (minimum first), see if it is "reachable"
	for (i = 0; i < listCount; i++) {
		//find the i'th minimum objective function value and the partition (idx) in which it occurs		
		minObj = matFunc.min(objectiveFunction, listCount, &idx, ((i==0)?NULL:&lastMinObj)); 

		if (idx == listCount-1) { //the last partition (i.e. one big cluster) is best, so return a threshold that will always produce such a result
			threshold = 1.0 * std::numeric_limits<double>::max();
			finalPartitionIdx = listCount - 1;
			break;
		} else {
			//compute threshold needed to have an icr smaller than it in the idx'th partition
			threshold = sclib::incrementDouble(icr[idx]); // + std::numeric_limits<double>::epsilon();

			//see if this is really the first partition in the list where icr falls below the threshold
			//if so, we finished optimization;
			//otherwise, we sadly cannot "reach" this minimum objective function partition with the icr stopping criterion and have to look for the 2ns (3rd, 4th, ...) best solution
			success = true;
			for (j = idx+1; j < listCount; j++) { //the remaining partitions mustn't exceed the thresholod, otherwise this optimal partition isn't reachable by icr
				if (icr[j] < threshold) {
					success = false;
					break;
				}
			}
			if (success == true) {
				finalPartitionIdx = idx;
				break;
			}
		}

		lastMinObj = minObj;
	}

	MFree_1D(icr);

	return threshold;
}

//====================================================================================================================
//  This method returns, for a given partition-list (that depicts the succession of a clustering process), the WCD 
//  penalty factor that, for the given values of an objective function per list element that must be minimized, 
//  yields the best result when WCD is used as a termination criterion for clustering. In finalPartitionIdx, the index 
//  of the partition that would be choosen as final by WCD with the computed penalty factor is returned.
//====================================================================================================================
double SC_SpeakerClusterer::tuneWCD(SC_Partition *pPartitionList, double *objectiveFunction, unsigned long int &finalPartitionIdx) {
	SC_MatrixFunctions matFunc;
	SC_Partition *pHook = pPartitionList;
	double p, penaltyFactor = 1.0, *x, *y, wcd, minObj, lastMinObj, greater, smaller;
	unsigned long int i, j, idx = 0, listCount = sclib::getListCount(pPartitionList);
	bool success;

	//get the wcd for each partition
	MArray_1D(x, listCount, double, "SC_SpeakerClusterer.tuneWCD: x");
	MArray_1D(y, listCount, double, "SC_SpeakerClusterer.tuneWCD: y");
	while (pHook != NULL) {
		if (pHook != pPartitionList) { //if this is not the first partition, the clusters have no workdata, so copy it from the first partition and destroy it afterwards
			pHook->getSpeechFrames(pPartitionList);
		}
		wcd = this->pDist->withinClusterDispersion(pHook, &x[idx], &y[idx]);
		if (pHook != pPartitionList) {
			pHook->releaseWorkData();
		}
		idx++;
		pHook = pHook->Next;
	}	
	//sclib::vectorOut("x.txt", x, listCount, false, this->pTweak);
	//sclib::vectorOut("y.txt", y, listCount, false, this->pTweak);

	//for each partition in the order of their objective function value (minimum first), see if it is "reachable"
	for (i = 0; i < listCount; i++) {
		//find the i'th minimum objective function value and the partition (idx) in which it occurs		
		minObj = matFunc.min(objectiveFunction, listCount, &idx, ((i==0)?NULL:&lastMinObj)); 

		//compare the idx'th partitions' wcd with each other's wcd - there must be a common penalty factor that 
		//makes the idx'th wcd minimal among all, or we can't reach this optimal partition
		success = true;
		greater = std::numeric_limits<double>::max() * -1.0;
		smaller = std::numeric_limits<double>::max();
		for (j = 0; j < listCount; j++) {
			if (j != idx) {
				p = (x[idx]-x[j]) / (y[j]-y[idx]); //for idx'th partition to have a smaller wcd as the j'th one, the following inequality must hold: x_idx-x_j < p*(y_j-y_idx)
				if (y[j]-y[idx] >= 0.0) { //care for sign-change in the inequality due to division by a negative number
					p = sclib::incrementDouble(p); //p + std::numeric_limits<double>::epsilon();
					greater = sclib::max(greater, p);
				} else {
					p = sclib::decrementDouble(p); //p - std::numeric_limits<double>::epsilon();
					smaller = sclib::min(smaller, p);
				}
				if (smaller < greater) {
					success = false; //contradiction in the set of unequalities => optimal partition unreachable!
					break;
				}
			}
		}
		if (success == true) {
			finalPartitionIdx = idx;
			if (smaller<std::numeric_limits<double>::max() && greater>std::numeric_limits<double>::max()*-1.0) {
				penaltyFactor = greater + (smaller-greater)/2.0; //the middle between them
			} else if (smaller < std::numeric_limits<double>::max()) {
				penaltyFactor = smaller;
			} else if (greater > std::numeric_limits<double>::max()*-1.0) {
				penaltyFactor = greater;
			} else {
				penaltyFactor = 0.0; //should never happen
			}
			break;
		}

		lastMinObj = minObj;
	}

	MFree_1D(x);
	MFree_1D(y);

	return penaltyFactor;
}

//====================================================================================================================
// Perform speaker clustering manually in an interactive manner, return the final partition
//====================================================================================================================
SC_Partition* SC_SpeakerClusterer::interactiveClustering(SC_Partition* &pFirstPartition) {
	SC_Partition *pFinalPartition, *pCurrentPartition, *pLastPartition;
	SC_Cluster *pHook, *pFirst, *pSecond, *pNew;
	bool finished, firstTime, restart;
	int id[2], origId[2], x, special;
	unsigned long int mergeCounter, actualClusterCount;

	unsigned int oldMm = this->pTweak->cluster.mergeMode;
	this->pTweak->cluster.mergeMode = sclib::mergeNone;

	do {
		int map[3][14] = {{401, 1301, 1401, 701, 1201, 601, 101, 501, 801, 301, 1001, 1101, 901, 201},
											{1301, 101, 1401, 1101, 1001, 801, 201, 301, 1201, 401, 901, 701, 501, 601},
											{1301, 601, 501, 801, 1101, 301, 701, 201, 401, 1401, 1201, 101, 901, 1001}};
		finished = false;
		restart = false;
		pCurrentPartition = pFirstPartition;
		pFinalPartition = pFirstPartition;
		mergeCounter = 0;
		actualClusterCount = sclib::getListCount(pCurrentPartition->getClusters());

		printf("\n\nWelcome to interactive clustering: \n");
		printf("Use special mappings for experiment (-1, 3, 2, 1)? ");
		scanf("%d", &special);
		
		do {
			printf("\ncurrent clusters:");
			pHook = pCurrentPartition->getClusters();
			x = 0;
			while (pHook != NULL) {
				printf("\nid: %d\tname: %s", x++, pHook->getSpeakerName());
				pHook = pHook->Next;
			}

			printf("\nnext mergence ('id, id' or '-1' to finish, '-2' to restart): ");
			firstTime = true;
			do { //repeat until correct input is found
				if (firstTime == false) {
					printf("\nonce again, couldn't recognize format: ");
				}
				scanf("%d, %d", &id[0], &id[1]);
				if (id[0] == -1) {
					finished = true;
					break;
				}
				if (id[0] == -2) {
					finished = true;
					restart = true;
					break;
				}
				firstTime = false;
			} while ((special==-1 && (id[0]<0 || (unsigned long int)(id[1])>sclib::getListCount(pCurrentPartition->getClusters()) || id[0]>=id[1])) ||
							 (special>0 && special<4 && (id[0]<1 || id[1]<1 || id[0]>14 || id[1]>14)));

			if (finished != true) { //create new partition, delete parents, add child
				if (special > 0) {
					origId[0] = id[0];
					origId[1] = id[1];
					id[0] = map[special-1][id[0]-1];
					id[0] = pCurrentPartition->getClusterIdxByID(id[0]);
					id[1] = map[special-1][id[1]-1];
					id[1] = pCurrentPartition->getClusterIdxByID(id[1]);
					if (id[0] > id[1]) {
						x = id[0];
						id[0] = id[1];
						id[1] = x;
						map[special-1][origId[0]-1] = map[special-1][origId[1]-1];
					} else if (id[0] < id[1]) {
						map[special-1][origId[1]-1] = map[special-1][origId[0]-1]; //after the mergence the merged cluster has the id of the first one, and the 2nd one belongs to it
					}
				}	

				if (id[0] != id[1]) {
					pCurrentPartition->Next = new SC_Partition(pCurrentPartition->getClusters(), 0.0, mergeCounter, false);
					if (pCurrentPartition != pFirstPartition) {
						pCurrentPartition->releaseWorkData();
					}
					pLastPartition = pCurrentPartition;
					pCurrentPartition = pCurrentPartition->Next;
					pFirst = sclib::getListWithIndex(pCurrentPartition->getClusters(), id[0]);
					pSecond = sclib::getListWithIndex(pCurrentPartition->getClusters(), id[1]);
					pLastPartition->setFather(sclib::getListWithIndex(pLastPartition->getClusters(), id[0]));
					pLastPartition->setMother(sclib::getListWithIndex(pLastPartition->getClusters(), id[1]));
					pNew = pFirst->mergeClusters(pSecond, mergeCounter++);
					for (x = 0; x < 2; x++) {
						pCurrentPartition->setClusters(sclib::removeFromList(pCurrentPartition->getClusters(), id[x]-x, actualClusterCount-x));
					}
					MFree_0D(pFirst);
					MFree_0D(pSecond);
					pNew->Next = pCurrentPartition->getClusters();
					pCurrentPartition->setClusters(pNew);
					pCurrentPartition->setChild(pNew);
					actualClusterCount--;
					pFinalPartition = pCurrentPartition;
				}
			}

			if (actualClusterCount < 2) {
				finished = true;
			}
		} while (finished == false);
		
		if (restart == true) {
			sclib::destructLinkedList(pFirstPartition->Next);
		}
	} while (restart == true);

	this->pTweak->cluster.mergeMode = oldMm;

	return pFinalPartition;
}

//====================================================================================================================
// Perform #repetitions runs of random clustering and print the results to the file; "random" here means (thats not
// so obvious because its a 2-stage process, first choosing cluster-count/-sizes and then assignung clusters to 
// individual segments, and, when compared to a human, a human will not act totally radnom on the first task): first
// a number of singel segments for a cluster is choosen according to the given distribution, then segments are 
// assigned to this cluster randomly. this process is iterated until no unassigned segments remain. the distribution 
// is expected as an array where the content at index i gives the probability that a cluster of size i is present in
// the clustering; maxClusterSize, then, is the maximal index into the array (yes, 0-size is included but should be
// 0); the repetitions exist to find a true random result by performing a kind of monte carlo simulation experiment 
// and averaging the random scores reached in each radnom run
//====================================================================================================================
void SC_SpeakerClusterer::randomClustering(const char *fileName, SC_GroundTruth *pGT, SC_Cluster *pClusterList, double *clusterSizeDistribution, unsigned int maxClusterSize, unsigned int repetitions) {
	SC_Partition *pCurrentPartition;
	SC_Cluster *pHook, *pFirst, *pSecond, *pNew;
	long int idx, id[2];
	unsigned int clusterSize;
	unsigned long int mergeCounter, s1, s2;
	std::vector<long int> remainingSegments;
	SC_SpeakerScore_Clustering *pScore = new SC_SpeakerScore_Clustering(this->pTweak, pGT, NULL, NULL, false);
	unsigned int clusters, correctHull, lines;
	int x;
	double lastPercentage = 0.0;
	std::map<std::pair<long int, long int>, unsigned long int> correctMap;
	std::map<std::pair<long int, long int>, unsigned long int>::iterator i;
	char *secondFileName = sclib::exchangeFileExtension(fileName, ".speakers.txt");

	unsigned int oldMm = this->pTweak->cluster.mergeMode;
	this->pTweak->cluster.mergeMode = sclib::mergeNone;

	//initialize the file with the speaker mappings (displaying for each pair of clusters if it was joined correctly)
	pFirst = pClusterList;
	while (pFirst != NULL) {
		pSecond = pFirst->Next;
		while (pSecond != NULL) {
			if (pGT->getSamplesSpeakerID(pFirst->getSegmentList(0, 1), sclib::modeGroundtruth) == pGT->getSamplesSpeakerID(pSecond->getSegmentList(0, 1), sclib::modeGroundtruth)) { //compare the gt-speaker-id of the segment-start-sample in each pair of segments from the two clusters to be merged
				correctMap[std::make_pair(pFirst->getSpeakerID(), pSecond->getSpeakerID())] = 0;
				break;
			}
			pSecond = pSecond->Next;
		}
		pFirst = pFirst->Next;
	}
	for (i = correctMap.begin(); i != correctMap.end(); ++i) {
		sclib::tupelOut(secondFileName, i->first.first, i->first.second, this->pTweak, "_", "\t");
	}
	sclib::stringOut(secondFileName, "", this->pTweak, "\n");

	clusterSizeDistribution[0] = 0.0; //make sure that clusters of size 0 have a probability of 0.0

	printf("\nSimulating random clustering: ");
	sclib::stringOut(fileName, "#clusters\t#correct(hull)\t#lines\toverallRecall\t\toverallPrecision\t\tmissclassRate\t\tavgClusterPurity\t\tavgClassPurity\t\tpurity\t\trandIndex\t\tbbnMetric\t\tder", this->pTweak, "\n");

	for (unsigned int r = 0; r < repetitions; r++) {
		pCurrentPartition = new SC_Partition(pClusterList, 0.0, 0, false, true);
		mergeCounter = 0;

		lines = 0;
		correctHull = 0;
		clusters = 0;
		correctMap.clear();

		//initialize the correctMap in the same order as the header for the file was printed
		pFirst = pClusterList;
		while (pFirst != NULL) {
			remainingSegments.push_back(pFirst->getSpeakerID());
			pSecond = pFirst->Next;
			while (pSecond != NULL) {
				if (pGT->getSamplesSpeakerID(pFirst->getSegmentList(0, 1), sclib::modeGroundtruth) == pGT->getSamplesSpeakerID(pSecond->getSegmentList(0, 1), sclib::modeGroundtruth)) { //compare the gt-speaker-id of the segment-start-sample in each pair of segments from the two clusters to be merged
					correctMap[std::make_pair(pFirst->getSpeakerID(), pSecond->getSpeakerID())] = 0;
					break;
				}
				pSecond = pSecond->Next;
			}
			pFirst = pFirst->Next;
		}

		do {
			//pick a size for the cluster to be created, according to the given distribution of cluster-sizes
			clusterSize = sclib::min((unsigned int)(remainingSegments.size()), sclib::drawIndexFromDistribution(clusterSizeDistribution, maxClusterSize+1));
			clusters++;
			
			//pick a first segment for this cluster
			idx = sclib::round(sclib::rand(0.0, remainingSegments.size()-1.0));
			id[0] = pCurrentPartition->getClusterIdxByID(remainingSegments[idx]);
			remainingSegments.erase(remainingSegments.begin()+idx);
			pFirst = sclib::getListWithIndex(pCurrentPartition->getClusters(), id[0]);
			
			//subsequently assign more segments to the cluster until the final size is reached
			if (clusterSize > 0) {
				for (unsigned int c = 0; c < clusterSize-1; c++) {
					//pick second segment
					idx = sclib::round(sclib::rand(0.0, remainingSegments.size()-1.0));
					id[1] = pCurrentPartition->getClusterIdxByID(remainingSegments[idx]);
					remainingSegments.erase(remainingSegments.begin()+idx);
					pSecond = sclib::getListWithIndex(pCurrentPartition->getClusters(), id[1]);

					//the first cluster must have the smaller index!
					if (id[0] > id[1]) {
						idx = id[0];
						id[0] = id[1];
						id[1] = idx;

						pHook = pFirst;
						pFirst = pSecond;
						pSecond = pHook;
					}

					//check if this mergence might be correct
					//TODO: this assumes that each final cluster is composed of at least 2 segments; it doesn't count single segments!
					for (s1 = 0; s1 < pFirst->getSegmentCount(); s1++) {
						for (s2 = 0; s2 < pSecond->getSegmentCount(); s2++) {
							if (pGT->getSamplesSpeakerID(pFirst->getSegmentList(s1, 1), sclib::modeGroundtruth) == pGT->getSamplesSpeakerID(pSecond->getSegmentList(s2, 1), sclib::modeGroundtruth)) { //compare the gt-speaker-id of the segment-start-sample in each pair of segments from the two clusters to be merged
								if (pFirst->getSegmentList(s1, 0) < pSecond->getSegmentList(s2, 0)) {
									correctMap[std::make_pair(pFirst->getSegmentList(s1, 0), pSecond->getSegmentList(s2, 0))]++;
								} else {
									correctMap[std::make_pair(pSecond->getSegmentList(s2, 0), pFirst->getSegmentList(s1, 0))]++;
								}
								correctHull++;
								break;
							}
						}
					}
					
					//merge it with the rest of the cluster, update the lnked list
					pNew = pFirst->mergeClusters(pSecond, mergeCounter++);
					for (x = 0; x < 2; x++) {
						pCurrentPartition->setClusters(sclib::removeFromList(pCurrentPartition->getClusters(), id[x]-x, sclib::getListCount(pCurrentPartition->getClusters())));
					}
					MFree_0D(pFirst);
					MFree_0D(pSecond);
					pNew->Next = pCurrentPartition->getClusters();
					pCurrentPartition->setClusters(pNew);
				
					//the complete cluster (up to this stage) is pointed to by pFirst
					pFirst = pNew;
					id[0] = 0; //the first cluster is now at index 0
					lines++;
				}
			}
		} while (remainingSegments.size() > 0);

		pCurrentPartition->storeIDs(pGT, true, true, true);
		pScore->setPartitionList(pCurrentPartition);
		pScore->setFinalPartition(pCurrentPartition);
		pScore->calcScores();

		sclib::scalarOut(fileName, clusters, this->pTweak, true, "\t");
		sclib::scalarOut(fileName, correctHull, this->pTweak, true, "\t");
		sclib::scalarOut(fileName, lines, this->pTweak, true,"\t");
		sclib::scalarOut(fileName, pScore->getOverallRecall(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getOverallRecall(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getOverallPrecision(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getOverallPrecision(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getMissclassificationRate(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getMissclassificationRate(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getAverageClusterPurity(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getAverageClusterPurity(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getAverageClassPurity(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getAverageClassPurity(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getPurity(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getPurity(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getRandIndex(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getRandIndex(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getBBNmetric(sclib::ceSample), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getBBNmetric(sclib::ceSegment), this->pTweak, true, "\t");
		sclib::scalarOut(fileName, pScore->getDER(), this->pTweak, true, "\n");
		for (i = correctMap.begin(); i != correctMap.end(); ++i) {
			sclib::scalarOut(secondFileName, i->second, this->pTweak, false, "\t");
		}
		sclib::stringOut(secondFileName, "", this->pTweak, "\n");

		MFree_0D(pCurrentPartition);
		lastPercentage = sclib::printPercentage(repetitions-1, r, lastPercentage, 0.0, r==0);
	}
	lastPercentage = sclib::printPercentage(repetitions-1, repetitions-1, lastPercentage, 0.0, false);

	this->pTweak->cluster.mergeMode = oldMm;
	
	MFree_1D(secondFileName);
	MFree_0D(pScore);

	return;
}
