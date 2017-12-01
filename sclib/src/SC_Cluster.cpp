/**************************************************************************/
/*	A Cluster is meant (in the end) to represent all Segments of one			*/
/*	Speaker.																															*/
/*																																				*/
/*    Responsibility:																											*/
/*      - encapsulates a cluster of speech; therefor it holds	(links to)	*/
/*				- models for speech and background															*/
/*				- the original speech-frames																		*/
/*				- a list which orders the frames in the timeline of the video		*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 17.05.2004																								*/
/**************************************************************************/

#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include "SC_Cluster.h"
#include "SC_Aux.h"
#include "SC_MatrixFunctions.h"
#include "SC_ModelHandler.h"
#include "SC_FeatureHandler.h"
#include <SV_Error.h>
#include <SV_DataIO.h>

//====================================================================================================================
//	constructor for a new cluster
//====================================================================================================================
SC_Cluster::SC_Cluster(SC_TweakableParameters* pTweak, SV_Data* pSpeech, unsigned long int* speechStart, unsigned long int* speechEnd, SC_Model* pMergedModel, SC_Model* pBackgroundModels, SC_Model* pSpeakerModels, unsigned long int segmentCount, long int speakerID, const char* speakerName, bool justLink) {
	this->Next = NULL;
	this->speakerID = speakerID;
	this->segmentCount = (speechEnd > 0 || pSpeech != NULL) ? segmentCount : 0;
	this->pTweak = pTweak;
  this->classSig = 196275;
  this->IDconfidenceScore = 1.0; //in the beginning, the given speaker-id is for sure...
	this->ownerID = sclib::noSpeaker; //no owner until scoring...
	this->ownerIsCorrect = false; 
	this->parentID[0] = sclib::noSpeaker;
	this->parentID[1] = sclib::noSpeaker;
	this->bestFitID = sclib::noSpeaker;
	this->mergenceCounter = 0;

  deepCopy(pMergedModel, pSpeakerModels, pBackgroundModels, pSpeech, justLink);
	
	if (strcmp(speakerName, "") != 0) {
		this->speakerName = new char[strlen(speakerName) + 1]; 
		strcpy(this->speakerName, speakerName);
	} else {
		this->speakerName = NULL;
	}

	if (this->segmentCount > 0) {
		MArray_2D(this->segmentList, (long int)(this->segmentCount), 3, unsigned long int, "SC_Cluster.SC_Cluster: segmentList");
		for (unsigned long int x = 0; x < this->segmentCount; x++) {
			this->segmentList[x][0] = speakerID;
			this->segmentList[x][1] = speechStart[x];
			this->segmentList[x][2] = speechEnd[x];
		}
	} else {
		this->segmentList = NULL;
	}
}

//====================================================================================================================
//	constructor for a new cluster created by the mergence-function (it gets a segmentList insted of speech-start and 
//  -end)
//====================================================================================================================
SC_Cluster::SC_Cluster(SC_TweakableParameters* pTweak, SV_Data* pSpeech, unsigned long int** segmentList, SC_Model* pMergedModel, SC_Model* pBackgroundModels, SC_Model* pSpeakerModels, unsigned long int segmentCount, long int speakerID, const char* speakerName, bool justLink) {
	this->Next = NULL;
	this->speakerID = speakerID;
	this->segmentCount = segmentCount;
  this->segmentList = segmentList;
	this->pTweak = pTweak;
  this->classSig = 196275;
  this->IDconfidenceScore = 1.0; //in the beginning, the given speaker-id is for sure...
	this->ownerID = sclib::noSpeaker; //no owner until scoring...
	this->ownerIsCorrect = false; 
	this->parentID[0] = sclib::noSpeaker;
	this->parentID[1] = sclib::noSpeaker;
	this->bestFitID = sclib::noSpeaker;
	this->mergenceCounter = 0;

  deepCopy(pMergedModel, pSpeakerModels, pBackgroundModels, pSpeech, justLink);

	if (strcmp(speakerName, "") != 0) {
		this->speakerName = new char[strlen(speakerName) + 1]; 
		strcpy(this->speakerName, speakerName);
	} else {
		this->speakerName = NULL;
	}
}

//====================================================================================================================
//	copy-constructor
//  ATTENTION: by default (and in contrary to both other constructors), the justLink-option is set =true
//  if noWorkData==true && justLink==false, the speech-frames and models aren't disposed at all
//====================================================================================================================
SC_Cluster::SC_Cluster(const SC_Cluster& pParent, bool justLink, bool noWorkData) {
	this->Next = pParent.Next;
	this->speakerID = pParent.speakerID;
	this->segmentCount = pParent.segmentCount;
  this->pTweak = pParent.pTweak;
  this->classSig = 196275;
  this->IDconfidenceScore = pParent.IDconfidenceScore;
	this->ownerID = pParent.ownerID;
	this->ownerIsCorrect = pParent.ownerIsCorrect;
	this->parentID[0] = pParent.parentID[0];
	this->parentID[1] = pParent.parentID[1];
	this->bestFitID = pParent.bestFitID;
	this->mergenceCounter = pParent.mergenceCounter;

	MArray_2D(this->segmentList, (long int)(this->segmentCount), 3, unsigned long int, "SC_Cluster.SC_Cluster: segmentList");
  for (unsigned long int x = 0; x < this->segmentCount; x++) {
    this->segmentList[x][0] = pParent.segmentList[x][0];
    this->segmentList[x][1] = pParent.segmentList[x][1];
    this->segmentList[x][2] = pParent.segmentList[x][2];
  }	

	if (justLink==false && noWorkData==true) {
		this->pSpeakerModels = NULL;
		this->pBackgroundModels = NULL;
		this->pMergedModel = NULL;
		this->pSpeechFrames = NULL;
	} else {
		deepCopy(pParent.pMergedModel, pParent.pSpeakerModels, pParent.pBackgroundModels, pParent.pSpeechFrames, justLink);
	}

	if (pParent.speakerName != NULL && strcmp(pParent.speakerName, "") != 0) {
		this->speakerName = new char[strlen(pParent.speakerName) + 1]; 
		strcpy(this->speakerName, pParent.speakerName);
	} else {
		this->speakerName = NULL;
	}
}

//====================================================================================================================
//	creates copies (!!!) of the (linked linked lists of) parameters for the cluster's data instead of just linking it
//  if the justLink parameter is =true, copys of the speechframes are made (so they are new objects which can have 
//  their own next-pointer), but they internally just link to the parent's data matrix
//====================================================================================================================
void SC_Cluster::deepCopy(SC_Model *pMergedModel, SC_Model *pSpeakerModels, SC_Model *pBackgroundModels, SV_Data *pSpeechFrames, bool justLink) {
  //copy speech frames
	if (pSpeechFrames != NULL) {
		this->pSpeechFrames = deepCopySpeechFrames(pSpeechFrames, justLink);
	} else {
		this->pSpeechFrames = NULL;
	}

	//copy merged model
  if (pMergedModel != NULL) {
    this->pMergedModel = deepCopyModels(pMergedModel, 1);
  } else {
    this->pMergedModel = NULL;
  }

  //copy background-models
  if (pBackgroundModels != NULL) {
    this->pBackgroundModels = deepCopyModels(pBackgroundModels, this->segmentCount);
  } else {
    this->pBackgroundModels = NULL;
  }

  //copy speaker models
  if (pSpeakerModels != NULL) {
    this->pSpeakerModels = deepCopyModels(pSpeakerModels, this->segmentCount, justLink);
  } else {
    this->pSpeakerModels = NULL;
  }
  
  return;
}

//====================================================================================================================
//  deepCopy of linked lists of models (speaker/background/merged); 
//	needed by deepCopy and (that's the reason for the outsourced method) the mergence-function
//  if listCount>0, only listCount elements are copied
//====================================================================================================================
SC_Model* SC_Cluster::deepCopyModels(SC_Model* pModelList, unsigned int listCount, bool justLink) {
	unsigned int count = 0;
	SC_Model *pNewHook = NULL, *pOldHook = pModelList, *pRes = NULL;
  SC_ModelHandler modelHandler(this->pTweak);
  
  if (pModelList != NULL) {
    pNewHook = modelHandler.copyModel(pModelList, justLink);

    pRes = pNewHook;
    while (pNewHook->Next != NULL && (count < listCount || listCount == 0)) {
      pNewHook->Next = modelHandler.copyModel((SC_Model*)(pOldHook->Next), justLink);
      pOldHook = (SC_Model*)(pOldHook->Next);
      pNewHook = (SC_Model*)(pNewHook->Next);
			count++;
    }
  }

  return pRes;
}

//====================================================================================================================
//  deepCopy (or copy with pure data-linkage) of linked lists of speech-frames
//	needed by deepCopy and (that's the reason for the outsourced method) the mergence-function
//====================================================================================================================
SV_Data* SC_Cluster::deepCopySpeechFrames(SV_Data* pSpeechFrames, bool justLink) {
  SV_Data *pNewDataHook = new SV_Data(*pSpeechFrames, justLink), *pOldDataHook = pSpeechFrames, *pRes;
  
  //copy speech frames
  pRes = pNewDataHook;
  while (pOldDataHook->Next != NULL) {
    pNewDataHook->Next = new SV_Data(*pOldDataHook->Next, justLink);
    pOldDataHook = pOldDataHook->Next;
    pNewDataHook = pNewDataHook->Next;
  }  

  return pRes;
}

//====================================================================================================================
//	destructor
//  all the data (models, speech-fames, segment-list) get destructed here; if the speech-frame's data-matrix was just
//  linked to it's partent, the SV_Data destructor handles this by just killing the link and leaving the original data 
//  as it was (so always remember to keep a pointer to the original, not justLinked cluster-objects to free them at the
//  end!)
//====================================================================================================================
SC_Cluster::~SC_Cluster() {
  releaseWorkData();

  MFree_1D(this->speakerName);
	MFree_2D(this->segmentList);
  this->classSig = 0;
}

//====================================================================================================================
//	steps trough the linked list of clusters and destructs all features, models and clusters
//====================================================================================================================
void SC_Cluster::killAllLinkedData(SC_Cluster* pToKill) {
	SC_Cluster *pHook, *pLast;

	pHook = pToKill;
	if (pHook != NULL) {
		do {
			pHook->killLinkedBackgroundModels(pHook);
			pHook->killLinkedSpeakerModels(pHook);
			pHook->killLinkedFeatureVectors(pHook);

      pLast = pHook;
			pHook = pHook->Next;

			MFree_0D(pLast);
		} while (pHook != NULL);
	}

  return;
}

//====================================================================================================================
//	destructs all the basic data-structures: (background-)models and speechFrames
//  this saves much of memory when all computations are done and the object is just holded to store the results in a 
//  partition
//====================================================================================================================
void SC_Cluster::releaseWorkData(void) {
  killLinkedBackgroundModels(this);
	killLinkedSpeakerModels(this);
  killLinkedFeatureVectors(this);
  MFree_0D(this->pMergedModel);
 
  return;
}

//====================================================================================================================
//	steps trough the linked list of feature vectors and destructs them; handling of linked data-matrices is done
//  by the SV_Data-destructor
//====================================================================================================================
void SC_Cluster::killLinkedFeatureVectors(SC_Cluster* pToKill) {
  sclib::destructLinkedList(pToKill->pSpeechFrames);
  pToKill->pSpeechFrames = NULL;

	return;
}

//====================================================================================================================
//	steps trough the linked list of background models and destructs them
//====================================================================================================================
void SC_Cluster::killLinkedBackgroundModels(SC_Cluster* pToKill) {
	sclib::destructLinkedList(pToKill->pBackgroundModels);
  pToKill->pBackgroundModels = NULL;

	return;
}

//====================================================================================================================
//	steps trough the linked list of background models and destructs them
//====================================================================================================================
void SC_Cluster::killLinkedSpeakerModels(SC_Cluster* pToKill) {
	sclib::destructLinkedList(pToKill->pSpeakerModels);
  pToKill->pSpeakerModels = NULL;

	return;
}

//====================================================================================================================
//	merge 'this' with 'pSecond' to a new cluster
//  the models of the new clusters are standalone copies; the speech-frames are also new objects, but the data-
//  matrix just points to the original one, so don't delete the parents!
//====================================================================================================================
SC_Cluster* SC_Cluster::mergeClusters(SC_Cluster* pSecond, unsigned int mergeCounter) {
	unsigned long int** newSegmentList;
  unsigned long int   newSegmentCount = this->segmentCount + pSecond->segmentCount;
  SC_Cluster* pNewCluster;
  SC_Model* pCombinedModel;
	char* speakerName;
	SV_Data* pFeatureHook = sclib::getListWithIndex(this->pSpeechFrames, this->segmentCount-1);
	SC_Model* pBackgroundHook	= sclib::getListWithIndex(this->pBackgroundModels, this->segmentCount-1);
	SC_Model* pSpeechHook = sclib::getListWithIndex(this->pSpeakerModels, this->segmentCount-1);
  SC_ModelHandler *pHandler = new SC_ModelHandler(this->pTweak);
 
  //connect the models and speech-frames for the child
	if (pFeatureHook != NULL) {
		pFeatureHook->Next = pSecond->pSpeechFrames;
	} else {
		pFeatureHook = pSecond->pSpeechFrames;
	}
	if (pBackgroundHook != NULL) {
		pBackgroundHook->Next = pSecond->pBackgroundModels;
	} else {
		pBackgroundHook = pSecond->pBackgroundModels;
	}
	if (pSpeechHook != NULL) {
		pSpeechHook->Next = pSecond->pSpeakerModels;
	} else {
		pSpeechHook = pSecond->pSpeakerModels;
	}
  
	//create the merged model, i.e. the one that represents the whole cluster in one model
  switch (this->pTweak->cluster.mergeMode) {
		case sclib::mergeNone:
			pCombinedModel = NULL;
			break;

		case sclib::mergeAddUp:
      pCombinedModel = pHandler->combineModels(this->getMergedModel(), pSecond->getMergedModel());
      break;

		case sclib::mergeRetrain:
      pCombinedModel = pHandler->combineModels(this->getMergedModel(), pSecond->getMergedModel(), this->pSpeechFrames, this->segmentCount+pSecond->getSegmentCount(), this->pBackgroundModels);
      break;

		default:
      REPORT_ERROR(SVLIB_BadArg, "Specified model merge mode unknown");
      break;
  }

	MArray_2D(newSegmentList, (long int)(newSegmentCount), 3, unsigned long int, "SC_Cluster.mergeClusters: newSegmentList");
  for (unsigned long int x = 0; x < newSegmentCount; x++) {
    if (x < this->segmentCount) {
	    newSegmentList[x][0] = this->segmentList[x][0];
      newSegmentList[x][1] = this->segmentList[x][1];
      newSegmentList[x][2] = this->segmentList[x][2];
    } else {
	    newSegmentList[x][0] = pSecond->segmentList[x - this->segmentCount][0];
      newSegmentList[x][1] = pSecond->segmentList[x - this->segmentCount][1];
      newSegmentList[x][2] = pSecond->segmentList[x - this->segmentCount][2];
    }
  }

	if (mergeCounter > 0) {
		speakerName = new char[strlen(this->speakerName) + strlen(pSecond->speakerName) + ((mergeCounter / 10) + 1) + 6];
		sprintf(speakerName, "(%d: %s %s)", mergeCounter, this->speakerName, pSecond->speakerName);
	} else {
		speakerName = new char[strlen(this->speakerName) + strlen(pSecond->speakerName) + 4];
		sprintf(speakerName, "(%s %s)", this->speakerName, pSecond->speakerName);
	}

	pNewCluster = new SC_Cluster(this->pTweak, this->pSpeechFrames, newSegmentList, pCombinedModel, this->pBackgroundModels, this->pSpeakerModels, this->segmentCount + pSecond->segmentCount, this->speakerID, speakerName, true);
	pNewCluster->setParentID(0, this->getSpeakerID());
	pNewCluster->setParentID(1, pSecond->getSpeakerID());
	pNewCluster->setMergenceCount(this->mergenceCounter + pSecond->getMergenceCount() + 1);
	MFree_1D(speakerName);
  MFree_0D(pCombinedModel);
  
  //remove the connections from the parents
	if (pFeatureHook != NULL) {
		pFeatureHook->Next = NULL;
	}
	if (pBackgroundHook != NULL) {
		pBackgroundHook->Next = NULL;
	}
	if (pSpeechHook != NULL) {
		pSpeechHook->Next = NULL;
	}

  MFree_0D(pHandler);

	return pNewCluster;
}

//====================================================================================================================
//	Compute covariance matrix of all feature vectors of this cluster
//====================================================================================================================
double** SC_Cluster::getCovar(bool diagonal) {
  SV_Data *pCompleteData = this->pSpeechFrames->MergeData(this->segmentCount);
  SC_MatrixFunctions matFunc;
  double** cov;
	
	cov = matFunc.cov(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col, NULL, 0, 0, 0, 0, diagonal);
  MFree_0D(pCompleteData);

  return cov;
}

//====================================================================================================================
//	If this cluster contains trajectorized features that have replaced original features during training, they are
//  unwinded and the total feature count is returned in T
//====================================================================================================================
double** SC_Cluster::getCovar(int &T, int &D, bool diagonal, bool unwindTrajectories) {
  SV_Data *pCompleteData = this->pSpeechFrames->MergeData(this->segmentCount);
  SC_MatrixFunctions matFunc;
  double** cov;

	if (unwindTrajectories==true && this->pTweak->modelTime.replaceTrainingData==true && pCompleteData->Hdr.Signature[2]>0) { //we have trajectories here that have replaced original features... unwind them
		SC_FeatureHandler handler(this->pTweak, false);
		SV_Data *pUnwinded = handler.unwindTrajectories(pCompleteData, pCompleteData->Hdr.Signature[3], pCompleteData->Hdr.Signature[4]);
		//SV_Data *pNorm = handler.loadFeature(this->pTweak->modelTime.normalizationFile);
		//handler.unNormalize(pUnwinded, pNorm, true);
		//MFree_0D(pNorm);
		MFree_0D(pCompleteData);
		pCompleteData = pUnwinded;
	}
	
	T = pCompleteData->Row;
	D = pCompleteData->Col;
	cov = matFunc.cov(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col, NULL, 0, 0, 0, 0, diagonal);
  MFree_0D(pCompleteData);

  return cov;
}

//====================================================================================================================
//	Return the covariance matrix of a this cluster as if it was merged with pSecond
//====================================================================================================================
double** SC_Cluster::getCombinedCovar(SC_Cluster *pSecond, int &T, int &D, bool diagonal, bool unwindTrajectories) {
	SV_Data *pCompleteData, *pHook = this->getSpeechFrames(this->getSegmentCount()-1), *pSave = pHook->Next;
  SC_MatrixFunctions matFunc;
  double** cov;

	pHook->Next = pSecond->getSpeechFrames();
	pCompleteData = this->getSpeechFrames()->MergeData(this->getSegmentCount()+pSecond->getSegmentCount());
	pHook->Next = pSave;

	if (unwindTrajectories==true && this->pTweak->modelTime.replaceTrainingData==true && pCompleteData->Hdr.Signature[2]>0) { //we have trajectories here that have replaced original features... unwind them
		SC_FeatureHandler handler(this->pTweak, false);
		SV_Data *pUnwinded = handler.unwindTrajectories(pCompleteData, pCompleteData->Hdr.Signature[3], pCompleteData->Hdr.Signature[4]);
		//SV_Data *pNorm = handler.loadFeature(this->pTweak->modelTime.normalizationFile);
		//handler.unNormalize(pUnwinded, pNorm, true);
		//MFree_0D(pNorm);
		MFree_0D(pCompleteData);
		pCompleteData = pUnwinded;
	}

	T = pCompleteData->Row;
	D = pCompleteData->Col;
	cov = matFunc.cov(pCompleteData->Mat, pCompleteData->Row, pCompleteData->Col, NULL, 0, 0, 0, 0, diagonal);
	MFree_0D(pCompleteData);

  return cov;
}

//====================================================================================================================
//	Needed by sclib::destructLinkedList()
//====================================================================================================================
bool SC_Cluster::Valid(void) {
  return (this->classSig == 196275) ? true: false;
}

//====================================================================================================================
//	Get the Number of samples for the segment with index 'segment'; if 'segment' =-1, return the number of all
//  samples for the whole cluster
//====================================================================================================================
unsigned long int SC_Cluster::getSampleCount(long int segment) {
  unsigned long int count = 0;

  for (unsigned long int i = 0; i < this->segmentCount; i++) {
    if (segment == -1 || segment == i) {
      count += this->segmentList[i][2] - this->segmentList[i][1] + 1; //end - start
    }
  }

  return count;
}

//====================================================================================================================
//	Get the Number of feature vectors for the segment with index 'segment'; if 'segment' =-1, return the number of all
//  feature vectors for the whole cluster
//====================================================================================================================
unsigned long int SC_Cluster::getFeatureVectorCount(long int segment) {
  unsigned long int count = 0;
  SV_Data *pHook = this->pSpeechFrames;

  for (unsigned long int i = 0; i < this->segmentCount; i++) {
    if (segment == -1 || segment == i) {
      count += pHook->Row;
    }
    pHook = pHook->Next;
  }

  return count;
}

//====================================================================================================================
//	If this cluster contains segments in it's segment-list, which are direct neighbores (i.e. the end+1 of the first 
//  one is the beginning of the next one), they are united to one complete single segment covering both ranges (needed 
//  so that the segment-based scoring works properly).
//  It is a bit difficult to decide whether two segments ae direct neighbors, because not border-frames, but border-
//  samples are stored in the segmentList, and because the original frames could overlap, the start-sample of the 
//  direct successor might lie shortly before the end-sample of it's predecessor. We workaround this like that:
//  A segment is the direct successor of another one, if it's start lies between start and end+1 of the other one
//  The number of united segments is returned
//====================================================================================================================
unsigned long int SC_Cluster::uniteNeighboringSegments(void) {
  unsigned long int x, y, z, segmentEnd, count = 0;
  unsigned long int **newList; //Array of the form: "segmentNr | segmentStart | segmentStop"

  //unite all connected segments to a single big one, make the parts invalid by setting all components =0
  for (x = 0; x < this->segmentCount; x++) {
    segmentEnd = this->segmentList[x][2];
    
    for (y = 0; y < this->segmentCount; y++) {
      if (x != y && !(this->segmentList[y][0] == 0 && this->segmentList[y][1] == 0 && this->segmentList[y][2] == 0)) { //don't look at segments previously made invalid; don't compare equal segments
        if (this->segmentList[y][1] >= this->segmentList[x][1] && this->segmentList[y][1] <= segmentEnd+1) {
          segmentEnd = this->segmentList[y][2];
          for (z = 0; z < 3; z++) {
            this->segmentList[y][z] = 0; //make this list-entry invalid
          }
          count++; //remember the change
        }
      }
    }
    
    this->segmentList[x][2] = segmentEnd;
  }

  //delete the entries in the list which have been made invalid by copying only the not-invalid entries into the new list
  if (count > 0) {
    MArray_2D(newList, (long int)(this->segmentCount - count), 3, unsigned long int, "SC_Cluster.uniteNeighboringSegments: newList");
    
    x = 0;
    for (y = 0; y < this->segmentCount; y++) {
      if (!(this->segmentList[y][0] == 0 && this->segmentList[y][1] == 0 && this->segmentList[y][2] == 0)) {
        for (z = 0; z < 3; z++) {
          newList[x][z] = this->segmentList[y][z];
        }
        x++;
      }
    }

    MFree_2D(this->segmentList);
    this->segmentList = newList;
    this->segmentCount -= count;
  }

  return count;
}

//====================================================================================================================
//	Set new speaker models
//====================================================================================================================
void SC_Cluster::setSpeakerModels(SC_Model* pNewModels) {
  killLinkedSpeakerModels(this);
  this->pSpeakerModels = deepCopyModels(pNewModels, this->segmentCount);

  return;
}

//====================================================================================================================
//	Set a new merged models
//====================================================================================================================
void SC_Cluster::setMergedModel(SC_Model* pNewModel) {
  MFree_0D(this->pMergedModel);
  this->pSpeakerModels = deepCopyModels(pNewModel, 1);

  return;
}

//====================================================================================================================
//	Set a new (linked list) background models
//====================================================================================================================
void SC_Cluster::setBackgroundModel(SC_Model* pNewModels) {
  killLinkedBackgroundModels(this);
  this->pBackgroundModels = deepCopyModels(pNewModels, this->segmentCount);

  return;
}

//====================================================================================================================
//  Assumes that the segments (i.e. the linked-list structure) remain the same or are changed accordingly before;
//  the downside is that afterwards the pNewSpeechFrames pointer is destroyed!
//====================================================================================================================
void SC_Cluster::setSpeechFrames(SV_Data *&pNewSpeechFrames) {
	SV_Data *pLastSegmentsNext;
	
	if (this->pSpeechFrames == NULL) {
		this->pSpeechFrames = pNewSpeechFrames;
	} else {
		if (this->segmentCount > 0) {
			pLastSegmentsNext = sclib::getListWithIndex(this->pSpeechFrames, this->segmentCount-1)->Next;
			
			//avoid killing the old pointer so that possibly existing backward link structures are preserved; instead, fill it with the new data without copying much
			if (this->pSpeechFrames->getJustLinked() == false) {
				MFree_2D(this->pSpeechFrames->Mat);
			} else {
				this->pSpeechFrames->Mat = NULL;
			}
			this->pSpeechFrames->Hdr = pNewSpeechFrames->Hdr;
			this->pSpeechFrames->Row = pNewSpeechFrames->Row;
			this->pSpeechFrames->Col = pNewSpeechFrames->Col;
			this->pSpeechFrames->Mat = pNewSpeechFrames->Mat;
			this->pSpeechFrames->setJustLinked(false); //yes, false!

			//now the data of the first element is at its place, care for the rest of the list
			sclib::destructLinkedList(this->pSpeechFrames->Next, sclib::max(0, this->segmentCount-2)); //ok, old stuff is dead now
			this->pSpeechFrames->Next = pNewSpeechFrames->Next; //ok, new stuff is plugged in now
			sclib::getListWithIndex(this->pSpeechFrames, this->segmentCount-1)->Next = pLastSegmentsNext; //ok, the forward connection is also reconstructed now

			//handle garbage of the first pNewSpeechFrames-element, which has been dissembled and merged into the first pSpeechFrames element 
			pNewSpeechFrames->setJustLinked(true); //the matrix now belongs to this->pSpeechFrames
			MFree_0D(pNewSpeechFrames);
			pNewSpeechFrames = this->pSpeechFrames;
		} else {
			this->pSpeechFrames = pNewSpeechFrames; //this is a situation we don't want because some things are inconsistent and undefined now, but we do the most sane thing
		}
	}

	return;
}

//====================================================================================================================
//	retrain the merged-model of this cluster using all available features; explicitly make it a model of the desired 
//  type (this way, this method can also be used to convert/switch between different model types)
//  if modelOrder=0, the best order is guessed by the modelHandler
//====================================================================================================================
void SC_Cluster::retrainMergedModel(unsigned long int modelType, unsigned short int modelOrder) {
  SC_ModelHandler *pHandler = new SC_ModelHandler(this->pTweak);
  unsigned short int bestOrder = (modelOrder > 0) ? modelOrder : pHandler->guessModelOrder(this->pSpeechFrames, this->pBackgroundModels, modelType, 1, this->pTweak->modelHandler.maxSpeakerModelOrder, this->segmentCount);
  SC_Model *pNewModel = pHandler->buildModel(this->pSpeechFrames, this->pBackgroundModels, bestOrder, modelType, this->segmentCount);
  
  if (pNewModel != NULL) {
    MFree_0D(this->pMergedModel);
    this->pMergedModel = pNewModel;
  }

  MFree_0D(pHandler);

  return;
}

//====================================================================================================================
// (re-)create the merged model, i.e. the single model representing all segments in this cluster, using the method 
// stated in the tweakable parameters (i.e. addup, rebuild or do nothing)
//====================================================================================================================
void SC_Cluster::reBuildMergedModel(void) {
	SC_ModelHandler modeller(this->pTweak, true);
	SC_Model *pHook, *pMerged = NULL, *pTmp;
	
	switch (this->pTweak->cluster.mergeMode) {
		case sclib::mergeNone: 
			MFree_0D(this->pMergedModel);
			break;

		case sclib::mergeRetrain:
			MFree_0D(this->pMergedModel);
			this->pMergedModel = modeller.buildModel(this->pSpeechFrames, this->pBackgroundModels, sclib::modeForeground, this->segmentCount);
			break;

		case sclib::mergeAddUp: //subsequently merge 2 models by adding up
			MFree_0D(this->pMergedModel);
			pMerged = this->pSpeakerModels;
			pHook = (this->pSpeakerModels != NULL) ? (SC_Model*)(this->pSpeakerModels->Next) : NULL;
			while (pHook != NULL) { //subsequently merge the next model in the list with the already comvined models by addign components up
				pTmp = pMerged->combineModels(pHook);
				if (pMerged != this->pSpeakerModels) {
					MFree_0D(pMerged);
				}
				pMerged = pTmp;
				pHook = (SC_Model*)(pHook->Next);
			}
			if (pMerged != this->pSpeakerModels) { //if there was more than one model in the list
				this->pMergedModel = pMerged;
			}
			break;

		default:
			REPORT_ERROR(SVLIB_BadArg, "given cluster-mergemode is unknown or unsupported");
			break;
	}

	return;
}

//====================================================================================================================
//	add a segment to this cluster for which no model exists (so merging can't be applied); if wished, the cluster's
//  merged-model gets retrained; if segmentId==0, segmentCount-1 is assumed as the id
//====================================================================================================================
void SC_Cluster::addSegment(SV_Data *pSpeech, unsigned long int speechStart, unsigned long int speechEnd, unsigned long int segmentId, bool justLink, bool retrainModel) {
  SV_Data *pHook; 
  unsigned long int **newSegmentList = NULL;
  
  //push the new speech frames in the linked list
	pHook = sclib::getLastInList(this->pSpeechFrames);
  pHook->Next = deepCopySpeechFrames(pSpeech, justLink);

  //increase the count of the segments of this cluster
  this->segmentCount++;

  //add the new segment-borders to the segment-list
   MArray_2D(newSegmentList, (long int)(this->segmentCount), 3, unsigned long int, "SC_Cluster.addSegment: newSegmentList");
  for (unsigned long int x = 0; x < this->segmentCount-1; x++) {
    newSegmentList[x][0] = this->segmentList[x][0]; //x;
    newSegmentList[x][1] = this->segmentList[x][1];
    newSegmentList[x][2] = this->segmentList[x][2];
  }
	newSegmentList[this->segmentCount-1][0] = (segmentId > 0) ? segmentId : this->segmentCount-1;
  newSegmentList[this->segmentCount-1][1] = speechStart;
  newSegmentList[this->segmentCount-1][2] = speechEnd;
  MFree_2D(this->segmentList);
  this->segmentList = newSegmentList;

  //retrain speaker-model, if wished
  if (retrainModel == true) {
    retrainMergedModel(this->pMergedModel->Hdr.ModelType, 0);
  }

  return;
}

//====================================================================================================================
//	return the linked list of features beginning with the specified entry of the list; return NULL if the list is 
//  shorter than specified
//====================================================================================================================
SV_Data* SC_Cluster::getSpeechFrames(long int segment) {
  SV_Data *pHook = this->pSpeechFrames;
  long int count = 0;
  
  while (segment < count && pHook->Next != NULL) {
    pHook = pHook->Next;
    count++;
  }

  if (segment < count) {
    pHook = NULL;
  }

  return pHook;
}

//====================================================================================================================
//  Store a copy of the given name as the new speaker name
//====================================================================================================================
void SC_Cluster::setSpeakerName(char* newName) {
	MFree_1D(this->speakerName); 

	MArray_1D(this->speakerName, strlen(newName)+1, char, "SC_Cluster.setSpeakerName: this->speakerName");
	sprintf(this->speakerName, "%s\0", newName);

	return;
}

//====================================================================================================================
//  Return a new list of |segStart|segEnd| tuples, ordered by segStart ascending; memory has to be released by the 
//  caller
//====================================================================================================================
unsigned long int** SC_Cluster::getOrderedSegments(void) {
	unsigned long int **orderedSegmentList, y;

	MArray_2D(orderedSegmentList, (long int)(this->segmentCount), 2, unsigned long int, "SC_Cluster.getOrderedSegments: orderedSegmentList");

	for (y = 0; y < this->segmentCount; y++) {
		orderedSegmentList[y][0] = this->segmentList[y][1];
		orderedSegmentList[y][1] = this->segmentList[y][2];
	}

	sclib::quickSort(orderedSegmentList, 0, this->segmentCount-1, 2, 0);

	return orderedSegmentList;
}

//====================================================================================================================
//  File-I/O for this class, so that the sate of a cluster-object can be saved to a file
//  Linked data gets saved with this file and is then non-linked after loading, so memory-requirements may rise...
//  Everything (also: linked list of clusters) gets saved except the pTweak member
//====================================================================================================================
bool SC_Cluster::save(const char *fileName, unsigned int listNr) {
	bool res = true;
	char ext[sclib::bufferSize], *fName;
	unsigned long int len, x, result;
	fstream clFile;
	SC_FeatureHandler featureHandler(this->pTweak);
	SC_Model *pModelHook;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes;
	io.getCurrentDatatypeSizes(codeSizes);

	sprintf(ext, ".cluster_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	clFile.open(fName, ios::out|ios::binary);  //truncate
	MFree_1D(fName);

	//machine header (to achieve machine independance)
	io.writeMachineHeader(&clFile, codeSizes);

	//simple types
	io.writeScalar(&clFile, this->classSig);
	io.writeScalar(&clFile, this->IDconfidenceScore);
	io.writeScalar(&clFile, this->speakerID);
	io.writeScalar(&clFile, this->ownerID);
	io.writeScalar(&clFile, this->ownerIsCorrect);
	io.writeArray(&clFile, this->parentID, 2);
	io.writeScalar(&clFile, this->bestFitID);
	io.writeScalar(&clFile, this->segmentCount);
	io.writeMatrix(&clFile, this->segmentList, this->segmentCount, 3);
	io.writeScalar(&clFile, this->mergenceCounter);
	
	//speaker name
	if (this->speakerName != NULL) {
		len = (unsigned long int)(strlen(this->speakerName));
		io.writeScalar(&clFile, len);
		io.writeArray(&clFile, this->speakerName, len);

	} else {
		len = 0;
		io.writeScalar(&clFile, len);
	}

	//speech features
	sprintf(ext, ".speechFrames_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	if (!featureHandler.saveFeatureList(fName, this->pSpeechFrames)) {
		REPORT_ERROR(SVLIB_Fail, "Saving SC_Cluster failed due to speech frames!");
		res = false;
	}
	MFree_1D(fName);

	//speaker model(s)
	pModelHook = this->pSpeakerModels;
	sprintf(ext, ".spkModels_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	len = sclib::getListCount(this->pSpeakerModels); //write to the cluster's file the nr of speaker-models in the linked list
	io.writeScalar(&clFile, len);
	for (x = 0; x < len; x++) {
		io.writeScalar(&clFile, pModelHook->Hdr.ModelType);
		if (x == 0) {
			pModelHook->OpenFile(fName, WRITE_MODEL);
		} else {
			pModelHook->OpenFile(fName, APPEND_MODEL);
		}
		result = pModelHook->SaveModel();
		pModelHook->CloseFile();
		if (result <= 0) {
			REPORT_ERROR(SVLIB_Fail, "Saving SC_Cluster failed due to speaker models!");
			res = false;
		}
		pModelHook = (SC_Model*)(pModelHook->Next);
	}
	MFree_1D(fName);


	//merged model
	if (this->pMergedModel == NULL) {
		io.writeScalar(&clFile, 0); //there is no model with type 0, so indicate this way that there has been none to save
	} else {
		io.writeScalar(&clFile, this->pMergedModel->Hdr.ModelType);
		sprintf(ext, ".mergedModel_%d\0", listNr);
		fName = sclib::exchangeFileExtension(fileName, ext);
		this->pMergedModel->OpenFile(fName, WRITE_MODEL);
		MFree_1D(fName);
		result = this->pMergedModel->SaveModel();
		this->pMergedModel->CloseFile();
		if (result <= 0) {
			REPORT_ERROR(SVLIB_Fail, "Saving SC_Cluster failed due to merged model!");
			res = false;
		}
	}

	//background model(s)
	pModelHook = this->pBackgroundModels;
	sprintf(ext, ".bgModels_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	len = sclib::getListCount(this->pBackgroundModels); //write to the cluster's file the nr of bg-models in the linked list
	io.writeScalar(&clFile, len);
	for (x = 0; x < len; x++) {
		io.writeScalar(&clFile, pModelHook->Hdr.ModelType);
		if (x == 0) {
			pModelHook->OpenFile(fName, WRITE_MODEL);
		} else {
			pModelHook->OpenFile(fName, APPEND_MODEL);
		}
		result = pModelHook->SaveModel();
		pModelHook->CloseFile();
		if (result <= 0) {
			REPORT_ERROR(SVLIB_Fail, "Saving SC_Cluster failed due to background models!");
			res = false;
		}
		pModelHook = (SC_Model*)(pModelHook->Next);
	}
	MFree_1D(fName);

	//care for the rest of the clusters in the linked list; write as the last entry in this file the nr of following clusters in this list
	len = sclib::getListCount(this);
	io.writeScalar(&clFile, len);
  if (clFile.good() != TRUE) { //close file here 'cause the next cluster is written to another file, so we don't need it open anymore
		REPORT_ERROR(SVLIB_Fail, "Saving SC_Cluster failed!");
		res = false;
	}
	clFile.close();
	if (len > 1) {
		if (!this->Next->save(fileName, listNr+1)) {
			REPORT_ERROR(SVLIB_Fail, "Saving linked list of SC_Cluster failed!");
			res = false;
		}
	}

	return res;
}

//====================================================================================================================
//  File-I/O for this class, so that the sate of a cluster-object can be loaded from a file
//  Linked data gets saved with this file and is then non-linked after loading, so memory-requirements may rise...
//  Everything (also: linked list of clusters) gets loaded except the pTweak member
//====================================================================================================================
SC_Cluster* SC_Cluster::load(const char *fileName, unsigned int listNr) {
	unsigned long int x, len;
	long int modelType, lastPos = 0;
	char ext[sclib::bufferSize], *fName;
	int bytes;
	SC_Cluster *pCluster = this;
	SC_Model *pModelHook = NULL, *pNewModel;
	SC_FeatureHandler featureHandler(this->pTweak);
	SC_ModelHandler modelHandler(this->pTweak, false);
	fstream clFile;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

	sprintf(ext, ".cluster_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	clFile.open(fName, ios::in|ios::binary);  //read
	if (clFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Loading SC_Cluster failed!");
		pCluster = NULL;		
	}
	MFree_1D(fName);

	//machine header (to achieve machine independance)
	bytes = io.readMachineHeader(&clFile, fileSizes, true);
	if (bytes > 0) {
		io.consumeBytes(&clFile, bytes);
	} else {
		bytes = 0;
	}

	//simple types
	io.readScalar(&clFile, this->classSig, codeSizes, fileSizes);
	io.readScalar(&clFile, this->IDconfidenceScore, codeSizes, fileSizes);
	io.readScalar(&clFile, this->speakerID, codeSizes, fileSizes);
	io.readScalar(&clFile, this->ownerID, codeSizes, fileSizes);
	io.readScalar(&clFile, this->ownerIsCorrect, codeSizes, fileSizes);
	io.readArray(&clFile, this->parentID, 2, codeSizes, fileSizes);
	io.readScalar(&clFile, this->bestFitID, codeSizes, fileSizes);
	io.readScalar(&clFile, this->segmentCount, codeSizes, fileSizes);
	MFree_2D(this->segmentList);
	MArray_2D(this->segmentList, (long int)(this->segmentCount), 3, unsigned long int, "SC_Cluster.load: segmentList");
	io.readMatrix(&clFile, this->segmentList, this->segmentCount, 3, codeSizes, fileSizes);
	io.readScalar(&clFile, this->mergenceCounter, codeSizes, fileSizes);

	//speaker name
	MFree_1D(this->speakerName);
	io.readScalar(&clFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->speakerName, len+1, char, "SC_Cluster.load: speakerName");
		io.readArray(&clFile, this->speakerName, len, codeSizes, fileSizes);
		this->speakerName[len] = '\0';
	}

	//speech features
	sclib::destructLinkedList(this->pSpeechFrames);
	sprintf(ext, ".speechFrames_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	this->pSpeechFrames = featureHandler.loadFeature(fName);
	MFree_1D(fName);
	if (this->pSpeechFrames == NULL) {
		REPORT_ERROR(SVLIB_Fail, "Loading SC_Cluster failed while loading speech frames!");
		pCluster = NULL;
	}

	//speaker model(s)
	sclib::destructLinkedList(this->pSpeakerModels);
	io.readScalar(&clFile, len, codeSizes, fileSizes);
	sprintf(ext, ".spkModels_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	lastPos = 0;
	for (x = 0; x < len; x++) {
		io.readScalar(&clFile, modelType, codeSizes, fileSizes);
		pNewModel = modelHandler.createRawModel(modelType);
		pNewModel->OpenFile(fName, READ_MODEL);
		pNewModel->getDFile()->seekg(lastPos);
		if (pNewModel->LoadModel() == NULL) {
			REPORT_ERROR(SVLIB_Fail, "Loading SC_Cluster failed while loading speaker models!");
			pCluster = NULL;
		}
		lastPos = pNewModel->getDFile()->tellg();
		pNewModel->CloseFile();
		if (x == 0) {
			this->pSpeakerModels = pNewModel;
			pModelHook = this->pSpeakerModels;
		} else {
			pModelHook->Next = pNewModel;
			pModelHook = (SC_Model*)(pModelHook->Next);
		}
	}
	MFree_1D(fName);

	//merged model
	MFree_0D(this->pMergedModel);
	io.readScalar(&clFile, modelType, codeSizes, fileSizes);
	if (modelType > 0) {
		sprintf(ext, ".mergedModel_%d\0", listNr);
		fName = sclib::exchangeFileExtension(fileName, ext);
		this->pMergedModel = modelHandler.createRawModel(modelType);
		this->pMergedModel->OpenFile(fName, READ_MODEL);
		MFree_1D(fName);
		if (this->pMergedModel->LoadModel() == NULL) {
			REPORT_ERROR(SVLIB_Fail, "Loading SC_Cluster failed while loading merged model!");
			pCluster = NULL;
		}
		this->pMergedModel->CloseFile();
	} else {
		this->pMergedModel = NULL;
	}

	//background model(s)
	sclib::destructLinkedList(this->pBackgroundModels);
	io.readScalar(&clFile, len, codeSizes, fileSizes);
	sprintf(ext, ".bgModels_%d\0", listNr);
	fName = sclib::exchangeFileExtension(fileName, ext);
	lastPos = 0;
	for (x = 0; x < len; x++) {
		io.readScalar(&clFile, modelType, codeSizes, fileSizes);
		pNewModel = modelHandler.createRawModel(modelType);
		pNewModel->OpenFile(fName, READ_MODEL);
		pNewModel->getDFile()->seekg(lastPos);
		if (pNewModel->LoadModel() == NULL) {
			REPORT_ERROR(SVLIB_Fail, "Loading SC_Cluster failed while loading background models!");
			pCluster = NULL;
		}
		lastPos = pNewModel->getDFile()->tellg();
		pNewModel->CloseFile();
		if (x == 0) {
			this->pBackgroundModels = pNewModel;
			pModelHook = this->pBackgroundModels;
		} else {
			pModelHook->Next = pNewModel;
			pModelHook = (SC_Model*)(pModelHook->Next);
		}
	}
	MFree_1D(fName);

	//care for the rest of the clusters in the linked list
	io.readScalar(&clFile, len, codeSizes, fileSizes);
	if (clFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Loading SC_Cluster failed!");
		pCluster = NULL;
	}
	clFile.close();
	if (len > 1) {
		this->Next = new SC_Cluster(this->pTweak);
		if (!this->Next->load(fileName, listNr+1)) {
			REPORT_ERROR(SVLIB_Fail, "Loading linked list of SC_Cluster failed!");
			pCluster = NULL;
		}
	}

	return pCluster;
}

//====================================================================================================================
// Dump cluster's content in ASCII 
//====================================================================================================================
ostream& operator<< (ostream& OutS, SC_Cluster& cluster) {
	unsigned long int x;
	SV_Data *pSpeechHook = cluster.pSpeechFrames;
	SC_Model *pModelHook = cluster.pBackgroundModels;

	OutS.setf(ios::fixed|ios::basefield);
  OutS.precision(5);

	OutS << "Name and ID: " << cluster.speakerName << ", " << cluster.speakerID << endl;
	OutS << "Segment Count: " << cluster.segmentCount << endl;
	OutS << "Segments: Nr | start-sample | end-sample" << endl;
	for (x = 0; x < cluster.segmentCount; x++) {
		OutS << "\t" << cluster.segmentList[x][0] << "\t" << cluster.segmentList[x][1] << "\t" << cluster.segmentList[x][2] << endl;
	}
	OutS << "ID Confidence Score: " << cluster.IDconfidenceScore << endl;
	OutS << "Owner ID: " << cluster.ownerID << endl;
	OutS << "Owner is correct: " << cluster.ownerIsCorrect << endl;
	OutS << "Parent IDs: " << cluster.parentID[0] << ", " << cluster.parentID[1] << endl;
	OutS << "Best-Fit ID: " << cluster.bestFitID << endl;
	OutS << "Mergence Counter: " << cluster.mergenceCounter << endl;
	OutS << "Speech Frames:" << endl;
	while (pSpeechHook != NULL) {
		OutS << *pSpeechHook;
		pSpeechHook = pSpeechHook->Next;
	}
	OutS << "Speaker Models:" << endl;
	pModelHook = cluster.pSpeakerModels;
	while (pModelHook != NULL) {
		OutS << *pModelHook;
		pModelHook = (SC_Model*)(pModelHook->Next);

	}
	OutS << "Merged Model:" << endl;
	if (cluster.pMergedModel != NULL) {
		OutS << *cluster.pMergedModel;
	}
	OutS << "Background Models:" << endl;
	while (pModelHook != NULL) {
		OutS << *pModelHook;
		pModelHook = (SC_Model*)(pModelHook->Next);
	}

	return(OutS);
}
