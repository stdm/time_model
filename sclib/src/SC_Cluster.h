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

#ifndef __SC_Cluster_H__
#define __SC_Cluster_H__

#include "SC_Api.h"
#include "SC_Model.h"
#include "SC_TweakableParameters.h"
#include <SV_Data.h>
#include <SV_Error.h>

class SCLIB_API SC_Cluster {

	private :

	protected:

		//====================================================================================================================
		//	creates copies (!!!) of the (linked linked lists of) parameters for the cluster's data instead of just linking it
    //  if the justLink parameter is =true, copys of the speechframes are made (so they are new objects which can have 
    //  their own next-pointer), but they internally just link to the parent's data matrix
		//====================================================================================================================
    void deepCopy(SC_Model *pMergedModel, SC_Model *pSpeakerModel, SC_Model *pBackgroundModel, SV_Data *pSpeechFrames, bool justLink = false);
    SC_Model* deepCopyModels(SC_Model* pModelList, unsigned int listCount = 0, bool justLink = false);
    SV_Data* deepCopySpeechFrames(SV_Data* pSpeechFrames, bool justLink);

		//====================================================================================================================
		//	destructs the linkedlist of background-models/speech-frames in one single cluster
		//====================================================================================================================
		void killLinkedFeatureVectors(SC_Cluster* pToKill);
		void killLinkedBackgroundModels(SC_Cluster* pToKill);
		void killLinkedSpeakerModels(SC_Cluster* pToKill);

		//====================================================================================================================
		//	the parameters which build up a cluster representing a speaker
		//====================================================================================================================
		char *speakerName;				         //to name a speaker (as hypothesized)
		long int speakerID;			           //to identify a speaker (as hypothesized)
		unsigned long int segmentCount;		 //how many (hopefully) speaker-homogenious segments contributed to this cluster?
		unsigned long int **segmentList;	 //array of the form: "segmentNr | segmentStart | segmentStop"
		SV_Data *pSpeechFrames;			       //the original speech-features corresponding to each segment, formed as a linked list
		SC_Model *pBackgroundModels;	     //the corresponding background-models for each segment, formed as a linked list
		SC_Model *pSpeakerModels;					 //the original speaker models corresponding with each segment as a linked list
		SC_Model *pMergedModel;						 //the merged model for the whole cluster; NULL if merging is switched off
    double IDconfidenceScore;          //the (normalized) (log) likelihood (ratio) that this cluster belongs to the claimed speaker (set by the speaker-identificator algorithm)
		long int ownerID;									 //set by the speaker-clustering-scoring class: gt-speaker-id of the speaker that uttered the most samples of the segments assembled in this cluster; ==sclib::noSpeaker if not set
		bool ownerIsCorrect;							 //set by the speaker-clustering-scoring class: true, if this is the biggest cluster of that owner (otherwise it is just fitting)
		long int parentID[2];							 //speaker-id's (hypothesized) of this cluster's two parents, if it originates from merging; otherwise it is ==sclib::noSpeaker
		long int bestFitID;								 //set by the speaker-id-class to remember which was the cluster this one was identified as belonging to; ==sclib::noSpeaker if this cluster is considered an impostor or if no id took place
		unsigned long int mergenceCounter; //tells of how many mergences this cluster was the result of

		SC_TweakableParameters *pTweak;		//for tweakable parameters
    long classSig;                    //needed by Valid()-method

		void setMergenceCount(unsigned long int count) {this->mergenceCounter = count; return;}

	public :

		//====================================================================================================================
		//	constructors
    //  ATTENTION: by default (and in contrary to both other constructors), the justLink-option of the copy-constructor 
    //  is set ==true
		//  if noWorkData==true && justLink==false, the speech-frames and models aren't disposed at all
  	//====================================================================================================================
	  SC_Cluster(SC_TweakableParameters* pTweak, SV_Data* pSpeech = NULL, unsigned long int* speechStart = 0, unsigned long int* speechEnd = 0, SC_Model* pMergedModel = NULL, SC_Model* pBackgroundModels = NULL, SC_Model* pSpeakerModels = NULL, unsigned long int segmentCount = 1, long int speakerID = sclib::noSpeaker, const char* speakerName = "", bool justLink = false);
		SC_Cluster(SC_TweakableParameters* pTweak, SV_Data* pSpeech, unsigned long int** segmentList, SC_Model* pMergedModel, SC_Model* pBackgroundModels, SC_Model* pSpeakerModels, unsigned long int segmentCount = 1, long int speakerID = sclib::noSpeaker, const char* speakerName = "", bool justLink = false);
    SC_Cluster(const SC_Cluster& pParent, bool justLink = true, bool noWorkData = false);

		//====================================================================================================================
		//	destructor
    //  all the data (models, speech-fames, segment-list) get destructed here; if the speech-frame's data-matrix was just
    //  linked to it's partent, the SV_Data destructor handles this by just killing the link and leaving the original data 
    //  as it was (so always remember to keep a pointer to the original, not justLinked cluster-objects to free them at the
    //  end!)
		//====================================================================================================================
		virtual ~SC_Cluster();
    
    bool Valid(void);
     
    //====================================================================================================================
		//	destructs all the basic data-structures: (background-)models and speechFrames
    //  this saves much of memory when all computations are done and the object is just holded to store the results in a 
    //  partition
		//====================================================================================================================
		void releaseWorkData(void);

    //====================================================================================================================
		//	steps trough a linked list of clusters and destructs them properly
		//====================================================================================================================
		void killAllLinkedData(SC_Cluster* pToKill);

		//====================================================================================================================
		//	give access to protected members
		//====================================================================================================================
		unsigned long int	getSegmentCount(void) {return this->segmentCount;}
		unsigned long int getMergenceCount(void) {return this->mergenceCounter;}
		unsigned long int** getSegmentList(void) {return this->segmentList;}
		unsigned long int* getSegmentList(unsigned long int index) {return this->segmentList[index];}
		unsigned long int getSegmentList(unsigned long int index, unsigned long int field) {return this->segmentList[index][field];}
		long int getSpeakerID(void)	{return this->speakerID;}
		char* getSpeakerName(void) {return this->speakerName;}    
		SC_Model* getSpeakerModels(void) {return this->pSpeakerModels;}
		SC_Model* getMergedModel(void) {return (this->pMergedModel==NULL && this->segmentCount==1) ? this->pSpeakerModels : this->pMergedModel;} //if there's only one segment, return the speaker model list, which has only one element naturally describing the whole cluster
		SC_Model* getBackgroundModels(unsigned long int idx = 0) {return sclib::getListWithIndex(this->pBackgroundModels, idx);}
    SV_Data* getSpeechFrames(long int segment = 0);
    unsigned int getFeatureDim(void) {return (unsigned int)(this->pSpeechFrames->Col);}
    unsigned long int getSampleCount(long int segment = -1);
    unsigned long int getFeatureVectorCount(long int segment = -1);
    double getIDconfidenceScore(void) {return this->IDconfidenceScore;}
		long int getOwnerID(void) {return this->ownerID;}
		bool getOwnerIsCorrect(void) {return this->ownerIsCorrect;}
		long int getParentID(int idx) {return (idx >= 0 && idx < 2) ? this->parentID[idx] : sclib::noSpeaker;}
		long int getBestFitID(void) {return this->bestFitID;}

		void setSpeakerName(char* newName); //creates a copy of the given name
		void setSpeakerID(long int newID) {this->speakerID = newID; return;}
    void setSpeakerModels(SC_Model* pNewModels);
		void setMergedModel(SC_Model* pNewModel);
    void setBackgroundModel(SC_Model* pNewModels);
    void setIDconfidenceScore(double score) {this->IDconfidenceScore = score; return;}
		void setOwnerID(long int speakerId) {this->ownerID = speakerId; return;}
		void setOwnerIsCorrect(bool correct) {this->ownerIsCorrect = correct; return;}
		void setParentID(int idx, long int id) {if (idx >= 0 && idx < 2) {this->parentID[idx] = id;} return;}
		void setBestFitID(long int id) {this->bestFitID = id; return;}
		void setSpeechFrames(SV_Data *&pNewSpeechFrames); //assumes that the segments (i.e. the linked-list structure) remain the same or are changed accordingly before; the downside is that afterwards the pNewSpeechFrames pointer is destroyed!

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
    unsigned long int uniteNeighboringSegments(void);

		//====================================================================================================================
		//	Compute covariance matrix of all feature vectors of this cluster
		//====================================================================================================================
    double** getCovar(bool diagonal = false);
    double** getCovar(int &T, int &D, bool diagonal = false, bool unwindTrajectories = true);

		//====================================================================================================================
		//	Return the covariance matrix of a this cluster as if it was merged with pSecond
		//====================================================================================================================
		double** getCombinedCovar(SC_Cluster *pSecond, int &T, int &D, bool diagonal = false, bool unwindTrajectories = true);

		//====================================================================================================================
		//	To form a linked list of clusters
		//====================================================================================================================
		SC_Cluster* Next;

		//====================================================================================================================
		//	Merge 'this' with 'pSecond' to a new cluster
		//  the models of the new clusters are standalone copies; the speech-frames are also new objects, but the data-
		//  matrix just points to the original one, so don't delete the parents!
		//====================================================================================================================
		SC_Cluster* mergeClusters(SC_Cluster* pSecond, unsigned int mergeCounter = 0);

		//====================================================================================================================
		//	Retrain the merged-model of this cluster using all available features; explicitly make it a model of the desired 
    //  type (this way, this method can also be used to convert/switch between different model types
    //  if modelOrder=0, the best order is guessed by the modelHandler
		//====================================================================================================================
    void retrainMergedModel(unsigned long int modelType, unsigned short int modelOrder = 0);

	  //====================================================================================================================
	  // (Re-)create the merged model, i.e. the single model representing all segments in this cluster, using the method 
		// stated in the tweakable parameters (i.e. addup, rebuild or do nothing)
	  //====================================================================================================================
		void reBuildMergedModel(void);

		//====================================================================================================================
		//	Add a segment to this cluster for which no model exists (so meerging can't be applied); if wished, the cluster's
    //  speaker-model gets retrained; if segmentId==0, segmentCount-1 is assumed as the id
		//====================================================================================================================
    void addSegment(SV_Data *pSpeech, unsigned long int speechStart, unsigned long int speechEnd, unsigned long int segmentId = 0, bool justLink = false, bool retrainModel = false);

		//====================================================================================================================
		//  Return a new list of |segStart|segEnd| tuples, ordered by segStart ascending; memory has to be released by the 
		//  caller
		//====================================================================================================================
		unsigned long int** getOrderedSegments(void); 

		//====================================================================================================================
		// File-I/O for this class, so that the sate of a cluster-object can be saved and loaded from/to a file
		// Linked data gets saved with this file and is then non-linked after loading, so memory-requirements may rise...
		//====================================================================================================================
		virtual bool save(const char *fileName, unsigned int listNr = 1);
		virtual SC_Cluster* load(const char *fileName, unsigned int listNr = 1);

	  //====================================================================================================================
	  // Methods for debug-output
	  //====================================================================================================================
	  friend  ostream& operator<<(ostream& os, SC_Cluster& cluster);
};

#endif
