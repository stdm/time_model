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

#ifndef __SC_Partition_H__
#define __SC_Partition_H__

#include "SC_Api.h"
#include "SC_Cluster.h"
#include "SC_GroundTruth.h"

class SCLIB_API SC_Partition {

	private :

	protected:

    SC_Cluster *pClusters;
    double globalCriterion;
    long classSig;
		SC_Cluster *pFather, *pMother, *pChild; //pointers to the origins (father & mother) and results (child) of merging: father and mother in this partition have a child in the successive partition; the child in this partition is the cluster that resulted from merging father and mother in the last partition => the first partition in a linked list has no child, the last one no parents!
		unsigned long int partitionNr; //index of the partition in the list created successively by the hierarchical clustering process; first partition (all small clusters) has partitionNr==0

	public :

		//====================================================================================================================
		//	constructor
    //  if justLink=true, the new objects gets just a pointer to the given list of clusters; otherwise, it makes a copy
    //  and stores it. if noWorkData==true (and justLink==false), the speech-frames and models in the clusters aren't 
		//  disposed at all (like calling the constructor and realeaseWorkData() directly after each other)
		//====================================================================================================================
	  SC_Partition(SC_Cluster* pClusters, double globalCriterion, unsigned long int partitionNr, bool justLink = false, bool noWorkData = false);
    SC_Partition(const SC_Partition& pParent, bool justLink = false);

    //====================================================================================================================
		//	destructor
    //  the clusters get destructed, regardless if they where copys or just linked
		//====================================================================================================================
		virtual ~SC_Partition();

    bool Valid(void);

    //====================================================================================================================
		//	provide access to protected members
		//====================================================================================================================
    SC_Cluster* getClusters(void) {return this->pClusters;}
    SC_Cluster* getClusterByName(const char* speakerName);
    SC_Cluster* getClusterByID(unsigned long int speakerID);
    int getClusterIdxByID(unsigned long int speakerID);
    double getGlobalCriterion(void) {return this->globalCriterion;}
		SC_Cluster* getFather(void) {return this->pFather;}
		SC_Cluster* getMother(void) {return this->pMother;}
		SC_Cluster* getChild(void) {return this->pChild;}
		    
    void setGlobalCriterion(double value) {this->globalCriterion = value; return;}
    void setClusters(SC_Cluster *pFirst) {this->pClusters = pFirst; return;}
		void setFather(SC_Cluster *pFather) {this->pFather = pFather; return;}
		void setMother(SC_Cluster *pMother) {this->pMother = pMother; return;}
		void setChild(SC_Cluster *pChild) {this->pChild = pChild; return;}

    //====================================================================================================================
		//	to form a linked list of succesive partitions during a clustering
		//====================================================================================================================
    SC_Partition *Next;

    //====================================================================================================================
		// to store the ID-information found in the cluster-objects in the frameList of the groundtruth-object
    // that means: remove all previous speaker-boundary/speech-seg-start/speech-seg-end/speaker-id labels and set them 
    // according to the information stored in the clusters of this partition
		//====================================================================================================================
    void storeIDs(SC_GroundTruth *pGT, bool removeOldSpeakerBoundarys = true, bool removeOldSpeakerIDs = true, bool removeOldSpeechSegments = true);

    //====================================================================================================================
		//	return the uncertainty diameter due to clustering, i.e. two framesteps of the used features
		//====================================================================================================================
		unsigned long int getUncertaintyDiameter(void);

    //====================================================================================================================
		//	destructs all the basic data-structures ((background-)models and speechFrames) of all clusters in this partition
    //  this saves much of memory when all computations are done and the object is just holded to store the results in a 
    //  partition
		//====================================================================================================================
		void releaseWorkData(void);

    //====================================================================================================================
		//	if this is not the first partition in a succession of partitions created during clustering, it doesn't hold any
		//  work data, e.g. speech frames; but the first partition of the process does; this method returns, for a given id
		//  of one of its clusters, a linked list of *new* objects (the objects are new and need destruction afterwards, the 
		//  data itself is just linked to save space) representing all frames for this cluster, which can be found in this
		//  first partition via the segment-boundaries, which don't get changed during clustering
		//====================================================================================================================
		SV_Data* getSpeechFramesOfCluster(unsigned long int speakerID, SC_Partition *pFirstPartition);

    //====================================================================================================================
		//	as above, for all clusters in the current partition; the speech frames are directly stored in the corresponding
		//  cluster objects
		//====================================================================================================================
		void getSpeechFrames(SC_Partition *pFirstPartition);
};

#endif
