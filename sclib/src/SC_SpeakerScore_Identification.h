/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the speaker       */
/*        identification according to the ground truth                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 01.05.2006																								*/
/**************************************************************************/

#ifndef __SC_SpeakerScore_Identification_H__
#define __SC_SpeakerScore_Identification_H__

#include "SC_SpeakerScore.h"
#include "SC_Cluster.h"
#include "SC_Partition.h"

class SCLIB_API SC_SpeakerScore_Identification : public SC_SpeakerScore {
	private:

protected:

		SC_Partition *pCandidatePartition, *pSpeakerPartition; //hold the given linked lists of clusters in partition-containers to have more elegant ways of access (per id/name)
		SC_SpeakerScore::SC_IndexMapping *pSpeakerMapping; //to remember the mapping between clusters in the speaker- and cacdidate-list the id-process generated

    //====================================================================================================================
    //	These score-variables get filled by calcScores(), before they are initialized with 0
    //====================================================================================================================
    double accuracy;
		
    //====================================================================================================================
		//  Adds a mapping between candidates and speakers
		//====================================================================================================================
		void addSpeakerMapping(long int speakerId, long int candidateId);

    //====================================================================================================================
		//  To make operator<<() kind of virtual...
		//====================================================================================================================
		virtual ostream& output(ostream& OutS); 

	public:

    //====================================================================================================================
    //	The constructor
    //  the final partition must be a pointer to an onject within the partition-list in order to get destructed!
    //====================================================================================================================
    SC_SpeakerScore_Identification(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT, SC_Cluster *pCandidateList, SC_Cluster *pSpeakerList);

    //====================================================================================================================
    //	The destructor
    //  destructs all the linked partitions, too!!!
    //====================================================================================================================
    virtual ~SC_SpeakerScore_Identification();

		//====================================================================================================================
		//	Fills the internal score-variables by computing their values according to the frameList in pGT, so that the
		//  get*()-Functions return reasonable values (before calling calcScores(), they return all 0)
		//  "start" and "end" refer to sample-numbers so that the area of the frameList for which scores shall be computed can
		//  be specified; this way, scores can be calculated only for parts of the video/corpus, e.g. for a scene.
		//  In algorithmicUncertaintyDiameter a value [in samples!] can be given which describes the precision with which the
		//  specific algorithm responsible for the results (and only the algorithm, not the gt...) can predict the place of 
		//  event-on- and -offsets
		//====================================================================================================================
    virtual void calcScores(unsigned long int start = 0, unsigned long int end = 0, unsigned long int algorithmicUncertaintyDiameter = 0); 
    
    //====================================================================================================================
    //	These get*()-functions give access to the results computed by calcScores()
    //====================================================================================================================
    virtual double getRecognitionRate(void) {return this->accuracy;}
		virtual double getAccuracy(void) {return this->accuracy;}

    //====================================================================================================================
    //	Functions to set certain values during the clustering-process to prepare score-calculation
    //====================================================================================================================
		void setCandidateList(SC_Cluster *pCandidates);
		void setSpeakerList(SC_Cluster *pSpeakers);

    //====================================================================================================================
    //	Returns an or-concatenated list of the audio-types this class is responsible for scoring 
		//  (e.g. sclib::atSpeech|sclib::atNoise) or sclib::noType if there is no such type (e.g. in case of speaker id)
		//  This is useful e.g. to now which types can be fed into class2idx()
    //====================================================================================================================
		long int responsibility(void) {return sclib::noType;}
};

#endif
