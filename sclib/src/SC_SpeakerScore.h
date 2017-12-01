/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Abstarct base class for speaker-related scores (because there   */
/*        is some common functionality among them													*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.12.2006																								*/
/**************************************************************************/

#ifndef __SC_SpeakerScore_H__
#define __SC_SpeakerScore_H__

#include "SC_Api.h"
#include "SC_Score.h"
#include "SC_Partition.h"
#include "SC_GroundTruth.h"
#include "SC_TweakableParameters.h"
#include <SV_Error.h>

class SCLIB_API SC_SpeakerScore : public SC_Score {
	private :

  protected :

    class SC_IndexMapping {
      public:
        long int idx1;
        long int idx2;
        SC_IndexMapping *Next;
				SC_IndexMapping(long int index1, long int index2) : idx1(index1), idx2(index2), Next(NULL) {};
        int Valid(void) {return 1;}
    };

    SC_SpeakerScore::SC_IndexMapping *pClassMapping; //maps gt-speaker-ids to col-indixes into the scatter-matrix
		SC_SpeakerScore::SC_IndexMapping *pClusterMapping[4]; //maps hypo-(cluster-)speaker-ids to row-indexes into the scatter-matrix (may differ between CEs)

		unsigned long int ca; //nr. of rows/columns (classes, i.e. gt speakers) in the scatter matrix

    //====================================================================================================================
		//  Returns true if the given segment-borders do not appear in any of the clusters in the exlude-list
		//====================================================================================================================
		bool segmentIsUnknwon(long int segmentStart, long int segmentEnd, SC_Cluster *pExludeList);

		//====================================================================================================================
		//	Determines and returns the owner (gt-speaker-)id of the given cluster; if the linked list of cluster-objects in
		//  pSegmentExludeList is not empty, only segments from pClutser are considered for determining the owner that do not
		//  occur in one of the clusters of that list (needed to do id-scoring); in differentCa, a class-count different to
		//  the member this->ca can be given to assure correct array sizes if called by a method that temporarily manages
		//  more classes (as createMergeHistory() in SC_SpeakerScore_Clustering does)
		//====================================================================================================================
		long int determineClusterOwner(SC_Cluster *pCluster, unsigned long int algorithmicUncertaintyDiameter = 0, int countedEntity = sclib::ceSample, SC_Cluster *pSegmentExcludeList = NULL, unsigned long int differentCa = 0);

		//====================================================================================================================
		//	Returns the number of segments that belong to clutsers with the given owner-id
		//====================================================================================================================
		unsigned long int getNrOfSegmentsWithOwner(SC_Cluster *pClusters, long int ownerID);

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
		double calcDER(unsigned long int softBoundaryDiameter);

  public :

    //====================================================================================================================
    //	The constructor
    //====================================================================================================================
    SC_SpeakerScore(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT);

    //====================================================================================================================
    //	The destructor
    //====================================================================================================================
    virtual ~SC_SpeakerScore();

		//====================================================================================================================
		//	Fills the internal score-variables by computing their values according to the frameList in pGT, so that the
		//  get*()-Functions return reasonable values (before calling calcScores(), they return all 0)
		//  "start" and "end" refer to sample-numbers so that the area of the frameList for which scores shall be computed can
		//  be specified; this way, scores can be calculated only for parts of the video/corpus, e.g. for a scene.
		//  In algorithmicUncertaintyDiameter a value [in samples!] can be given which describes the precision with which the
		//  specific algorithm responsible for the results (and only the algorithm, not the gt...) can predict the place of 
		//  event-on- and -offsets
		//====================================================================================================================
    virtual void calcScores(unsigned long int start = 0, unsigned long int end = 0, unsigned long int algorithmicUncertaintyDiameter = 0) = 0;

    //====================================================================================================================
    //	Convert class-tags (gt-speaker-id's) to indices into the scatter-matrix/result-vectors or vice versa; 
		//  SVLIB_Fail is returned if the mapping can't be established
    //====================================================================================================================
		virtual long int class2idx(long int classTag);
		virtual long int idx2class(unsigned long int classIdx, char* clusterName = NULL, bool shortName = false);

    //====================================================================================================================
    //	Returns an or-concatenated list of the audio-types this class is responsible for scoring 
		//  (e.g. sclib::atSpeech|sclib::atNoise) or sclib::noType if there is no such type (e.g. in case of speaker id)
		//  This is useful e.g. to now which types can be fed into class2idx()
    //====================================================================================================================
		long int responsibility(void) {return sclib::noType;}
};

#endif
