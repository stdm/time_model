/**************************************************************************/
/*    Responsibility:																											*/
/*      - Performs the speaker idetification of a linked list of features */
/*        against a set of speaker-clusters                               */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 26.04.2006																								*/
/**************************************************************************/

#ifndef __SC_SpeakerIdentificator_H__
#define __SC_SpeakerIdentificator_H__

#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Cluster.h"
#include "SC_MixtureModel.h"
#include "SC_SpeakerScore_Identification.h"

class SCLIB_API SC_SpeakerIdentificator {
	
  private:

  protected:
	
    //====================================================================================================================
		// There are some tweakable parameters in the SC_Lib library; they can be centraly managed in this class.
		//====================================================================================================================
    SC_TweakableParameters *pTweak;
    
    SC_MixtureModel *pBackgroundSpeakerModel;
    bool verbose;
  
    //====================================================================================================================
		// Performs score normalization (accounting for unknown impostors and imperfect features/models) under the "distance 
    // normalization" paradigm published in "A Monte-Carlo Method for Score Normalization in Automatic Speaker 
    // Verification using Kullback-Leibler Distances" by Ben, Blouet, Bimbot on ICASSP 2002
    // DNorm is comparable (in terms of recognition performance) to ZNorm and TNorm but dosn't need extra training data
    // or speaker models to learn the normalization parameters.
		//====================================================================================================================
    double dNorm(double clientScore, SC_Model *pClientModel, SC_MixtureModel *pWorldModel);

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_SpeakerIdentificator(SC_TweakableParameters *pTweak, SC_MixtureModel *pBackgroundSpeakerModel, bool verbose = true);
		virtual ~SC_SpeakerIdentificator();

    //====================================================================================================================
    // Identify to which cluster (if any) each feature-set in the linked list of candidates belongs; the candidates are 
    // also organized in "clusters" to enhance the features with segment-border information (no speaker-models are needed
    // therein!); the speakerList of clusters is assumed to contain *real* clusters as obtained by the speaker-clusterer 
    // object
		// updateClusters should be set ==true so that correct speaker-information is stroed in the groundtruth
		// TODO: - threshold-finding to exclude "impostors"; 0.0 is just the ideal threshold when the world-model scores 
		//         higher than any speaker-model...
		//       - at the moment, identification is performed based on the merged speaker model; integrate id based on other
		//         linkage methods, such as average, complete or single linkage based on single segment models
    //====================================================================================================================
		SC_SpeakerScore_Identification* identifySpeakers(SC_GroundTruth *pGT, SC_Cluster *pCandidateList, SC_Cluster *pSpeakerList, double impostorThreshold = 0.0, bool updateClusters = true, bool retrainClusters = false, long int progressCount = 0, long int progressMax = 0, bool (*doProceed)(long int &curr, long int max, unsigned int cnt) = NULL);
};

#endif
