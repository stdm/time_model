/**************************************************************************/
/*    Responsibility:																											*/
/*		  - provides the possbility to build speaker- or noise models of    */
/*        all implemented types                                           */
/*      - helps to save/load models to/from a file                        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.02.2006																								*/
/**************************************************************************/

#ifndef __SC_ModelHandler_H__
#define __SC_ModelHandler_H__

#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_Cluster.h"
#include "SC_Corpus.h"
#include <SV_Data.h>

class SCLIB_API SC_ModelHandler {
	private :

  protected :

    SC_TweakableParameters *pTweak;
    bool verbose;

  public :
		
    SC_ModelHandler(SC_TweakableParameters *pTweak, bool verbose = true);
    virtual ~SC_ModelHandler();

    //====================================================================================================================
    //	Change the parameter-container to extract features with differing parameter-settings
    //====================================================================================================================
    void setTweak(SC_TweakableParameters *pTweak) {this->pTweak = pTweak; return;}

 		//====================================================================================================================
		// Speaker Modelling with Noise-Compensation
		// explicitModels is a table whichs rows have the form 'sceneNr|segmentNr|fileName of explicit background model'. If 
		// explicitModelCount >0, it is parsed, and for a segment specified by scene- and segmentNr (0 is a wildcard in this 
    // case: 0 as sceneNr means all scenes, dito for segmentNr), the background-model will be loaded from the specified 
		// file instead of being estimated.
    // segment-boundarys are sample-based
    // returned are the speaker-models (together with their respective features and segments-boundarys packed into 
    // cluster-objects) of all valid (long enough) segments as the linked list 'pClusters'; the features and segments-
    // boundarys of all invalid segments are returned in the linked list 'pInvalidClusters' so that overall each single
    // speech-segment of the given segment is found in one of the clusters afterwards
		//====================================================================================================================
	  int buildSpeakerModels(SC_Corpus *pCorpus, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data* pFeatures, SC_Cluster* &pClusters, SC_Cluster* &pInvalidClusters, char*** explicitModels	= NULL, unsigned short int explicitModelCount = 0);

    //====================================================================================================================
    //	Used by buildSpeakerModels(): gets the features, a background-model and a name/id for this segment and constructs 
    //  the model, which is returned in a cluster-object (for this reason, the additionalInfo is needed); debug-outputting 
    //  is also handled; layer can be sclib::modeForeground or sclib::modeBackground to determine whether a speech- or 
    //  background model should be built; segment-start and -end are just for the purpose of storing the start- and 
    //  endsamples of the features in the cluster object.
    //====================================================================================================================
    SC_Cluster* buildModel(SV_Data *pFeatures, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Model* pBackgroundModel, unsigned long int layer, char* name = "", unsigned long int ID = 0);
    
    //====================================================================================================================
    //	Used by buildSpeakerModels(): gets the features and maybe a background-model for this segment and constructs 
    //  the model, which is returned; debug-outputting is also handled
    //  layer can be sclib::modeForeground or sclib::modeBackground to determine whether a speech- or background model 
    //  should be built
    //====================================================================================================================
    SC_Model* buildModel(SV_Data *pFeatures, SC_Model* pBackgroundModel, unsigned long int layer, unsigned long int segmentsToMerge);

    //====================================================================================================================
    //	Generates a Model for a given segment of feature-vectors of specific type
    //  segmentStart end segmentEnd (in samples) mus be start- and endsample of the given features!!!
    //====================================================================================================================
		SC_Model* buildModel(SC_GroundTruth* pGT, SV_Data* pFeatures, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Model *pBackground, unsigned long int types, unsigned long int typesNot, unsigned short int modelOrder, unsigned short int modelType);

    //====================================================================================================================
    //	Creates a new object of the desired type
    //====================================================================================================================
    SC_Model* createRawModel(unsigned long int modelType, SC_Model *pBackgroundModel = NULL, unsigned short int modelOrder = 0, unsigned short int dim = 1);

    //====================================================================================================================
    //	Generates and saves a Model for a given segment of feature-vectors of specific type; returns nr. of bytes written
    //====================================================================================================================
    int saveModel(const char *fileName, SC_GroundTruth* pGT, SV_Data* pFeatures, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Model *pBackground, unsigned long int types, unsigned long int typesNot, unsigned short int modelOrder, unsigned short int modelType);

    //====================================================================================================================
		//	Generates a Model and returns it
		//  in contrast to the abovementioned buildModel()-method, this one can merge all the frames in the linked list 
		//  pFeature together and builds the model instead of copying frames of desired type together; if the given 
		//  modelOrder==0, a search for the best order is conducted
    //====================================================================================================================
    SC_Model* buildModel(SV_Data *pFeatures, SC_Model* pBackground, unsigned short int modelOrder, unsigned long int modelType, unsigned long int segmentsToMerge);

    //====================================================================================================================
    //	Generates a Model and saves it while returning nr. of bytes written
    //  in contrast to the abovementioned buildModel()-method, this one just merges all the frames in the linked list 
    //  pFeature together  and  builds the model instead of copying frames of desired type together
    //====================================================================================================================
    int saveModel(const char *fileName, SV_Data *pFeatures, SC_Model* pBackground, unsigned short int modelOrder, unsigned short int modelType);
    
    //====================================================================================================================
    //	save a model to a file; handles report-printing as well; uses debug-dir as a prefix to the filename only if wanted
    //====================================================================================================================
    int saveModel(const char *fileName, SC_Model *pModel, bool useDebugDir = true);

    //====================================================================================================================
    //	load a model from a file according to the tweakable parameters
    //====================================================================================================================
    SC_Model* loadModel(const char *fileName, unsigned long int modelType);

    //====================================================================================================================
    //	Creates a copy of the given parent-model of its own kind
    //====================================================================================================================
    SC_Model* copyModel(SC_Model *pParentModel, bool justLink = false);

    //====================================================================================================================
    //	Combine 2 models by adding their "components"
    //  In contrats to the methods in the different models themselfes, this function "knows" which model can be combined
    //  with which one
    //====================================================================================================================
    SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);

    //====================================================================================================================
    //	Combine 2 models by retraining a new one on the concatenated training data of both parants
    //  In contrats to the methods in the different models themselfes, this function "knows" which model can be combined
    //  with which one
    //====================================================================================================================
    SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);

    //====================================================================================================================
    // Computes the SNR between given signal- and noise features. The features are meant to contain frame-Energys in 
    // col 0 (e.g. EnergyZCR features, or MFCC). If they contain other featuires, this function can also calculate a SNR 
    // over 0 to maxColIdx cols, but maybe the resukt isn't meaningful then.
    //====================================================================================================================
    double getSNR(SV_Data* pSignalEnergy, SV_Data* pNoiseEnergy, unsigned int maxColIdx = 0);

    //====================================================================================================================
    // For mixture models: Guess the optimal model-order depending on the guessMode in the tweakable parameters:
    //  - Do a search over the space of available modelorders, building a rough model (only few EM steps)  
    //    for each and returning the order which maximized BIC criterion
    //  ...or...
    //  - Following a heuristic, which maps the length of the segment to the # of gaussians
    //====================================================================================================================
    unsigned short int guessModelOrder(SV_Data* pFeatures, SC_Model* pBackground, unsigned long int modelType, unsigned short int minOrder = 1, unsigned short int maxOrder = 2048, unsigned long int segmentsToMerge = 1);

    //====================================================================================================================
    //	Returns true if the given model is a child of SC_MixtureModel
    //====================================================================================================================
    bool isMixtureModel(SC_Model* pModel);
    bool isMixtureModel(long int modelType);

    //====================================================================================================================
    //	Returns true if the given model type returns averaged scores (i.e. scores divided by the number of test patterns)
    //====================================================================================================================
		bool averagesItsScores(long int modelType);

    //====================================================================================================================
    //	This methods returns the score for the given data and the given model; the speciality is the last parameter can be
    //  used to force a model to compute un-normalized (i.e. no likelihood rations as in the GMM-UBM) scores without 
    //  knowledge of the particular abilities of the model; this way, all special knowlege about the model-internas 
    //  remains within this class; as a sugar, the score is retuned as a double instead of an SV_Data object; if 
		//  averageScores is false, the result will not be divided by the number of testpatterns as is usually (but not al-
		//  ways) done by the models.
    //====================================================================================================================
    double testModel(SC_Model* pModel, SV_Data* pData, unsigned long int segmentsToMerge, bool forceNoNormalization = false, SC_Model* pBackgroundModel = NULL, bool averageScores = true);
};

#endif

