/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle Wei-Ho (MIX2MAX) Tsais test    */
/*        data                                                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 31.08.2005																								*/
/**************************************************************************/

#ifndef __SC_Corpus_Wesley_H__
#define __SC_Corpus_Wesley_H__

#include "SC_Corpus.h"
#include "SC_ModelHandler.h"
#include "SC_Cluster.h"

class SCLIB_API SC_Corpus_Wesley : public SC_Corpus {
	
  private:

  protected:
	
    SC_ModelHandler *pModeller;
    bool convertFeatures;
    unsigned long int modelOrigin;

    //====================================================================================================================
    //	load or train models, return linked list of clusters
		//  if the postfix is given, the trained models get saved immediately as if saved by savedModel()
    //====================================================================================================================
    SC_Cluster* constructModels(const char* targetListFileName, const char* corpusFileName, unsigned long int layer = sclib::modeForeground, const char* savePostFix = "");
    SC_Cluster* constructModel(const char* corpusFileName, const char* targetID, SC_Model* pBackgroundModel = NULL, const char* modelFileName = NULL, unsigned long int layer = sclib::modeForeground, const char* savePostFix = "");

    //====================================================================================================================
    //	collect all the data fpr one singer and return it as a cluster-object
    //====================================================================================================================
    SC_Cluster* constructFeatures(const char* corpusFileName, const char* targetID);

    //====================================================================================================================
		// functions to load a model saved by Wesleys algorithms
		//====================================================================================================================
    SC_Model* loadWesleyGMM(const char* fileName); //is always of type sclib::mtGMM_new
    SC_Model* loadWesleyMIXMAX(const char* fileName);
 
  public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Corpus_Wesley(SC_TweakableParameters* pTweak, unsigned long int sampleRate, bool mfcc2lfbe = false, unsigned long int modelOrigin = sclib::modeSClib);
		virtual ~SC_Corpus_Wesley();

    //====================================================================================================================
    //	builds a model for each ID mentioned in the targetListFile and found in the corpusFile; returns a linked list of
    //  clusters containing the models and the data;
		//  if savePostFix != '', the models are not only trained, but also saved immediately after creation (saving via 
		//  saveModels() needs all models trained before saving the first, so there is lot of unsafed work if training is 
		//  slow...)
    //====================================================================================================================
    SC_Cluster* trainModels(const char* targetListFileName, const char* corpusFileName, unsigned long int layer = sclib::modeForeground, const char* savePostFix = "");
    SC_Cluster* trainModel(const char* corpusFileName, const char* targetID, SC_Model* pBackgroundModel = NULL, unsigned long int layer = sclib::modeForeground, const char* savePostFix = "");

    //====================================================================================================================
    //	fill the cluster-objects as in trainModels(), but load the models from saved files instead of rebuilding them
    //  if loading a model fails, it is rebuilded from the feature vectors as in trainModels()
		//  if savePostFix != '', the models are not only trained, but also saved immediately after creation 
    //====================================================================================================================
    SC_Cluster* loadModels(const char* targetListFileName, const char* corpusFileName, unsigned long int layer = sclib::modeForeground);
    SC_Cluster* loadModel(const char* corpusFileName, const char* targetID, const char* modelFileName, SC_Model* pBackgroundModel = NULL, unsigned long int layer = sclib::modeForeground);
    
    //====================================================================================================================
    //	write the models to binary files; the filenames get the given postfix (should include the '.')
    //====================================================================================================================
    long saveModels(SC_Cluster* pModels, const char* postfix = ".gmm");
    long saveModel(SC_Cluster* pModelWithInfo, const char* postfix = ".gmm", bool verbose = true);

    //====================================================================================================================
    //	test the prebuild models against the files test-corpus
    //====================================================================================================================
    void testModels(const char* targetListFileName, const char* corpusFileName, SC_Cluster* trainedModels, const char* resultFileName);

    //====================================================================================================================
		// change the sampleRate
		//====================================================================================================================
    void setSampleRate(unsigned long int sampleRate);

    //====================================================================================================================
    //	alter the modelOrigin-member; return true, if new value is valid, otherwise return false and don't change anything
    //====================================================================================================================
    bool setModelOrigin(unsigned long int newOrigin);

    //====================================================================================================================
    //	convert MFCC features to log-FbE features by invoking an inverse dct on the MFCCs
    //====================================================================================================================
    void mfcc2lfbe(SV_Data* pMFCC);

    //====================================================================================================================
		// just that the interface is fullfilled... this function doesn't to anything in this corpus
		//====================================================================================================================
    SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries = false) {return NULL;};
};

#endif

