/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle Ralph's MPEG7 video-files,    */
/*        which have essentially no ground-truth available (rest like     */
/*        SCiVo)                                                          */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 11.05.2006																								*/
/**************************************************************************/

#ifndef __SC_Corpus_MPEG7_H__
#define __SC_Corpus_MPEG7_H__

#include "SC_Corpus.h"
#ifdef SC_USE_JNI
	#include <jni.h>
#endif

class SCLIB_API SC_Corpus_MPEG7 : public SC_Corpus {
	
  private:

  protected:
	
  public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Corpus_MPEG7(SC_TweakableParameters* pTweak, double videoFrameRate, const char *audioFileName, const char *sceneFile = "");
#ifdef SC_USE_JNI
		SC_Corpus_MPEG7(SC_TweakableParameters* pTweak, double videoFrameRate, JNIEnv *env, jobject jStreamObject);
#endif
    SC_Corpus_MPEG7(SC_TweakableParameters* pTweak, double videoFrameRate, const char *audioFileName, int *sceneList, int sceneListLength);
		virtual ~SC_Corpus_MPEG7();

    //====================================================================================================================
		// provides the possibility to load the desired samples for conformity reasons (it is a straight-forward call of
    // the corresponding method in SC_SegmentationHandler)
    // the given borders reamin unchanged in this method
		//====================================================================================================================
    SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries = false);

		//====================================================================================================================
    // write the given feature-set to a file formated to fit Ralph's needs (videoFrameStart, videoFrameEnd, features, ...)
		//====================================================================================================================
    void featureOut(const char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pData, bool printHeader = false, bool printVideoFrames = true, const char *featureName = NULL, bool soleHeader = false);

		//====================================================================================================================
		// another method to output features for input on java side: this time a CSV file is written, containing in the first
		// row the feature-name of each column (as a column name for db-tables), and then per row the features for one video-
		// frame, regardless of what was the audioframe-size (mean of all audio-frames beginning in the video-frame is 
		// written). the first column contains videoframe-numbers; pData is meant to be an array of datasets each containing 
		// a different feature (as returned by SC_FeatureHandler.extractFeatures())
		//====================================================================================================================
		virtual void lowLevelFeaturesOut(const char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pData, unsigned long int selectionMask = 0xFFFFFFFF);

		//====================================================================================================================
    // output all highlevel-features (=algorithmic results) for the whole video at once alike to the low-level features
		//====================================================================================================================
		virtual void highLevelFeaturesOut(const char *fileName, long int selectionMask = sclib::atSilence|sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atBackground|sclib::atMusic|sclib::atAction|sclib::atUndefined, bool selectSpkID = false);

		//====================================================================================================================
    // parses the given XML file and returns an array of video-frame-based shot-boundaries along with the size of the 
		// array
		//====================================================================================================================
		static int* getShotList(char *xmlFileName, double videoFrameRate, int &listLength);
};

#endif
