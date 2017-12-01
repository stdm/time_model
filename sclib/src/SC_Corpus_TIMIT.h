/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle the TIMIT corpus              */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.04.2006																								*/
/**************************************************************************/

#ifndef __SC_Corpus_TIMIT_H__
#define __SC_Corpus_TIMIT_H__

#include "SC_Corpus.h"
#include "SC_Signal.h"
#include "SC_GroundTruth_TIMIT.h"
#include <SV_Data.h>

class SCLIB_API SC_Corpus_TIMIT : public SC_Corpus {
	
  private:

  protected:
	  
    bool verbose; //controls some console-output

  public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Corpus_TIMIT(SC_TweakableParameters* pTweak, const char *corpusFileName, bool verbose = true, bool setSceneBoundaries = true, bool setAllSpeakerBoundaries = false, unsigned int sbSkipCount = 0, bool sceneBoundariesOnlyAtSpeakerBoundaries = false);
		virtual ~SC_Corpus_TIMIT();

		//====================================================================================================================
		// provides the possibility to load the desired samples even if they are located in different files
		// due to the huge size of the TIMIT corpus it is not recommended to load it into memory in one single part, but
		// to operate on it scene-wise ("scenes" are here "files").
		// if the given segment-borders nearly fit the real borders of a file (nearly=within fewer samples than 1 GT-internal
		// frameSize), the segment-borders are aligned to the file-borders (that's the reason for the reference parameters)
		//====================================================================================================================
    SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries = false);

    //====================================================================================================================
    // don't use this method!!! reason:
    // at the end of each file there are some frames missing because they don't fit completely into the loaded signal
    // part (only the samples of one file are loaded at one moment). as more and more files are loaded and features are 
    // extracted, the frame-numbers get more and more out of sync with the way frames are mapped to samples (and vice 
    // versa) in the ground-truth object.
    // so better use scene-based analysis as in SCiVo or the TIMIT data (scene-boundarys are set at each file's first 
    // sample) for recognition/dectection involving the groundtruth, and then maybe copy samples together (using the
    // copyFramesTogether()-method of the ground-truth object) to build bigger models etc.
    // nevertheless, this method's implementation can serve as an example for how to concatenate the feature-sets after
    // all ground-truth work has been done.
    //====================================================================================================================
		// there might be a memory problem when first loading a huge part of the corpus by calling loadSignal and the doing
    // feature extraction (using the SC_FeatureHandler class) on this large sample-array (feature extraction even 
    // involves min. 1 copy of the comlete samples...).
    // so this method provides a way of specifiying which features should be extracted on which part of the signal 
    // (can be located in different files as above), and the features get extracted file-by-file and are then concatenated
		//====================================================================================================================
    //SV_Data** extractFeatures(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int featureTypes);

    void setVerbose(bool verbose) {this->verbose = verbose; return;}

		//====================================================================================================================
    // Return the fileName of the audio file containg the given startSample
		//====================================================================================================================
		virtual char* getCurrentAudioFileName(unsigned long int startSample) {unsigned long int offset, sampleCount; return ((SC_GroundTruth_TIMIT*)this->getGT())->whichFile(startSample, offset, sampleCount);}
};

#endif
