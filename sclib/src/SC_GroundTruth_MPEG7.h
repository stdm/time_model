/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth do accomplish its targets in the    */
/*        context of MPEG7-Corpus files (Video's without! annotation)     */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 11.05.2006																								*/
/**************************************************************************/

#ifndef __SC_GroundTruth_MPEG7_H__
#define __SC_GroundTruth_MPEG7_H__

#include "SC_Api.h"
#include "SC_GroundTruth.h"
#ifdef SC_USE_JNI
	#include <jni.h>
#endif

class SCLIB_API SC_GroundTruth_MPEG7 : public SC_GroundTruth {

  private:

	protected :

		//====================================================================================================================
		// If a cutlist is available from videana/mediana (either fileName or list!)
		//====================================================================================================================
    char *sceneFileName;
		int *sceneList;
		int sceneListLength;
		int postEndCutCount; //tells how many cuts (if any) from the list lie behind the last frame reported by the adio-decoder (which is typically less than reported by the video-decoder, somehow)

		//====================================================================================================================
		// Methods to read in the results of the extern audio-/video-segmentation-module from a file and initializes the 
		// frameList
		//====================================================================================================================
		bool readSceneFile(const char* sceneFileName);

		//====================================================================================================================
		// Methods to read in the results of the extern audio-/video-segmentation-module from an integer-array of cuts/scenes
		// and initializes the frameList; the return value indicates success (true) or error (false)
		//====================================================================================================================
		bool readSceneList(int* sceneList, int sceneListLength);

 		//====================================================================================================================
		// Read Groundtruth from file(s) and initialize internal data structures (frameList etc.)
		//====================================================================================================================
    virtual bool readGroundTruth(void);

    unsigned long int videoFrame2FLI(unsigned long int frame, unsigned int alignment = sclib::alignmentStart);
    unsigned long int FLI2videoFrame(unsigned long int frameListIndex, unsigned int alignment = sclib::alignmentStart);

		//====================================================================================================================
		// Some classdata to remember important parameters of the signal, the video and the feature-vectors
		//====================================================================================================================
		double videoFrameSize; //in samples

	public :

    //====================================================================================================================
		// Constructor, Destructor
		//====================================================================================================================
    SC_GroundTruth_MPEG7(SC_TweakableParameters *pTweak, double videoFrameRate, const char* audioFileName, const char *sceneFileName = "");
#ifdef SC_USE_JNI
		SC_GroundTruth_MPEG7(SC_TweakableParameters *pTweak, double videoFrameRate, JNIEnv *env, jobject jStreamObject);
#endif
    SC_GroundTruth_MPEG7(SC_TweakableParameters *pTweak, double videoFrameRate, const char* audioFileName, int *sceneList, int sceneListLength);
		virtual ~SC_GroundTruth_MPEG7();

 		//====================================================================================================================
		// Methods to do some debugging output
		//====================================================================================================================
    virtual	void segmentStatisticsOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd);
		virtual ostream& output(ostream& OutS, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types = sclib::noType, unsigned long int typesNot = sclib::noType, int origin = sclib::modeHypothesized, bool suppressHeader = false); 
		virtual ostream& output(ostream& OutS, unsigned long int frameCount, unsigned long int firstFrameNr, long int *segmentationResult, long int *speakerIDs, double **probabilities, bool suppressHeader = false, const char *separator = "\t", long int selectionMask = sclib::atSilence|sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atBackground|sclib::atMusic|sclib::atAction|sclib::atUndefined, bool selectSpkID = false); 

 		//====================================================================================================================
		// check the ground-truth (or hypothesized results) for consistency; this includes the check for mutually exclusive
		// tupels of audio-types (e.g. no internal frame can be speech and noise together); returns true if everything is ok
		// violations of rules are reported in a file if debug output SCLIB_BD_GT_INCOSISTENCIES is turned on
		//====================================================================================================================
		virtual bool checkConsistency(unsigned long int segmentStart, unsigned long int segmentEnd, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// File-I/O for this class, so that the sate of a groundtruth-object can be saved and loaded from/to a file
		//====================================================================================================================
		virtual bool save(const char *fileName);
		virtual SC_GroundTruth* load(const char *fileName);

		//====================================================================================================================
		// Method to convert the findings of the algorithms (hypo-col of frameList, speaker-ids and probabiliy-list) to 
		// video-frame based lists that can be used by videana; memory for the new lists is allocated herein
		// The return value is the length of the 3 return-parameters
		//====================================================================================================================
		unsigned long int getResults(long int* &segmentationResults, long int* &speakerIDs, double** &probabilities);

		//====================================================================================================================
		//	output all analyzed parts of the frameList (really ALL parts are assumed to have been analyzed in the context of 
		//  this class & task!); mimic the behaviour of the base class by simultaneous using getResults(), which converts FLIs 
		//  to videoFrames
		//====================================================================================================================
		virtual void analyzedOut(const char* fileName, sclib::OpenMode mode = ios_base::out|ios_base::app) {
			analyzedOutEx(fileName, mode);
		}
    virtual void analyzedOutEx(const char* fileName, sclib::OpenMode mode = ios_base::out|ios_base::app, long int selectionMask = sclib::atSilence|sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atBackground|sclib::atMusic|sclib::atAction|sclib::atUndefined, bool selectSpkID = false);

		//====================================================================================================================
		// Reads a cutlist from the specified file and converts it to an int-array containing the end of a cut per entry; the 
		// length if this array is returned
		// (those cutList-arrays are subject to the audioSegmentation()-function in SC_Lib that is called by the java 
		//  frontend; this method is static because it doesn't need any object knowledge but needs to be called before 
		//  audioSegmentation() if only the file-based cutList is available)
		//====================================================================================================================
		static unsigned long int cutListFile2cutListArray(const char* fileName, int* &cutList);

		//====================================================================================================================
		// returns the video-framerate (in frames per second)
		//====================================================================================================================
		virtual double getVideoFrameRate(void);

		//====================================================================================================================
		// return the number of shots found by the video-decoder after the end of the file as reported by the audio-decoder
		//====================================================================================================================
		int getPostEndCutCount(void) {return this->postEndCutCount;}
};

#endif
