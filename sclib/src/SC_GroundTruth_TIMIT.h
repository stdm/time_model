/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth to accomplish its targets in the    */
/*        context of TIMIT-Corpus files (LDC's pure speech corpus)        */
/*      - reads and organizes and provides access to TIMIT annotation/    */
/*        ground-truth files                                              */
/*      - makes parts of the corpus available as if they are one big file */
/*      - operates on lists of TIMIT-files specifiying the "parts"        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 07.04.2006																								*/
/**************************************************************************/

#ifndef __SC_GroundTruth_TIMIT_H__
#define __SC_GroundTruth_TIMIT_H__

#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_GroundTruth.h"

class SCLIB_API SC_GroundTruth_TIMIT : public SC_GroundTruth {

  private:

	protected :

    class SC_FileList {
      public:
        char fileName[sclib::bufferSize];
        unsigned long int startSample;
        unsigned long int endSample;
        SC_FileList *Next;
        int Valid(void) {return 1;}
    };

    bool verbose; //to control some console output
		bool setSceneBoundariesAtFileBoundaries; //if true, artificial scene-boudaries are included in the groundtruth at the beginning of a new file
		bool setSpeakerBoundariesAtFileBoundaries; //if true (non-standard), a speaker-boundary is introduced at each beginning of a new file (first speech occurence). might be useful for speaker recognition experiments on TIMIT without shuffeling the corpus files...
		bool setSceneBoundariesAtSpeakerBoundaries; //if true (non-standard), artificial scene-boundaries are included in the groundtruth at the beginning of a new speaker segment
		unsigned int sbSkipCount; //if setSpeakerBoundariesAtFileBoundaries==true (and setSceneBoundariesAtFileBoundaries==false), this parameter is applied: it tells how many of those speaker boundaries at file boundaries should be skipped (i.e. not be set) if the files are speaker-homogenious; this way, longer speaker-homogenious segments can be created without copying files together

    //====================================================================================================================
		// The framelist: 
    //  - The 0. component holds an entry for every audioFrame describing it's content with the bitflag-constants 
    //    defined above, as stated in the ground-truth. 
    //  - The 1. component holds the same, but as hypothesized by the algorithms of this lib. 
    //  - The 2. component holds the ground-truth Speaker-ID for this frame (if available)
    //  - The 3. component holds the hypothesized Speaker-ID as computed by the classification algorithms
    //  - The 4. component holds the phone-type (only one at a time) as stated by the ground-truth
    //  - The 5. component holds the phone-type as hypothesized by the algorithms
		//====================================================================================================================
		//long int** frameList; //already defined in the base class; comment appears here only becaus of slightly changed meaning ( 4th+5th component)

    //====================================================================================================================
		// Remember which file's content was mapped to which part of the frameList
		//====================================================================================================================
		SC_GroundTruth_TIMIT::SC_FileList *fileList;
    void addFile(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd);

		//====================================================================================================================
		// Initialize the frameList with it's new 4. column
		//====================================================================================================================
		virtual void initFrameList(void);

		//====================================================================================================================
		// Extracts the information regarding audio-types from the given phone; also returns a precise (numerical) 
    // representation of the phone-type, if a variable to alter is given (!= NULL)
		//====================================================================================================================
    unsigned long int phone2audioType(const char* phone, unsigned long int &phoneType);

    //====================================================================================================================
		// This class holds the groundtruth for a whole corpus as if it where a single file; this function reads the 
    // parameters of all files in the corpus to measure it's overall size and create a mapping between the frames in 
    // the framelist and the corresponding files storing the actual signals
		//====================================================================================================================
 		bool measureCorpusSize(const char *corpusFileName);

    //====================================================================================================================
		// Read Groundtruth from file(s) and initialize internal data structures (frameList etc.)
		//====================================================================================================================
    virtual bool readGroundTruth(void);

	public :

    //====================================================================================================================
		// Constructor, Destructor
		//====================================================================================================================
    SC_GroundTruth_TIMIT(SC_TweakableParameters *pTweak, const char* corpusFileName, bool verbose = true, bool setSceneBoundariesAtFileBoundaries = true, bool setSpeakerBoundariesAtFileBoundaries = false, unsigned int sbSkipCount = 0, bool sceneBoundariesAtSpeakerBoundaries = false);
		virtual ~SC_GroundTruth_TIMIT();

    //====================================================================================================================
    // Insert a phone-label into the framelist
    //====================================================================================================================
    void setPhone(unsigned long int segmentStart, unsigned long int segmentEnd, int phoneType, int origin = sclib::modeHypothesized);

		//====================================================================================================================
    // Returns a list with all phones occuring in the specified segment in order of appearence
    //====================================================================================================================
    int* getPhones(unsigned long int segmentStart, unsigned long int segmentEnd, int &listSize, int origin = sclib::modeHypothesized);

		//====================================================================================================================
    // Adds a new column to the feature-set (new 0th column) containing phone-labels for each frame; if replicate==true,
		// frames containing more than one phone (the majority) will be replicated to form one example for each of its phones;
		// otherwise, each frame will be decorated with the phone-label that is valid for most of the time of the frame.
    //====================================================================================================================
		void addPhoneLabels(SV_Data* &pData, unsigned long int segmentStart, unsigned long int segmentEnd, int origin = sclib::modeHypothesized, bool replicate = false);

		//====================================================================================================================
    // Return the string-representation of a phones to a given phone-label; one entry is per FLI
    //====================================================================================================================
		const char* phoneType2string(int phoneType);

		//====================================================================================================================
		// Returns for a given phone-label it's audioType; reports incosistencies in the matchings between phone-labels 
		// and the corresponding string-prepresentations (found in phoneType2string() and phone2audioType(.) for both 
		// directions
		//====================================================================================================================
		unsigned long int phone2audioType(unsigned long int phoneType);

    //====================================================================================================================
    // Get the filename out of the internal list holding the desired start-sample; the sample-nr out of the frameList 
    // corresponding with the first sample in this file is also returned as offset; the sample-count in this file is 
    // returned in sampleCount
    //====================================================================================================================
    char* whichFile(unsigned long int segmentStart, unsigned long int &offset, unsigned long int &sampleCount);

    void setVerbose(bool verbose) {this->verbose = verbose; return;}

		//====================================================================================================================
		// File-I/O for this class, so that the sate of a groundtruth-object can be saved and loaded from/to a file
		//====================================================================================================================
		virtual bool save(const char *fileName);
		virtual SC_GroundTruth* load(const char *fileName);

		//====================================================================================================================
		// Many experiments on TIMIT assign each audio file of a speaker (there are 10 overall) to either of two utterances: 
		// based on a lexicographic ordering of all 10 filenames per speaker, the first 8 are regarded as utterance 1, the
		// remaining 2 are utterance 2. For the given audioFilenName and speakername, this method returns the utterance number
		// returns sclib::noSpeaker on error
		//====================================================================================================================
		int getUtteranceNumber(const char *audioFileName, const char *speakerName);

		//====================================================================================================================
		// Extracts and returns the letter F or M in front of the TIMIT speaker name (last folder name)
		//====================================================================================================================
		char getGenderFromFilePath(const char* audioFileName);
};

#endif
