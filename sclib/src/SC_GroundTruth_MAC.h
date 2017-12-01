/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth to accomplish its targets in the    */
/*        context of MAC-Corpus files (thilo's movie sound audio          */
/*        classification corpus)                                          */
/*      - makes parts of the corpus available as if they are one big file */
/*      - operates on lists of MAC-files specifiying the "parts"          */
/*      - each file contains only one audio-type (music/background/       */
/*        noisy speech/pure speech/breath/action) and maybe silence-      */
/*        regions                                                         */
/*      - the groundtruth about the audio-type is extracted from the      */
/*        name of the first subdirectory under the corpus' base path;     */
/*        there is no groundtruth about silence regions!                  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 12.06.2006																								*/
/**************************************************************************/

#ifndef __SC_GroundTruth_MAC_H__
#define __SC_GroundTruth_MAC_H__

#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_GroundTruth.h"

class SCLIB_API SC_GroundTruth_MAC : public SC_GroundTruth {

  private:

	protected :

    class SC_FileList {
      public:
        char fileName[sclib::bufferSize];
        unsigned long int startSample;
        unsigned long int endSample;
        long int types; //audio-type contained in this file (possible because of the constarined of only one type (and silence) per file in this corpus)
        SC_FileList *Next;
        int Valid(void) {return 1;}
    };

    bool verbose; //to control some console output

    //====================================================================================================================
		// Remember which files content was mapped to which part of the frameList
		//====================================================================================================================
		SC_GroundTruth_MAC::SC_FileList *fileList;
    void addFile(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd);

    //====================================================================================================================
		// This class holds the groundtruth for a whole corpus as if it where a aignle file; this function reads the 
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
    SC_GroundTruth_MAC(SC_TweakableParameters *pTweak, const char* corpusFileName, bool verbose = true);
		virtual ~SC_GroundTruth_MAC();

    //====================================================================================================================
    // get the filename out of the internal list holding the desired start-sample; the sample-nr out of the frameList 
    // corresponding with the first sample in this file is also returned as offset; the sample-count in this file is 
    // returned in sampleCount
    //====================================================================================================================
    char* whichFile(unsigned long int segmentStart, unsigned long int &offset, unsigned long int &sampleCount, long int &types);

    void setVerbose(bool verbose) {this->verbose = verbose; return;}

		//====================================================================================================================
		// File-I/O for this class, so that the sate of a groundtruth-object can be saved and loaded from/to a file
		//====================================================================================================================
		virtual bool save(const char *fileName);
		virtual SC_GroundTruth* load(const char *fileName);
};

#endif
