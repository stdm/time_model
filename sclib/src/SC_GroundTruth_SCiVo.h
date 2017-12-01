/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth do accomplish its targets in the    */
/*        context of SCiVo-Corpus files (Video's and annotation)          */
/*      - reads and organizes and provides access to SCiVo annotation     */
/*        ground-truth files                                              */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.02.2004																								*/
/**************************************************************************/

#ifndef __SC_GroundTruth_SCiVo_H__
#define __SC_GroundTruth_SCiVo_H__

#include "SC_Api.h"
#include "SC_GroundTruth.h"

class SCLIB_API SC_GroundTruth_SCiVo : public SC_GroundTruth {

  private:

	protected :

		//====================================================================================================================
		// Methods to read in the results of the extern audio-/video-segmentation-module and initialize the frameList
		//====================================================================================================================
		long int readSceneFile(const char* sceneFileName);
		long int readSegmentFile(const char* segmentFileName);

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
    char *segmentFileName;
    char *sceneFileName;

	public :

    //====================================================================================================================
		// Constructor, Destructor
		//====================================================================================================================
    SC_GroundTruth_SCiVo(SC_TweakableParameters *pTweak, double videoFrameRate, const char* audioFileName, const char* sceneFileName, const char* segmentFileName);
		virtual ~SC_GroundTruth_SCiVo();

		//====================================================================================================================
		// Methods to do some debugging output
		//====================================================================================================================
    virtual	void segmentStatisticsOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd);
		virtual ostream& output(ostream& OutS, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types = sclib::noType, unsigned long int typesNot = sclib::noType, int origin = sclib::modeHypothesized, bool suppressHeader = false); 

		//====================================================================================================================
		// File-I/O for this class, so that the sate of a groundtruth-object can be saved and loaded from/to a file
		//====================================================================================================================
		virtual bool save(const char *fileName);
		virtual SC_GroundTruth* load(const char *fileName);

		//====================================================================================================================
		// returns the video-framerate (in frames per second)
		//====================================================================================================================
		virtual double getVideoFrameRate(void);
};

#endif
