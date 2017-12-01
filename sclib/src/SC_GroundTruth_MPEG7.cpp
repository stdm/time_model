/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth do accomplish its targets in the    */
/*        context of MPEG7-Corpus files (Video's without! annotation)     */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 11.05.2006																								*/
/**************************************************************************/

#include <list>
#include <math.h>
#include <limits>
#include <string.h>
#include <iomanip>
#include "SC_GroundTruth_MPEG7.h"
#include "SC_Aux.h"
#include "SC_Signal_MPEG.h"
#include "SC_Signal_jWAVE.h"
#include <SV_DataIO.h>

/**************************************************************************/
/*       	                Constructor, Destructor     									  */
/**************************************************************************/

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth_MPEG7::SC_GroundTruth_MPEG7(SC_TweakableParameters *pTweak, double videoFrameRate, const char* audioFileName, const char* sceneFileName) : SC_GroundTruth(pTweak, audioFileName) {
  //some init-stuff
  this->gtType = sclib::gtMPEG7;
	this->uncertaintyRegion = 0;
	this->postEndCutCount = 0;

	//set video-framerate; try to extract it from the signal-parameters if this is stored in a SC_Signal_MPEG class
	if ((this->pSignalPrototype->getSignalType() == sclib::stMP2 ||
		   this->pSignalPrototype->getSignalType() == sclib::stMP2 ||
			 this->pSignalPrototype->getSignalType() == sclib::stMPEG ||
			 this->pSignalPrototype->getSignalType() == sclib::stOther) &&
	  	 ((SC_Signal_MPEG*)this->pSignalPrototype)->getHeader()->videoFrameRate > 0.0) {
		this->videoFrameSize = (double)(this->pSignalPrototype->SigPar.SRate) / ((SC_Signal_MPEG*)this->pSignalPrototype)->getHeader()->videoFrameRate; //length of a video-frame in samples
	} else {
		this->videoFrameSize = (double)(this->pSignalPrototype->SigPar.SRate) / videoFrameRate; //length of a video-frame in samples
	}

	this->pConverter->setVideoFrameSize(this->videoFrameSize);
  this->sceneFileName = new char[strlen(sceneFileName)+1];
	strcpy(this->sceneFileName, sceneFileName);
	this->sceneList = NULL;
	this->sceneListLength = 0;

  //initialize the framelist
  readGroundTruth();
}

#ifdef SC_USE_JNI
//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth_MPEG7::SC_GroundTruth_MPEG7(SC_TweakableParameters *pTweak, double videoFrameRate, JNIEnv *env, jobject jStreamObject) : SC_GroundTruth(pTweak, "") {
	//some init-stuff
  this->gtType = sclib::gtMPEG7;
	this->uncertaintyRegion = 0;
	this->postEndCutCount = 0;

	//do what wasn't done in the base class constructore due to not loading a sound from file
  this->audioFileName = new char[8];
  sprintf(this->audioFileName, "%s", "jStream");
  this->pSignalPrototype = new SC_Signal_jWAVE();
	((SC_Signal_jWAVE*)this->pSignalPrototype)->setStream(env, jStreamObject);
	this->pSignalPrototype->LoadSignal("jStream"); //for jStream source, load and cache the complete signal in the groundtruth class!
  this->audioSampleRate = this->pSignalPrototype->SigPar.SRate;
	this->pConverter->setAudioSampleRate(this->audioSampleRate);
  this->audioSampleCount = this->pSignalPrototype->getSampleCount();
  this->internalFrameSize = this->pConverter->ms2sample(this->pTweak->groundTruth.internalFrameSize); //in samples!!!
  this->internalFrameCount = (unsigned long int)(ceil((double)(this->audioSampleCount) / (double)(this->internalFrameSize)));
  initFrameList();

	//set video-framerate
	this->videoFrameSize = (double)(this->pSignalPrototype->SigPar.SRate) / videoFrameRate; //length of a video-frame in samples

	this->pConverter->setVideoFrameSize(this->videoFrameSize);
  this->sceneFileName = new char[2];
	sprintf(this->sceneFileName, "%s", "");
	this->sceneList = NULL;
	this->sceneListLength = 0;

  //initialize the framelist
  readGroundTruth();
}
#endif

//====================================================================================================================
//	constructor that receives an integer-array as a scene/cut-list  instead of a file containing the video-frame 
//  nubers where changes occured
//====================================================================================================================
SC_GroundTruth_MPEG7::SC_GroundTruth_MPEG7(SC_TweakableParameters *pTweak, double videoFrameRate, const char* audioFileName, int *sceneList, int sceneListLength) : SC_GroundTruth(pTweak, audioFileName) {
  //some init-stuff
  this->gtType = sclib::gtMPEG7;
  //this->videoFrameSize = (double)(this->pSignalPrototype->SigPar.SRate) / videoFrameRate; //length of a video-frame in samples
	this->uncertaintyRegion = 0;
	this->postEndCutCount = 0;

	//set video-framerate; try to extract it from the signal-parameters if this is stored in a SC_Signal_MPEG class
	if ((this->pSignalPrototype->getSignalType() == sclib::stMP2 ||
		   this->pSignalPrototype->getSignalType() == sclib::stMP2 ||
			 this->pSignalPrototype->getSignalType() == sclib::stMPEG ||
			 this->pSignalPrototype->getSignalType() == sclib::stOther) &&
	  	 ((SC_Signal_MPEG*)this->pSignalPrototype)->getHeader()->videoFrameRate > 0.0) {
		this->videoFrameSize = (double)(this->pSignalPrototype->SigPar.SRate) / ((SC_Signal_MPEG*)this->pSignalPrototype)->getHeader()->videoFrameRate; //length of a video-frame in samples
	} else {
		this->videoFrameSize = (double)(this->pSignalPrototype->SigPar.SRate) / videoFrameRate; //length of a video-frame in samples
	}

	this->pConverter->setVideoFrameSize(this->videoFrameSize);
	this->sceneFileName = new char[2];;
	sprintf(this->sceneFileName, "\0");
	this->sceneList = sceneList;
	this->sceneListLength = sceneListLength;

  //initialize the framelist
  readGroundTruth();
	
	this->sceneList = NULL; //we only need this as a member because readGroundTruth is parameter-free, but don't want to save it twice (in addition to the frameList, where it's values have been copied to)
	this->sceneListLength = 0;
}

//====================================================================================================================
//	destructor 
//====================================================================================================================
SC_GroundTruth_MPEG7::~SC_GroundTruth_MPEG7() {
  MFree_1D(this->sceneFileName);
}

//====================================================================================================================
//	initialize the framelist with the data of the scene-/segment-List
//====================================================================================================================
bool SC_GroundTruth_MPEG7::readGroundTruth(void) {
	bool success = true, generateScenes = true;
	long int res = SVLIB_Ok;

  //initialization of the frameList with zeros for the labels and -1 for the speaker-id's has been already done by initFrameList()
	//so care for scene/shots now
	if (strncmp(this->sceneFileName, "", sclib::bufferSize) != 0) {
		success = readSceneFile(this->sceneFileName);
		generateScenes = !success;
	}	else if (this->sceneList != NULL && this->sceneListLength > 0) {
		success = readSceneList(this->sceneList, this->sceneListLength);
		generateScenes = !success;
	}
	
	if (generateScenes == true) {
		//set an artificial scene start at the first frame
		setFrame(0, sclib::atSceneBoundary|sclib::atShotBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
		setFrame(0, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);

		//generate pseudo scene boundarys if wished
		if (this->pTweak->groundTruth.pseudoSceneLength > 0) {
			generatePseudoSceneBoundarys(this->pConverter->ms2sample(this->pTweak->groundTruth.pseudoSceneLength));
		}
	}
  
  return success;
}

//====================================================================================================================
//	read scene-list from file to internal data-structure
//====================================================================================================================
bool SC_GroundTruth_MPEG7::readSceneFile(const char* filename) {
	int *cutList = NULL;
	int length = cutListFile2cutListArray(filename, cutList);
	bool res = (length > 0) ? readSceneList(cutList, length) : false;
	
	MFree_1D(cutList);

	return res;
}

//====================================================================================================================
// Reads a cutlist from the specified file and converts it to an int-array containing the end of a cut per entry; the 
// length if this array is returned
// (those cutList-arrays are subject to the audioSegmentation()-function in SC_Lib that is called by the java 
//  frontend; this method os static because it doesn't need any object knowledge but needs to be called before 
//  audioSegmentation() if only the file-based cutList is available)
//====================================================================================================================
unsigned long int SC_GroundTruth_MPEG7::cutListFile2cutListArray(const char* fileName, int* &cutList) {
	FILE* inFile = fopen(fileName, "r");
	char start[sclib::bufferSize], end[sclib::bufferSize], length[sclib::bufferSize], type[sclib::bufferSize];
	unsigned long int count = 0;
	int scanStatus;

	if (inFile != NULL) {
		//count cuts in the file
		while (!feof(inFile)) {
			scanStatus = fscanf(inFile, "%s\t%s\t%s\t%s\n", start, end, length, type);
			if (sclib::isNum(start[0]) == true) {
				count++;
			}
		}

		//alocate memory for the cutList
		MArray_1D(cutList, count, int, "SC_GroundTruth_MPEG7.cutListFile2cutListArray: cutList");
		count = 0;
	  
		//extract the information form the file
		fseek(inFile, 0, SEEK_SET);
		count = 0;
		while (!feof(inFile)) {
			scanStatus = fscanf(inFile, "%s\t%s\t%s\t%s\n", start, end, length, type);
			
			if (sclib::isNum(start[0]) == true) {
				cutList[count++] = atoi(end); //for hard cuts, start=end; for gradual transitions, use the end als the first frame of the new scene/shot as in the scivo-groundtruth
			}
		}

		fclose(inFile);
	}

  return count;
}

//====================================================================================================================
// Methods to read in the results of the extern audio-/video-segmentation-module from an integer-array of cuts/scenes
// and initializes the frameList; the return value indicates success (true) or error (false)
//====================================================================================================================
bool SC_GroundTruth_MPEG7::readSceneList(int* sceneList, int sceneListLength) {
	unsigned long int sceneBoundary = 0, lastBoundary;
  bool failure = false;

	//set an artificial scene start at the first frame
	setFrame(0, sclib::atSceneBoundary|sclib::atShotBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
	setFrame(0, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);

	//extract the information form the file
	for (int x = 0; x < sceneListLength; x++) {
    lastBoundary = sceneBoundary;
		sceneBoundary = (unsigned long int)(sceneList[x]);
    sceneBoundary = videoFrame2FLI(sclib::max(0, sceneBoundary + this->pTweak->groundTruth.videoFrameMachineOffset), sclib::alignmentStart);
    
    //make sure the scene-list is correct and all boundarys are in consecutive order
    if (!(sceneBoundary > lastBoundary || sceneBoundary == 0)) {
      printf("\n Error in ground truth data: Scene boundary at video-frame %d !> %d", FLI2videoFrame(sceneBoundary, sclib::alignmentStart), FLI2videoFrame(lastBoundary, sclib::alignmentStart));
      failure = true;
    }  

		//remember cuts after the end of the multimedia file as reported by the audio-decoder (which somehow ofton reports too less frames)
		if (sceneBoundary >= this->internalFrameCount) {
			this->postEndCutCount++;
		}

    //copy the information into frameList
    if (sceneBoundary < this->internalFrameCount) {
      setFrame(sceneBoundary, sclib::atSceneBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth); //we read cut-detection lists, not scene-detection lists...
      setFrame(sceneBoundary, sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
      setFrame(sceneBoundary, sclib::atSceneBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
      setFrame(sceneBoundary, sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
    }
	}

  return (!failure);
}

//====================================================================================================================
//	convert video-frame-nr into internal frame-nr
//====================================================================================================================
unsigned long int SC_GroundTruth_MPEG7::videoFrame2FLI(unsigned long int frame, unsigned int alignment) {
  return this->pConverter->videoFrame2audioFrame(frame, this->internalFrameSize, this->internalFrameSize, alignment);
}

//====================================================================================================================
//	convert internal frame-nr into video-frame-nr
//====================================================================================================================
unsigned long int SC_GroundTruth_MPEG7::FLI2videoFrame(unsigned long int frameListIndex, unsigned int alignment) {
  return this->pConverter->audioFrame2videoFrame(frameListIndex, this->internalFrameSize, this->internalFrameSize, alignment);
}

//====================================================================================================================
//  To make operator<<() kind of virtual...
//====================================================================================================================
ostream& SC_GroundTruth_MPEG7::output(ostream& OutS, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types, unsigned long int typesNot, int origin, bool suppressHeader) {
  long int id, width = 10, precision = 3;
  unsigned long int x, start = sample2FLI(segmentStart), end = sample2FLI(segmentEnd);
  bool correct;
  char hypothesizedSpeakerName[sclib::bufferSize];
  
  OutS << setiosflags(ios_base::left|ios_base::fixed|ios_base::showpoint);

  if (suppressHeader == false) {
		OutS << this->audioFileName << endl;
		OutS << "FrameRate\t" << getVideoFrameRate() << endl;

		OutS << "FRAME#" << "\t"
				 << "VIDEO#" << "\t"
			   << SC_GroundTruth::getAudioTypeName(sclib::atNoise, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atNoise, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atSpeech, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeech, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atPause, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atPause, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atVoiced, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atVoiced, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atUnvoiced, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atUnvoiced, true) << "_p\t"
				 //<< "MALE" << "\t"
				 //<< "FEMALE" << "\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atSilence, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSilence, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atShort, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atShort, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atPureSpeech, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atPureSpeech, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atNoisySpeech, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atNoisySpeech, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atBackground, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atBackground, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atMusic, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atMusic, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atAction, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atAction, true) << "_p\t"
				 //<< "pBREATH" << "\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atUndefined, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atUndefined, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atSceneBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSceneBoundary, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atShotBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atShotBoundary, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atSpeakerBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeakerBoundary, true) << "_p\t"
				 << SC_GroundTruth::getAudioTypeName(sclib::atNoiseBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atNoiseBoundary, true) << "_p\t"
				 << "ART_B" << "\t"
         << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentStart, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentStart, true) << "_p\t"
         << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentEnd, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentEnd, true) << "_p\t"
				 << "SPK_ID" << "\t" << "SPK_ID_p" << "\t"
				 << "SPK"
				 << endl;
  }

  for (x = start; x <= ((end < this->internalFrameCount) ? end : this->internalFrameCount - 1); x++) {
    if (testFrame(x, types, false, typesNot, true, origin) == true) {
      id = getSpeakerGIDFromHID(this->frameList[x][3], correct); //the mapping has to be established earlier
      if (testFrame(x, sclib::atSpeech, false, sclib::noType, true, origin) == true && id >= 0) {
        if (correct == true) {
          sprintf(hypothesizedSpeakerName, "%s\0",  this->speakerNames[this->frameList[x][2]]);
        } else {
          sprintf(hypothesizedSpeakerName, "%s*\0",  this->speakerNames[this->frameList[x][2]]); //add an asterisc to the name if the hypothesized speaker is only fitting to the gt-speaker and not correct (i.e. not on the biggest cluster)
        }
      } else {
        sprintf(hypothesizedSpeakerName, "%s\0", "-");
      }

			OutS << x << "\t"
					 << FLI2videoFrame(x, sclib::alignmentStart) << "\t"
			     << ((this->frameList[x][1] & sclib::atNoise) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atNoise)] : 1.0) << "\t"
					 << ((this->frameList[x][1] & sclib::atSpeech) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atSpeech)] : 1.0) << "\t"
					 <<	((this->frameList[x][1] & sclib::atPause) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atPause)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atVoiced) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atVoiced)] : 1.0) << "\t"
					 << ((this->frameList[x][1] & sclib::atUnvoiced) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atUnvoiced)] : 1.0) << "\t"
           //<< ((this->frameList[x][1] & sclib::atMaleVoice) > 0) << " (" << ((this->frameList[x][0] & sclib::atMaleVoice) > 0) << ")\t"
           //<< ((this->frameList[x][1] & sclib::atFemaleVoice) > 0) << " (" << ((this->frameList[x][0] & sclib::atFemaleVoice) > 0) << ")\t"
					 << ((this->frameList[x][1] & sclib::atSilence) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atSilence)] : 1.0) << "\t"
					 << ((this->frameList[x][1] & sclib::atShort) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atShort)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atPureSpeech) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atPureSpeech)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atNoisySpeech) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atNoisySpeech)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atBackground) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atBackground)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atMusic) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atMusic)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atAction) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atAction)] : 1.0) << "\t"
           //<< ((this->frameList[x][1] & sclib::atBreath) > 0) << " (" << ((this->frameList[x][0] & sclib::atBreath) > 0) << ")\t"
           << ((this->frameList[x][1] & sclib::atUndefined) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atUndefined)] : 1.0) << "\t"
					 << ((this->frameList[x][1] & sclib::atSceneBoundary) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atSceneBoundary)] : 1.0) << "\t"
					 << ((this->frameList[x][1] & sclib::atShotBoundary) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atShotBoundary)] : 1.0) << "\t"
					 << ((this->frameList[x][1] & sclib::atSpeakerBoundary) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atSpeakerBoundary)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atNoiseBoundary) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atNoiseBoundary)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atArtificialBoundary) > 0) << " (" << ((this->frameList[x][0] & sclib::atArtificialBoundary) > 0) << ")\t"
					 << ((this->frameList[x][1] & sclib::atSpeechSegmentStart) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atSpeechSegmentStart)] : 1.0) << "\t"
           << ((this->frameList[x][1] & sclib::atSpeechSegmentEnd) > 0) << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][sclib::bitPosition(sclib::atSpeechSegmentEnd)] : 1.0) << "\t"
					 << this->frameList[x][3] << "\t" << setprecision(precision) << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? this->probabilityList[x][getProbabilityListDim()-1] : 1.0) << "\t"
           << (testFrame(x, sclib::atSpeech, false, sclib::noType, true, origin) && id >= 0 ? this->speakerNames[id] : "-")
					 << endl;
    }
	}

	return OutS;
}

//====================================================================================================================
// for a given segment, this function writes the ratio of brutto speech length (all frames labeled sclib::atSpeech) to
// netto speech length (all sclib::atSpeech-frames without sclib::atPause|sclib::atSilence|sclib::atUnvoiced) into a 
// file; also tells how many seconds of sppech where too short to analyze
//====================================================================================================================
void SC_GroundTruth_MPEG7::segmentStatisticsOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd) {
	unsigned long int bruttoLength = 0, nettoLength = 0, shortLength = 0, frame;
  long int speechStart, speechEnd;
	char *outLine	= new char[1024];
	
  for (unsigned long int y = segmentStart; y <= segmentEnd; y++) {
		getNextSegment(y, speechStart, speechEnd, sclib::atSpeech);
		if ((speechStart != sclib::noSegment) && (speechEnd != sclib::noSegment) && (speechEnd <= (long)(segmentEnd))) {

      //calculate the brutto- and netto-length of the speech segments
      for (frame = sample2FLI((unsigned long)speechStart); frame <= sample2FLI((unsigned long)speechEnd); frame++) {
		    if (testFrame(frame, sclib::atSpeech) == true) {
			    bruttoLength++;
			    if (!(testFrame(frame, sclib::atPause) == true) && !(testFrame(frame, sclib::atUnvoiced) == true)) {
				    nettoLength++;
			    }
          if (testFrame(frame, sclib::atShort) == true) {
            shortLength++;
          }
		    }
	    }
			
      //print it to a file
      sprintf(outLine, "start: %12u;\t end: %12u;\t brutto: %8.3f s; \tnetto: %8.3f s; \tratio: %7.3f %%; \ttoo short: %8.3f s\0", this->pConverter->sample2videoFrame(segmentStart, sclib::alignmentStart), this->pConverter->sample2videoFrame(segmentEnd, sclib::alignmentEnd), (double)this->pConverter->audioFrame2ms(bruttoLength, this->internalFrameSize, this->internalFrameSize, sclib::alignmentStart)/1000.0, (double)this->pConverter->audioFrame2ms(nettoLength, this->internalFrameSize, this->internalFrameSize, sclib::alignmentStart)/1000.0, ((double)nettoLength/(double)bruttoLength)*100, (double)this->pConverter->audioFrame2ms(shortLength, this->internalFrameSize, this->internalFrameSize, sclib::alignmentStart)/1000.0);
	    sclib::scalarOut(fileName, outLine, this->pTweak);
	      
      y = speechEnd;
		} //segment-boundarys are valid
		else {break;}
	} //for all frames in this scene

	MFree_1D(outLine);
	
	return;
}

//====================================================================================================================
// check the ground-truth (or hypothesized results) for consistency; this includes the check for mutually exclusive
// tupels of audio-types (e.g. no internal frame can be speech and noise together); returns true if everything is ok
// violations of rules are reported in a file or on screen (if fileName != "")
//====================================================================================================================
bool SC_GroundTruth_MPEG7::checkConsistency(unsigned long int segmentStart, unsigned long int segmentEnd, int origin) {
	unsigned long int y, startFrame, endFrame;
	bool res = true, existedBefore = false;
	ostream *OutS = &cerr;
	fstream fileOut;
	char *fName = NULL, postfixedFilename[sclib::bufferSize];
	long int t, problems, idx = (origin = sclib::modeGroundtruth) ? 0 : 1;
	int nrOfTypes = sizeof(this->frameList[0][0])*8 - 1;
	char typeNameBuffer[sclib::bufferSize];

	//prepare file-output instead of console output, if wished
	if (sclib::bitTest(pTweak->debug.debugMode, sclib::dbConsistencyCheck) == true) {
		sprintf(postfixedFilename, "%s%s%s\0", "inconsistencies", ((origin == sclib::modeGroundtruth) ? "_gt" : "_hypo"), ".txt");
		fName = new char[strlen(this->pTweak->debug.debugDir) + strlen(postfixedFilename) + 1];
		sprintf(fName, "%s%s\0", this->pTweak->debug.debugDir, postfixedFilename);
		existedBefore = sclib::fileExists(fName); //remember if this file existed before creation (to decidde if it can be deleted afterwards if nothing was written...)
		fileOut.open(fName, ios_base::out|ios_base::app);
		OutS = &fileOut;
	}

	startFrame = sample2FLI(segmentStart);
	endFrame = sample2FLI((segmentEnd < this->audioSampleCount) ? segmentEnd : this->audioSampleCount-1);

	if (origin == sclib::modeHypothesized) { //no groundtruth for MPEG7-files, so also no check
		for (y = startFrame; y < endFrame; y++) {
			//check for mutual exclusiveness
			t = 1;
			for (int i = 0; i < nrOfTypes; i++) { //biggest type: 2^nrOfTypes
				if (testFrame(y, t, false, sclib::noType, false, origin) == true) {
					if (areMutualExclusiveTypes(t, this->frameList[y][idx], problems) == true) { //test for mutual exclusiveness of the current type t (which is present in the current frame y) and any other type present in y
						*OutS << "FLI " << y << " (sample " << FLI2sample(y)  << "): " 
									<< "Mutual exclusiveness error with " << getAudioTypeName(t, true)
									<< " and " << getAudioTypesNames(problems, typeNameBuffer, ", ", true) << "\n";
						res = false;
					}
				}
				t = t << 1;
			}

			//v/uv result available?
			if (this->pTweak->segmentationHandler.vUvDetectorMode != sclib::algorithm_nothing && origin != sclib::modeGroundtruth) {
				if (testFrame(y, sclib::atSpeech, false, sclib::atVoiced|sclib::atUnvoiced|sclib::atPause, false, origin) == true) {
					*OutS << "FLI " << y << "(sample " << FLI2sample(y)  << "): No voiced/unvoiced-type classification result present\n";
					res = false;
				}
			}

			//ATC result available?
			if (this->pTweak->segmentationHandler.audioTypeMode != sclib::algorithm_nothing && origin != sclib::modeGroundtruth) {
				if (testFrame(y, (SC_GroundTruth::detailedAudioTypes|sclib::atSilence), false, sclib::noType, false, origin) == false) {
					*OutS << "FLI " << y << "(sample " << FLI2sample(y)  << "): No audio-type classification result present\n";
					res = false;
				}
			}
		}
	}

	if (sclib::bitTest(this->pTweak->debug.debugMode, sclib::dbConsistencyCheck) == true) {
		fileOut.close();
		if (res == true && existedBefore == false) {
			remove(fName); //try to delete the file if it didn't exist before and nothing has been written here
		}
	}

	MFree_1D(fName);

	return res;
}

//====================================================================================================================
// File-I/O for this class, so that the sate of a groundtruth-object can be saved to a file
// All members except the pSignalPrototype and the pTweak are saved
//====================================================================================================================
bool SC_GroundTruth_MPEG7::save(const char *fileName) {
	unsigned long int len;
	long int x;
 	bool res = true;
	int bytes;
	fstream gtFile;
	SC_GroundTruth::SC_SpeakerMapping *pHook = this->pSpeakerMapping;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes;
	io.getCurrentDatatypeSizes(codeSizes);

	gtFile.open(fileName, ios::out|ios::binary);  //truncate

	//simple types
	bytes = io.writeMachineHeader(&gtFile, codeSizes);
	bytes += io.writeScalar(&gtFile, this->gtType);
	bytes += io.writeScalar(&gtFile, this->uncertaintyRegion);
	bytes += io.writeScalar(&gtFile, this->internalFrameSize);
	bytes += io.writeScalar(&gtFile, this->internalFrameCount);
	bytes += io.writeScalar(&gtFile, this->audioSampleRate);
	bytes += io.writeScalar(&gtFile, this->audioSampleCount);
	bytes += io.writeScalar(&gtFile, this->postEndCutCount);
	bytes += io.writeMatrix(&gtFile, this->frameList, this->internalFrameCount, 4);

	//types where length-information is necessary before the actual data
	if (this->audioFileName != NULL) {
		len = (unsigned long int)(strlen(this->audioFileName));
		bytes += io.writeScalar(&gtFile, len);
		bytes += io.writeArray(&gtFile, this->audioFileName, len);
	} else {
		len = 0;
		bytes += io.writeScalar(&gtFile, len);
	}

	if (this->sceneFileName != NULL) {
		len = (unsigned long int)(strlen(this->sceneFileName));
		gtFile.write((char*)(&(len)), sizeof(unsigned long int));
		gtFile.write((char*)(this->sceneFileName), len * sizeof(char));
	} else {
		len = 0;
		gtFile.write((char*)(&(len)), sizeof(unsigned long int));
	}

	if (this->probabilityList != NULL) {
		len = getProbabilityListDim();
		bytes += io.writeScalar(&gtFile, len);
		bytes += io.writeMatrix(&gtFile, this->probabilityList, this->internalFrameCount, len);
	} else {
		len = 0;
		bytes += io.writeScalar(&gtFile, len);
	}

	len = sclib::maxSpeakers;
	bytes += io.writeScalar(&gtFile, len);
	for (x = 0; x < sclib::maxSpeakers; x++) {
		if (this->speakerNames[x] != NULL) {
			len = (unsigned long int)(strlen(this->speakerNames[x]));
			bytes += io.writeScalar(&gtFile, len);
			bytes += io.writeArray(&gtFile, this->speakerNames[x], len);
		} else {
			len = 0;
			bytes += io.writeScalar(&gtFile, len);
		}
	}
	
	len = sclib::getListCount(this->pSpeakerMapping);
	bytes += io.writeScalar(&gtFile, len);
	for (x = 0; x < (long int)(len); x++) {
		bytes += io.writeScalar(&gtFile, pHook->correct);
		bytes += io.writeScalar(&gtFile, pHook->groundTruthID);
		bytes += io.writeScalar(&gtFile, pHook->hypothesizedID);
		pHook = pHook->Next;
	}

  if (gtFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Saving SC_GroundTruth_MPEG7 Failed!");
		res = false;
	}

	gtFile.close();

  return res;
}

//====================================================================================================================
// File-I/O for this class, so that the sate of a groundtruth-object can be loaded from a file
// Everything gets loaded except the pSignalPrototype and the pTweak member
// 'this' is returned or NULL in case of error
//====================================================================================================================
SC_GroundTruth* SC_GroundTruth_MPEG7::load(const char *fileName) {
	unsigned long int len, len2;
	long int x, groundTruthID, hypothesizedID;
	bool correct;
	int bytes;
	fstream gtFile;
	SC_GroundTruth *pGT = this;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

	gtFile.open(fileName, ios::in|ios::binary);  //read

	//simple types
	bytes = io.readMachineHeader(&gtFile, fileSizes, true);
	if (bytes > 0) {
		io.consumeBytes(&gtFile, bytes);
	} else {
		bytes = 0;
	}
	bytes += io.readScalar(&gtFile, this->gtType, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->uncertaintyRegion, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->internalFrameSize, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->internalFrameCount, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->audioSampleRate, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->audioSampleCount, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->postEndCutCount, codeSizes, fileSizes);
	MFree_2D(this->frameList);
  MArray_2D(this->frameList, (long int)(this->internalFrameCount), 4, long int, "SC_GroundTruth.load: frameList"); 
	bytes += io.readMatrix(&gtFile, this->frameList, this->internalFrameCount, 4, codeSizes, fileSizes);

	//types where length-information is necessary before the actual data
	MFree_1D(this->audioFileName);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->audioFileName, len+1, char, "SC_GroundTruth.load: audioFileName");
		bytes += io.readArray(&gtFile, this->audioFileName, len, codeSizes, fileSizes);
		this->audioFileName[len] = '\0';
	}

	MFree_1D(this->sceneFileName);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->sceneFileName, len+1, char, "SC_GroundTruth.load: audioFileName");
		bytes += io.readArray(&gtFile, this->sceneFileName, len, codeSizes, fileSizes);
		this->audioFileName[len] = '\0';
	}

	MFree_2D(this->probabilityList);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_2D(this->probabilityList, (long int)(this->internalFrameCount), len, double, "SC_GroundTruth.load: probabilityList"); 
		bytes += io.readMatrix(&gtFile, this->probabilityList, this->internalFrameCount, len, codeSizes, fileSizes);
	}

	for (x = 0; x < sclib::maxSpeakers; x++) {
		MFree_1D(this->speakerNames[x]);
	}
	MFree_1D(this->speakerNames);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len != sclib::maxSpeakers) {
		REPORT_ERROR(SVLIB_BadData, "New max. speaker-name-count differs from old one!");
	}
	MArray_1D(this->speakerNames, len, char*, "SC_GroundTruth.load: speakerNames");
	for (x = 0; x < (long int)(len); x++) {
		bytes += io.readScalar(&gtFile, len2, codeSizes, fileSizes);
		if (len2 > 0) {
			MArray_1D(this->speakerNames[x], len2+1, char, "SC_GroundTruth.load: speakerNames[x]");
			bytes += io.readArray(&gtFile, this->speakerNames[x], len2, codeSizes, fileSizes);
			this->speakerNames[x][len2] = '\0';
		} else {
			this->speakerNames[x] = NULL;
		}
	}

	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	for (x = 0; x < (long int)(len); x++) {
		bytes += io.readScalar(&gtFile, correct, codeSizes, fileSizes);
		bytes += io.readScalar(&gtFile, groundTruthID, codeSizes, fileSizes);
		bytes += io.readScalar(&gtFile, hypothesizedID, codeSizes, fileSizes);
		addSpeakerMapping(groundTruthID, hypothesizedID, correct);
	}

	if (gtFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Loading SC_GroundTruth_MPEG7 Failed!");
		pGT = NULL;
	}

	gtFile.close();

	return pGT;
}

//====================================================================================================================
// Method to convert the findings of the algorithms (hypo-col of frameList, speaker-ids and probabiliy-list) to 
// video-frame based lists that can be used by videana; memory for the new lists is allocated herein
// The return value is the length of the 3 return-parameters
//====================================================================================================================
unsigned long int SC_GroundTruth_MPEG7::getResults(long int* &segmentationResults, long int* &speakerIDs, double** &probabilities) {
	unsigned long int videoFrameCount = sclib::round(this->audioSampleCount / this->videoFrameSize);
	long int currentType, segStart, segEnd, vfStart, vfEnd, bits = sizeof(this->frameList[0][0])*8 - 1; //-1 because the frameList is signed;
	unsigned long int dim = getProbabilityListDim();
	double maxProbability;

	//allocate memory
	MFree_1D(segmentationResults);
	MArray_1D(segmentationResults, videoFrameCount, long int, "SC_GroundTruth_MPEG7.getResults: segmentationResults");
	MFree_1D(speakerIDs);
	//if (this->pTweak->speakerClusterer.doClustering == true) { //otherwise there are no IDs
		MArray_1D(speakerIDs, videoFrameCount, long int, "SC_GroundTruth_MPEG7.getResults: segmentatpeakerIDs");
	//}
	MFree_2D(probabilities);
	if (this->pTweak->groundTruth.storeProbabilityInformation == true) { //otherwise there are no probabilities
		MArray_2D(probabilities, (long int)(videoFrameCount), dim, double, "SC_GroundTruth_MPEG7.getResults: probabilities");
	}

	for (unsigned long int vf = 0; vf < videoFrameCount; vf++) {
		vfStart = sclib::min(this->audioSampleCount-1, this->pConverter->videoFrame2sample(vf, sclib::alignmentStart)); //necessary to always compute it new here to not accumulate rounding errors due to odd frame-rates (NTSC)...
		vfEnd = sclib::min(this->audioSampleCount-1, this->pConverter->videoFrame2sample(vf, sclib::alignmentEnd));

		//care for segmentation results
		segmentationResults[vf] = sclib::noType;
		for (short int b = 0; b < bits; b++) { //loop over all types
			currentType = sclib::bit(b);
			if (testSegment(vfStart, vfEnd, true, currentType) > 0) { 
				//set a type as appearing in this frame if it appears min. 1 time rather than appearing min. 50% of the time
				//thus we avoid video-frames having no labels at all because they appear only less than 50%, boundaries for example wouldn't appear then at all
				//this way we may get video-frames having mutually exclusive labels, but that's ok 'cause that's what really happened in the signal
				segmentationResults[vf] |= currentType;
			}
		}

		//care for speaker-ids
		if (speakerIDs != NULL) {
			speakerIDs[vf] = getMajorSpeakerID(vfStart, vfEnd); //get most often occuring speaker-id (or sclib::noSpeaker in case of no speaker) for this frame
		} else {
			speakerIDs[vf] = sclib::noSpeaker; //videana expects a valid speakerID-array, so fill it with a valid value in case of no ids
		}

		//care for probabilities
		if (probabilities != NULL) {
			segStart = sclib::min(this->internalFrameCount-1, sample2FLI(vfStart));
			segEnd = sclib::min(this->internalFrameCount-1, sample2FLI(vfEnd));
			for (unsigned long int c = 0; c < dim; c++) { //loop over all columns of the probabilities-matrix
				maxProbability = std::numeric_limits<double>::max() * (-1);
				for (int i = segStart; i <= segEnd; i++) {
					if (this->probabilityList[i][c] > maxProbability) {
						maxProbability = this->probabilityList[i][c];
					}
				}
				probabilities[vf][c] = maxProbability;
			}
		}
	} //for vf

	return videoFrameCount;
}

//====================================================================================================================
//	output all analyzed parts of the frameList (really ALL parts are assumed to have been analyzed in the context of 
//  this class & task!); mimic the behaviour of the base class by simultaneous using getResults(), which converts FLIs 
//  to videoFrames
//====================================================================================================================
void SC_GroundTruth_MPEG7::analyzedOutEx(const char* fileName, sclib::OpenMode mode, long int selectionMask, bool selectSpkID) {
  //long int sceneStart, sceneEnd;
  //unsigned long int sceneNr = 1;
	fstream fileOut;
	char *fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
  //bool printedSomething = false;
	unsigned long int videoFrames;
	long int *segmentationResults = NULL, *speakerIDs = NULL;
	double** probabilities = NULL;

	sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, mode);
	MFree_1D(fName);

  /*for (unsigned long int y = 0; y < this->audioSampleCount; y++) {
    getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary);
    if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
      if ((sceneNr < pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & pTweak->general.sceneSelection) && (pTweak->general.sceneSelection != 0))) {
        sceneNr++; 
        y = sceneEnd; 
        continue;
      }*/

			videoFrames = getResults(segmentationResults, speakerIDs, probabilities);
      
			output(fileOut, videoFrames, 0, segmentationResults, speakerIDs, probabilities, false, "\t", selectionMask, selectSpkID);
			
			MFree_1D(segmentationResults);
			MFree_1D(speakerIDs);
			MFree_2D(probabilities);

			/*printedSomething = true;
      sceneNr++;
      y = sceneEnd;
      if (sceneNr > pTweak->general.lastScene) {
				break;
			}
    }
  }*/

  fileOut.close();

	return;
}

//====================================================================================================================
//  To make operator<<() kind of virtual in a 2nd form: with the results returned by getResults()
//====================================================================================================================
ostream& SC_GroundTruth_MPEG7::output(ostream& OutS, unsigned long int frameCount, unsigned long int firstFrameNr, long int *segmentationResult, long int *speakerIDs, double **probabilities, bool suppressHeader, const char *separator, long int selectionMask, bool selectSpkID) {
  long int id;
  unsigned long int x;
	bool correct, speakers = speakerIDs != NULL && selectSpkID == true;
  
  OutS << setiosflags(ios_base::left|ios_base::fixed|ios_base::showpoint);

  if (suppressHeader == false) {
		OutS << this->audioFileName << endl;
		OutS << "FrameRate\t" << getVideoFrameRate() << endl;

		//OutS << "FRAME#" << "\t";
    OutS << "videoFrame" << separator;
		if (sclib::bitTest(selectionMask, sclib::atNoise) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atNoise, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atNoise, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeech) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atSpeech, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atSpeech, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atPause) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atPause, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atPause, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atVoiced) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atVoiced, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atVoiced, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atUnvoiced) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atUnvoiced, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atUnvoiced, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atMaleVoice) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atMaleVoice, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atMaleVoice, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atFemaleVoice) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atFemaleVoice, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atFemaleVoice, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSilence) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atSilence, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atSilence, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atShort) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atShort, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atShort, true) << "_p\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atPureSpeech) == true) {
		OutS << SC_GroundTruth::getAudioTypeName(sclib::atPureSpeech, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atPureSpeech, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atNoisySpeech) == true) {
		OutS << SC_GroundTruth::getAudioTypeName(sclib::atNoisySpeech, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atNoisySpeech, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atBackground) == true) {
		OutS << SC_GroundTruth::getAudioTypeName(sclib::atBackground, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atBackground, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atMusic) == true) {
		OutS << SC_GroundTruth::getAudioTypeName(sclib::atMusic, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atMusic, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atAction) == true) {
		OutS << SC_GroundTruth::getAudioTypeName(sclib::atAction, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atAction, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atBreath) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atBreath, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atBreath, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atUndefined) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atUndefined, true) << separator << SC_GroundTruth::getAudioTypeName(sclib::atUndefined, true) << "_p" << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSceneBoundary) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atSceneBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSceneBoundary, true) << "_p\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atShotBoundary) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atShotBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atShotBoundary, true) << "_p\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeakerBoundary) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atSpeakerBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeakerBoundary, true) << "_p\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atNoiseBoundary) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atNoiseBoundary, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atNoiseBoundary, true) << "_p\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atArtificialBoundary) == true) {
			OutS << "ART_B" << "\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeechSegmentStart) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentStart, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentStart, true) << "_p\t";
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeechSegmentEnd) == true) {
			OutS << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentEnd, true) << "\t" << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentEnd, true) << "_p\t";
		}
		if (speakers == true) {
			OutS << "SPK_ID" << separator << "SPK_ID_p" << separator;
			OutS << "SPK";
		}
		OutS << endl;
  }

  for (x = 0; x < frameCount; x++) {
		if (speakers == true) {
			id = getSpeakerGIDFromHID(speakerIDs[x], correct); //the mapping has to be established earlier
		}

		//OutS << videoFrame2audioFrame(firstFrameNr+x, this->internalFrameSize, this->internalFrameSize) << "\t";
		OutS << firstFrameNr + x << separator;
		if (sclib::bitTest(selectionMask, sclib::atNoise) == true) {
			OutS << ((segmentationResult[x] & sclib::atNoise) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atNoise)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeech) == true) {
			OutS << ((segmentationResult[x] & sclib::atSpeech) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atSpeech)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atPause) == true) {
			OutS <<	((segmentationResult[x] & sclib::atPause) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atPause)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atVoiced) == true) {
			OutS << ((segmentationResult[x] & sclib::atVoiced) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atVoiced)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atUnvoiced) == true) {
			OutS << ((segmentationResult[x] & sclib::atUnvoiced) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atUnvoiced)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atMaleVoice) == true) {
			OutS << ((segmentationResult[x] & sclib::atMaleVoice) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atMaleVoice)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atFemaleVoice) == true) {
			OutS << ((segmentationResult[x] & sclib::atFemaleVoice) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atFemaleVoice)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSilence) == true) {
			OutS << ((segmentationResult[x] & sclib::atSilence) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atSilence)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atShort) == true) {
			OutS << ((segmentationResult[x] & sclib::atShort) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atShort)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atPureSpeech) == true) {
			OutS << ((segmentationResult[x] & sclib::atPureSpeech) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atPureSpeech)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atNoisySpeech) == true) {
			OutS << ((segmentationResult[x] & sclib::atNoisySpeech) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atNoisySpeech)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atBackground) == true) {
			OutS << ((segmentationResult[x] & sclib::atBackground) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atBackground)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atMusic) == true) {
			OutS << ((segmentationResult[x] & sclib::atMusic) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atMusic)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atAction) == true) {
			OutS << ((segmentationResult[x] & sclib::atAction) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atAction)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atBreath) == true) {
			OutS << ((segmentationResult[x] & sclib::atBreath) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atBreath)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atUndefined) == true) {
			OutS << ((segmentationResult[x] & sclib::atUndefined) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atUndefined)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSceneBoundary) == true) {
			OutS << ((segmentationResult[x] & sclib::atSceneBoundary) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atSceneBoundary)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atShotBoundary) == true) {
			OutS << ((segmentationResult[x] & sclib::atShotBoundary) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atShotBoundary)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeakerBoundary) == true) {
			OutS << ((segmentationResult[x] & sclib::atSpeakerBoundary) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atSpeakerBoundary)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atNoiseBoundary) == true) {
			OutS << ((segmentationResult[x] & sclib::atNoiseBoundary) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atNoiseBoundary)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atArtificialBoundary) == true) {
			OutS << ((segmentationResult[x] & sclib::atArtificialBoundary) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atArtificialBoundary)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeechSegmentStart) == true) {
			OutS << ((segmentationResult[x] & sclib::atSpeechSegmentStart) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atSpeechSegmentStart)] : 1.0) << separator;
		}
		if (sclib::bitTest(selectionMask, sclib::atSpeechSegmentEnd) == true) {
			OutS << ((segmentationResult[x] & sclib::atSpeechSegmentEnd) > 0) << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][sclib::bitPosition(sclib::atSpeechSegmentEnd)] : 1.0) << separator;
		}
		if (speakers == true) {
			OutS << speakerIDs[x] << separator << ((this->pTweak->groundTruth.storeProbabilityInformation == true) ? probabilities[x][getProbabilityListDim()-1] : 1.0) << separator;
			OutS << (((segmentationResult[x] & sclib::atSpeech) > 0 && id >= 0) ? this->speakerNames[id] : "-");
		}
		OutS << endl;
	}

	return OutS;
}

//====================================================================================================================
// returns the video-framerate (in frames per second)
//====================================================================================================================
double SC_GroundTruth_MPEG7::getVideoFrameRate(void) {
	return (double)(this->audioSampleRate) / this->videoFrameSize;
}
