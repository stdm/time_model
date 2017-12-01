/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth do accomplish its targets in the    */
/*        context of SCiVo-Corpus files (Video's and annotation)          */
/*      - reads and organizes and provides acces to SCiVo annotation/     */
/*        ground-truth files                                              */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.02.2004																								*/
/**************************************************************************/

#include <list>
#include <math.h>
#include <limits>
#include <string.h>
#include <iomanip>
#include "SC_GroundTruth_SCiVo.h"
#include "SC_Aux.h"
#include "SC_Signal_MPEG.h"
#include <SV_DataIO.h>

/**************************************************************************/
/*       	                Constructor, Destructor     									  */
/**************************************************************************/

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth_SCiVo::SC_GroundTruth_SCiVo(SC_TweakableParameters *pTweak, double videoFrameRate, const char* audioFileName, const char* sceneFileName, const char* segmentFileName) : SC_GroundTruth(pTweak, audioFileName) {
  //some init-stuff
  this->gtType = sclib::gtSCiVo;

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
	this->uncertaintyRegion = sclib::round(2 * this->videoFrameSize * 3); //twice because of +-

  this->segmentFileName = new char[strlen(segmentFileName)+1];
	strcpy(this->segmentFileName, segmentFileName);

  this->sceneFileName = new char[strlen(sceneFileName)+1];
	strcpy(this->sceneFileName, sceneFileName);

  //read the ground-truth files
  readGroundTruth();
}

//====================================================================================================================
//	destructor 
//====================================================================================================================
SC_GroundTruth_SCiVo::~SC_GroundTruth_SCiVo() {
  MFree_1D(this->segmentFileName);
  MFree_1D(this->sceneFileName);
}

//====================================================================================================================
//	initialize the framelist with the data of the scene-/segment-List
//====================================================================================================================
bool SC_GroundTruth_SCiVo::readGroundTruth(void) {
  bool generateScenes = false, success = true;
  long int res;

	//initialization of the frameList with zeros for the labels and -1 for the speaker-id's has been already done by initFrameList()

  //read the sceneList, copy it's information into the frameList
  res = readSceneFile(this->sceneFileName);
  if (res == SVLIB_BadData) {
    REPORT_ERROR(SVLIB_BadData, "Error in scene ground-truth data\n");
    success = false;
  } else if (res == 0) {
     //a scenefile isn't necessary, we can generate pseudo scenes instead :-)
    generateScenes = true;
  }
  
	//copy the information of segmentList to frameList
  res = readSegmentFile(this->segmentFileName);

  if (res == SVLIB_BadData) {
    REPORT_ERROR(SVLIB_BadData, "\nError in audio ground truth data\n");
    success = false;
  } else if (res == 0) {
    REPORT_ERROR(SVLIB_BadData, "\nNo Data in Segmentation-File\n");
    success = false;
  }

  //generate pseudo scene boundarys if there was no scene file
  if (generateScenes == true) {
    generatePseudoSceneBoundarys(this->pConverter->ms2sample(this->pTweak->groundTruth.pseudoSceneLength));
  }

  return success;
}

//====================================================================================================================
//	read scene-list from file to internal data-structure
//====================================================================================================================
long int SC_GroundTruth_SCiVo::readSceneFile(const char* filename) {
	FILE* inFile = fopen(filename, "r");
  char* buffer = new char[sclib::bufferSize];
	unsigned long int count = 0, sceneBoundary = 0, lastBoundary;
  bool failure = false;

	if (inFile == NULL) {return 0;}

	//extract the information form the file
	while (!feof(inFile)) {
    sclib::readline(inFile, buffer, sclib::bufferSize);
    if (sclib::isNum(buffer[0]) == true) {
      count++;
      lastBoundary = sceneBoundary;
			sceneBoundary = sclib::getNextIntFromString(buffer, sclib::bufferSize);
      sceneBoundary = videoFrame2FLI(sclib::max(0, sceneBoundary + this->pTweak->groundTruth.videoFrameMachineOffset), sclib::alignmentStart);
      
      //make sure the scene-list is correct and all boundaries are in consecutive order
      //assert(sceneBoundary > lastBoundary || sceneBoundary == 0);
      if (!(sceneBoundary > lastBoundary || sceneBoundary == 0)) {
        printf("\n Error ground truth data: Scene boundary at video-frame %d !> %d", FLI2videoFrame(sceneBoundary, sclib::alignmentStart), FLI2videoFrame(lastBoundary, sclib::alignmentStart));
        failure = true;
      }  

      //copy the information into frameList
      if (sceneBoundary < this->internalFrameCount) {
        setFrame(sceneBoundary, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
        setFrame(sceneBoundary, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
      }
		}
	}

  fclose(inFile);
	MFree_1D(buffer);

  return (failure == true) ? SVLIB_BadData : (long int)count;
}

//====================================================================================================================
//	read segmentation-list from file to internal data-structure
//====================================================================================================================
long int SC_GroundTruth_SCiVo::readSegmentFile(const char* filename) {
	FILE *inFile = fopen(filename, "r");
  char *buffer = new char[sclib::bufferSize], *comment = NULL, *speakerName = NULL, *oldSpeakerName = NULL;
	unsigned long int x, count = 0, segmentBoundary = 0, startFrame = 0, endFrame, oldBoundary = 0;
	int speakerID, oldSpeakerID = sclib::noSpeaker;
	int segmentType, oldSegmentType, olderSegmentType;
  bool generateSpeakerBoundarys = true, failure = false;
	const int segSpeech = 0, segNonSpeech = 1; //values as used in the segment-files

	if (inFile == NULL) {return 0;}

	//extract the information from the file
	while (!feof(inFile)) {
    sclib::readline(inFile, buffer, sclib::bufferSize);
    if (sclib::isNum(buffer[0]) == true) {

			//extract one line, fill the right varaiables, do value conversion
      segmentBoundary = sclib::getNextIntFromString(buffer, sclib::bufferSize); 
			segmentType = sclib::getNextIntFromString(buffer, sclib::bufferSize);
      segmentType = (segmentType == segSpeech) ? sclib::atSpeech : sclib::atNoise;
      comment = sclib::extractComment(buffer, sclib::bufferSize);
      
			//initializations, if this is the first loop
			if (count == 0) { 
        oldSegmentType = segmentType;
        olderSegmentType = segmentType;
      }
			
			//make sure the audio-list is correct and all boundarys are in consecutive order
      if (!(segmentBoundary > oldBoundary || segmentBoundary == 0)) { 
        printf("\n Error in ground truth data: Segment boundary at video-frame %d !> %d", segmentBoundary, oldBoundary);
        failure = true;
      }

      //compute the real endframe
      endFrame = videoFrame2FLI(sclib::max(0, segmentBoundary + this->pTweak->groundTruth.videoFrameMachineOffset), sclib::alignmentEnd); //endFrame is the start of the current segment and the end of the last one
      if (endFrame >= this->internalFrameCount) {
        endFrame = this->internalFrameCount - 1;
      }

      //is there a new speaker-name for the current segment given?
      if (sclib::bitTest(segmentType, sclib::atSpeech) == true) {
        if ((comment != NULL) && (strncmp(comment, "\0", sclib::bufferSize) > 0)) {
          MFree_1D(oldSpeakerName);
          oldSpeakerName = speakerName; //not just the speaker-name of the last segment (which could be "" in case of no change or non-speech) but the last valid speaker name given
          speakerName = comment;
          insertSpeakerName(speakerName, sclib::bufferSize);
					if ((oldSpeakerName == NULL) || (strncmp(oldSpeakerName, speakerName, sclib::bufferSize) != 0)) { //different names => different speakers
						setFrame(endFrame, sclib::atSpeakerBoundary|((oldSpeakerName == NULL) ? sclib::atArtificialBoundary : sclib::noType), false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
						setFrame(endFrame, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
						generateSpeakerBoundarys = false;
 					}
        } else {
          MFree_1D(comment);
        }
      } else {
        MFree_1D(comment);
      }

      //set the last segment according to the  oldSegmentType, set the oldSpeakerID if it was a speech-segment
      speakerID = getSpeakerIDfromName(speakerName, sclib::bufferSize);
		  for (x = startFrame; x < endFrame; x++) {
        setFrame(x, oldSegmentType, false, oldSpeakerID, sclib::modeLabelAdd, sclib::modeGroundtruth);
        setFrame(x, oldSegmentType, false, oldSpeakerID, sclib::modeLabelAdd, sclib::modeHypothesized);
				//artificaially create groundtruth for sub-classes of speech and noise
				setFrame(x, (sclib::bitTest(oldSegmentType, sclib::atSpeech) == true) ? sclib::atNoisySpeech : sclib::atBackground, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
		  }

      //assume a speaker-boundary direct after a scene-boundary if the speech is continous
      if ((sclib::bitTest(oldSegmentType, sclib::atSpeech) == true) && (testFrame(startFrame, sclib::atSceneBoundary) == true) && (sclib::bitTest(segmentType, sclib::atSpeech) == true)) {
        setFrame(startFrame, sclib::atSpeakerBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
        setFrame(startFrame, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
      }

      //prepare for the next loop
		  startFrame = endFrame;
		  olderSegmentType = oldSegmentType;
      oldSegmentType = segmentType;
			oldSpeakerID = speakerID;
      oldBoundary = segmentBoundary;

			count++;
		}
	}
  MFree_1D(speakerName);
  MFree_1D(oldSpeakerName);
	MFree_1D(buffer);
	fclose(inFile);

  //generate artificial speaker changes after each non-speech segment if none where given in the ground-truth
  if (generateSpeakerBoundarys == true) {
		//assume a speaker-boundary within each new speech-segment after a nonspeech-segment 
    for (x = 1; x < this->internalFrameCount; x++) {
      if ((testFrame(x, sclib::atSpeech) == true) && (testFrame(x-1, sclib::atSpeech) == false)) {
			  setFrame(x, sclib::atSpeakerBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
        setFrame(x, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
		  }
    }
  }

  //set the very last segment of the file
  for (x = startFrame; x < this->internalFrameCount; x++) {
		setFrame(x, oldSegmentType, false, ((sclib::bitTest(oldSegmentType, sclib::atSpeech) == true) ? speakerID : sclib::noSpeaker), sclib::modeLabelAdd, sclib::modeGroundtruth);
    setFrame(x, oldSegmentType, false, ((sclib::bitTest(oldSegmentType, sclib::atSpeech) == true) ? speakerID : sclib::noSpeaker), sclib::modeLabelAdd, sclib::modeHypothesized);
		//artificaially create groundtruth for sub-classes of speech and noise
		setFrame(x, (sclib::bitTest(oldSegmentType, sclib::atSpeech) == true) ? sclib::atNoisySpeech : sclib::atBackground, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
	}

  //set a speaker-boundary if the first frame of the file is a speech frame
  if (testFrame(0, sclib::atSpeech) == true) {
		setFrame(0, sclib::atSpeakerBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
    setFrame(0, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
  }

  //set a speaker-boundary at the beginning of each first speech-segment in a new scene
  generateSpeakerBoundarys = false;
  for (x = 1; x < this->internalFrameCount; x++) { //x=1,because we just tested the first frame in the block above!
    if (testFrame(x, sclib::atSceneBoundary) == true) {
      generateSpeakerBoundarys = true;
    }
    if ((testFrame(x, sclib::atSpeech) == true) && (generateSpeakerBoundarys == true)) {
      setFrame(x, sclib::atSpeakerBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
      setFrame(x, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
      generateSpeakerBoundarys = false;
		}
  }

  return (failure == true) ? SVLIB_BadData : (long int)count;
}

//====================================================================================================================
//	convert video-frame-nr into internal frame-nr
//====================================================================================================================
unsigned long int SC_GroundTruth_SCiVo::videoFrame2FLI(unsigned long int frame, unsigned int alignment) {
  return this->pConverter->videoFrame2audioFrame(frame, this->internalFrameSize, this->internalFrameSize, alignment);
}

//====================================================================================================================
//	convert internal frame-nr into video-frame-nr
//====================================================================================================================
unsigned long int SC_GroundTruth_SCiVo::FLI2videoFrame(unsigned long int frameListIndex, unsigned int alignment) {
  return this->pConverter->audioFrame2videoFrame(frameListIndex, this->internalFrameSize, this->internalFrameSize, alignment);
}

//====================================================================================================================
//  To make operator<<() kind of virtual...
//====================================================================================================================
ostream& SC_GroundTruth_SCiVo::output(ostream& OutS, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types, unsigned long int typesNot, int origin, bool suppressHeader) {
  long int id, width = 3;
  unsigned long int x, start = sample2FLI(segmentStart), end = sample2FLI(segmentEnd);
  bool correct;
  char hypothesizedSpeakerName[sclib::bufferSize];

  OutS << setiosflags(ios_base::left|ios_base::fixed|ios_base::showpoint);

  if (suppressHeader == false) {
		OutS << setw(10)				<< "FRAME#"
         << setw(10)				<< "VIDEO#"
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atNoise, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atSpeech, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atPause, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atVoiced, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atUnvoiced, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atMaleVoice, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atFemaleVoice, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atSilence, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atShort, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atPureSpeech, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atNoisySpeech, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atBackground, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atMusic, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atAction, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atBreath, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atUndefined, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atSceneBoundary, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atShotBoundary, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atSpeakerBoundary, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atNoiseBoundary, true)
				 << setw(width+2)		<< SC_GroundTruth::getAudioTypeName(sclib::atArtificialBoundary, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentStart, true)
				 << setw(2*width+2) << SC_GroundTruth::getAudioTypeName(sclib::atSpeechSegmentEnd, true)
				 << setw(2*width+2)	<< "SPK_ID"
				 << setw(2*width+2) << "SPK"
				 << endl;
  }

  for (x = start; x <= ((end < this->internalFrameCount) ? end : this->internalFrameCount - 1); x++) {
    if (testFrame(x, types, false, typesNot, false, origin) == true) {
      id = getSpeakerGIDFromHID(this->frameList[x][3], correct); //the mapping has to be established earlier
      if (testFrame(x, sclib::atSpeech, false, sclib::noType, true, origin) == true && id >= 0) {
        if (correct == true) {
          sprintf(hypothesizedSpeakerName, "%s\0",  this->speakerNames[id]);
        } else {
          sprintf(hypothesizedSpeakerName, "%s*\0",  this->speakerNames[id]); //add an asterisc to the name if the hypothesized speaker is only fitting to the gt-speaker and not correct (i.e. not on the biggest cluster)
        }
      } else {
        sprintf(hypothesizedSpeakerName, "%s\0", "-");
      }

			OutS << setw(10) << x
           << setw(10) << FLI2videoFrame(x, sclib::alignmentStart)
			     << setw(width) << ((this->frameList[x][1] & sclib::atNoise) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atNoise) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atSpeech) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atSpeech) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atPause) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atPause) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atVoiced) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atVoiced) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atUnvoiced) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atUnvoiced) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atMaleVoice) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atMaleVoice) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atFemaleVoice) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atFemaleVoice) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atSilence) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atSilence) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atShort) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atShort) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atPureSpeech) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atPureSpeech) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atNoisySpeech) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atNoisySpeech) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atBackground) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atBackground) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atMusic) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atMusic) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atAction) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atAction) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atBreath) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atBreath) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atUndefined) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atUndefined) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atSceneBoundary) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atSceneBoundary) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atShotBoundary) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atShotBoundary) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atSpeakerBoundary) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atSpeakerBoundary) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atNoiseBoundary) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atNoiseBoundary) > 0) << ")"
           << "(" << setw(width) << ((this->frameList[x][0] & sclib::atArtificialBoundary) > 0) << ")"
					 << setw(width) << ((this->frameList[x][1] & sclib::atSpeechSegmentStart) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atSpeechSegmentStart) > 0) << ")"
           << setw(width) << ((this->frameList[x][1] & sclib::atSpeechSegmentEnd) > 0) << "(" << setw(width) << ((this->frameList[x][0] & sclib::atSpeechSegmentEnd) > 0) << ")"
					 << setw(2*width+2) << this->frameList[x][3]
           << setw(15) << hypothesizedSpeakerName << "(" << (testFrame(x, sclib::atSpeech, false, sclib::noType, true, sclib::modeGroundtruth) && this->frameList[x][2] >= 0 ? this->speakerNames[this->frameList[x][2]] : "-") << ")"
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
void SC_GroundTruth_SCiVo::segmentStatisticsOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd) {
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
// File-I/O for this class, so that the sate of a groundtruth-object can be saved to a file
// All members except the pSignalPrototype and the pTweak are saved
//====================================================================================================================
bool SC_GroundTruth_SCiVo::save(const char *fileName) {
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
	bytes += io.writeMatrix(&gtFile, this->frameList, this->internalFrameCount, 4);
	bytes += io.writeScalar(&gtFile, this->videoFrameSize);

	//types where length-information is necessary before the actual data
	if (this->audioFileName != NULL) {
		len = (unsigned long int)(strlen(this->audioFileName));
		bytes += io.writeScalar(&gtFile, len);
		bytes += io.writeArray(&gtFile, this->audioFileName, len);
	} else {
		len = 0;
		bytes += io.writeScalar(&gtFile, len);
	}

	if (this->segmentFileName != NULL) {
		len = (unsigned long int)(strlen(this->segmentFileName));
		bytes += io.writeScalar(&gtFile, len);
		bytes += io.writeArray(&gtFile, this->segmentFileName, len);
	} else {
		len = 0;
		bytes += io.writeScalar(&gtFile, len);
	}

	if (this->sceneFileName != NULL) {
		len = (unsigned long int)(strlen(this->sceneFileName));
		bytes += io.writeScalar(&gtFile, len);
		bytes += io.writeArray(&gtFile, this->sceneFileName, len);
	} else {
		len = 0;
		bytes += io.writeScalar(&gtFile, len);
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
		REPORT_ERROR(SVLIB_Fail, "Saving SC_GroundTruth_SCiVo Failed!");
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
SC_GroundTruth* SC_GroundTruth_SCiVo::load(const char *fileName) {
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
	MFree_2D(this->frameList);
  MArray_2D(this->frameList, (long int)(this->internalFrameCount), 4, long int, "SC_GroundTruth.load: frameList"); 
	bytes += io.readMatrix(&gtFile, this->frameList, this->internalFrameCount, 4, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->videoFrameSize, codeSizes, fileSizes);

	//types where length-information is necessary before the actual data
	MFree_1D(this->audioFileName);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->audioFileName, len+1, char, "SC_GroundTruth.load: audioFileName");
		bytes += io.readArray(&gtFile, this->audioFileName, len, codeSizes, fileSizes);
		this->audioFileName[len] = '\0';
	}

	MFree_1D(this->segmentFileName);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->segmentFileName, len+1, char, "SC_GroundTruth.load: audioFileName");
		bytes += io.readArray(&gtFile, this->segmentFileName, len, codeSizes, fileSizes);
		this->segmentFileName[len] = '\0';
	}

	MFree_1D(this->sceneFileName);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->sceneFileName, len+1, char, "SC_GroundTruth.load: audioFileName");
		bytes += io.readArray(&gtFile, this->sceneFileName, len, codeSizes, fileSizes);
		this->sceneFileName[len] = '\0';
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
		REPORT_ERROR(SVLIB_Fail, "Loading SC_GroundTruth_SCiVo Failed!");
		pGT = NULL;
	}

	gtFile.close();

	return pGT;
}

//====================================================================================================================
// returns the video-framerate (in frames per second)
//====================================================================================================================
double SC_GroundTruth_SCiVo::getVideoFrameRate(void) {
	return (double)(this->audioSampleRate) / this->videoFrameSize;
}
