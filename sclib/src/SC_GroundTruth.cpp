/**************************************************************************/
/*    Responsibility:																											*/
/*		  - organizes file-io (audio-stream, scene-list, segement-list)			*/
/*			- organizes timeline in terms of video-frames											*/
/*			- maps seconds and audio-samples to video-frames									*/
/*			- privides access to the audio-segments														*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.02.2004																								*/
/**************************************************************************/

#include <list>
#include <map>
#include <math.h>
#include <limits>
#include <string.h>
#include <iomanip>
#include <assert.h>
#include "SC_GroundTruth.h"
#include "SC_Aux.h"
#include "SC_SignalHandler.h"
#include <SV_DataIO.h>

/**************************************************************************/
/*       	                Constructor, Destructor     									  */
/**************************************************************************/

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth::SC_GroundTruth(SC_TweakableParameters *pTweak, const char* audioFileName) {
	SC_SignalHandler *pLoader = NULL;
  
  //some init-stuff
  this->gtType = sclib::gtStandard;
  this->pTweak = pTweak;
	this->frameList = NULL;
	MArray_1D(this->speakerNames, sclib::maxSpeakers, char*, "SC_GroundTruth: speakerNames");
  for (unsigned int x = 0; x < sclib::maxSpeakers; x++) {
    this->speakerNames[x] = NULL;
  }
  this->pSpeakerMapping = NULL;
	this->uncertaintyRegion = 0;
	this->probabilityList = NULL;
	this->pConverter = new SC_Conversion();

  //maybe a derived class (as the TIMIT-groundtruth) has to handle this a little bit different...
  if (audioFileName != NULL && strncmp(audioFileName, "", sclib::bufferSize) != 0) {
    this->audioFileName = new char[strlen(audioFileName)+1];
	  strcpy(this->audioFileName, audioFileName);

    //open the signal to have a prototype of it's parameters
    pLoader = new SC_SignalHandler(this->pTweak, sclib::stGuess);
    this->pSignalPrototype = pLoader->openSignal(audioFileName, this->pTweak->signalHandler.forceSampleRate);
    MFree_0D(pLoader);

    this->audioSampleRate = this->pSignalPrototype->SigPar.SRate;
		this->pConverter->setAudioSampleRate(this->audioSampleRate);
    this->audioSampleCount = this->pSignalPrototype->getSampleCount();
	  this->internalFrameSize = this->pConverter->ms2sample(this->pTweak->groundTruth.internalFrameSize); //in samples!!!
    this->internalFrameCount = (unsigned long int)(ceil((double)(this->audioSampleCount) / (double)(this->internalFrameSize)));

    //init the frame-list
    initFrameList();
  }
}

//====================================================================================================================
//	destructor 
//====================================================================================================================
SC_GroundTruth::~SC_GroundTruth() {
	MFree_2D(this->frameList);
	MFree_1D(this->audioFileName);
  MFree_0D(this->pSignalPrototype);
	MFree_2D(this->probabilityList);
	MFree_0D(this->pConverter);
  
  sclib::destructLinkedList(this->pSpeakerMapping);
  
  for (unsigned int x = 0; x < sclib::maxSpeakers; x++) {
    MFree_1D(this->speakerNames[x]);
  }
	MFree_1D(this->speakerNames);
}

//====================================================================================================================
//	initialize the framelist
//====================================================================================================================
void SC_GroundTruth::initFrameList(void) {
  unsigned long int probabilityListDim;
	
	MArray_2D(this->frameList, (long int)(this->internalFrameCount), 4, long int, "initFrameList: frameList"); 

	probabilityListDim = getProbabilityListDim();
	if (this->pTweak->groundTruth.storeProbabilityInformation == true) {
		MArray_2D(this->probabilityList, (long int)(this->internalFrameCount), probabilityListDim, double, "initFrameList: probabilityList"); 
	}

	//init frameList with zeros for the labels and -1 for the speaker-id's
  for (unsigned long int y = 0; y < this->internalFrameCount; y++) {
		this->frameList[y][0] = sclib::noType;  //ground-truth frame-attributes
    this->frameList[y][1] = sclib::noType;  //hypothesized frame-attributes
		this->frameList[y][2] = sclib::noSpeaker; //ground-truth speaker-id
    this->frameList[y][3] = sclib::noSpeaker; //hypothesized speaker-id
		if (this->pTweak->groundTruth.storeProbabilityInformation == true) {
			for (unsigned int x = 0; x < probabilityListDim; x++) {
				this->probabilityList[y][x] = 1.0;
			}
		}
	}

	return;
}

//====================================================================================================================
//  Inserts a speaker-name into the speakerNames[]-array and returns it's speaker-id (index into the array). If the
//  name already exists, the existing index is returned; if anything goes wrong (e.g. speakerName is NULL or '\0' or 
//  MAX_SPEAKERS speakers are reached), -1 (sclib::noSpeaker) is returned
//====================================================================================================================
int SC_GroundTruth::insertSpeakerName(const char* speakerName, int size) {
  int speakerID = sclib::noSpeaker;
  char* newName;
  
  if ((speakerName != NULL) && (strncmp(speakerName, "\0", size) >= 0)) {
    for (int x = 0; x < sclib::maxSpeakers; x++) {
      if (this->speakerNames[x] != NULL) { //check if this name already exists; if so, return it's index
        if (strncmp(this->speakerNames[x], speakerName, size) == 0) {
          speakerID = x;
          break;
        }
      } else { //insert as a new name, return the next free index
        newName = new char[size+1];
				sprintf(newName, "%s", speakerName);
        this->speakerNames[x] = newName;
        speakerID = x;
        break;
      }
    }
  }

  return speakerID;
}

//====================================================================================================================
//  Returns the index (= ground-truth speaker-id) of this speaker-name in into the speakerNames[]-array; 
//  returns -1 (sclib::noSpeaker) if the given speaker-name is not found
//====================================================================================================================
int SC_GroundTruth::getSpeakerIDfromName(const char* speakerName, int size) {
  int speakerID = sclib::noSpeaker;
  
  if ((speakerName != NULL) && (strncmp(speakerName, "\0", size) >= 0)) {
    for (int x = 0; x < sclib::maxSpeakers; x++) {
      if (this->speakerNames[x] != NULL && strncmp(this->speakerNames[x], speakerName, size) == 0) {
        speakerID = x;
        break;
      }
    }
  }

  return speakerID;
}

//====================================================================================================================
//  Returns the speaker-name for a given ground-truth speaker-id; returns "" if the given speaker-id is unknown
//====================================================================================================================
const char* SC_GroundTruth::getSpeakerName(int speakerGID) {
	if (speakerGID >= 0 && speakerGID < sclib::maxSpeakers) {
		return this->speakerNames[speakerGID];
	} else {
		return "";
	}
}

//====================================================================================================================
// returns the numer of different names in the speakerNames[]-array; this corresponds to the nr. of different 
// speakers if ground-truthing is carryed out the right way
// counts only those speakers which appear in the scenes to evaluate: >= firstScene, <= lastScene, in sceneSelection
// count only those whose segments where long enough to model/cluster
// count only those segments between start & end (if end > 0, otherwise all)
//====================================================================================================================
unsigned int SC_GroundTruth::getRealSpeakerCount(bool considerSegmentLengthThresholds, unsigned long int start, unsigned long int end) {
  unsigned int count = 0;
  unsigned long int x, y, speakerID, sceneNr = 1, msPerGaussian, modelThreshold, clusterThreshold, segLength, stop, surroundingSegLength;
  long int sceneStart, sceneEnd, speechStart, speechEnd;
  std::list<unsigned long int> speakerIDs;

	msPerGaussian = this->pConverter->ms2sample(this->pTweak->modelHandler.msPerGaussian); //this corresponds with complete speaker-segments, i.e. all single speech-segments between consecutive speaker-boundaries
	modelThreshold = this->pConverter->ms2sample(this->pTweak->general.shortSpeechThreshold); //this corresponds with single speech-segments
	clusterThreshold = this->pConverter->ms2sample(this->pTweak->speakerClusterer.speechSegLengthThreshold); //this also corresponds with complete speaker-segments

	stop = (end > 0 && end < getAudioSampleCount()) ? end : getAudioSampleCount()-1;

  for (y = 0; y <= stop; y++) { //go from the beginning to see all scenes, even if speakers before "start" are not counted
    //go through all scenes
		getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward, sclib::modeGroundtruth); 
    if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment && sceneStart <= (long int)(stop)) {
      //is this scene to evaluate?
      if ((sceneNr <= this->pTweak->general.lastScene) &&(sceneNr >= this->pTweak->general.firstScene) && ((sclib::bit(sceneNr) & this->pTweak->general.sceneSelection) || (this->pTweak->general.sceneSelection == 0))) {
        for (x = (unsigned long)sceneStart; x <= (unsigned long)sceneEnd; x++) {
          //go to all speech-segments
					getNextSegment(x, speechStart, speechEnd, sclib::atSpeech, sclib::searchForward, sclib::modeGroundtruth, true, true, sclib::atSpeakerBoundary|sclib::atSceneBoundary); 
          if ((speechStart != sclib::noSegment) && (speechEnd != sclib::noSegment) && (speechStart <= sceneEnd)) {
						//is the speech-segment in the range to evaluate?
						if (sclib::intersect(start, stop, speechStart, speechEnd) > 0) {
							//is it long enough?
							segLength = (unsigned long int)(speechEnd - speechStart + 1);
							surroundingSegLength = getSurroundingSpeakersSegmentSize(speechStart, speechEnd);
							if ((considerSegmentLengthThresholds == false) || (
									(segLength > modelThreshold || modelThreshold == 0) && 
									(surroundingSegLength > clusterThreshold || clusterThreshold == 0) &&
									(surroundingSegLength > msPerGaussian))) {
								//is only one speaker to evaluate?
								if (strncmp(this->pTweak->modelHandler.onlyThisSpeaker, "", sclib::bufferSize) == 0 || 
										isSpeakerHomogenious(speechStart, speechEnd, this->pTweak->modelHandler.onlyThisSpeaker) == true) {
									//insert each occuring speakerID in the list
									speakerID = getSpeakerID(sample2FLI(speechStart), sclib::modeGroundtruth);
									speakerIDs.push_back(speakerID);
								}
							}
						}
          } else {
            break;
          }
          x = speechEnd;
        }
      }
      sceneNr++;
      y = sceneEnd;
      if (sceneNr > this->pTweak->general.lastScene) {
				break;
			}
    } else {
			break;
		}
  }

  //find out how many differnet speaker-ID's there are in the list
  speakerIDs.sort();
  speakerIDs.unique();
  count = (unsigned int)speakerIDs.size();

  return count;
}

//====================================================================================================================
// Returns the number of all samples having the sclib::atSpeech label within the speaker-boundary the given segment is 
// belonging to (needed if the segment-length-threshold for clustering should be evaluated: A new cluster contains all 
// speech belonging to one speaker-segment, so not the length of a single speech-segment has to be compared to the
// segLengthThreshold, but all speech in it's speakers segment).
// Here, origin==sclib::modeGroundtruth makes perfect sense as a default value 'cause the method is meant to be used
// to forecast the number of speakers/ce's to be recognized by the algorithms depending on the groundtruth
// ATTENTION: Relys on the fact that a new cluster always contains the speech of one speaker's segment, i.e. it is 
//            created by SC_ModelHandler.buildSpeakerModels()
//====================================================================================================================
unsigned long int SC_GroundTruth::getSurroundingSpeakersSegmentSize(unsigned long int segmentStart, unsigned long int segmentEnd, int origin) {
	unsigned long int res = 0;
	long int boundaryStart, boundaryEnd;

	getNextBoundary(segmentStart, boundaryStart, boundaryEnd, sclib::atSpeakerBoundary, sclib::searchMiddle, origin);
	if (boundaryStart != sclib::noSegment && boundaryEnd != sclib::noSegment) {
		res = testSegment(boundaryStart, boundaryEnd, false, sclib::atSpeech, false, sclib::noType, false, origin);
	} else {
		res = segmentEnd - segmentStart + 1; //this should never happen, though...
	}

	return res;
}

//====================================================================================================================
// counts the number of scenes to evaluate according to firstScene, lastScene and sceneSelection
//====================================================================================================================
unsigned long int SC_GroundTruth::getRealSceneCount(void) {
  long int sceneStart, sceneEnd;
	unsigned long int count = 0, sceneNr = 1;

  for (unsigned long int y = 0; y <= FLI2sample(this->internalFrameCount-1, sclib::alignmentEnd); y++) {
    //go through all scenes
		getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary, sclib::searchForward, sclib::modeGroundtruth); 
    if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
      //is this scene to evaluate?
      if ((sceneNr <= this->pTweak->general.lastScene) &&(sceneNr >= this->pTweak->general.firstScene) && ((sclib::bit(sceneNr) & this->pTweak->general.sceneSelection) || (this->pTweak->general.sceneSelection == 0))) {
				count++;
			}
      sceneNr++;
      y = sceneEnd;
      if (sceneNr > this->pTweak->general.lastScene) {break;}
    } else {break;}
  }

	return count;
}

//====================================================================================================================
// method to generate pseudo scene boundarys at suitable positions
// this is helpful if no scene-file is available to read real detected scene boundarys, but if one wishes to part the
// audio into little pieces anyway (to reduce computational load and memory requirements)
// this function analyzes the framelist and generates scene-boundarys approximately all 'approximateSceneLength' 
// ms, but cares for not cutting homogenious speech segments into pieces. care also that each scene has some noise 
// (max. 250ms after a extension because of a connected speech segment) for background estimation at it's start and 
// end; approximateSceneLength is meant to be in samples
// ATTENTION: The generated boundarys are set for both the groundtruth- and hypothesized data
//====================================================================================================================
unsigned long int SC_GroundTruth::generatePseudoSceneBoundarys(unsigned long int approximateSceneLength) {
  unsigned long int sceneLength = sample2FLI(approximateSceneLength);
  unsigned long int sceneBoundary = sceneLength, sceneCount = 0;
	long int speechStart, speechEnd;

  //remove all previously set sclib::atSceneBoundary labels
  setSegment(0, FLI2sample(this->internalFrameCount, sclib::alignmentEnd), sclib::atSceneBoundary, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeGroundtruth);
  setSegment(0, FLI2sample(this->internalFrameCount, sclib::alignmentEnd), sclib::atSceneBoundary, false, sclib::noSpeaker, sclib::modeLabelRemove, sclib::modeHypothesized);

	//boundary at the first frame
  setFrame(0, sclib::atSceneBoundary|sclib::atShotBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
  setFrame(0, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);

	if (sceneLength > 0) {
		while (sceneBoundary < this->internalFrameCount) {
			//don't part an existing, connected sclib::atSpeech segment; step to it's end
			if (testFrame(sceneBoundary, sclib::atSpeech) == true) {
				getNextSegment(FLI2sample(sceneBoundary, sclib::alignmentStart), speechStart, speechEnd, sclib::atSpeech, sclib::searchMiddle);
				sceneBoundary = speechEnd + 1; //there mus be a segment because the current frame has the speech-label; therefore, no check for valid borders here...
	      
				//search for the next speech-segment start to give this "scene" max. 250ms of noise for background estimation at it's end.
				getNextSegment(sceneBoundary, speechStart, speechEnd, sclib::atSpeech, sclib::searchForward); //note that sceneBoundary is sample-based here due to the bove line...
				if (speechStart != sclib::noSegment && speechEnd != sclib::noSegment) {
					sceneBoundary += sclib::max((speechStart-sceneBoundary)/2, this->pConverter->ms2sample(250, sclib::alignmentEnd));
				}

				sceneBoundary = sample2FLI(sceneBoundary);
			}

			//set the found boundary
			if (sceneBoundary < this->internalFrameCount) {
				setFrame(sceneBoundary, sclib::atSceneBoundary|sclib::atShotBoundary|sclib::atArtificialBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
				setFrame(sceneBoundary, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
				sceneBoundary += sceneLength;
				sceneCount++;
			}
		}
	}

	return sceneCount;
}

/**************************************************************************/
/*				       methods for dealing with the frameList										*/
/**************************************************************************/

//====================================================================================================================
//	Method to get the start and end of a boundary (scene, speaker, ...) in terms of samples out of the frameList
//  the behaviour of this method is dependant on the flag 'direction':
//    - sclib::searchForward:  starting from 'initPosition', the next start- and endpoint of a segment of the 
//                             desired type is returned (the standart case)
//    - sclib::searchBackward: searching from 'initPosition' backward, the last segement completely before the 
//                             'initPosition' is returned
//    - sclib::searchMiddle:   it is searched for a startpoint before- and a endpoint after 'initPosition'; if a 
//                             boundary lies on the initposition, the following segment (starting at 'initPosition') 
//                             is returned
//	ATTENTION:	boundaries in the frameList are only marked by their beginning, the start of the next one means the 
//              end+1 of the last one!
//====================================================================================================================
void	SC_GroundTruth::getNextBoundary(unsigned long int initPosition, long int &segmentStart, long int &segmentEnd, long int type, unsigned int direction, int origin) {
	long int x = 0;

	segmentStart = sclib::noSegment;
	segmentEnd = sclib::noSegment;

  if (type == sclib::atSpeakerBoundary || type == sclib::atSceneBoundary || type == sclib::atNoiseBoundary || type == sclib::atShotBoundary || type == sclib::atArtificialBoundary) {
    switch (direction) {
      case sclib::searchWithin: //the same as sclib::searchForward, so no break here!!!
      case sclib::searchForward: {
	      for (x = sample2FLI(initPosition); x < (long)this->internalFrameCount; x++) { //get boundary-start
		      if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
			      segmentStart = x;
			      break;
		      }
	      }
	      for (x = segmentStart + 1; x < (long)this->internalFrameCount; x++) { //get boundary-end
		      if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
			      segmentEnd = x - 1;
			      break;
		      } else if (x == this->internalFrameCount - 1) {
			      segmentEnd = x;
			      break;
		      }
	      }
        if (segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
					segmentStart = sclib::max(initPosition, FLI2sample(segmentStart, sclib::alignmentStart));
					segmentEnd = sclib::max(segmentStart, FLI2sample(segmentEnd, sclib::alignmentEnd));
        }
        break;
      } //sclib::searchForward
      
      case sclib::searchBackward: {
	      for (x = sample2FLI(initPosition); x >= 0; x--) { //get boundary-end
		      if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
			      segmentEnd = x - 1;
			      break;
		      }
	      }
	      for (x = segmentEnd - 1; x >= 0; x--) { //get boundary-start
		      if ((testFrame(x, type, true, sclib::noType, true, origin) == true) || (x == 0)) {
			      segmentStart = x;
			      break;
		      }
	      }
        if (segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
					segmentEnd = sclib::min(initPosition, FLI2sample(segmentEnd, sclib::alignmentEnd));
					segmentStart = sclib::min(segmentStart, FLI2sample(segmentStart, sclib::alignmentStart));
        }
        break;
      } //sclib::searchBackward

      case sclib::searchMiddle: {
	      for (x = sample2FLI(initPosition); x >= 0; x--) { //get boundary-start
		      if ((testFrame(x, type, true, sclib::noType, true, origin) == true) || (x == 0)) {
			      segmentStart = x;
			      break;
		      }
	      }
 	      for (x = segmentStart + 1; x < (long)this->internalFrameCount; x++) { //get boundary-end
		      if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
			      segmentEnd = x - 1;
			      break;
		      } else if (x == this->internalFrameCount - 1) {
			      segmentEnd = x;
			      break;
		      }
        }
        if (segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
          segmentStart = FLI2sample(segmentStart, sclib::alignmentStart);
          segmentEnd = FLI2sample(segmentEnd, sclib::alignmentEnd);
        }
        break;
      } //sclib::searchMiddle

      default: break;
    }
  } else { //type is a known boundary-type 
    REPORT_ERROR(SVLIB_BadArg, "Boundary-type unknown");
  }
	
  return;
}

//====================================================================================================================
//	Method to get the start and end of a segment (speech, noise, silence, ...) in terms of samples out of the 
//  framelist. The behaviour of this method is dependant on the flag 'direction':
//    - sclib::searchForward:  starting from 'initPosition', the next start- and endpoint of a segment of the 
//                             desired type is returned (the standart case)
//    - sclib::searchBackward: searching from 'initPosition' backward, the last segement completely before the 
//                             'initPosition' is returned
//    - sclib::searchMiddle:   it is searched for a startpoint before- and a endpoint after 'initPosition'; if a 
//                             boundary lies on the initposition, the following segment (starting at 
//                             'initPosition') is returned
//    - sclib::searchWithin:   as sclib::searchForward, but if the startFrame already has the desired 
//                             segmentType, it is returned as the segmentStart negelecting the fact that the
//                             segment started earlier
//  If the speaker-id changes within a segment and speakerBordersSegment==true, the changepoint is also considered as 
//  the segment-end; the same is for other boundaries (speaker, scene or shot, may be OR-concatenated in the OR sense)
//  when boundaryBordersSegment==true;
//	ATTENTION:	segments are marked by continous labeling; e.g. if one frame is labeled as silence and the next is 
//              not, this means the end of this segment!
//====================================================================================================================
void SC_GroundTruth::getNextSegment(unsigned long int initPosition, long int &segmentStart, long int &segmentEnd, long int type, unsigned int direction, int origin, bool speakerBordersSegment, bool boundaryBordersSegment, unsigned long int boundaryType) {
	long int speakerID;

	segmentStart = sclib::noSegment;
  segmentEnd = sclib::noSegment;

  switch (direction) {
  
    case sclib::searchForward: {
	    for (unsigned long int x = sample2FLI(initPosition); x < this->internalFrameCount; x++) {
				//get segment-start
		    if (segmentStart == sclib::noSegment) { 
					if (testFrame(x, type, true, sclib::noType, true, origin) == true) { //frame must have the wanted type
						speakerID = getSpeakerID(x, origin);
						if (x == 0) { //frame may be the first one...
							segmentStart = x;
						} else if (testFrame(x-1, type, true, sclib::noType, true, origin) == false) { //...or the type just changed to the desired one...
							segmentStart = x;
						} else if (boundaryBordersSegment == true && testFrame(x, boundaryType, false, sclib::noType, true, origin) == true) { //...or it is the beginning of a new boundary...
							segmentStart = x;
						} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x-1, origin)) { //... or the beginning of a new speaker's segment
					    segmentStart = x;
						}
			    }
		    }
				
				//get segment-end
		    if (segmentStart != sclib::noSegment) { 
					if (testFrame(x, type, true, sclib::noType, false, origin) == true) { //frame must have the desired type
						if (x == this->internalFrameCount - 1) { //frame may be the very last one...
							segmentEnd = x;
						} else if (testFrame(x+1, type, true, sclib::noType, true, origin) == false) { //...or frame has the type, but it's successor hasn't...
							segmentEnd = x;
						} else if (boundaryBordersSegment == true && testFrame(x+1, boundaryType, false, sclib::noType, true, origin) == true) { //...or it is the beginning of a new boundary...
							segmentEnd = x; 
						} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x+1, origin)) { //...or frame's successor has a new speaker
							segmentEnd = x;
						}
					}
		    }

		    if ((segmentStart != sclib::noSegment) && (segmentEnd != sclib::noSegment)) {
					segmentStart = sclib::max(initPosition, FLI2sample(segmentStart, sclib::alignmentStart));
					segmentEnd = sclib::max(segmentStart, FLI2sample(segmentEnd, sclib::alignmentEnd));
          break;
		    }
	    } //for x...
      break;
    } //sclib::searchForward

    case sclib::searchBackward: {
	    for (unsigned long int x = sample2FLI(initPosition); x >= 0; x--) {
				//get segment-end
		    if (segmentEnd == sclib::noSegment) {
					if (testFrame(x, type, true, sclib::noType, false, origin) == true) {
						speakerID = getSpeakerID(x, origin);
						if (x == this->internalFrameCount - 1) {
							segmentEnd = x;
						} else if (testFrame(x+1, type, true, sclib::noType, false, origin) == false) {
							segmentEnd = x;
						} else if (boundaryBordersSegment == true && testFrame(x+1, boundaryType, false, sclib::noType, true, origin) == true) {
							segmentEnd = x;
						} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x+1, origin)) {
							segmentEnd = x;
						}
					}
		    }

				//get segment-start
		    if (segmentEnd != sclib::noSegment) { 
					if (testFrame(x, type, true, sclib::noType, false, origin) == true) {
						if (x == 0) {
							segmentStart = x;
						} else if (testFrame(x-1, type, true, sclib::noType, true, origin) == false) {
							segmentStart = x;
						} else if (boundaryBordersSegment == true && testFrame(x, boundaryType, false, sclib::noType, true, origin) == true) {
							segmentStart = x;
						} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x, origin)) {
							segmentStart = x;
						}
					}
		    }

		    if ((segmentStart != sclib::noSegment) && (segmentEnd != sclib::noSegment)) {
					segmentEnd = sclib::min(initPosition, FLI2sample(segmentEnd, sclib::alignmentEnd));
					segmentStart = sclib::min(segmentEnd, FLI2sample(segmentStart, sclib::alignmentStart));
			    break;
		    }
	    } //for x...
      break;
    } //sclib::searchBackward

    case sclib::searchMiddle: {
			//get segment-start
	    for (unsigned long int x = sample2FLI(initPosition); x > 0; x--) { 
				if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
					speakerID = getSpeakerID(x, origin);					
					if (x == 0) {
						segmentStart = FLI2sample(x, sclib::alignmentStart);
	          break;
					} else if (testFrame(x-1, type, true, sclib::noType, true, origin) == false) {
						segmentStart = FLI2sample(x, sclib::alignmentStart);
		        break;
					} else if (boundaryBordersSegment == true && testFrame(x, boundaryType, false, sclib::noType, true, origin) == true) {
						segmentStart = FLI2sample(x, sclib::alignmentStart);
			      break;
					} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x-1, origin)) {
						segmentStart = FLI2sample(x, sclib::alignmentStart);
				    break;
					}
			  }
	    }

			//get segment-end
      for (unsigned long int x = sample2FLI(initPosition + 1); x < this->internalFrameCount; x++) { 
				if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
					if (x == this->internalFrameCount - 1) {
						segmentEnd = FLI2sample(x, sclib::alignmentEnd);
						break;
					} else if (testFrame(x+1, type, true, sclib::noType, true, origin) == false) {
						segmentEnd = FLI2sample(x, sclib::alignmentEnd);
						break;
					} else if (boundaryBordersSegment == true && testFrame(x+1, boundaryType, false, sclib::noType, true, origin) == true) {
						segmentEnd = FLI2sample(x, sclib::alignmentEnd);
						break;
					} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x+1, origin)) {
						segmentEnd = FLI2sample(x, sclib::alignmentEnd);
						break;
					}
				}
		  }
      break;
    } //sclib::searchMiddle

    case sclib::searchWithin: {
	    for (unsigned long int x = sample2FLI(initPosition); x < this->internalFrameCount; x++) {
				//get segment-start
		    if (segmentStart == sclib::noSegment) { 
			    if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
				    segmentStart = x;
						speakerID = getSpeakerID(x, origin);
			    }
		    }

				//get segment-end
		    if (segmentStart != sclib::noSegment) { 
					if (testFrame(x, type, true, sclib::noType, true, origin) == true) {
						if (x == this->internalFrameCount - 1) {
							segmentEnd = x;
						} else if (testFrame(x+1, type, true, sclib::noType, true, origin) == false) {
							segmentEnd = x;
						} else if (boundaryBordersSegment == true && testFrame(x+1, boundaryType, false, sclib::noType, true, origin) == true) {
							segmentEnd = x;
						} else if (speakerBordersSegment == true && speakerID != getSpeakerID(x+1, origin)) {
							segmentEnd = x;
						}
					}
		    }

		    if ((segmentStart != sclib::noSegment) && (segmentEnd != sclib::noSegment)) {
					segmentStart = sclib::max(initPosition, FLI2sample(segmentStart, sclib::alignmentStart));
					segmentEnd = sclib::max(segmentStart, FLI2sample(segmentEnd, sclib::alignmentEnd));
			    break;
		    }
	    } //for x...
      break;
    } //sclib::searchWithin

    default: break;
  } //switch
	
  return;
}

//====================================================================================================================
//  If uniteTypes==true, this method just calls getNextSegment() and returns it's result; new functionality is added 
//  for the other case: Then, the closest segment borders (minimum distance between appropriate segment-border and 
//  initPosition) having one of the types in 'types' is returned together with it's type; it is assumed that the given
//  types only appear mutually exclusive, a fact that can also be checked for if checkMutualTypeExclusiveness==true
//====================================================================================================================
unsigned long int SC_GroundTruth::getClosestSegment(unsigned long int initPosition, long int &segmentStart, long int &segmentEnd, long int types, bool uniteTypes, unsigned int direction, int origin, bool speakerBordersSegment, bool boundaryBordersSegment, unsigned long int boundaryType, bool checkMutualTypeExclusiveness) {
	unsigned long int foundType = sclib::noType;
	long int segStart[sizeof(this->frameList[0][0])*8], segEnd[sizeof(this->frameList[0][0])*8], min = std::numeric_limits<long int>::max();
	long int currentType, x, y, bits = sizeof(this->frameList[0][0])*8 - 1; //-1 because the frameList is signed
	char buffer[sclib::bufferSize] = "\0";

	segmentStart = sclib::noSegment;
	segmentEnd = sclib::noSegment;

	if (uniteTypes == true) {
		getNextSegment(initPosition, segmentStart, segmentEnd, types, direction, origin, speakerBordersSegment, boundaryBordersSegment, boundaryType);
		if (segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
			foundType = types;
		}
	} else {
		//assume mutual exclusiveness of the given types and find and return the segment of this type that lies closest to the initPosition
		for (x = 0; x < bits; x++) {
			segStart[x] = sclib::noSegment; //just an initialization for the loop below over y<x
			segEnd[x] = sclib::noSegment; 
      currentType = sclib::bit(x);
      if (sclib::bitTest(types, currentType)) {
			  getNextSegment(initPosition, segStart[x], segEnd[x], currentType, direction, origin, speakerBordersSegment, boundaryBordersSegment, boundaryType);

			  //remember segment closest to the initPosition; the definition of close differs with search-direction: if searched forward or
			  //within, it is bound to the minimum distance of the segmentStart to the initPosition; if searched backward, it is bound to the
			  //minimum distance of the segmentEnd 
			  if (segStart[x] != sclib::noSegment && segEnd[x] != sclib::noSegment) {
				  if (direction == sclib::searchForward || direction == sclib::searchMiddle || direction ==sclib::searchWithin) {
					  if (abs(sclib::max(segStart[x], initPosition)-sclib::min(segStart[x], initPosition))+1 < min) {
						  min = abs(sclib::max(segStart[x], initPosition) - sclib::min(segStart[x], initPosition)) + 1;
						  segmentStart = segStart[x];
						  segmentEnd = segEnd[x];
						  foundType = currentType;
					  }
				  } else {
					  if (segEnd[x]-(long)(initPosition)+1 < min) {
						  min = segEnd[x] - initPosition + 1;
						  segmentStart = segStart[x];
						  segmentEnd = segEnd[x];
						  foundType = currentType;
					  }
				  }

					//check assumption of mutual exclusiveness of wanted types; warn if it is violated
					if (checkMutualTypeExclusiveness == true) {
						for (y = 0; y < x; y++) {
							if ((segStart[y] != sclib::noSegment && segEnd[y] != sclib::noSegment) && //compare only with those types for which segment-borders where found
									((segStart[y]  >= segStart[x] && segStart[y] <= segEnd[x]) || //start of seg. y lies in the range of seg. x
										(segEnd[y]	 >= segStart[x] && segEnd[y]	 <= segEnd[x]) || //end   of seg. y lies in the range of seg. x
										(segStart[x] >= segStart[y] && segStart[x] <= segEnd[y]) || //start of seg. x lies in the range of seg. y
										(segEnd[x]	 >= segStart[y] && segEnd[x]	 <= segEnd[y]))) { //end   of seg. x lies in the range of seg. y
								sprintf(buffer, "%s%d%s%d%s%d%s%d%s%\0", "SC_GroundTruh.getClosestSegment(): Violation of mutual exclusiveness of types between segments x (", segStart[x], " - ", segEnd[x], ") and y (", segStart[y], " - ", segEnd[y], ")!");
								REPORT_ERROR(SVLIB_BadData, buffer);
							}
						}
					} //mutual exclusiveness check
			  } //segment is valid

      } //type is wanted
		} //loop over all possible types
	} //andTypes==false

	return foundType;
}

//====================================================================================================================
//  returns the next occasion of the given type in the given direction; helpful for e.g. sclib::atSpeechSegment*
//====================================================================================================================
long int SC_GroundTruth::getNextOccasion(unsigned long int initPosition, long int type, unsigned int direction, int origin, unsigned int alignment) {
	long int occasion = sclib::noSegment;

  switch (direction) {
    case sclib::searchForward:
		case sclib::searchMiddle:
		case sclib::searchWithin:
	    for (unsigned long int x = sample2FLI(initPosition); x < this->internalFrameCount; x++) {
				if (testFrame(x, type, true, sclib::noType, true, origin) == true) { //frame must have the wanted type
					occasion = x;
					break;
		    }
	    }
			if (occasion != sclib::noSegment) {
				occasion = sclib::max(initPosition, FLI2sample(occasion, alignment));
	    }
      break;

    case sclib::searchBackward:
	    for (unsigned long int x = sample2FLI(initPosition); x >= 0; x--) {
				if (testFrame(x, type, true, sclib::noType, false, origin) == true) {
					occasion = x;
					break;
				}
	    }
			if (occasion != sclib::noSegment) {
				occasion = sclib::min(initPosition, FLI2sample(occasion, alignment));
	    }
			break;

    default: break;
  } //switch
	
  return occasion;
}

//====================================================================================================================
//	Method to label a segment (given by start- and end-point in samples) with specific marker(s) (types) by processing
//  each frame in it with setFrame().
//====================================================================================================================
void SC_GroundTruth::setSegment(unsigned long int segmentStart, unsigned long int segmentEnd, long int types, bool uniteTypes, int speakerID, int action, int origin, bool forceNoSpeaker, bool playByTheRules) {
  unsigned long int segStart = sample2FLI(segmentStart); //convert samples to FLIs, use only FLIs below
  unsigned long int segEnd = sample2FLI(segmentEnd);
	bool unite = (uniteTypes==true || sclib::isPowerOfTwo(types)==true) ? true : false; //this saves computational time in case of a single type in setFrame()

  if (origin != sclib::modeHypothesized && origin != sclib::modeGroundtruth) {
    REPORT_ERROR(SVLIB_BadArg, "Decide whether to use groundtruth- or hypothesized data!");
  }

  if (segStart >= this->internalFrameCount)	{
    segStart = this->internalFrameCount - 1;
  }
	if (segEnd >= this->internalFrameCount)		{
    segEnd = this->internalFrameCount - 1;
  }

	for (unsigned long int x = segStart; x <= segEnd; x++) {
		setFrame(x, types, unite, speakerID, action, origin, forceNoSpeaker, playByTheRules);
	}

	return;
}

//====================================================================================================================
//	Method to label a segment (given by start- and end-point in samples) with specific marker(s) (type(s));
//  only those frames in the segment get processed by setFrame() which have the desired types; the number of actually 
//  processed (by setFrame()) frames is returned
//====================================================================================================================
unsigned long int SC_GroundTruth::setSegmentIf(unsigned long int segmentStart, unsigned long int segmentEnd, long int ifTypes, bool uniteIfTypes, long int ifTypesNot, bool uniteIfTypesNot, long int types, bool uniteTypes, int speakerID, int action, int origin, bool forceNoSpeaker, bool playByTheRules) {
  unsigned long int res = 0;
  unsigned long int segStart = sample2FLI(segmentStart); //convert samples to FLIs, use only FLIs below
  unsigned long int segEnd = sample2FLI(segmentEnd);

  if (origin != sclib::modeHypothesized && origin != sclib::modeGroundtruth) {
    REPORT_ERROR(SVLIB_BadArg, "Decide whether to use groundtruth- or hypothesized data!");
  }

  if (segStart >= this->internalFrameCount)	{
    segStart = this->internalFrameCount - 1;
  }
	if (segEnd >= this->internalFrameCount)		{
    segEnd = this->internalFrameCount - 1;
  }

	for (unsigned long int x = segStart; x <= segEnd; x++) {
		if (testFrame(x, ifTypes, uniteIfTypes, ifTypesNot, uniteIfTypesNot, origin) == true) {
			res++;
			setFrame(x, types, uniteTypes, speakerID, action, origin, forceNoSpeaker, playByTheRules);
		}
	}

	return res;
}

//====================================================================================================================
//	Method to label the single frame FLI with specific marker(s) (type(s))
//   - types: signle audioType-constant (SCLIB_AT_*) or OR-concatenated list of them
//   - uniteTypes: if true and the labelsshould be removed, they are removed only if they all apply to the frame, 
//                 otherwise each one is removed separately
//	 - speakerID: the frame is given this speaker-id if it contains SPEECH; if speakerID==sclib::noSpeaker, [...]
//   - forceNoSpeaker: [...] has to be true to set (or better: remove) it anyhow
//   - action: sclib::modeLabelRemove or sclib::modeLabelAdd
//   - orighin: sclib::modeGroundtruth or sclib::modeHypothesized
//   - playByTheRules: if true, a set of heuristics is applied to make the action not violate any rule/assumption 
//                     associated with the organisation of labels in the framelist; e.g. for label-adding it is 
//                     checked then if the given types are consistent with mutual exclusiveness rules one by one, and 
//                     already attached labels get removed if problems exist.
//====================================================================================================================
void SC_GroundTruth::setFrame(unsigned long int FLI, long int types, bool uniteTypes, int speakerID, int action, int origin, bool forceNoSpeaker, bool playByTheRules) {
  unsigned short idxAT = (origin == sclib::modeHypothesized) ? 1 : 0;
  unsigned short idxID = (origin == sclib::modeHypothesized) ? 3 : 2;
	long int t = 1, nrOfTypes = (sizeof(this->frameList[0][0])*8) - 1;
	long int problems;
  
  if (FLI >= this->internalFrameCount)	{
    FLI = this->internalFrameCount - 1;
  }

	if (action == sclib::modeLabelAdd) {

		if (playByTheRules == true) { //check rules, and change previous content of the frame if it doesn't fit to the new type(s)
			for (int i = 0; i < nrOfTypes; i++) { //care for each type in types separately so that trespassings within types are detected and corrected; biggestType: 2^nrOfTypes
				if (sclib::bitTest(types, t) == true) {
					if (areMutualExclusiveTypes(t, this->frameList[FLI][idxAT], problems) == true) {
						this->frameList[FLI][idxAT] ^= problems; //remove problematic types
					}
					this->frameList[FLI][idxAT] |= t;
				}
				t = t << 1;
			}
		} else {
			//care for the types themselfs
			this->frameList[FLI][idxAT] |= types;
		}
		
		//care for the speaker-id; here, no special rules except that it is only set for SPEECH segments apply
    if ((speakerID >= 0 || forceNoSpeaker == true) && (testFrame(FLI, sclib::atSpeech, false, sclib::noType, true, origin) == true)) {
      this->frameList[FLI][idxID] = speakerID;
    }

  } else if (action == sclib::modeLabelRemove) { //TODO

		if (uniteTypes == false) {
			for (int i = 0; i < nrOfTypes; i++) { //remove each type separately
				if (sclib::bitTest(types, t) == true) {
					if (testFrame(FLI, t, false, sclib::noType, true, origin) == true) {
						this->frameList[FLI][idxAT] ^= t;
					}
				}
				t = t << 1;
			}
		} else {
			if (testFrame(FLI, types, true, sclib::noType, true, origin) == true) {
				this->frameList[FLI][idxAT] ^= types;
			}
		}

		if ((speakerID >= 0 || forceNoSpeaker == true) && (testFrame(FLI, sclib::atSpeech, false, sclib::noType, true, origin) == true)) {
      this->frameList[FLI][idxID] = speakerID;
    }

  }

	return;
}

//====================================================================================================================
//  returns the next occasion of the given types/no-types combination (regarding the unite* flags as in testFrame())
//  in the given direction; helpful for e.g. sclib::atSpeechSegment*
//====================================================================================================================
long int SC_GroundTruth::getNextFrame(unsigned long int initFLI, long int types, bool uniteTypes, unsigned int typesNot, bool uniteTypesNot, unsigned int direction, int origin) {
	long int occasion = sclib::noSegment;

  switch (direction) {
    case sclib::searchForward:
		case sclib::searchMiddle:
		case sclib::searchWithin:
	    for (unsigned long int x = initFLI; x < this->internalFrameCount; x++) {
				if (testFrame(x, types, uniteTypes, typesNot, uniteTypesNot, origin) == true) {
					occasion = x;
					break;
		    }
	    }
      break;

    case sclib::searchBackward:
	    for (unsigned long int x = initFLI; x >= 0; x--) {
				if (testFrame(x, types, uniteTypes, typesNot, uniteTypesNot, origin) == true) {
					occasion = x;
					break;
				}
	    }
			break;

    default: break;
  } //switch
	
  return occasion;	
}

//====================================================================================================================
//	set the probability for the occurence of the given type(s) in the framelist from segmentStart to segmentEnd; 
//  if type==sclib::noType, the probability for the speakerID-decision is set
//====================================================================================================================
void SC_GroundTruth::setProbability(unsigned long int segmentStart, unsigned long int segmentEnd, long int types, double probability) {
  unsigned long int segStart = sample2FLI(segmentStart); //convert samples to FLIs, use only FLIs below
  unsigned long int segEnd = sample2FLI(segmentEnd);
	int idx, nrOfTypes = (sizeof(this->frameList[0][0])*8) - 1;
	long int t = 1;

	if (this->pTweak->groundTruth.storeProbabilityInformation == true) {
		if (segStart >= this->internalFrameCount)	{
			segStart = this->internalFrameCount - 1;
		}
		if (segEnd >= this->internalFrameCount)		{
			segEnd = this->internalFrameCount - 1;
		}

		if (sclib::isPowerOfTwo(types) == true) {

			idx = (types == sclib::noType) ? getProbabilityListDim()-1 : sclib::bitPosition(types);
			for (unsigned long int y = segStart; y <= segEnd; y++) {
				this->probabilityList[y][idx] = probability;
			}

		} else {

			for (idx = 0; idx < nrOfTypes; idx++) {
				if (sclib::bitTest(types, t) == true) {
					for (unsigned long int y = segStart; y <= segEnd; y++) {
						this->probabilityList[y][idx] = probability;
					}
				}
				t = t << 1;
			}

		}
	}

	return;
}

//====================================================================================================================
//	Method that labels very short silence-segments (within speech) as 'pause' (instead of silence) and longer 
//	sclib::atSilence segments as sclib::atSilence alone (without the tag sclib::atSpeech)
//	thereby, 'very short' is defined by the parameter 'threshold', which is expressed in samples
//====================================================================================================================
void SC_GroundTruth::silence2pause(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned int threshold, int origin) {
	long int silenceStart, silenceEnd;
	unsigned long int silenceLength, frame;

  //new code that acknowledges the fact that now SPEECH and SILENCE are mutually exclusive
  for (unsigned long int x = segmentStart; x <= segmentEnd; x++) {
		getNextSegment(x, silenceStart, silenceEnd, sclib::atSilence, sclib::searchForward, origin);
		if (silenceStart != sclib::noSegment && silenceEnd != sclib::noSegment && silenceStart < (long int)(segmentEnd)) {
			silenceLength = silenceEnd - silenceStart;
			if (silenceLength < threshold && silenceStart > 0 && silenceEnd < (long int)(this->audioSampleCount-1)) { //very short sclib::atSilence segment not at the borders of the file
        if (testFrame(sample2FLI(silenceStart-1), sclib::atSpeech, false, sclib::noType, false, origin) == true &&
            testFrame(sample2FLI(silenceEnd+1), sclib::atSpeech, false, sclib::noType, false, origin) == true) {
				  for (frame = sample2FLI((unsigned long int)silenceStart); frame <= sample2FLI((unsigned long int)silenceEnd); frame++) {
					  setFrame(frame, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelRemove, origin); //delete sclib::atSilence-mark
            setFrame(frame, sclib::atPause|sclib::atSpeech|sclib::atNoisySpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, origin); //add	sclib::atPause-mark
				  }
        }
			}
			x = silenceEnd; //normally x=silenceEnd+1, but the loop increments x for us
		} else {
			break;
		}
	}

	return;
}

//====================================================================================================================
//	This function marks speech segments (given by their boundarys in samples), which are shorter than 
//  'segmentLengthThreshold' samples, with the sclib::atShort tag
//====================================================================================================================
void SC_GroundTruth::markShortSpeechSegments(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int segmentLengthThreshold) {
  long int speechStart, speechEnd;

  for (unsigned long int y = segmentStart; y <= segmentEnd; y++) {
		this->getNextSegment(y, speechStart, speechEnd, sclib::atSpeech);
		
    if ((speechStart != sclib::noSegment) && (speechEnd != sclib::noSegment) && (speechStart <= (long)(segmentEnd))) {
      if (speechEnd > (long)(segmentEnd)) {
        speechEnd = (long)(segmentEnd);
      }
      if ((unsigned long int)(speechEnd - speechStart + 1) < segmentLengthThreshold) {
				this->setSegment(speechStart, speechEnd, sclib::atShort);
			} 
			y = speechEnd; //normally y=segmentEnd+1, but the next loop increments y for us automatically		
    } else {
      break;
    }
  }

  return;
}

//====================================================================================================================
//	Due to the windowed signal-analysis used during feature-extraction, the last audio-frame in a segment given by
//	start and end (in samples) isn't just end-start: As the frame-length increases, more and more frames at the 
//  segment-end doesn't completely fit into this segment, but overhang into the next one (they start in the desired 
//  segment,  but their end is beyond). So this frames meybe aren't computed during feature-extraction. 
//
//	This function gives the frame-count in a given scene while incorporating this knowledge
//====================================================================================================================
unsigned long int	SC_GroundTruth::getAudioFrameCountInSegment(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned int audioFrameSize, unsigned int audioFrameStep) {
	unsigned long int sampleCount = segmentEnd - segmentStart + 1;
	unsigned long int	frameCount = sclib::getRowCount(sampleCount, audioFrameSize, audioFrameStep); //(sampleCount / audioFrameStep) - (audioFrameSize / audioFrameStep) + 1; //should be ok so, i have veryfied it many times...

	return frameCount;
}

//====================================================================================================================
//	This function gives the last frame-nr in a segment specified by start and end (in samples) while incorporating the 
//  above knowledge.
//====================================================================================================================
unsigned long int SC_GroundTruth::getLastAudioFrameNrInSegment(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned int audioFrameSize, unsigned int audioFrameStep) {
	unsigned long int lastFrame = this->pConverter->sample2audioFrame(segmentStart, audioFrameSize, audioFrameStep) + getAudioFrameCountInSegment(segmentStart, segmentEnd, audioFrameSize, audioFrameStep);
	
	return lastFrame;
}

//====================================================================================================================
// Speaker-specific segments in the frameList are modeled only by a flag at the beginning of a new speaker's segment,
// so the proposed end of the last speech-segment is the beginning of the new one - 1. To get the real end of the last 
// segment, this function counts backwards from the proposed end on till it finds the first sclib::atSpeech-frame. 
// This is then the real segment-end; input- and output-values are sample-based
//====================================================================================================================
unsigned long int	SC_GroundTruth::getSpeechEnd(unsigned long int proposedSpeechEnd, int origin) {
	unsigned long int speechEnd = sample2FLI(proposedSpeechEnd);

  while ((testFrame(speechEnd, sclib::atSpeech, false, sclib::noType, true, origin) == false) && (speechEnd > 0)) {
		speechEnd--;
	}
		
	return FLI2sample(speechEnd, sclib::alignmentEnd);
}

//====================================================================================================================
//	method to test whether a given type exists at least once in the frameList between start and end (or, if end is bad 
//  [<start, >frameCount], internalFrameCount) given in terms of samples.
//====================================================================================================================
bool SC_GroundTruth::existsSegmentType(unsigned long int segmentStart, unsigned long int segmentEnd, long int type, int origin) {
  unsigned long int end = (((segmentEnd > segmentStart) && (sample2FLI(segmentEnd) < this->internalFrameCount)) ? segmentEnd : FLI2sample(this->internalFrameCount, sclib::alignmentEnd));

  for (unsigned long int x = sample2FLI(segmentStart); x < sample2FLI(end); x++) {
    if (testFrame(x, type, false, sclib::noType, true, origin) == true) {
      return true;
    }
	}
	
  return false;
}

//====================================================================================================================
// tests all frames in the given sample-borders and returns the number of samples meeting the requirements.
// andFrames means: all frames have to meet the requirements, otherwise 0 is returned; for other parameters, see 
// method testFrame() for their respective meaning.
//====================================================================================================================
unsigned long int SC_GroundTruth::testSegment(unsigned long int segmentStart, unsigned long int segmentEnd, bool andFrames, unsigned int types, bool uniteTypes, unsigned int typesNot, bool uniteTypesNot, int origin) {
  unsigned long int res = 0, cnt = 0, FLIstart = sample2FLI(segmentStart), FLIend = sclib::min(sample2FLI(segmentEnd), this->internalFrameCount-1);

  for (unsigned long int frame = FLIstart; frame <= FLIend; frame++) {
    if (testFrame(frame, types, uniteTypes, typesNot, uniteTypesNot, origin) == true) {
      res++;
    } else {
      if (andFrames == true) {
        res = 0;
        break;
      }
    }
    cnt++;
  }

  return sclib::min(res*this->internalFrameSize, segmentEnd-segmentStart+1); //(cnt > 0) ? res / cnt : 0;
}

//====================================================================================================================
// can be used instead of testSegment() when - after each internal frame, which's size and concept is hidden outside 
// this class - the result if the frame is positive or negative according to the given criteria is needed; this result
// is returned in terms of nr. of positive samples in the positiveSamples parameter (>0 => FLI was positive), and 
// positiveSamplePosition holds the startSample of this FLI; as long as the function-return-value is true, the tested 
// FLI was between segmentStart and segmentEnd
//====================================================================================================================
bool SC_GroundTruth::testSegmentCallback(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned short int &positiveSamples, unsigned long int &positiveSamplePosition, unsigned int types, bool uniteTypes, unsigned int typesNot, bool uniteTypesNot, int origin) {
  static unsigned long int lastStart = 0;
	static unsigned long int lastEnd = 0;
	static unsigned long int frame = 0;
	bool res;

	positiveSamplePosition = 0;

	if (lastStart == segmentStart && lastEnd == segmentEnd) {
		if (frame <= sclib::min(sample2FLI(segmentEnd), this->internalFrameCount-1)) {
			if (testFrame(frame, types, uniteTypes, typesNot, uniteTypesNot, origin) == true) {
				positiveSamplePosition = this->pConverter->audioFrame2sample(frame, this->internalFrameSize, this->internalFrameSize, sclib::alignmentStart);
				positiveSamples = (unsigned short int)(sclib::min(segmentEnd, this->pConverter->audioFrame2sample(frame, this->internalFrameSize, this->internalFrameSize, sclib::alignmentEnd)) -  positiveSamplePosition + 1);
			} else {
				positiveSamples = 0;
			}
			frame++;
			res = true;
		} else {
			positiveSamples = 0;
			res = false;
		}
	} else {
		lastStart = segmentStart;
		lastEnd = segmentEnd;
		frame = sample2FLI(segmentStart);
		if (testFrame(frame, types, uniteTypes, typesNot, uniteTypesNot, origin) == true) {
			positiveSamplePosition = sclib::max(segmentStart, this->pConverter->audioFrame2sample(frame, this->internalFrameSize, this->internalFrameSize, sclib::alignmentStart));
			positiveSamples = (unsigned short int)(sclib::min(segmentEnd, this->pConverter->audioFrame2sample(frame, this->internalFrameSize, this->internalFrameSize, sclib::alignmentEnd)) -  positiveSamplePosition + 1);
		} else {
			positiveSamples = 0;
		}
		frame++;
		res = true;
	}

  return res;
}

//====================================================================================================================
// Returns true if the questionable frame (internal frame, this method can't be called from outside!) has labels 
//  according to this rules:
//  - types:
//    - uniteTypes==true : every label in types (concatenated with OR: e.g. sclib::atSpeech|sclib::atPause) must apply
//                         to the frame
//    - uniteTypes==false: min. one of the labels in types must apply to the frame
//  - typesNot:
//    - uniteTypesNot==t.: all of the labels in typesNot together must not apply to the frame (a single one is ok, though)
//    - uniteTypesNot==f.: no single one of the labels in typesNot must apply to the frame
//====================================================================================================================
bool SC_GroundTruth::testFrame(unsigned long int FLI, unsigned int types, bool uniteTypes, unsigned int typesNot, bool uniteTypesNot, int origin) {
  unsigned short idx = (origin == sclib::modeHypothesized) ? 1 : 0;
  bool result = false;
  
  if (origin != sclib::modeHypothesized && origin != sclib::modeGroundtruth) {
    REPORT_ERROR(SVLIB_BadArg, "Decide whether to use groundtruth- or hypothesized data!");
  }

  //type is present in this frame or irrelevant
  if ((types == sclib::noType) || ((uniteTypes == false) && (
      ((types & sclib::atPureSpeech)					  && (this->frameList[FLI][idx] & sclib::atPureSpeech)) || 
      ((types & sclib::atNoisySpeech)						&& (this->frameList[FLI][idx] & sclib::atNoisySpeech)) || 
      ((types & sclib::atBackground)				    && (this->frameList[FLI][idx] & sclib::atBackground)) || 
      ((types & sclib::atMusic)				          && (this->frameList[FLI][idx] & sclib::atMusic)) || 
      ((types & sclib::atAction)				        && (this->frameList[FLI][idx] & sclib::atAction)) || 
      ((types & sclib::atBreath)								&& (this->frameList[FLI][idx] & sclib::atBreath)) || 
      ((types & sclib::atUndefined)							&& (this->frameList[FLI][idx] & sclib::atUndefined)) || 
			((types & sclib::atNoise)								  && (this->frameList[FLI][idx] & sclib::atNoise)) || 
      ((types & sclib::atSpeech)							  && (this->frameList[FLI][idx] & sclib::atSpeech)) ||
      ((types & sclib::atPause)								  && (this->frameList[FLI][idx] & sclib::atPause)) ||
      ((types & sclib::atVoiced)  						  && (this->frameList[FLI][idx] & sclib::atVoiced)) ||
      ((types & sclib::atUnvoiced)						  && (this->frameList[FLI][idx] & sclib::atUnvoiced)) ||
      ((types & sclib::atMaleVoice)						  && (this->frameList[FLI][idx] & sclib::atMaleVoice)) ||
      ((types & sclib::atFemaleVoice)					  && (this->frameList[FLI][idx] & sclib::atFemaleVoice)) ||
      ((types & sclib::atSilence)							  && (this->frameList[FLI][idx] & sclib::atSilence)) ||
      ((types & sclib::atShort)           		  && (this->frameList[FLI][idx] & sclib::atShort)) ||
      ((types & sclib::atSceneBoundary)					&& (this->frameList[FLI][idx] & sclib::atSceneBoundary)) ||
      ((types & sclib::atSpeakerBoundary)				&& (this->frameList[FLI][idx] & sclib::atSpeakerBoundary)) ||
      ((types & sclib::atShotBoundary)		      && (this->frameList[FLI][idx] & sclib::atShotBoundary)) ||
      ((types & sclib::atNoiseBoundary)         && (this->frameList[FLI][idx] & sclib::atNoiseBoundary)) ||
      ((types & sclib::atArtificialBoundary)    && (this->frameList[FLI][idx] & sclib::atArtificialBoundary)) ||
			((types & sclib::atSpeechSegmentStart)		&& (this->frameList[FLI][idx] & sclib::atSpeechSegmentStart)) ||
			((types & sclib::atSpeechSegmentEnd)			&& (this->frameList[FLI][idx] & sclib::atSpeechSegmentEnd))
     )) || ((uniteTypes == true) && (
      (!(types & sclib::atPureSpeech)					  || ((types & sclib::atPureSpeech)				    && (this->frameList[FLI][idx] & sclib::atPureSpeech))) && 
      (!(types & sclib::atNoisySpeech)					|| ((types & sclib::atNoisySpeech)				  && (this->frameList[FLI][idx] & sclib::atNoisySpeech))) && 
      (!(types & sclib::atBackground)		      	|| ((types & sclib::atBackground)		        && (this->frameList[FLI][idx] & sclib::atBackground))) && 
      (!(types & sclib::atMusic)			          || ((types & sclib::atMusic)		            && (this->frameList[FLI][idx] & sclib::atMusic))) && 
      (!(types & sclib::atAction)				        || ((types & sclib::atAction)			          && (this->frameList[FLI][idx] & sclib::atAction))) && 
      (!(types & sclib::atBreath)								|| ((types & sclib::atBreath)							  && (this->frameList[FLI][idx] & sclib::atBreath))) && 
      (!(types & sclib::atUndefined)						|| ((types & sclib::atUndefined)			 		  && (this->frameList[FLI][idx] & sclib::atUndefined))) && 
			(!(types & sclib::atNoise)                || ((types & sclib::atNoise)							  && (this->frameList[FLI][idx] & sclib::atNoise))) && 
      (!(types & sclib::atSpeech)               || ((types & sclib::atSpeech)							  && (this->frameList[FLI][idx] & sclib::atSpeech))) &&
      (!(types & sclib::atPause)                || ((types & sclib::atPause)							  && (this->frameList[FLI][idx] & sclib::atPause))) &&
      (!(types & sclib::atVoiced)               || ((types & sclib::atVoiced)  						  && (this->frameList[FLI][idx] & sclib::atVoiced))) &&
      (!(types & sclib::atUnvoiced)             || ((types & sclib::atUnvoiced)						  && (this->frameList[FLI][idx] & sclib::atUnvoiced))) &&
      (!(types & sclib::atMaleVoice)            || ((types & sclib::atMaleVoice)					  && (this->frameList[FLI][idx] & sclib::atMaleVoice))) &&
      (!(types & sclib::atFemaleVoice)          || ((types & sclib::atFemaleVoice)				  && (this->frameList[FLI][idx] & sclib::atFemaleVoice))) &&
      (!(types & sclib::atSilence)              || ((types & sclib::atSilence)						  && (this->frameList[FLI][idx] & sclib::atSilence))) &&
      (!(types & sclib::atShort)                || ((types & sclib::atShort)           		  && (this->frameList[FLI][idx] & sclib::atShort))) &&
      (!(types & sclib::atSceneBoundary)        || ((types & sclib::atSceneBoundary)			  && (this->frameList[FLI][idx] & sclib::atSceneBoundary))) &&
      (!(types & sclib::atSpeakerBoundary)      || ((types & sclib::atSpeakerBoundary)		  && (this->frameList[FLI][idx] & sclib::atSpeakerBoundary))) &&
      (!(types & sclib::atShotBoundary)         || ((types & sclib::atShotBoundary)					&& (this->frameList[FLI][idx] & sclib::atShotBoundary))) &&
      (!(types & sclib::atNoiseBoundary)        || ((types & sclib::atNoiseBoundary)        && (this->frameList[FLI][idx] & sclib::atNoiseBoundary))) &&
      (!(types & sclib::atArtificialBoundary)   || ((types & sclib::atArtificialBoundary)   && (this->frameList[FLI][idx] & sclib::atArtificialBoundary))) &&
			(!(types & sclib::atSpeechSegmentStart)		|| ((types & sclib::atSpeechSegmentStart)		&& (this->frameList[FLI][idx] & sclib::atSpeechSegmentStart))) &&
			(!(types & sclib::atSpeechSegmentEnd)			|| ((types & sclib::atSpeechSegmentEnd)			&& (this->frameList[FLI][idx] & sclib::atSpeechSegmentEnd)))
     ))) {

    if (typesNot == sclib::noType) { //there are no negative types given, so the result is true after passing the positive type test

      result = true;

    } else {
      if (uniteTypesNot == false) { //switch between the methods to interpret the typesNot-list
        
        //if one of the labels in typesNot also appears in the framelist at position 'index', the result is false
        if (!(((typesNot & sclib::atPureSpeech)					  && (this->frameList[FLI][idx] & sclib::atPureSpeech)) || 
              ((typesNot & sclib::atNoisySpeech)					&& (this->frameList[FLI][idx] & sclib::atNoisySpeech)) || 
              ((typesNot & sclib::atBackground)			      && (this->frameList[FLI][idx] & sclib::atBackground)) || 
              ((typesNot & sclib::atMusic)			          && (this->frameList[FLI][idx] & sclib::atMusic)) || 
              ((typesNot & sclib::atAction)				        && (this->frameList[FLI][idx] & sclib::atAction)) || 
              ((typesNot & sclib::atBreath)								&& (this->frameList[FLI][idx] & sclib::atBreath)) || 
              ((typesNot & sclib::atUndefined)						&& (this->frameList[FLI][idx] & sclib::atUndefined)) || 
              ((typesNot & sclib::atNoise)								&& (this->frameList[FLI][idx] & sclib::atNoise)) || 
              ((typesNot & sclib::atSpeech)							  && (this->frameList[FLI][idx] & sclib::atSpeech)) ||
              ((typesNot & sclib::atPause)								&& (this->frameList[FLI][idx] & sclib::atPause)) ||
              ((typesNot & sclib::atVoiced)  						  && (this->frameList[FLI][idx] & sclib::atVoiced)) ||
              ((typesNot & sclib::atUnvoiced)						  && (this->frameList[FLI][idx] & sclib::atUnvoiced)) ||
              ((typesNot & sclib::atMaleVoice)					  && (this->frameList[FLI][idx] & sclib::atMaleVoice)) ||
              ((typesNot & sclib::atFemaleVoice)				  && (this->frameList[FLI][idx] & sclib::atFemaleVoice)) ||
              ((typesNot & sclib::atSilence)							&& (this->frameList[FLI][idx] & sclib::atSilence)) ||
              ((typesNot & sclib::atShort)           		  && (this->frameList[FLI][idx] & sclib::atShort)) ||
              ((typesNot & sclib::atSceneBoundary)			  && (this->frameList[FLI][idx] & sclib::atSceneBoundary)) ||
              ((typesNot & sclib::atSpeakerBoundary)		  && (this->frameList[FLI][idx] & sclib::atSpeakerBoundary)) ||
              ((typesNot & sclib::atShotBoundary)					&& (this->frameList[FLI][idx] & sclib::atShotBoundary)) ||
              ((typesNot & sclib::atNoiseBoundary)        && (this->frameList[FLI][idx] & sclib::atNoiseBoundary)) ||
              ((typesNot & sclib::atArtificialBoundary)   && (this->frameList[FLI][idx] & sclib::atArtificialBoundary)) ||
			        ((typesNot & sclib::atSpeechSegmentStart)		&& (this->frameList[FLI][idx] & sclib::atSpeechSegmentStart)) ||
			        ((typesNot & sclib::atSpeechSegmentEnd)			&& (this->frameList[FLI][idx] & sclib::atSpeechSegmentEnd)))
           ) {
          result = true;
         }

      } else {

        //all of the labels in typesNot together must not apply to the frame (a single one is ok, though)
        //none of the labels in typesNot appears in the frameList at position 'index'
        if (
            (!(typesNot & sclib::atPureSpeech)					  || ((typesNot & sclib::atPureSpeech)						&& !(this->frameList[FLI][idx] & sclib::atPureSpeech))) &&
            (!(typesNot & sclib::atNoisySpeech)					  || ((typesNot & sclib::atNoisySpeech)						&& !(this->frameList[FLI][idx] & sclib::atNoisySpeech))) &&
            (!(typesNot & sclib::atBackground)			      || ((typesNot & sclib::atBackground)			     	&& !(this->frameList[FLI][idx] & sclib::atBackground))) &&
            (!(typesNot & sclib::atMusic)			            || ((typesNot & sclib::atMusic)				          && !(this->frameList[FLI][idx] & sclib::atMusic))) &&
            (!(typesNot & sclib::atAction)			          || ((typesNot & sclib::atAction)				        && !(this->frameList[FLI][idx] & sclib::atAction))) &&
            (!(typesNot & sclib::atBreath)								|| ((typesNot & sclib::atBreath)								&& !(this->frameList[FLI][idx] & sclib::atBreath))) &&
            (!(typesNot & sclib::atUndefined)						  || ((typesNot & sclib::atUndefined)							&& !(this->frameList[FLI][idx] & sclib::atUndefined))) &&
            (!(typesNot & sclib::atNoise)					  		  || ((typesNot & sclib::atNoise)								  && !(this->frameList[FLI][idx] & sclib::atNoise))) &&
            (!(typesNot & sclib::atSpeech)								|| ((typesNot & sclib::atSpeech)								&& !(this->frameList[FLI][idx] & sclib::atSpeech))) &&
            (!(typesNot & sclib::atPause)								  || ((typesNot & sclib::atPause)								  && !(this->frameList[FLI][idx] & sclib::atPause))) &&
            (!(typesNot & sclib::atVoiced)  							|| ((typesNot & sclib::atVoiced)  							&& !(this->frameList[FLI][idx] & sclib::atVoiced))) &&
            (!(typesNot & sclib::atUnvoiced)							|| ((typesNot & sclib::atUnvoiced)							&& !(this->frameList[FLI][idx] & sclib::atUnvoiced))) &&
            (!(typesNot & sclib::atMaleVoice)							|| ((typesNot & sclib::atMaleVoice)							&& !(this->frameList[FLI][idx] & sclib::atMaleVoice))) &&
            (!(typesNot & sclib::atFemaleVoice)						|| ((typesNot & sclib::atFemaleVoice)						&& !(this->frameList[FLI][idx] & sclib::atFemaleVoice))) &&
            (!(typesNot & sclib::atSilence)							  || ((typesNot & sclib::atSilence)							  && !(this->frameList[FLI][idx] & sclib::atSilence))) &&
            (!(typesNot & sclib::atShort)								  || ((typesNot & sclib::atShort)								  && !(this->frameList[FLI][idx] & sclib::atShort))) &&
			      (!(typesNot & sclib::atSceneBoundary)					|| ((typesNot & sclib::atSceneBoundary)					&& !(this->frameList[FLI][idx] & sclib::atSceneBoundary))) &&
            (!(typesNot & sclib::atSpeakerBoundary)				|| ((typesNot & sclib::atSpeakerBoundary)				&& !(this->frameList[FLI][idx] & sclib::atSpeakerBoundary))) &&
            (!(typesNot & sclib::atShotBoundary)			    || ((typesNot & sclib::atShotBoundary)			    && !(this->frameList[FLI][idx] & sclib::atShotBoundary))) &&
            (!(typesNot & sclib::atNoiseBoundary)					|| ((typesNot & sclib::atNoiseBoundary)         && !(this->frameList[FLI][idx] & sclib::atNoiseBoundary))) &&
            (!(typesNot & sclib::atArtificialBoundary)    || ((typesNot & sclib::atArtificialBoundary)    && !(this->frameList[FLI][idx] & sclib::atArtificialBoundary))) &&
            (!(typesNot & sclib::atSpeechSegmentStart)		|| ((typesNot & sclib::atSpeechSegmentStart)		&& !(this->frameList[FLI][idx] & sclib::atSpeechSegmentStart))) &&
			      (!(typesNot & sclib::atSpeechSegmentEnd)			|| ((typesNot & sclib::atSpeechSegmentEnd)			&& !(this->frameList[FLI][idx] & sclib::atSpeechSegmentEnd)))
           ) {  
          result = true;
        }

      } //andTypesNot = true
    } //typesNot != sclib::noType
  } //types found

  return result;
}

//====================================================================================================================
// Return the speaker-id for the given internal frame
//====================================================================================================================
long int SC_GroundTruth::getSpeakerID(unsigned long int FLI, int origin) {
  unsigned short idx = (origin == sclib::modeHypothesized) ? 3 : 2;

  if (origin != sclib::modeHypothesized && origin != sclib::modeGroundtruth) {
    REPORT_ERROR(SVLIB_BadArg, "Decide whether to use groundtruth- or hypothesized data!");
  }
  if (FLI >= this->internalFrameCount) {
    FLI = this->internalFrameCount - 1;
  }

  return this->frameList[FLI][idx];
}

//====================================================================================================================
// Return the speaker-id for the given sample-nr
//====================================================================================================================
long int SC_GroundTruth::getSamplesSpeakerID(unsigned long int sample, int origin) {
  return getSpeakerID(sample2FLI(sample), origin);
}

//====================================================================================================================
// Test if the given segment (given by boundarys in samples) is speaker-homogenious (all frames are labeled with the 
// same speaker-id)
//====================================================================================================================
bool SC_GroundTruth::isSpeakerHomogenious(long int segmentStart, long int segmentEnd, int origin) {
  unsigned long int startFrame = sample2FLI((unsigned long)segmentStart), endFrame = sample2FLI((unsigned long)segmentEnd);
	long int oldSpeakerID = sclib::noSpeaker, currentSpeakerID = sclib::noSpeaker;
	bool result = true;
    
  for (unsigned long int frame = startFrame; frame <= endFrame; frame++) {
    if (testFrame(frame, sclib::atSpeech, false, sclib::noType, true, origin) == true) { //there may be parts during a segment without speech and therefore also with no valid speaker label (see function below): don't evaluate them for this purpose
			currentSpeakerID = getSpeakerID(frame, origin);
      if (oldSpeakerID != sclib::noSpeaker && oldSpeakerID != currentSpeakerID) {
        result = false;
        break;
      }
			oldSpeakerID = currentSpeakerID;
    }
  }

  return result;
}

//====================================================================================================================
// Test if the given segment (given by boundarys in samples) is speaker-homogenious (all frames are labeled with the 
// same speaker-id) and if it's speaker is the one specified by speaker-name; does only make sense on ground-truth 
// data, because only there speaker names may be available
//====================================================================================================================
bool SC_GroundTruth::isSpeakerHomogenious(long int segmentStart, long int segmentEnd, const char* speakerName) {
  bool result = true;
  unsigned long start = sample2FLI((unsigned long int)segmentStart), end = sample2FLI((unsigned long int)segmentEnd);
  long int oldSpeakerID, speakerID = getSpeakerID(start, sclib::modeGroundtruth), wantedSpeakerID = getSpeakerIDfromName(speakerName, sclib::bufferSize);
  
  if (wantedSpeakerID == sclib::noSpeaker || wantedSpeakerID != speakerID) {
    result = false;
  } else {
    for (unsigned long int frame = start+1; frame <= end; frame++) {
      oldSpeakerID = speakerID;
      speakerID = getSpeakerID(frame, sclib::modeGroundtruth);
      if (speakerID != oldSpeakerID && speakerID != sclib::noSpeaker && oldSpeakerID != sclib::noSpeaker) { //there may be parts during a segment without speech and therefore also with no valid speaker label
        result = false;
        break;
      }
    }
  }

  return result;
}

//====================================================================================================================
// gives for a given sampleNr the scene in which it lies
// if lastSample and lastScene are given, the search starts from lastSample, assuming scene-nr lastScene at this 
// position
//====================================================================================================================
unsigned long int SC_GroundTruth::sample2scene(unsigned long int sample, unsigned long int lastSample, unsigned long int lastScene){
	unsigned long int y = 0, scene = 1;
	long int segStart, segEnd;

	//use cached previous results, if available
	if (lastSample > 0 && lastScene > 0) {
		y = lastSample;
		scene = lastScene;
	}

	while (y < sample) {
		getNextBoundary(y, segStart, segEnd, sclib::atSceneBoundary, sclib::searchMiddle);
		if (segStart == sclib::noSegment || segEnd == sclib::noSegment || sclib::intersect(segStart, segEnd, sample, sample) > 0) {
			break;
		} else {
			scene++;
			y = segEnd + 1;
		}
	}

	return scene;
}

//====================================================================================================================
// gives for a given sampleNr the shot in which it lies
// if lastSample and lastShot are given, the search starts from lastSample, assuming shot-nr lastShot at this position
//====================================================================================================================
unsigned long int SC_GroundTruth::sample2shot(unsigned long int sample, unsigned long int lastSample, unsigned long int lastShot){
	unsigned long int y = 0, shot = 1;
	long int segStart, segEnd;

	//use cached previous results, if available
	if (lastSample > 0 && lastShot > 0) {
		y = lastSample;
		shot = lastShot;
	}

	while (y < sample) {
		getNextBoundary(y, segStart, segEnd, sclib::atShotBoundary, sclib::searchMiddle);
		if (segStart == sclib::noSegment || segEnd == sclib::noSegment || sclib::intersect(segStart, segEnd, sample, sample) > 0) {
			break;
		} else {
			shot++;
			y = segEnd + 1;
		}
	}

	return shot;
}

//====================================================================================================================
// convert sample-nr to internal audioframe-nr
//====================================================================================================================
unsigned long int SC_GroundTruth::sample2FLI(unsigned long int sample) {
  return this->pConverter->sample2audioFrame(sample, this->internalFrameSize, this->internalFrameSize, sclib::alignmentStart);
}

//====================================================================================================================
// convert internal audioframe-nr to sample-nr
//====================================================================================================================
unsigned long int SC_GroundTruth::FLI2sample(unsigned long int frameListIndex, unsigned int alignment) {
	unsigned long int sample = this->pConverter->audioFrame2sample(frameListIndex, this->internalFrameSize, this->internalFrameSize, alignment);

	return sclib::min(sample, this->audioSampleCount-1); //necessary because the last internal frames can exceed the number available samples
}

//====================================================================================================================
//	print the frameList from the given start- to end-sample into a file
//====================================================================================================================
void SC_GroundTruth::frameListOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types, unsigned long int typesNot, sclib::OpenMode mode, int origin, bool suppressHeader) {
	fstream fileOut;
	char *fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];

	sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, mode);
	MFree_1D(fName);
	output(fileOut, segmentStart, segmentEnd, types, typesNot, origin, suppressHeader);
	fileOut.close();

	return;
}

//====================================================================================================================
//	output all analyzed parts of the frameList (acording to selected scenes in the tweakable parameters)
//====================================================================================================================
void SC_GroundTruth::analyzedOut(const char* fileName, sclib::OpenMode mode) {
  long int sceneStart, sceneEnd;
  unsigned long int sceneNr = 1;
	fstream fileOut;
	char *fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
  bool printedSomething = false;

	sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, mode);
	MFree_1D(fName);

  for (unsigned long int y = 0; y < this->audioSampleCount; y++) {
    getNextBoundary(y, sceneStart, sceneEnd, sclib::atSceneBoundary);
    if (sceneStart != sclib::noSegment && sceneEnd != sclib::noSegment) {
      if ((sceneNr < pTweak->general.firstScene) || (!(sclib::bit(sceneNr) & pTweak->general.sceneSelection) && (pTweak->general.sceneSelection != 0))) {
        sceneNr++; 
        y = sceneEnd; 
        continue;
      }

      output(fileOut, sceneStart, sceneEnd, sclib::noType, sclib::noType, sclib::modeHypothesized, printedSomething);
      printedSomething = true;

      sceneNr++;
      y = sceneEnd;
      if (sceneNr > pTweak->general.lastScene) {
				break;
			}
    }
  }

  fileOut.close();

	return;
}

//====================================================================================================================
//  To make operator<<() kind of virtual...
//====================================================================================================================
ostream& SC_GroundTruth::output(ostream& OutS, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types, unsigned long int typesNot, int origin, bool suppressHeader) {
  long int id, width = 3;
  unsigned long int x, start = sample2FLI(segmentStart), end = sample2FLI(segmentEnd);
  bool correct;
  char hypothesizedSpeakerName[sclib::bufferSize];

  OutS << setiosflags(ios_base::left|ios_base::fixed|ios_base::showpoint);

  if (suppressHeader == false) {
		OutS << setw(10)				<< "FRAME#"
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
				 << setw(2*15+2)		<< "SPK"
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
//  operator<<(), so that output can be directed to different streams
//====================================================================================================================
ostream& operator<< (ostream& OutS, SC_GroundTruth& pGT) {
  return pGT.output(OutS, 0, pGT.getAudioSampleCount());
}

//====================================================================================================================
// for a given segment, this function writes the ratio of brutto speech length (all frames labeled sclib::atSpeech) to
// netto speech length (all sclib::atSpeech-frames without sclib::atPause|sclib::atSilence|sclib::atUnvoiced) into a 
// file; also tells how many seconds of speech where too short to analyze
//====================================================================================================================
void SC_GroundTruth::segmentStatisticsOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd) {
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
      sprintf(outLine, "start: %12u;\t end: %12u;\t brutto: %8.3f s; \tnetto: %8.3f s; \tratio: %7.3f %%; \ttoo short: %8.3f s\0", this->pConverter->sample2ms(segmentStart), this->pConverter->sample2ms(segmentEnd), (double)this->pConverter->audioFrame2ms(bruttoLength, this->internalFrameSize, this->internalFrameSize, sclib::alignmentEnd)/1000.0, (double)this->pConverter->audioFrame2ms(nettoLength, this->internalFrameSize, this->internalFrameSize, sclib::alignmentEnd)/1000.0, ((double)nettoLength/(double)bruttoLength)*100, (double)this->pConverter->audioFrame2ms(shortLength, this->internalFrameSize, this->internalFrameSize, sclib::alignmentEnd)/1000.0);
	    sclib::scalarOut(fileName, outLine, this->pTweak);
	      
      y = speechEnd;
		} //segment-boundarys are valid
		else {break;}
	} //for all frames in this scene

	MFree_1D(outLine);
	
	return;
}

//====================================================================================================================
// The speech segment length in this segment (given by borders in samples) will be set to the given parameters: 
// Longer segments will be split, shorter segments will be erased. All values are sample-based 
//====================================================================================================================
void SC_GroundTruth::setSpeechSegLength(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int minLength, unsigned long int maxLength) {
	setMaxSpeechSegLength(segmentStart, segmentEnd, maxLength);
	setMinSpeechSegLength(segmentStart, segmentEnd, minLength);

	return;
}

//====================================================================================================================
// The maximum speech segment length in the specified scene will be set to maxLength samples; longer segments
// will be split, with speaker-boundarys before each newly created segment!
//====================================================================================================================
void SC_GroundTruth::setMaxSpeechSegLength(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int maxLength, int origin) {
	unsigned long int x, speechLength, frame;
	long int speechStart, speechEnd;

	for (x = segmentStart; x <= segmentEnd; x++) {
		getNextSegment(x, speechStart, speechEnd, sclib::atSpeech, sclib::searchForward, origin);
		if (speechStart != sclib::noSegment && speechEnd != sclib::noSegment) {
			speechLength = speechEnd - speechStart;
			if (speechLength > maxLength) {
        frame = sample2FLI(speechStart + maxLength);
        setFrame(frame+1, sclib::atSpeech, false, sclib::noSpeaker, sclib::modeLabelRemove, origin);
				if (testFrame(frame+2, sclib::atSpeech, false, sclib::noType, true, origin) == true) {
          setFrame(frame+2, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, origin);
        }
				x = FLI2sample(frame + 1, sclib::alignmentEnd);
			} else {
				x = speechEnd;
			}
		} else {
			break;
		}
	}

	return;
}

//====================================================================================================================
// The minimum speech segment length in the specified scene will be set to minLength samples; shorter segments
// will be erased
//====================================================================================================================
void SC_GroundTruth::setMinSpeechSegLength(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int minLength) {
	unsigned long int x, speechLength;
	long int speechStart, speechEnd;

	for (x = segmentStart; x <= segmentEnd; x++) {
		getNextSegment(x, speechStart, speechEnd, sclib::atSpeech);
		if (speechStart != sclib::noSegment && speechEnd != sclib::noSegment) {
			speechLength = speechEnd - speechStart;
			if (speechLength < minLength) {
				setSegment(speechStart, speechEnd, sclib::atSpeech, false, sclib::noSpeaker, sclib::modeLabelRemove);
				setSegment(speechStart, speechEnd, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelRemove);
				setSegment(speechStart, speechEnd, sclib::atShort);
			}
			x = speechEnd;
		} else {
			break;
		}
	}

	return;
}

//====================================================================================================================
// Copys all frames between start and end together in one new SV_Data-object, which has one of the frame-types
// defined by 'types' (types may be multiple types added together with OR: e.g. sclib::atSpeech || sclib::atNoise || 
// sclib::atPause) and definitely not one of the frame-types defined by 'typesNot'.
// If needed, one can specifiy the parameters additionalCols and startCol: The new SV_Data-object will then have 
// additionalCols more columns than it's parent, with the original data copied from startCol on.
// If setStartSample is true and additionalCols is at least 1 and startCol is at least 1, the number of the first 
// sample belonging to this frame (overall) is inserted in the zeroth col of each vector. If startCol or 
// additionalCols is 0 but setting startSample is wished anyway, a second SV_Data-object is created with only 1 column
// and rows corresponding to the rows of the returned object; this 2nd oject will have the sample-numbers in its 
// single column and will be accessible via the ->Next pointer of the first SV_Data-object.
// 'start' and 'end' must be real sample-numbers, not relative to 'offset'!!!
// The offset means which sample is the first one in the first frame (first row) of the pComplete object
//====================================================================================================================
SV_Data* SC_GroundTruth::copyFramesTogether(SV_Data* pComplete, unsigned long int offset, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types, unsigned long int typesNot, unsigned int additionalCols, unsigned int startCol, bool setStartSample) {
  unsigned long int frame, x, z = 0, dim, segmentLength = 0, start, end;
	double srRatio = pComplete->Hdr.sampleRate / (double)(this->pSignalPrototype->SigPar.SRate);
	unsigned long int frameSize = sclib::round(pComplete->Hdr.frameSize / srRatio), frameStep = sclib::round(pComplete->Hdr.frameStep / srRatio); //frame-parameters in terms of the corpus' sample-rate (as offset, segStart & segEnd are supposed to be)
  SV_Data* pCollection = NULL, *pStartSamples = NULL;
  unsigned long int startFrame = this->pConverter->sample2audioFrame(segmentStart-offset, frameSize, frameStep);
  unsigned long int endFrame = startFrame + getAudioFrameCountInSegment(segmentStart, segmentEnd, frameSize, frameStep);

	if (endFrame > (unsigned long int)pComplete->Row) {
    endFrame = pComplete->Row;
  }

  //get length of actual segment
  for (frame = startFrame; frame < endFrame; frame++) {
    start = offset + this->pConverter->audioFrame2sample(frame, frameSize, frameStep, sclib::alignmentStart);
    end = offset + this->pConverter->audioFrame2sample(frame, frameSize, frameStep, sclib::alignmentEnd);
    if (testSegment(start, end, true, types, false, typesNot, false) > 0) {
      segmentLength++;    
    }
  }

  if (segmentLength > 0) {
    dim = pComplete->Col;
	  pCollection = new SV_Data(segmentLength, dim + additionalCols);
    pCollection->Hdr = pComplete->Hdr;
		//for (z = 0; z < 8; z++) {
		//	pCollection->Hdr.Signature[z] = pComplete->Hdr.Signature[z];
		//}
		if (setStartSample == true && (additionalCols < 1 || startCol < 1)) { //create an additional data-object to store the start-samples of the frames correspondign to the rows in pCollection if no additional space is there to stoe it in the 0th column of pCollection itself but shpuld be done anyway
			pStartSamples = new SV_Data(segmentLength, 1);
			pStartSamples->Hdr = pComplete->Hdr;
			pCollection->Next = pStartSamples;
		}

    //copy relevant frames
		z = 0;
    for (frame = startFrame; frame < endFrame; frame++) {
      if (testSegment(offset + this->pConverter->audioFrame2sample(frame, frameSize, frameStep, sclib::alignmentStart), 
                      offset + this->pConverter->audioFrame2sample(frame, frameSize, frameStep, sclib::alignmentEnd),
                      true, types, false, typesNot, false) > 0) {
			  if (setStartSample == true) {
					if ((additionalCols >= 1) && (startCol >= 1)) {
						pCollection->Mat[z][0] = (float)(offset + this->pConverter->audioFrame2sample(frame, frameSize, frameStep, sclib::alignmentStart));
						for (x = 0; x < dim; x++) {
							pCollection->Mat[z][x+startCol] = pComplete->Mat[frame][x];
						}
					} else {
						pStartSamples->Mat[z][0] = (float)(offset + this->pConverter->audioFrame2sample(frame, frameSize, frameStep, sclib::alignmentStart));
						for (x = 0; x < dim; x++) {
							pCollection->Mat[z][x+startCol] = pComplete->Mat[frame][x];
						}
					}
			  } else {
				  for (x = 0; x < dim; x++) {
					  pCollection->Mat[z][x+startCol]	=	pComplete->Mat[frame][x];
				  }
			  }
        z++;
      } //if testSegment()
    } //for frame
  } //segmentLength > 0

  return pCollection;
}

//====================================================================================================================
//  Read an ASCII-file containing a list when to load which explicit background model. It returns the nr of rows of 
//  the table explicitModels, which has the form desired  by SC_SpeechProcessing::buildSpeakerModels():
//    - explicitModels[x][0][0]: sceneNumber or 0 (then the modell applies to all scenes)
//    - explicitModels[x][1][0]: validSegmentNumber or 0 (then the modell applies for all valid segments in the scene)
//    - explicitModels[x][2]   : the full path/filename to the model-file to load
//  Here's a constrained: the scene-/segment-numbers must be <=127 and the filename has to be <=255 characters in 
//  length.
//====================================================================================================================
unsigned short int SC_GroundTruth::readExplicitModelList(const char* modelListFileName, char*** &explicitModels) {
	FILE* inFile;
  char* buffer	= new char[sclib::bufferSize];
	unsigned short int count = 0, max = 0;
	int value;

  if (modelListFileName == NULL || strcmp(modelListFileName, "") == 0) {
    MFree_1D(buffer);
    return 0;
  }

  inFile = fopen(modelListFileName, "r");
	if (inFile == NULL) {
    MFree_1D(buffer);
    return 0;
  }

	//get count of relevant lines
	while (!feof(inFile)) {
    sclib::readline(inFile, buffer, sclib::bufferSize);
    if (sclib::isNum(buffer[0]) == true) {
			max++;
		}
	}
	fseek(inFile, 0, SEEK_SET);

	MArray_3D(explicitModels, max, 3, sclib::bufferSize, char, "readExplicitModelList: explicitModels"); 
  
	//extract the information form the file
	while (!feof(inFile)) {
    sclib::readline(inFile, buffer, sclib::bufferSize);
    if (sclib::isNum(buffer[0]) == true) {
      value = sclib::getNextIntFromString(buffer, sclib::bufferSize);
      assert((value >= 0) && (value < 128));
			explicitModels[count][0][0] = (char)(value); //scene-nr

      value = sclib::getNextIntFromString(buffer, sclib::bufferSize);
      assert((value >= 0) && (value < 128));
			explicitModels[count][1][0] = (char)(value); //valid-segment-nr

      sclib::lTrim(buffer);
			sprintf(explicitModels[count][2], "%s\0", buffer); //path to model-file
			count++;
		}
	}

	fclose(inFile);
	MFree_1D(buffer);
	return count;
}

//====================================================================================================================
// establishes a mapping between speaker-id's from the groundtruth and from the clustering-process; if the mapping
// already exists, just the 'correct'-flag is altered
//====================================================================================================================
void SC_GroundTruth::addSpeakerMapping(long int groundTruthID, long int hypothesizedID, bool correct) {
  SC_GroundTruth::SC_SpeakerMapping *sm = new SC_GroundTruth::SC_SpeakerMapping(), *pHook;

  sm->groundTruthID = groundTruthID;
  sm->hypothesizedID = hypothesizedID;
  sm->correct = correct;
  sm->Next = NULL;
  
  if (this->pSpeakerMapping == NULL) { //this is the first entry
    this->pSpeakerMapping = sm;
  } else {
    pHook = this->pSpeakerMapping;
    while (pHook->Next != NULL) {
      if (pHook->groundTruthID == groundTruthID && pHook->hypothesizedID == hypothesizedID) { //don't add a new mapping if it already exists; just alter the 'correct'-flag
        pHook->correct = correct;
        MFree_0D(sm);
        break;
      }
      pHook = pHook->Next;
    }
    if (pHook->groundTruthID == groundTruthID && pHook->hypothesizedID == hypothesizedID) { //the last entry is missed by the above loop... don't add a new mapping if it already exists; just alter the 'correct'-flag
      pHook->correct = correct;
      MFree_0D(sm);
    }
    if (sm != NULL) { //add a new mapping if it wasn't already found
      pHook->Next = sm; 
    }
  }
  
  return;
}

//====================================================================================================================
// returns the corresponding groundTruth-speakerID, if a speakermapping was established earlier
//====================================================================================================================
long int SC_GroundTruth::getSpeakerGIDFromHID(long int hypothesizedID, bool &correct) {
  SC_GroundTruth::SC_SpeakerMapping *pHook = this->pSpeakerMapping;
  long ID = sclib::noSpeaker;

  correct = false;

  //find the index
  while (pHook != NULL) {
    if (pHook->hypothesizedID == hypothesizedID) {
      ID = pHook->groundTruthID;
      correct = pHook->correct;
      break;
    }
    pHook = pHook->Next;
  }

  return ID;
}

//====================================================================================================================
// returns the corresponding hypothesized speakerID, if a speakermapping was established earlier
//====================================================================================================================
long int SC_GroundTruth::getSpeakerHIDFromGID(long int groundTruthID, bool &correct) {
  SC_GroundTruth::SC_SpeakerMapping *pHook = this->pSpeakerMapping;
  long ID = sclib::noSpeaker;

  correct = false;

  //find the index
  while (pHook != NULL) {
    if (pHook->groundTruthID == groundTruthID) {
      ID = pHook->hypothesizedID;
      correct = pHook->correct;
      break;
    }
    pHook = pHook->Next;
  }

  return ID;
}

//====================================================================================================================
// removes all speaker-mappings
//====================================================================================================================
void SC_GroundTruth::removeAllSpeakerMappings(void) {
  sclib::destructLinkedList(this->pSpeakerMapping);

  return;
}

//====================================================================================================================
// The most common audio-type in the given segment is returned, as well as (in the percent parameter) the percentage 
// with which it occurs; in the statistics-parameter, an array is returned in which each cell tells how many FLIs of 
// the type of pow(2, index) are present in the segment (the last (32th) entry holds the sum of the previous ones); 
// wantedTypes specifies which audio-types should be considered, default is all. 
// This method exploits that the audio-types are just bitflags
//====================================================================================================================
unsigned long int SC_GroundTruth::getPrevailingAudioType(unsigned long int segmentStart, unsigned long int segmentEnd, double &percent, unsigned long int *&statistics, unsigned long int wantedTypes, int origin) {
	unsigned long int audioType, type, max = 0, bits = sizeof(this->frameList[0][0])*8-1;
	unsigned long int frame, x, startFrame, endFrame;

	startFrame = sample2FLI(segmentStart);
	endFrame = sample2FLI(segmentEnd);

	//init return values
	MFree_1D(statistics);
	MArray_1D(statistics, bits+1, unsigned long int, "SC_GroundTruth.getPrevailingAudioType: statistics");
	for (x = 0; x <= bits; x++) {
		statistics[x] = 0;
	}
  percent = 0.0;
	audioType = sclib::noType;

  //create a statistic: which (wanted) audioType is present how ofton in the current segment, with special care for soft boundaries
  for (frame = startFrame; frame <= endFrame; frame++) {
    type = 1;
    for (x = 0; x < bits; x++) { //check for each audioType if it is wanted and if it is present in the frame (there are 32 possible audioTypes due to a 32 sclib::bit variable at the moment...)
      if (sclib::bitTest(wantedTypes, type) && (testFrame(frame, type, false, sclib::noType, true, origin) == true)) {
        statistics[x] += this->internalFrameSize;
			}
      type *= 2;
    }
  }

  //find the winner
  type = bits;
  for (x = 0; x < bits; x++) {
		statistics[bits] += statistics[x]; //copy the sum of the previous elements into the last one
    if (statistics[x] > max) {
      type = x;
      max = statistics[x];
    }
  }

  //calculate result
  if (type < bits) {
    audioType = (int)(pow(2.0, (double)(type)));
		percent = (statistics[bits] > 0) ? (double)(statistics[type]) / (double)(statistics[bits]) : 0.0; //statistics[bits] is the sum
	}

  return audioType;
}

//====================================================================================================================
//	Returns the ID of the speaker that appears most ofton in the given segment or sclib::noSpeaker if there is none
//====================================================================================================================
long int SC_GroundTruth::getMajorSpeakerID(unsigned long int segmentStart, unsigned long int segmentEnd, int origin) {
	unsigned long int FLI, startFrame, endFrame;
	long int max = 0, winner = sclib::noSpeaker;
	int idxSpk = (origin == sclib::modeHypothesized) ? 3 : 2;
	std::map<long int, int> speakerCount;

	if (segmentStart != sclib::noSegment && segmentEnd != sclib::noSegment) {
		startFrame = sample2FLI(segmentStart);
		endFrame = sample2FLI(segmentEnd);
		
		for (FLI = startFrame; FLI <= endFrame; FLI++) {
			if (this->frameList[FLI][idxSpk] != sclib::noSpeaker) {
				speakerCount[this->frameList[FLI][idxSpk]] += 1;
			}
		}

		for (std::map<long int, int>::iterator i = speakerCount.begin(); i != speakerCount.end(); i++) {
			if (i->second > max) {
				winner = i->first;
				max = i->second;
			}
		}
	}

	return winner;
}

//====================================================================================================================
// Go through the given segment and return:
//  - overallGT: #samples having the type(s) according to groundTruth
//  - overallHypo: #samples having the type(s) according to algorithmic results
//  - overallContradiction: #samples that have the type(s) according to origin but not according to the other column
//  - overallAgreement: #samples havin the type(s) in both GT and hypo columns
//====================================================================================================================
void SC_GroundTruth::getTypeStatistics(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int &overallGT, unsigned long int &overallHypo, unsigned long int &overallErrors, unsigned long int &overallAgreement, unsigned long int types, bool uniteTypes, int origin){
	unsigned long int frame, segStart = sample2FLI(segmentStart), segEnd = sample2FLI(segmentEnd);
	bool hypoHit, gtHit;

	overallGT = 0;
	overallHypo = 0;
	overallAgreement = 0;
	overallErrors = 0;

	for (frame = segStart; frame < segEnd; frame++) {
		gtHit = testFrame(frame, types, uniteTypes, sclib::noType, false, sclib::modeGroundtruth);
    hypoHit = testFrame(frame, types, uniteTypes, sclib::noType, false, sclib::modeHypothesized);
		if (gtHit == true) {
			overallGT++;
		}
		if (hypoHit == true) {
			overallHypo++;
		}
		if (gtHit == true && hypoHit == true) {
			overallAgreement++;
		} else if ((origin == sclib::modeGroundtruth && gtHit == true && hypoHit == false) ||
			         (origin == sclib::modeHypothesized && gtHit == false && hypoHit == true)) {
			overallErrors++;
		}
	}

	//convert FLIs to samples
	overallGT *= this->internalFrameSize;
	overallHypo *= this->internalFrameSize;
	overallAgreement *= this->internalFrameSize;
	overallErrors *= this->internalFrameSize;

	return;
}

//====================================================================================================================
// check the ground-truth (or hypothesized results) for consistency; this includes the check for mutually exclusive
// tupels of audio-types (e.g. no internal frame can be speech and noise together); returns true if everything is ok
// violations of rules are reported in a file or on screen (fileName != "")
//====================================================================================================================
bool SC_GroundTruth::checkConsistency(unsigned long int segmentStart, unsigned long int segmentEnd, int origin) {
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
// Returns the total width of the area in which an event may really have been occured of it is told in the 
// groundtruth to have been occured at exact position X (the middle of that region);
//====================================================================================================================
unsigned long int SC_GroundTruth::getUncertaintyRegionWidth(bool includeGroundtruthUncertainty) {
	long int sum = 0, maxFeatureFrameStep = 0;

	sum = 2 * this->internalFrameSize; //due to quantization effects during result storing, double because of +-

	if (includeGroundtruthUncertainty == true) {
		sum += this->uncertaintyRegion; //uncertainty during (hand made) ground-truthing
	}

	return sum;
}

//====================================================================================================================
// Returns true if the first parameter represents a subtype of the second one (e.g. MUSIC is a subtype of NOISE), 
// false if it isn't (e.g. because the first parameter is never a subtype)
//====================================================================================================================
bool SC_GroundTruth::isSubType(long int questionableType, long int type) {
	switch (questionableType) {
		case sclib::atNoise:
			return false;
		case sclib::atPureSpeech:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atNoisySpeech:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atBackground:
			return (type == sclib::atNoise) ? true : false;
		case sclib::atMusic:
			return (type == sclib::atNoise) ? true : false;
		case sclib::atAction:
			return (type == sclib::atNoise) ? true : false;
		case sclib::atBreath:
			return (type == sclib::atNoise) ? true : false;
		case sclib::atUndefined:
			return (type == sclib::atNoise) ? true : false;
		case sclib::atSpeech:
			return false;
		case sclib::atPause:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atVoiced:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atUnvoiced:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atMaleVoice:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atFemaleVoice:
			return (type == sclib::atSpeech) ? true : false;
		case sclib::atSilence:
			return false;
		case sclib::atShort:
			return false;
		case sclib::atSceneBoundary:
			return false;
		case sclib::atSpeakerBoundary:
			return false;
		case sclib::atNoiseBoundary:
			return false;
		case sclib::atShotBoundary:
			return false;
		case sclib::atSpeechSegmentStart:
			return false;
		case sclib::atSpeechSegmentEnd:
			return false;
		case sclib::atArtificialBoundary:
			return false;
		default:
			return false;
	}
}

//====================================================================================================================
// Returns true if the first parameter is semantically the opposite of the second one (e.g. SPEECH and NOISE, meaning 
// the absence of one means the presence of the other), false if it isn't (e.g. BREATH and MUSIC, they don't have such 
// a struct relationship although the< are mutually exclusive)
//====================================================================================================================
bool SC_GroundTruth::isOppositeType(long int questionableType, long int type) {
	switch (questionableType) {
		case sclib::atNoise:
      return (type == sclib::atSpeech) ? true : false;
		case sclib::atPureSpeech:
			return (type == sclib::atNoisySpeech) ? true : false;
		case sclib::atNoisySpeech:
			return (type == sclib::atPureSpeech) ? true : false;
		case sclib::atBackground:
			return false;
		case sclib::atMusic:
			return false;
		case sclib::atAction:
			return false;
		case sclib::atBreath:
			return false;
		case sclib::atUndefined:
			return false;
		case sclib::atSpeech:
			return  (type == sclib::atNoise) ? true : false;;
		case sclib::atPause:
			return false;
		case sclib::atVoiced:
			return (type == sclib::atUnvoiced) ? true : false;
		case sclib::atUnvoiced:
			return (type == sclib::atVoiced) ? true : false;
		case sclib::atMaleVoice:
			return (type == sclib::atFemaleVoice) ? true : false;
		case sclib::atFemaleVoice:
			return (type == sclib::atMaleVoice) ? true : false;
		case sclib::atSilence:
      return (type == sclib::atSpeech || type == sclib::atNoise) ? true : false;
		case sclib::atShort:
			return false;
		case sclib::atSceneBoundary:
			return false;
		case sclib::atSpeakerBoundary:
			return false;
		case sclib::atNoiseBoundary:
			return false;
		case sclib::atShotBoundary:
			return false;
		case sclib::atSpeechSegmentStart:
			return false;
		case sclib::atSpeechSegmentEnd:
			return false;
		case sclib::atArtificialBoundary:
			return false;
		default:
			return false;
	}
}

//====================================================================================================================
// Returns true if the first parameter and the second parameters (can be OR-concatenated) are mutually exclusive (e.g. 
// VOICED and UNVOICED speech); in the parameter problem a OR-concatenated list of those types actally violating the
// mutual exclusiveness is given
//====================================================================================================================
bool SC_GroundTruth::areMutualExclusiveTypes(long int questionableType, long int types, long int &problem) {
	bool res = false;

	problem = sclib::noType;

	switch (questionableType) {
		case sclib::atNoise:
			//the term after '&' represents all the types that must not coincide with the questionable type; 
			//the '&'-operation with types results in all those audio-types that actually occur in types (thay are "problematic"); 
			//if this is > 0 (i.e. != sclib::noType), questionableType and types are mutually exclusive
			problem |= types & (sclib::atSilence|SC_GroundTruth::speechRelatedTypes);
			return problem != sclib::noType;
		case sclib::atPureSpeech:
			problem |= types & ((SC_GroundTruth::detailedAudioTypes^sclib::atPureSpeech)|sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atNoisySpeech:
			problem = types & ((SC_GroundTruth::detailedAudioTypes^sclib::atNoisySpeech)|sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atBackground:
			problem |= types & ((SC_GroundTruth::detailedAudioTypes^sclib::atBackground)|sclib::atSilence|SC_GroundTruth::speechRelatedTypes);
			return problem != sclib::noType;
		case sclib::atAction:
			problem |= types & ((SC_GroundTruth::detailedAudioTypes^sclib::atAction)|sclib::atSilence|SC_GroundTruth::speechRelatedTypes);
			return problem != sclib::noType;
		case sclib::atBreath:
			problem |= types & ((SC_GroundTruth::detailedAudioTypes^sclib::atBreath)|sclib::atSilence|SC_GroundTruth::speechRelatedTypes);
			return problem != sclib::noType;
		case sclib::atUndefined:
			problem |= types & ((SC_GroundTruth::detailedAudioTypes^sclib::atUndefined)|sclib::atSilence);
			return problem != sclib::noType;
		case sclib::atSpeech:
			problem |= types & (sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atPause:
			return false; //still unclear/open if it can also be used for interrupted noise segments (currently it is only used for interrupted speech segments)
		case sclib::atVoiced:
			problem |= types & (sclib::atUnvoiced|sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atUnvoiced:
			problem |= types & (sclib::atVoiced|sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atMaleVoice:
			problem |= types & (sclib::atFemaleVoice|sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atFemaleVoice:
			problem |= types & (sclib::atMaleVoice|sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atSilence:
			problem |= types & ((SC_GroundTruth::speechRelatedTypes^sclib::atSpeakerBoundary)|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atShort: //still unclear/open if it can also be used for short noise segments (currently it is only used for short speech segments)
			return false;
		case sclib::atSceneBoundary:
			return false;
		case sclib::atSpeakerBoundary:
			problem |= types & (sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atNoiseBoundary:
			problem |= types & (sclib::atSilence|SC_GroundTruth::speechRelatedTypes);
			return problem != sclib::noType;
		case sclib::atShotBoundary:
			return false;
		case sclib::atSpeechSegmentStart:
			problem |= types & (sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atSpeechSegmentEnd:
			problem |= types & (sclib::atSilence|SC_GroundTruth::noiseRelatedTypes);
			return problem != sclib::noType;
		case sclib::atArtificialBoundary:
			return false;
		default:
			return false;
	}

	return res;
}

//====================================================================================================================
// Returns a string constant representing the Name (or abbreviation) of/for the given type
//====================================================================================================================
const char* SC_GroundTruth::getAudioTypeName(long int audioType, bool shortName) {
	switch (audioType) {
		case sclib::atNoise:
      return (shortName == false) ? "Noise" : "NOIZ";
		case sclib::atPureSpeech:
      return (shortName == false) ? "Pure speech" : "pSPCH";
		case sclib::atNoisySpeech:
      return (shortName == false) ? "Noisy speech" : "nSPCH";
		case sclib::atBackground:
      return (shortName == false) ? "Background noise" : "BG";
		case sclib::atMusic:
      return (shortName == false) ? "Music" : "MUSIC";
		case sclib::atAction:
      return (shortName == false) ? "Action sound" : "ACTN";
		case sclib::atBreath:
      return (shortName == false) ? "Breating sound" : "BREATH";
		case sclib::atUndefined:
      return (shortName == false) ? "Undefined sound" : "UNDEF";
		case sclib::atSpeech:
      return (shortName == false) ? "Speech" : "SPCH";
		case sclib::atPause:
      return (shortName == false) ? "Pause" : "PAUSE";
		case sclib::atVoiced:
      return (shortName == false) ? "Voiced speech" : "VOICED";
		case sclib::atUnvoiced:
      return (shortName == false) ? "Unvoiced speech" : "UNVOIC";
		case sclib::atMaleVoice:
      return (shortName == false) ? "Male speech" : "MALE";
		case sclib::atFemaleVoice:
      return (shortName == false) ? "Female speech" : "FEMALE";
		case sclib::atSilence:
      return (shortName == false) ? "Silence" : "SILNC";
		case sclib::atShort:
      return (shortName == false) ? "Short segment" : "SHORT";
		case sclib::atSceneBoundary:
      return (shortName == false) ? "Scene boundary" : "SCENE_B";
		case sclib::atSpeakerBoundary:
      return (shortName == false) ? "Speaker boundary" : "SPK_B";
		case sclib::atNoiseBoundary:
      return (shortName == false) ? "Noise boundary" : "NOIZ_B";
		case sclib::atShotBoundary:
      return (shortName == false) ? "Shot boundary" : "SHOT_B";
		case sclib::atSpeechSegmentStart:
      return (shortName == false) ? "Speaker's segment start" : "START";
		case sclib::atSpeechSegmentEnd:
      return (shortName == false) ? "Speaker's segment end" : "END";
		case sclib::atArtificialBoundary:
      return (shortName == false) ? "Artificial boundary" : "ART_B";
		default:
			return "";
	}	
}

//====================================================================================================================
// Fills the buffer (and returns a pointer to it) with the concatenated (separated by the given separator) names of 
// the types present in audioTypes (may be OR-concatenated)
//====================================================================================================================
char* SC_GroundTruth::getAudioTypesNames(long int audioTypes, char *buffer, const char *separator, bool shortNames) {
	int nrOfTypes = sizeof(audioTypes)*8 - 1;
	long int t = 1;
	bool firstEntry = true;
	char tmpBuffer[sclib::bufferSize];

	sprintf(buffer, "\0");

	for (int i = 0; i < nrOfTypes; i++) {
		if (sclib::bitTest(audioTypes, t) == true) {
			if (firstEntry == true) { //don't write the separator (which will otherwise be the last character in this string) for the first entry
				sprintf(buffer, "%s", getAudioTypeName(t, shortNames));
				firstEntry = false;
			} else { //write new entries before the previous one(s)
				sprintf(tmpBuffer, "%s%s%s", getAudioTypeName(t, shortNames), separator, buffer);
				sprintf(buffer, "%s", tmpBuffer);
			}
		}
		t = t << 1;
	}

	return buffer;
}

//====================================================================================================================
// File-I/O for this class, so that the sate of a groundtruth-object can be saved to a file
// All members except the pSignalPrototype and the pTweak are saved
//====================================================================================================================
bool SC_GroundTruth::save(const char *fileName) {
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

	//types where length-information is necessary before the actual data
	if (this->audioFileName != NULL) {
		len = (unsigned long int)(strlen(this->audioFileName));
		bytes += io.writeScalar(&gtFile, len);
		bytes += io.writeArray(&gtFile, this->audioFileName, len);
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
		REPORT_ERROR(SVLIB_Fail, "Saving SC_GroundTruth Failed!");
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
SC_GroundTruth* SC_GroundTruth::load(const char *fileName) {
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

	//types where length-information is necessary before the actual data
	MFree_1D(this->audioFileName);
	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	if (len > 0) {
		MArray_1D(this->audioFileName, len+1, char, "SC_GroundTruth.load: audioFileName");
		bytes += io.readArray(&gtFile, this->audioFileName, len, codeSizes, fileSizes);
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
		REPORT_ERROR(SVLIB_Fail, "Loading SC_GroundTruth Failed!");
		pGT = NULL;
	}

	gtFile.close();

	return pGT;
}

//====================================================================================================================
// This method returns a (gt) speaker-id each time it is called until it returns false when all speaker-ids are 
// returned; with a nested loop, over all speakers, therein over all segments, therein doing loading, feature 
// extraction and modeling, speaker models for all speakers in the corpus can be built.
//====================================================================================================================
bool SC_GroundTruth::getSpeakersCallback(long int &speakerId) {
	static long int lastSpeakerId = sclib::noSpeaker;
	bool res = false;

	speakerId = sclib::noSpeaker; //initialization in case of no further speaker

	//this method relies on the fact that new speakers are added subsequently to the speakerNames array in order of their appearance.
	//if an entry in the array has the value NULL, this means that no further speakers are present.
	if (lastSpeakerId != sclib::noSpeaker) {
		if (lastSpeakerId+1 < sclib::maxSpeakers && this->speakerNames[lastSpeakerId+1] != NULL) {
			speakerId = lastSpeakerId + 1;
			lastSpeakerId = speakerId;
			res = true;
		}
	} else {
		if (this->speakerNames[0] != NULL) {
			speakerId = 0;
			lastSpeakerId = 0;
			res = true;
		}
	}

	return res;
}

//====================================================================================================================
// This method gives, for a given (gt) speaker-id, all segments belonging to this speaker; with a nested loop, over 
// all speakers, therein over all segments, therein doing loading, feature extraction and modeling, speaker models for 
// all speakers in the corpus can be built.
//====================================================================================================================
bool SC_GroundTruth::getSpeakersSegmentsCallback(long int speakerId, long int &segmentStart, long int &segmentEnd) {
	static long int lastSpeakerId = sclib::noSpeaker;
	static long int lastSegmentEnd = sclib::noSegment;
	bool res;
	unsigned long int initPosition;
	long int currentSpeakerId;

	if (speakerId == sclib::noSpeaker) { //a trivial case at the beginning
		lastSpeakerId = sclib::noSpeaker;
		lastSegmentEnd = sclib::noSegment;
		return false;
	}

	initPosition = (lastSegmentEnd == sclib::noSegment || lastSpeakerId != speakerId) ? 0 : lastSegmentEnd + 1;
	do { //search the speaker-homogenious segments from the last one on to find the next one of the given speaker
		getNextBoundary(initPosition, segmentStart, segmentEnd, sclib::atSpeakerBoundary, sclib::searchForward, sclib::modeGroundtruth);
		currentSpeakerId = getMajorSpeakerID(segmentStart, segmentEnd, sclib::modeGroundtruth); //checks inside if the segment is valid, otherwise returns sclib::noSpeaker
		initPosition = segmentEnd + 1; //for the next loop
	} while (currentSpeakerId != speakerId && currentSpeakerId != sclib::noSpeaker); //this is only true if the segment is also valid, see above

	if (currentSpeakerId == speakerId) { //we found another valid segment of the given speaker
		lastSpeakerId = speakerId;
		lastSegmentEnd = segmentEnd;
		res = true;
	} else { //there are no further segments of the given speaker
		segmentStart = sclib::noSegment;
		segmentEnd = sclib::noSegment;
		lastSpeakerId = sclib::noSpeaker;
		lastSegmentEnd = sclib::noSegment;
		res = false;
	}

	return res;
}

//====================================================================================================================
//  Calculate the Diarization error Rate (DER) according to Huang, Marchertet, Visweswariah, Potamianos, "The IBM RT07 
//  Evaluation System for Speaker Diarization on Lecutre Meetings", 2007:
//    "In accordance to NIST scoring, results are reported in terms of diarization error rate (DER). DER is calculated 
//     by first finding the optimal one-to-one mapping between reference speakers and the hypothesized ones, and then 
//     computing the percentage of time that is wrongly assigned according to the optimal mapping. DER includes 
//     speaker error time, missed speaker time, and false alarm speaker time, thus also taking SAD errors into 
//     account" (SAD=speech activity detection)
//  It is placed here to save lots of computational time when the getNext* methods are used as compared to just
//  traversing the internal data strucuture once herein
//====================================================================================================================
double SC_GroundTruth::calcDER(unsigned long int softBoundaryDiameter) {
	unsigned long int wrong = 0, sbdFrames = this->sample2FLI(softBoundaryDiameter);
	long int gtId, hypoId, start, end, speechStart, speechEnd;
	std::map<long int, std::pair<long int, bool> > idMap;
	std::map<long int, std::pair<long int, bool> >::iterator pMap;
	std::pair<long int, bool> gtIdCorrect;
  SC_GroundTruth::SC_SpeakerMapping *pHook = this->pSpeakerMapping;

	//create a faster way to access groundtruth speaker id's from given hypthesized speaker ids 
	//before calling getSpeakerGIDFromHID() (and each time traversing a linked list) internalFrameCount times
  while (pHook != NULL) {
		gtIdCorrect = std::make_pair(pHook->groundTruthID, pHook->correct);
		idMap[pHook->hypothesizedID] = gtIdCorrect;
    pHook = pHook->Next;
  }

	//get boundaries of first speech segment (and first non-speech segment as well)
	speechStart = getNextFrame(0, sclib::atSpeechSegmentStart, true, sclib::noType, true, sclib::searchForward, sclib::modeHypothesized); //get the boundaries for the first speech segment: segment-start...
	speechEnd = getNextFrame(speechStart, sclib::atSpeechSegmentEnd, true, sclib::noType, true, sclib::searchForward, sclib::modeHypothesized); //...and the respective segment-end
	start = 0; //this are the boundaries of the first nonspeech-segment, then: it just preceeds the first speech segment
	end = speechStart - 1;

	//count all wrong (in the 3-fold sense above) frames
	for (unsigned long int f = 0; f < this->internalFrameCount; f++) {
		if ((long int)(f) > end) { //if the current speech- or nonspeech-segment is exhausted, look for the next one:
			if (end < speechEnd) { //we where in a nonspeech-segment, so we already now the boundaries of the next speech segment
				start = speechStart;
				end = speechEnd;
			} else { //the current speech segment is over, so get the boundaries of the next one but first account for the noise-segment in between
				start = speechEnd + 1;
				speechStart = getNextFrame(f, sclib::atSpeechSegmentStart, true, sclib::noType, true, sclib::searchForward, sclib::modeHypothesized); //get the boundaries for the first speech segment: segment-start...
				if (speechStart == sclib::noSegment) {
					end = this->internalFrameCount - 1;
				} else {
					speechEnd = getNextFrame(speechStart+1, sclib::atSpeechSegmentEnd, true, sclib::noType, true, sclib::searchForward, sclib::modeHypothesized); //...and the respective segment-end
					end = speechStart - 1;
				}
			}
		}
		
		if (sclib::isBetween(start+sbdFrames, f, end-sbdFrames) == true) { //account for soft boundaries by not counting any wrong decisions made in the diameter around noise/speech segment boundaries
			gtId = getSpeakerID(f, sclib::modeGroundtruth);
			hypoId = getSpeakerID(f, sclib::modeHypothesized);
			pMap = idMap.find(hypoId);
			if (pMap == idMap.end()) { //we must distinguish the case of not finding this hypoId from the normal case because in the case of not finding, we don't want the standard return value of std::map (which is <0,false>), but <sclib::noSpeaker,false>
				gtIdCorrect = std::make_pair(sclib::noSpeaker, false);
			} else {
				gtIdCorrect = std::make_pair(pMap->second.first, pMap->second.second);
			}
			
			if ((hypoId!=sclib::noSpeaker && gtId==sclib::noSpeaker) ||
					(hypoId==sclib::noSpeaker && gtId!=sclib::noSpeaker) ||
					(gtIdCorrect.first!=gtId) ||
					(gtIdCorrect.first!=sclib::noSpeaker && gtIdCorrect.first==gtId && gtIdCorrect.second == false)) {
				wrong++;
			}
		}
	}

	return ((double)(wrong)*(double)(this->internalFrameSize)) / (double)(this->audioSampleCount);
}
