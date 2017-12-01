/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth to accomplish its targets in the    */
/*        context of TIMIT-Corpus files (LDC's pure speech corpus)        */
/*      - reads and organizes and provides access to TIMIT annotation/     */
/*        ground-truth files                                              */
/*      - makes parts of the corpus available as if they are one big file */
/*      - operates on lists of TIMIT-files specifiying the "parts"        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 07.04.2006																								*/
/**************************************************************************/

#include <math.h>
#include <string.h>
#include <iomanip>
#include <map>
#include <algorithm>
#include <string>
#include "SC_GroundTruth_TIMIT.h"
#include "SC_SignalHandler.h"
#include "SC_Classifier.h"
#include <SV_DataIO.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth_TIMIT::SC_GroundTruth_TIMIT(SC_TweakableParameters *pTweak, const char* corpusFileName, bool verbose, bool setSceneBoundariesAtFileBoundaries, bool setSpeakerBoundariesAtFileBoundaries, unsigned int sbSkipCount, bool setSceneBoundariesAtSpeakerBoundaries) : SC_GroundTruth(pTweak, NULL) {
  //some init-stuff
  this->gtType = sclib::gtTIMIT;
  this->audioFileName = new char[strlen(corpusFileName)+1]; //yes, a missuse of this member, but seems ok ;-)
	strcpy(this->audioFileName, corpusFileName);
  this->fileList = NULL;
  this->verbose = verbose;
	this->uncertaintyRegion = 0;
	this->setSceneBoundariesAtFileBoundaries = setSceneBoundariesAtFileBoundaries;
	this->setSpeakerBoundariesAtFileBoundaries = setSpeakerBoundariesAtFileBoundaries;
	this->sbSkipCount = sbSkipCount;
	this->setSceneBoundariesAtSpeakerBoundaries = setSceneBoundariesAtSpeakerBoundaries;

  if (measureCorpusSize(corpusFileName) == true) {
    initFrameList();

    //read the ground-truth files
    readGroundTruth();
  }
}

//====================================================================================================================
//	destructor 
//====================================================================================================================
SC_GroundTruth_TIMIT::~SC_GroundTruth_TIMIT() {
  sclib::destructLinkedList(this->fileList);
}

//====================================================================================================================
//	initialize the framelist
//====================================================================================================================
void SC_GroundTruth_TIMIT::initFrameList(void) {
  MArray_2D(this->frameList, (long int)(this->internalFrameCount), 6, long int, "initFrameList: frameList"); 
	if (this->pTweak->groundTruth.storeProbabilityInformation == true) {
		MArray_2D(this->probabilityList, (long int)(this->internalFrameCount), getProbabilityListDim(), double, "initFrameList: probabilityList"); 
	}

	//init frameList with zeros for the labels and -1 for the speaker-id's
  for (unsigned long int y = 0; y < this->internalFrameCount; y++) {
		this->frameList[y][0] = sclib::noType;  //ground-truth frame-attributes
    this->frameList[y][1] = sclib::noType;  //hypothesized frame-attributes
		this->frameList[y][2] = sclib::noSpeaker; //ground-truth speaker-id
    this->frameList[y][3] = sclib::noSpeaker; //hypothesized speaker-id
    this->frameList[y][4] = sclib::noPhone; //ground-truth phone-type 
    this->frameList[y][5] = sclib::noPhone; //hypothesized phone-type 
		if (this->pTweak->groundTruth.storeProbabilityInformation == true) {
			for (unsigned int x = 0; x < getProbabilityListDim(); x++) {
				this->probabilityList[y][x] = 1.0;
			}
		}
	}

	return;
}

//====================================================================================================================
// This class holds the groundtruth for a whole corpus as if it where a single file; this function reads the 
// parameters of all files in the corpus to measure it's overall size and create a mapping between the frames in 
// the framelist and the corresponding files storing the actual signals
//====================================================================================================================
bool SC_GroundTruth_TIMIT::measureCorpusSize(const char *corpusFileName) {
  unsigned long int overallSampleCount = 0;
	bool res = true;
  int count = 0;
	char line[sclib::bufferSize], *waveformFileName, errStr[sclib::bufferSize], *completeFileName = NULL, *baseDir = NULL;
  FILE *corpusFile;

	SC_SignalHandler signalHandler(this->pTweak, ((sclib::like(corpusFileName, "%conTIMIT%")==true)?sclib::stWave:sclib::stNIST)); //conTIMIT data has RIFF/Wave files, TIMIT is NIST format
  SC_Signal *pSignal = NULL;

  if (this->verbose == true) {
    printf("Measuring TIMIT corpus size for corpus file %s: ", corpusFileName);
  }

  baseDir = sclib::extractPath(corpusFileName);

  if ((corpusFile = fopen(corpusFileName, "r")) != NULL) {
		while (!feof(corpusFile)) {
			//get the filename, discard what may follow after whitespaces
			if (sclib::readline(corpusFile, line, sclib::bufferSize) > 1) {
				waveformFileName = sclib::getNextStringFromString(line, sclib::bufferSize);

				//construct complete waveformFileName
				MArray_1D(completeFileName, strlen(baseDir) + strlen(waveformFileName) + 1, char, "SC_GroundTruth_TIMIT.measureCorpusSize: completeFileName");
				sprintf(completeFileName, "%s%s\0", baseDir, waveformFileName);
				MFree_1D(waveformFileName);

				//open signal and store it's partameters
				pSignal = signalHandler.openSignal(completeFileName, this->pTweak->signalHandler.forceSampleRate);
				if (pSignal == NULL) {
					sprintf(errStr, "Waveform-file %s could not be opened", completeFileName);
					REPORT_ERROR(SVLIB_FileErr, errStr);
					res = false;
				} else {
					addFile(completeFileName, overallSampleCount, overallSampleCount+pSignal->getSampleCount()-1); //remember which files are mapped to what sample-positions
					overallSampleCount += pSignal->getSampleCount();
					if (count == 0) {
						this->pSignalPrototype = pSignal;
					} else {
						if (this->pSignalPrototype->SigPar.Encode != pSignal->SigPar.Encode || this->pSignalPrototype->SigPar.NChannel != pSignal->SigPar.NChannel || this->pSignalPrototype->SigPar.SRate != pSignal->SigPar.SRate) { // || this->pSignalPrototype->SigPar.StByte != pSignal->SigPar.StByte) {
							sprintf(errStr, "Corpus-File %s has different signal parameters than the rest", completeFileName);
							REPORT_ERROR(SVLIB_FileErr, errStr);
							res = false;
						}
						MFree_0D(pSignal);
					}
				}

				MFree_1D(completeFileName);
				count++;
			} //if corpus line not empty
    }

    fclose(corpusFile);

    this->audioSampleRate = this->pSignalPrototype->SigPar.SRate;
    this->audioSampleCount = overallSampleCount;
		this->pConverter->setAudioSampleRate(this->audioSampleRate);
	  this->internalFrameSize = this->pConverter->ms2sample(this->pTweak->groundTruth.internalFrameSize); //in samples!!!
    this->internalFrameCount = (unsigned long int)(ceil((double)(this->audioSampleCount) / (double)(this->internalFrameSize)));

  } else {
    res = false;
  }

  MFree_1D(baseDir);

  if (this->verbose == true) {
		printf("done (size: %0.1f min, %d files)!\n", this->pConverter->sample2ms(overallSampleCount, this->pSignalPrototype->SigPar.SRate)/1000.0/60.0, count);
  }

  return res;
}

//====================================================================================================================
// Read Groundtruth from file(s) and initialize internal data structures (frameList etc.)
//====================================================================================================================
bool SC_GroundTruth_TIMIT::readGroundTruth(void) {
  unsigned long int start = 0, end = 0, offset = 0, types, phoneType, oldEnd, fileLength;
  bool res = true, firstSpeech;
  int spkID = sclib::noSpeaker, oldSpkID;
	char line[sclib::bufferSize], *waveformFileName, errStr[sclib::bufferSize], *tmp = NULL, *completeFileName = NULL, *baseDir = NULL, *fileName, *pathName;
	char sex, speaker[sclib::bufferSize], phone[sclib::bufferSize], *lastSlash;
  FILE *corpusFile, *phoneFile;
  SC_SignalHandler signalHandler(this->pTweak, sclib::stNIST);
  SC_Signal *pSignal = NULL;
  SC_GroundTruth_TIMIT::SC_FileList *pFileListHook = this->fileList;
	double lastPercentage = 0.0;
	double sampleRateRatio = (this->pTweak->signalHandler.forceSampleRate > 0) ? (double)(this->pTweak->signalHandler.forceSampleRate)/16000.0 : 1.0;
	int speakerFileCount;
	bool isConTimit = false, changeBetweenFiles;
	float sec;

  if (this->verbose == true) {
    printf("Loading ground truth for corpus file %s: ", this->audioFileName);
  }

  baseDir = sclib::extractPath(this->audioFileName);
	if (sclib::like(baseDir, "%conTIMIT%") == true) { //we can also handle the conTIMIT dataset from Kotti, Benetos, Kotropoulos, "Computationally Efficient and Robust BIC-Based Speaker Segmentation", 2008 with this class; it has different groundtruth, though
		isConTimit = true;
		setSegment(0, this->audioSampleCount-1, sclib::atSpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false); //only speech in conTIMIT
		setSegment(0, this->audioSampleCount-1, sclib::atSpeech, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
	}

  if ((corpusFile = fopen(this->audioFileName, "r")) != NULL) {
    while (!feof(corpusFile)) {
			if (sclib::readline(corpusFile, line, sclib::bufferSize) > 1) {
				waveformFileName = sclib::getNextStringFromString(line, sclib::bufferSize);
				
				if (isConTimit == true) { //conTIMIT: concatenated TIMIT files with scd groundtruth; see http://poseidon.csd.auth.gr/LAB_RESEARCH/Latest/data or Kotti et al., "Computationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
					changeBetweenFiles = (sclib::getNextIntFromString(line, sclib::bufferSize) == 0) ? false : true; //conTIMIT corpus files may indicate by an integer ==0 after a given filename that there is NO speaker change between this file and the one in the next line

					//construct complete filename
					fileName = sclib::extractFileName(waveformFileName);
					pathName = sclib::extractPath(waveformFileName, true);
					MArray_1D(tmp, strlen(baseDir) + strlen(pathName) + strlen("../groundtruth/") + strlen(fileName) + 1, char, "SC_GroundTruth_TIMIT.readGroundTruth: tmp");
					sprintf(tmp, "%s%s%s%s\0", baseDir, pathName, "../groundtruth/", fileName); //groundtruth for contimit resides in a different directory
					completeFileName = sclib::exchangeFileExtension(tmp, ".txt");
					MFree_1D(tmp);
					MFree_1D(fileName);
					MFree_1D(pathName);

					//read the speaker change ground truth of conTIMIT
					offset = pFileListHook->startSample;
					if ((phoneFile = fopen(completeFileName, "r")) != NULL) {
						//set scene-boundaries at file-bouhdaries, because we do not have any way to be sure about the speakerchange indicated at each new file-beginning by the groundtruth
						//=> conTIMIT data has to processed scene-wise normally!!!
						if (this->setSceneBoundariesAtFileBoundaries == true) {
							setSegment(offset, offset, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false);
							setSegment(offset, offset, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
						}
						setSegment(pFileListHook->startSample, pFileListHook->endSample, sclib::atSpeech|sclib::atVoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false);
						setSegment(pFileListHook->startSample, pFileListHook->endSample, sclib::atSpeech|sclib::atVoiced, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
						if (changeBetweenFiles == true && pFileListHook->endSample < this->audioSampleCount-1) { //typically the speaker changes from file to file if this is not the last file in the corpus
							setSegment(pFileListHook->endSample+1, pFileListHook->endSample+1, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false);
							setSegment(pFileListHook->endSample+1, pFileListHook->endSample+1, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
							if (this->setSceneBoundariesAtSpeakerBoundaries == true) {
								setSegment(offset, offset, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false);
								setSegment(offset, offset, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
							}
						}
						while (fscanf(phoneFile, "%f ", &sec) != EOF) { 
							start = offset + this->pConverter->ms2sample(sclib::round(sec*1000.0)); //conTIMIT groundtruth is given in terms of seconds
							if (start > 0) { //this would only mark the beginning of the first segment, not an actual change point
								setSegment(start, start, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false);
								setSegment(start, start, sclib::atSpeakerBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
							}
						}
						fclose(phoneFile);
					} else {
						sprintf(errStr, "Segmentation-file %s could not be opened", completeFileName);
						REPORT_ERROR(SVLIB_FileErr, errStr);
						res = false;
					}

				} else {
					//construct complete complete filename
					MArray_1D(tmp, strlen(baseDir) + strlen(waveformFileName) + 1, char, "SC_GroundTruth_TIMIT.readGroundTruth: tmp");
					sprintf(tmp, "%s%s\0", baseDir, waveformFileName);
					completeFileName = sclib::exchangeFileExtension(tmp, ".PHN");
					MFree_1D(tmp);

					//extract speaker-name and sex
					tmp = sclib::extractPath(completeFileName, false);
					lastSlash = strrchr(tmp, '/'); //separates the last directory name which contains information on the speaker's sex and name
					sex = lastSlash[1]; //[0] is '/', the next character specifies males 'M' and females 'F'
					sprintf(speaker, "%s", lastSlash+2); //following the sex-specifier is the speaker's "name", usually 4 characters long
					MFree_1D(tmp);
					oldSpkID = spkID; //remember the last speaker to set speaker boundaries correctly
					spkID = insertSpeakerName(speaker, (int)(strlen(speaker)));

					if (oldSpkID == spkID) { //how many files of this speaker have been read in a row?
						speakerFileCount++;
					} else {
						speakerFileCount = 1;
					}

					//open phone-file
					offset = pFileListHook->startSample;
					if ((phoneFile = fopen(completeFileName, "r")) != NULL) {
						oldEnd = 0;
						firstSpeech = true; //first speech part in each file
						while (fscanf(phoneFile, "%d %d %s", &start, &end, phone) != EOF) {
							start = sclib::round(start * sampleRateRatio); //to compensate for differences between the original sample-boundaries in the gt-files (based on 16kHz data) and the forced samplerate of the loaded signal
							end = sclib::round(end * sampleRateRatio);
							
							end -= 1; //the endings in the phone-file are exclusive, while herein they are inclusive, thus decrement by 1

							//correct for errors in the .PHN files 
							if (oldEnd == 0 && start > 0) { //the transcript doesn't start with offset 0, so include silence part, scene/shot-boundary and silence-phone
								if (this->setSceneBoundariesAtFileBoundaries == true) {
									setSegment(offset, offset, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
									setSegment(offset, offset, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
								}
								setPhone(offset, offset+start-1, sclib::phone_h_sharp, sclib::modeGroundtruth);
								setSegment(offset, offset+start-1, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, true); //also remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
								setPhone(offset, offset+start-1, sclib::phone_h_sharp, sclib::modeHypothesized);
								setSegment(offset, offset+start-1, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //also remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
							}
							if (end > pFileListHook->endSample-pFileListHook->startSample+1) { //the end is too long, so cut it down
								end = pFileListHook->endSample-pFileListHook->startSample+1;
							}
							if (oldEnd > 0 && start > 0 && start != oldEnd+1) { //there are gaps between consecutive labels, so fill them with silence, too
								setPhone(offset+oldEnd+1, offset+start-1, sclib::phone_h_sharp, sclib::modeGroundtruth);
								setSegment(offset+oldEnd+1, offset+start-1, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, true); //also remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
								setPhone(offset+oldEnd+1, offset+start-1, sclib::phone_h_sharp, sclib::modeHypothesized);
								setSegment(offset+oldEnd+1, offset+start-1, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //also remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
							}

							//determine audio-types
							types = phone2audioType(phone, phoneType);
							if (sclib::bitTest(types, sclib::atSpeech) == true) {
								types |= ((sex == 'F' || sex == 'f') ? sclib::atFemaleVoice : sclib::atMaleVoice);
							}

							//error during phone-translation?
							if (phoneType == sclib::noType) {
								printf("Phone %s in transcription-file %s not recognized!\n", phone, completeFileName);
							}

							//store the findings in the frameList
							setPhone(offset+start, offset+end, phoneType, sclib::modeGroundtruth);
							setSegment(offset+start, offset+end, types, false, spkID, sclib::modeLabelAdd, sclib::modeGroundtruth, false, true); //also remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
							setPhone(offset+start, offset+end, phoneType, sclib::modeHypothesized);
							setSegment(offset+start, offset+end, types, false, spkID, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //also remove previously set types that may exist here (at the borders between two phones) due to frame-quantization

							//care for scene boundaries (at the beginning of new files, but maybe not each new file...)
							if ((start == 0) && //it is a beginning of a file
								  ( (this->setSceneBoundariesAtFileBoundaries == true) ||  //we set a scene boundary at every beginning OR...
									  ( (this->setSceneBoundariesAtSpeakerBoundaries == true) && //we set it a file boundaries where the speaker id changes
									    ((oldSpkID != spkID) || ((this->sbSkipCount>0) && (speakerFileCount%(this->sbSkipCount+1)==1))) //this makes sure that a speaker boundary will be set at the first speech occurence in this file
								 )) ) {
								setSegment(offset+start, offset+start, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
								setSegment(offset+start, offset+start, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
							}

							//care for speaker boundaries
							//explanation for the "(this->sbSkipCount>0)?speakerFileCount%(this->sbSkipCount+1)==1:true"-term: 
							//if some consecutive speaker boundaries should be skipped (in order to create longer speaker segments), 
							//set a boundarie at the beginning of each (this->sbSkipCount+2)th file, i.e. (this->sbSkipCount+1) files are merged; 
							//example: you want 2 files being merged, then set this->sbSkipCount==1, then boundaries will be set at the beginning of the 1st, 3rd 6th, 9th, ... file of the same speaker.
							if ((sclib::bitTest(types, sclib::atSpeech) == true) && 
								  (firstSpeech == true) && 
									( (oldSpkID != spkID) || 
									  (this->setSceneBoundariesAtFileBoundaries==false && this->setSpeakerBoundariesAtFileBoundaries==true && ((this->sbSkipCount>0)?(speakerFileCount%(this->sbSkipCount+1)==1):true))
								 )) {
								firstSpeech = false;
								setSegment(offset+start, offset+start, sclib::atSpeakerBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeGroundtruth);
								setSegment(offset+start, offset+start, sclib::atSpeakerBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeHypothesized);
								if (oldSpkID == sclib::noSpeaker) {
									setSegment(offset+start, offset+start, sclib::atArtificialBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeGroundtruth);
									setSegment(offset+start, offset+start, sclib::atArtificialBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeHypothesized);
								}
							}

							oldEnd = end;
						}

						//correct for errors in the .PHN files 
						//the transcript is too short, so fill the rest with silence phones/audio-types
						fileLength = pFileListHook->endSample - pFileListHook->startSample + 1;
						if (end < fileLength-1) {
							setPhone(offset+end+1, offset+fileLength-1, sclib::phone_h_sharp, sclib::modeGroundtruth);
							setSegment(offset+end+1, offset+fileLength-1, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, true); //remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
							setPhone(offset+end+1, offset+fileLength-1, sclib::phone_h_sharp, sclib::modeHypothesized);
							setSegment(offset+end+1, offset+fileLength-1, sclib::atSilence, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, true); //remove previously set types that may exist here (at the borders between two phones) due to frame-quantization
						}

						fclose(phoneFile);
					} else { //if gt-files can't be found, assume that everythin gin the files is speech and care for scene/speaker-boundaries
						printf("\nphone-file '%s' could not be opened, assuming only speech in corresponding TIMIT audio file", completeFileName);
						start = pFileListHook->startSample; //this->pConverter->ms2sample(this->pConverter->sample2ms(pFileListHook->startSample)); //this way, boundaries are set at the beginnings of ms, not at a sample in between as with the oure file sample counts
						end = pFileListHook->endSample; //this->pConverter->ms2sample(this->pConverter->sample2ms(pFileListHook->endSample), sclib::alignmentEnd);
						setSegment(start, end, sclib::atSpeech, false, spkID, sclib::modeLabelAdd, sclib::modeGroundtruth);
						setSegment(start, end, sclib::atSpeech, false, spkID, sclib::modeLabelAdd, sclib::modeHypothesized);
						if ((this->setSceneBoundariesAtFileBoundaries == true) ||  //we set a scene boundary at every beginning OR...
								  ( (this->setSceneBoundariesAtSpeakerBoundaries == true) && //we set it a file boundaries where the speaker id changes
								    ((oldSpkID != spkID) || ((this->sbSkipCount>0) && (speakerFileCount%(this->sbSkipCount+1)==1))) //this makes sure that a speaker boundary will be set at the first speech occurence in this file
							  ) ) {
							setSegment(start, start, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
							setSegment(start, start, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
						}
						if ((oldSpkID != spkID) || 
								(this->setSceneBoundariesAtFileBoundaries==false && this->setSpeakerBoundariesAtFileBoundaries==true && ((this->sbSkipCount>0)?(speakerFileCount%(this->sbSkipCount+1)==1):true))
							 ) {
							setSegment(start, start, sclib::atSpeakerBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeGroundtruth);
							setSegment(start, start, sclib::atSpeakerBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeHypothesized);
							if (oldSpkID == sclib::noSpeaker) {
								setSegment(start, start, sclib::atArtificialBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeGroundtruth);
								setSegment(start, start, sclib::atArtificialBoundary, false, spkID, sclib::modeLabelAdd, sclib::modeHypothesized);
							}
						}
						//sprintf(errStr, "Phone-file %s could not be opened", completeFileName);
						//REPORT_ERROR(SVLIB_FileErr, errStr);
						res = false;
					}
				} //not conTIMIT, but real TIMIT corpus file

				if (this->verbose == true) {
					lastPercentage = sclib::printPercentage(this->audioSampleCount, offset, lastPercentage, 10.0, !offset);
				}

				MFree_1D(waveformFileName);
				MFree_1D(completeFileName);
				pFileListHook = pFileListHook->Next;
			} //corpus file line not emtpy
		}  //while corpusFile

    fclose(corpusFile);

		if (this->setSceneBoundariesAtFileBoundaries == false) { //include a first scene in every case
			setSegment(0, 0, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
			setSegment(0, 0, sclib::atSceneBoundary|sclib::atShotBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
		}

  } else {
    res = false;
  }

  MFree_1D(baseDir);

  if (this->verbose == true) {
		sclib::printPercentage(1, 1, 0.0, 0.0, false);
    printf("done!\n");
  }

  return res;
}

//====================================================================================================================
// Extracts the information regarding audio-types from the given phone; also returns a precise (numerical) 
// representation of the phone-type, if a variable to alter is given (!= NULL)
//====================================================================================================================
unsigned long int SC_GroundTruth_TIMIT::phone2audioType(const char* phone, unsigned long int &phoneType) {
  unsigned long int phoneNum = sclib::noType, type = sclib::noType;

	if (sclib::like(phone, "b") == true) {
    phoneNum = sclib::phone_b;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "bcl") == true) {
    phoneNum = sclib::phone_bcl;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "d") == true) {
    phoneNum = sclib::phone_d;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "dcl") == true) {
    phoneNum = sclib::phone_dcl;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "g") == true) {
    phoneNum = sclib::phone_g;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "gcl") == true) {
    phoneNum = sclib::phone_gcl;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "p") == true) {
    phoneNum = sclib::phone_p;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "pcl") == true) {
    phoneNum = sclib::phone_pcl;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "t") == true) {
    phoneNum = sclib::phone_t;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "tck") == true) {
    phoneNum = sclib::phone_tck;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "k") == true) {
    phoneNum = sclib::phone_k;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "kcl") == true) {
    phoneNum = sclib::phone_kcl;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "dx") == true) {
    phoneNum = sclib::phone_dx;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech; //was Uv
  } else if (sclib::like(phone, "tcl") == true) {
    phoneNum = sclib::phone_tcl;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "q") == true) {
    phoneNum = sclib::phone_q;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech; //was V
  } else if (sclib::like(phone, "jh") == true) {
    phoneNum = sclib::phone_jh;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ch") == true) {
    phoneNum = sclib::phone_ch;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "s") == true) {
    phoneNum = sclib::phone_s;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "sh") == true) {
    phoneNum = sclib::phone_sh;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "z") == true) {
    phoneNum = sclib::phone_z;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "zh") == true) {
    phoneNum = sclib::phone_zh;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "f") == true) {
    phoneNum = sclib::phone_f;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "th") == true) {
    phoneNum = sclib::phone_th;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "v") == true) {
    phoneNum = sclib::phone_v;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "dh") == true) {
    phoneNum = sclib::phone_dh;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "m") == true) {
    phoneNum = sclib::phone_m;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech; //was V
  } else if (sclib::like(phone, "n") == true) {
    phoneNum = sclib::phone_n;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech; //was V
  } else if (sclib::like(phone, "ng") == true) {
    phoneNum = sclib::phone_ng;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech; // was V
  } else if (sclib::like(phone, "em") == true) {
    phoneNum = sclib::phone_em;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech; //was Uv
  } else if (sclib::like(phone, "en") == true) {
    phoneNum = sclib::phone_en;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech; // was Uv
  } else if (sclib::like(phone, "eng") == true) {
    phoneNum = sclib::phone_eng;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech; //was Uv
  } else if (sclib::like(phone, "nx") == true) {
    phoneNum = sclib::phone_nx;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech; //was Uv
  } else if (sclib::like(phone, "l") == true) {
    phoneNum = sclib::phone_l;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "r") == true) {
    phoneNum = sclib::phone_r;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "w") == true) {
    phoneNum = sclib::phone_w;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "y") == true) {
    phoneNum = sclib::phone_y;//bing
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "hh") == true) {
    phoneNum = sclib::phone_hh;
    type = sclib::atSpeech|sclib::atUnvoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "hv") == true) {
    phoneNum = sclib::phone_hv;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "el") == true) {
    phoneNum = sclib::phone_el;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech; //was Uv
  } else if (sclib::like(phone, "iy") == true) {
    phoneNum = sclib::phone_iy;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ih") == true) {
    phoneNum = sclib::phone_ih;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "eh") == true) {
    phoneNum = sclib::phone_eh;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ey") == true) {
    phoneNum = sclib::phone_ey;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ae") == true) {
    phoneNum = sclib::phone_ae;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "aa") == true) {
    phoneNum = sclib::phone_aa;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "aw") == true) {
    phoneNum = sclib::phone_aw;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ay") == true) {
    phoneNum = sclib::phone_ay;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ah") == true) {
    phoneNum = sclib::phone_ah;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ao") == true) {
    phoneNum = sclib::phone_ao;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "oy") == true) {
    phoneNum = sclib::phone_oy;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ow") == true) {
    phoneNum = sclib::phone_ow;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "uh") == true) {
    phoneNum = sclib::phone_uh;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "uw") == true) {
    phoneNum = sclib::phone_uw;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ux") == true) {
    phoneNum = sclib::phone_ux;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "er") == true) {
    phoneNum = sclib::phone_er;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ax") == true) {
    phoneNum = sclib::phone_ax;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ix") == true) {
    phoneNum = sclib::phone_ix;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "axr") == true) {
    phoneNum = sclib::phone_axr;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "ax-h") == true) {
    phoneNum = sclib::phone_ax_h;
    type = sclib::atSpeech|sclib::atVoiced|sclib::atPureSpeech;
  } else if (sclib::like(phone, "pau") == true) {
    phoneNum = sclib::phone_pau;
    type = sclib::atSpeech|sclib::atPause|sclib::atPureSpeech;
  } else if (sclib::like(phone, "epi") == true) {
    phoneNum = sclib::phone_epi;
    type = sclib::atSpeech|sclib::atPause|sclib::atPureSpeech;
  } else if (sclib::like(phone, "h#") == true) {
    phoneNum = sclib::phone_h_sharp;
    type = sclib::atSilence;
  }

  if (&phoneType != NULL) {
    phoneType = phoneNum;
  }

  return type;
}

//====================================================================================================================
// Returns for a given phone-label it's audioType; reports incosistencies in the matchings between phone-labels 
// and the corresponding string-prepresentations (found in phoneType2string() and phone2audioType(.) for both 
// directions
//====================================================================================================================
unsigned long int SC_GroundTruth_TIMIT::phone2audioType(unsigned long int phoneType) {
	unsigned long int audioType = sclib::noType, returnedPhoneType;
	const char* phoneString;
	
	if (phoneType != sclib::noPhone) {
		phoneString = phoneType2string(phoneType);
		audioType = phone2audioType((char*)(phoneString), returnedPhoneType);
		if (returnedPhoneType != phoneType) {
			REPORT_ERROR(SVLIB_BadData, "Matching between phone-type-labels and phone-string-representations is inconsistent!");
		}
	}

	return audioType;
}

//====================================================================================================================
// Return the string-representation of a phone to a given phone-label //Bing
//====================================================================================================================
const char* SC_GroundTruth_TIMIT::phoneType2string(int phoneType) {
   if (phoneType == sclib::phone_b) return "b";
	 else if (phoneType == sclib::phone_bcl) return "bcl";  
   else if (phoneType ==  sclib::phone_d) return "d";    
   else if (phoneType == sclib::phone_dcl) return "dcl";  
   else if (phoneType == sclib::phone_g) return "g";   
   else if (phoneType ==  sclib::phone_gcl) return "gcl";   
   else if (phoneType ==  sclib::phone_p) return "p";  
   else if (phoneType == sclib::phone_pcl) return "pcl";   
   else if (phoneType == sclib::phone_t) return "t";   
   else if (phoneType ==  sclib::phone_tck) return "tck";    
	 else if (phoneType == sclib::phone_k) return "k";   
	 else if (phoneType == sclib::phone_kcl) return "kcl";    
   else if (phoneType ==  sclib::phone_dx) return "dx";    
   else if (phoneType ==  sclib::phone_tcl) return "tcl";  
   else if (phoneType ==  sclib::phone_q) return "q";   
   else if (phoneType ==  sclib::phone_jh) return "jh";
   else if (phoneType == sclib::phone_ch) return "ch";   
   else if (phoneType == sclib::phone_s) return "s";   
   else if (phoneType == sclib::phone_sh) return "sh";   
   else if (phoneType == sclib::phone_z) return "z";    
   else if (phoneType ==  sclib::phone_zh) return "zh";   
   else if (phoneType ==  sclib::phone_f) return "f";    
   else if (phoneType == sclib::phone_th) return "th";   
   else if (phoneType == sclib::phone_v) return "v";    
   else if (phoneType == sclib::phone_dh) return "dh";   
   else if (phoneType == sclib::phone_m) return "m";   
   else if (phoneType ==  sclib::phone_n) return "n";    
   else if (phoneType == sclib::phone_ng) return"ng";    
   else if (phoneType ==  sclib::phone_em) return "em";   
   else if (phoneType ==  sclib::phone_en) return "en";    
	 else if (phoneType ==  sclib::phone_eng) return "eng";    
   else if (phoneType == sclib::phone_nx) return "nx";  
   else if (phoneType ==  sclib::phone_l) return "l";   
   else if (phoneType ==  sclib::phone_r) return "r";   
	 else if (phoneType ==  sclib::phone_w) return "w";   
	 else if (phoneType == sclib::phone_y) return "y";   
	 else if (phoneType ==  sclib::phone_hh) return "hh";   
	 else if (phoneType ==  sclib::phone_hv) return "hv";   
	 else if (phoneType ==  sclib::phone_el) return "el";    
	 else if (phoneType ==  sclib::phone_iy) return "iy";   
   else if (phoneType ==  sclib::phone_ih) return "ih";  
	 else if (phoneType == sclib::phone_eh) return "eh";
	 else if (phoneType == sclib::phone_ey) return "ey";   
	 else if (phoneType == sclib::phone_ae) return "ae";   
	 else if (phoneType ==  sclib::phone_aa) return "aa";    
	 else if (phoneType == sclib::phone_aw) return"aw";   
	 else if (phoneType ==  sclib::phone_ay) return "ay";   
	 else if (phoneType == sclib::phone_ah) return "ah";   
	 else if (phoneType == sclib::phone_ao) return "ao";    
	 else if (phoneType ==  sclib::phone_oy) return "oy";   
	 else if (phoneType == sclib::phone_ow) return "ow";    
	 else if (phoneType == sclib::phone_uh) return "uh";   
	 else if (phoneType == sclib::phone_uw) return "uw";  
   else if (phoneType ==  sclib::phone_ux) return "ux";   
   else if (phoneType == sclib::phone_er) return "er";    
	 else if (phoneType == sclib::phone_ax) return "ax";   
	 else if (phoneType ==  sclib::phone_ix) return "ix";   
   else if (phoneType == sclib::phone_axr) return "axr";   
   else if (phoneType == sclib::phone_ax_h) return "ax-h";    
   else if (phoneType == sclib::phone_pau) return "pau";  
   else if (phoneType ==  sclib::phone_epi) return "epi";    
	 else if (phoneType == sclib::phone_h_sharp) return "h#";
	 else return "";
}

//====================================================================================================================
//	insert a phone-label into the framelist
//====================================================================================================================
void SC_GroundTruth_TIMIT::setPhone(unsigned long int segmentStart, unsigned long int segmentEnd, int phoneType, int origin) {
  unsigned short int idxAT = (origin == sclib::modeHypothesized) ? 5 : 4;
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
		this->frameList[x][idxAT] = phoneType;
	}

	return;
}

//====================================================================================================================
// Returns a list with all phones occuring in the specified segment in order of appearence; one entry is per FLI
//====================================================================================================================
int* SC_GroundTruth_TIMIT::getPhones(unsigned long int segmentStart, unsigned long int segmentEnd, int &listSize, int origin){
  unsigned short int idxAT = (origin == sclib::modeHypothesized) ? 5 : 4;
  unsigned long int segStart = sample2FLI(segmentStart); //convert samples to FLIs, use only FLIs below
  unsigned long int segEnd = sample2FLI(segmentEnd);
	int *phoneList, y = 0;
	
	listSize = segEnd - segStart + 1;
	MArray_1D(phoneList, listSize, int, "SC_GroundTruth_TIMIT.getPhone: phoneList");

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
		phoneList[y++] = this->frameList[x][idxAT];
	}

	return phoneList;
}

//====================================================================================================================
// Adds a new column to the feature-set (new 0th column) containing phone-labels for each frame; if replicate==true,
// frames containing more than one phone (the majority) will be replicated to form one example for each of its phones;
// otherwise, each frame will be decorated with the phone-label that is valid for most of the time of the frame.
//====================================================================================================================
void SC_GroundTruth_TIMIT::addPhoneLabels(SV_Data* &pData, unsigned long int segmentStart, unsigned long int segmentEnd, int origin, bool replicate) {
	unsigned long int start, end;
	int *phoneList, listSize, d, l, t, t2 = 0, finalRowCount = 0, phone, maxCount;
	float **newMatrix, **tmp;
	std::map<int, int> distinctPhones;
	std::map<int, int>::iterator i;


	if (replicate == true) {
		//estimate how many phones are present in each frame => final row count, because each frame has to be replicated distinctPhoneTypes times
		for (t = 0; t < pData->Row; t++) {
			start = segmentStart + t*pData->Hdr.frameStep; //sample-numers of start and end of this frame in the whole multimedia object
			end = start + pData->Hdr.frameSize;

			phoneList = getPhones(start, end, listSize, origin); //get phones present in this frame
			for (l = 0; l < listSize; l++) {
				distinctPhones[phoneList[l]]++; //count which phone appears how ofton
			}
			finalRowCount += (int)(distinctPhones.size()); //number of distinct phones

			distinctPhones.clear(); //reset this hash
			MFree_1D(phoneList);
		}
	} else {
		finalRowCount = pData->Row;
	}

	//populate the new data matrix with added phone label column and replicated vectors for each phone occurence
	MArray_2D(newMatrix, finalRowCount, pData->Col+1, float, "SC_GroundTruth_TIMIT.addPhoneLabels: newMatrix");
	for (t = 0; t < pData->Row; t++) {
		start = segmentStart + t*pData->Hdr.frameStep;
		end = start + pData->Hdr.frameSize;

		phoneList = getPhones(start, end, listSize, origin);

		if (replicate == true) {
			for (l = 0; l < listSize; l++) {
				distinctPhones[phoneList[l]]++;
				if (distinctPhones[phoneList[l]] == 1) { //if this phone oocured for the first time in the current frame
					newMatrix[t2][0] = (float)(phoneList[l]); //add phone label as new column 0
					for (d = 0; d < pData->Col; d++) { //ad rest of feature values
						newMatrix[t2][d+1] = pData->Mat[t][d];
					}
					t2++; //increase row number after an added row
				}
			}
		} else {
			//find the phone occuring most ofton
			for (l = 0; l < listSize; l++) {
				distinctPhones[phoneList[l]]++;
			}
			maxCount = 0;
			for (i = distinctPhones.begin(); i != distinctPhones.end(); ++i) {
				if (i->second > maxCount) {
					phone = i->first;
					maxCount = i->second;
				}
			}

			newMatrix[t2][0] = (float)(phone); //add phone label as new column 0
			for (d = 0; d < pData->Col; d++) { //ad rest of feature values
				newMatrix[t2][d+1] = pData->Mat[t][d];
			}
			t2++; //increase row number after an added row
		}

		distinctPhones.clear();
		MFree_1D(phoneList);		
	}
	
	//change the pData object to have the new matrix with its parameters as its content
	pData->Row = finalRowCount;
	pData->Col++;
	tmp = pData->Mat;
	pData->Mat = newMatrix;
	if (pData->getJustLinked() == false) {
		MFree_2D(tmp);
	} else {
		pData->setJustLinked(true);
	}

	return;
}

//====================================================================================================================
// Add a file with it's sample-borders to the internal list
//====================================================================================================================
void SC_GroundTruth_TIMIT::addFile(const char* fileName, unsigned long int startSample, unsigned long int endSample) {
  SC_GroundTruth_TIMIT::SC_FileList *pItem = new SC_GroundTruth_TIMIT::SC_FileList(), *pHook;

  pItem->endSample = endSample;
  sprintf(pItem->fileName, "%s\0", fileName);
  pItem->Next = NULL;
  pItem->startSample = startSample;
  
  if (this->fileList == NULL) {
    this->fileList = pItem;
  } else {
    pHook = this->fileList;
    while (pHook->Next != NULL) {
      pHook = pHook->Next;
    }
    pHook->Next = pItem;
  }

  return;
}


//====================================================================================================================
// Get the filename out of the internal list holding the desired start-sample; the sample-nr out of the frameList 
// corresponding with the first sample in this file is also returned as offset; the sample-count in this file is 
// returned in sampleCount
//====================================================================================================================
char* SC_GroundTruth_TIMIT::whichFile(unsigned long int segmentStart, unsigned long int &offset, unsigned long int &sampleCount) {
  SC_GroundTruth_TIMIT::SC_FileList *pHook = this->fileList;
  char *fileName = NULL;

  offset = 0;
  sampleCount = 0;

  while (pHook != NULL) {
    if (segmentStart >= pHook->startSample && segmentStart <= pHook->endSample) { //start is within this file
      MArray_1D(fileName, strlen(pHook->fileName)+1, char, "SC_GroundTruth_TIMIT.whichFile: fileName");
      sprintf(fileName, "%s\0", pHook->fileName);
      offset = pHook->startSample;
      sampleCount = pHook->endSample - pHook->startSample + 1;
      break;
    }
    pHook = pHook->Next;
  }

  return fileName; 
}

//====================================================================================================================
// File-I/O for this class, so that the sate of a groundtruth-object can be saved to a file
// All members except the pSignalPrototype and the pTweak are saved
//====================================================================================================================
bool SC_GroundTruth_TIMIT::save(const char *fileName) {
	unsigned long int len, len2;
	long int x;
 	bool res = true;
	int bytes;
	fstream gtFile;
	SC_GroundTruth::SC_SpeakerMapping *pSpeakerHook = this->pSpeakerMapping;
	SC_GroundTruth_TIMIT::SC_FileList *pFileHook = this->fileList;
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
	bytes += io.writeMatrix(&gtFile, this->frameList, this->internalFrameCount, 6);
	bytes += io.writeScalar(&gtFile, this->verbose);
	bytes += io.writeScalar(&gtFile, this->setSceneBoundariesAtFileBoundaries);
	bytes += io.writeScalar(&gtFile, this->setSpeakerBoundariesAtFileBoundaries);
	bytes += io.writeScalar(&gtFile, this->sbSkipCount);
	bytes += io.writeScalar(&gtFile, this->setSceneBoundariesAtSpeakerBoundaries);

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
		bytes += io.writeScalar(&gtFile, pSpeakerHook->correct);
		bytes += io.writeScalar(&gtFile, pSpeakerHook->groundTruthID);
		bytes += io.writeScalar(&gtFile, pSpeakerHook->hypothesizedID);
		pSpeakerHook = pSpeakerHook->Next;
	}

	len = sclib::getListCount(this->fileList);
	len2 = sclib::bufferSize;
	bytes += io.writeScalar(&gtFile, len);
	for (x = 0; x < (long int)(len); x++) {
		bytes += io.writeScalar(&gtFile, len2);
		bytes += io.writeArray(&gtFile, pFileHook->fileName, sclib::bufferSize);
		bytes += io.writeScalar(&gtFile, pFileHook->startSample);
		bytes += io.writeScalar(&gtFile, pFileHook->endSample);
		pFileHook = pFileHook->Next;
	}

  if (gtFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Saving SC_GroundTruth_TIMIT Failed!");
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
SC_GroundTruth* SC_GroundTruth_TIMIT::load(const char *fileName) {
	unsigned long int len, len2, startSample, endSample;
	long int x, groundTruthID, hypothesizedID;
	bool correct;
	char *name = NULL;
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
  MArray_2D(this->frameList, (long int)(this->internalFrameCount), 6, long int, "SC_GroundTruth.load: frameList"); 
	bytes += io.readMatrix(&gtFile, this->frameList, this->internalFrameCount, 6, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->verbose, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->setSceneBoundariesAtFileBoundaries, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->setSpeakerBoundariesAtFileBoundaries, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->sbSkipCount, codeSizes, fileSizes);
	bytes += io.readScalar(&gtFile, this->setSceneBoundariesAtSpeakerBoundaries, codeSizes, fileSizes);

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

	bytes += io.readScalar(&gtFile, len, codeSizes, fileSizes);
	for (x = 0; x < (long int)(len); x++) {
		bytes += io.readScalar(&gtFile, len2, codeSizes, fileSizes);
		MArray_1D(name, len2+1, char, "SC_GroundTruth_TIMIT.load: name");
		bytes += io.readArray(&gtFile, name, len2, codeSizes, fileSizes);
		name[len2] = '\0';
		bytes += io.readScalar(&gtFile, startSample, codeSizes, fileSizes);
		bytes += io.readScalar(&gtFile, endSample, codeSizes, fileSizes);
		addFile(name, startSample, endSample);
		MFree_1D(name);
	}

	if (gtFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Loading SC_GroundTruth_TIMIT Failed!");
		pGT = NULL;
	}

	gtFile.close();

	return pGT;
}

//====================================================================================================================
// Many experiments on TIMIT assign each audio file of a speaker (there are 10 overall) to either of two utterances: 
// based on a lexicographic ordering of all 10 filenames per speaker, the first 8 are regarded as utterance 1, the
// remaining 2 are utterance 2. For the given audioFilenName and speakername, this method returns the utterance number
// returns sclib::noSpeaker on error
//====================================================================================================================
int SC_GroundTruth_TIMIT::getUtteranceNumber(const char *audioFileName, const char *speakerName) {
	int speakerID = this->getSpeakerIDfromName(speakerName, sclib::bufferSize);
	
	if (speakerID != sclib::noSpeaker) {
		//form a list of all files of the wanted speaker
		std::vector<std::string> audioFilesOfSpeaker;
		SC_GroundTruth_TIMIT::SC_FileList *pHook = this->fileList;
		while (pHook != NULL) {
			int currentSpeaker = this->getMajorSpeakerID(pHook->startSample, pHook->endSample, sclib::modeGroundtruth);
			if (currentSpeaker == speakerID) {
				audioFilesOfSpeaker.push_back(std::string(pHook->fileName));
			}
			pHook = pHook->Next;
		}
		
		//sanity check: if there are more than 10 files for this speaker, this method makes no sense because basic assumptions are not met
		if (audioFilesOfSpeaker.size() == 10) {
			//sort that list
			std::sort(audioFilesOfSpeaker.begin(), audioFilesOfSpeaker.end());

			//find the position of the given audioFile in that list
			std::string wantedFile = std::string(audioFileName);
			unsigned int position = 0;
			bool found = false;
			for(unsigned int i = 0; i < audioFilesOfSpeaker.size(); i++) {
				position++;
				if (audioFilesOfSpeaker[i].compare(wantedFile) == 0) {
					found = true;
					break;
				}
			}
			if (found) {
				return (position <= 8) ? 1 : 2;
			}
		} //10 files for this speaker
	}

	return sclib::noSpeaker;
}

//====================================================================================================================
// Extracts and returns the letter F or M in front of the TIMIT speaker name (last folder name)
//====================================================================================================================
char SC_GroundTruth_TIMIT::getGenderFromFilePath(const char* audioFileName) {
	//extract speaker-name and sex
	char *tmp = sclib::extractPath(audioFileName, false);
	char *lastSlash = strrchr(tmp, '/'); //separates the last directory name which contains information on the speaker's sex and name
	char sex = lastSlash[1]; //[0] is '/', the next character specifies males 'M' and females 'F'
	return sex;
}
