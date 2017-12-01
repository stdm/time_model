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

#include <math.h>
#include <string.h>
#include <iomanip>
#include "SC_GroundTruth_MAC.h"
#include "SC_SignalHandler.h"
#include <SV_DataIO.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth_MAC::SC_GroundTruth_MAC(SC_TweakableParameters *pTweak, const char* corpusFileName, bool verbose) : SC_GroundTruth(pTweak, NULL) {
  //some init-stuff
  this->gtType = sclib::gtMAC;
  this->audioFileName = new char[strlen(corpusFileName)+1]; //yes, a missuse of this member, but seems ok ;-
	strcpy(this->audioFileName, corpusFileName);
  this->fileList = NULL;
  this->verbose = verbose;
	this->uncertaintyRegion = 0;

  if (measureCorpusSize(corpusFileName) == true) {
    initFrameList();

    //read the ground-truth files
    readGroundTruth();
  }
}

//====================================================================================================================
//	destructor 
//====================================================================================================================
SC_GroundTruth_MAC::~SC_GroundTruth_MAC() {
  sclib::destructLinkedList(this->fileList);
}

//====================================================================================================================
// This class holds the groundtruth for a whole corpus as if it where a aignle file; this function reads the 
// parameters of all files in the corpus to measure it's overall size and create a mapping between the frames in 
// the framelist and the corresponding files storing the actual signals
//====================================================================================================================
bool SC_GroundTruth_MAC::measureCorpusSize(const char *corpusFileName) {
  unsigned long int overallSampleCount = 0, lineCount;
  bool res = true;
  int count = 0;
  char waveformFileName[sclib::bufferSize], errStr[sclib::bufferSize], *completeFileName = NULL, *baseDir = NULL;
  FILE *corpusFile;
  SC_SignalHandler signalHandler(this->pTweak, sclib::stWave);
  SC_Signal *pSignal = NULL;
	double lastPercentage = 0.0;

  if (this->verbose == true) {
    printf("Measuring MAC corpus size for corpus file %s: ", corpusFileName);
  }

  lineCount = sclib::countLines(corpusFileName);
  baseDir = sclib::extractPath(corpusFileName);

  if ((corpusFile = fopen(corpusFileName, "r")) != NULL) {
    while (fscanf(corpusFile, "%s", waveformFileName) != EOF) {
      //construct complete waveformFileName
      MArray_1D(completeFileName, strlen(baseDir) + strlen(waveformFileName) + 1, char, "SC_GroundTruth_MAC.measureCorpusSize: completeFileName");
      sprintf(completeFileName, "%s%s\0", baseDir, waveformFileName);

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
          if (this->pSignalPrototype->SigPar.Encode != pSignal->SigPar.Encode || this->pSignalPrototype->SigPar.NChannel != pSignal->SigPar.NChannel || this->pSignalPrototype->SigPar.SRate != pSignal->SigPar.SRate || this->pSignalPrototype->SigPar.StByte != pSignal->SigPar.StByte) {
            sprintf(errStr, "Corpus-File %s has different signal parameters than the rest", completeFileName);
            REPORT_ERROR(SVLIB_FileErr, errStr);
            res = false;
          }
          MFree_0D(pSignal);
        }
      }

      if (this->verbose == true) {
				lastPercentage = sclib::printPercentage(lineCount, count, lastPercentage, 5.0, !count);
      }

      MFree_1D(completeFileName);
      count++;
    }

    fclose(corpusFile);

    this->audioSampleRate = this->pSignalPrototype->SigPar.SRate;
		this->pConverter->setAudioSampleRate(this->audioSampleRate);
    this->audioSampleCount = overallSampleCount;
	  this->internalFrameSize = this->pConverter->ms2sample(this->pTweak->groundTruth.internalFrameSize); //in samples!!!
    this->internalFrameCount = (unsigned long int)(ceil((double)(this->audioSampleCount) / (double)(this->internalFrameSize)));

  } else {
    res = false;
  }

  MFree_1D(baseDir);

  if (this->verbose == true) {
		sclib::printPercentage(1, 1, 0.0, 0.0, false);
    printf("done! Size: %0.1f min.\n",  this->pConverter->sample2ms(this->audioSampleCount)/60000.0);
  }

  return res;
}

//====================================================================================================================
// Read Groundtruth from file(s) and initialize internal data structures (frameList etc.)
// This is rather siple here because all relevant information has already been extracted in the measureCorpusSize()
// function
//====================================================================================================================
bool SC_GroundTruth_MAC::readGroundTruth(void) {
  bool res = true;
  SC_GroundTruth_MAC::SC_FileList *pFileListHook = this->fileList;

  if (this->verbose == true) {
    printf("Loading ground truth for corpus file %s: ", this->audioFileName);
  }

  if (this->fileList != NULL) {
    pFileListHook = this->fileList;
    while (pFileListHook != NULL) {
			setSegment(pFileListHook->startSample, pFileListHook->startSample, sclib::atSceneBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth, false, false);
			setSegment(pFileListHook->startSample, pFileListHook->startSample, sclib::atSceneBoundary, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized, false, false);
      setSegment(pFileListHook->startSample, pFileListHook->endSample, pFileListHook->types, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeGroundtruth);
      setSegment(pFileListHook->startSample, pFileListHook->endSample, pFileListHook->types, false, sclib::noSpeaker, sclib::modeLabelAdd, sclib::modeHypothesized);
      pFileListHook = pFileListHook->Next;
    }
  } else {
    res = false;
  }

  if (this->verbose == true) {
    printf("done!\n");
  }

  return res;
}

//====================================================================================================================
// add a file with it's sample-borders to the internal list
//====================================================================================================================
void SC_GroundTruth_MAC::addFile(const char* fileName, unsigned long int startSample, unsigned long int endSample) {
  long int types;
  SC_GroundTruth_MAC::SC_FileList *pItem = new SC_GroundTruth_MAC::SC_FileList(), *pHook;

  //extract audio-type from path
  if (strstr(fileName, "/ACTION/") != NULL) {
    types = sclib::atAction|sclib::atNoise;
  } else if (strstr(fileName, "/BACKGROUND/") != NULL) {
    types = sclib::atBackground|sclib::atNoise;
  } else if (strstr(fileName, "/MUSIC/") != NULL) {
    types = sclib::atMusic|sclib::atNoise;
  } else if (strstr(fileName, "/BREATH/") != NULL) {
    types = sclib::atBreath|sclib::atNoise;
  } else if (strstr(fileName, "/NOISYSPEECH/") != NULL) {
    types = sclib::atNoisySpeech|sclib::atSpeech;
  } else if (strstr(fileName, "/PURESPEECH/") != NULL) {
    types = sclib::atPureSpeech|sclib::atSpeech;
  } else {
		if (this->verbose == true) {
			printf(" file %s has undefined audio-type", fileName);
		}
    types = sclib::atUndefined;
  }

  pItem->endSample = endSample;
  sprintf(pItem->fileName, "%s\0", fileName);
  pItem->Next = NULL;
  pItem->startSample = startSample;
  pItem->types = types;

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
// get the filename out of the internal list holding the desired start-sample; the sample-nr out of the frameList 
// corresponding with the first sample in this file is also returned as offset; the sample-count in this file is 
// returned in sampleCount
//====================================================================================================================
char* SC_GroundTruth_MAC::whichFile(unsigned long int segmentStart, unsigned long int &offset, unsigned long int &sampleCount, long int &types) {
  SC_GroundTruth_MAC::SC_FileList *pHook = this->fileList;
  char *fileName = NULL;

  offset = 0;
  sampleCount = 0;
  types = sclib::noType;

  while (pHook != NULL) {
    if (segmentStart >= pHook->startSample && segmentStart <= pHook->endSample) { //start is within this file
      MArray_1D(fileName, strlen(pHook->fileName)+1, char, "SC_GroundTruth_MAC.whichFile: fileName");
      sprintf(fileName, "%s\0", pHook->fileName);
      offset = pHook->startSample;
      sampleCount = pHook->endSample - pHook->startSample + 1;
      types = pHook->types;
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
bool SC_GroundTruth_MAC::save(const char *fileName) {
	unsigned long int len, len2;
	long int x;
 	bool res = true;
	int bytes;
	fstream gtFile;
	SC_GroundTruth::SC_SpeakerMapping *pSpeakerHook = this->pSpeakerMapping;
	SC_GroundTruth_MAC::SC_FileList *pFileHook = this->fileList;
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
	bytes += io.writeScalar(&gtFile, this->verbose);

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
		REPORT_ERROR(SVLIB_Fail, "Saving SC_GroundTruth_MAC Failed!");
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
SC_GroundTruth* SC_GroundTruth_MAC::load(const char *fileName) {
	unsigned long int len, len2, startSample, endSample;
	long int x, groundTruthID, hypothesizedID;
	bool correct;
	char *name = NULL;
	int bytes;
	fstream gtFile;
	SC_GroundTruth_MAC *pGT = this;
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
	bytes += io.readScalar(&gtFile, this->verbose, codeSizes, fileSizes);

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
		MArray_1D(name, len2+1, char, "SC_GroundTruth_MAC.load: name");
		bytes += io.readArray(&gtFile, name, len2, codeSizes, fileSizes);
		name[len2] = '\0';
		bytes += io.readScalar(&gtFile, startSample, codeSizes, fileSizes);
		bytes += io.readScalar(&gtFile, endSample, codeSizes, fileSizes);
		addFile(name, startSample, endSample);
		MFree_1D(name);
	}

	if (gtFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Loading SC_GroundTruth_MAC Failed!");
		pGT = NULL;
	}

	gtFile.close();

	return pGT;
}
