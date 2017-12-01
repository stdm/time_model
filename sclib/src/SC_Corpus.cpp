/**************************************************************************/
/*    Responsibility:																											*/
/*      - base class for objects to encapsulates all algorithms and       */
/*        methods to handle a specific corpus of data and it's specific   */
/*        tasks. because the corpora are highly different and common      */
/*        taks and methods are seldom, there are only few methods in the  */
/*        base class                                                      */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.04.2006																								*/
/**************************************************************************/

#include "SC_Corpus.h"
#include "SC_Aux.h"
#include "SC_FeatureHandler.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Corpus::SC_Corpus(SC_TweakableParameters* pTweak) {
  this->pTweak = pTweak;
  this->pGT = NULL;
	this->pSignalHandler = new SC_SignalHandler(this->pTweak);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Corpus::~SC_Corpus() {
  MFree_0D(this->pGT);
	MFree_0D(this->pSignalHandler);
}

//====================================================================================================================
// write the given feature-set to a file formated to fit Ralph's needs (videoFrameStart, videoFrameEnd, features, ...)
//====================================================================================================================
void SC_Corpus::featureOut(char *fileName, char *featureName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pData, bool printHeader) {
	fstream fileOut;
	char *fName;
	unsigned long int start, end;
	double srRatio = (double)(this->pGT->getAudioSampleRate()) / (double)(pData->Hdr.sampleRate);

  fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
  sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, ios_base::out|ios_base::app);

  if (printHeader == true) {
		fileOut << "Feature type: " << featureName << endl;
    fileOut << "Rows: " << pData->Row << endl;
    fileOut << "Columns: " << pData->Col << endl;
    fileOut << "Frame-size [samples]: " << pData->Hdr.frameSize << endl;
    fileOut << "Frame-step [samples]: " << pData->Hdr.frameStep << endl;
    fileOut << "Start [samples]: " << segmentStart << endl;
    fileOut << "End [samples]: " << segmentEnd << endl;
		fileOut << "(Sample rate: " << pData->Hdr.sampleRate << "/s)" << endl;
    fileOut << endl;
		fileOut << "Start sample\tEnd sample\tFeature values" << endl;
  }
 
  for (int y = 0; y < pData->Row; y++) {
		start = (unsigned long int)(srRatio * y*pData->Hdr.frameStep);
		end = (unsigned long int)(srRatio * (y*pData->Hdr.frameStep + pData->Hdr.frameSize)) - 1;
    fileOut << segmentStart+start << "\t" << segmentStart+end << "\t";

    for (int x = 0; x < pData->Col; x++) {
      fileOut << pData->Mat[y][x];
      if (x < pData->Col-1) {
        fileOut << "\t";
      }
    }
    fileOut << endl;
  }

  fileOut.close();
	MFree_1D(fName);

	return;
}

//====================================================================================================================
// outputs the wanted low-level-features of all voiced speech frames in the given segment, preceded by each frame's 
// gt-speaker-id
//====================================================================================================================
void SC_Corpus::speakerFeaturesOut(char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pData, unsigned long int selectionMask) {
	int x, y, z, rows = 0;
	char *fName, *tmp, featureName[sclib::bufferSize];
	unsigned long int start, end, frameSize, frameStep;
	const char *separator = "\t";
	fstream fileOut;
	SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);

  fName = new char[strlen(this->pTweak->debug.debugDir) + strlen(fileName) + 1];
  sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, ios_base::out|ios_base::app);
	MFree_1D(fName);

	//write header (feature-names or replacement if can't be deduced)
  fileOut << setiosflags(ios_base::left|ios_base::fixed|ios_base::showpoint);
	fileOut << this->pGT->getAudioFileName() << endl;
	fileOut << "SpeakerID"; //title of first column
	for (z = 0; z < (int)(pFeatureHandler->getFeatureCount()); z++) {
		if (pData[z] != NULL && sclib::bitTest(selectionMask, sclib::bit(z)) == true) {
			if (rows == 0) { //check for equal rowcount
				rows = pData[z]->Row;
				frameSize = pData[z]->Hdr.frameSize;
				frameStep = pData[z]->Hdr.frameStep;
			} else if (rows != pData[z]->Row) {
				REPORT_ERROR(SVLIB_BadArg, "Datasets need to have equal rowcount for outputting");
			}

			tmp = pFeatureHandler->getFeatureName(pData[z]->Hdr.ID);
			for (y = 0; y < pData[z]->Col; y++) {
				if (strncmp(tmp, "", sclib::bufferSize) == 0) {
					sprintf(featureName, "%d_%d", z, y); 
				} else if (pData[z]->Col == 1) {
					sprintf(featureName, "%s", tmp); 
				} else {
					sprintf(featureName, "%s_%d", tmp, y); 
				}
				fileOut << separator << featureName;
			}
			MFree_1D(tmp);
		}
	}
	fileOut << endl;

	//write speaker-id and audio-frames of voiced speech frames in the given segment
	for (y = 0; y < rows; y++) {
		start = segmentStart + (y*frameStep);
		end = start + frameSize;
		if (this->pGT->testSegment(start, end, true, sclib::atSpeech|sclib::atVoiced, true, sclib::atPause, false, sclib::modeHypothesized) > 0) {
			fileOut << this->pGT->getSamplesSpeakerID(start, sclib::modeGroundtruth); //first column: speaker-id
			for (z = 0; z < (int)(pFeatureHandler->getFeatureCount()); z++) {
				if (pData[z] != NULL && sclib::bitTest(selectionMask, sclib::bit(z)) == true) {
					for (x = 0; x < pData[z]->Col; x++) {
						fileOut << separator << pData[z]->Mat[y][x];
					}
				}
			}	
			fileOut << endl;
		}
	}
	
  fileOut.close();
	MFree_0D(pFeatureHandler);

	return;
}

//====================================================================================================================
// Return the fileName of the audio file containg the given startSample
//====================================================================================================================
char* SC_Corpus::getCurrentAudioFileName(unsigned long int startSample) {
	char *fileName;

	MArray_1D(fileName, sclib::bufferSize, char, "SC_Corpus.getCurentAudioFileName: fileName");
	sprintf(fileName, "%s", this->getGT()->getAudioFileName());

	return fileName;
}
