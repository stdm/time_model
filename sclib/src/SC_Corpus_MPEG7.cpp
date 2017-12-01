/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle Ralph's MPEG7 video-files,    */
/*        which have essentially no ground-truth available (rest like     */
/*        SCiVo)                                                          */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 11.05.2006																								*/
/**************************************************************************/

#include <iostream>
#include <list>
#include "SC_Corpus_MPEG7.h"
#include "SC_SignalHandler.h"
#include "SC_GroundTruth_MPEG7.h"
#include "SC_Aux.h"
#include "SC_FeatureHandler.h"
#include "SC_MatrixFunctions.h"
#include "SC_Signal_jWAVE.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Corpus_MPEG7::SC_Corpus_MPEG7(SC_TweakableParameters* pTweak, double videoFrameRate, const char *audioFileName, const char *sceneFile) : SC_Corpus(pTweak) {
  this->pGT = new SC_GroundTruth_MPEG7(this->pTweak, videoFrameRate, audioFileName, sceneFile);
}

#ifdef SC_USE_JNI
//====================================================================================================================
//	constructor opening stream from java
//====================================================================================================================
SC_Corpus_MPEG7::SC_Corpus_MPEG7(SC_TweakableParameters* pTweak, double videoFrameRate, JNIEnv *env, jobject jStreamObject) : SC_Corpus(pTweak) {
  this->pGT = new SC_GroundTruth_MPEG7(this->pTweak, videoFrameRate, env, jStreamObject);
}
#endif

//====================================================================================================================
//	constructor that receives an integer-array as a scene/cut-list  instead of a file containing the video-frame 
//  nubers where changes occured
//====================================================================================================================
SC_Corpus_MPEG7::SC_Corpus_MPEG7(SC_TweakableParameters* pTweak, double videoFrameRate, const char *audioFileName, int *sceneList, int sceneListLength) : SC_Corpus(pTweak) {
	this->pGT = new SC_GroundTruth_MPEG7(this->pTweak, videoFrameRate, audioFileName, sceneList, sceneListLength);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Corpus_MPEG7::~SC_Corpus_MPEG7() {

}

//====================================================================================================================
// provides the possibility to load the desired samples for conformity reasons (it is a straight-forward call of
// the corresponding method in SC_SegmentationHandler)
// the given borders reamin unchanged in this method
//====================================================================================================================
SC_Signal* SC_Corpus_MPEG7::loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries) {
  SC_Signal *pSignal = NULL;
	short *samples = NULL, *buf;
	unsigned long int len;

	if (this->pGT->getSignalPrototype()->getSignalType() != sclib::stJWave) {	
		pSignal = this->pSignalHandler->loadSignal(this->pGT->getAudioFileName(), segmentStart, sclib::min(segmentEnd, this->pGT->getAudioSampleCount()-1));
	} else { //in the case of a jStream source, all the samples are already cached in the groundtruth object and just need to be copied out
		pSignal = new SC_Signal_jWAVE((SC_Signal_jWAVE*)(this->pGT->getSignalPrototype()));	
		len = sclib::min(segmentEnd, this->getGT()->getAudioSampleCount()-1) - segmentStart + 1;
		MArray_1D(samples, len, short, "SC_Corpus_MPEG7.loadSignal: samples");
		buf = this->pGT->getSignalPrototype()->GetBuf_L();
		for (unsigned long int t = 0; t < len; t++) {
			samples[t] = buf[segmentStart+t];
		}
		pSignal->setBuf_L(samples, len);
	}

  return pSignal;
}

//====================================================================================================================
// write the given feature-set to a file formated to fit Ralph's needs (videoFrameStart, videoFrameEnd, features, ...)
//====================================================================================================================
void SC_Corpus_MPEG7::featureOut(const char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pData, bool printHeader, bool printVideoFrames, const char *featureName, bool soleHeader) {
	fstream fileOut;
	char *fName, *featName = const_cast<char*>(featureName);
	unsigned long int videoFrameStart, videoFrameEnd, start, end;
	double srRatio = (double)(this->pGT->getAudioSampleRate()) / (double)(pData->Hdr.sampleRate);
	SC_FeatureHandler handler(this->pTweak, false);

  fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
  sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, ios_base::out|ios_base::app);

  if (printHeader == true) {
		if (featName == NULL) {
			featName = handler.getFeatureName(pData->Hdr.ID);
		}
		fileOut << "Feature type: " << featName << endl;
		if (featName != featureName) {
			MFree_1D(featName);
		}
		if (soleHeader == false) {
			fileOut << "Rows: " << pData->Row << endl;
		}
    fileOut << "Columns: " << pData->Col << endl;
    fileOut << "Frame-size [samples]: " << pData->Hdr.frameSize << endl;
    fileOut << "Frame-step [samples]: " << pData->Hdr.frameStep << endl;
		if (soleHeader == false) {
			fileOut << "Start [samples] (videoFrame): " << segmentStart << " (" << this->pGT->getConverter()->sample2videoFrame(segmentStart) << ")" << endl;
			fileOut << "End [samples] (videoFrame): " << segmentEnd << " (" << this->pGT->getConverter()->sample2videoFrame(segmentEnd) << ")" << endl;
		}
		fileOut << "Sample rate [Hz]: " << pData->Hdr.sampleRate << endl;
    fileOut << endl;
		if (printVideoFrames == true) {
			fileOut << "VF start\tVF end\tFeature values" << endl;
		} else {
			fileOut << "Feature values" << endl;
		}
  }
 
  for (int y = 0; y < pData->Row; y++) {
		if (printVideoFrames == true) {
			start = (unsigned long int)(srRatio * y*pData->Hdr.frameStep);
			end = (unsigned long int)(srRatio * (y*pData->Hdr.frameStep + pData->Hdr.frameSize)) - 1;
			videoFrameStart = this->pGT->getConverter()->sample2videoFrame(segmentStart + start); //first sample of current audioFrame -> videoFrame
			videoFrameEnd = this->pGT->getConverter()->sample2videoFrame(segmentStart + end); //last sample of current audioFrame -> videoFrame
			fileOut << videoFrameStart << "\t" << videoFrameEnd << "\t";
		}

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
// another method to output features for input on java side: this time a CSV file is written, containing in the first
// row the feature-name of each column (as a column name for db-tables), and then per row the features for one video-
// frame, regardless of what was the audioframe-size (mean of all audio-frames beginning in the video-frame is 
// written). the first column contains videoframe-numbers; pData is meant to be an array of datasets each containing 
// a different feature (as returned by SC_FeatureHandler.extractFeatures())
//====================================================================================================================
void SC_Corpus_MPEG7::lowLevelFeaturesOut(const char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pData, unsigned long int selectionMask) {
	int x, y, z, rows = 0;
	char *fName, *tmp, featureName[sclib::bufferSize];
	unsigned long int audioFrameStart, audioFrameEnd, videoFrame, videoFrameStart, videoFrameEnd;
	double *featureVec = NULL, srRatio = 0.0;
	const char *separator = "\t";
	fstream fileOut;
	SC_FeatureHandler *pFeatureHandler = new SC_FeatureHandler(this->pTweak);
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();

  fName = new char[strlen(this->pTweak->debug.debugDir) + strlen(fileName) + 1];
  sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
	fileOut.open(fName, ios_base::out|ios_base::app);
	MFree_1D(fName);

	//write header (feature-names or replacement if can't be deduced)
  fileOut << setiosflags(ios_base::left|ios_base::fixed|ios_base::showpoint);
	fileOut << this->pGT->getAudioFileName() << endl;
	fileOut <<  "FrameRate" << separator << ((SC_GroundTruth_MPEG7*)this->pGT)->getVideoFrameRate() << endl;

	fileOut << "videoFrame"; //this is the first column to be written
	for (z = 0; z < (int)(pFeatureHandler->getFeatureCount()); z++) {
		if (pData[z] != NULL && sclib::bitTest(selectionMask, sclib::bit(z)) == true) {
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

	//write mean of audio-frames per video-frame in the given segment
	videoFrameStart = this->pGT->getConverter()->sample2videoFrame(segmentStart, sclib::alignmentStart);
	videoFrameEnd = this->pGT->getConverter()->sample2videoFrame(segmentEnd, sclib::alignmentEnd);
	for (videoFrame = videoFrameStart; videoFrame <= videoFrameEnd; videoFrame++) {
		fileOut << videoFrame;

		for (z = 0; z < (int)(pFeatureHandler->getFeatureCount()); z++) {
			if (pData[z] != NULL && sclib::bitTest(selectionMask, sclib::bit(z)) == true) {
				srRatio = (double)(this->pGT->getAudioSampleRate()) / (double)(pData[z]->Hdr.sampleRate);
				audioFrameStart = this->pGT->getConverter()->videoFrame2sample(videoFrame, sclib::alignmentStart); //first sample of the videoFrame
				audioFrameStart = (unsigned long int)((audioFrameStart - segmentStart > (srRatio*pData[z]->Hdr.frameSize)) ? audioFrameStart - segmentStart - (srRatio*pData[z]->Hdr.frameSize) : 0); //startsample of a possible first frame with segmentStart as reference (not video-start); only the end of this frame has to be in the videoFrame, not necessarily the beginning
				audioFrameStart = (unsigned long int)ceil((double)(audioFrameStart) / (double)(srRatio*pData[z]->Hdr.frameStep)); //index into pData[z]->Mat[]
				audioFrameEnd = this->pGT->getConverter()->videoFrame2sample(videoFrame, sclib::alignmentEnd); //last sample of the videoFrame
				audioFrameEnd = (unsigned long int)floor((double)(audioFrameEnd - segmentStart) / (double)(srRatio*pData[z]->Hdr.frameStep)); //index into pData[z]->Mat[]

				featureVec = pFunc->mean(pData[z]->Mat, pData[z]->Row, pData[z]->Col, audioFrameStart, sclib::min(audioFrameEnd+1, pData[z]->Row));
				for (x = 0; x < pData[z]->Col; x++) {
					fileOut << separator << featureVec[x];
				}
				MFree_1D(featureVec);
			}
		}		

		fileOut << endl;
	}
	
  fileOut.close();
	MFree_0D(pFeatureHandler);
	MFree_0D(pFunc);

	return;
}

//====================================================================================================================
// output all highlevel-features (=algorithmic results) for the whole video at once alike to the low-level features
//====================================================================================================================
void SC_Corpus_MPEG7::highLevelFeaturesOut(const char *fileName, long int selectionMask, bool selectSpkID) {
	((SC_GroundTruth_MPEG7*)this->pGT)->analyzedOutEx(fileName, ios_base::out|ios_base::app, selectionMask, selectSpkID);
	return;
}

//====================================================================================================================
// parses the given XML file and returns an array of video-frame-based shot-boundaries along with the size of the 
// array
//====================================================================================================================
int* SC_Corpus_MPEG7::getShotList(char *xmlFileName, double videoFrameRate, int &listLength) {
	int startPos, i, *shotList = NULL, h = 0, m = 0, s = 0, frac = 0, total = 0, videoFrame;
	std::list<int> shotContainer;
	FILE* inFile = fopen(xmlFileName, "r");
	char buffer[sclib::bufferSize];
	bool inShot = false;
	SC_Conversion converter(16000, 16000.0/videoFrameRate); //we can assume andy samplerate here

	if (inFile != NULL) {
		while (!feof(inFile)) {
			sclib::readline(inFile, buffer, sclib::bufferSize);
			if (sclib::like(buffer, "%<VideoSegment id=\"$shot\">%") == true) {
				inShot = true; //start of a shot boundary segment
			} else if (sclib::like(buffer, "%</VideoSegment>%") == true) {
				inShot = false; //end of this segment
			} else if (inShot == true) {
				if (sclib::like(buffer, "%<MediaTimePoint>%") == true) { //the actual time of the boundary is in this line
					//a line has the form "                      <MediaTimePoint>T00:00:13:2F25</MediaTimePoint>"
					//how to interpret the timestamp (from http://www-nlpir.nist.gov/projects/tv2003/common.shot.ref/time.elements):
					//mediaTimePointType: a datatype specifying a time stamp of the media using Gregorian date and day time 
					//without specifying the TZD (see also TimePoint = datatype).
          //
          //YYYY-MM-DDThh:mm:ss:nnnFNNN 
          //
          //The following lexicals are used for digits of the corresponding date/time elements: 
          //
          //Y: Year, can be a variable number of digits, 
          //M: Month, 
          //D:Day, 
          //h: hour, 
          //m: minute, 
          //s: second, 
          //n: number of fractions, nnn can be any number between 0 and NNN-1 (NNN and with it nnn can have an 
					//   arbitrary number of digits).
          //N: number of fractions of one second which are counted by nnn. NNN can have a arbitrary number of digits and is not limited to =
          //   three.
          //Also delimiters for the time specification (T) and the number of fractions of one second are used (F).
					startPos = 0;
					startPos = sclib::getNextIntFromString(buffer, sclib::bufferSize, h, startPos, ":F<");
					startPos = sclib::getNextIntFromString(buffer, sclib::bufferSize, m, startPos, ":F<");
					startPos = sclib::getNextIntFromString(buffer, sclib::bufferSize, s, startPos, ":F<");
					startPos = sclib::getNextIntFromString(buffer, sclib::bufferSize, frac, startPos, ":F<");
					startPos = sclib::getNextIntFromString(buffer, sclib::bufferSize, total, startPos, ":F<");
					videoFrame = converter.ms2videoFrame(sclib::round((double)(frac*1000)/(double)(total)) + s*1000 + m*60000 + h*3600000);
					shotContainer.push_back(videoFrame);
				}
			}
		}

		fclose(inFile);
	}

	//copy content of linked-list container to the desired array
	if (shotContainer.size() > 0) {
		listLength = (int)(shotContainer.size());
		MArray_1D(shotList, listLength, int, "SC_Corpus_MPEG7.getShotList: shotList");
		for (i = 0; i < listLength; i++) {
			shotList[i] = shotContainer.front();
			shotContainer.pop_front();
		}
	} else {
		listLength = 0;
	}

	return shotList;	
}
