/**************************************************************************/
/*    Responsibility: 																										*/
/*      - return an SC_Signal object able to read/write the specified     */
/*        audio-type                                                      */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#include "SC_Aux.h"
#include "SC_SignalHandler.h"
#include "SC_Signal_WAVE.h"
#include "SC_Signal_NIST.h"
#include "SC_Signal_jWAVE.h"
#include "SC_Resampling.h"
#include <SV_Error.h>
#ifdef SC_USE_JNI
	#include <jni.h>
#endif

//====================================================================================================================
// default constructor
//====================================================================================================================
SC_SignalHandler::SC_SignalHandler(SC_TweakableParameters *pTweak, unsigned short int signalType) {
  this->pTweak = pTweak;
  this->signalType = signalType;
	this->pMPEG = NULL;
}

//====================================================================================================================
// destructor 
//====================================================================================================================
SC_SignalHandler::~SC_SignalHandler() {
	MFree_0D(pMPEG);
}

//====================================================================================================================
//	guess the signal-type from the filename's extension or the files content,if the extension is ambigous
//====================================================================================================================
unsigned short int SC_SignalHandler::guessSignalType(const char *fileName) {
  unsigned short int sType = sclib::stOther; //make "other" the default-type, if it can not be estimated from the filename (because the underlying ffmpeg lib can autodetect & load most things...)
	char *ext = sclib::uCase(sclib::extractExtension(fileName)), buffer[5];
	fstream signalFile;

  if (ext != NULL) {
    if (strncmp(ext, "WAV", 3) == 0) {
			//open signal file and read first 4 byte: for MS-RIFF-WAV files they should be 'RIFF', for NIST SPHERE files it should be 'NIST'
			if (sclib::fileExists(fileName) == true) {
				signalFile.open(fileName, ios::in|ios::binary);
					signalFile.read((char*)(buffer), 4 * sizeof(char));
				signalFile.close();
				buffer[4] = '\0';
				if (strncmp(buffer, "RIFF", 4) == 0) {
					sType = sclib::stWave;
				} else if (strncmp(buffer, "NIST", 4) == 0) {
					sType = sclib::stNIST;
				}
			}
    } else if (strncmp(ext, "MP2", 3) == 0) {
      sType = sclib::stMP2;
    } else if (strncmp(ext, "MPG", 3) == 0 || strncmp(ext, "MPEG", 4) == 0) {
      sType = sclib::stMPEG;
    } else if (strncmp(ext, "MP3", 3) == 0) {
      sType = sclib::stMP3;
    } else if (strncmp(ext, "JSTREAM", 7) == 0) {
			sType = sclib::stJWave;
		}

    MFree_1D(ext);
	} else {
		char buf[sclib::bufferSize];
		sprintf(buf, "%s", fileName);
		sclib::uCase(buf);
		if (strncmp(fileName, "JSTREAM", 7) == 0) {
			sType = sclib::stJWave;
		}
	}

  return sType;
}

//====================================================================================================================
//	just opens the signal so that it's parameters can be extracted; if forceSampleRate>0, the signal (here, at least, 
//  the parameters in the header) will be resampled to the given samplerate regardless of it's original sampleRate
//====================================================================================================================
SC_Signal* SC_SignalHandler::openSignal(const char* fileName, int forceSampleRate, void* jniEnvironment, void* jStreamObject) {
  unsigned short int sType = (this->signalType == sclib::stGuess) ? guessSignalType(fileName) : this->signalType;
  char *fName, *prefix;
  SC_Signal *pSignal = NULL;
	bool enough = false;
  
  //install, if wished, a prefix before the debugDir
  fName = sclib::extractFileName(fileName);
  prefix = sclib::exchangeFileExtension(fName, "");
  MFree_1D(fName);
  this->pTweak->setDebugPrefix(prefix);
  MFree_1D(prefix);

	//load signal-parameters; try to load it with a more general class if a first attempt fails (max. 2 loops here...)
	while (enough == false) {
		switch (sType) {
			case sclib::stWave: {
				pSignal = new SC_Signal_WAVE(fileName);
				if (pSignal->isOk() == false) {
					sType = sclib::stOther; //use more general class on error
				} else {
					enough = true;
				}
				break;
			}
			case sclib::stMP2: case sclib::stMPEG: case sclib::stMP3: case sclib::stOther: {
				pSignal = new SC_Signal_MPEG(fileName, sType, true, (forceSampleRate > 0) ? forceSampleRate : this->pTweak->signalMpeg.outputSampleRate, this->pTweak->signalMpeg.outputChannelCount, this->pTweak->signalMpeg.hqResampling, this->pTweak->signalMpeg.fastSeeking);
				enough = true; //we have no more options on failure if we already used the SC_Signal_MPEG class
				break;
			}
			case sclib::stNIST: {
				pSignal = new SC_Signal_NIST(fileName);
				if (pSignal->isOk() == false) {
					sType = sclib::stOther; //use more general class on error
					MFree_0D(pSignal);
				} else {
					enough = true;
				}
				break;
			}
			case sclib::stJWave: {
#ifdef SC_USE_JNI
				pSignal = new SC_Signal_jWAVE((JNIEnv*)(jniEnvironment), (jobject)(jStreamObject));
#endif
				enough = true;
				break;
			}
		}
	}

	if (pSignal == NULL || pSignal->isOk() == false) {
		REPORT_ERROR(SVLIB_BadData, "Opening the signal failed ultimately");
		MFree_0D(pSignal);
	} else {
		if (forceSampleRate > 0 && forceSampleRate != pSignal->SigPar.SRate) { //"resample" parameters if wished
			pSignal->setSampleCount(sclib::round(pSignal->getSampleCount() * (forceSampleRate / (double)(pSignal->SigPar.SRate))));
			pSignal->SigPar.SRate = forceSampleRate;
		}
	}

  return pSignal;
}

//====================================================================================================================
// Meta-loader of SC_Signal_MPEG based signals: Instead of just opening the signal with a new (fresh-state) object,
// this method uses a permanent SC_Signal_MPEG object that loads all content and then links it to the given pSignal or
// sets it NULL in case of error; this way, the seeking problem of ffmpeg can be "workarounded" for the prominent case
// that the next consecutive call of the load-method just wants the segment that directly succeeds the previously
// opened one
//====================================================================================================================
void SC_SignalHandler::loadMPEG(unsigned long int start, unsigned long int end, SC_Signal_MPEG* &pSignal) {
	char fileName[sclib::bufferSize];
	
	if (pSignal != NULL) { //it is mandatory that the signal could be opened, so it has correct format and the file could be read
		//(re-)create the permanent copy if it doesn't exist or the source to load changed
		if (this->pMPEG == NULL || strcmp(this->pMPEG->getFileName(), pSignal->getFileName()) != 0) {
			MFree_0D(this->pMPEG);
			this->pMPEG = new SC_Signal_MPEG(pSignal);
			if (this->pMPEG->readHeader(pSignal->getFileName()) == false) {
				MFree_0D(this->pMPEG);
				MFree_0D(pSignal);
				return;
			}
		}
		
		//load the samples in the permanent copy (this can make use of object-internal caching if a segment 
		//immediately following the last loaded segment is going to be loaded) and then 
		//copy it to the new signal, forgetting the buffer in the permanent one
		sprintf(fileName, "%s\0", pMPEG->getFileName()); //need to give a copy, because during readHeader it is freed and filled with the given one, which happen to be the same pointers in the below case...
		if (this->pMPEG->LoadSignal(fileName, start, end) > 0) {
			pSignal->setBuf_L(this->pMPEG->GetBuf_L(), this->pMPEG->GetLen());
			this->pMPEG->forgetBuf_L();
		} else {
			printf(" <nothing loaded> ");
		}
	}

	return;
}

//====================================================================================================================
//	load the speech-signal between the given segment-boundarys (in samples)
//  install also new debug-prefix, if wished, in the tweakable parameters, according to the fileName
//  if forceSampleRate>0, the signal will be resampled to the given samplerate regardless of it's original sampleRate
//====================================================================================================================
SC_Signal* SC_SignalHandler::loadSignal(const char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, int forceSampleRate) {
  unsigned long int originalStart, originalEnd;
	double srRatio;
	SC_Signal *pResampledSignal = NULL, *pSignal = openSignal(fileName, 0); //intentionally open the signal with its orignal sampleRate to compute the original segment-boundaries (the ones given are based on the forcedSampleRate, but LoadSignal() needs them based on the signals orignal rate)

	//TODO: jStream sources need different handling because there no seeking and no re-opening is possible 
	if (pSignal != NULL) {
		if (segmentEnd > 0 && segmentStart <= segmentEnd) {
			if (pSignal->getSignalType() == sclib::stMP2 || pSignal->getSignalType() == sclib::stMP3 || pSignal->getSignalType() == sclib::stMPEG || pSignal->getSignalType() == sclib::stOther) {
				//SC_Signal_MPEG::LoadSignal() is already able to handle given segment-boundaries based on the wished (rather the original) sampleRate
				loadMPEG(segmentStart, segmentEnd, ((SC_Signal_MPEG*&)pSignal));  //special treatment for SC_Signal_MPEG based signal objects to overcome ffmpeg's problems with seeking
			} else {
				//all other LoadSignal() methods need the boundaries as described above
				if (forceSampleRate > 0 && forceSampleRate != pSignal->SigPar.SRate) {
					srRatio = (double)(forceSampleRate) / (double)(pSignal->SigPar.SRate);
					originalStart = sclib::round(segmentStart / srRatio);
					originalEnd = sclib::round(segmentEnd / srRatio);
				} else {
					originalStart = segmentStart;
					originalEnd = segmentEnd;
				}
				pSignal->LoadSignal(fileName, originalStart, originalEnd);
			}
		} else {
			pSignal->LoadSignal(fileName);
		}
	}

	if (forceSampleRate > 0 && pSignal->SigPar.SRate != forceSampleRate) {
		pResampledSignal = resample(pSignal, (double)(forceSampleRate), true);
		MFree_0D(pSignal);
		pSignal = pResampledSignal;
	}

  return pSignal;
}

//====================================================================================================================
// write an audio-file only incorporating samples belonging to frames between start and end of specific type
// if pSigPar already contains loaded samples of the correct segment-length, it is not tried to load the samples again
// this way, the problem with loading samples from corpora with more than one audio source file (e.g. TIMIT/MAC) is 
// avoided by relocating the signal-loading work to outside this class, where access may be possible to the 
// SC_Corpus_*.loadSignal() method.
// if markGTboundary represents a groundtruth-based boundary-type (speaker/background/artificial) and this type is 
// encountered in the samples to write, it is marked with a short clicking noise in the resulting signal (to be 
// visible in a spectrum analyzer)
//====================================================================================================================
long SC_SignalHandler::storeSignal(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SC_Signal *pSigPar, SC_GroundTruth *pGT, unsigned long int type, unsigned long int typesNot, int origin, unsigned long int markGTboundary, bool uniteTypes, bool uniteTypesNot) {
	unsigned long int x, sampleCount; 
  unsigned short int sType = pSigPar->getSignalType();
	long res = SVLIB_Fail;
	short *samples = NULL, *hook = NULL;
  SC_Signal	*pSignal = NULL;
	char *fName;
	bool justLinked;
	
	if (segmentEnd > pGT->getAudioSampleCount()) {
		segmentEnd = pGT->getAudioSampleCount();
	}

	//load signal-parameters
	switch (sType) {
		case sclib::stWave: {
			pSignal = new SC_Signal_WAVE((SC_Signal_WAVE*)pSigPar);
			break;
		}
		case sclib::stMP2: case sclib::stMPEG: case sclib::stMP3: case sclib::stOther: {
			pSignal = new SC_Signal_MPEG((SC_Signal_MPEG*)pSigPar);
			break;
		}
		case sclib::stNIST: {
			pSignal = new SC_Signal_NIST((SC_Signal_NIST*)pSigPar);
			break;
		}
	}

	//load interesting part of the signal
	if (pSigPar->GetBuf_L() != NULL && pSigPar->GetLen() == segmentEnd-segmentStart+1) {
		pSignal->setBuf_L(pSigPar->GetBuf_L(), pSigPar->GetLen());
		justLinked = true;
	} else {
		pSignal->LoadSignal(pSigPar->getFileName(), segmentStart, segmentEnd);
		justLinked = false;
	}
		
	//copy the desired samples together
	MArray_1D(samples, segmentEnd-segmentStart+1, short, "SC_SignalHandler.storeSignal: samples"); //make it as big as might be necessary (ofton too big, though...)
  sampleCount	= 0;
	hook	= pSignal->GetBuf_L();
	for (x = segmentStart; x <= segmentEnd; x++) {//andFrames was ==false below => never mind, only 1 FLI is under consideration here!!!
    if (pGT->testSegment(x, x, true, type, uniteTypes, typesNot, uniteTypesNot, origin) > 0) { //here we make #internalFrameSize times too many tests, but thats ok for the benefit of not doing any boundary-error here
			if (markGTboundary != sclib::noType && pGT->testSegment(x, x, true, markGTboundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) {
				samples[sampleCount++] = std::numeric_limits<short>::max(); //a "click" sound
			} else {
				samples[sampleCount++] = hook[x-segmentStart];
			}
    } else if (markGTboundary != sclib::noType && pGT->testSegment(x, x, true, markGTboundary, false, sclib::noType, false, sclib::modeGroundtruth) > 0) {
			samples[sampleCount++] = std::numeric_limits<short>::max(); //a "click" sound
		}
	}
	
	if (sampleCount > 0) {
		//write the new audio-file
		fName = new char[strlen(this->pTweak->debug.debugDir) + strlen(fileName) + 1];
		sprintf(fName, "%s%s\0", this->pTweak->debug.debugDir, fileName);
		if (justLinked == true) {
			pSignal->forgetBuf_L();
		}
		pSignal->setBuf_L(samples, sampleCount);
		res	= pSignal->SaveSignal(fName);	
		
    MFree_1D(fName);
		MFree_0D(pSignal);
	} else {
		MFree_0D(samples); //needs only to be freed explicitly if this is not done during destruction of pSignal
	}

	return res;
}

//====================================================================================================================
// write an audio-file containing the samples in the given column of a linked list of sv_data-objects
//====================================================================================================================
long SC_SignalHandler::storeSignal(const char* fileName, SV_Data* pSamples, unsigned short int col, SC_Signal *pSigPar) {
  unsigned long int signalLength = 0, offset = 0;
  long res, i;
  short *signal = NULL;
  SV_Data *pHook = pSamples;

  //determine the full sample count
  while (pHook != NULL) {
    signalLength += pHook->Row;
    pHook = pHook->Next;
  }
  
  //copy all samples together in one array of shorts
  MArray_1D(signal, signalLength, short, "SC_SignalHandler.storeSignal: signal");
  pHook = pSamples;
  while (pHook != NULL) {
    for (i = 0; i < pHook->Row; i++) {
      signal[offset + i] = (short)pHook->Mat[i][col];
    }
    offset += pHook->Row;
    pHook = pHook->Next;
  }
  
  res = storeSignal(fileName, signal, signalLength, pSigPar); //signal is freed in this function
  
  return res;
}

//====================================================================================================================
// write an audio-file containing the samples in the given given row (all rows for row==-1) of a sv_data-object
//====================================================================================================================
long SC_SignalHandler::storeSignal(const char* fileName, SV_Data* pSamples, SC_Signal *pSigPar, int row) {
  unsigned long int offset = 0;
  long res, i, start, end;
  short *signal = NULL;
  SV_Data *pHook = pSamples;
	char *fname = NULL, postfix[sclib::bufferSize];
  
	start = (row >= 0) ? row : 0;
	end = (row >= 0) ? start+1 : pSamples->Row;

	for (int s = start; s < end; s++) {
		MArray_1D(signal, pSamples->Col, short, "SC_SignalHandler.storeSignal: signal"); //copy all samples together in one array of shorts
		for (i = 0; i < pSamples->Col; i++) {
			signal[i] = (short)pSamples->Mat[s][i];
		}
 
		sprintf(postfix, "_%d", s);
		fname = sclib::addPostfixToFilename(fileName, postfix);
	  res = storeSignal(fname, signal, pSamples->Col, pSigPar); //signal is freed in this function
		MFree_1D(fname);
	}
  
  return res;
}

//====================================================================================================================
// write an audio-file containing the samples in the pSignal-buffer of given length
// the signal-array gets destructed by this function!!!
//====================================================================================================================
long SC_SignalHandler::storeSignal(const char* fileName, short* &signal, unsigned long int signalLength, SC_Signal *pSigPar) {
  unsigned short int sType = pSigPar->getSignalType();
  char *fName = new char[strlen(this->pTweak->debug.debugDir) + strlen(fileName) + 1];
  long res;
  SC_Signal	*pNewSignal = NULL;
	sprintf(fName, "%s%s\0", this->pTweak->debug.debugDir, fileName);

  switch (sType) {
    case sclib::stWave: {
      pNewSignal = new SC_Signal_WAVE((SC_Signal_WAVE*)pSigPar);
      break;
    }
		case sclib::stMP2: case sclib::stMPEG: case sclib::stMP3: case sclib::stOther: {
      pNewSignal = new SC_Signal_MPEG((SC_Signal_MPEG*)pSigPar);
      break;
    }
    case sclib::stNIST: {
      pNewSignal = new SC_Signal_NIST((SC_Signal_NIST*)pSigPar);
      break;
    }
  }
	
	pNewSignal->setBuf_L(signal, signalLength);

  res	= pNewSignal->SaveSignal(fName);

  MFree_1D(fName);
  MFree_0D(pNewSignal);

	return res;
}

//====================================================================================================================
//	just a wrapper around the corresponding functions in SC_Signal, so that this class needs not to be exported
//  (gives compiler warnings because the base-class from SV_Lib is not exported... ugly, but works ;-)
//  in addition to the SC_Signal method, this one can be forced to use the debug-dir from pTweak
//====================================================================================================================
long SC_SignalHandler::saveSignal(const char* fileName, SC_Signal *pSignal, bool useDebugDir) {
  char *temp = NULL;
  long res = SVLIB_Fail;

  if (useDebugDir == true) {
    temp = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
    sprintf(temp, "%s%s\0", pTweak->debug.debugDir, fileName);
  } else {
    temp = const_cast<char*>(fileName);
  }
  
  res = pSignal->SaveSignal(temp);

  if (temp != fileName) {
    MFree_1D(temp);
  }

  return res;
}

//====================================================================================================================
//	just a wrapper around the corresponding functions in SC_Signal, so that this class needs not to be exported
//  (gives compiler warnings because the base-class from SV_Lib is not exported... ugly, but works ;-)
//====================================================================================================================
bool SC_SignalHandler::implantSamples(SC_Signal* pSignal, unsigned long startSample, short* samples, unsigned long length) {
  return pSignal->implantSamples(startSample, samples, length);
}

//====================================================================================================================
//	return a new signal object with the given signal resampled to the newSampleRate, with possibly done
//  lowpass-filtering to avoid aliasing; if no resampling is necessary (because current and wished sampleRate already 
//  are the same, a copy of pSignal is returned)
//====================================================================================================================
SC_Signal* SC_SignalHandler::resample(SC_Signal *pSignal, double newSampleRate, bool doFiltering) {
  unsigned long int newLength;
  short *newSignal = NULL;
  SC_Resampling *pResampler;
  SC_Signal *pResampledSignal;
  
	if ((double)(pSignal->SigPar.SRate) != newSampleRate) {
		pResampler = new SC_Resampling(this->pTweak);
		newSignal = pResampler->resample(pSignal->GetBuf_L(), pSignal->GetLen(), pSignal->SigPar.SRate, newSampleRate, newLength, NULL, true, doFiltering);

		pResampledSignal = openSignal(pSignal->getFileName(), 0);
		pResampledSignal->SigPar.SRate = (int)(newSampleRate);
		pResampledSignal->setBuf_L(newSignal, newLength);
		pResampledSignal->setSampleCount(sclib::round(pResampledSignal->getSampleCount() * (newSampleRate / (double)(pSignal->SigPar.SRate))));

		MFree_0D(pResampler);
	} else { //just create a copy in case no resampling is necessary
		newLength = pSignal->GetLen();
		
		if (newLength > 0) { 
			MArray_1D(newSignal, newLength, short, "SC_SignalHandler.resample: newSignal");
			for (int i = 0; i < pSignal->GetLen(); i++) {
				newSignal[i] = pSignal->GetBuf_L()[i];
			}
		}
		
		pResampledSignal = openSignal(pSignal->getFileName(), 0);
		pResampledSignal->setBuf_L(newSignal, newLength);
	}

  return pResampledSignal;
}

//====================================================================================================================
//	wrapper around SC_Signal's setBuf_L() so that this does not need to be dll-exported...
//====================================================================================================================
void SC_SignalHandler::exchangeSignalBuffer(SC_Signal *pSignal, short *newSignal, unsigned long int newSignalLength) {
	pSignal->setBuf_L(newSignal, newSignalLength);

	return;
}
