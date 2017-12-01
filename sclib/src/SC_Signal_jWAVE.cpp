/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Signal_jWAVE to load wav files via Java-streams							*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.01.2009																								*/
/**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <fstream>
#include "SC_Aux.h"
#include "SC_Signal_jWAVE.h"
#include "SC_Signal_WAVE.h"
#include <SV_Error.h>

//====================================================================================================================
// default constructor
//====================================================================================================================
SC_Signal_jWAVE::SC_Signal_jWAVE() : SC_Signal(), headerLength(44) {
	this->signalType = sclib::stJWave;
	this->io.getCurrentDatatypeSizes(this->sizeInCode);
	createWavFileSizes(this->sizeInFile);

	//some standrad values if this class is just used for write out some hearable samples
  this->SigPar.NChannel = 1;
	this->SigPar.SRate = 16000;
	this->SigPar.StByte = this->headerLength;
	this->SigPar.Encode	= sclib::encodingPCM;
	this->header.bitsPerSample = 16;
	this->header.format = 1;
	this->header.sampleRate = 16000;
	this->header.channelCount = 1;
}

//====================================================================================================================
// this constructor analyses the wav-header without loading samples into memory
//====================================================================================================================
#ifdef SC_USE_JNI
SC_Signal_jWAVE::SC_Signal_jWAVE(JNIEnv *env, jobject jStreamObject) : SC_Signal(), headerLength(44) {
	this->signalType = sclib::stWave;
	this->io.getCurrentDatatypeSizes(this->sizeInCode);
	createWavFileSizes(this->sizeInFile);
	this->jio.setEnvironment(env, jStreamObject);

	if (this->jio.isInitialized() != true) {
 		printf("Load WAVE data failed (jStream not initialized)!");
		this->failure = true;
	} else {
		//read header
		if (fillHeaderInfo() == 0) {
 			printf("Load WAVE data failed (can't retrieve header)!");
			this->failure = true;
		}

		//Copy header info to derived header
		this->SigPar.NChannel = this->header.channelCount;
		this->SigPar.SRate = this->header.sampleRate;
		this->SigPar.StByte = this->headerLength;
		if (this->header.format != 1) {
			this->SigPar.Encode = 0;   
		} else {	
			this->SigPar.Encode	= sclib::encodingPCM;
		}
	}
}
#endif

//====================================================================================================================
// this is a copy-constructor, which copys the header-info of a given class
//====================================================================================================================
SC_Signal_jWAVE::SC_Signal_jWAVE(SC_Signal_jWAVE* oldSignal) : SC_Signal(oldSignal), headerLength(44) {
	this->io.getCurrentDatatypeSizes(this->sizeInCode);
	createWavFileSizes(this->sizeInFile);
#ifdef SC_USE_JNI
	this->jio.setEnvironment(oldSignal->jio.getEnvironment(), oldSignal->jio.getStream());
#endif
	
	//fill header
	this->header = oldSignal->header;
  	
	//Copy header info to derived header
	this->SigPar.StByte = this->headerLength;
	if (oldSignal->header.format != 1) {
		this->SigPar.Encode = 0;   
	} else {	
	  this->SigPar.Encode	= sclib::encodingPCM;
	}
}

//====================================================================================================================
// destructor 
//====================================================================================================================
SC_Signal_jWAVE::~SC_Signal_jWAVE() {

}

#ifdef SC_USE_JNI
//====================================================================================================================
// sets the stream object from which to read
//====================================================================================================================
void SC_Signal_jWAVE::setStream(JNIEnv *env, jobject jStreamObject) {
	this->jio.setEnvironment(env, jStreamObject);

	return;
}
#endif

//====================================================================================================================
// create a size-descriptor containing the data-type-sizes (and endianness) in RIFF-WAV compliant files
//====================================================================================================================
void SC_Signal_jWAVE::createWavFileSizes(SV_DataIO::SV_DatatypeSizes &sizes) {
	sizes.boolSize = sizeof(BYTE);
	sizes.charSize = sizeof(BYTE);
	sizes.shortSize = sizeof(WORD);
	sizes.intSize = sizeof(DWORD);
	sizes.endianness = SVLIB_LITTLE_ENDIAN;
	
	//all other data types don't occur in wav-files (at least not in the parts parseable so far by this class...)
	//setting their sizes to 0 will raise errors if the system nontheless happens to use them
	sizes.longSize = 0;
	sizes.longLongSize = 0;
	sizes.floatSize = 0;
	sizes.doubleSize = 0;
	sizes.longDoubleSize = 0;
	
	return;
}

//====================================================================================================================
//  Load comlete WAVE format signal into memory              
//====================================================================================================================
long SC_Signal_jWAVE::LoadSignal(const char *FName) {
	this->failure = false;
  MFree_1D(this->Buf_L);
	this->Len = 0;

#ifdef SC_USE_JNI
	if (this->jio.isInitialized() != true) {
 		printf("Load WAVE data failed (jStream not initialized)!");
		this->failure = true;
	} else {
		/*read header*/
		if (fillHeaderInfo()==0) {
 			printf("Load WAVE data failed (can't retrieve header)!");
			this->failure = true;
		} else {
			/*Copy header info*/
			this->SigPar.NChannel = this->header.channelCount;
			this->SigPar.SRate = this->header.sampleRate;
			this->SigPar.StByte = this->headerLength;
			if (this->header.format != 1) {
				this->SigPar.Encode = 0;   
			} else {	
				this->SigPar.Encode	= sclib::encodingPCM;
			}

			/*allocate mem*/
			MArray_1D(this->Buf_L, this->sampleCount, short, "LoadSignal: Buf_L"); 
			this->Len = this->sampleCount;

			/*Read WAVE data into memory*/
			if (Load_WAVE(this->Buf_L, 0, this->Len) == -1) {
				MFree_1D(this->Buf_L);
				this->Len = 0;
 				printf("Load WAVE data failed (can't retrieve data)!");   
				this->failure = true;
			}
		}
	}
#endif

	return(this->Len);
};

//====================================================================================================================
//  Load WAVE format signal from start-sample to end-sample into memory
//  this method expects the header already to be filled by the constructor!              
//====================================================================================================================
long SC_Signal_jWAVE::LoadSignal(const char *FName, unsigned long int start, unsigned long int end) {
  this->failure = false;
  MFree_1D(this->Buf_L);
	this->Len = 0;

#ifdef SC_USE_JNI
	if (this->jio.isInitialized() != true) {
 		printf("Load WAVE data failed (jStream not initialized)!");
		this->failure = true;
	} else {
		//allocate mem
		MFree_1D(this->Buf_L);
		MArray_1D(this->Buf_L, end-start+1, short, "LoadSignal: Buf_L"); 
		this->Len = end-start+1;

		//Read WAVE data into memory
		if (Load_WAVE(Buf_L, start, end) == -1) {
			MFree_1D(this->Buf_L);
			this->Len = 0;
 			printf("Load WAVE data failed (can't retrieve data)!");   
			this->failure = true;
		}
	}
#endif

	return(this->Len);
};

//====================================================================================================================
//  This procedure can be called to load speech samples from a RIFF WAVE format (Microsoft Windows Waveform  audio) 
//  data file. if success, return the number of samples in SigBuf, otherwise return(-1).                               
//                                                       
//  Speech buf allocated from outside                    
//====================================================================================================================
int SC_Signal_jWAVE::Load_WAVE(short *Speech, unsigned long int start, unsigned long int end) {
  unsigned long read = 0;

	//Check if we can use this speech data
	if (this->header.channelCount != 1)		{return(SVLIB_Fail);} //data must be mono
	if (this->header.format != 1)					{return(SVLIB_Fail);} //		-"-			 PCM-format
	if (this->header.bitsPerSample != 16)	{return(SVLIB_Fail);}	//		-"-			 16 sclib::bit

#ifdef SC_USE_JNI
	//load samples (data is always in little endian (intel) format)
	read = this->jio.readArray(this->Buf_L, this->Len, this->sizeInCode, this->sizeInFile) / this->sizeInFile.shortSize;
#endif

	if (read != end-start) {
		return(SVLIB_Fail);
  }	else {
		return(end-start);
	}
};

//====================================================================================================================
//  This procedure read a few important fields in WAVE header. if success, return (1), otherwise return(0)     
//====================================================================================================================
int SC_Signal_jWAVE::fillHeaderInfo() {
#ifdef SC_USE_JNI
	char id[4];	//four bytes to hold signatures
	WORD block_align;	//our 16 sclib::bit format info values
	DWORD format_length, avg_bytes_sec;	//our 32 sclib::bit format info values
	unsigned char endianness = this->io.getEndianness();
#endif

	this->failure = false;
  MFree_1D(this->fileName);
  this->fileName = new char[8];
  sprintf(this->fileName, "jStream");

#ifdef SC_USE_JNI
	this->jio.readArray(id, 4, this->sizeInCode, this->sizeInFile); //read signature "RIFF"
	if (strncmp(id, "RIFF", 4) != 0) {return 0;}

	this->jio.readScalar(this->header.fileSize, this->sizeInCode, this->sizeInFile); //read in 32bit size value

	this->jio.readArray(id, 4, this->sizeInCode, this->sizeInFile); //read signature "WAVE"
	if (strncmp(id,"WAVE", 4) != 0) {return 0;}

	this->jio.readArray(id, 4, this->sizeInCode, this->sizeInFile); //read signature "fmt " (start of format chunk)
	if(strncmp(id,"fmt ", 4) != 0) {return 0;}

	this->jio.readScalar(format_length, this->sizeInCode, this->sizeInFile); //read length of format chunk (must be 16 for recognized format)
	if (format_length != 16) {return 0;}

	this->jio.readScalar(this->header.format, this->sizeInCode, this->sizeInFile); //read format type

	this->jio.readScalar(this->header.channelCount, this->sizeInCode, this->sizeInFile); //1 mono, 2 stereo

	this->jio.readScalar(this->header.sampleRate, this->sizeInCode, this->sizeInFile); //like 44100, 22050, etc...

	this->jio.readScalar(avg_bytes_sec, this->sizeInCode, this->sizeInFile); //some un-needed data...
	this->jio.readScalar(block_align, this->sizeInCode, this->sizeInFile);

	this->jio.readScalar(this->header.bitsPerSample, this->sizeInCode, this->sizeInFile); //8 sclib::bit or 16 sclib::bit file?

	this->jio.readArray(id, 4, this->sizeInCode, this->sizeInFile); //read signature 'data'
	if (strncmp(id,"data", 4) != 0) {return 0;}

	this->jio.readScalar(this->header.dataLength, this->sizeInCode, this->sizeInFile); //how many bytes of sound data we have

	//calculate sampleCount:
	this->sampleCount = this->header.dataLength / (this->header.channelCount * (this->header.bitsPerSample / 8));
#endif

	return(1);
}

//====================================================================================================================
//  Generate a new wav-file with the header of this (adpted to the new data) and the samples stored in pSamples
//====================================================================================================================
long SC_Signal_jWAVE::SaveSignal(const char *FName) {
  long int res = 0;
  SC_Signal_WAVE *pWaveWriter = new SC_Signal_WAVE();

  //create a WAV-file header only with the necessary information to create a new WAV file
  pWaveWriter->setHeader(0, 1, this->SigPar.NChannel, this->SigPar.SRate, 16, 0); //use the SigPar.* parameters here (instead of header.*) that could have been modified from outside (together with resampling)

	SC_Signal_WAVE *pWaveHook = pWaveWriter;
	SC_Signal_jWAVE *pHook = this;
	while (pHook != NULL) {
		pWaveHook->setBuf_L(pHook->GetBuf_L(), pHook->GetLen());
		pHook = (SC_Signal_jWAVE*)(pHook->Next);
		if (pHook != NULL) {
			pWaveHook->Next = new SC_Signal_WAVE();
			pWaveHook = (SC_Signal_WAVE*)(pWaveHook->Next);
		}
	}	

  //pWaveWriter->setBuf_L(this->Buf_L, this->Len);
  res = pWaveWriter->SaveSignal(FName);
  //pWaveWriter->forgetBuf_L();

	pWaveHook = pWaveWriter;
	while (pWaveHook != NULL) {
		pWaveHook->forgetBuf_L();
		pWaveHook = (SC_Signal_WAVE*)(pWaveHook->Next);
	}

  return res;
}
