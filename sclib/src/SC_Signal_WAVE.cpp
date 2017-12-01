/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Signal to load data with RIFF WAVE header (Microsoft Windows */
/*        Audio Waveform format, PCM)																			*/
/*																																				*/
/*    Responsibility:																											*/
/*      - analyses MS-Windows WAV filetype (RIFF WAVE) header and loads		*/
/*        samples according to it																					*/
/*																																				*/
/*    Limitations:																												*/
/*      - only PCM format is supported (uncompressed)											*/
/*      - only 16 sclib::bit bitrate, mono channel is supported									*/
/*			- no byte-order-swapping at the moment, so only use on intel			*/
/*				machines!																												*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.02.2004																								*/
/**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <fstream>
#include "SC_Aux.h"
#include "SC_Signal_WAVE.h"
#include <SV_Error.h>

//====================================================================================================================
// default constructor
//====================================================================================================================
SC_Signal_WAVE::SC_Signal_WAVE() : SC_Signal() {
	this->signalType = sclib::stWave;
	this->io.getCurrentDatatypeSizes(this->sizeInCode);
	createWavFileSizes(this->sizeInFile);

	//some standrad values if this class is just used for write out some hearable samples
  this->SigPar.NChannel = 1;
	this->SigPar.SRate = 16000;
	this->SigPar.StByte = SCLIB_WAVE_HEADER_LENGTH;
	this->SigPar.Encode	= sclib::encodingPCM;
	this->header.bitsPerSample = 16;
	this->header.format = 1;
	this->header.sampleRate = 16000;
	this->header.channelCount = 1;
}

//====================================================================================================================
// this constructor analyses the wav-header without loading samples into memory
//====================================================================================================================
SC_Signal_WAVE::SC_Signal_WAVE(const char* fileName) : SC_Signal() {
	this->signalType = sclib::stWave;
	this->io.getCurrentDatatypeSizes(this->sizeInCode);
	createWavFileSizes(this->sizeInFile);

	//read header
	if (fillHeaderInfo(fileName) == 0) {
 		printf("Load WAVE data failed (can't retrieve header)!");
		this->failure = true;
	};
	
	//Copy header info to derived header
  this->SigPar.NChannel = this->header.channelCount;
	this->SigPar.SRate = this->header.sampleRate;
	this->SigPar.StByte = SCLIB_WAVE_HEADER_LENGTH;

	if (this->header.format != 1) {
		this->SigPar.Encode = 0;   
  } else {	
	  this->SigPar.Encode	= sclib::encodingPCM;
	}
}

//====================================================================================================================
// this is a copy-constructor, which copys the header-info of a given class
//====================================================================================================================
SC_Signal_WAVE::SC_Signal_WAVE(SC_Signal_WAVE* oldSignal) : SC_Signal(oldSignal) {
	this->io.getCurrentDatatypeSizes(this->sizeInCode);
	createWavFileSizes(this->sizeInFile);
	
	//fill header
	this->header = oldSignal->header;
  	
	//Copy header info to derived header
	this->SigPar.StByte = SCLIB_WAVE_HEADER_LENGTH;
	if (oldSignal->header.format != 1) {
		this->SigPar.Encode = 0;   
	} else {	
	  this->SigPar.Encode	= sclib::encodingPCM;
	}
}

//====================================================================================================================
// destructor 
//====================================================================================================================
SC_Signal_WAVE::~SC_Signal_WAVE() {

}

//====================================================================================================================
// create a size-descriptor containing the data-type-sizes (and endianness) in RIFF-WAV compliant files
//====================================================================================================================
void SC_Signal_WAVE::createWavFileSizes(SV_DataIO::SV_DatatypeSizes &sizes) {
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
long SC_Signal_WAVE::LoadSignal(const char *FName) {
	this->failure = false;
  MFree_1D(this->Buf_L);
	this->Len = 0;

  //read header
	if (fillHeaderInfo(FName)==0) {
 		printf("Load WAVE data failed (can't retrieve header)!");
		this->failure = true;
	} else {
		//Copy header info
		this->SigPar.NChannel = this->header.channelCount;
		this->SigPar.SRate = this->header.sampleRate;
		this->SigPar.StByte = SCLIB_WAVE_HEADER_LENGTH;
		if (this->header.format != 1) {
			this->SigPar.Encode = 0;   
		} else {	
			this->SigPar.Encode	= sclib::encodingPCM;
		}

		//allocate mem
		MFree_1D(this->Buf_L);
		MArray_1D(this->Buf_L, this->sampleCount, short, "LoadSignal: Buf_L"); 
		this->Len = this->sampleCount;

		//Read WAVE data into memory
		if (Load_WAVE(FName, this->Buf_L, 0, this->Len) == -1) {
			MFree_1D(this->Buf_L);
			this->Len = 0;
 			printf("Load WAVE data failed (can't retrieve data)!");   
			this->failure = true;
		}
	}

	return(this->Len);
};

//====================================================================================================================
//  Load WAVE format signal from start-sample to end-sample into memory
//  this method expects the header already to be filled by the constructor!              
//====================================================================================================================
long SC_Signal_WAVE::LoadSignal(const char *FName, unsigned long int start, unsigned long int end) {
  this->failure = false;
	
	/*allocate mem*/
  MFree_1D(this->Buf_L);
	MArray_1D(this->Buf_L, end-start+1, short, "LoadSignal: Buf_L"); 
	this->Len = end-start+1;

	/*Read WAVE data into memory*/
	if (Load_WAVE(FName, Buf_L, start, end) == -1) {
		MFree_1D(this->Buf_L);
		this->Len = 0;
 		printf("Load WAVE data failed (can't retrieve data)!");   
		this->failure = true;
	}

	return(this->Len);
};

//====================================================================================================================
//  This procedure can be called to load speech samples from a RIFF WAVE format (Microsoft Windows Waveform  audio) 
//  data file. if success, return the number of samples in SigBuf, otherwise return(-1).                               
//                                                       
//  Speech buf allocated from outside                    
//====================================================================================================================
int SC_Signal_WAVE::Load_WAVE(const char *FName, short *Speech, unsigned long int start, unsigned long int end) {
	fstream inFile;
  unsigned long read = 0;

	//Check if we can use this speech data
	if (this->header.channelCount != 1)		{return(SVLIB_Fail);} //data must be mono
	if (this->header.format != 1)					{return(SVLIB_Fail);} //		-"-			 PCM-format
	if (this->header.bitsPerSample != 16)	{return(SVLIB_Fail);}	//		-"-			 16 sclib::bit
	
	//Open Speech Data File
	inFile.open(FName, ios::in|ios::binary);
	if (inFile.fail()) {
    return(SVLIB_Fail);
  }
	inFile.seekg(this->SigPar.StByte + (start * this->sizeInFile.shortSize), ios_base::beg);

	//load samples (data is always in little endian (intel) format)
	read = this->io.readArray(&inFile, this->Buf_L, this->Len, this->sizeInCode, this->sizeInFile) / this->sizeInFile.shortSize;

	inFile.close();

	if (read != end-start+1) {
		return(SVLIB_Fail);
  }	else {
		return(end-start+1);
	}
};

//====================================================================================================================
//  This procedure read a few important fields in WAVE header. if success, return (1), otherwise return(0)     
//====================================================================================================================
int SC_Signal_WAVE::fillHeaderInfo (const char *FName) {
	fstream inFile;
	char id[4];	//four bytes to hold signatures
	WORD block_align;	//our 16 sclib::bit format info values
	DWORD format_length, avg_bytes_sec;	//our 32 sclib::bit format info values
	unsigned char endianness = this->io.getEndianness();
  
	this->failure = false;
  MFree_1D(this->fileName);
  this->fileName = new char[strlen(FName) + 1];
  strcpy(this->fileName, FName);
  this->fileName[strlen(FName)] = '\0';

	inFile.open(FName, ios::in|ios::binary);
	if (inFile.fail()) {
		//REPORT_ERROR(SVLIB_FileErr, "fillHeaderInfo");
		return 0;
	} 

	this->io.readArray(&inFile, id, 4, this->sizeInCode, this->sizeInFile); //read signature "RIFF"
	if (strncmp(id, "RIFF", 4) != 0) {return 0;}

	this->io.readScalar(&inFile, this->header.fileSize, this->sizeInCode, this->sizeInFile); //read in 32bit size value

	this->io.readArray(&inFile, id, 4, this->sizeInCode, this->sizeInFile); //read signature "WAVE"
	if (strncmp(id,"WAVE", 4) != 0) {return 0;}

	this->io.readArray(&inFile, id, 4, this->sizeInCode, this->sizeInFile); //read signature "fmt " (start of format chunk)
	if(strncmp(id,"fmt ", 4) != 0) {return 0;}

	this->io.readScalar(&inFile, format_length, this->sizeInCode, this->sizeInFile); //read length of format chunk (must be 16 for recognized format)
	if (format_length != 16) {return 0;}

	this->io.readScalar(&inFile, this->header.format, this->sizeInCode, this->sizeInFile); //read format type

	this->io.readScalar(&inFile, this->header.channelCount, this->sizeInCode, this->sizeInFile); //1 mono, 2 stereo

	this->io.readScalar(&inFile, this->header.sampleRate, this->sizeInCode, this->sizeInFile); //like 44100, 22050, etc...

	this->io.readScalar(&inFile, avg_bytes_sec, this->sizeInCode, this->sizeInFile); //some un-needed data...
	this->io.readScalar(&inFile, block_align, this->sizeInCode, this->sizeInFile);

	this->io.readScalar(&inFile, this->header.bitsPerSample, this->sizeInCode, this->sizeInFile); //8 sclib::bit or 16 sclib::bit file?

	this->io.readArray(&inFile, id, 4, this->sizeInCode, this->sizeInFile); //read signature 'data'
	if (strncmp(id,"data", 4) != 0) {return 0;}

	this->io.readScalar(&inFile, this->header.dataLength, this->sizeInCode, this->sizeInFile); //how many bytes of sound data we have

	inFile.close();

	//calculate sampleCount:
	this->sampleCount = this->header.dataLength / (this->header.channelCount * (this->header.bitsPerSample / 8));

	return(1);
}

//====================================================================================================================
//  Generate a new wav-file with the header of this (adpted to the new data) and the samples stored in pSamples
//  Regards the linked-list nature of SC_Signal!
//====================================================================================================================
long SC_Signal_WAVE::SaveSignal(const char *FName) {
	long res;

	//could have been changed form outside...
	this->header.channelCount = this->SigPar.NChannel;
	this->header.sampleRate = this->SigPar.SRate;

	//determine nr. of samples in possibly linked list
	SC_Signal_WAVE *pHook = this;
	int overallSampleCount = 0;
	while (pHook != NULL) {
		overallSampleCount += pHook->GetLen();
		pHook = (SC_Signal_WAVE*)(pHook->Next);
	}
	
	pHook = this;
	res = storeHeader(FName, overallSampleCount);
	while (pHook != NULL) {
		res += storeSamples(FName, pHook->GetBuf_L(), res);
		pHook = (SC_Signal_WAVE*)(pHook->Next);
	}

	return res;
}

//====================================================================================================================
//  Save the current header to a file
//====================================================================================================================
long SC_Signal_WAVE::storeHeader(const char* fileName, int numSamples) {
	fstream OutFile;
	int res = 0;
	int sampleCount = (numSamples==0) ? this->GetLen() : numSamples;

  OutFile.open(fileName, ios_base::out|ios_base::binary);

	res += this->io.writeArray(&OutFile, "RIFF", 4, this->sizeInCode, this->sizeInFile);
	
	res += this->io.writeScalar(&OutFile, (DWORD)(SCLIB_WAVE_HEADER_LENGTH - 8 + (GetLen()*this->sizeInFile.shortSize)), this->sizeInCode, this->sizeInFile); //riff chunk length (filesize - 8)

	res += this->io.writeArray(&OutFile, "WAVE", 4, this->sizeInCode, this->sizeInFile);
	
	res += this->io.writeArray(&OutFile, "fmt ", 4, this->sizeInCode, this->sizeInFile);
	
	res += this->io.writeScalar(&OutFile, (DWORD)(16), this->sizeInCode, this->sizeInFile); //format chunk length

	res += this->io.writeScalar(&OutFile, (WORD)(this->header.format), this->sizeInCode, this->sizeInFile); //format tag/channels
	res += this->io.writeScalar(&OutFile, (WORD)(this->header.channelCount), this->sizeInCode, this->sizeInFile);
	
	res += this->io.writeScalar(&OutFile, (DWORD)(this->header.sampleRate), this->sizeInCode, this->sizeInFile); //sample-rate

	res += this->io.writeScalar(&OutFile, (DWORD)( this->header.sampleRate*this->header.channelCount*(this->header.bitsPerSample/8)), this->sizeInCode, this->sizeInFile); //sample-rate*channels*(Bits/8)

	res += this->io.writeScalar(&OutFile, (WORD)((this->header.bitsPerSample/8)*this->header.channelCount), this->sizeInCode, this->sizeInFile); //block-align/bits per sample
	res += this->io.writeScalar(&OutFile, (WORD)(this->header.bitsPerSample), this->sizeInCode, this->sizeInFile);
	
	res += this->io.writeArray(&OutFile, "data ", 4, this->sizeInCode, this->sizeInFile);
	
	res += this->io.writeScalar(&OutFile, (DWORD)(sampleCount*this->sizeInFile.shortSize), this->sizeInCode, this->sizeInFile); //actual length of raw data

	if (OutFile.fail() || res != SCLIB_WAVE_HEADER_LENGTH) {
    return SVLIB_Fail;
  }
	OutFile.close();
	
	return res;
}

//====================================================================================================================
//  Save the samples to a file
//====================================================================================================================
long SC_Signal_WAVE::storeSamples(const char* fileName, short* pSamples, int seekTo) {
	unsigned char endianness = this->io.getEndianness();
	int res;
	fstream OutFile;

  OutFile.open(fileName, ios_base::out|ios_base::app|ios_base::binary); //ios::out|ios::binary
	OutFile.seekg(seekTo); //SCLIB_WAVE_HEADER_LENGTH
	res = this->io.writeArray(&OutFile, this->Buf_L, this->Len, this->sizeInCode, this->sizeInFile);

  if (OutFile.fail() || res != this->Len*this->sizeInFile.shortSize) {
    return SVLIB_Fail;
  }
	OutFile.close();

	return (res);
}

//====================================================================================================================
//  fill the header; fileSize and dataSize are only set if != 0
//====================================================================================================================
void SC_Signal_WAVE::setHeader(unsigned long fileSize, unsigned short format, unsigned short channelCount, unsigned long sampleRate, unsigned short bitsPerSample, unsigned long dataLength) {
  if (fileSize > 0) {
    this->header.fileSize = fileSize;
  }
  this->header.format = format;
  this->header.channelCount = channelCount;
	this->SigPar.NChannel = channelCount;
  this->header.sampleRate = sampleRate;
	this->SigPar.SRate = sampleRate;
  this->header.bitsPerSample = bitsPerSample;
  if (dataLength > 0) {
    this->header.dataLength = dataLength;
  }
}
