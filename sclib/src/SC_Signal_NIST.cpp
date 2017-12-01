/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Signal to load SPHERE headered(NIST) files with some comfort */
/*                                                                        */
/*    Also based heavily on Jialong He's SV_Signal_NIST class             */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 07.04.2006																								*/
/**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "SC_Signal_NIST.h"
#include "SC_Signal_WAVE.h"
#include "SC_Aux.h"
#include <SV_Error.h>

#define SCLIB_NIST_HEADER_LENGTH 1024

//====================================================================================================================
// default constructor
//====================================================================================================================
SC_Signal_NIST::SC_Signal_NIST() : SC_Signal() {
  this->fileName = NULL;
	this->signalType = sclib::stNIST;
}

//====================================================================================================================
// this constructor analyses the wav-header without loading samples into memory
//====================================================================================================================
SC_Signal_NIST::SC_Signal_NIST(const char* fileName) : SC_Signal() {
	if (loadHeader(fileName, this->header) == 0) { //read header
 		printf("Load NIST data failed (can't retrieve header)!");   
		this->failure = true;
  } else { //success -> copy to internal data-structures
    //copy fileName
    MFree_1D(this->fileName);
    this->fileName = new char[strlen(fileName) + 1];
    strcpy(this->fileName, fileName);
    this->fileName[strlen(fileName)] = '\0';

	  //Copy header info to derived header
    this->SigPar.Encode = this->header.sample_coding;
	  this->SigPar.SRate = this->header.sample_rate;
    this->SigPar.NChannel = this->header.channel_count;
   	this->SigPar.StByte = SCLIB_NIST_HEADER_LENGTH; //sart-byte
  
    //sample-count
    this->sampleCount = header.sample_count;

		//signal type
		this->signalType = sclib::stNIST;
  }
}

//====================================================================================================================
// this is a copy-constructor, which copys the header-info of a given class
//====================================================================================================================
SC_Signal_NIST::SC_Signal_NIST(SC_Signal_NIST* oldSignal) : SC_Signal(oldSignal) {
	this->header = oldSignal->header;
}

//====================================================================================================================
// destructor 
//====================================================================================================================
SC_Signal_NIST::~SC_Signal_NIST() {

}

//====================================================================================================================
// decompress SHORTEN data (this function can be found in de_short.cpp in the svlib libraray)
//====================================================================================================================
extern int de_short(FILE *filei, char *OutBuf, long BufSize);

//====================================================================================================================
//  Load complte NIST format signal into memory
//====================================================================================================================
long SC_Signal_NIST::LoadSignal (const char *FName) {
	this->failure = false;
	
	if (loadHeader(FName, this->header) == 0) { //Read Sphere header info

    printf("Load NIST headers failed!");   
		MFree_1D(this->Buf_L);
		this->Len = 0;
		this->failure = true;

  } else {

    //fill meta-data: fileName & sampleCount
    MFree_1D(this->fileName);
    this->fileName = new char[strlen(FName) + 1];
    strcpy(this->fileName, FName);
    this->fileName[strlen(FName)] = '\0';
	  this->sampleCount = this->header.sample_count;

	  //Copy header info 
	  SigPar.NChannel = this->header.channel_count;
	  SigPar.SRate = this->header.sample_rate;
	  SigPar.StByte = SCLIB_NIST_HEADER_LENGTH;
	  SigPar.Encode = this->header.sample_coding;

	  //Read complete NIST data into memory
	  if (LoadNIST(FName, 0, this->Len, this->header, this->Buf_L, this->Len) == SVLIB_Fail) {
		  MFree_1D(this->fileName);
      this->sampleCount = 0;
 	    printf("Load NIST data failed!");   
			this->failure = true;
	  }
  }

	return(this->Len);
}

//====================================================================================================================
//  Load NIST format signal from start-sample to end-sample into memory
//====================================================================================================================
long SC_Signal_NIST::LoadSignal(const char *FName, unsigned long int start, unsigned long int end) {
	this->failure = false;
	
	if (loadHeader(FName, this->header) ==0) { //Read Sphere header info

    printf("Load NIST headers failed!");   
		MFree_1D(Buf_L);
		this->Len = 0;
		this->failure  = true;

  } else {

    //fill meta-data: fileName & sampleCount
    MFree_1D(this->fileName);
    this->fileName = new char[strlen(FName) + 1];
    strcpy(this->fileName, FName);
    this->fileName[strlen(FName)] = '\0';
	  this->sampleCount = this->header.sample_count;

	  //Copy header info 
	  SigPar.NChannel = this->header.channel_count;
	  SigPar.SRate = this->header.sample_rate;
	  SigPar.StByte = SCLIB_NIST_HEADER_LENGTH;
	  SigPar.Encode = this->header.sample_coding;

	  //Read complete NIST data into memory
	  if (LoadNIST(FName, start, end, this->header, this->Buf_L, this->Len) == SVLIB_Fail) {
		  MFree_1D(this->fileName);
      this->sampleCount = 0;
 	    printf("Load NIST data failed!");   
	  }
  }

	return(this->Len);
}

//====================================================================================================================
//  This procedure can be called to load speech samples from NIST's SPHERE format data file.
//  if sucess, return the number of samples in SigBuf, otherwise, return(SVLIB_Fail).
//  Speech buf is allocated here!!!
//====================================================================================================================
int SC_Signal_NIST::LoadNIST(const char *FName, unsigned long int start, unsigned long int end, SC_Signal_NIST::SC_HeaderTypeSPHERE &HInfo, short* &buffer, int &length) {
  FILE *fid;
	fstream fid2;
	char *pByte;
  long read = 0, Cnt, i = 0;
  short *tempBuf = NULL; //buffer to temporary hold the complete signal
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;

	this->io.getCurrentDatatypeSizes(codeSizes);
	fileSizes.shortSize = (unsigned char)(this->header.sample_n_bytes); //the only two informations needed are the sizeOf(short) and the file's endianness
	fileSizes.endianness = (this->header.SwapByte == false) ? codeSizes.endianness : (codeSizes.endianness == SVLIB_LITTLE_ENDIAN) ? SVLIB_BIG_ENDIAN : SVLIB_LITTLE_ENDIAN;

  //reset everything
  MFree_1D(buffer);
  length = 0;

	//Check if we can use this speech data
	if (end > (unsigned long)HInfo.sample_count || start > end || HInfo.channel_count != 1) {
		return(SVLIB_Fail);
  }

	//Open Speech Data File
	if (HInfo.sample_coding == PCM_2) {
		fid2.open(FName, ios::in|ios::binary);
		if (fid2.fail()) {
			return SVLIB_Fail;
		}
		fid2.seekg(SCLIB_NIST_HEADER_LENGTH, ios_base::beg);
	} else {
	  fid = fopen(FName, "rb");
	  if (fid==NULL) return(SVLIB_Fail);
	  fseek(fid, SCLIB_NIST_HEADER_LENGTH, SEEK_SET);
	}

  //we can't load just the desired portion because the signal maybe needs to be decompressed first
  //so we load the complete signal, convert it to plain samples and then cut the desired part out of it
  MArray_1D(tempBuf, HInfo.sample_count, short, "SC_Signal_NIST.LoadNIST: tempBuf");

	//depending sample_coding, load samples
	switch (HInfo.sample_coding) {
	  case PCM_2 :
			read = this->io.readArray(&fid2, tempBuf, HInfo.sample_count, codeSizes, fileSizes) / fileSizes.shortSize; //at least this we can alter...
	    break;

	  case PCM_2_Shorten :
	    read = de_short(fid, (char*)tempBuf, HInfo.sample_count*2);
      read /= 2;   /* returned is bytes not samples */
		  if (HInfo.SwapByte) {
				this->io.swapBytes((char*)tempBuf, 2, HInfo.sample_count*2);
		    //swab((char*)tempBuf, (char*)tempBuf, HInfo.sample_count*2);
		  }
		  break;   // de_short already swaped depend on byte order.

	  case ULAW :
			pByte = (char*)tempBuf;
	    read = (long)fread(pByte, 1, HInfo.sample_count, fid);
      for (Cnt=read-1; Cnt<=0; Cnt--) {
			  tempBuf[Cnt] = (short)ulaw2linear(pByte[Cnt]);
		  }
		  break;

	  case ULAW_Shorten :
		  MArray_1D(pByte, HInfo.sample_count*2, char, "TmpBuf");
	    read = de_short(fid, pByte, HInfo.sample_count*2);
      for (Cnt=0; Cnt<read; Cnt++) tempBuf[Cnt] = ulaw2linear(pByte[Cnt]); 
		  MFree_1D(pByte);
	    break;

	  default: ;
	}  // switch

	if (HInfo.sample_coding == PCM_2) {
		fid2.close();
	} else {
		fclose(fid);
	}

	if (read != HInfo.sample_count) {
    MFree_1D(tempBuf);
    return(SVLIB_Fail);
  }	else {
    if (start == 0 && end == HInfo.sample_count) { //no need to copy here, we want the complete signal
      buffer = tempBuf;
      length = HInfo.sample_count;
    } else { //copy the desired part out of tempbuffer and free that thing
      MArray_1D(buffer, end-start+1, short, "SC_Signal_NIST.LoadNIST: buffer");
      for (Cnt = start; Cnt <= (long)end; Cnt++) { //TODO: something is wrong here with the borders...
        buffer[i++] = tempBuf[Cnt];
      }
      length = end - start + 1;
			if (i != length) {
				printf("Wrong number of samples read!");
				this->failure = true;
			}
      MFree_1D(tempBuf);
    }
    return length;
  }
}

//====================================================================================================================
//  This procedure read a few important fields in SPHERE header. if success, return (1), otherwise return(0) 
//====================================================================================================================
int SC_Signal_NIST::loadHeader(const char *FName, SC_Signal_NIST::SC_HeaderTypeSPHERE &HInfo) {
  ifstream InFile;
	char header[SCLIB_NIST_HEADER_LENGTH], *Field;
	short Int2;
	long Int4;
	char StrV1[100], StrV2[100], StrV3[100];

	this->failure = false;

  //Read header into memory
	InFile.open(FName, ios::in|ios::binary);
	InFile.read(header, SCLIB_NIST_HEADER_LENGTH);
	if (InFile.fail()) {
		//REPORT_ERROR(SVLIB_FileErr, "HeaderInfo"); 	
		return 0;
	}

	InFile.close();

	//analysis header contents
	header[1023] = 0;  //make header a valid string
	if (strncmp(header, "NIST_1A", 7) != 0) return(0);

	//channel_count
	Field = strstr(header, "channel_count");
	if (Field == NULL) return(0);

	sscanf(Field, "%s %s %hd", StrV1, StrV2, &Int2);
	HInfo.channel_count = Int2;

	//sample_count
	Field = strstr(header, "sample_count");
	if (Field == NULL) return (0);

	sscanf(Field, "%s %s %ld", StrV1, StrV2, &Int4);
	HInfo.sample_count = Int4;

	//sample_rate
	Field = strstr(header, "sample_rate");
	if (Field == NULL) return (0);

	sscanf(Field, "%s %s %ld", StrV1, StrV2, &Int4);
	HInfo.sample_rate = Int4;

	//sample_n_bytes
	Field = strstr(header, "sample_n_bytes");
	if (Field == NULL) return(0);

	sscanf(Field, "%s %s %hd", StrV1, StrV2, &Int2);
	HInfo.sample_n_bytes = Int2;

	//sample_byte_format
	Field = strstr(header, "sample_byte_format");
	if (Field != NULL) {
	  sscanf (Field, "%s %s %s", StrV1, StrV2, StrV3);
	  if ( (IsBigEndian() && strcmp(StrV3, "10")==0) || (!IsBigEndian() && strcmp(StrV3, "01")==0)) {
      HInfo.SwapByte = 0;
    }	else {
      HInfo.SwapByte = 1;
    }
	} else {
    HInfo.SwapByte = 0; //suppose native format
  }

  //sample_coding
	Field = strstr(header, "sample_coding");
	if (Field != NULL) {
		sscanf(Field, "%s %s %s", StrV1, StrV2, StrV3);
		
		if (strcmp(StrV3, "pcm")==0) //determine encoding string contents (modify this part to recognize more types)
			HInfo.sample_coding = PCM_2;
		else if (strstr(StrV3, "pcm") && strstr(StrV3, "shorten"))
			HInfo.sample_coding = PCM_2_Shorten;
	  else if (strcmp(StrV3, "ulaw") == 0 || strcmp(StrV3, "mu-law") == 0)
			HInfo.sample_coding = ULAW;
		else if (strstr(StrV3, "ulaw") && strstr(StrV3, "shorten"))
			HInfo.sample_coding = ULAW_Shorten;
		else
			HInfo.sample_coding = UNKNOWN;
	} else {
    HInfo.sample_coding = PCM_2;  //default is two byte PCM
  }

	return(1); //Have all need info
}

//====================================================================================================================
//  Generate a new wav (!!!) file with samples stored in Buf_L
//  because we don't want to encode a new SPHERE file here, we just create a new WAV file instead
//====================================================================================================================
long SC_Signal_NIST::SaveSignal(const char *FName) {
  long int res = 0;
  SC_Signal_WAVE *pWaveWriter = new SC_Signal_WAVE();

  //create a WAV-file header only with the necessary information to create a new WAV file
  pWaveWriter->setHeader(0, 1, this->SigPar.NChannel, this->SigPar.SRate, 16, 0); //use the SigPar.* parameters here (instead of header.*) that could have been modified from outside (together with resampling)

	SC_Signal_WAVE *pWaveHook = pWaveWriter;
	SC_Signal_NIST *pHook = this;
	while (pHook != NULL) {
		pWaveHook->setBuf_L(pHook->GetBuf_L(), pHook->GetLen());
		pHook = (SC_Signal_NIST*)(pHook->Next);
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
