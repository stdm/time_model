//************************************************************************
//    Loading NIST format speech signal into memory
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "SV_Signal_NIST.h"
#include "SV_Error.h"

using namespace std;

//block by thilo: _swab() is not available in gcc, but swab() is deprectaded in msvc++
#ifdef __GNUC__
	#ifndef _swab
		#define _swab swab
	#endif
#endif
//end by thilo


static char SV_LibID[] = "Copyright (c) by Jialong He";

//--------------------------------------
// decompress SHORTEN data
//--------------------------------------
int de_short(FILE *filei, char *OutBuf, long BufSize);


//==========================================
// default constructor
//==========================================
SV_Signal_NIST::SV_Signal_NIST() : SV_Signal(){


}



//==========================================
// destructor 
//==========================================
SV_Signal_NIST::~SV_Signal_NIST() {


}

//=======================================================*/
//  Load NIST format signal into memory                  */
//=======================================================*/
long SV_Signal_NIST::LoadSignal (char *FName) {

	HeaderType HInfo;

	//------------------------------------
	// Read Sphere header info
	//------------------------------------ 
	if (SPHERE_HeaderInfo (FName, &HInfo)==0) {
 	    REPORT_ERROR(SVLIB_BadData, "Load NIST data failed!");   
	};

	//------------------------------------
	// Copy header info and allocate mem
	//------------------------------------
	SigPar.NChannel  = HInfo.channel_count;
	SigPar.SRate     = HInfo.sample_rate;
	SigPar.StByte    = HLen;
	SigPar.Encode	 = HInfo.sample_coding;

	if (HInfo.sample_count != Len) {
		MFree_1D(Buf_L);
		MArray_1D(Buf_L, HInfo.sample_count, short, "Buf_L"); 
		Len = HInfo.sample_count;
	}

	//------------------------------------
	// Read NIST data into memory
	//------------------------------------
	if (LoadNIST(FName, Buf_L, HInfo) == -1) {
		MFree_1D(Buf_L);
		Len = 0;
 	    REPORT_ERROR(SVLIB_BadData, "Load NIST data failed!");   
	}

	return(Len);
};


/*=======================================================*/
/*  This procedure can be called to load speech samples  */
/*  from NIST's SPHERE format data file.                 */
/*                                                       */
/*  if sucess, return the number of samples in SigBuf    */
/*  otherwise, return(-1).                               */
/*                                                       */
/*  Speech buf allocated from outside                    */
/*=======================================================*/
int SV_Signal_NIST::LoadNIST(char *FName, short *Speech, HeaderType &HInfo) {

    FILE *fid;
	char  *pByte;
    long  Readed=0, Cnt;

	/*--------------------------------------*/
	/* Check if we can use this speech data */
	/*--------------------------------------*/
	if (HInfo.channel_count != 1) {
		return(-1);
	};

	/*--------------------------------------*/
	/* Open Speech Data File                */
	/*--------------------------------------*/
	fid = fopen(FName, "rb");
	if (fid==NULL) return(-1);
	fseek(fid, HLen, SEEK_SET);

	/*------------------------------------------*/
	/*  depending sample_coding, load samples   */
	/*------------------------------------------*/
	switch (HInfo.sample_coding) {
	  case PCM_2 :
	       Readed = (long)fread(Speech, 2, HInfo.sample_count, fid);
		   if (HInfo.SwapByte) {
		       _swab((char*)Speech, (char*)Speech, HInfo.sample_count*2); //by thilo: swab() to _swab() (POSIX function is deprecated since vc2005, new one is ISO C++ conformant)
		   }
	       break;

	  case PCM_2_Shorten :
	       Readed = de_short(fid, (char*)Speech, HInfo.sample_count*2);
           Readed /= 2;   /* returned is bytes not samples */
		   if (HInfo.SwapByte) {
		       _swab((char*)Speech, (char*)Speech, HInfo.sample_count*2); //by thilo: swab() to _swab() (POSIX function is deprecated since vc2005, new one is ISO C++ conformant)
		   }
		   break;   // de_short already swaped depend on byte order.

	  case ULAW :
			pByte = (char*) Speech;
	       Readed = (long)fread(pByte, 1, HInfo.sample_count, fid);
           for (Cnt=Readed-1; Cnt<=0; Cnt--) {
			   Speech[Cnt] = (short)ulaw2linear(pByte[Cnt]);
		   }
		   
		   break;

	  case ULAW_Shorten :
		   MArray_1D(pByte, HInfo.sample_count*2, char, "TmpBuf");
	       Readed = de_short(fid, pByte, HInfo.sample_count*2);
           for (Cnt=0; Cnt<Readed; Cnt++) Speech[Cnt] = ulaw2linear(pByte[Cnt]); 
		   MFree_1D(pByte);
	       break;

	  default: ;
	}  // switch


	fclose(fid);

	if (Readed != HInfo.sample_count) {
		return(-1);
    }
	else return(HInfo.sample_count);

};


/*=======================================================*/
/*  This procedure read a few inportant fields in SPHERE */
/*  header. if success, return (1), otherwise return(0)  */
/*=======================================================*/
int SV_Signal_NIST::SPHERE_HeaderInfo (char *FName, HeaderType *HInfo) {

	ifstream InFile;

	char header[HLen], *Field;
	short   Int2;
	long    Int4;
	char    StrV1[100], StrV2[100], StrV3[100];

    /*----------------------------*/
    /* Read header into memory    */
    /*----------------------------*/
	InFile.open(FName, ios::in|ios::binary);
	InFile.read(header, HLen);
	if (InFile.fail()) {
		REPORT_ERROR(SVLIB_FileErr, "HeaderInfo"); 	
	}

	InFile.close();

	/*----------------------------*/
	/* analysis header contents   */
	/*----------------------------*/
	header[1023] = 0;   /* make header a valid string */
	if (strncmp(header, "NIST_1A", 7) != 0) return(0);

	/*---- channel_count -----*/
	Field = strstr(header, "channel_count");
	if (Field == NULL) return(0);

	sscanf(Field, "%s %s %hd", StrV1, StrV2, &Int2);
	HInfo->channel_count = Int2;

	/*---- sample_count -----*/
	Field = strstr(header, "sample_count");
	if (Field == NULL) return (0);

	sscanf(Field, "%s %s %ld", StrV1, StrV2, &Int4);
	HInfo->sample_count = Int4;

	/*---- sample_rate -----*/
	Field = strstr(header, "sample_rate");
	if (Field == NULL) return (0);

	sscanf(Field, "%s %s %ld", StrV1, StrV2, &Int4);
	HInfo->sample_rate = Int4;

	/*---- sample_n_bytes -----*/
	Field = strstr(header, "sample_n_bytes");
	if (Field == NULL) return(0);

	sscanf(Field, "%s %s %hd", StrV1, StrV2, &Int2);
	HInfo->sample_n_bytes = Int2;

	/*---- sample_byte_format -----*/
	Field = strstr(header, "sample_byte_format");
	if (Field != NULL) {
	    sscanf (Field, "%s %s %s", StrV1, StrV2, StrV3);
	    if ( (IsBigEndian() && strcmp(StrV3, "10")==0) || 
			(!IsBigEndian() && strcmp(StrV3, "01")==0)) {HInfo->SwapByte = 0;}
		else  {HInfo->SwapByte = 1;}
	}
	else HInfo->SwapByte = 0;  /* suppose native format */

     /*---- sample_coding -----*/
	Field = strstr(header, "sample_coding");
	if (Field != NULL) {

		sscanf(Field, "%s %s %s", StrV1, StrV2, StrV3);
		/*-------------------------------------------*/
		/* determine encoding string contents        */
		/* modify this part to recognize more types  */
		/*-------------------------------------------*/
		if (strcmp(StrV3, "pcm")==0)
			HInfo->sample_coding = PCM_2;
		else if (strstr(StrV3, "pcm") && strstr(StrV3, "shorten"))
			HInfo->sample_coding = PCM_2_Shorten;
	    else if (strcmp(StrV3, "ulaw") == 0 || strcmp(StrV3, "mu-law") == 0)
			HInfo->sample_coding = ULAW;
		else if (strstr(StrV3, "ulaw") && strstr(StrV3, "shorten"))
			HInfo->sample_coding = ULAW_Shorten;
		else
			HInfo->sample_coding = UNKNOWN;
	}
    else  {HInfo->sample_coding = PCM_2;}  /* default is two byte PCM*/

	return(1);  /* Have all need info */
}




