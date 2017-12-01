//########################################################################
//  
// A C++ class library for automatic speech recognition and 
// speaker recognition (identification and verification). 
//
// This class library is provided "as is" without any express 
// or implied warranty of any kind with respect to this software. 
// In particular the author shall not be liable for any direct, 
// indirect, special, incidental or consequential damages arising 
// in any way from use of the software.
//
//
// Author   : Jialong He,  Copyright (C), all rights reserved. 
// Date     : May, 1999. 
//
// Contact  : Jialong_He@bigfoot.com, Jialong_He@homemail.com
// Web Page : www.bigfoot.com/~Jialong_He
//########################################################################

//************************************************************************
//    Base class for loading speech signal into memory.
//    This class can load raw speech data.
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Signal_H__
#define __SV_Signal_H__

#include "SV_General.h"
#include "SV_DataIO.h" //by thilo to do machine- and os-independant io

//------------------------------- 
// Parameter Block for Signal
//------------------------------- 
typedef struct {
	int SRate;                     // Sampling Rate (samples per second)  
    int NChannel;				   // Number of channel (1, or 2)
	int Encode;                    // 1: short int, 2: ulaw 
	int StByte;					   // Starting Byte to read (for ignore Hdr) 
} SIG_PAR;

/*------------------------------------*/
/* Current recognized coding format:  */
/*------------------------------------*/
#define UNKNOWN          0
#define PCM_2            1
#define PCM_2_Shorten    2
#define ULAW             3
#define ULAW_Shorten     4

//===================================================================
//  SV_Signal class, base class 
//===================================================================
class SV_Signal {

private :
	long ClassSig;					// Class signature

protected :
	short *Buf_L;                  // point to Left Channel Buffer 
	short *Buf_R;                  // point to Right Channel Buffer
	int   Len;                     // Number of samples 
	SV_DataIO io;									 // by thilo to do machine- and os-independant io

public :
	//------------------------------- 
	// public data members
	//------------------------------- 
	SIG_PAR SigPar;

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Signal();
	virtual ~SV_Signal();

	//------------------------------------
	// public methods
	//------------------------------------
	unsigned char linear2ulaw(int sample);    // linear --> ulaw
	int ulaw2linear(unsigned char ulawbyte);  // ulaw --> linear
	unsigned char linear2alaw(int sample);		
	int alaw2linear(unsigned char alawbyte);
	unsigned char alaw2ulaw(unsigned char alawbyte);
	unsigned char ulaw2alaw(unsigned char ulawbyte);


  int  Valid (void);                        // test this class signature

  int  IsBigEndian (void);
  void SwapByteOrder (void);
	void EndPoint(short *Buf, int Len, int *Start, int *End);
  virtual long FileSize(char *FName);       // only work for disk file (not pipe)
	
	//------------------------------------
	//  Provide access to internal buffers  
	//------------------------------------
	long  GetLen(void)  {return (Len);}
	short *GetBuf_L(void) {return (Buf_L);}
	short *GetBuf_R(void) {return (Buf_R);}

	//------------------------------------
	//  Main method to be used 
	//------------------------------------
	virtual long LoadSignal(char *FName);
	virtual long SaveSignal(char *FName);

};   // class SV_Signal


#endif

