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
//    Derived from SV_Signal to load data with SPHRE header (NIST format).
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Signal_NIST_H__
#define __SV_Signal_NIST_H__

#include "SV_General.h"
#include "SV_Signal.h"

#define HLen 1024      /* NIST header length */

/*----------------------------*/
/* Important Header Info      */
/*----------------------------*/
typedef struct {
   short channel_count;
   short sample_n_bytes;
   long  sample_count;
   long  sample_rate;
   short SwapByte;        /* 1: machine's byte order is different from data. */
   short sample_coding;

} HeaderType;


//===================================================================
//  This is the base class of loading speech file into memory.
//  It can be used to derive a class to load signal in specific format 
//  such as NIST format.
//===================================================================
class SV_Signal_NIST : public SV_Signal {


private :

protected :
	int SPHERE_HeaderInfo (char *FName, HeaderType *HInfo);
	int LoadNIST(char *FName, short *Speech, HeaderType &HInfo);

public :

	//------------------------------- 
	// Constructor / Destructor
	//------------------------------- 
	SV_Signal_NIST();
	virtual ~SV_Signal_NIST();

	//------------------------------- 
	// public methods
	//------------------------------- 
	virtual long LoadSignal(char *FName);

};


#endif
