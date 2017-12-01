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
//    Extract MFCC vector sequence from a given signal.
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Feature_MFCC_H__
#define __SV_Feature_MFCC_H__

#include "SV_General.h"
#include "SV_Feature.h"
#include "SV_Data.h"

//===================================================================
//  This class is derived from SV_Feature to extract MFCC from 
//  a whole utterance. 
//===================================================================
class SV_Feature_MFCC : public SV_Feature {

private :

protected :
	int *Mel_Index(void);  // Mel-Scale index for samples

public :

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Feature_MFCC();
	virtual ~SV_Feature_MFCC();

	//-------------------------------
	// override base class method
	// return MFCC vector sequence
	//------------------------------- 
	virtual SV_Data *ExtractFeature(void);

};   // class SV_Signal

#endif // SV_Featuer_MFCC_H
