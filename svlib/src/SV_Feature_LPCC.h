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
//    Extract LPCC vector sequence from a given signal.
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Feature_LPCC_H__
#define __SV_Feature_LPCC_H__

#include "SV_General.h"
#include "SV_Feature.h"
#include "SV_Data.h"

//-------------------------------
// Constant for FType
//-------------------------------
#define		Coef_LPC		1		// get LPC coefficients
#define		Coef_LPCC		2		// DEFAULT get LPCC, 


class SV_Feature_LPCC : public SV_Feature {

private :
			

protected :

public :

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Feature_LPCC();
	virtual ~SV_Feature_LPCC();
	
	//------------------------------------- 
	// override base class method, 
	// return LPCC or LPC vector sequence
	//------------------------------------- 
	virtual SV_Data *ExtractFeature(void);

};   // class SV_Feature_LPCC


#endif   // SV_Feature_LPCC_H
