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
//    Calculate pitch contour.
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Feature_Pitch_H__
#define __SV_Feature_Pitch_H__

#include "SV_General.h"
#include "SV_Feature.h"


//===================================================================
//  This class is derived from SV_Feature to extract Pitch Period
//===================================================================
class SV_Feature_Pitch : public SV_Feature {

private :
	

protected :

	float Pitch_Detector (float *PulseVal, float *PulsePos, int PulseNum);
	void Fill_PitchMatrix(void);		

public :

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Feature_Pitch();
	virtual ~SV_Feature_Pitch();

	//------------------------------- 
	// return pitch contour
	//-------------------------------
	virtual SV_Data *ExtractFeature(void);
	

};   // class SV_Feature_Pitch











#endif // MFCC_H
