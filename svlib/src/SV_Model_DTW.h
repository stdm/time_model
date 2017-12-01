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
//    Calculate distance between two patterns by
//    Dynamic Time Wrapping (DTW)
//
//
//    Author  : Jialong HE
//    Date    : March 21, 1999
//************************************************************************
#ifndef __SV_Model_DTW_H__
#define __SV_Model_DTW_H__

#include <iostream>
#include "SV_General.h"
#include "SV_Model.h"
#include "SV_Data.h"

//===================================================================
//  DTW class
//===================================================================
class SV_Model_DTW : public SV_Model {

private :

protected :
	SV_Data *RefData;    // linked list of reference patterns
	int NumRef;          // number of reference patterns

public :

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Model_DTW();
	virtual ~SV_Model_DTW();

	//------------------------------- 
	// public methods
	//-------------------------------
	double DTW_Comp(SV_Data &First, SV_Data& Second);
	double dtw_dist(float **TestFrame, int TestLen, float **RefFrame, int RefLen, int Dim);
	virtual int TrainModel(SV_Data *TrainData); //by thilo: changed return-value from void to int (needed in derived class to indicate error)

};

#endif
