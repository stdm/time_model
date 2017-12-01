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
//    Design FIR filter with window method. Filtering a signal
//    with LOWPASS, HIGHPASS, BANDPASS, BANDSTOP  filter.
//
//    Author  : Jialong HE
//    Date    : March 21, 1999
//
//************************************************************************
#ifndef __GN_Filter_H__
#define __GN_Filter_H__

#include "SV_General.h"

//------------------------
// Filter type constants
//------------------------
#define LOW_PASS     1
#define HIGH_PASS    2
#define BAND_PASS    3
#define BAND_STOP    4

//===================================================================
//  This class provides FIR design and filtering methods
//  Design is based on window the simple method. 
//  
//===================================================================
class GN_Filter {

private :

protected:

	double *Coef;    // hold FIR filter's coef.
	int   Order;     // number of coef.
	int   FType;
	int   fir_dsgn(int Len, double CutLow, double CutHigh);
	void  remez(int numband, double bands[], double des[], double weight[], int numtaps = 127, int type = 1);

public :

	//------------------------------- 
	// Constructor/Destructor
	//-------------------------------
	GN_Filter();
	virtual ~GN_Filter();

	//------------------------------- 
	// Access to protected members
	//-------------------------------
	int  Get_FType(void)   {return(FType);};  // Filter Type
	int  Get_Order(void)   {return(Order);};  // No.of coef.
	double *Get_Coef(void) {return(Coef);};	  // coef. buf				

	//------------------------------------------------------ 
	// Design a FIR filter, suppose sampling rate is 2 Hz, 
	// Nyquist Freq is 1 Hz,  Cutoff<1.
	//------------------------------------------------------
	int  LowPass(double CutOff, int NCoef = 127);
	int  HighPass(double CutOff, int NCoef = 127);
	int  BandPass(double LowCut, double HighCut, int NCoef = 127);
	int  BandStop(double LowCut, double HighCut, int NCoef = 127);


	//----------------------------------------- 
	// Filtering Sig with current filter 
	//----------------------------------------- 
	void FIR_Filter(float *Sig, int SampleNum);

};   // class GN_Filter


#endif
