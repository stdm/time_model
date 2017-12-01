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
//    Generate typical waveforms.
//
//
//    Author  : Jialong HE
//    Date    : March 24, 1999
//************************************************************************
#ifndef __GN_WaveGen_H__
#define __GN_WaveGen_H__

#include "SV_General.h"


//===================================================================
//  This class is a waveform generator.  
//==============================s====================================
class GN_WaveGen {

private :

protected:

	double *Wave;        // hold wave 
	int   Length;        // *Wave length

	void   AlloWaveBuf(int Len);

public :

	//------------------------------- 
	// Constructor / Destructor  
	//-------------------------------
	GN_WaveGen();
	virtual ~GN_WaveGen();
	
	//------------------------------- 
	// access to protected members  
	//-------------------------------
	int     GetLength(void) {return (Length);};
	double* GetWave(void) {return (Wave);};

	//----------------------------------- 
	// Typical waveform generators
	// return pointer of buffer (*Wave)
	//-----------------------------------
	double *Square_Wave(int HighDur, int LowDur, int Len);
	double *Sawtooth_Wave(int RampDur, int LowDur, int Len);
	double *Tone_Wave(double Cycle, double Phase, int Len);
	double *NoiseUni_Wave(int Len);
	double *NoiseGaus_Wave(double Mean, double Std, int Len);
	double *ModTone_Wave(double ModCyl, double CarCyl, int Len);

};   // class GN_WaveGen


#endif
