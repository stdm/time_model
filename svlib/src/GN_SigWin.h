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
//    Signal window functions. 
//
//
//    Author  : Jialong HE
//    Date    : March 24, 1999
//************************************************************************
#ifndef __GN_SigWin_H__
#define __GN_SigWin_H__

#include "SV_General.h"
	
#define  HAMMING	1
#define  HANNING	2
#define  TRIANGLE	3
#define  RECT		4
#define  BLACKMAN	5
#define  HARRIS		6


//===================================================================
//  This class provides commonly used signal windows.
//  call each method, return butter contains desired window function. 
//===================================================================
class GN_SigWin {

protected:

	double *SigWin;         // hold window function
	int   Length;           // SigWin length 
	int   WinType;          // window type

	void Alloc(int Len);

public :

	//------------------------------- 
	// Constructor/Destructor
	//------------------------------- 
	GN_SigWin ();
	virtual ~GN_SigWin ();

	//------------------------------- 
	// access to protected members
	//------------------------------- 
	int   GetLength(void) {return (Length);};
	double* GetSigWin(void) {return (SigWin);};

	//--------------------------------------------- 
	// design window functions, 
	// reture the pointer to the window function. 
	//---------------------------------------------
	double *Hamming(int Len);
	double *Hanning(int Len);
	double *Triangle(int Len);
	double *Blackman(int Len);
	double *Harris(int Len);
	
	//------------------------------------------------------- 
	// Multiplying Sig with current window functions
	// return a pointer to the window function.
	//------------------------------------------------------- 
	double *ApplyWindow(double *Sig, int Len, int WinType = HAMMING);

};   // class GN_Window

#endif   // GN_Window
