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
//    This class calculates radix-2 FFT/IFFT, DFT/IDFT, DCT/IDCT  
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//
//************************************************************************
#ifndef __GN_FFT_H__
#define __GN_FFT_H__

#include "SV_General.h"


//===================================================================
//  This class provides fast Fourier transform (FFT, IFFT)
//  discrete Fourier transform (any length) (DFT, IDFT)
//  cosine transform (DCT, IDCT) and real input FFT
//   
//===================================================================
class GN_FFT {

private :

public :

	//------------------------------- 
	// Constructor/Destructor
	//------------------------------- 
	GN_FFT ();
	virtual ~GN_FFT ();

	//------------------------------- 
	// FFT/IFFT, FFT for real number
	// Len must be 2**m, e.g., 1024, 2048
	//------------------------------- 
    void fft(COMPLEX *In_Out, int Len);		
    void ifft(COMPLEX *In_Out, int Len);
    void rfft(double *In, COMPLEX *Out, int Len); // *Out is Len/2  

	//------------------------------- 
	// Fast Discrete Cosine Transform
	// Len must be 2**m,
	//------------------------------- 
	void dct(double *In_Out, int Len);  
	void idct(double *In_Out, int Len);  
	
	//-------------------------------------------
	// Discrete Fourier Transform of any length
	//-------------------------------------------
    void dft(COMPLEX *In, COMPLEX *Out, int Len);
    void idft(COMPLEX *In, COMPLEX *Out, int Len);
	
	//-------------------------------------------
	// Mag_angle(): Convert to Polar (Mag/Phase)
	// Real_Imag(): Convert to Real/Imag 
	//-------------------------------------------
    void Mag_Angle(COMPLEX *In_Out, int Len);      
    void Real_Imag(COMPLEX *In_Out, int Len);

};   // class GN_FFT

#endif   // GN_FFT
