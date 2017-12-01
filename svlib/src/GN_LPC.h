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
//    LPC analysis using autocorrelation method.
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//
//************************************************************************
#ifndef __GN_LPC_H__
#define __GN_LPC_H__

//===================================================================
//  Interface definiation for LPC class
//===================================================================
class GN_LPC {

private :

protected :

	int LPC_Order;
    double AutoCorrelation(float *Signal, int Samples, double *AutoCorr, int Order); //by thilo: changed return type from void to double to return original CorCoeff[0] so to be able to reverse the normalization applied here
    void DUrbin(double *CorrCoef, double *LPCArray, double *Reflect);

public :
	//-------------------------------
	// Constructor/Destructor
	//-------------------------------
	GN_LPC();
	virtual ~GN_LPC();
	
	//-------------------------------------------------------------
	// Calculate LPC coefficients
	// ***NOTE*** LPC coefficients should be (1, a1, a2, a3, ...)
	// but the first "1" is not returned in LpcCoef.
	//--------------------------------------------------------------
	double CalcLPC (float* Sig, int Len, double* LpcCoef, int LpcOrder); //by thilo: changed return value from void to double to return the LPC filter's gain

	//---------------------------------------------------
	// Convert LPC coef to Cepstral coefficients
	//---------------------------------------------------
	void Lpc2Cep (double* LpcCoef, int LpcOrder, double* CepCoef, int CepOrder);
	
	//==============================================================
	//  line Spectrum Frequency <-> LPC coefficients
	//  ***NOTE*** LsfCoef has (LpcOrder + 1)
	//  LsfCoef should be (0, b1, b2, b3, ...), where bi are angles
	//  divided by Pi, but the first zero (b0) is not returned/used.
	//==============================================================
    void Lpc2Lsf (double* LpcCoef, double *LsfCoef, int LpcOrder);
	void Lsf2Lpc (double* LsfCoef, double *LpcCoef, int LpcOrder);

	//----------------------------------------------------
	// residual signal (inverse filtering signal)
	//----------------------------------------------------
	void Residual(float* Sig, float* Rst, int Len, double* LpcCoef, int LpcOrder);

};   // class GN_LPC


#endif



