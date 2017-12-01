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
//    This is the base class of all feature extraction classes.
//    It contains commonly used functions.
//    Specific features can be derived from this class (e.g. MFCC, LPCC)
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Feature_H__
#define __SV_Feature_H__

#include "SV_General.h"
#include "SV_Signal.h"
#include "SV_Data.h"

//-------------------------------- 
// Analysis parameters 
//-------------------------------- 
typedef struct {
	//======================
	// LPCC related
	//======================
	int		WinSz;				// default 256 samples
	int		StpSz;				// default 128 samples
	long	FType;				// feature sets, bit ORed;
	int		LPC_Order;			// LPC order
	int		LPCC_Order;			// number of LPCC
	double	Alpha;				// preemphesis, default=0.97
	int 	HammingWin;			// 1: Hamming window, 0: rect window (no window)
	int		RmvSilence;			// 1: remove silence part, 0: intact

	//======================
	// extra parameter of MFCC 
	//======================
	int		SRate;	            // sampling rate of Sig;
	int		MFCC_Order;			// number of MFCC
	int		NFilter;			// number of filters in Filter bank 
	int		FFTSz;				// FFT length 
	int		DEnergy;			// 1: output delta value for first MFCC (energy)

	//======================
	// new para for Pitch 
	//======================
	int		Smooth;				// 1: smooth pitch contour with 3-point median filter

} PARA_TYPE;

//===================================================================
//  Base class of feature extraction algorithms.
//===================================================================
class SV_Feature {

private :
	long ClassSig;

protected :

	float *Sig;   // signal to be analysized
	long  Len;    // Length of *Sig

public :
	//------------------------------- 
	// parameters for analysis method
	//------------------------------- 
	PARA_TYPE Para;

	//------------------------------- 
	// Constructor/Destructor 
	//------------------------------- 
	SV_Feature();
	virtual ~SV_Feature();

	//------------------------------- 
	// Load signal into *Sig
	//------------------------------- 
	void CopySignal(SV_Signal& SigObj);
	void CopySignal(short *SigBuf, int BufLen);

  //====================================================================================================================
	//	by thilo: manipulate internal signal pointer to link to external memory (it can be selected  if the memory 
	//  currently pointed at should be freed before setting to a new location)
  //====================================================================================================================
	void setSignal(float* signal, unsigned long int len, bool freeOldSignal = true);

	//------------------------------------------------ 
	// return 1,  if Sig loaded, otherwise return 0 
	//------------------------------- ----------------
	int  IsSigLoaded(void);
	int  Valid (void);

	//---------------------------
	// access to protected members
	//---------------------------
	float *GetSig(void) {return(Sig);}
	long  GetLen(void) {return(Len);}
	
	//-------------------------------------------------------- 
	//  Horizontally concate two patterns to form 
	//  a wider feature. Internal allocate memory for new data
	//-------------------------------------------------------- 
	SV_Data *Concat(SV_Data &First, SV_Data &Second);

	//-------------------------------------------------------- 
	//  Obtain dynamic feaure sets (Delta), 
	//  if applying two times, obtaining Delta-Delta feature
	//  Internal allocate memory for new data
	//-------------------------------------------------------- 
	SV_Data *Dynamic(SV_Data &DataSet);

	//-------------------------------------------------------- 
	//  Cepstral mean normalization, 
	//  retrun mean vector of cepstral sequence
	//-------------------------------------------------------- 
	void CMN(SV_Data *DataSet);
	
	//-------------------------------------------------------- 
	//  Preprocessing *Sig
	//-------------------------------------------------------- 
	void ApplyHammingWin(double *Sig, int Num);
	void PreEmphasize(double Alpha = 0.97);
	void FIR_filtering(double LowCut, double HighCut); // normlized Nyquest Freq 0-1 Hz
	long RemoveSilence(double Ratio=0.3, int AnaSize = 32);

	//-------------------------------------------------------- 
	//  This method must be overrided by derived class 
	//  to extract different features.
	//-------------------------------------------------------- 
	virtual SV_Data *ExtractFeature(void);

};   // class SV_Feature


#endif   // SV_Feature
