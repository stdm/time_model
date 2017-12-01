//************************************************************************
//    
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <math.h>
#include "SV_Feature_MFCC.h"
#include "GN_FFT.h"
#include "SV_Error.h"
#define Norm 10                 /* scale factor */

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Feature_MFCC::SV_Feature_MFCC() : SV_Feature(){

	Para.WinSz		= 256;			// default 256 samples
	Para.StpSz		= 128;			// default 128 samples
	Para.MFCC_Order = 16;			// number of MFCC, D_MFCC, DD_MFCC
	Para.NFilter    = 60;			// Number of filters in filter bank
	Para.FFTSz	    = 1024;			// FFT size
	Para.Alpha		= 0.97;			// preemphesis, default=0.97
	Para.HammingWin = 1;			// Applying Hamming window
	Para.RmvSilence = 0;			// remove silence part, def. Not
	Para.DEnergy	= 1;			// use Delta engergy
	
}

//==========================================
// default desstructor
//==========================================
SV_Feature_MFCC::~SV_Feature_MFCC() {


}

//================================================================
// Mel-Scale defination is
//  
//   mel = 2595*log10(1 + f/700); 
//
//  or f = 700*(exp( mel/2595*log(10) ) - 1);
//    
//  n filter's sample range is IndArray[n] -> IndArray[n+1]
// 
//================================================================
int *SV_Feature_MFCC::Mel_Index(void) {  // Mel-Scale index for samples
	
	int		*IndArray;
	double MelMax, MelStp, Freq;
	
	MArray_1D(IndArray, Para.NFilter+1, int, "Index Array");
	MelMax = 2595*log10(1 + Para.SRate / 2.0/700.0);  // half sampling rate
	MelStp = MelMax / Para.NFilter;  // Mel-Freq for each filter 

	for (int Cnt=0; Cnt<=Para.NFilter; Cnt++) {
		Freq = 700*(exp( Cnt*MelStp / 2595*log(10.0) ) - 1.0);
		//--------------------------------------------
		// when Cnt=Para.NFilter, Freq = SRate/2
		// IndArray is Para.FFTSz/2
		//--------------------------------------------
		IndArray[Cnt] = int(Freq * Para.FFTSz / Para.SRate);
	}

	return(IndArray);
}

//==========================================
// This is the engine of deriving features
//==========================================
SV_Data *SV_Feature_MFCC::ExtractFeature(void) {

	SV_Data *DataSet = NULL;
	GN_FFT	Trans;
	COMPLEX *Spec;

	float	*SigBuf;
	double  *Segment;
	double  SUM, *LogSpec, **Coef;
	long	SigLen;
	int		*MelInd, Col, Cnt, Start, End, FrmCnt;
	int RowCnt, ColCnt;
	
	
	//--------------------------------
	// Check if Signal is loaded
	//--------------------------------
	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}
	else return (NULL);	
	//--------------------------------------------------
	// Processing Signal in Sig  
	//--------------------------------------------------
	if (Para.RmvSilence) {
		RemoveSilence(0.3, 32);
		SigLen = GetLen();  // new SigLen
	}

	//--------------------------------
	// Allocate memory for DataSet
	//--------------------------------
	DataSet = new SV_Data;
	if (DataSet==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}
	DataSet->Row = (SigLen - Para.WinSz + Para.StpSz) / Para.StpSz;
	DataSet->Col = Para.MFCC_Order;
	DataSet->Alloc();

	//-----------------------------------
	// check FFT size
	//-----------------------------------
	if (Para.WinSz > Para.FFTSz) {
		REPORT_ERROR(SVLIB_BadArg, "FFT length should be larger than WinSz");
	}

	MelInd = Mel_Index();		// get mel-scale index

	MArray_1D(Segment, Para.FFTSz, double, "MFCC:Segment");
	MArray_1D(Spec, (Para.FFTSz/2), COMPLEX, "MFCC:Spec");
	MArray_1D(LogSpec, (Para.NFilter), double, "MFCC:LogSpec");

	/*-----------------------------------------------------*/
	/* Calculate cosine transform coefficients 
	/*-----------------------------------------------------*/
	MArray_2D(Coef, Para.MFCC_Order, Para.NFilter, double, "Cosine Coeff");
    for (RowCnt=0; RowCnt<Para.MFCC_Order; RowCnt++) {
		for (ColCnt=0; ColCnt<Para.NFilter; ColCnt++) {
			Coef[RowCnt][ColCnt] = cos((RowCnt) * (ColCnt-0.5) * 3.141592653589793 / Para.NFilter);
		}
	}

	//--------------------------------
	// Calculate MFCC for each frame
	//--------------------------------
	for (FrmCnt=0; FrmCnt<DataSet->Row; FrmCnt++) {

		//---------------------------------------------------
		// Not destroy content in SigBuf, make a copy
		//---------------------------------------------------
		for (Col=0; Col<Para.WinSz; Col++) {
			Segment[Col] = SigBuf[FrmCnt*Para.StpSz + Col];
		}

		for (Col=Para.WinSz; Col<Para.FFTSz; Col++) {
			Segment[Col] = 0.0;  // append zeros
		}
	     
		if (Para.HammingWin) {
			ApplyHammingWin(Segment, Para.WinSz);
		}

		Trans.rfft(Segment, Spec, Para.FFTSz);   // FFT for real signal 

		//-----------------------------
		// set DC to zero
		//-----------------------------
		Spec[0].real = 0.0;
		Spec[0].imag = 0.0;

		//--------------------------------
		//  Sum spec within a filter bank
		//-------------------------------
		for (Cnt=0; Cnt<Para.NFilter; Cnt++) {
			SUM	  = 0.0;
			Start = MelInd[Cnt];
			End   = MelInd[Cnt+1];
			for (Col=Start; Col<End; Col++) {
				SUM +=  sqrt(Spec[Col].real * Spec[Col].real +
							Spec[Col].imag * Spec[Col].imag);
			}

			if ((Start != End) && (SUM!=0.0)) {
				SUM /= End - Start;          /* Average magnitude */
				LogSpec[Cnt] = log10(SUM);   // log spectrum 
			}
			else {LogSpec[Cnt] = 0.0;}

		}   /* for Cnt */
  
		/*----------------------------*/
		/*  do cosine transform       */
		/*----------------------------*/
		for (Cnt=0; Cnt<Para.MFCC_Order; Cnt++) {
			SUM = 0.0;
			for (Col=0; Col<Para.NFilter; Col++) {
				SUM += LogSpec[Col] * Coef[Cnt][Col];
			}
			DataSet->Mat[FrmCnt][Cnt] = float(SUM / Norm);
		} /* for Cnt */	


	}	//	for (FrmCnt)

	//----------------------------------------------------------
	// first MFCC is engery, normally absolute energy is no use
	// but delta energy useful
	//----------------------------------------------------------
	if (Para.DEnergy) {
		
		for (FrmCnt=DataSet->Row-1; FrmCnt>=1; FrmCnt--) {
			DataSet->Mat[FrmCnt][0] = DataSet->Mat[FrmCnt][0] - DataSet->Mat[FrmCnt-1][0];
		}
	}

	MFree_1D(Segment);
	MFree_1D(Spec);
	MFree_1D(LogSpec);
	MFree_1D(MelInd);
	MFree_2D(Coef);
	return(DataSet);
}

