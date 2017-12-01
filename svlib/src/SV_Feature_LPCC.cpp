//************************************************************************
//    
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <math.h>
#include <stdlib.h>
#include "GN_LPC.h"
#include "SV_Feature_LPCC.h"
#include "SV_Error.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Feature_LPCC::SV_Feature_LPCC() : SV_Feature(){

	Para.WinSz = 256;				// default 256 samples
	Para.StpSz = 128;				// default 128 samples
	Para.FType = Coef_LPCC;				// feature sets, bit ORed;
	Para.LPC_Order  = 14;			// LPC order
	Para.LPCC_Order = 16;			// number of LPCC, D_LPCC, DD_LPCC
	Para.Alpha = 0.97;				// preemphesis, default=0.97
	Para.HammingWin = 1;			// Applying Hamming window
	Para.RmvSilence = 0;			// remove silence part, def. Not
}

//==========================================
// default desstructor
//==========================================
SV_Feature_LPCC::~SV_Feature_LPCC() {


}

//==========================================
// This is the engine of deriving features
//==========================================
SV_Data *SV_Feature_LPCC::ExtractFeature(void) {

	SV_Data *DataSet = NULL;
	float	*SigBuf, *Segment;
	long	SigLen;
	int		Col;

	//--------------------------------
	// Check if Signal is loaded
	//--------------------------------
	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}
	else return (NULL);	

	DataSet = new SV_Data;
	if (DataSet==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}

	//--------------------------------------------------
	// Processing Signal in Sig  
	//--------------------------------------------------
	if (Para.RmvSilence) {
		RemoveSilence(0.3, 32);
		SigLen = GetLen();  // new SigLen
	}

	PreEmphasize(Para.Alpha);
	//--------------------------------------------------
	// find out how many Row and Col for feature vectors  
	//--------------------------------------------------
	DataSet->Row = (SigLen - Para.WinSz + Para.StpSz) / Para.StpSz;
	if (Para.FType == Coef_LPC) {
		DataSet->Col  = Para.LPC_Order;
	}
	else { // (Para.FType == LPCC) 
		DataSet->Col  = Para.LPCC_Order;
	}

	DataSet->Alloc();	// allocate memory for Mat

	//--------------------------------------------------
	// Prepare run LPC analysis  
	//--------------------------------------------------
	GN_LPC  LPC_Eng;	// LPC analysis engine
	double  *LpcBuf, *CepBuf;
	
	MArray_1D(LpcBuf, Para.LPC_Order, double, "LpcBuf");
	MArray_1D(CepBuf, Para.LPCC_Order, double, "CepBuf");
	MArray_1D(Segment, Para.WinSz, float, "Segment");

    double Factor = 2.0 * 3.1415926 / (Para.WinSz -1.0);
	//--------------------------------------------------
	// Extract features from each frame  
	//--------------------------------------------------
	for (int FrmCnt=0; FrmCnt<DataSet->Row; FrmCnt++) {

		//---------------------------------------------------
		// Not destroy content in SigBuf, make a copy
		//---------------------------------------------------
		for (Col=0; Col<Para.WinSz; Col++) {
			Segment[Col] = SigBuf[FrmCnt*Para.StpSz + Col];
   		    if (Para.HammingWin) {
				Segment[Col] = float(Segment[Col] * (0.54 - 0.46*cos(Factor * Col)));		
			}

		}


		//--------------------------------------
		// Calculte LPC coef for one segment
		//--------------------------------------
		LPC_Eng.CalcLPC (Segment, Para.WinSz, LpcBuf, Para.LPC_Order);
		
		//--------------------------------------
		// Copy to DataSet
		//--------------------------------------
		
		switch (Para.FType) {
			case Coef_LPC:
				for (Col=0; Col<DataSet->Col; Col++) {
					DataSet->Mat[FrmCnt][Col] = float(LpcBuf[Col]);					
				}
				break;
			
			case Coef_LPCC:
				LPC_Eng.Lpc2Cep (LpcBuf, Para.LPC_Order, CepBuf, Para.LPCC_Order);
				for (Col=0; Col<DataSet->Col; Col++) {
					DataSet->Mat[FrmCnt][Col] = float(CepBuf[Col]);					
				}
				break;

			default: REPORT_ERROR(SVLIB_BadArg, "ExtractFeature");

		}  // switch


	}	// for (FrmCnt)

	MFree_1D(Segment);
	MFree_1D(LpcBuf);
	MFree_1D(CepBuf);
	return(DataSet);
}


