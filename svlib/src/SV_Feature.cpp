//************************************************************************7
//    Implementation Feature class
//
//
//    Author  : Jialong HE
//    Date    : April 28, 1999
//************************************************************************
#include <math.h>
#include <fstream>
#include <limits> //by thilo
#include "SV_Feature.h"
#include "SV_Data.h"
#include "SV_Signal.h"
#include "SV_Error.h"
#include "SV_General.h"
#include "GN_Filter.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Feature::SV_Feature() {

	Sig   = NULL;
	Len   = 0;
	ClassSig = 1995131;
}

//==========================================
// destructor
//==========================================
SV_Feature::~SV_Feature() {
	
	if (Sig != NULL) {
		MFree_1D(Sig);
		Len = 0;	
	}

	ClassSig = 0;
}

//=======================================================*/
//  Test if current class is valid        
//  return 1: yes, it is valid             
//  return 0: no, it is not valid, maybe deleted.
//=======================================================*/
int SV_Feature::Valid (void) {
	if (ClassSig == 1995131) {return (1);}
	else {return(0);}
}

//================================================
// copy signal from short to internal float
//
// TODO: May need set other parameters
//================================================
void SV_Feature::CopySignal(short *SigBuf, int BufLen) {

	//--------------------------------------------
	// Error check
	//--------------------------------------------
	if (SigBuf == NULL || BufLen == 0) {
       REPORT_ERROR(SVLIB_BadArg, "CopySignal");   
	}

	//-------------------------------------------------
	// if new signal has different length than before
	//-------------------------------------------------
	if (BufLen != Len) {
		MFree_1D(Sig);
		MArray_1D(Sig, BufLen, float, "CopySignal");
		Len = BufLen;
	}	

	for (int Cnt=0; Cnt<Len; Cnt++) {
		Sig[Cnt] = float(SigBuf[Cnt]);
	}
}

//================================================
// copy signal from SV_Signal to internal float
//
// TODO: May need set other parameters
//================================================
void SV_Feature::CopySignal(SV_Signal& SigObj) {
	long BufLen;
	short *SigBuf;

	BufLen	= SigObj.GetLen();
	SigBuf	= SigObj.GetBuf_L();

	Para.SRate = SigObj.SigPar.SRate;
  CopySignal(SigBuf, BufLen);
}

//================================================
// Test if signal is loaded into Feature object
//================================================
int SV_Feature::IsSigLoaded(void) { 

	if (Sig != NULL && Len > 0) {
		return (1);
	}
	else {return (0);}

}

//================================================
// Preemphasize (highpass) internal signal
//================================================
void SV_Feature::PreEmphasize(double Alpha) {
	for (int Cnt = Len-1; Cnt>0; Cnt--) {
		Sig[Cnt] = Sig[Cnt] - float(Alpha) * Sig[Cnt - 1];
	}
}

//================================================
// Multiplying signal with hamming window
//================================================
void SV_Feature::ApplyHammingWin(double *Sig, int Num) {

    double Factor = 2.0 * 3.1415926 / (Num -1.0);

	for (int Cnt = Num-1; Cnt>0; Cnt--) {
		Sig[Cnt] = Sig[Cnt] * (0.54 - 0.46*cos(Factor * Cnt));
	}

}
//================================================
// Bandpass filtering internal signal
//================================================
void SV_Feature::FIR_filtering(double LowCut, double HighCut) {

	GN_Filter FIR;

	if (LowCut>0 && LowCut<HighCut && HighCut < 1) {
		FIR.BandPass(LowCut, HighCut);		// design FIR filter
		FIR.FIR_Filter(Sig, Len);			// filtering data
	}

}

//================================================
// Remove silence parts based energy threshold
//================================================
long SV_Feature::RemoveSilence(double Ratio, int AnaSize) {

	double *Energy;
	double MinVal, MeanVal, Threshold;
	long NSeg, Cnt, SegNum, SaveInd, AnaInd;
	
	NSeg = Len / AnaSize + 2;
	MArray_1D(Energy, NSeg, double, "RemoveSilence");

	//----------------------------------
	// Calculate Energy Contour
	//----------------------------------
	SegNum   = 0; 
	AnaInd   = 0;
	while (AnaInd + AnaSize < Len ) {
		//----------------------------------
		// Calculate energy of this segment
		//----------------------------------
		Energy[SegNum] = 0.0;
		for (Cnt=0; Cnt<AnaSize; Cnt++) {
			Energy[SegNum] += fabs(Sig[AnaInd + Cnt]);
		}
		
		AnaInd  += AnaSize;
		SegNum++;		// next segment
	}

	//------------------------------------------
	// Find Min() and Mean() of Energy contour  
	//------------------------------------------
	MinVal = Energy[0]; MeanVal = 0.0;
	for (Cnt=0; Cnt<SegNum; Cnt++) {
		if 	 (Energy[Cnt] < MinVal) {
			MinVal = Energy[Cnt]; 
		}
		MeanVal += Energy[Cnt];
	}
	MeanVal /= SegNum;

	//-------------------------------------------------
	// Experienced Threshold value
	//-------------------------------------------------
	Threshold = MinVal + Ratio * (MeanVal - MinVal);

	//-------------------------------------------------
	// Copy segments that have energy higher than Thres
	//-------------------------------------------------
	SaveInd  = 0;                 /* starting save samples */
	for (int SegCnt=0; SegCnt<SegNum; SegCnt++) {

		//-----------------------------------
		// keep this segment
		//-----------------------------------
		if (Energy[SegCnt] > Threshold) {

			AnaInd = SegCnt*AnaSize;
			for (Cnt=0; Cnt<AnaSize; Cnt++) {
				Sig[SaveInd + Cnt] = Sig[AnaInd + Cnt];
			}
			SaveInd += AnaSize;
		}	
	}

	Len = SaveInd;  // new length of signal	
	MFree_1D(Energy);
	return(SaveInd);
}

//===================================================
// Concatenate two feature sets to form a longer one 
// return NULL if can not form a new SV_Data
//===================================================
SV_Data *SV_Feature::Concat(SV_Data &First, SV_Data &Second) {

	SV_Data *NewData = NULL;
	int Row, Col;

	if (First.Row != 0 && First.Row == Second.Row) {
		NewData = new SV_Data;
		NewData->Hdr = First.Hdr;  // take first record's Header
		NewData->Row = First.Row; 
		NewData->Col = First.Col + Second.Col; 
		NewData->Alloc();  // allocate memory for Mat
		//---------------------------
		// Concatenate two data sets
		//---------------------------
		for (Row=0; Row<NewData->Row; Row++) {
			for (Col=0; Col<First.Col; Col++) {
				NewData->Mat[Row][Col] = First.Mat[Row][Col]; 
			}

			for (Col=0; Col<Second.Col; Col++) {
				NewData->Mat[Row][First.Col + Col] = Second.Mat[Row][Col]; 
			}
		}

	}  // if 
	return (NewData);
}	

//===================================================
// Remove mean of cepstral sequence 
//===================================================
void  SV_Feature::CMN(SV_Data *DataSet) {

	double *Mean;
	int Row, Col;

	MArray_1D(Mean, DataSet->Col, double, "CMN: Mean");
	
	//-----------------------------------------
	// Calculate cepstral mean
	//-----------------------------------------
	for (Col=0; Col<DataSet->Col; Col++) {
		Mean[Col] = 0;
		for (Row=0; Row<DataSet->Row; Row++) {
			Mean[Col] += DataSet->Mat[Row][Col];
		}	
		Mean[Col] /= DataSet->Row;
	}

	//-----------------------------------------
	// Substract mean from DataSet->Mat
	//-----------------------------------------
	for (Row=0; Row<DataSet->Row; Row++) {
		for (Col=0; Col<DataSet->Col; Col++) {
			DataSet->Mat[Row][Col] -= float(Mean[Col]);
		}	
	}

	MFree_1D(Mean);
}

//===================================================
// Provide only a base function, override by derived  
// classes such as MFCC, LPCC etc.
//===================================================
SV_Data *SV_Feature::ExtractFeature(void) {
	return (NULL);
}



//===================================================
// Calculate dynamic (delta) feature.  
// if input is delta feature, it gives delta-delta feature
//
//
//  See Rabinar's book Eq. 3.91
//===================================================
SV_Data *SV_Feature::Dynamic(SV_Data &DataSet) {

	SV_Data *NewData = NULL;
	int Row, Col, SubCnt;

	if (DataSet.Row>0 && DataSet.Col>0) {
		//---------------------------------------------------
		// Allocate memory for NewData
		//---------------------------------------------------
		NewData = new SV_Data;
		if (NewData == NULL) {
			REPORT_ERROR(SVLIB_NoMem, "Dynamic");
		}

		//---------------------------------------------------
		// NewData has the same size as DataSet
		//---------------------------------------------------
		NewData->Row = DataSet.Row;
		NewData->Col = DataSet.Col;
		NewData->Alloc();

		//--------------------------------------
		// For 1st row, simply make a copy
		//--------------------------------------
		for (Col=0; Col<DataSet.Col; Col++) {
			NewData->Mat[0][Col] = DataSet.Mat[0][Col];	
		}
		
		//---------------------------------------------------
		// for other row, calculate linear regression
		//---------------------------------------------------
		for (Row=1; Row<DataSet.Row; Row++) {
			for (Col=0; Col<DataSet.Col; Col++) {

				NewData->Mat[Row][Col] = 0.0;
				if (Row>=3 && Row<DataSet.Row-3) {
					//------------------------------
					// use frame (-3, 3), 7 frames 
					//------------------------------
					for (SubCnt = -3; SubCnt <= 3; SubCnt++) {
						NewData->Mat[Row][Col] += SubCnt * DataSet.Mat[Row + SubCnt ][Col];
					}
				}

				else {
					//------------------------------
					// use two frames -1, 0 
					//------------------------------
					NewData->Mat[Row][Col] = DataSet.Mat[Row][Col] - DataSet.Mat[Row -1][Col];
				}		
			}  // for Col
		}  // for Row
	}	// if (DataSet.Row>0 && DataSet.Col>0)

	return (NewData);
}

//====================================================================================================================
//	by thilo: manipulate internal signal pointer to link to external memory (it can be selected if the memory 
//  currently pointed at should be freed before setting to a new location)
//====================================================================================================================
void SV_Feature::setSignal(float* signal, unsigned long int len, bool freeOldSignal) {

	if (freeOldSignal == true) {
		MFree_1D(this->Sig);
	}

	this->Sig = signal;
	this->Len = len;

	return;
}
