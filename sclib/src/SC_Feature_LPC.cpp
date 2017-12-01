//************************************************************************
//    Calculte LPC coefficients.
//
//    Author  : Jun Zhou
//    Date    : March 22, 2006
//************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SC_Feature_LPC.h"
#include "SC_Signal.h"
#include "SC_Aux.h"
#include "SC_Transform.h"
#include <GN_LPC.h>
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_LPC::SC_Feature_LPC(int sampleRate, int frameSize, int frameStep, unsigned int window, double preemphasis, unsigned short LPCorder, bool computeGain) : SV_Feature() {
	this->Para.StpSz = frameStep; //by thilo
	this->Para.WinSz = frameSize; //by thilo
	this->Para.Alpha = preemphasis; // by thilo
	this->Para.LPC_Order = LPCorder; //by thilo: changed the parameter from extra protected variable to this one inside the Para-struct
	this->Para.SRate = sampleRate; // by thilo
	this->Para.HammingWin = window; //by thilo
	this->computeGain = computeGain; //by thilo
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Feature_LPC::~SC_Feature_LPC() {

}

//====================================================================================================================
//	extract LPC coefficients per frame
//====================================================================================================================
SV_Data *SC_Feature_LPC::ExtractFeature(void) {
	SV_Data *DataSet = NULL;
	float	*Segment;
  float *Signal;
  long SigLen;
  SC_Transform transformer; //by thilo
  double *wnd, gain; //by thilo
	int frameCount; //by thilo
		
	if (IsSigLoaded()) { 
		Signal = GetSig();
		SigLen = GetLen();
	}	else {
		return (NULL);
	}

	frameCount = (int)(sclib::getRowCount(SigLen, Para.WinSz, Para.StpSz)); //(SigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1; //(SigLen - Para.WinSz) / Para.StpSz + 1; //by thilo: mathematically equivalent, but this way it looks like in the other classes...
	if (frameCount > 0) {
		// Preeamphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		DataSet = new SV_Data;
		DataSet->Row = frameCount;
		DataSet->Col = this->Para.LPC_Order + (this->computeGain==true?1:0); //by thilo
		DataSet->Alloc();

		if (DataSet==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		} 

		DataSet->Hdr.frameSize = this->Para.WinSz; //by thilo
		DataSet->Hdr.frameStep = this->Para.StpSz; //by thilo
		DataSet->Hdr.sampleRate = this->Para.SRate; //by thilo
		DataSet->Hdr.ID = sclib::featureLPC; //by thilo
		DataSet->Hdr.Signature[1] = (this->computeGain == true) ? 1 : 0; //by thilo: remember if the last column is the gain rather than a last LPC coefficient

		//--------------------------------------------------
		// Prepare run LPC analysis  
		//--------------------------------------------------
		GN_LPC  LPC_Eng;	// LPC analysis engine
		double  *LpcBuf;
		
		MArray_1D(LpcBuf, this->Para.LPC_Order, double, "LpcBuf");
		MArray_1D(Segment,Para.WinSz, float, "Segment");
		
		//const double Factor = 2.0 * 3.1415926 / (Para.WinSz -1.0); //by thilo: seems we need windows...
		//double windows_coeff;
		wnd = transformer.window(this->Para.WinSz, this->Para.HammingWin); //by thilo: this way, the window-type can be choosen
		
		//--------------------------------------------------
		// Extract features from each frame  
		//--------------------------------------------------
		for (int FrmCnt=0; FrmCnt<DataSet->Row; FrmCnt++) {

			//---------------------------------------------------
			// content of Signal copy to Segment
			//---------------------------------------------------
			for (int Col=0; Col<Para.WinSz; Col++) {
				Segment[Col] = (float)(Signal[FrmCnt*Para.StpSz + Col] * wnd[Col]); //by thilo: added multiplication with window
				//windows_coeff =  0.54 - 0.46*cos(Factor * Col);
  			//Segment[Col] *= (float)(windows_coeff);	//by thilo: seems we need windows...
			}

			//--------------------------------------
			// Calculte LPC coef for one segment
			//--------------------------------------
			gain = LPC_Eng.CalcLPC (Segment, Para.WinSz, LpcBuf, this->Para.LPC_Order); //by thilo: catch new gain return value
			
			//--------------------------------------
			// Copy to Matrix
			//--------------------------------------
			for (int Col=0; Col<this->Para.LPC_Order; Col++) {
				DataSet->Mat[FrmCnt][Col] = float(LpcBuf[Col]);					
			}
			if (this->computeGain == true) { //block by thilo
				DataSet->Mat[FrmCnt][this->Para.LPC_Order] = (float)(gain);
			}

		}	// end of for(FrmCnt)

		MFree_1D(Segment);
		MFree_1D(LpcBuf);
		MFree_1D(wnd); //by thilo
	}
	
	return(DataSet);	
}
