/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for extracting short time Energy							            */
/*			- Copy from SC_Feature_EnergyZCR										              */
/*    Author  : Bing Shi													                        */
/*    Date    : 11.04.2006																								*/
/**************************************************************************/

#include <math.h>
#include "SC_Aux.h"
#include "SC_Feature_STE.h"
#include "SC_TweakableParameters.h"
#include <GN_FFT.h>
#include <SV_Error.h>

//==========================================
// constructor
//==========================================
SC_Feature_STE::SC_Feature_STE(int sampleRate, int frameSize, int frameStep, double preEmphasize, bool useButterworth, bool scaleResult) : SV_Feature(){
	this->Para.WinSz = frameSize;
	this->Para.StpSz = frameStep;
  this->Para.Alpha = preEmphasize;
	this->Para.SRate = sampleRate;
	this->useButterworth = useButterworth;
	this->scaleResult = scaleResult;
}

//==========================================
// default destructor
//==========================================
SC_Feature_STE::~SC_Feature_STE() {

}

//==========================================
// Calculate Energy and Zero Crossing Rate  
//==========================================
SV_Data* SC_Feature_STE::ExtractFeature(void) {
	SV_Data *pData = NULL;
	double CurrE; 
  short ThisSample;
	int FrameCnt, SampleCnt, sampleNr;
	float *SigBuf;
	long SigLen;
	double scaleFactor = (this->Para.WinSz * 1000.0) / (double)(this->Para.SRate); //needed to convert frame-energy to energy per ms
	int frameCount;

	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}	else {
		return (NULL);
	}

	frameCount = (int)(sclib::getRowCount(SigLen, Para.WinSz, Para.StpSz)); //(SigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1; //(SigLen - Para.WinSz + Para.StpSz) / Para.StpSz;
	if (frameCount > 0) {
		//apply special hardcoded butterworth filter, if wished
		if (this->useButterworth == true) {
			if (this->Para.SRate = 8000) {
				butterworthBandPass(SigBuf, SigLen);
			} else {
				REPORT_ERROR(SVLIB_BadArg, "The hardcoded Butterworth bandpass can only be used with 8kHz signals");
				FIR_filtering((200.0/(double)(this->Para.SRate))*2.0, (3400.0/(double)(this->Para.SRate))*2.0);
			}
		}

		// Preeamphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}
		pData->Row = frameCount;
		pData->Col = 1;
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureSTE;
		pData->Alloc();

		for (FrameCnt = 0; FrameCnt < pData->Row; FrameCnt++) {
			CurrE = 0.0;
	   
			sampleNr = FrameCnt * this->Para.StpSz;
			for (SampleCnt = 1; SampleCnt < this->Para.WinSz; SampleCnt++) {
				ThisSample = (short)SigBuf[sampleNr + SampleCnt];
				
				CurrE += (double)(ThisSample*ThisSample);  // changed from abs(thisSample) to square(thisSample) on 21.02.2006 by thilo	
			}

			if (this->scaleResult == true) {
				CurrE /= scaleFactor; //per-frame => per-ms
			}
			pData->Mat[FrameCnt][0] = (CurrE > 0.0) ? (float)(10.0*log10(CurrE)) : (float)(10.0*-87.3);  //(float)CurrE; //amplitude in db according to http://en.wikipedia.org/wiki/Decibel
		}  // for FrameCnt
	}
	 
	return pData;
}

//====================================================================================================================
//	A hardcoded implementation of the Butterworth bandpass filter as specified in "Cepstrum-Based Pitch Detection 
//  Using a New Statistical V/UV Classification Algorithm" (Ahmadi, Spanias 1999): 9th order, 8kHz sample-rate, pass-
//  band 200-3400Hz; the filter was generated using http://www-users.cs.york.ac.uk/~fisher/mkfilter 
//====================================================================================================================
void SC_Feature_STE::butterworthBandPass(float* signal, long int length) {
  int i; 
	long int t;
	const double GAIN = 6.354485727e+00;
	double xv[19], yv[19]; //in the originally created code, xv has size NZEROS+1 and yv NPOLES+1 with NZEROS=NPOLES=18

	//Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
	//Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 9 -a 2.5000000000e-02 4.2500000000e-01 -l
	//initializations are done by the compiler via the static keyword in the originally generated code...
	for (i = 0; i < 18+1; i++) {
		xv[i] = 0.0;
		yv[i] = 0.0;
	}

	for (t = 0; t < length; t++) {
		for (i = 0; i < 18; i++) {
			xv[i] = xv[i+1];
			yv[i] = yv[i+1];
		}

		xv[18] = signal[t] / GAIN; //originally, this line came beofer the initialisation of yv[0]-yv[17], but this doesn't matter 'case they are both independent of xv[18]
		
		yv[18] =   (xv[18] - xv[0]) + 9 * (xv[2] - xv[16]) + 36 * (xv[14] - xv[4])
								+ 84 * (xv[6] - xv[12]) + 126 * (xv[10] - xv[8])
								+ (  0.0247650485 * yv[0])  + (  0.0453054927 * yv[1])
								+ ( -0.2694106370 * yv[2])  + ( -0.5045292753 * yv[3])
								+ (  1.3237268119 * yv[4])  + (  2.5142646221 * yv[5])
								+ ( -3.8888805689 * yv[6])  + ( -7.3557006394 * yv[7])
								+ (  7.6304375952 * yv[8])  + ( 13.8987072240 * yv[9])
								+ (-10.5788125720 * yv[10]) + (-17.5060272390 * yv[11])
								+ ( 10.6404877430 * yv[12]) + ( 14.5198001250 * yv[13])
								+ ( -7.7217913562 * yv[14]) + ( -7.3763465895 * yv[15])
								+ (  3.7745859324 * yv[16]) + (  1.8294126096 * yv[17]);
		
		signal[t] = (float)(yv[18]);
	}
  
	return;
}
