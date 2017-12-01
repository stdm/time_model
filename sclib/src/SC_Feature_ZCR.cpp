/**************************************************************************/
/*    Class to extract Zero Crossing Rate Feature                         */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 12.02.2006																								*/
/**************************************************************************/

#include "SC_Feature_ZCR.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_ZCR::SC_Feature_ZCR(int sampleRate, int frameSize, int frameStep, double preemphasis, bool useChebyshev, bool scaleResults) : SV_Feature() {
	this->Para.StpSz = frameStep;
	this->Para.WinSz = frameSize;
	this->Para.Alpha = preemphasis;
	this->Para.SRate = sampleRate;
	this->useChebyshev = useChebyshev;
	this->scaleResults = scaleResults;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Feature_ZCR::~SC_Feature_ZCR() {

}

//====================================================================================================================
//	extract Zero Crossing Rate per frame, normalized by 2*(frameSize-1)
//====================================================================================================================
SV_Data *SC_Feature_ZCR::ExtractFeature(void) {
	float	*SigBuf;
	long SigLen;
	short	ThisSample, LastSample;
	SV_Data *pData = NULL;	
	int frameCount, tempZCR;
	double scaleFactor = (this->Para.WinSz * 1000.0) / (double)(this->Para.SRate); //needed to convert frame-zcr to zcr per ms

	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}	else {
		return (NULL);
	}

	// number of frames in the signal
	frameCount = (int)(sclib::getRowCount(SigLen, Para.WinSz, Para.StpSz)); //(SigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1;
	if (frameCount > 0) {
		//apply special hardcoded chebyshev filter, if wished
		if (this->useChebyshev == true) {
			if (this->Para.SRate = 8000) {
				chebyshevHighPass(SigBuf, SigLen);
			} else { //use standrad FIR filter with the chebysev's passband instead
				REPORT_ERROR(SVLIB_BadArg, "The hardcoded Chebyshev highpass can only be used with 8kHz signals");
				FIR_filtering((100.0/(double)(this->Para.SRate))*2.0, (3999.0/(double)(this->Para.SRate))*2.0);
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
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureZCR;

		for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex) {
			tempZCR = 0;
			
			for (int i = 1; i < this->Para.WinSz; ++i) {
				ThisSample = (short)SigBuf[frameIndex * this->Para.StpSz + i];
				LastSample = (short)SigBuf[frameIndex * this->Para.StpSz + i - 1];
				
				if (((ThisSample>=0) && (LastSample<0)) || ((ThisSample<0) && (LastSample>=0))) {
					tempZCR++; //sign changed, a zero crossing found 
				}
			}
			
			pData->Mat[frameIndex][0] = ((float)(tempZCR)); // / (2*(Para.WinSz-1)); //by thilo: killed old scale factor, don't see any rationale behind the /2... (update: it was due to a missinterpretation of the formula in li/kuo/zhangs paper on audio type classification!)
			if (this->scaleResults == true) {
				pData->Mat[frameIndex][0] /= (float)(scaleFactor); //per-frame => per-ms
			}
		}  // for frameIndex
	}
	
	return pData;
}

//====================================================================================================================
//	A hardcoded implementation of the Chebyshev highpass filter as specified in "Cepstrum-Based Pitch Detection 
//  Using a New Statistical V/UV Classification Algorithm" (Ahmadi, Spanias 1999): 9th order, 8kHz sample-rate, pass-
//  band 100-4000Hz; the filter was generated using http://www-users.cs.york.ac.uk/~fisher/mkfilter 
//====================================================================================================================
void SC_Feature_ZCR::chebyshevHighPass(float* signal, long int length) {
  int i; 
	long int t;
	const double GAIN = 1.355851745e+00;
	double xv[10], yv[10]; //in the originally created code, xv has size NZEROS+1 and yv NPOLES+1 with NZEROS=NPOLES=9

	//Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
	//Command line: /www/usr/fisher/helpers/mkfilter -Bu -Bp -o 9 -a 2.5000000000e-02 4.2500000000e-01 -l
	//initializations are done by the compiler via the static keyword in the originally generated code...
	for (i = 0; i < 9+1; i++) {
		xv[i] = 0.0;
		yv[i] = 0.0;
	}

	for (t = 0; t < length; t++) {
		for (i = 0; i < 9; i++) {
			xv[i] = xv[i+1];
			yv[i] = yv[i+1];
		}

		xv[9] = signal[t] / GAIN; //originally, this line came beofer the initialisation of yv[0]-yv[8], but this doesn't matter 'case they are both independent of xv[9]
		
    yv[9] =   (xv[9] - xv[0]) + 9 * (xv[1] - xv[8]) + 36 * (xv[7] - xv[2])
              + 84 * (xv[3] - xv[6]) + 126 * (xv[5] - xv[4])
              + (  0.5406311218 * yv[0]) + ( -5.2066128453 * yv[1])
              + ( 22.2849043680 * yv[2]) + (-55.6429164670 * yv[3])
              + ( 89.3279148920 * yv[4]) + (-95.6244500340 * yv[5])
              + ( 68.2626559800 * yv[6]) + (-31.3372242870 * yv[7])
              + (  8.3950972663 * yv[8]);
		
		signal[t] = (float)(yv[9]);
	}

	return;
}
