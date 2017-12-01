#include "SC_Feature_NFR.h"
#include "SC_Transform.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_NFR::SC_Feature_NFR(int sampleRate, int frameSize, double preemphasis, double nfrThreshold, int FFTsize, unsigned short window) : SV_Feature() {
	this->Para.StpSz = frameSize;
  this->Para.WinSz = frameSize;
	this->Para.Alpha = preemphasis;
	this->nfrThreshold = nfrThreshold; //depends on the SNR; 0.3 is a good choice to detect music
  this->Para.FFTSz = FFTsize;
  this->Para.HammingWin = window;
	this->Para.SRate = sampleRate;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Feature_NFR::~SC_Feature_NFR() {

}

//====================================================================================================================
//	extract the Noise Frame Ratio feature per frame; this means: it is not really a ratio, just a measure if a frame 
//  is considered as noise (1) or not (0); the ratio has to computed later for a given sub-clip
//====================================================================================================================
SV_Data *SC_Feature_NFR::ExtractFeature(void) {
  double *twoFrames = NULL, *normalizedCorrelation = NULL, maximum;
	float *SigBuf;
	long SigLen;
	int frameCount, noiseFrames = 0;
  unsigned long int i;
	SV_Data *pData = NULL;	
  SC_Transform *pTrans = NULL;

	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}	else {
		return (NULL);
	}

	// number of frames in the signal
	frameCount = (int)(sclib::getRowCount(SigLen, Para.WinSz, Para.StpSz)); //(SigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1; //(SigLen - Para.WinSz + Para.StpSz) / Para.StpSz;
	
	if (frameCount > 1) {
		// Preeamphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin);

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
		pData->Hdr.ID = sclib::featureNFR;

		MArray_1D(twoFrames, 2*this->Para.WinSz, double, "SC_Feature_BandPeriodicity.ExtractFeature:twoFrames");

		for (unsigned long int frameNr = 1; frameNr < (unsigned long int)(frameCount); frameNr++) {
			// make a copy of a frame into frameArray
			for (i = 0; i < (unsigned long int)(this->Para.WinSz); i++) {
				twoFrames[i] = SigBuf[(frameNr-1)*this->Para.StpSz+i]; //previous frame
				twoFrames[i + this->Para.WinSz] = SigBuf[frameNr*this->Para.StpSz+i]; //current frame
			}

			//extract maximum peak of normalized correlation between the two frames
 			normalizedCorrelation = pTrans->normalizedCorrelation(twoFrames, this->Para.WinSz);
			maximum = 0.0;
			for (i = 0; i < (unsigned long int)(this->Para.WinSz-1); i++) {  //the correlation-result has one element less because of no correlation-computation for 0-shift (would be always =1)
				if (maximum < normalizedCorrelation[i]) {
					maximum = normalizedCorrelation[i];
				}
			}
			pData->Mat[frameNr][0] = (maximum < this->nfrThreshold) ? (float)(1.0) : (float)(0.0);

			MFree_1D(normalizedCorrelation);
		}
		pData->Mat[0][0] = pData->Mat[1][0]; //there is no meaningful value for the first frame because it has no predecessor

		MFree_1D(twoFrames);
		MFree_0D(pTrans);
	}
	
	return pData;
}
