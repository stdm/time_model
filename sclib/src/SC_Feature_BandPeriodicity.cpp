#include "SC_Feature_BandPeriodicity.h"
#include "SC_Transform.h"
#include "SC_Aux.h"
#include <SV_Error.h>
#include <GN_Filter.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_BandPeriodicity::SC_Feature_BandPeriodicity(int sampleRate, int frameSize, int frameStep, double preemphasis, unsigned short int FFTsize, unsigned short int window) : SV_Feature() {
	this->Para.StpSz = frameStep;
  this->Para.WinSz = frameSize;
	this->Para.Alpha = preemphasis;
  this->Para.FFTSz = FFTsize;
  this->Para.HammingWin = window; 
	this->Para.SRate = sampleRate; 
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Feature_BandPeriodicity::~SC_Feature_BandPeriodicity() {

}

//====================================================================================================================
//	returns a new array (must be destructed on caller side!) containing a bandpass-filtered version of the loaded 
//  signal; low and high are normalized to (0, 1)
//====================================================================================================================
float* SC_Feature_BandPeriodicity::filterSignal(float* signal, unsigned long int len, double lowCut, double highCut) {
	float *filteredSignal;
	
  MArray_1D(filteredSignal, len, float, "SC_Feature_BandPeriodicity.filterSignal: filteredSignal");
	
  for (unsigned long int i = 0; i < len; i++) {
    filteredSignal[i] = signal[i];
  }
	
  GN_Filter filter = GN_Filter();
	filter.BandPass(lowCut, highCut);
	filter.FIR_Filter(filteredSignal, len);
	
  return filteredSignal;
}

//====================================================================================================================
//	extract Band Periodicity feature per frame: the maximum peak in the normalized correlation fuction between the 
//  current and the previous frame per frequency sub-band (0-1/8, 1/8-1/4, 1/4-1/2, 1/2-1 normalized frequency)
//====================================================================================================================
SV_Data *SC_Feature_BandPeriodicity::ExtractFeature(void) {
  double **subBands = NULL, *twoFrames = NULL, *normalizedCorrelation = NULL, maximum;
  float *SigBuf, *filteredSignal = NULL;
	int SigLen, frameCount;
  int i; 
	SV_Data	*pData = NULL;
  SC_Transform *pTrans = NULL;
  
	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}	else {
		return (NULL);
	}

	// number of frames in the signal
	frameCount = (int)(sclib::getRowCount(SigLen, Para.WinSz, Para.StpSz)); //(SigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1;

	if (frameCount > 1) {
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}
		pData->Row = frameCount;
		pData->Col = 4;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureBandPeriodicity;

		pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin);

		MArray_1D(twoFrames, 2*this->Para.WinSz, double, "SC_Feature_BandPeriodicity.ExtractFeature:twoFrames");
		MArray_2D(subBands, 4, 2, double, "SC_Feature_BandPeriodicity.ExtractFeature: subBands");
		subBands[0][0] = 0.0;   subBands[0][1] = 0.125;
		subBands[1][0] = 0.125; subBands[1][1] = 0.25;
		subBands[2][0] = 0.25;  subBands[2][1] = 0.5;
		subBands[3][0] = 0.5;   subBands[3][1] = 1.0;

		//for all subbands
		for (unsigned long int subBand = 0; subBand < 4; subBand++) {
			filteredSignal = filterSignal(SigBuf, SigLen, subBands[subBand][0], subBands[subBand][1]);

			//for all frames
			for (int frameNr = 1; frameNr < frameCount; frameNr++) {
				// make a copy of a frame into frameArray
				for (i = 0; i < this->Para.WinSz; i++) {
					twoFrames[i] = filteredSignal[(frameNr-1)*this->Para.StpSz+i]; //previous frame
					twoFrames[i + this->Para.WinSz] = filteredSignal[frameNr*this->Para.StpSz+i]; //current frame
				}

				//extract maximum peak of normalized correlation between the two frames
 				normalizedCorrelation = pTrans->normalizedCorrelation(twoFrames, this->Para.WinSz);
				maximum = 0.0;
				for (i = 0; i < this->Para.WinSz-1; i++) { //the correlation-result has one element less because of no correlation-computation for 0-shift (would be always =1)
					if (maximum < normalizedCorrelation[i]) {
						maximum = normalizedCorrelation[i];
					}
				}
				pData->Mat[frameNr][subBand] = (float)(maximum);
				if (!sclib::isFinite(pData->Mat[frameNr][subBand])) {
					frameNr = frameNr;
				}

				MFree_1D(normalizedCorrelation);
			}
	  
			pData->Mat[0][subBand] = pData->Mat[1][subBand]; //there is no (other) meaningful value for the first frame because it has no predecessor

			MFree_1D(filteredSignal);
		}
		MFree_2D(subBands);
		MFree_1D(twoFrames);
		MFree_0D(pTrans);
	}

	return(pData);
}
