#include <math.h>
#include "SC_Feature_BrightnessBandwidth.h"
#include "SC_Transform.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_BrightnessBandwidth::SC_Feature_BrightnessBandwidth(int sampleRate, int frameSize, int frameStep, int fftLength, unsigned short window, double preemphasis) : SV_Feature() {
	this->Para.SRate = sampleRate;
	this->Para.StpSz = frameStep;
	this->Para.WinSz = frameSize;
	this->Para.Alpha = preemphasis;
	this->Para.FFTSz = fftLength;
  this->Para.HammingWin = window;
}

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_BrightnessBandwidth::~SC_Feature_BrightnessBandwidth() {

}

//====================================================================================================================
//	extract brightness & bandwidth per frame
//====================================================================================================================
SV_Data *SC_Feature_BrightnessBandwidth::ExtractFeature(void) {
	SV_Data	*pData = NULL;
	double *frame, *powerSpectrum;
  float	*sigBuf;
	double brightness, bandwidth, frequency, energy, frequencyWeightedEnergy;
  long sigLen;
	int frameCount, fftCoeffCount;
	SC_Transform *pTrans;
	
	if (IsSigLoaded()) { 
		sigBuf = GetSig();
		sigLen = GetLen();
	}	else {
		return (NULL);
	}

	// number of frames in the signal
	frameCount = (int)(sclib::getRowCount(sigLen, Para.WinSz, Para.StpSz)); //(sigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1;
	if (frameCount > 0) {
		// number of the fft-coeffezients 
		fftCoeffCount = Para.FFTSz / 2 + 1;

		// creates the fft
		pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin);
		
			// Preeamphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}
		pData->Row = frameCount;
		pData->Col = 2;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureBrightnessBandwidth;

		MArray_1D(frame, this->Para.WinSz, double, "SC_Feature_Brigthness_Bandwidth.ExtractFeature: frameArray");
		
		// for all frames
		for (int frameNr = 0; frameNr < frameCount; ++frameNr) {
			// make a copy of a frame into frameArray
			for (int i=0; i<Para.WinSz; ++i) {
				frame[i] = sigBuf[frameNr*Para.StpSz+i];
			}

			// execute a fft on the frame
			powerSpectrum = pTrans->powerSpectrum(frame, this->Para.WinSz);

			//*****************************************************************
			// Calculation of the brightness in a frame
			//*****************************************************************
			energy = 0;
			frequencyWeightedEnergy = 0;
		
			// iterate over the fft-coefficients
			for (int i=0; i<fftCoeffCount; ++i) {
				energy += powerSpectrum[i];
				frequency = (((double)Para.SRate) / (2*fftCoeffCount)) * (i+1);
				frequencyWeightedEnergy += frequency*powerSpectrum[i];
			}

			// brigthness
			brightness = (energy != 0.0) ? frequencyWeightedEnergy / energy : 0.0;

			//*****************************************************************
			// Calculation of the bandwidth in a frame
			//*****************************************************************
			frequencyWeightedEnergy = 0;
		
			// iterate over the fft-coefficients
			for (int i=0; i<fftCoeffCount; ++i) {
				frequency = (((double)Para.SRate) / (2*fftCoeffCount)) * (i+1);
				frequencyWeightedEnergy += (frequency-brightness)*(frequency-brightness) * powerSpectrum[i];
			}
			
			// bandwidth
			bandwidth = (energy != 0.0) ? sqrt(frequencyWeightedEnergy / energy) : 0.0;

			MFree_1D(powerSpectrum);

			//*****************************************************************
			// Store the result in the matrix
			//*****************************************************************
			pData->Mat[frameNr][0] = (float) brightness;
			pData->Mat[frameNr][1] = (float) bandwidth;
		} // endfor: frameNr

		MFree_1D(frame);
		delete pTrans;
	}

	return(pData);
}

