#include <math.h>
#include <float.h>
#include <limits>
#include "SC_Feature_SpectrumFlux.h"
#include "SC_Transform.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_SpectrumFlux::SC_Feature_SpectrumFlux(int sampleRate, int frameSize, int frameStep, int fftLength, unsigned short window, double preemphasis) : SV_Feature() {
	this->Para.StpSz = frameStep;
  this->Para.WinSz = frameSize;
	this->Para.FFTSz = fftLength;
	this->Para.Alpha = preemphasis;
  this->Para.HammingWin = window; 
	this->Para.SRate = sampleRate;
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Feature_SpectrumFlux::~SC_Feature_SpectrumFlux() {

}

//====================================================================================================================
//	exptracts Spectrum Flux per frame: the avergage difference per frequency-bin in magnitude spectra of two succesive 
//  frames
//====================================================================================================================
SV_Data *SC_Feature_SpectrumFlux::ExtractFeature(void) {
	SV_Data	*pData = NULL;
	double spectrumFlux, temp;
	double smallDelta = numeric_limits<double>::epsilon(); //TODO: was DBL_EPSILON before...
  double *frame = NULL, *magnitudeSpectrum = NULL, *lastSpectrum = NULL;
	float *sig;
	int frameCount, spectrumCount = Para.FFTSz / 2 + 1, len;
	SC_Transform *pTrans = NULL;

	if (IsSigLoaded()) { 
		sig = GetSig();
		len = GetLen();
	}	else {
		return (NULL);
	}

  // number of frames in the signal
	frameCount = (int)(sclib::getRowCount(len, Para.WinSz, Para.StpSz)); //(len / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1;

	if (frameCount > 1) {
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		MArray_1D(frame, Para.StpSz, double, "SC_Feature_SpectrumFlux.initSpectrumValues: frame");
		pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin);

		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}
		pData->Row=frameCount;
		pData->Col=1;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureSpectrumFlux;

		// for all frames in the signal
		for (int frameIdx=0; frameIdx<frameCount; ++frameIdx) {
			// make a copy of the current frame
			for (int i=0; i<Para.StpSz; ++i) {
				frame[i] = sig[frameIdx*Para.StpSz+i];
			}

			// generate fourier-transform
			magnitudeSpectrum = pTrans->magnitudeSpectrum(frame, Para.WinSz);
		
			if (frameIdx > 0) {
				spectrumFlux = 0.0;
				for (int k = 0; k<spectrumCount; ++k) { //TODO: k=1..spectrumCount?!?
					temp = log(magnitudeSpectrum[k] + smallDelta) - log(lastSpectrum[k] + smallDelta);
					spectrumFlux += temp * temp;
				}
				pData->Mat[frameIdx][0] = (float)(spectrumFlux / (double)(spectrumCount-1));
				if (!sclib::isFinite(pData->Mat[frameIdx][0])) {
					frameIdx = frameIdx;
				}
			}

			MFree_1D(lastSpectrum);
			lastSpectrum = magnitudeSpectrum;
		}
		MFree_1D(lastSpectrum);
		pData->Mat[0][0] = pData->Mat[1][0]; //set the SF-value of the first frame to that of the second one, because the first one has no predecessor

		MFree_1D(frame);
		MFree_0D(pTrans);
	}

	return(pData);
}
