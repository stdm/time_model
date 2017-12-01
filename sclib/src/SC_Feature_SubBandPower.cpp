#include <math.h>
#include "SC_Aux.h"
#include "SC_Feature_SubBandPower.h"
#include "SC_Transform.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Feature_SubBandPower::SC_Feature_SubBandPower(int sampleRate, int frameSize, int frameStep, int fftLength, unsigned short window, double preemphasis) : SV_Feature() {
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
SC_Feature_SubBandPower::~SC_Feature_SubBandPower() {

}

//====================================================================================================================
//  sum up the entries in the array from "from" to "to"
//====================================================================================================================
double SC_Feature_SubBandPower::sumArray(double *values, int from, int to) {
	double temp = 0;
	for (int i=from; i<to; ++i) {
		temp += values[i];
	}
	return temp;
}

//====================================================================================================================
//	extracts per frame: 1. Short Time Energy (derived from log of integral of powerspectrum)
//                      2. Ratio of sub-band 1 (integral of 0   - 1/8 powerspectrum) to exp(STE)
//                      3. Ratio of sub-band 2 (integral of 1/8 - 1/4 powerspectrum) to exp(STE)
//                      4. Ratio of sub-band 3 (integral of 1/4 - 1/2 powerspectrum) to exp(STE)
//                      5. Ratio of sub-band 4 (integral of 1/2 - 1   powerspectrum) to exp(STE)
//====================================================================================================================
SV_Data *SC_Feature_SubBandPower::ExtractFeature(void) {
	SV_Data *pData = NULL;
	double shortTimeEnergy, subBandEnergyOne, subBandEnergyTwo, subBandEnergyThree, subBandEnergyFour;
	double *frame, *powerSpectrum;
  float *sigBuf;
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
		//fftCoeffCount = Para.StpSz / 2 +1;
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
		pData->Col = 5;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureSubbandPower;

		MArray_1D(frame, this->Para.WinSz, double, "SC_Feature_SubBandPower.ExtractFeature: frame");
		
		// for all frames
		for (int frameNr = 0; frameNr < frameCount; ++frameNr) {
			// make a copy of a frame into frameArray
			for (int i=0; i<Para.WinSz; ++i) {
				frame[i] = sigBuf[frameNr*Para.StpSz+i];
			}

			// execute a fft on the frame
			powerSpectrum = pTrans->powerSpectrum(frame, this->Para.WinSz);		

			// calculate the whole energy of the frame
			shortTimeEnergy = (sumArray(powerSpectrum, 0, fftCoeffCount));  // this is the exp(Short-Time-Energy)

			// calculate the energy of the first subband [0;fftCoeffCount/8]
			subBandEnergyOne = (shortTimeEnergy != 0.0) ? (sumArray(powerSpectrum, 0, fftCoeffCount/8)) / shortTimeEnergy : 0.0;

			// calculate the energy of the second subband [fftCoeffCount/8;fftCoeffCount/4]
			subBandEnergyTwo = (shortTimeEnergy != 0.0) ? (sumArray(powerSpectrum, fftCoeffCount/8+1, fftCoeffCount/4)) / shortTimeEnergy : 0.0;

			// calculate the energy of the third subband [fftCoeffCount/4;fftCoeffCount/2]
			subBandEnergyThree = (shortTimeEnergy != 0.0) ? (sumArray(powerSpectrum, fftCoeffCount/4+1, fftCoeffCount/2)) / shortTimeEnergy : 0.0;

			// calculate the energy of the forth subband [fftCoeffCount/2;fftCoeffCount]
			subBandEnergyFour = (shortTimeEnergy != 0.0) ? (sumArray(powerSpectrum, fftCoeffCount/2+1, fftCoeffCount)) / shortTimeEnergy : 0.0;

			MFree_1D(powerSpectrum);

			pData->Mat[frameNr][0] = (shortTimeEnergy > 0.0) ? (float)log(shortTimeEnergy) : (float)(-87.3); //-87.336544750553102 is the log of the smallest non-neg. float number
			pData->Mat[frameNr][1] = (float) subBandEnergyOne;
			pData->Mat[frameNr][2] = (float) subBandEnergyTwo;
			pData->Mat[frameNr][3] = (float) subBandEnergyThree;
			pData->Mat[frameNr][4] = (float) subBandEnergyFour;

			//if (!sclib::isFinite(pData->Mat[frameNr][0]) || !sclib::isFinite(pData->Mat[frameNr][1]) || !sclib::isFinite(pData->Mat[frameNr][2]) || !sclib::isFinite(pData->Mat[frameNr][3]) || !sclib::isFinite(pData->Mat[frameNr][4])) {
			//  frameNr = frameNr;
			//}
		}

		MFree_1D(frame);
		MFree_0D(pTrans);
	}

	return pData;
}

