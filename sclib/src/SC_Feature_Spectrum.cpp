/**************************************************************************/
/*    Copied from:																												*/
/*      - This is a altered version of SC_Feature_FbE               			*/
/*																																				*/
/*    Responsibility:																											*/
/*      - Class for extracting Powerspectrum features as used in the  		*/
/*				SC_Enhancement class.                                           */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 30.03.2005																								*/
/**************************************************************************/

#include <math.h>
#include "SC_Aux.h"
#include "SC_Feature_Spectrum.h"
#include "SC_Transform.h"
#include <SV_Error.h>

//====================================================================================================================
// constructor
// for use with SC_Enhancement, it should hold:
//  - FFTlength = frameLength
//  - StpSz     = frameLength / 2
//  - therefore, frameLength must be a power of 2
// but that is checked for there...
//====================================================================================================================
SC_Feature_Spectrum::SC_Feature_Spectrum(int sampleRate, int frameLength, int frameStep, double preemphasize, int FFTsize, unsigned short window, bool logarithmize, bool createPhase) : SV_Feature(){
  this->Para.WinSz = frameLength;
	this->Para.StpSz = frameStep;
	this->Para.FFTSz = FFTsize;
	this->Para.SRate = sampleRate;
  this->Para.Alpha = preemphasize;
  this->Para.HammingWin = window;
	this->logarithmize = logarithmize;
	this->createPhase = createPhase;
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_Spectrum::~SC_Feature_Spectrum() {
	
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data* SC_Feature_Spectrum::ExtractFeature(void) {
	SV_Data *pData = NULL, *pPhase;
  SC_Transform *pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin);
 	float *sigBuf = NULL;
  double *frame = NULL;
  double *powerSpectrum, *phaseSpectrum = NULL;
  long int sigLen, frameCnt, sampleCnt;
  int frameCount;

	if (IsSigLoaded()) { 
		sigBuf = GetSig();
		sigLen = GetLen();
	}	else {
		return (NULL);
	}

	frameCount = (int)(sclib::getRowCount(sigLen, Para.WinSz, Para.StpSz)); //(sigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1;
	if (frameCount > 0) {
		// Preemphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}

		pData->Row = frameCount;
		pData->Col = this->Para.FFTSz/2 +1; //this is the length of the resulting power spectrum
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureSpectrum;
		pData->Hdr.Signature[1] = 0;

		if (this->createPhase == true) {
			pPhase = new SV_Data(pData->Row, this->Para.FFTSz); //also extract the pase spectrum and return it in a 2nd data container, identifiable by the 1 in Hdr.Signature[1]
			pPhase->Hdr = pData->Hdr;
			pPhase->Hdr.Signature[1] = 1;
			pData->Next = pPhase;
		}

		MArray_1D(frame, this->Para.WinSz, double, "SC_Feature_Spectrum.ExtractFeature: frame");

		for (frameCnt = 0; frameCnt < pData->Row; frameCnt++) {
			for (sampleCnt = 0; sampleCnt < this->Para.WinSz; sampleCnt++) {
				frame[sampleCnt] = sigBuf[frameCnt * this->Para.StpSz + sampleCnt];
			}
			//sclib::vectorOut("frame.txt", frame, this->Para.WinSz, true);

			if (this->createPhase == true) {
				powerSpectrum = pTrans->powerSpectrum(frame, this->Para.WinSz, phaseSpectrum);
			} else {
				powerSpectrum = pTrans->powerSpectrum(frame, this->Para.WinSz);
			}
			//sclib::vectorOut("power.txt", powerSpectrum, this->Para.FFTSz/2+1, true);
			//sclib::vectorOut("phase.txt", phaseSpectrum, this->Para.FFTSz, true);

			//MFree_0D(frame);
			//pTrans->powerPhase2ft(powerSpectrum, phaseSpectrum, false, this->Para.WinSz, false, false);
			//frame = pTrans->ifft(powerSpectrum, phaseSpectrum, this->Para.WinSz);
			//pTrans->tapering(frame, this->Para.WinSz, this->Para.HammingWin, 1.0, false, true);
			//sclib::vectorOut("frame2.txt", frame, this->Para.WinSz, true);

			for (sampleCnt = 0; sampleCnt < this->Para.FFTSz; sampleCnt++) {
				if (sampleCnt < this->Para.FFTSz/2 +1) {
					if (this->logarithmize == true) {
						pData->Mat[frameCnt][sampleCnt] = (float)((powerSpectrum[sampleCnt] <= 0.0) ? 0.0 : log(powerSpectrum[sampleCnt]));
					} else {
						pData->Mat[frameCnt][sampleCnt] = (float)(powerSpectrum[sampleCnt]);
					}
				}
				if (this->createPhase == true) {
					pPhase->Mat[frameCnt][sampleCnt] = (float)(phaseSpectrum[sampleCnt]);
				}
			}

			MFree_1D(phaseSpectrum);
			MFree_1D(powerSpectrum);
		}

		MFree_1D(frame);
	}
	MFree_0D(pTrans);

  return(pData);
}

//====================================================================================================================
// another engine...
// the signal is already framed (each row in pFrames = one frame)... compute a new dataset containing their spectra
//====================================================================================================================
SV_Data* SC_Feature_Spectrum::ExtractFeature(SV_Data *pFrames) {
	SV_Data *pData = NULL, *pPhase;
  SC_Transform *pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin);
  double *frame = NULL;
  double *powerSpectrum, *phaseSpectrum = NULL;
  long int frameCnt, sampleCnt;
  
	pData = new SV_Data;
	if (pData==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}

	pData->Row = pFrames->Row;
	pData->Col = this->Para.FFTSz/2 +1; //this is the length of the resulting power spectrum
	pData->Alloc();
  pData->Hdr.frameSize = this->Para.WinSz;
  pData->Hdr.frameStep = this->Para.StpSz;
  pData->Hdr.sampleRate = this->Para.SRate;
	pData->Hdr.ID = sclib::featureSpectrum;
	pData->Hdr.Signature[1] = 0;

	if (this->createPhase == true) {
		pPhase = new SV_Data(*pData, true); //also extract the pase spectrum and return it in a 2nd data container, identifiable by the 1 in Hdr.Signature[1]
		pPhase->Alloc();
		pData->Next = pPhase;
		pPhase->Hdr.Signature[1] = 1;
	}

  MArray_1D(frame, this->Para.WinSz, double, "SC_Feature_Spectrum.ExtractFeature: frame");

  for (frameCnt = 0; frameCnt < pData->Row; frameCnt++) {
    for (sampleCnt = 0; sampleCnt < this->Para.WinSz; sampleCnt++) {
      frame[sampleCnt] = pFrames->Mat[frameCnt][sampleCnt];
    }
  
		if (this->createPhase == true) {
			powerSpectrum = pTrans->powerSpectrum(frame, this->Para.WinSz, phaseSpectrum);
		} else {
			powerSpectrum = pTrans->powerSpectrum(frame, this->Para.WinSz);
		}

    for (sampleCnt = 0; sampleCnt < this->Para.FFTSz/2 +1; sampleCnt++) {
			if (this->logarithmize == true) {
				pData->Mat[frameCnt][sampleCnt] = (float)((powerSpectrum[sampleCnt] <= 0.0) ? 0.0 : log(powerSpectrum[sampleCnt]));
			} else {
				pData->Mat[frameCnt][sampleCnt] = (float)(powerSpectrum[sampleCnt]);
			}
			if (this->createPhase == true) {
				pPhase->Mat[frameCnt][sampleCnt] = (float)(phaseSpectrum[sampleCnt]);
			}
    }

		MFree_1D(phaseSpectrum);
    MFree_1D(powerSpectrum);
  }

  MFree_0D(pTrans);
  MFree_1D(frame);

  return(pData);
}
