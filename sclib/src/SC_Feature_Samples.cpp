/**************************************************************************/
/*    Responsibility:																											*/
/*      - Provides the signal samples in a SV_Data container and chopped  */
/*        into frames                                                     */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 29.01.2009																								*/
/**************************************************************************/

#include "SC_Aux.h"
#include "SC_Feature_Samples.h"
#include <SV_Error.h>

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_Samples::SC_Feature_Samples(int sampleRate, int frameLength, int frameStep) : SV_Feature(){
  this->Para.WinSz = frameLength;
	this->Para.StpSz = frameStep;
	this->Para.SRate = sampleRate;
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_Samples::~SC_Feature_Samples() {
	
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data* SC_Feature_Samples::ExtractFeature(void) {
	SV_Data *pData = NULL;
 	float *sigBuf = NULL;
  long int sigLen, frameCnt, sampleCnt;
  int frameCount;

	if (IsSigLoaded()) { 
		sigBuf = GetSig();
		sigLen = GetLen();
	}	else {
		return NULL;
	}

	frameCount = (int)(sclib::getRowCount(sigLen, Para.WinSz, Para.StpSz));
	if (frameCount > 0) {
		pData = new SV_Data;
		if (pData == NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}

		pData->Row = frameCount;
		pData->Col = this->Para.WinSz;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureSamples;

		for (frameCnt = 0; frameCnt < pData->Row; frameCnt++) {
			for (sampleCnt = 0; sampleCnt < this->Para.WinSz; sampleCnt++) {
				pData->Mat[frameCnt][sampleCnt] = sigBuf[frameCnt * this->Para.StpSz + sampleCnt];
			}
		}
	}

  return pData;
}
