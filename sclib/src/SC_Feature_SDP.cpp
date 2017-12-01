/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for extracting Symmetrized Dot Patterns per frame         */
/*																																				*/
/*    Author  : Bing Shi                           												*/
/*    Date    : 11.04.2006																								*/
/**************************************************************************/

#include "SC_Feature_SDP.h"
#include "SC_SDP.h"
#include <SV_Error.h>
//#include "SC_Signature.h"
//#include "SC_Model_PictureCentroids.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_SDP::SC_Feature_SDP(int sampleRate, int frameLength, int frameStep, double preemphasize, int m, int lag, int color, int n, int pictureSize, int tau) : SV_Feature() {
	this->Para.WinSz = frameLength;
	this->Para.StpSz = frameStep;
	this->Para.SRate = sampleRate;
	this->Para.Alpha = preemphasize;

	this->m = m;
	this->lag = lag;
	this->color = color;
	this->n = n;
	this->pictureSize = pictureSize;
	this->tau = tau;
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_SDP::~SC_Feature_SDP() {
	
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data* SC_Feature_SDP::ExtractFeature(void) {
	SV_Data *pData = NULL;
 	float *sigBuf = NULL;
  double *frame;
  long int sigLen, frameCnt, sampleCnt;
	SC_SDP *pSDP = new SC_SDP(this->m, this->lag, this->color, this->n, this->pictureSize, this->tau);
	unsigned long int **bmp;
	int frameCount;
  
	if (IsSigLoaded()) { 
		sigBuf = GetSig();
		sigLen = GetLen();
	}	else {
		return (NULL);
	}

	frameCount = (int)(sclib::getRowCount(sigLen, this->Para.WinSz, this->Para.StpSz));
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
		pData->Col = this->pictureSize * this->pictureSize;
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureSDP;
		pData->Alloc();
		
		MArray_1D(frame, this->Para.WinSz, double, "SC_Feature_SDP.ExtractFeature: frame");

		char* fileName = new char[sclib::bufferSize];
		
		for (frameCnt = 0; frameCnt < pData->Row; frameCnt++) {
			for (sampleCnt = 0; sampleCnt < this->Para.WinSz; sampleCnt++) {
				frame[sampleCnt] = sigBuf[frameCnt * this->Para.StpSz + sampleCnt];
			}
			bmp = pSDP->createSDP(frame, this->Para.WinSz, false, true);
			pSDP->sdp2vector(bmp, pData->Mat[frameCnt]);
			
			/*
			char fileName[sclib::bufferSize];
			sprintf(fileName,"%s%5d%s\0", "SDP", frameCnt, ".bmp");	
			pSDP->saveBitMap(fileName, bmp, false);
			*/
			
			MFree_2D(bmp);
		}

		/*for(int i = 0; i < 2; i++){
			sprintf(fileName, "%s%d%s\0", "../TESTSDP/2/Small_Pictures/k", i+1, ".bmp");
			bmp = pSDP->loadBitMap(fileName);
			//sprintf(fileName, "%s%d%s\0", "D:/testset_shi/o", i, ".bmp");
			//pSDP->saveBitMap(fileName, bmp,false);
			pSDP->sdp2vector(bmp, pData->Mat[i]);
			MFree_2D(bmp);
		}*/

		MFree_1D(fileName);
		MFree_1D(frame);
	}
	MFree_0D(pSDP);

  return pData;
}
