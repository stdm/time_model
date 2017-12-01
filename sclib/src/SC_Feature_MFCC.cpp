/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for extracting MFCC's via different algorithms: svlib's		*/
/*        original implementation or sclib's new one (via SC_Feature_FbE) */
/*																																				*/
/*    Author  : Thilo Stadelmann															            */
/*    Date    : 27.02.2004																								*/
/**************************************************************************/

#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_Feature_MFCC.h"
#include "SC_Feature_FbE.h"
#include <SV_Feature_MFCC.h>
#ifdef SC_USE_MATLAB
	#include "SC_Matlab.h"
#endif

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_MFCC::SC_Feature_MFCC(int sampleRate, int frameSize, int frameStep, int MFCCorder, int filterBankSize, int fftSize, unsigned short int window, double preemphasize, unsigned char coeffSelection[8], int dEnergy, bool CMN, bool addDeltas, bool addDeltaDeltas, int method, unsigned char sclib_frequencyScale, unsigned char sclib_smoothing, double sclib_minFilterBankFrequency, double sclib_maxFilterBankFrequency) : SV_Feature() {
	this->Para.WinSz = frameSize;
	this->Para.StpSz = frameStep;
	this->Para.MFCC_Order = MFCCorder;
	this->Para.NFilter = filterBankSize;
	this->Para.FFTSz = fftSize;
	this->Para.HammingWin = window;
	this->Para.RmvSilence = 0; //we don't use this one, it will break up the coherencey of frames with ground truth data because frames will be omitted, not just marked
	this->Para.DEnergy = dEnergy;
	this->Para.SRate = sampleRate;
	this->Para.Alpha = preemphasize;

	this->method = method;
	for (char i = 0; i < 8; i++) {
		this->coeffSelection[i] = (coeffSelection!=NULL)?coeffSelection[i]: 0xFF;
	}
	this->CMN = CMN;
	this->addDeltas = addDeltas;
	this->addDeltaDeltas = addDeltaDeltas;
	this->svlib_removeSilence = svlib_removeSilence;
	this->sclib_frequencyScale = sclib_frequencyScale;
	this->sclib_smoothing = sclib_smoothing;
	this->sclib_minFilterBankFrequency = sclib_minFilterBankFrequency;
	this->sclib_maxFilterBankFrequency = sclib_maxFilterBankFrequency;
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_MFCC::~SC_Feature_MFCC() {
	
}

//====================================================================================================================
// method to select specific coeeficents out of the feature-vector
// selectionString is a bit-flag-map, where each setted bit corresponds to a col of pData which shouldn't be discarded
// the first bit in selectionString[0] corresponds with the first column, the last bit in selectionstring[7] with the
// 64th column.
//====================================================================================================================
void SC_Feature_MFCC::selectCoefficients(SV_Data* &pData, unsigned char selectionString[8]) {
	int newCoeff = 0, oldCoeff = pData->Col;
	float **oldMat = NULL; 

	//get new col-count and init matrix accordingly
	for (int x = 0; x < pData->Col; x++) {
		//if (sclib::bitTest(selectionMap, sclib::bit(x+1)) == true) {
		if (sclib::bitTest(selectionString[x/8], sclib::bit(x%8)) == true) { //the x'th bit in selectionString is the x%8'th bit in selectionString[x/8]
			newCoeff++;
			//printf("%d, ", x+1);
		}
	}

	if (newCoeff == pData->Col) {
		return;
	}

	oldMat = pData->Mat; //save pointer to old matrix...
	pData->setJustLinked(true); //...ensure that it doesn't get killed by Alloc()...
	pData->Col = newCoeff;
	pData->Alloc(); //...and create new data matrix inside old SV_Data object

	//copy selected features
	for (int y = 0; y < pData->Row; y++) {
		newCoeff = 0;
		for (int x = 0; x < oldCoeff; x++) {
			//if (sclib::bitTest(selectionMap, sclib::bit(x+1)) == true) {
			if (sclib::bitTest(selectionString[x/8], sclib::bit(x%8)) == true) { //the x'th bit in selectionString is the x%8'th bit in selectionString[x/8]
				pData->Mat[y][newCoeff++] = oldMat[y][x];
			}
		}
	}

	//destroy old matrix, return new one
	MFree_2D(oldMat);

	return;
}

//==========================================
// This is the engine of deriving features
//==========================================
SV_Data *SC_Feature_MFCC::ExtractFeature(void) {
	SV_Data *pMFCCs = NULL, *pDeltas = NULL, *pDeltaDeltas = NULL, *pTemp = NULL;
	SV_Feature *pExtractor = NULL;

	switch (this->method) {
		case sclib::modeSClib: {
			pExtractor = new SC_Feature_FbE(this->Para.SRate, this->Para.WinSz, this->Para.StpSz, this->Para.NFilter, this->Para.FFTSz, this->Para.HammingWin, this->Para.Alpha, this->sclib_minFilterBankFrequency, this->sclib_maxFilterBankFrequency, this->Para.MFCC_Order, (this->Para.DEnergy!=0)?true:false, this->sclib_frequencyScale, sclib::resultCepstrum, this->sclib_smoothing, 1.0);

			//extract mfccs
			pExtractor->setSignal(this->Sig, this->Len, false);
			pMFCCs = pExtractor->ExtractFeature();
			pExtractor->setSignal(NULL, 0, false);

			//set header
			if (pMFCCs != NULL) {
				pMFCCs->Hdr.Signature[1] = sclib::modeSClib; //encode which extractor was used
			}
			break;
		}

		case sclib::modeSVlib: {
			pExtractor = new SV_Feature_MFCC();
			
			//set parameters
			pExtractor->Para.SRate = this->Para.SRate;
			pExtractor->Para.WinSz = this->Para.WinSz;
			pExtractor->Para.StpSz = this->Para.StpSz;
			pExtractor->Para.MFCC_Order = this->Para.MFCC_Order;
			pExtractor->Para.NFilter = this->Para.NFilter;
			pExtractor->Para.FFTSz = this->Para.FFTSz;
			pExtractor->Para.HammingWin = this->Para.HammingWin;
			pExtractor->Para.RmvSilence = this->Para.RmvSilence;
			pExtractor->Para.DEnergy = this->Para.DEnergy;
			pExtractor->Para.Alpha = this->Para.Alpha; //this is not used inside SV_Feature_MFCC, so we need to do preemphasis here
			
			//extract mfccs
			pExtractor->setSignal(this->Sig, this->Len, false);
			pExtractor->PreEmphasize(pExtractor->Para.Alpha); //as said before, this is not done in ExtractFeature() in this class
			pMFCCs = pExtractor->ExtractFeature();
			pExtractor->setSignal(NULL, 0, false);

			//set header (svlib version doesn't do this)
			if (pMFCCs != NULL) {
				pMFCCs->Hdr.ID = sclib::featureMFCC;
				pMFCCs->Hdr.frameSize = this->Para.WinSz;
				pMFCCs->Hdr.frameStep = this->Para.StpSz;
				pMFCCs->Hdr.sampleRate = this->Para.SRate;
				pMFCCs->Hdr.Signature[1] = sclib::modeSVlib; //encode which extractor was used
			}
			break;
		}

		case sclib::modeSlaney: {
#ifdef SC_USE_MATLAB
			//initialize Matlab
			SC_Matlab mHandler;
			if (!mHandler.initialize()) {
				return NULL;
			}

			//prepare input arguments in Matlab format
	    double *dSig = NULL;
			MArray_1D(dSig, this->Len, double, "SC_Feature_MFCC.ExtractFeature: dSig");
			for (int x = 0; x < this->Len; x++) {
				dSig[x] = this->Sig[x] / 32768.0; //normalize to [-1,1) as done in Matlab's wavread(); TODO: this assumes 16bit/sample original audio data
			}
			MatlabRealMatrix input = mHandler.vector2matlab(dSig, this->Len, false);
			MFree_1D(dSig);

			MatlabRealMatrix samplingRate = mHandler.scalar2matlab(this->Para.SRate);
			MatlabRealMatrix frameRate = mHandler.scalar2matlab(this->Para.SRate/this->Para.StpSz); //in Matlab, the windowStep (in samples) is computed as windowsStep=samplingRate/frameRate

			//call Matlab function
			//uses (hard coded): preemphasis with alpha=0.97, a hamming window, fftSize of 512, 40 filters, lowest frequency of 133.333Hz, results in 13 MFCCs
			MatlabRealMatrix ceps = mHandler.sc_mlfMfcc(input, samplingRate, frameRate);
			pMFCCs = mHandler.matlab2svdata(ceps, true);
			pMFCCs->Hdr.ID = sclib::featureMFCC;
			pMFCCs->Hdr.frameSize = 256; //hard coded in the auditory toolbox
			pMFCCs->Hdr.frameStep = this->Para.StpSz;
			pMFCCs->Hdr.sampleRate = this->Para.SRate;
			pMFCCs->Hdr.Signature[1] = sclib::modeSlaney; //encode which extractor was used

			//collect output and clean up memory
			mHandler.freeMatrix(input);
			mHandler.freeMatrix(samplingRate);
			mHandler.freeMatrix(frameRate);
			mHandler.freeMatrix(ceps);

			//terminate Matlab: done in lib-constructor (SC_Lib)
			//mHandler.terminate();
#else
			REPORT_ERROR(SVLIB_BadArg, "Choosen MFCC extractor by Slaney not available!\n");
			return NULL;
#endif
			break;
		}

		default:
			REPORT_ERROR(SVLIB_BadArg, "Choosen MFCC extractor unknown!\n");
			return NULL;
	}

	if (pMFCCs != NULL) {
		//CMN
		if (this->CMN == true) {
			pExtractor->CMN(pMFCCs);
		}

		//select specific coefficients, if wished
		selectCoefficients(pMFCCs, this->coeffSelection);

		//add deltas, delta-deltas
		if (this->addDeltas == true) {
			pDeltas = pExtractor->Dynamic(*pMFCCs);
			if (this->addDeltaDeltas == true) {
				pDeltaDeltas = pExtractor->Dynamic(*pDeltas);
			}
			pTemp = pExtractor->Concat(*pMFCCs, *pDeltas);
			MFree_0D(pMFCCs);
			MFree_0D(pDeltas);
			pMFCCs = pTemp;
			if (pDeltaDeltas != NULL) {
				pTemp = pExtractor->Concat(*pMFCCs, *pDeltaDeltas);
				MFree_0D(pMFCCs);
				MFree_0D(pDeltaDeltas);
				pMFCCs = pTemp;
			}
		}
	}

	MFree_0D(pExtractor);

	return pMFCCs;
}
