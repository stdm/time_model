/**************************************************************************/
/*    Responsibility:																											*/
/*		  - provides the possbility to extract all implemented features     */
/*        in one framework                                                */
/*      - helps to save/load features to a file                           */
/*      - helps aggregating features per sub-clip or to a combined vector */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.02.2006																								*/
/**************************************************************************/

#include <vector>
#include <map>
#include "SC_FeatureHandler.h"
#include "SC_Aux.h"
#include "SC_Feature_MFCC.h"
#include "SC_Feature_FbE.h"
#include "SC_Feature_Spectrum.h"
#include "SC_Feature_BandPeriodicity.h"
#include "SC_Feature_BrightnessBandwidth.h"
#include "SC_Feature_NFR.h"
#include "SC_Feature_SpectrumFlux.h"
#include "SC_Feature_SubBandPower.h"
#include "SC_Feature_ZCR.h"
#include "SC_Feature_STE.h"//Bing
#include "SC_Feature_LPC.h" //jun
#include "SC_Feature_LPCresidual.h" //jun
#include "SC_Feature_SDP.h" // Bing
#include "SC_Feature_LSP.h"
#include "SC_Feature_Pitch.h"
#include "SC_Feature_Formant.h"
#include "SC_Feature_Samples.h"
#include "SC_Conversion.h"
#include "SC_SignalHandler.h"
#include "SC_Signal_WAVE.h"
#include "SC_MatrixFunctions.h"
#include "SC_Clusterer.h"
#include "SC_Synthesis.h"
#include "SC_MD5.h"
#include <SV_DataIO.h>

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_FeatureHandler::SC_FeatureHandler(SC_TweakableParameters *pTweak, bool verbose) {
  this->pTweak = pTweak;

	//commented because maybe we only need methods that don't need tweakable parameters...
  //if (this->pTweak == NULL) {
  //  REPORT_ERROR(SVLIB_BadArg, "Need tweakable parameters");
  //}
  
	this->featureCount = 22; //Bing   // nan   +1 //jun +2 // basti +2 // Bing +1 for SC_FEATURE_SDP
	this->verbose = verbose;
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_FeatureHandler::~SC_FeatureHandler() {

}

//====================================================================================================================
//	return a (new) string containing the feature's name for the given feature-nr (a sclib::feature* constant or a 
//  concatenation of those)
//====================================================================================================================
char* SC_FeatureHandler::getFeatureName(unsigned long int featureNr) {
  char *name = new char[sclib::bufferSize*4];
	unsigned long int singleFeature = 1;

	sprintf(name, "");

	for (unsigned int i = 0; i < this->featureCount; i++) {
		if (sclib::bitTest(featureNr, singleFeature) == true) {
			switch (singleFeature) {
				case sclib::featureBandPeriodicity: {
					sprintf(name, "%s%s|", name, "BandPeriodicity");
					break;
				}
				case sclib::featureBrightnessBandwidth: {
					sprintf(name, "%s%s|", name, "BrightnessAndBandwidth");
					break;
				}
				case sclib::featureFbE: {
					sprintf(name, "%s%s|", name, "FilterBankEnergy");
					break;
				}
				case sclib::featureMFCC: {
					sprintf(name, "%s%s|", name, "MFCC");
					break;
				}
				case sclib::featureNFR: {
					sprintf(name, "%s%s|", name, "NoiseFrameRatio");
					break;
				}
				case sclib::featureSpectrum: {
					sprintf(name, "%s%s|", name, "Spectrum");
					break;
				}
				case sclib::featureSpectrumFlux: {
					sprintf(name, "%s%s|", name, "SpectrumFlux");
					break;
				}
				case sclib::featureSubbandPower: {
					sprintf(name, "%s%s|", name, "SubbandPower");
					break;
				}
				case sclib::featureZCR: {
					sprintf(name, "%s%s|", name, "ZCR");
					break;
				}
				case sclib::featureSTE: {
					sprintf(name, "%s%s|", name, "ShortTimeEnergy");
					break;
				}
				case sclib::featureSDP: { // by Bing
					sprintf(name, "%s%s|", name, "SDP");
					break;
				}
				case sclib::featureLPC: {
					sprintf(name, "%s%s|", name, "LPC");
					break;
				}
				case sclib::featureLPCresidual: {
					sprintf(name, "%s%s|", name, "LPCresidual");
					break;
				}
				case sclib::featurePitch: {
					sprintf(name, "%s%s|", name, "Pitch");
					break;
				}
				case sclib::featureLSP: {
					sprintf(name, "%s%s|", name, "LSP");
					break;
				}
				case sclib::featureFormant: {
					sprintf(name, "%s%s|", name, "Formant");
					break;
				}
				case sclib::featureSamples: {
					sprintf(name, "%s%s|", name, "Samples");
					break;
				}
				default: {
					sprintf(name, "%s%s|", name, "");
					break;
				}
			} //switch feature
		} //if feature included
		singleFeature *= 2;
	} //for all features

	name[strlen(name)-1] = '\0'; //remove last '|'

  return name;
}

//====================================================================================================================
//	This method is called by extractFeatured() to extract a single feature with a given extractor and care for proper
//  debug output and progress report
//====================================================================================================================
SV_Data* SC_FeatureHandler::extractFeature(SC_Corpus *pCorpus, SC_Signal *pSignal, SV_Feature* pExtractor, SC_TweakableParameters::SC_FeaturePar *pParameters, unsigned long int segmentStart, unsigned long int segmentEnd) {
	double lf = 0.0, hf = 0.0;
	SV_Data *pFeature = NULL;
	char *fileName = NULL, *featureName = NULL;
	SC_SignalHandler sigHandler(this->pTweak);
	SC_Signal *pFinalSignal = NULL;

	//show progress, if wished
	/*if (this->verbose == true) {
		char *tmp = getFeatureName(pParameters->featureType);
		printf("(%s)", tmp);
		MFree_1D(tmp);
	}*/

	//check sampleRate, do conversion, if needed
	if (pParameters != NULL && sclib::round(pParameters->sampleRate) != pSignal->SigPar.SRate && sclib::round(pParameters->sampleRate) > 0.0) {
		pFinalSignal = sigHandler.resample(pSignal, pParameters->sampleRate, true);
	} else {
		pFinalSignal = pSignal;
	}

	//normalize cutoff frequencies
	if (pParameters != NULL) {
		lf = (pParameters->lowCut / (double)(pFinalSignal->SigPar.SRate)) * 2.0;
		hf = (pParameters->highCut / (double)(pFinalSignal->SigPar.SRate)) * 2.0; 
	}

	//filter signal and extract feature
	pExtractor->CopySignal(*pFinalSignal);
	if (pFinalSignal != pSignal) {
		MFree_0D(pFinalSignal);
	}
	if (lf > 0.0 || (hf > 0.0 && hf < 1.0)) { //only filter if there is something to do (both hf&lf==0 have a special meaning: do nothing!)
		pExtractor->FIR_filtering(lf, hf);
	}
	pFeature = pExtractor->ExtractFeature();

	//do debug output, if wished
	if (pFeature != NULL && this->pTweak->debug.debugMode & sclib::dbFeatures) {
		featureName = getFeatureName(sclib::bit(pFeature->Hdr.ID));
		fileName = sclib::addPostfixToFilename("feature.txt", featureName);
		pCorpus->featureOut(fileName, featureName, segmentStart, segmentEnd, pFeature, true);
		MFree_1D(featureName);
		MFree_1D(fileName);
	}

	//show progress, if wished
	if (this->verbose == true) {
		printf(".");
	}

	return pFeature;
}

//====================================================================================================================
//  this method receives an already allocated pFeatures array of proper (this->featureCount) size and fills it with 
//  the respective selected feature-sets according to the parameters in pTweak->feature*. if an array-entry already 
//  includes a feature-set, the new one is added to the end of the linked list it represents.
//====================================================================================================================
void	SC_FeatureHandler::extractFeatures(SV_Data** &pFeatures, SC_Corpus *pCorpus, SC_Signal *pSignal, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int featureTypes) {
	int idx, frameSize, frameStep;
	SV_Feature* pExtractor = NULL;
	SV_DataIO IO;
  double actualNormLength, expectedNormLength;
	SV_Signal *pTemp = NULL;
	SC_SignalHandler *pSigHandler = NULL;
	SV_Data *pFeature;

  actualNormLength = (double)(pSignal->GetLen()) / (double)(pSignal->SigPar.SRate);
  expectedNormLength = (double)(segmentEnd-segmentStart+1) / (double)(pCorpus->getGT()->getSignalPrototype()->SigPar.SRate);
  if (actualNormLength < expectedNormLength) {REPORT_ERROR(SVLIB_BadData, "\nSignal needs to be loaded before feature-extraction!\n");}

  //extract MFCCs
  idx = sclib::bitPosition(sclib::featureMFCC);
	if (sclib::bitTest(featureTypes, sclib::featureMFCC) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureMfcc.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureMfcc.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
		pExtractor = new SC_Feature_MFCC(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureMfcc.MFCCorder, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.coeffSelection, this->pTweak->featureMfcc.dEnergy, this->pTweak->featureMfcc.CMN, this->pTweak->featureMfcc.addDeltas, this->pTweak->featureMfcc.addDeltaDeltas, this->pTweak->featureMfcc.method, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_smoothing, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureMfcc, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureMfcc.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }
  
  //Extract Log-Filterbank Energy
  idx = sclib::bitPosition(sclib::featureFbE);
  if (sclib::bitTest(featureTypes, sclib::featureFbE) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureFbe.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureFbe.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_FbE(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureFbe.filterBankSize, this->pTweak->featureFbe.FFTsize, this->pTweak->featureFbe.window, this->pTweak->featureFbe.preEmphasizeFactor, this->pTweak->featureFbe.minFilterBankFrequency, this->pTweak->featureFbe.maxFilterBankFrequency, this->pTweak->featureFbe.MFCCorder, this->pTweak->featureFbe.dEnergy, this->pTweak->featureFbe.frequencyScale, this->pTweak->featureFbe.resultType, this->pTweak->featureFbe.smoothing, 1.0);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureFbe, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureFbe.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract powerSpectrum
  idx = sclib::bitPosition(sclib::featureSpectrum);
  if (sclib::bitTest(featureTypes, sclib::featureSpectrum) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSpectrum.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSpectrum.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_Spectrum(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureSpectrum.preEmphasizeFactor, this->pTweak->featureSpectrum.FFTsize, this->pTweak->featureSpectrum.window, this->pTweak->featureSpectrum.logarithmize, this->pTweak->featureSpectrum.createPhase);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureSpectrum, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureSpectrum.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract Band-Periodicity
  idx = sclib::bitPosition(sclib::featureBandPeriodicity);
  if (sclib::bitTest(featureTypes, sclib::featureBandPeriodicity) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureBandPeriodicity.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureBandPeriodicity.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_BandPeriodicity(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureBandPeriodicity.preEmphasizeFactor, this->pTweak->featureBandPeriodicity.FFTsize, this->pTweak->featureBandPeriodicity.window);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureBandPeriodicity, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureBandPeriodicity.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract Brightness&Bandwidth
  idx = sclib::bitPosition(sclib::featureBrightnessBandwidth);
  if (sclib::bitTest(featureTypes, sclib::featureBrightnessBandwidth) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureBrightnessBandwidth.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureBrightnessBandwidth.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_BrightnessBandwidth(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureBrightnessBandwidth.FFTsize, this->pTweak->featureBrightnessBandwidth.window, this->pTweak->featureBrightnessBandwidth.preEmphasizeFactor);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureBrightnessBandwidth, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureBrightnessBandwidth.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract Noise Frame (Ratio)
  idx = sclib::bitPosition(sclib::featureNFR);
  if (sclib::bitTest(featureTypes, sclib::featureNFR) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureNfr.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureNfr.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_NFR(pSignal->SigPar.SRate, frameSize, this->pTweak->featureNfr.preEmphasizeFactor, this->pTweak->featureNfr.NFRthreshold, this->pTweak->featureNfr.FFTsize, this->pTweak->featureNfr.window);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureNfr, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureNfr.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract Spectrum Flux
  idx = sclib::bitPosition(sclib::featureSpectrumFlux);
  if (sclib::bitTest(featureTypes, sclib::featureSpectrumFlux) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSpectrumFlux.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSpectrumFlux.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_SpectrumFlux(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureSpectrumFlux.FFTsize, this->pTweak->featureSpectrumFlux.window, this->pTweak->featureSpectrumFlux.preEmphasizeFactor);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureSpectrumFlux, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureSpectrumFlux.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract Subband Power
  idx = sclib::bitPosition(sclib::featureSubbandPower);
  if (sclib::bitTest(featureTypes, sclib::featureSubbandPower) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSubBandPower.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSubBandPower.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_SubBandPower(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureSubBandPower.FFTsize, this->pTweak->featureSubBandPower.window, this->pTweak->featureSubBandPower.preEmphasizeFactor);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureSubBandPower, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureSubBandPower.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //Extract Zero Crossing Rate
  idx = sclib::bitPosition(sclib::featureZCR);
  if (sclib::bitTest(featureTypes, sclib::featureZCR) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureZcr.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
	  frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureZcr.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
	  pExtractor = new SC_Feature_ZCR(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureZcr.preEmphasizeFactor, this->pTweak->featureZcr.useChebyshev, this->pTweak->featureZcr.scaleResult);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureZcr, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureZcr.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }
	
  //by Bing
	//Extract STE
	idx = sclib::bitPosition(sclib::featureSTE);
  if (sclib::bitTest(featureTypes, sclib::featureSTE) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSte.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
		frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSte.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
		pExtractor = new SC_Feature_STE(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureSte.preEmphasizeFactor, this->pTweak->featureSte.useButterworth, this->pTweak->featureSte.scaleResult);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureSte, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureSte.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }
  
  //by Bing
  //Extract sclib::featureSDP
	idx = sclib::bitPosition(sclib::featureSDP);
  if (sclib::bitTest(featureTypes,  sclib::featureSDP) == true) {
		frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSdp.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
		frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSdp.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
		pExtractor = new SC_Feature_SDP(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureSdp.preEmphasizeFactor, this->pTweak->featureSdp.m, this->pTweak->featureSdp.lag, this->pTweak->featureSdp.color, this->pTweak->featureSdp.n, this->pTweak->featureSdp.pictureSize, this->pTweak->featureSdp.tau);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureSdp, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureSdp.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //by Jun
	//Extract LPC
	idx = sclib::bitPosition(sclib::featureLPC);
  if (sclib::bitTest(featureTypes, sclib::featureLPC) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureLpc.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples //by thilo
    frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureLpc.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate); //by thilo
		pExtractor = new SC_Feature_LPC(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureLpc.window, this->pTweak->featureLpc.preEmphasizeFactor, this->pTweak->featureLpc.LPCorder, this->pTweak->featureLpc.computeGain);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureLpc, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureLpc.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

  //by Jun
  //Extract LPCresidual
  idx = sclib::bitPosition(sclib::featureLPCresidual);
  if (sclib::bitTest(featureTypes, sclib::featureLPCresidual) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureLpcResidual.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples //by thilo
		frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureLpcResidual.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate); //by thilo
   	pExtractor = new SC_Feature_LPCresidual(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureLpcResidual.preEmphasizeFactor, this->pTweak->featureLpcResidual.order);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureLpcResidual, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureLpcResidual.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

	//extract Pitch information per frame
  //can also be done later as a side product of the v/uv detection algorithms
	idx = sclib::bitPosition(sclib::featurePitch);
  if (sclib::bitTest(featureTypes, sclib::featurePitch) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featurePitch.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
    frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featurePitch.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
    pExtractor = new SC_Feature_Pitch(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featurePitch.method, this->pTweak->featurePitch.sqrt, this->pTweak->featurePitch.esps_cand_thresh, this->pTweak->featurePitch.esps_lag_weight, this->pTweak->featurePitch.esps_freq_weight, this->pTweak->featurePitch.esps_trans_cost, this->pTweak->featurePitch.esps_trans_amp, this->pTweak->featurePitch.esps_trans_spec, this->pTweak->featurePitch.esps_voice_bias, this->pTweak->featurePitch.esps_double_cost, this->pTweak->featurePitch.esps_min_f0, this->pTweak->featurePitch.esps_max_f0, this->pTweak->featurePitch.esps_n_cands, this->pTweak->featurePitch.esps_wind_dur, false); //don't make the pitch extractor verbose, it talks too much...
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featurePitch, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featurePitch.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

	//Extract LSPs
	idx = sclib::bitPosition(sclib::featureLSP);
  if (sclib::bitTest(featureTypes, sclib::featureLSP) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureLsp.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
    frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureLsp.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
    pExtractor = new SC_Feature_LSP(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureLsp.window, this->pTweak->featureLsp.preEmphasizeFactor, this->pTweak->featureLsp.LPCorder, this->pTweak->featureLsp.method, this->pTweak->featureLsp.delta, this->pTweak->featureLsp.bisections, this->pTweak->featureLsp.minSeparation, this->pTweak->featureLsp.maxLoops);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureLsp, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureLsp.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
	}

	//extract Formant information per frame
	idx = sclib::bitPosition(sclib::featureFormant);
  if (sclib::bitTest(featureTypes, sclib::featureFormant) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureFormant.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
    frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureFormant.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
    pExtractor = new SC_Feature_Formant(pSignal->SigPar.SRate, frameSize, frameStep, this->pTweak->featureFormant.esps_lpc_ord, this->pTweak->featureFormant.esps_lpc_type, this->pTweak->featureFormant.esps_w_type, this->pTweak->featureFormant.esps_ds_freq, this->pTweak->featureFormant.esps_wdur, this->pTweak->featureFormant.esps_nom_f1, this->pTweak->featureFormant.esps_cor_wdur, this->pTweak->featureFormant.esps_frame_int, this->pTweak->featureFormant.preEmphasizeFactor, this->pTweak->featureFormant.esps_nform, false); //the formant tracker is too verbose...
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureFormant, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureFormant.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

	//extract samples per frame
	idx = sclib::bitPosition(sclib::featureSamples);
  if (sclib::bitTest(featureTypes, sclib::featureSamples) == true) {
    frameSize = pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSamples.frameSize, sclib::alignmentStart, pSignal->SigPar.SRate); //ms -> samples
    frameStep	= pCorpus->getGT()->getConverter()->ms2sample(this->pTweak->featureSamples.frameStep, sclib::alignmentStart, pSignal->SigPar.SRate);
    pExtractor = new SC_Feature_Samples(pSignal->SigPar.SRate, frameSize, frameStep);
		pFeature = extractFeature(pCorpus, pSignal, pExtractor, &this->pTweak->featureSamples, segmentStart, segmentEnd);
		if (pFeature != NULL) {
			pFeature->Hdr.Signature[0] = this->pTweak->featureSamples.client;
			sclib::addToList(pFeatures[idx], pFeature);
		}
		MFree_0D(pExtractor);
  }

	return;
}

//====================================================================================================================
//	Method to control feature-extraction procedure: An this->featureCount-dimensional array of feature-vector classes 
//  is returned, whith the following content:
//    - pFeatures[0]:  MFCCs or NULL
//    - pFeatures[1]:  NULL (was: Energy&ZCR)
//    - pFeatures[2]:  FbE's or NULL
//    - pFeatures[3]:  Spectrum or NULL
//    - pFeatures[4]:  Band Periodicity or NULL
//    - pFeatures[5]:  Brightness&Bandwidth or NULL
//    - pFeatures[6]:  NFR or NULL
//    - pFeatures[7]:  Spectrum Flux or NULL
//    - pFeatures[8]:  Subband Power or NULL
//    - pFeatures[9]:  ZCR or NULL
//    - pFeatures[10]: STE or NULL                  //Bing
//		- pFeatures[11]: CepstralPeak or NULL         //Bing
//    - pFeatures[12]: Wavelet Energy Distribution  //Nan
//    - pFeatures[13]: LPC or NULL;                 //Jun
//    - pFeatures[14]: LPCresidual or NULL          //Jun
//    - pFeatures[15]: AAP or NULL                  //Basti
//    - pFeatures[16]: BFDAC or NULL                //Basti
//    - pFeatures[17]: SDP or NULL                  //Bing
//    - pFeatures[18]: Pitch or NULL
//    - pFeatures[19]: LSPs or NULL
//    - pFeatures[20]: Formants or NULL
//    - pFeatures[21]: Samples or NULL
//  NULL means that this class of feature-vectors wasn't selected by the featureTypes parameter; each entry might 
//  include a linked list of features of the same class, but extracted using different sets of parameters. in the 
//  standard case, only the features for the standard parameter sets (pTweak->feature*) are extracted, but additional 
//  parametersets can be given in the linked list pSpecialParameters.
//====================================================================================================================
SV_Data** SC_FeatureHandler::extractFeatures(SC_Corpus *pCorpus, SC_Signal *pSignal, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int featureTypes, SC_TweakableParameters::SC_FeaturePar *pSpecialParameters, bool extractOnlySpecialFeatures) {
	SV_Data **pFeatures = NULL;
	unsigned long int i;
	SC_TweakableParameters::SC_FeaturePar **pOriginalParameters = NULL, *pParameterHook;

  MArray_1D(pFeatures, this->featureCount, SV_Data*, "SC_FeatureHandler.extractFeatures: pFeatures");
	MArray_1D(pOriginalParameters, this->featureCount, SC_TweakableParameters::SC_FeaturePar*, "SC_FeatureHandler.extractFeatures: pOriginalParameters");
	for (i = 0; i < this->featureCount; i++) {
		pFeatures[i] = NULL;
		pOriginalParameters[i] = NULL;
	}

	//1. extract all the wanted features using standard parameters
	if (extractOnlySpecialFeatures == false || pSpecialParameters == NULL) {
		extractFeatures(pFeatures, pCorpus, pSignal, segmentStart, segmentEnd, featureTypes);
	}

	//2. care for each special parameter set separately, if there are any; only extract featuresets a 2nd time if parameters differ (cared for in extractWithSpecialParameters()!)
	pParameterHook = pSpecialParameters;
	while (pParameterHook != NULL) {
		if (sclib::bitTest(featureTypes, pParameterHook->featureType) == true) { //only extract features that are wanted (normally, providing a special parameter set and wanting that feature should coincide, but who knows...)
			switch (pParameterHook->featureType) {
				case sclib::featureBandPeriodicity:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureBandPeriodicityPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureBandPeriodicityPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureBandPeriodicity), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureBrightnessBandwidth:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureBrightnessBandwidthPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureBrightnessBandwidthPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureBrightnessBandwidth), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureFbE:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureFbEPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureFbEPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureFbe), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureLPC:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureLPCPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureLPCPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureLpc), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureLPCresidual:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureLPCresidualPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureLPCresidualPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureLpcResidual), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureLSP:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureLSPPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureLSPPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureLsp), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureMFCC:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureMFCCPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureMFCCPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureMfcc), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureNFR:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureNFRPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureNFRPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureNfr), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featurePitch:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeaturePitchPar*)(pParameterHook), (SC_TweakableParameters::SC_FeaturePitchPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featurePitch), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureSDP:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureSDPPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureSDPPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureSdp), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureSpectrum:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureSpectrumPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureSpectrumPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureSpectrum), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureSpectrumFlux:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureSpectrumFluxPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureSpectrumFluxPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureSpectrumFlux), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureSTE:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureSTEPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureSTEPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureSte), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureSubbandPower:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureSubBandPowerPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureSubBandPowerPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureSubBandPower), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureZCR:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureZCRPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureZCRPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureZcr), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureFormant:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureFormantPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureFormantPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureFormant), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				case sclib::featureSamples:
					extractWithSpecialParameters(pFeatures, pParameterHook->featureType, (SC_TweakableParameters::SC_FeatureSamplesPar*)(pParameterHook), (SC_TweakableParameters::SC_FeatureSamplesPar*)(pOriginalParameters[sclib::bitPosition(pParameterHook->featureType)]), &(this->pTweak->featureSamples), pCorpus, pSignal, segmentStart, segmentEnd, extractOnlySpecialFeatures);
					break;
				default:
					REPORT_ERROR(SVLIB_BadArg, "Feature type unknown in extraction using special parameter sets");
					break;
			}
		}

		pParameterHook = pParameterHook->Next;
	}

	//clean up
	for (i = 0; i < this->featureCount; i++) {
		MFree_0D(pOriginalParameters[i]);
	}
	MFree_1D(pOriginalParameters);

	return pFeatures;
}

//====================================================================================================================
//	Save an array of (linked lists of) features as returned by extractFeatures()
//====================================================================================================================
bool SC_FeatureHandler::saveFeatures(const char *fileName, SV_Data** pFeatures, unsigned long int featureTypes) {
 	int res;
  bool success = true;
  SV_DataIO IO;
	SV_Data *pHook;
  
  //save feature to file, if filename is != ""
	if (strlen(fileName) > 0) {
    IO.OpenFile(fileName, WRITE_REC);
		for (unsigned long int i = 0; i < this->featureCount; i++) {
			pHook = pFeatures[i];
			while (pHook != NULL) {
		    res = IO.PutDataRec(*pHook);
				if (res == 0) {
          success = false;
        }
				pHook = pHook->Next;
      }
		}
    IO.CloseFile();
	}

  return success;
}

//====================================================================================================================
//	Save a single feature-set (this function is able to append data to an existing feature-file); linked lists of 
//  feature-sets are not supported
//  appends if the file already exists
//====================================================================================================================
bool SC_FeatureHandler::saveFeature(const char *fileName, SV_Data *pFeature) {
 	int res;
  bool success = true;
  SV_DataIO IO;
  
  //save feature to file, if filename is != ""
	if (strlen(fileName) > 0) {
    IO.OpenFile(fileName, WRITE_REC);
		res = IO.PutDataRec(*pFeature);
    if (res == 0) {
      success = false;
    }
    IO.CloseFile();
  }

  return success;  
}

//====================================================================================================================
//	Save a linked list of features to a file 
//====================================================================================================================
bool SC_FeatureHandler::saveFeatureList(const char *fileName, SV_Data *pList) {
 	int res;
  bool success = true;
  SV_DataIO IO;
  SV_Data *pHook = pList;
  
  //save feature to file, if filename is != ""
	if (strlen(fileName) > 0) {
    IO.OpenFile(fileName, WRITE_REC);
    while (pHook != NULL) {
      res = IO.PutDataRec(*pHook);
      if (res == 0) {
        success = false;
        break;
      }
      pHook = pHook->Next;
    }
    IO.CloseFile();
  }

  return success;  
}

//====================================================================================================================
//	Load all features from the file and return them in one big (merged) SV_Data object instead as a linked list
//====================================================================================================================
SV_Data* SC_FeatureHandler::loadAndMergeFeatures(const char *fileName, long int firstCol, long int lastCol, bool pivot) {
  long int rows, x, y, col, firstColumn, lastColumn;
  int cols;
  SV_DataIO IO;
  SV_Data *pData = NULL, *pTemp;

  if (sclib::fileExists(fileName) == true) {
    IO.OpenFile(fileName, READ_REC);
		rows = IO.getAllRecRowCount(cols);
    if (cols > 0 && rows > 0) {
      firstColumn = (firstCol >= 0 && firstCol < cols) ? firstCol : 0;
      lastColumn = (lastCol >= 0 && lastCol < cols) ? lastCol : cols-1;
      cols = lastColumn - firstColumn + 1;
			if (pivot == false) {
				pData = new SV_Data(rows, cols);
			} else {
				pData = new SV_Data(cols, rows);
			}
      rows = 0;
      pTemp = IO.GetDataRec();
      while (pTemp != NULL) {
        if (rows == 0) { //copy header
          pData->Hdr.ByteOrder = pTemp->Hdr.ByteOrder;
          pData->Hdr.frameSize = pTemp->Hdr.frameSize;
          pData->Hdr.frameStep = pTemp->Hdr.frameStep;
          pData->Hdr.ID = pTemp->Hdr.ID;
        } else if (pTemp->Hdr.ByteOrder != pData->Hdr.ByteOrder || //check header
                   pTemp->Hdr.frameSize != pData->Hdr.frameSize ||
                   pTemp->Hdr.frameStep != pData->Hdr.frameStep ||
                   pTemp->Hdr.ID != pData->Hdr.ID) {
          REPORT_ERROR(SVLIB_BadData, "Merging different Features is undefined!");
        }

        for (y = 0; y < pTemp->Row; y++) {
					col = 0;
          for (x = firstColumn; x <= lastColumn; x++) {
						if (pivot == false) {
							pData->Mat[rows + y][col++] = pTemp->Mat[y][x];
						} else {
							pData->Mat[col++][rows + y] = pTemp->Mat[y][x];
						}
          }
        }
        rows += pTemp->Row;
        MFree_0D(pTemp);
        pTemp = IO.GetDataRec();
      }
    }
    IO.CloseFile();
  }

  return pData;
}

//====================================================================================================================
//	Load all features from the file and return them in one big (merged) SVMproblem object instead as a linked list
//  that can be used for one-class SVM training (i.e. labels always positive, no weights)
//====================================================================================================================
SC_SVMproblem* SC_FeatureHandler::loadAndMergeSvmFeatures(const char *fileName, long int firstCol, long int lastCol) {
  long int rows, x, y, firstColumn, lastColumn;
  int cols;
  SV_DataIO IO;
	SV_Data *pTemp, *pTemplate;
  SC_SVMproblem *pProblem = new SC_SVMproblem;
	SC_SVMnode *pNodes = NULL;

  if (sclib::fileExists(fileName) == true) {
    IO.OpenFile(fileName, READ_REC);
		rows = IO.getAllRecRowCount(cols);
    if (cols > 0 && rows > 0) {
      firstColumn = (firstCol >= 0 && firstCol < cols) ? firstCol : 0;
      lastColumn = (lastCol >= 0 && lastCol < cols) ? lastCol : cols-1;
      cols = lastColumn - firstColumn + 1;
      pProblem->l = rows;
			MArray_1D(pProblem->y, pProblem->l, double, "SC_FeatureHandler.loadAndMergeSvmFeatures: pProblem->y");
			MArray_1D(pProblem->x, pProblem->l, SC_SVMnode*, "SC_FeatureHandler.loadAndMergeSvmFeatures: pProblem->x");
			pProblem->W = NULL;
			pTemplate = new SV_Data();
	    rows = 0;
      pTemp = IO.GetDataRec();
      while (pTemp != NULL) {
        if (rows == 0) { //copy header
          pTemplate->Hdr.ByteOrder = pTemp->Hdr.ByteOrder;
          pTemplate->Hdr.frameSize = pTemp->Hdr.frameSize;
          pTemplate->Hdr.frameStep = pTemp->Hdr.frameStep;
          pTemplate->Hdr.ID = pTemp->Hdr.ID;
        } else if (pTemp->Hdr.ByteOrder != pTemplate->Hdr.ByteOrder || //check header
                   pTemp->Hdr.frameSize != pTemplate->Hdr.frameSize ||
                   pTemp->Hdr.frameStep != pTemplate->Hdr.frameStep ||
                   pTemp->Hdr.ID != pTemplate->Hdr.ID) {
          REPORT_ERROR(SVLIB_BadData, "Merging different Features is undefined!");
        }

        for (y = 0; y < pTemp->Row; y++) {
					pProblem->y[rows+y] = (double)(sclib::labelPositive);
				  MArray_1D(pNodes, cols+1, SC_SVMnode, "SC_FeatureHandler.loadAndMergeSvmFeatures: pNodes");
					for (x = 0; x < cols; x++) {
						pNodes[x].index = x + 1;
						pNodes[x].value = pTemp->Mat[y][x+firstColumn];
					}
					pNodes[cols].index = -1; //this tells the SVM that the vector is finished
					pProblem->x[rows+y] = pNodes;
        }
        rows += pTemp->Row;
        MFree_0D(pTemp);
        pTemp = IO.GetDataRec();
      }
			MFree_0D(pTemplate);
		}
    IO.CloseFile();
  }

  return pProblem;
}

//====================================================================================================================
//	Merges a linked list of feature-vectors to one big matrix not in main memory but on hard-disk using a temporary
//  file; the original list is replaced by the merged dataset to save memory (at all time there is only one version
//  of the data in memory)
//  if a filename is given explicitly, the data is stored therein and NOT deleted after merging (debug-dir isn't used)
//====================================================================================================================
bool SC_FeatureHandler::mergeOnDisk(SV_Data* &pDataList, const char *fileName) {
	bool res = true;
	const char tmpName[] = "featureList.tmp\0";
	char tmp[sclib::bufferSize];

	if (fileName == NULL || strncmp(fileName, "", 1) == 0) {
		sprintf(tmp, "%s%s", this->pTweak->debug.debugDir, tmpName); //use temporary file
	} else {
		sprintf(tmp, "%s", fileName); //use given fileName
	}

	res = saveFeatureList(tmp, pDataList);

	if (res != false) {
		sclib::destructLinkedList(pDataList);
		pDataList = loadAndMergeFeatures(tmp);
		if (pDataList == NULL) {
			res = false;
		}
	}

	if (fileName == NULL || strncmp(fileName, "", 1) == 0) { //only delete the file if internal temporary filename was used
		if (remove(tmp) != 0) {
			res = false;	
		}
	}

	return res;
}

//====================================================================================================================
//	Load an array of features as returned by extractFeatures(); if there are more than one dataset of each featureType
//  in the file, the corresponding entry in the feature-array is a linked list
//====================================================================================================================
SV_Data** SC_FeatureHandler::loadFeatures(const char *fileName, unsigned long int featureTypes) {
	SV_DataIO IO;
  SV_Data** pFeatures = NULL;
	unsigned long int feature;

  MArray_1D(pFeatures, this->featureCount, SV_Data*, "SC_FeatureHandler.extractFeatures: pFeatures");

  //if exists, read features from file
  if (sclib::fileExists(fileName) == true) {
		for (unsigned long int i = 0; i < this->featureCount; i++) {
			feature = sclib::bit(i);
			if (featureTypes == 0 || sclib::bitTest(featureTypes, feature) == true) {
				IO.OpenFile(fileName, READ_REC); //always open&close because it is read till the end in order to get all parts of the possibly linked list for this id!
				pFeatures[i] = IO.GetAllRec((int)(feature)); //TODO: cast to int save?
				IO.CloseFile();
			} else {
				pFeatures[i] = NULL;
			}
    }
  }
	
  return pFeatures;
}

//====================================================================================================================
//	Load a single feature, if featureType > 0 (if there are more than 1 entries with given featureType, a linked list 
//  is returned; if featureType == 0, return a linked list of all featuresets in the file
//====================================================================================================================
SV_Data* SC_FeatureHandler::loadFeature(const char *fileName, unsigned long int featureType) {
	SV_DataIO IO;
  SV_Data* pFeature = NULL;
	unsigned long int feature;

	if (!sclib::isPowerOfTwo(featureType)) {
		REPORT_ERROR(SVLIB_BadArg, "Only a single feature (or all fatures by specifiying 0) can be loaded here");
	}

  //if exists, read features from file
  if (sclib::fileExists(fileName) == true) {
    IO.OpenFile(fileName, READ_REC);
    if (featureType == 0) { //load all features in this file, save as a linked list
        pFeature = IO.GetAllRec();
    } else { //load only the specified feature
		  for (unsigned long int i = 0; i < this->featureCount; i++) {
				feature = sclib::bit(i);
				if (featureType == feature) {
					pFeature = IO.GetAllRec((int)(feature)); //TODO: is cast to int save here?
					break; 
				}
      }
    }
    IO.CloseFile();
  }
	
  return pFeature;
}

//====================================================================================================================
//	Normalize a set of feature-vectors by rescaling each vector x to the interval [0..1] by using the min/max values 
//  for each dimension (best obtained from a big set of training-data): x' = (x-min)/(max-min)
//====================================================================================================================
SV_Data* SC_FeatureHandler::normalize(SV_Data *pFeatures, double *min, double *max, bool replace) {
	double *spread;
	unsigned long int t, d;
	SV_Data *pRes = pFeatures;

	MArray_1D(spread, pFeatures->Col, double, "SC_FeatureHandler.normalize: spread");
	for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
		spread[d] = max[d] - min[d];
	}

	if (replace == false) {
		pRes = new SV_Data(pFeatures->Row, pFeatures->Col);
		pRes->Hdr = pFeatures->Hdr;
	}

	for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
    for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
			pRes->Mat[t][d] = (float)(sclib::getBetween(min[d], (pFeatures->Mat[t][d]-min[d])/spread[d], max[d])); //scale to [0..1]  
    } 
  }

	MFree_1D(spread);

  return pRes;
}

//====================================================================================================================
//	Does normalization as above; the min thereby comes from the first row of pNorm while the max comes from the second 
//  row
//====================================================================================================================
SV_Data* SC_FeatureHandler::normalize(SV_Data *pFeatures, SV_Data *pNorm, bool replace) {
	double *spread;
	unsigned long int t, d;
	SV_Data *pRes = pFeatures;

	MArray_1D(spread, pFeatures->Col, double, "SC_FeatureHandler.normalize: spread");
	for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
		spread[d] = pNorm->Mat[1][d] - pNorm->Mat[0][d];
	}
	
	if (replace == false) {
		pRes = new SV_Data(pFeatures->Row, pFeatures->Col);
		pRes->Hdr = pFeatures->Hdr;
	}

	for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
    for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
			pRes->Mat[t][d] = sclib::getBetween(pNorm->Mat[0][d], (pFeatures->Mat[t][d]-pNorm->Mat[0][d])/spread[d], pNorm->Mat[1][d]); //scale to [0..1]  
    }
  }

	MFree_1D(spread);

  return pRes;
}

//====================================================================================================================
//	normalizes the length of each column vector to unity (1)
//====================================================================================================================
void SC_FeatureHandler::normalizeColumns(SV_Data *pFeatures) {
	unsigned long int t, d;
	double sum;

  for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
		sum = 0.0;
		for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
			sum += pFeatures->Mat[t][d];
    }

		for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
			pFeatures->Mat[t][d] = (float)(pFeatures->Mat[t][d] / sum);
    }
	}

	return;
}

//====================================================================================================================
//	normalizes the length of each row vector to unity (1)
//====================================================================================================================
void SC_FeatureHandler::normalizeRows(SV_Data *pFeatures) {
	unsigned long int t, d;
	double sum;

  for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
		sum = 0.0;
		for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
			sum += pFeatures->Mat[t][d];
    }

		for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
			pFeatures->Mat[t][d] = (float)(pFeatures->Mat[t][d] / sum);
    }
	}

	return;
}

//====================================================================================================================
//	Removes the effect of normalization
//====================================================================================================================
SV_Data* SC_FeatureHandler::unNormalize(SV_Data *pFeatures, SV_Data *pNorm, bool replace) {
	double *spread;
	unsigned long int t, d;
	SV_Data *pRes = pFeatures;

	MArray_1D(spread, pFeatures->Col, double, "SC_FeatureHandler.unNormalize: spread");
	for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
		spread[d] = pNorm->Mat[1][d] - pNorm->Mat[0][d];
	}
	
	if (replace == false) {
		pRes = new SV_Data(pFeatures->Row, pFeatures->Col);
		pRes->Hdr = pFeatures->Hdr;
	}

	for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
    for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
			pRes->Mat[t][d] = (float)((pFeatures->Mat[t][d]*spread[d]) + pNorm->Mat[0][d]);
    }
  }

	MFree_1D(spread);

  return pRes;
}

//====================================================================================================================
//	Create a normalization matrix as expected by normalize()
//====================================================================================================================
SV_Data* SC_FeatureHandler::createNormalizationMatrix(SV_Data *pFeatures) {
	double *min, *max;
	SV_Data *pNorm = NULL;
  SC_MatrixFunctions *pFunc = NULL;
	
	if (pFeatures != NULL && pFeatures->Col > 0 && pFeatures->Row > 0) {
		pNorm = new SV_Data(2, pFeatures->Col);
		pFunc = new SC_MatrixFunctions();

		min = pFunc->min(pFeatures->Mat, pFeatures->Row, pFeatures->Col);
		max = pFunc->max(pFeatures->Mat, pFeatures->Row, pFeatures->Col);
		
		for (int x = 0; x < pNorm->Col; x++) {
			pNorm->Mat[0][x] = (float)(min[x]);
			pNorm->Mat[1][x] = (float)(max[x]);
		}

		MFree_1D(min);
		MFree_1D(max);
	  MFree_0D(pFunc);
	}

	return pNorm;
}

//====================================================================================================================
//	Standardize a set of feature-vectors by rescaling each vector x by the mean m and standard deviation sd of this 
//  feature's distribution (best obtained from a big set of training-data): x' = (x-m)/sd
//====================================================================================================================
SV_Data* SC_FeatureHandler::standardize(SV_Data *pFeatures, double *mean, double *sd, bool replace) {
	SV_Data *pRes = pFeatures;

	if (replace == false) {
		pRes = new SV_Data(pFeatures->Row, pFeatures->Col);
		pRes->Hdr = pFeatures->Hdr;
	}

	for (unsigned long int t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
    for (unsigned long int d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
      pRes->Mat[t][d] = (float)(sclib::zTransform(pFeatures->Mat[t][d], mean[d], sd[d]));
    }
  }

  return pRes;
}

//====================================================================================================================
//	Does standardization as above; the mean thereby comes from the first row of pNorm while the sd comes from the 
//  second row
//====================================================================================================================
SV_Data* SC_FeatureHandler::standardize(SV_Data *pFeatures, SV_Data *pStd, bool replace) {
	SV_Data *pRes = pFeatures;

	if (replace == false) {
		pRes = new SV_Data(pFeatures->Row, pFeatures->Col);
		pRes->Hdr = pFeatures->Hdr;
	}
  
	for (unsigned long int t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
    for (unsigned long int d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
      pRes->Mat[t][d] = (float)(sclib::zTransform(pFeatures->Mat[t][d], pStd->Mat[0][d], pStd->Mat[1][d]));
    }
  }

  return pRes;
}

//====================================================================================================================
//	Create a standardization matrix as expected by standardize()
//====================================================================================================================
SV_Data* SC_FeatureHandler::createStandardizationMatrix(SV_Data *pFeatures) {
	double *mean, *sd;
	SV_Data *pStd = NULL;
  SC_MatrixFunctions *pFunc = NULL;
	
	if (pFeatures != NULL && pFeatures->Col > 0 && pFeatures->Row > 0) {
		pStd = new SV_Data(2, pFeatures->Col);
		pFunc = new SC_MatrixFunctions();

		mean = pFunc->mean(pFeatures->Mat, pFeatures->Row, pFeatures->Col);
		sd = pFunc->std(pFeatures->Mat, pFeatures->Row, pFeatures->Col, mean);
		
		for (int x = 0; x < pStd->Col; x++) {
			pStd->Mat[0][x] = (float)(mean[x]);
			pStd->Mat[1][x] = (float)(sd[x]);
		}

		MFree_1D(mean);
		MFree_1D(sd);
	  MFree_0D(pFunc);
	}

	return pStd;
}

//====================================================================================================================
//	Concatenate the columns of all non-NULL array-entries of features (normalization should be done before or 
//  afterwards to ensure that the scales of the features are similar; otherwise, concatenation in not meaningful) in 
//  the array going from 0 to featureCount-1 (featureCount is necessary is a parameter to be able to handle feature-
//  sets of other cardinality as returned by extractFeatures(), as is the case in the LZL ATC algorithm)
//  If selectionMap != 0, only those features are concatenated that have their index set as a bit-flag in this map. 
//  If individualColSelectionMap != NULL, for each previously selected (non-NULL, selected in selectionMap or 
//  selectionMap == 0) featureset it is lookad at it's index-position: Only these columns are concatenated which's 
//  bitflag (first features has index 1!) is set to 1
//  The featureMap provides, indexed by the featureType, start- and end-index of the columns of this featureType
//====================================================================================================================
SV_Data* SC_FeatureHandler::combineFeatureVectors(SV_Data** pFeatures, unsigned long int featureCount, unsigned long int selectionMap, unsigned long int *individualColSelectionMap) {
	std::map<unsigned long int, std::pair<unsigned int, unsigned int> > featureMap;
	//this version doesn't return the featureMap
	return combineFeatureVectors(pFeatures, featureCount, featureMap, selectionMap, individualColSelectionMap);
}
SV_Data* SC_FeatureHandler::combineFeatureVectors(SV_Data** pFeatures, unsigned long int featureCount, std::map<unsigned long int, std::pair<unsigned int, unsigned int> > &featureMap, unsigned long int selectionMap, unsigned long int *individualColSelectionMap) {
  unsigned long int colCount = 0, rowCount = 0, t, d, i;
  unsigned int frameSize = 0, frameStep = 0, sampleRate = 0;
  SV_Data *pResult = NULL;
	std::pair<unsigned int, unsigned int> startStop;
	bool first;
	unsigned int firstCol;

  //calculate final number of cols
  for (i = 0; i < featureCount; i++) {
		first = true;
    if (pFeatures[i] != NULL && (selectionMap == 0 || sclib::bitTest(selectionMap, sclib::bit(i)) == true)) {
      if (frameSize == 0 && frameStep == 0 && sampleRate == 0) { //check for equal frame/sampleRate parameters of all non-NULL features-sets
        frameSize = pFeatures[i]->Hdr.frameSize;
        frameStep = pFeatures[i]->Hdr.frameStep;
        sampleRate = pFeatures[i]->Hdr.sampleRate;
      } else if (pFeatures[i]->Hdr.frameSize != frameSize || pFeatures[i]->Hdr.frameStep != frameStep || pFeatures[i]->Hdr.sampleRate != sampleRate) {
        REPORT_ERROR(SVLIB_BadData, "\nFeatures should have equal frame-parameters for combination!\n");
      }
			if (individualColSelectionMap != NULL) { //only concat selected cols per feature-set
				for (d = 1; d <= (unsigned long int)(pFeatures[i]->Col); d++) {
					if (sclib::bitTest(individualColSelectionMap[i], sclib::bit(d)) == true) {
						if (first == true) {
							firstCol = colCount;
							first = false;
						}
						startStop = std::make_pair(firstCol, colCount);
						featureMap[sclib::bit(i)] = startStop; //insert the tupel several times, the last time the correct endCol is inserted; this is not the fasted way possible, but simple to implement and speed isn't critical here...
						colCount++;
					}
				}
			} else { //concat all cols
				startStop = std::make_pair(colCount, colCount+pFeatures[i]->Col-1);
				featureMap[sclib::bit(i)] = startStop;
				colCount += pFeatures[i]->Col;
			}
      if (rowCount == 0) { //check for equal set-length of all non-NULL sets
        rowCount = pFeatures[i]->Row;
      } else if (rowCount != pFeatures[i]->Row) {
        REPORT_ERROR(SVLIB_BadData, "\nFeatures can't be combined if they have differnt row-counts!\n");
      }
    }
  }

  if (rowCount > 0) {
    pResult = new SV_Data(rowCount, colCount);
    pResult->Hdr.frameSize = frameSize; //we checked above that all feature-sets have equal parameters... 
    pResult->Hdr.frameStep = frameStep;
    pResult->Hdr.sampleRate = sampleRate;
		pResult->Hdr.ID = 0;

    //concatenate cols
    for (t = 0; t < rowCount; t++) {
      colCount = 0;
      for (i = 0; i < featureCount; i++) {
        if (pFeatures[i] != NULL && (selectionMap == 0 || sclib::bitTest(selectionMap, sclib::bit(i)) == true)) {
					if (t == 0) {
						pResult->Hdr.ID |= pFeatures[i]->Hdr.ID;
					}
          for (d = 0; d < (unsigned long int)(pFeatures[i]->Col); d++) {
						if (individualColSelectionMap == NULL || sclib::bitTest(individualColSelectionMap[i], sclib::bit(d+1)) == true) {
	            pResult->Mat[t][colCount] = pFeatures[i]->Mat[t][d];
							colCount++;
						}
          }
        }
      }
    }
  }

  return pResult;
}


//====================================================================================================================
//	Can be called prior to combineFeatureVectors to bring all feature sets to a common frame-rate/-size (smallest 
//  frame-rate is used as the common ground)
//====================================================================================================================
void SC_FeatureHandler::equalizeFrameParameters(SV_Data** pFeatures, unsigned long int featureCount, unsigned long int selectionMap) {
	unsigned int i, rows = std::numeric_limits<unsigned int>::max(), idx;
	double iv[1] = {0.0}; //ignore-value-column for pitch feature
	SV_Data *pConverted;

  //find smallest row-count (i.e. smallest frame-rate)
  for (i = 0; i < featureCount; i++) {
    if (pFeatures[i] != NULL && (selectionMap == 0 || sclib::bitTest(selectionMap, sclib::bit(i)) == true)) {
			if ((unsigned int)(pFeatures[i]->Row) < rows) {
				idx = i;
				rows = pFeatures[i]->Row;
			}
    }
  }

	//convert all other selected sets to this frame parameters
  for (i = 0; i < featureCount; i++) {
    if (pFeatures[i] != NULL && (selectionMap == 0 || sclib::bitTest(selectionMap, sclib::bit(i)) == true)) {
			if (pFeatures[i]->Hdr.sampleRate != pFeatures[idx]->Hdr.sampleRate) { //convert sample rate by just simultaneously changing rate- and frame-parameters
				pFeatures[i]->Hdr.frameSize = sclib::round(pFeatures[i]->Hdr.frameSize * (double)(pFeatures[idx]->Hdr.sampleRate)/(double)(pFeatures[i]->Hdr.sampleRate));
				pFeatures[i]->Hdr.frameSize = sclib::round(pFeatures[i]->Hdr.frameStep * (double)(pFeatures[idx]->Hdr.sampleRate)/(double)(pFeatures[i]->Hdr.sampleRate));
				pFeatures[i]->Hdr.sampleRate = pFeatures[idx]->Hdr.sampleRate;
			}
			if (pFeatures[i]->Hdr.frameSize != pFeatures[idx]->Hdr.frameSize || pFeatures[i]->Hdr.frameStep != pFeatures[idx]->Hdr.frameStep) {
				pConverted = convertFrameRate(pFeatures[i], (rows-1)*pFeatures[idx]->Hdr.frameStep+pFeatures[idx]->Hdr.frameSize, pFeatures[idx]->Hdr.frameSize, pFeatures[idx]->Hdr.frameStep, (pFeatures[i]->Hdr.ID==sclib::featurePitch)?iv:NULL);
				pConverted->Next = pFeatures[i]->Next;
				MFree_0D(pFeatures[i]);
				pFeatures[i] = pConverted;
			}
			if (pFeatures[i]->Row != pFeatures[idx]->Row) {
				REPORT_ERROR(SVLIB_Fail, "still unequal row counts after equalization of frame parameters!");
			}
		}
	}

	return;
}

//====================================================================================================================
//	Load features in HTK (Hidden Markov Model Toolkit) format;
//  TODO: make it machine-independant by using SV_DataIO's fgeatures of reading/writing binary data
//        (but Wesley's files aren't written machine-independant anyway...)
//====================================================================================================================
SV_Data* SC_FeatureHandler::loadHTKfeatures(const char *fileName) {
  SC_HTK_HeaderType header;
  FILE *inFile;
  SV_Data *pFeatures = NULL;
  size_t res = 0;
  int x;
  float *row;
	char tmp[sclib::bufferSize];

  //can the file be opened?
  if ((inFile=fopen(fileName,"rb")) == NULL) {
		sprintf(tmp, "Load HTK data failed (could not open file '%s')!", fileName);
    REPORT_ERROR(SVLIB_BadData, tmp);
    return NULL;
  }
  
  //fill header
  res = fread(&(header.sampleCount), sizeof(int), 1, inFile);
  res += fread(&(header.sampleStep), sizeof(int), 1, inFile);
  res += fread(&(header.sampleSize), sizeof(short), 1, inFile);
  res += fread(&(header.sampleType), sizeof(short), 1, inFile);
  header.featureDim = header.sampleSize / sizeof(float); 
  //header.fileName = new char[strlen(fileName) + 1]; //ok, it is freed just below in this function, but if someone else needs a full header, he has a full example of how to fill it here ;-)
  //strcpy(header.fileName, fileName);
  //header.fileName[strlen(fileName)] = '\0';

  //get memory
  pFeatures = new SV_Data(header.sampleCount, header.featureDim);
  MArray_1D(row, header.featureDim, float, "SC_FeatureHandler.loatHTKfeatures: row");

  //can't be determined out of the file's header... must be set later!!!
  //pFeatures->Hdr.frameSize = 
  //pFeatures->Hdr.frameStep = 

  switch (header.sampleType) {
    case 6: {
      pFeatures->Hdr.ID = sclib::featureMFCC;
      break;
    }
    case 7: {
      pFeatures->Hdr.ID = sclib::featureFbE;
      break;
    }
    case 8: {
      pFeatures->Hdr.ID = sclib::featureFbE;
      break;
    }
  }

  //read data
  for (long int y = 0; y < header.sampleCount; y++){
	  res += fread(row, sizeof(float), header.featureDim, inFile);
    for (x = 0; x < header.featureDim; x++) {
      pFeatures->Mat[y][x] = row[x];
    }
  }

  //fre memory
  MFree_1D(row);  
  fclose(inFile);

  if (res != ((header.sampleCount * header.featureDim) + 4)) {
    REPORT_ERROR(SVLIB_FileErr, "Not all bytes could be read from the HTK featurefile");
  }

	return pFeatures;
}

//====================================================================================================================
//	Save features in HTK (Hidden Markov Model Tolkit) format return nr. of bytes written on success, SVLib_fail
//  (-1 normally...) on failure
//  TODO: make it machine-independant by using SV_DataIO's features of reading/writing binary data
//        (but Wesley's files aren't written machine-independant anyway...)
//====================================================================================================================
int SC_FeatureHandler::saveHTKfeatures(const char *fileName, SV_Data *pFeatures) {
  SC_HTK_HeaderType header;
  FILE *outFile;
  int x;
  float *row;
  int res = 0;
	char tmp[sclib::bufferSize];
  
  if ((outFile=fopen(fileName,"wb")) == NULL) {
		sprintf(tmp, "Save HTK data failed (could not open file '%s')!", fileName);
    REPORT_ERROR(SVLIB_BadData, tmp);
    return SVLIB_Fail;
  }
  //fill header
  header.sampleCount = pFeatures->Row;
  header.featureDim = pFeatures->Col;
  header.sampleStep = 0;
  header.sampleSize = pFeatures->Hdr.frameSize;
  header.sampleType = pFeatures->Hdr.frameStep;

  //write header
  res = (int)(fwrite(&(header.sampleCount), sizeof(int), 1, outFile));
  res += (int)(fwrite(&(header.sampleStep), sizeof(int), 1, outFile));
  res += (int)(fwrite(&(header.sampleSize), sizeof(short), 1, outFile));
  res += (int)(fwrite(&(header.sampleType), sizeof(short), 1, outFile));

  MArray_1D(row, header.featureDim, float, "LoadSignal: row");
  
  //write data
  for (long int y = 0; y < header.sampleCount; y++){
    for (x = 0; x < header.featureDim; x++) {
      row[x] = pFeatures->Mat[y][x];
    }
	  res += (int)(fwrite(row, sizeof(float), header.featureDim, outFile));
  }

  MFree_1D(row);  
  fclose(outFile);

  if (res != ((header.sampleCount * header.featureDim) + 4)) {
    REPORT_ERROR(SVLIB_Fail, "Save HTK data failed (not all bytes written)!");
  }

	return res;
}

//====================================================================================================================
//	Returns the overall number of rows in all the SV_Data objects in the (possibly) linked list of objects
//  counts only the rows in the first segmentsToCount objects or all if this is =0
//====================================================================================================================
unsigned long int SC_FeatureHandler::countFeaturesInList(SV_Data *pFeatureList, unsigned long int segmentsToCount) {
  unsigned long int res = 0, count = 0;
  SV_Data *pHook = pFeatureList;

  while (pHook != NULL && (segmentsToCount == 0 || count < segmentsToCount)) {
    res += pHook->Row;
    pHook = pHook->Next;
    count++;
  }

  return res;
}

//====================================================================================================================
//	Destruct array of (linked lists of) features as returned by extractFeatures() method
//====================================================================================================================
void SC_FeatureHandler::freeFeatures(SV_Data** &pFeatures) {
	if (pFeatures != NULL) {
		for (unsigned int x = 0; x < this->featureCount; x++) {
			if (pFeatures[x] != NULL) {
				sclib::destructLinkedList(pFeatures[x]);
			}
		}
		MFree_1D(pFeatures);
	}

  return;
}

//====================================================================================================================
//	Check if the values of the features are alright (finite)
//  The variable parameters provide information on where th (first) error occured in case of defective features
//====================================================================================================================
bool SC_FeatureHandler::checkFeatures(SV_Data **pFeatures, unsigned long int &feature, unsigned long int &row, unsigned long int &column) {
	if (pFeatures == NULL) {
		return false;
	}

	for (unsigned int i = 0; i < this->featureCount; i++) {
		if (pFeatures[i] != NULL) {
			if (!checkFeatures(pFeatures[i], row, column)) {
				feature = sclib::bitValue(i);
				return false;
			}
		}
	}

  return true;
}

bool SC_FeatureHandler::checkFeatures(SV_Data *pFeatures, unsigned long int &row, unsigned long int &column) {
  for (int y = 0; y < pFeatures->Row; y++) {
    for (int x = 0; x < pFeatures->Col; x++) {
      if (!sclib::isFinite(pFeatures->Mat[y][x])) {
				row = y; 
				column = x;
        return false;
      }
    }
  }

  return true;
}

//====================================================================================================================
//	Convert to another frame-size/-step by weighted averaging all given frames in the range of a new sized frame
//  originalSampleRange is the number of samples used to construct the originalFrameSet (used in order to deduce 
//  correct amount of new frames (there might be one missing in the original set that fits in the new one)). 
//  "ignoreValues" may point to an array of values resembling a valid row of pOriginalFrameSet, which will be excluded
//  from bulding the linear weighted mean (e.g. if there are values of special meaning, like 0.0 as a value describing 
//  that there is no pitch in an pitch-feature-set)
//====================================================================================================================
SV_Data* SC_FeatureHandler::convertFrameRate(SV_Data *pOriginalFrameSet, unsigned long int originalSampleRange, unsigned int newFrameSize, unsigned int newFrameStep, double *ignoreValues) {
	int newFrameCount = sclib::getRowCount(originalSampleRange, newFrameSize, newFrameStep);
	int startSample, endSample, startFrame, endFrame, frameStart, z;
	double weight, *weightSum, *mean;
	SC_Conversion converter(pOriginalFrameSet->Hdr.sampleRate);
	SV_Data *pNewFrameSet = NULL;

	//create new frame-set
	pNewFrameSet = new SV_Data(newFrameCount, pOriginalFrameSet->Col);
  pNewFrameSet->Hdr.frameSize = newFrameSize;
  pNewFrameSet->Hdr.frameStep = newFrameStep;
  pNewFrameSet->Hdr.sampleRate = pOriginalFrameSet->Hdr.sampleRate;
	pNewFrameSet->Hdr.ID = pOriginalFrameSet->Hdr.ID;
	for (z = 0; z < 8; z++) {
		pNewFrameSet->Hdr.Signature[z] = pOriginalFrameSet->Hdr.Signature[z];
	}

	MArray_1D(mean, pOriginalFrameSet->Col, double, "SC_FeatureHandler.convertFrameRate: mean");
	MArray_1D(weightSum, pOriginalFrameSet->Col, double, "SC_FeatureHandler.convertFrameRate: weightSum");
	for (z = 0; z < pOriginalFrameSet->Col; z++) {
		mean[z] = 0.0;
		weightSum[z] = 0.0;
	}

	for (int x = 0; x < pNewFrameSet->Row; x++) {
		//first and last sample inside the current new frame
		startSample = converter.audioFrame2sample(x, newFrameSize, newFrameStep, sclib::alignmentStart);
		endSample = converter.audioFrame2sample(x, newFrameSize, newFrameStep, sclib::alignmentEnd);

		//first and last old frame overlapping with the current new frame
		startFrame = converter.sample2audioFrame(startSample, pOriginalFrameSet->Hdr.frameSize, pOriginalFrameSet->Hdr.frameStep, sclib::alignmentEnd);
		endFrame = sclib::min(pOriginalFrameSet->Row-1, converter.sample2audioFrame(endSample, pOriginalFrameSet->Hdr.frameSize, pOriginalFrameSet->Hdr.frameStep, sclib::alignmentStart));
		for (int y = startFrame; y <= endFrame; y++) {
			frameStart = converter.audioFrame2sample(y, pOriginalFrameSet->Hdr.frameSize, pOriginalFrameSet->Hdr.frameStep, sclib::alignmentStart);
			weight = (double)(sclib::intersect(startSample, endSample, frameStart, frameStart+pOriginalFrameSet->Hdr.frameSize)) / (double)(newFrameSize); //fraction of samples of overlap between both frames related to overall amount of samples in new frame
			for (int z = 0; z < pOriginalFrameSet->Col; z++) {
				if (ignoreValues == NULL || pOriginalFrameSet->Mat[y][z] != ignoreValues[z]) {
					mean[z] += pOriginalFrameSet->Mat[y][z] * weight;
					weightSum[z] += weight; //due to the individual ignoreValues per column, each column can have an individual weight
				}
			}
		}

		//care for proper normalization (the intersection-regions of different old frames with the current new frame may have overlapped so that weightSum > 1.0) and copy the new frame
		for (int z = 0; z < pOriginalFrameSet->Col; z++) {
			if (weightSum[z] > 0.0) { //weightSum[z] == 0 iif all values for this column where ignored; then, not just only use 0.0 as the new feature value (because it's the default value of mean), but the corresponding ignoreValue itself
				pNewFrameSet->Mat[x][z] = (float)(mean[z] / weightSum[z]);
			} else {
				pNewFrameSet->Mat[x][z] = (float)(ignoreValues[z]);
			}
			mean[z] = 0.0;
			weightSum[z] = 0.0;
		}
	}

	MFree_1D(mean);
	MFree_1D(weightSum);

	return pNewFrameSet;
}

//====================================================================================================================
//	Each Algorithm class gets as input the array pFeatures of all extracted features, one array-entry for a feature.
//  It expects that if it needs sclib::featureXYZ (a feature ID) for processing, it will find it at position 
//  pFeatures[sclib::bit(sclib::featureXYZ]. But there is only room for one such feature set, what happens if 
//  different instances of the same feature where extracted using different parameters (each algorithm can specifiy 
//  special parameter sets for the features it uses in case they deviate from the standard püarameters)? they are 
//  deposited as linked lists at the specified location in the array; but, further on, how does the algorithm class 
//  now find "its" feature-set (with it's special parameters) in the linked list? it doesn't have to find it, it just
//  always takes the first element of the linked list. THIS method, called just before the call to the algorithm 
//  itself, cares for putting the right linked-list-elements at the frontal position, given a linked list of the 
//  special parameters an algorithm wants. so, after calling this method, all non-specific feature-sets (those 
//  extracted using the standrad parameters) reside in the front position of the linked lists in the pFeatures-array-
//  entries, except for those feature-sets for which this method got a special parameter-set; in this case, those 
//  feature-sets corresponding to the given parameter-set are pushed to the front (coherence between feature-sets and 
//  parameter-sets is established via the "client"-member of a parameter-set and the first character in the 
//  "Hdr.Signature"-enry of a feature-set, which should both hold a non-zero algorithm ID).
//====================================================================================================================
void SC_FeatureHandler::prepareFeatureSet(SV_Data **pFeatures, SC_TweakableParameters::SC_FeaturePar *pSpecialParameterList) {
	SC_TweakableParameters::SC_FeaturePar *pParameterHook;
	SV_Data *pFeatureHook, *pPreviousFeature;
	unsigned long int featureType;
	unsigned long int algorithm;

	for (unsigned long int i = 0; i < this->featureCount; i++) {
		featureType = sclib::bit(i);

		if (pFeatures[i] != NULL && pFeatures[i]->Next != NULL) { //to save time, do the following only if there is a linked list of features at the current position
			//determine for which client (std: sclib::algorithm_none) the current feature set has to be pushed forward
			algorithm = sclib::algorithm_nothing;
			pParameterHook = pSpecialParameterList;
			while (pParameterHook != NULL) {
				if (pParameterHook->featureType == featureType) {
					algorithm = pParameterHook->client;
					break;
				}
				pParameterHook = pParameterHook->Next;
			}

			//find the feature-set for the determined client, if it is NULL, nothing is to be done
			pPreviousFeature = NULL;
			pFeatureHook = pFeatures[i];
			while (pFeatureHook != NULL) {
				if (pFeatureHook->Hdr.Signature[0] == (char)(algorithm)) {
					break;
				}
				pPreviousFeature = pFeatureHook;
				pFeatureHook = pFeatureHook->Next;
			}

			//bring the found feature-set to the first position in the list (pPreviousFeature contains its direct successor, pFeatureHook the found set itself) 
			if (pFeatureHook != NULL && pFeatureHook != pFeatures[i]) {
				//remove the found feature-set from the list
				if (pPreviousFeature != NULL) {
					pPreviousFeature->Next = pFeatureHook->Next;
				}

				//install the found feature set as the new list-start
				pPreviousFeature = pFeatures[i];
				pFeatures[i] = pFeatureHook;
				pFeatures[i]->Next = pPreviousFeature;
			}
		}
	}

	return;
}

//====================================================================================================================
//  returns true if all wanted feature sets (and those with special parameters) are existent, false if one is NULL
//====================================================================================================================
bool SC_FeatureHandler::checkFeatureSet(SV_Data **pFeatures, unsigned long int featureTypes, SC_TweakableParameters::SC_FeaturePar *pSpecialParameterList) {
	SC_TweakableParameters::SC_FeaturePar *pParameterHook;
	SV_Data *pFeatureHook;
	unsigned long int featureType;
	bool found;
	SC_Conversion converter;

	for (unsigned long int i = 0; i < this->featureCount; i++) {
		featureType = sclib::bit(i);
		if (sclib::bitTest(featureTypes, featureType) == true) { //check all wanted features on existence
			if (pFeatures[i] == NULL) { //wanted feature doesn't exist at all => report failure
				return false;
			} else { //wanted feature with special parameters doesn't exist => report failure
				pParameterHook = pSpecialParameterList;
				while (pParameterHook != NULL) {
					if (pParameterHook->featureType == featureType) { //we want special parameters for this feature
						pFeatureHook = pFeatures[i];
						found = false;
						while (pFeatureHook != NULL) {
							if (pFeatureHook->Hdr.Signature[0] == pParameterHook->client) { //we found the feature set with special parameters
								found = true;
								break;
							}	
							pFeatureHook = pFeatureHook->Next;
						}
						if (found == false) { //last chance: check if special parameters equal standard parameters so that no extra effort was taken to extract it twice
							if (pParameterHook->frameSize == converter.sample2ms(pFeatures[i]->Hdr.frameSize, pFeatures[i]->Hdr.sampleRate) &&
								  pParameterHook->frameStep == converter.sample2ms(pFeatures[i]->Hdr.frameStep, pFeatures[i]->Hdr.sampleRate) &&
									(pParameterHook->sampleRate == pFeatures[i]->Hdr.sampleRate || pParameterHook->sampleRate == 0.0)) {
								found = true;
							}
						}
						if (found == false) { //we didn't find it!
							return false;
						}
					}
					pParameterHook = pParameterHook->Next;
				} //while				
			} //if special parameters are wanted
		} //if it is a wanted feature
	} //for all possible features

	return true;
}

//====================================================================================================================
//	randomizes the order of rows in the given frameset; in fact, not frames but blocks of size #preserveBlockSize 
//  frames are re-ordered, i.e. if e.g. preserveBlockSize==5, the first 5 frames will be consecutive also in the 
//  original, than comes the next block of 5 originally consecutive samples, and so forth
//====================================================================================================================
void SC_FeatureHandler::splice(SV_Data *pFeatures, unsigned int preserveBlockSize) {
	unsigned long int i, j, d, idx;
	double randomValue;
	vector<unsigned long int> inVec;
	SV_Data *pTemp = NULL;

	if (pFeatures != NULL && pFeatures->Row > 0) {
		//create a list of all possible indices into the feature-matrix
		for (i = 0; i < (unsigned long int)(pFeatures->Row); i += preserveBlockSize) {
			inVec.push_back(i);
		}

		//create a mapping from input-indices to output indices and reorder the set in a temp matrix...
		pTemp = new SV_Data(pFeatures->Row, pFeatures->Col);
		i = 0;
		do {
			randomValue = sclib::rand(0.0, (double)(inVec.size()-1)); //randomly select an *index* into the set of still available indices of the original problem
			idx = sclib::round(randomValue);

			for (j = 0; j < sclib::min(preserveBlockSize, pFeatures->Row-inVec[idx]); j++) {
				for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
					pTemp->Mat[i][d] = pFeatures->Mat[inVec[idx]+j][d]; //the inVec[idx]+j'th row will be stored at position i
				}
				i++;
			}
			inVec.erase(inVec.begin()+idx); //remove the just selected block-index from the list of available (still selectable) indices
		} while (inVec.empty() == false);

		//...and copy the content back
		for (i = 0; i < (unsigned long int)(pFeatures->Row); i++) {
			for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
				pFeatures->Mat[i][d] = pTemp->Mat[i][d];
			}
		}
		MFree_0D(pTemp);
	}

	return;
}

//====================================================================================================================
//	creates and returns a new set of features with #intermediateFrameCount frames between two successive frames of the
//  original frame-set, which interpolate between them via a sigmoid function of given steepness; blendStep 
//  corresponds with the preserveBlockSize-parameter of splice(), so that glueing frames could only be inserted all 
//  #blendStep frames.
//====================================================================================================================
SV_Data* SC_FeatureHandler::blend(SV_Data *pFeatures, unsigned int intermediateFrameCount, double steepness, unsigned int blendStep) {
	unsigned long int t, i, d, currentIdx = 0;
	SV_Data *pBlended = NULL;
	double oneAt, pointFiveAt, zeroAt, sig1, sig2;

	if (pFeatures != NULL && pFeatures->Row > 0) {
		pBlended = new SV_Data(((pFeatures->Row-1)/blendStep)*intermediateFrameCount + pFeatures->Row, pFeatures->Col);
		pBlended->Hdr = pFeatures->Hdr;

		for (t = 0; t < (unsigned long int)(pFeatures->Row); t++) {
			//copy the t'th frame
			for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
				pBlended->Mat[currentIdx][d] = pFeatures->Mat[t][d];
			}

			//add blending frames at the correct positions
			if (t+1>1 && (t+1)%blendStep==0 && t+1<(unsigned long int)(pFeatures->Row)) {
				//interpolate on #intermediateFrameCount frames between t'th and (t+1)th frame
				oneAt = (double)(currentIdx);
				pointFiveAt = (double)(currentIdx + (intermediateFrameCount+1.0)/2.0);
				zeroAt = (double)(currentIdx + intermediateFrameCount + 1);
				currentIdx++; //idx where to start first interpolated frame
				for (i = 0; i < intermediateFrameCount; i++) {
					sig1 = sclib::sigmoid((double)(currentIdx), zeroAt, pointFiveAt, oneAt, 0.000001, steepness);
					sig2 = 1.0 - sig1;
					for (d = 0; d < (unsigned long int)(pFeatures->Col); d++) {
						pBlended->Mat[currentIdx][d] = (float)(sig1*pFeatures->Mat[t][d] + sig2*pFeatures->Mat[t+1][d]);
					}
					currentIdx++;
				}
			} else {
				currentIdx++;
			}
		}

		//report error if index computation went awry
		if (currentIdx > (unsigned long int)(pBlended->Row)) {
			REPORT_ERROR(SVLIB_Fail, "Index error while blending features");
		}
	}

	return pBlended;
}

//====================================================================================================================
//	same algorithm as above, but different storage container: featurs is an array with pointer to individual array-
//  rows, one for each row (i.e. its *not* a matrix that can be deallocated by MFree_2D()!). the result is alike.
//====================================================================================================================
double** SC_FeatureHandler::blend(double **features, unsigned int length, unsigned int dim, unsigned int intermediateFrameCount, double steepness, unsigned int blendStep, unsigned int &resultingRows) {
	unsigned long int t, i, d, currentIdx = 0;
	double **blended = NULL, *row = NULL;
	double oneAt, pointFiveAt, zeroAt, sig1, sig2;

	if (features != NULL && length > 0) {
		resultingRows = (length-1)*intermediateFrameCount + length;
		MArray_1D(blended, resultingRows, double*, "SC_FeatureHandler.blend: blended");

		for (t = 0; t < (unsigned long int)(length); t++) {
			//copy the t'th frame
			MArray_1D(row, dim, double, "SC_Featurehandler.blend: row");
			blended[currentIdx] = row;
			for (d = 0; d < (unsigned long int)(dim); d++) {
				blended[currentIdx][d] = features[t][d];
			}

			//add blending frames at the correct positions
			if ((t+1>1 || blendStep==1) && (t+1)%blendStep==0 && t+1<(unsigned long int)(length)) {
				//interpolate on #intermediateFrameCount frames between t'th and (t+1)th frame
				oneAt = (double)(currentIdx);
				pointFiveAt = (double)(currentIdx + (intermediateFrameCount+1.0)/2.0);
				zeroAt = (double)(currentIdx + intermediateFrameCount + 1);
				currentIdx++; //idx where to start first interpolated frame
				for (i = 0; i < intermediateFrameCount; i++) {
					MArray_1D(row, dim, double, "SC_Featurehandler.blend: row");
					blended[currentIdx] = row;
					sig1 = sclib::sigmoid((double)(currentIdx), zeroAt, pointFiveAt, oneAt, 0.000001, steepness);
					sig2 = 1.0 - sig1;
					for (d = 0; d < (unsigned long int)(dim); d++) {
						blended[currentIdx][d] = (float)(sig1*features[t][d] + sig2*features[t+1][d]);
					}
					currentIdx++;
				}
			} else {
				currentIdx++;
			}
		}

		//report error if index computation went awry
		if (currentIdx > (unsigned long int)(resultingRows)) {
			REPORT_ERROR(SVLIB_Fail, "Index error while blending features");
		}
	}

	return blended;
}

//====================================================================================================================
//	takes a single feature set (ignoring further linked sets) and returns a new linked list of shorter feature-sets
//  (#splitCount elements) that have the same time-order as the frames in pFeatureSet and an overlap of 
//  #overlapInFrames frames between two consecutive sets
//====================================================================================================================
SV_Data* SC_FeatureHandler::splitFeatureSet(SV_Data *pFeatureSet, unsigned int splitCount, unsigned int overlapInFrames) {
	unsigned int s, i, j, setSize, setStep, currentSetSize;
	SV_Data *pFirst = NULL, *pHook, *pNew; 
	
	setSize = (unsigned int)(pFeatureSet->Row + overlapInFrames*(splitCount-1) + 0.5) / splitCount; //number of rows in each featureset; resulted from solving the equation in sclib::getRowCount() for windowSize using MuPAD light (with windowStep=windowSize-overlapInFrames)
	setStep = setSize - overlapInFrames;

	for (s = 0; s < splitCount; s++) {
		currentSetSize = (s == splitCount-1) ? pFeatureSet->Row-(s*setStep) : setSize; //make the last set so big (maybe little bigger than the rest in fact) that all the rest of the frames fit inside
		pNew = new SV_Data(currentSetSize, pFeatureSet->Col);

		for (i = 0; i < currentSetSize; i++) { //copy the feature-values
			for (j = 0; j < (unsigned int)(pFeatureSet->Col); j++) {
				pNew->Mat[i][j] = pFeatureSet->Mat[s*setStep + i][j];
			}
		}

		if (s == 0) { //link the list
			pFirst = pNew;
			pHook = pFirst;
		} else {
			pHook->Next = pNew;
			pHook = pHook->Next;
		}
	}

	return pFirst;
}

//====================================================================================================================
//	the feature selection algorithm presented in Kotti, Benetos, Kotropoulos, "Computationally Efficient and Robust 
//  BIC-Based Speaker Segmentation", 2008. Input is a linked list of feature-sets with correspondig boolean labels to
//  assign each list element to one of two classes (the performance measure uses the relationship if inter- and intra-
//  class scatter), and the number of feature columns to select. The algorithm uses a depth-first-search branch-and-
//  bound startegy to find the most effective feature set. Returned is an array of charcters, wehere each bit 
//  corresponds with the respective column of the original feature set (see mfcc-feature-class for more on this).
//  see also: "Classification, Parameter Estimation and State Estimation: An Engineering Approach Using MATLAB", 
//  chapter 6
//====================================================================================================================
unsigned char* SC_FeatureHandler::branchBoundFeatureSelection(SV_Data *pFeatureList, bool *labels, unsigned int resultingColumnCount) {
	unsigned char *selected = NULL;
	SC_MatrixFunctions matFunc;
	unsigned long int l, listCount = sclib::getListCount(pFeatureList), trueCount = 0, falseCount = 0;
	SV_Data *pHook;
	double *trueMean, *falseMean, *grossMean, *tmpMean, **withinScatter, **betweenScatter, *diff;
	int d, x, y, x_2, resultLength;
	char *colKey, *bestColKey;
	std::map<std::string, int> visited;
	double bestPerformance;

	if (listCount > 1) {
		//get class-dependent- and gross mean vectors
		trueMean = matFunc.zeros(pFeatureList->Col);
		falseMean = matFunc.zeros(pFeatureList->Col);
		grossMean = matFunc.zeros(pFeatureList->Col);
		pHook = pFeatureList;
		for (l = 0; l < listCount; l++) { //how many true/false samples do we have overall?
			if (labels[l] == true) {
				trueCount += pHook->Row;
			} else {
				falseCount += pHook->Row;
			}
			pHook = pHook->Next;
		}
		pHook = pFeatureList;
		for (l = 0; l < listCount; l++) {
			tmpMean = matFunc.mean(pHook->Mat, pHook->Row, pHook->Col);
			
			for (d = 0; d < pFeatureList->Col; d++) {
				grossMean[d] += tmpMean[d] * (double)(pHook->Row) / (double)(trueCount+falseCount); //renormalize the means
				if (labels[l] == true) {
					trueMean[d] += tmpMean[d] * (double)(pHook->Row) / (double)(trueCount);
				} else {
					falseMean[d] += tmpMean[d] * (double)(pHook->Row) / (double)(falseCount);
				}
			}

			pHook = pHook->Next;
		}

		//get within-class scatter matrix (ML-estimated covariance matrix)
		withinScatter = matFunc.zeros(pFeatureList->Col, pFeatureList->Col);
		betweenScatter = matFunc.zeros(pFeatureList->Col, pFeatureList->Col);
    MArray_1D(diff, pFeatureList->Col, double, "SC_FeatureHandler.branchBoundFeatureSelection: diff");
		pHook = pFeatureList;
		for (l = 0; l < listCount; l++) { //this is standard covariance computation copied from SC_MatrixFunctions, with the additional functionality of
			                                // (a) handling the linked list of feature vectors as one big matrix
			                                // (b) using the true- or false-mean according to the label of the current list-element
																			// (c) normalization by N (instead of division by N-1) is applied (this is an ML estimate of the covariance, not the unbiased estimat...)
			for (y = 0; y < pHook->Row; y++) {
				for (x = 0; x < pHook->Col; x++) {
					diff[x] = (double)(pHook->Mat[y][x]) - ((labels[l] == true)?trueMean[x]:falseMean[x]);
				}
				for (x = 0; x < pHook->Col; x++) {
					for (x_2 = 0; x_2 < pHook->Col; x_2++) {
						withinScatter[x][x_2] += (diff[x] * diff[x_2]) / (double)(trueCount+falseCount);
					}
				}
			}
			pHook = pHook->Next;
		}

		//get betwee-class scatter
		for (y = 0; y < 2; y++) { //loop over the two class-dependent means
			for (x = 0; x < pFeatureList->Col; x++) {
				diff[x] = ((y==0)?trueMean[x]:falseMean[x]) - grossMean[x];
			}
			for (x = 0; x < pFeatureList->Col; x++) {
				for (x_2 = 0; x_2 < pFeatureList->Col; x_2++) {
					betweenScatter[x][x_2] += (((y==0)?(double)(trueCount):(double)(falseCount)) * (diff[x] * diff[x_2])) / (double)(trueCount+falseCount);
				}
			}
		}
		MFree_1D(diff);
		MFree_1D(trueMean);
		MFree_1D(falseMean);
		MFree_1D(grossMean);

		//this helpful to avoid errors from failed matrix inversions down in branch() if there are zero-values inside the scatter matrix (i.e. uncorrelated features)
		matFunc.replaceZeros(withinScatter, pFeatureList->Col, pFeatureList->Col);
		matFunc.replaceZeros(betweenScatter, pFeatureList->Col, pFeatureList->Col);

		//branch & bound
		MArray_1D(colKey, pFeatureList->Col+1, char, "SC_FeatureHandler.branchBoundFeatureSelection: colKey");
		MArray_1D(bestColKey, pFeatureList->Col+1, char, "SC_FeatureHandler.branchBoundFeatureSelection: bestColKey");
		for (d = 0; d < pFeatureList->Col; d++) {
			colKey[d] = 'x'; //this is regarded by the branch() method as the sign the column d is "active" in the current feature set
		}
		colKey[pFeatureList->Col] = '\0';
		bestPerformance = std::numeric_limits<double>::max()*-1.0;
		branch(withinScatter, betweenScatter, colKey, pFeatureList->Col, 0, pFeatureList->Col-resultingColumnCount+1, visited, bestPerformance, bestColKey, &matFunc);
		
		//store result
		resultLength = sclib::max(1, (int)(ceil(sclib::ld((double)(pFeatureList->Col))))); //so many characters that there is 1 bit for each column in the original feature set
		MArray_1D(selected, resultLength, unsigned char, "SC_FeatureHandler.branchBoundFeatureSelection: selected");
		for (d = 0; d < resultLength; d++) {
			selected[d] = 0;
		}
		for (d = 0; d < pFeatureList->Col; d++) {
			if (bestColKey[d] == 'x') {
				selected[d/8] |= sclib::bit(d%8); //the d'th column corresponds with the d%8'th bit in selected[d/8] => set it!
			}
		}

		MFree_1D(colKey);
		MFree_1D(bestColKey);
		MFree_2D(withinScatter);
		MFree_2D(betweenScatter);
	}

	return selected;
}

//====================================================================================================================
//	this is called by branchBoundFeatureSelection() to do the recursive branching with backtracking; the scatter 
//  matrixes are used to compute the performance measure inside, see also branchBoundFeatureSelection(). the colKey
//  is an array with a cell per column of the original feature set, having a 'x' character at those locations where
//  the corresponding column belongs to the current sub-featureset. the result is the best performing colKey.
//====================================================================================================================
void SC_FeatureHandler::branch(double **withinScatter, double **betweenScatter, char *colKey, unsigned int originalDim, unsigned int level, unsigned int maxLevel, std::map<std::string, int> &visited, double &bestPerformance, char* &bestColKey, SC_MatrixFunctions *pMatFunc) {
	double performance, **subWithin, **subBetween;
	char *newColKey;
	unsigned int activeCounter;
	unsigned long int rows, cols;
	
	if (level < maxLevel) {
		if (visited[colKey] == 0) { //only process branches that have not been visited before (same column constellations will arise in several branches and are prohibited by this line from being processed more than once)
			visited[colKey]++; //now we visited this column combination
			
			//evaluate performance measure for current feature subset: tr(Sw^-1 * Sb)
			subWithin = pMatFunc->getSubMatrix(withinScatter, originalDim, originalDim, colKey, colKey, cols, rows, 'x'); //get sub-matrixes of currently selected feature-set
			subBetween = pMatFunc->getSubMatrix(betweenScatter, originalDim, originalDim, colKey, colKey, cols, rows, 'x');
			if (pMatFunc->inv(subWithin, cols) < 0) { //the ^-1 part
				REPORT_ERROR(SVLIB_BadData, "within class scatter matrix during branch&bound feature selection not invertible");
			}
			pMatFunc->mult(subWithin, subBetween, rows, cols, rows, cols, true); //the * part
			performance = pMatFunc->trace(subWithin, cols); //the tr() part
			MFree_2D(subWithin);
			MFree_2D(subBetween);
			
			if (performance > bestPerformance) { //only proceed if this branch is better then the best from the desired deepest level
				if (level == maxLevel-1){ //if this run represents a feature set of the dimensionality we want, we have found the best seen so far
					bestPerformance = performance; //remember performance value of this run if it is on the deepest level
					sprintf(bestColKey, "%s", colKey);
				} else {
					//do recursive branching by selecting a feature-subset with one col less, subsequently in the loop leaving out each col of the currently active set
					MArray_1D(newColKey, originalDim+1, char, "SC_FeatureHandler.branch: newColKey");
					for (unsigned int i = 0; i < originalDim-level; i++) { //for all currently active cols
						sprintf(newColKey, "%s", colKey); //prepare the new colKey
						activeCounter = 0;
						for (unsigned int c = 0; c < originalDim; c++) { //remove the i'th "active-marker" from the newColKey and branch for this new subset (all currently active cols without the i'th active col)
							if (newColKey[c] == 'x') { //(original) column nr. c is marked active in the current feature set
								if (activeCounter++ == i) { //we want to deactivate the i'th active column
									newColKey[c] = '_'; //we do deactivation by marking this column with an underscore in stead af an 'x'
									break;
								}
							}
						} //for c
						branch(withinScatter, betweenScatter, newColKey, originalDim, level+1, maxLevel, visited, bestPerformance, bestColKey, pMatFunc);
					} //for i
					MFree_1D(newColKey);
				} //not the last level
			} //if performance best in this level
		} //if not already visited
	} //if maxlevel not reached

	return;
}

//====================================================================================================================
//	return a subset of columns of the given feature set
//====================================================================================================================
SV_Data* SC_FeatureHandler::getSubSet(SV_Data *pFeatures, unsigned int startCol, unsigned int endCol) {
	SV_Data *pSubSet = NULL;

	if (pFeatures != NULL && endCol >= startCol && endCol < (unsigned int)(pFeatures->Col)) {
		pSubSet = new SV_Data(pFeatures->Row, endCol-startCol+1);
		pSubSet->Hdr = pFeatures->Hdr;
		for (int y = 0; y < pFeatures->Row; y++) {
			for (unsigned int x = startCol; x <= endCol; x++) {
				pSubSet->Mat[y][x-startCol] = pFeatures->Mat[y][x];
			}
		}
	}

	return pSubSet;
}

//====================================================================================================================
//	replaces outliers in the data by the median of the resprective dimension; regard 1.5*IQR (inter quartile range) as
//  mild outliers, 3.0*IQR as extreme outliers. http://in.answers.yahoo.com/question/index?qid=1006052718450
//  the percentage of outliers is returned
//====================================================================================================================
double SC_FeatureHandler::removeOutliers(SV_Data *pFeatures, const char mode) {
	double IQR, lowerBound, upperBound, factor = 3.0; //the standrad for sclib::outlierRemoveExtreme
	float *percentiles, *col;
	unsigned long int count = 0;

	if (mode != sclib::outlierRemoveNone) {
		SC_MatrixFunctions matFunc;
		if (mode == sclib::outlierRemoveMild) {
			factor = 1.5;
		}

		for (int d = 0; d < pFeatures->Col; d++) {
			col = matFunc.getCol(pFeatures->Mat, pFeatures->Row, pFeatures->Col, d);
			percentiles = sclib::percentiles(col, pFeatures->Row, true);

			IQR = percentiles[75] - percentiles[25]; //the inter quartile range
			lowerBound = percentiles[25] - factor*IQR;
			upperBound = percentiles[75] + factor*IQR;

			for (int t = 0; t < pFeatures->Row; t++) {
				if (pFeatures->Mat[t][d] < lowerBound || pFeatures->Mat[t][d] > upperBound) {
					pFeatures->Mat[t][d] = percentiles[50]; //replace outliers with the median value
					count++;
				}
			}

			MFree_1D(col);
			MFree_1D(percentiles);
		}
	}

	return (double)(count) / (double)(pFeatures->Row*pFeatures->Col);
}

//====================================================================================================================
//	set all values below minPitch to zero (important to clean pitch sampled from a GMM where zero isn't zero)
//====================================================================================================================
void SC_FeatureHandler::cleanPitch(SV_Data *pPitch, int column, double minPitch) {
	for (int y = 0; y < pPitch->Row; y++) {
		if (pPitch->Mat[y][column] < minPitch) {
			pPitch->Mat[y][column] = 0.0f;
		}
	}

	return;
}

//====================================================================================================================
//	concatenates #framesPerTrajectory frames to a feature trajectory in order to incorporate time-dependencies into a
//  feature vector; move #trajectoryStep frames forward to begin the next trajectory.
//  if removeTiming==true, the original data is clustered via #clusteringIterations of kMeans (k=templateCount; if 0, 
//  66% of #frames) and each vector is the replaced by the nearest template; then, successive identical templates are 
//  reduced to just one instance before the trajectories are formed. this way it is hoped to remove timing differencs 
//  in trajectories (e.g. faster or slower articulation of syllables) by retaining the original spectral content
//====================================================================================================================
SV_Data* SC_FeatureHandler::createTrajectories(SV_Data *pFeature, unsigned int framesPerTrajectory, unsigned int trajectoryStep, bool removeTiming, unsigned int templateCount, unsigned int clusteringIterations) {
	int x, y = 0, t, d, lastX = -1, finalRowCount = pFeature->Row;
	SV_Data *pTrajectories, *pTemplates, *pReplaced = pFeature;
	SC_Clusterer clusterer(this->pTweak, this->verbose);
	double shrinkFactor = 0.66; //sounds well in experiments: speech is speeded up but still well intellegible

	//remove timing information if wished
	if (removeTiming == true) {
		if (templateCount == 0) {
			templateCount = sclib::round(pFeature->Row * shrinkFactor);
		}
		pTemplates = clusterer.kMeans(templateCount, pFeature, clusteringIterations, false, 5, false);
		pReplaced = new SV_Data(pFeature->Row, pFeature->Col);
		pReplaced->Hdr = pFeature->Hdr;

		//pick the nearest template for each vector; if it is different from its predecesor's template, insert the template in the pReplaced set
		//ohterwise the vector at all; this way, pReplaced has as many rows as pFeature, but only #finalRowCount rows are filled afterwards
		for (t = 0; t < pFeature->Row; t++) {
			x = (int)(clusterer.pickNearest(pFeature->Mat[t], pTemplates));

			if (x != lastX) {
				for (d = 0; d < pFeature->Col; d++) {
					pReplaced->Mat[y][d] = pTemplates->Mat[x][d];
				}
				y++;
			}

			lastX = x;
		}
		finalRowCount = y;

		//SC_Synthesis synthesizer(this->pTweak);
		//char buffer[sclib::bufferSize];
		//sprintf(buffer, "../data/test/replaced_%f.wav", shrinkFactor);
		//pReplaced->Row = finalRowCount; //TODO
		//synthesizer.feature2wav(buffer, pReplaced);

		MFree_0D(pTemplates);
	}

	//create trajectory container; take pFeature as the template to copy header info from because it wasn't filled in pReplaced
	pTrajectories = new SV_Data(sclib::getRowCount(finalRowCount, framesPerTrajectory, trajectoryStep), framesPerTrajectory*pFeature->Col);
	pTrajectories->Hdr = pFeature->Hdr;
	pTrajectories->Hdr.Signature[2]++; //remember that this is trajectory-set by incrementing the 3rd signature-entry (is 0 elsewise); copy the rest
	pTrajectories->Hdr.Signature[3] = (char)(framesPerTrajectory); //remember those parameters to be able to unwind trajectories later on
	pTrajectories->Hdr.Signature[4] = (char)(trajectoryStep); 
	pTrajectories->Hdr.frameSize = pFeature->Hdr.frameSize + (framesPerTrajectory-1)*pFeature->Hdr.frameStep;
	pTrajectories->Hdr.frameStep = trajectoryStep * pFeature->Hdr.frameStep;

	//create the actual trajectories
	for (y = 0; y < pTrajectories->Row; y++) {
		x = 0;
		for (t = 0; t < (int)(framesPerTrajectory); t++) {
			for (d = 0; d < pReplaced->Col; d++) {
				pTrajectories->Mat[y][x++] = pReplaced->Mat[y*trajectoryStep + t][d];
			}
		}
	}

	//clean up
	if (pReplaced != pFeature) {
		MFree_0D(pReplaced);
	}

	return pTrajectories;
}

//====================================================================================================================
//	convert trajectories back to single frames; just what was removed when removing timing misses
//====================================================================================================================
SV_Data* SC_FeatureHandler::unwindTrajectories(SV_Data *pTrajectories, unsigned int framesPerTrajectory, unsigned int trajectoryStep) {
	SV_Data *pFrames = NULL;
	int rowCount, colCount, y = 0, x, d, t;

	if (pTrajectories != NULL) {
		if (pTrajectories->Hdr.Signature[2] == 0) {
			pFrames = new SV_Data(*pTrajectories);
		} else {
			colCount = pTrajectories->Col / framesPerTrajectory;
			rowCount = framesPerTrajectory + (pTrajectories->Row-1)*trajectoryStep;
			pFrames = new SV_Data(rowCount, colCount);
			pFrames->Hdr = pTrajectories->Hdr;
			pFrames->Hdr.Signature[2] = 0;
			pFrames->Hdr.Signature[3] = 0;
			pFrames->Hdr.Signature[4] = 0;
			pFrames->Hdr.frameStep = pTrajectories->Hdr.frameStep / trajectoryStep;
			pFrames->Hdr.frameSize = pTrajectories->Hdr.frameSize - (framesPerTrajectory-1)*pFrames->Hdr.frameStep;

			for (t = 0; t < pTrajectories->Row-1; t++) {
				x = 0;
				for (d = 0; d < (int)(pFrames->Col*trajectoryStep); d++) {
					if (d>0 && d%pFrames->Col==0) {
						x = 0;
						y++;
					}
					pFrames->Mat[y][x++] = pTrajectories->Mat[t][d];
				}
				y++;
			}
			x = 0;
			for (d = 0; d < pTrajectories->Col; d++) { //last trajectory
				if (d>0 && d%pFrames->Col==0) {
					x = 0;
					y++;
				}
				pFrames->Mat[y][x++] = pTrajectories->Mat[pTrajectories->Row-1][d];
			}
		}
	}

	return pFrames;
}

//====================================================================================================================
//	load data stored in libSVM format
//====================================================================================================================
SV_Data* SC_FeatureHandler::loadLibSVMfeatures(const char *fileName, int* &labels) {
	SV_Data *pFeatures;  
	unsigned long int lineCount, t = 0, dimCount = 0;
  FILE *file = NULL; 
	char buffer[sclib::bufferSize], *bufferPtr;
	int startPos, d, bufferLength;
	double value;
	
	if (sclib::fileExists(fileName) == true) {
		lineCount = sclib::countLines(fileName);
		file = fopen(fileName, "r");
			//count nr. of colons, that is: number of dimensions; assume that it is the full dumber of dimensions in this first vector
			bufferLength = sclib::readline(file, buffer, sclib::bufferSize);
			bufferPtr = buffer;
			while (bufferPtr[0] != '\0') { 
				if (*(bufferPtr++) == ':') {
					dimCount++;
				}
			}

			//create destination objects
			pFeatures = new SV_Data(lineCount-1, dimCount);
			MFree_1D(labels);
			MArray_1D(labels, lineCount, int, "SC_FeatureHandler.loadLibSVMdata: labels");

			//parse file and fill destinations
			while (buffer[0]=='1'  || buffer[0]=='0') {
				startPos = sclib::getNextIntFromString(buffer, bufferLength, labels[t], 0, " :\t\n\0\r");

				do {
					startPos = sclib::getNextIntFromString(buffer, bufferLength, d, startPos, " :\t\n\0\r");
					startPos = sclib::getNextDoubleFromString(buffer, bufferLength, value, startPos, " :\t\n\0\r");
					if (d > 0) {
						pFeatures->Mat[t][d-1] = (float)(value);
					} else {
						break;
					}
				} while (startPos < bufferLength);
				
				t++;
				bufferLength = sclib::readline(file, buffer, sclib::bufferSize); //get next line
			}
    fclose(file);
  }
	
  return pFeatures;
}

//====================================================================================================================
//	returns the MD5 checksum of the data matrix (not regarding the header or linked list nature or stuff)
//====================================================================================================================
char* SC_FeatureHandler::getChecksum(SV_Data *pFeature, unsigned long int segmentsToMerge) {
	SC_MD5 md5er;
	int error = 0, count = 0;
	char *checksum = NULL;
	SV_Data *pHook = pFeature;

	while ((error>=0) && (pHook!=NULL) && ((count<(int)(segmentsToMerge)) || (segmentsToMerge==0))) {
		md5er.Update((unsigned char*)(&pHook->Mat[0][0]), sizeof(float)*pHook->Row*pHook->Col, error);
		count++;
	}

	if (error >= 0) {
		checksum = md5er.Final(error);
		if (error < 0) {
			MFree_1D(checksum);
		}
	}

	return checksum;
}

//====================================================================================================================
//	shrink the given dataset by percent% via setting the row-count to a smaller value (the original data remains the 
//  same; it can be restored by setting pFeature->Row to its original value)
//====================================================================================================================
int SC_FeatureHandler::reduceByPercent(SV_Data *pFeature, double percent) {
	int originalRowCount = pFeature->Row;
	int newRowCount = sclib::round((double)(pFeature->Row) * (1.0-(percent/100.0)));

	pFeature->Row = newRowCount;	

	return originalRowCount;
}
