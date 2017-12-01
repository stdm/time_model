/**************************************************************************/
/*    Copied from:																												*/
/*      - This is a altered version of SV_Feature_MFCC	by Jialong He			*/
/*																																				*/
/*    Responsibility:																											*/
/*      - Class for extracting Filterbank Energys 										    */
/*																																				*/
/*    Background:                                                         */
/*      - The different frequency scales are due to the definition in 		*/
/*				"A comparative study of traditional and newly proposed          */
/*				 features for recognition of speech under stress", Bou-Ghazale, */
/*				Hansen, IEEE Trans. Speech & Audio Proc., July 2000							*/
/*																																				*/
/*    Author  : Jialong He / Thilo Stadelmann															*/
/*    Date    : 18.10.2004																								*/
/**************************************************************************/

#include "SC_Aux.h"
#include "SC_Feature_FbE.h"
#include "SC_Transform.h"
#include "SC_Signal_WAVE.h"
#include "SC_Conversion.h"
#include "SC_FeatureHandler.h"
#include <GN_Filter.h>

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_FbE::SC_Feature_FbE(int sampleRate, int frameSize, int frameStep, int filterBankSize, int FFTsize, int window, double preEmphasizeFactor, double lowestFrequency, double highestFrequency, int MFCCorder, bool dEnergy, unsigned char scale, unsigned char result, unsigned int smoothing, double taperingLength) : SV_Feature() {
	this->Para.WinSz = frameSize;
	this->Para.StpSz = frameStep;
	this->Para.SRate = sampleRate;
  this->Para.NFilter = filterBankSize;
	this->Para.FFTSz = FFTsize;
	this->Para.HammingWin = window;
	this->Para.RmvSilence = 0;
	this->Para.Alpha = preEmphasizeFactor;
  this->result = result;
  this->scale = scale;
  this->lowFrequency = lowestFrequency;
	this->highFrequency = (highestFrequency <= 0.0) ? this->Para.SRate/2 : highestFrequency;
  this->Para.MFCC_Order = MFCCorder;
  this->Para.DEnergy = (dEnergy == true) ? 1 : 0;
	this->taperingLength = taperingLength;

  //determine the smoothing-strength
  switch (smoothing) {
    case sclib::smoothLight: {
      this->smoothSpectrum = (unsigned int)floor(this->Para.FFTSz - this->Para.FFTSz/10.0); //cuts highest 10 percent of the frequencys of the powerSpectrum
      break;
    }
    case sclib::smoothMiddle: {
      this->smoothSpectrum = (unsigned int)floor(this->Para.FFTSz / 3.0); //holds lowest one third of the frequencys of the powerspectrum
      break;
    }
    case sclib::smoothHeavy: {
      this->smoothSpectrum = 6; //holds only the lowest 6 frequencys of the powerSpectrum
      break;
    }
    case sclib::smoothNone: {
      this->smoothSpectrum = 0;
      break;
    }
    default: {
      this->smoothSpectrum = 0;
      break;
    }
  }
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_FbE::~SC_Feature_FbE() {
	
}

//====================================================================================================================
// returns an array of the filter-parameters for filterbank spaced according to "mode" and the lowest and highest 
// frequency to cover (the upper and lower end are the center-frequencies of the neighboring filters):
// the first component holds the lower end frequency (Hz)
//     second                    center-frequency (Hz)
//     third                     upper end frequency (Hz)
//====================================================================================================================
double** SC_Feature_FbE::createSpacing(unsigned char mode, double fMinHz, double fMaxHz, int filterCount) {
  double **filterParam;
  double (*Hz2mode)(double f) = getFrequencyMappingHz2X(mode);
  double (*mode2Hz)(double f) = getFrequencyMappingX2Hz(mode);
  double fMin = (*Hz2mode)(fMinHz);
  double fMax = (*Hz2mode)(fMaxHz);
  double width = (fMax - fMin) / (filterCount + 1); //actually, half the width

  MArray_2D(filterParam, (long)(filterCount), 3, double, "SC_Feature_FbE.createLinearSpacing: filterParam");

  for (int x = 0; x < filterCount; x++) {
    filterParam[x][0] = (x > 0) ? filterParam[x-1][1] : fMinHz; //lower end
    filterParam[x][1] = (*mode2Hz)(fMin + (width * (x+1))); //center
    if (x > 0) { //upper end
      if (x < filterCount-1) {
        filterParam[x-1][2] = filterParam[x][1];
      } else {
        filterParam[x-1][2] = filterParam[x][1];
        filterParam[x][2] = fMaxHz;
      }
    }
  }

  return filterParam;
}

//====================================================================================================================
// returns a pointer to the frequency-convert-function Hz->mode
//====================================================================================================================
double (*SC_Feature_FbE::getFrequencyMappingHz2X(unsigned char mode))(double) {
  double (*fConvert)(double f);

  switch (mode) {
    case sclib::scaleMel: {
      fConvert = &Hz2mel; 
      break;
    }
    case sclib::scaleExpoLog: {
      fConvert = &Hz2expoLog;
      break;
    }
    case sclib::scaleModMel: {
      fConvert = &Hz2mMel;
      break;
    }
		case sclib::scaleBark: {
			fConvert = &Hz2bark;
			break;
		}
    default: {
      fConvert = &Hz2Hz;
      break;
    }
  }

  return fConvert;
}

//====================================================================================================================
// returns a pointer to the frequency-convert-function mode->Hz
//====================================================================================================================
double (*SC_Feature_FbE::getFrequencyMappingX2Hz(unsigned char mode))(double) {
  double (*fConvert)(double f);

  switch (mode) {
    case sclib::scaleMel: {
      fConvert = &mel2Hz; 
      break;
    }
    case sclib::scaleExpoLog: {
      fConvert = &expoLog2Hz;
      break;
    }
    case sclib::scaleModMel: {
      fConvert = &mMel2Hz;
      break;
    }
		case sclib::scaleBark: {
			fConvert = &bark2Hz;
			break;
		}
    default: {
      fConvert = &Hz2Hz;
      break;
    }
  }

  return fConvert;
}

//====================================================================================================================
// converts linear frequency (Hz) to linear frequeny (Hz) (the identity function :-)
//====================================================================================================================
double SC_Feature_FbE::Hz2Hz(double f) {
  return f;
}

//====================================================================================================================
// converts linear frequency (Hz) to mel
//====================================================================================================================
double SC_Feature_FbE::Hz2mel(double f) {
  return 2595.0 * log10(1.0 + (f / 700.0));
}

//====================================================================================================================
// converts mel to linear frequency (Hz)
//====================================================================================================================
double SC_Feature_FbE::mel2Hz(double mel) {
  return 700.0 * (pow((double)(10.0), (double)((mel / 2595.0))) - 1.0);
}

//====================================================================================================================
// converts linear frequency (Hz) to  modified mel
//====================================================================================================================
double SC_Feature_FbE::Hz2mMel(double f) {
  return 3070.0 * log10(1.0 + (f / 1000.0));
}

//====================================================================================================================
// converts modified mel to linear frequency (Hz)
//====================================================================================================================
double SC_Feature_FbE::mMel2Hz(double mMel) {
  return 1000.0 * (pow((double)(10.0), (double)((mMel / 3070.0))) - 1.0);
}

//====================================================================================================================
// converts linear frequency (Hz) to  expoLog
//====================================================================================================================
double SC_Feature_FbE::Hz2expoLog(double f) {
  if (f <= 2000.0) {
    return 700.0 * (pow((double)(10.0), (double)((f / 3988.0))) - 1.0);
  } else {
    return 2595.0 * log10(1.0 + (f / 700.0));
  }
}

//====================================================================================================================
// converts expoLog to linear frequency (Hz)
//====================================================================================================================
double SC_Feature_FbE::expoLog2Hz(double expoLog) {
  if (expoLog <= 1521.27615) {
    return 3988.0 * log10(1.0 + (expoLog / 700.0));
  } else {  
    return 700.0 * (pow((double)(10.0), (double)((expoLog / 2595.0))) - 1.0);
  }
}

//====================================================================================================================
// converts linear frequency (Hz) to bark according to H. Traunmüller (1990) "Analytical expressions for the tonotopic 
// sensory scale" J. Acoust. Soc. Am. 88: 97-100. 
//====================================================================================================================
double SC_Feature_FbE::Hz2bark(double f) {
	double z = -0.53;
	
	if (f > 0.0) {
		z = (26.81 / (1.0 + 1960 / f )) - 0.53; //also known as the critical band rate
	}

	if (z < 2.0) {
		z += 0.15 * (2.0 - z);
	} else if (z > 20.1) {
		z += 0.22 * (z - 20.1);
	}

	return z;
}

//====================================================================================================================
// converts bark to linear frequency (Hz)
// see e.g. http://www.ling.su.se/staff/hartmut/bark.htm
//====================================================================================================================
double SC_Feature_FbE::bark2Hz(double z) {
	double f;

	if (z < 2.0) {
		z = (z - 0.3) / 0.85;
	} else if (z > 20.1) {
		z = (z + 4.422) / 1.22;
	}

	if (z + 0.53 == 0.0) {
		f = 0.0;
	} else {
		f = 1960.0 / (26.81 / (z+0.53) - 1.0);
	}

	return f;
}

//====================================================================================================================
// returns, for a frequency in [Hz], its critical bandwidth in [Hz] according to the Bark scale; if usePublishedScale
// is true, not the bandwidth of exact the given frequency f is returned, but the bandwidth of the bark band this 
// frequency falls into as originally published is returned; additionally, the centerFrequancy xorresponding with this 
// band is returned.
//====================================================================================================================
double SC_Feature_FbE::getCriticalBandwidth(double f, double &centerFrequency, double &leftEdge, double &rightEdge, bool usePublishedScale) {
	double z, bw;

	if (usePublishedScale == false ) {
		z = Hz2bark(f);
		bw = 52548.0 / (z*z - 52.56*z + 690.39);
		centerFrequency = f;
		leftEdge = f - bw/2.0;
		rightEdge = f + bw/2.0;
	} else {
		if (f < 100.0) { //scale according to http://ccrma.stanford.edu/~jos/bbt/Bark_Frequency_Scale.html
			bw = 100.0;
			centerFrequency = 50.0;
			leftEdge = 0.0;
			rightEdge = 100.0;
		} else if (f < 200.0) {
			bw = 100.0;
			centerFrequency = 150.0;
			leftEdge = 100.0;
			rightEdge = 200.0;
		} else if (f < 300.0) {
			bw = 100.0;
			centerFrequency = 250.0;
			leftEdge = 200.0;
			rightEdge = 300.0;
		} else if (f < 400.0) {
			bw = 100.0;
			centerFrequency = 350.0;
			leftEdge = 300.0;
			rightEdge = 400.0;
		} else if (f < 510.0) {
			bw = 110.0;
			centerFrequency = 450.0;
			leftEdge = 400.0;
			rightEdge = 510.0;
		} else if (f < 630.0) {
			bw = 120.0;
			centerFrequency = 570.0;
			leftEdge = 510.0;
			rightEdge = 630.0;
		} else if (f < 770.0) {
			bw = 140.0;
			centerFrequency = 700.0;
			leftEdge = 630.0;
			rightEdge = 770.0;
		} else if (f < 920.0) {
			bw = 150.0;
			centerFrequency = 840.0;
			leftEdge = 770.0;
			rightEdge = 920.0;
		} else if (f < 1080.0) {
			bw = 160.0;
			centerFrequency = 1000.0;
			leftEdge = 920.0;
			rightEdge = 1080.0;
		} else if (f < 1270.0) {
			bw = 190.0;
			centerFrequency = 1170.0;
			leftEdge = 1080.0;
			rightEdge = 1270.0;
		} else if (f < 1480.0) {
			bw = 210.0;
			centerFrequency = 1370.0;
			leftEdge = 1270.0;
			rightEdge = 1480.0;
		} else if (f < 1720.0) {
			bw = 240.0;
			centerFrequency = 1600.0;
			leftEdge = 1480.0;
			rightEdge = 1720.0;
		} else if (f < 2000.0) {
			bw = 280.0;
			centerFrequency = 1850.0;
			leftEdge = 1720.0;
			rightEdge = 2000.0;
		} else if (f < 2320.0) {
			bw = 320.0;
			centerFrequency = 2150.0;
			leftEdge = 2000.0;
			rightEdge = 2320.0;
		} else if (f < 2700.0) {
			bw = 380.0;
			centerFrequency = 2500.0;
			leftEdge = 2320.0;
			rightEdge = 2700.0;
		} else if (f < 3150.0) {
			bw = 450.0;
			centerFrequency = 2900.0;
			leftEdge = 2700.0;
			rightEdge = 3150.0;
		} else if (f < 3700.0) {
			bw = 550.0;
			centerFrequency = 3400.0;
			leftEdge = 3150.0;
			rightEdge = 3700.0;
		} else if (f < 4400.0) {
			bw = 700.0;
			centerFrequency = 4000.0;
			leftEdge = 3700.0;
			rightEdge = 4400.0;
		} else if (f < 5300.0) {
			bw = 900.0;
			centerFrequency = 4800.0;
			leftEdge = 4400.0;
			rightEdge = 5300.0;
		} else if (f < 6400.0) {
			bw = 1100.0;
			centerFrequency = 5800.0;
			leftEdge = 5300.0;
			rightEdge = 6400.0;
		} else if (f < 7700.0) {
			bw = 1300.0;
			centerFrequency = 7000.0;
			leftEdge = 6400.0;
			rightEdge = 7700.0;
		} else if (f < 9500.0) {
			bw = 1800.0;
			centerFrequency = 8500.0;
			leftEdge = 7700.0;
			rightEdge = 9500.0;
		} else if (f < 12000.0) {
			bw = 2500.0;
			centerFrequency = 10500.0;
			leftEdge = 9500.0;
			rightEdge = 12000.0;
		} else if (f < 15500.0) {
			bw = 3500.0;
			centerFrequency = 13500.0;
			leftEdge = 12000.0;
			rightEdge = 15500.0;
		} else if (f < 20500.0) {
			bw = 5000.0;
			centerFrequency = 17500.0; //by thilo
			leftEdge = 15500.0;
			rightEdge = 20500.0;
		} else if (f < 27000.0) {
			bw = 100.0;
			centerFrequency = 23500.0; //by thilo
			leftEdge = 20500.0;
			rightEdge = 27000.0;
		} else {
			bw = 100.0;
			centerFrequency = 31500.0; //by thilo
			leftEdge = 27000.0;
			rightEdge = 50000.0; //by thilo
		}
	}

	return bw;
}

//====================================================================================================================
// returns, for a frequency in [Hz], its ERB in [Hz] according to B.C.J. Moore and B.R. Glasberg (1983) "Suggested 
// formulae for calculating auditory-filter bandwidths and excitation patterns" J. Acoust. Soc. Am. 74: 750-753. 
// the formula is valid for frequencies in the range 0.1-6.5kHz
//====================================================================================================================
double SC_Feature_FbE::getEquivalentRectangularBandwidth(double f) {
	double ERB = 0.0;

	if (f >= 100.0 && f <= 6500.0) {
		ERB = 11.17 * sclib::ln((f+312.0) / (f + 14675.0)) + 43.0;
	}

	return ERB;
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data *SC_Feature_FbE::ExtractFeature(void) {
	SV_Data *pData = NULL;
  SC_Transform *pTrans = new SC_Transform(this->Para.FFTSz, this->Para.HammingWin, this->taperingLength);
  double **filterParam = createSpacing(this->scale, this->lowFrequency, this->highFrequency, this->Para.NFilter);
 	float *sigBuf = NULL;
  double *frame = NULL;
  double *fbe = NULL;
  double *cepstrum = NULL;
  double *spectrum  = NULL;
  unsigned long int  sigLen, frameCnt, filterCnt, sampleCnt;
	int frameCount;
	SC_MatrixFunctions matFunc;

	if (IsSigLoaded()) { 
		sigBuf = GetSig();
		sigLen = GetLen();
	}	else {
		return (NULL);
	}

	frameCount = (int)(sclib::getRowCount(sigLen, Para.WinSz, Para.StpSz)); //(sigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1; //(sigLen - this->Para.WinSz + this->Para.StpSz) / this->Para.StpSz;
	if (frameCount > 0) {
		// Preeamphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}

		pData->Row = frameCount;
		pData->Col = (this->result == sclib::resultCepstrum) ? this->Para.MFCC_Order : this->Para.NFilter;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = (this->result == sclib::resultCepstrum) ? sclib::featureMFCC : sclib::featureFbE;

		MArray_1D(frame, this->Para.WinSz, double, "SC_Feature_FbE.ExtractFeature: frame");
		MArray_1D(fbe, this->Para.NFilter, double, "SC_Feature_FbE.ExtractFeature: features");

		for (frameCnt = 0; frameCnt < (unsigned long)(pData->Row); frameCnt++) {
			for (sampleCnt = 0; sampleCnt < (unsigned long)(this->Para.WinSz); sampleCnt++) {
				frame[sampleCnt] = sigBuf[frameCnt * this->Para.StpSz + sampleCnt];
			}

			//sclib::vectorOut("frame.txt", frame, this->Para.WinSz);

			//compute spectrum
			if (this->smoothSpectrum == 0) { //the spectrum shouldn't be smoothed
				spectrum = pTrans->magnitudeSpectrum(frame, this->Para.WinSz); //according to Slaney, Davis&Mermelstein, HTK and svlib we also use the magnitude-spectrum here... does it do harm => no, seems ok!
			} else { //if this->smoothSpectrum > 0, it's value determines the length of the smoothing window (the more the length, the less is the smoothing-effect!)
				spectrum = pTrans->smoothedMagnitudeSpectrum(frame, this->Para.WinSz, this->smoothSpectrum);
			}

			//sclib::vectorOut("spec.txt", spectrum, this->Para.WinSz/2+1);

			//create filterbank-energys
			for (filterCnt = 0; filterCnt < (unsigned int)(this->Para.NFilter); filterCnt++) {
				fbe[filterCnt] = evalFilter(spectrum, filterParam[filterCnt][0], filterParam[filterCnt][1], filterParam[filterCnt][2], this->Para.SRate, this->Para.FFTSz);
				assert(fbe[filterCnt] >= 0.0);
				if (this->result != sclib::resultLinear) { //TODO: tutorial on text-independent sv suggests 20*log() here for db-scale...
					fbe[filterCnt] = (fbe[filterCnt] == 0.0) ? -87.3 : log(fbe[filterCnt]); //TODO: Slaney uses log10 here... //-87.336544750553102 is the log of the smallest non-neg. float number
				}
			}

			//sclib::vectorOut("logFbe.txt", fbe, this->Para.NFilter);

			//decide whether to compute the cepstral coefficients instead of simple filterbank energy
			assert(this->Para.NFilter >= this->Para.MFCC_Order);
			if (this->result == sclib::resultCepstrum) {
				cepstrum = pTrans->dct(fbe, this->Para.NFilter, this->Para.MFCC_Order);
				for (filterCnt = 0; filterCnt < (unsigned int)(this->Para.MFCC_Order); filterCnt++) {
					pData->Mat[frameCnt][filterCnt] = (float)(cepstrum[filterCnt]);
				}

				/*
				sclib::vectorOut("mfcc.txt", cepstrum, this->Para.MFCC_Order);
				double* backCeps = pTrans->idct(cepstrum, this->Para.MFCC_Order);
				sclib::vectorOut("backCeps.txt", backCeps, this->Para.MFCC_Order);
				*/
	      
				MFree_1D(cepstrum);
			} else {
				for (filterCnt = 0; filterCnt < (unsigned int)(this->Para.NFilter); filterCnt++) {
					pData->Mat[frameCnt][filterCnt] = (float)(fbe[filterCnt]);
				}
			}

			MFree_1D(spectrum);
		}

		// for MFCC's, first MFCC is engery, normally absolute energy is no use but delta energy useful
		if (this->result == sclib::resultCepstrum && this->Para.DEnergy != 0) {
			//for (frameCnt = pData->Row-1; frameCnt >= 1; frameCnt--) {
			//	pData->Mat[frameCnt][0] = pData->Mat[frameCnt][0] - pData->Mat[frameCnt-1][0];
			//}
			matFunc.differentiate(pData->Mat, pData->Row, pData->Col, 0);
		}

		MFree_1D(frame);
		MFree_1D(fbe);
	}
	MFree_0D(pTrans);
	MFree_2D(filterParam);

  return(pData);
}

//====================================================================================================================
// evaluates the specified filter for the given spectrum. filterCenter is the filter's center frequency and 
// filterMin/-Max or the left and right point where the filter reaches zero, all as returned by createSpacing()
//
// this is taken from m.slaneys auditory toolbox for matlab (calculation of mfccFilterWeights) with some changes: 
//  - there, frequency = (idx * sampleRate) / fftSize, not frequency = (idx * sampleRate) / (2 * fftSize)
//  - in the rising branch, we ask for frequency >= filterMin &&..., there it is just frequency = filterMin &&..
//    (this misses indices!)
//====================================================================================================================
double SC_Feature_FbE::evalFilter(double* powerSpec, double filterMin, double filterCenter, double filterMax, double sampleRate, unsigned long int fftSize) {
  double frequency, converter = sampleRate / (double)(fftSize), filterValue = 0.0;
  double height = 2.0 / (filterMax - filterMin), normRising = filterCenter - filterMin, normFalling = filterMax - filterCenter;

  for (unsigned long int idx = (unsigned long)ceil(filterMin/converter); idx < fftSize/2; idx++) {
    frequency = idx * converter;
    if (frequency >= filterMin && frequency <= filterCenter) {
      filterValue += powerSpec[idx] * ((height * (frequency - filterMin)) / normRising);
    } else if (frequency > filterCenter && frequency < filterMax) {
      filterValue += powerSpec[idx] * ((height * (filterMax - frequency)) / normFalling);
    } else if (frequency > filterMax) {
      break;
    }
  }

  return filterValue;
}

//====================================================================================================================
// returns the filter coefficient for the given triangular filter and the proposed idx into a (power-)spectrum. just
// the same routines as in the other evalFilter() method.
//====================================================================================================================
double SC_Feature_FbE::evalFilter(unsigned int idx, double filterMin, double filterCenter, double filterMax, double sampleRate, unsigned long int fftSize) {
  double frequency, converter = sampleRate / (double)(fftSize), filterValue = 0.0;
  double height = 2.0 / (filterMax - filterMin), normRising = filterCenter - filterMin, normFalling = filterMax - filterCenter;

  frequency = idx * converter;
  if (frequency >= filterMin && frequency <= filterCenter) {
    filterValue = (height * (frequency - filterMin)) / normRising;
  } else if (frequency > filterCenter && frequency < filterMax) {
    filterValue = (height * (filterMax - frequency)) / normFalling;
  }

  return filterValue;
}

//====================================================================================================================
// prints the current filterbank-spacing to an ASCII-file
//====================================================================================================================
void SC_Feature_FbE::spacingOut(void) {
  spacingOut(this->scale, this->lowFrequency, this->highFrequency);

  return;
}

//====================================================================================================================
// prints the desired filterbank-spacing to an ASCII-file
// the file contains a matrix with as many rows as int.-frequencys till the maximum frequency in the desired scale.
// each col represents one filter
//====================================================================================================================
void SC_Feature_FbE::spacingOut(unsigned char mode, double fminHz, double fmaxHz) {
  double **filterParam = createSpacing(mode, fminHz, fmaxHz, this->Para.NFilter);
  double **scale;
  double delta_f = (double)(this->Para.SRate) / (double)(this->Para.FFTSz);
  long int idx, idx_min, idx_max, Idx = (long int)(floor(filterParam[this->Para.NFilter-1][2]));
  long int f;

  MArray_2D(scale, (long)(this->Para.NFilter), Idx, double, "SC_Feature_FbE.spacingOut: scale");

  for (idx = 0; idx < Idx; idx++) {
    for (f = 0; f < this->Para.NFilter; f++) {
      scale[f][idx] = 0.0;
    }
  }

  for (f = 0; f < this->Para.NFilter; f++) {
    idx_min = sclib::max((long int)filterParam[f][0], 0);
    idx_max = sclib::min((long int)filterParam[f][2], Idx - 1);


    for (idx = idx_min; idx <= idx_max; idx++) {
      if (idx <= filterParam[f][1]) {
        scale[f][idx] = (1.0 / (filterParam[f][1] - filterParam[f][0])) * (idx - idx_min);
      } else {
        scale[f][idx] = (-1.0 / (filterParam[f][2] - filterParam[f][1])) * (idx - filterParam[f][1]) + 1.0;
      }
    }
  }

  sclib::matrixOut("filterBank.txt", scale, this->Para.NFilter, Idx);

  MFree_2D(filterParam);
  MFree_2D(scale);

  return;
}
