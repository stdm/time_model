/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class that unites and encapsulates many algorithms for sound    */
/*        synthesis                                                       */
/*																																				*/
/*    Author  : Thilo Stadelmann															            */
/*    Date    : 12.06.2008																								*/
/**************************************************************************/

#include "SC_Synthesis.h"
#include "SC_Signal_WAVE.h"
#include "SC_MatrixFunctions.h"
#include "SC_FeatureHandler.h"
#include "SC_Feature_FbE.h"
#include <GN_WaveGen.h>

//====================================================================================================================
// constructor / destructor
//====================================================================================================================
SC_Synthesis::SC_Synthesis(SC_TweakableParameters *pTweak, unsigned int fftLength, double taperingLength, unsigned int taperingMode, unsigned int olaMaxIterations, double olaErrortarget, bool verbose) : olaMaxIterations(olaMaxIterations), olaErrorTarget(olaErrorTarget), verbose(verbose) {
	this->pTweak = pTweak;
	this->pTransform = new SC_Transform(fftLength, taperingMode, taperingLength, sclib::modeSClib);
}

SC_Synthesis::~SC_Synthesis(void) {
	MFree_0D(this->pTransform);
}

//=====================================================================================================================
// implementation of a convolutive (i.e. FIR) filter via the fft: convolution in the time-domain is multiplication in 
// the fft-donmain, this is why the inputs are the signal's and filter-kernel's (i.e. filter's impulse response's) 
// spectra in polar form (for convenience, can be given in rectangular form if spectraArePolar==false, then both mag 
// and phase array sa re meant to have length usedFftLen/2+1, else phase is meant to habve the full usedFftLen); the 
// original length figures must be given in order for the function to decides (and warn in case) if circular 
// convolution happened or everything went fine; the output signal has length 
// originalSignalLength+originalFilterLength-2; see dspguide book, pp. 183-184 and pp. 311-318.
// NULL is returned in case of error (e.g. circular convolution).
// if filterPhase==NULL all filterPhase[i] are assumed to be equal to 0.0 (i.e. zero phase) for convenience
//=====================================================================================================================
double* SC_Synthesis::fftFilter(double *signalMagnitude, double *signalPhase, unsigned int originalSignalLen, double *filterMagnitude, double *filterPhase, unsigned int originalFilterLen, unsigned int usedFftLen, unsigned int &resultLen, bool spectraArePolar) {
	double *newSignal = NULL, *newMagnitude, *newPhase, *tmp;
	unsigned int i, spectraLen = usedFftLen/2+1;
	unsigned long int oldFftLen = this->pTransform->getFftLen();
	
	resultLen = originalFilterLen+originalSignalLen-2;	
	if (resultLen <= usedFftLen) {
		this->pTransform->setFftLen(usedFftLen); //just to be sure...
	
		//do the filtering in the frequency domain
		MArray_1D(newMagnitude, (int)((spectraArePolar==true)?spectraLen:usedFftLen), double, "SC_Synthesis.fftFilter: newMagnitude");
		MArray_1D(newPhase, (int)(usedFftLen), double, "SC_Synthesis.fftFilter: newPhase");
		if (spectraArePolar == true) { //polar form
			for (i = 0; i < usedFftLen; i++) {
				if (i < spectraLen) {
					newMagnitude[i] = signalMagnitude[i] * filterMagnitude[i]; //dspguide, p.182
				}
				if (filterPhase != NULL) {
					newPhase[i] = signalPhase[i] + filterPhase[i]; 
				} else { //assume all filterPhase[i]==0.0 for filterPhase==NULL
					newPhase[i] = signalPhase[i];
				}
			}
		} else { //rectangular form
			for (i = 0; i < usedFftLen; i++) {
				if (filterPhase != NULL) {
					newMagnitude[i] = signalMagnitude[i]*filterMagnitude[i] - signalPhase[i]*filterPhase[i]; //dspguide, equation 9-1, p.182
					newPhase[i] = signalPhase[i]*filterMagnitude[i] + signalMagnitude[i]*filterPhase[i];
				} else { //assume all filterPhase[i]==0.0 for filterPhase==NULL
					newMagnitude[i] = signalMagnitude[i]*filterMagnitude[i];
					newPhase[i] = signalPhase[i]*filterMagnitude[i];
				}
			}
		}
		
		//transform back to time domain
		if (spectraArePolar == true) {
			this->pTransform->magnitudePhase2ft(newMagnitude, newPhase, false, false, false);
		}
		tmp = this->pTransform->ifft(newMagnitude, newPhase, usedFftLen);
		MFree_1D(newMagnitude);
		MFree_1D(newPhase);
		
		//store the result in a container of proper size
		if (usedFftLen > resultLen) {
			MArray_1D(newSignal, (int)(resultLen), double, "SC_Synthesis.fftFilter: newSignal"); //realloc the new signal because its shorter!
			for (i = 0; i < resultLen; i++) {
				newSignal[i] = tmp[i];
			}
			MFree_1D(tmp);	
		}  else {
			newSignal = tmp;
		}
		
		this->pTransform->setFftLen(oldFftLen); //restore previous state to avoid side effects
	} else {
		resultLen = 0;
	}
	
	return newSignal;
}

//=====================================================================================================================
// returns the filter coefficients for a filter with the desired frequency/phase response; see Steven W. Smith,
// "Digital Signal Processing - A Practical Guide for Engineers and Scientists", pp. 297-300. 
// The last 3 parameters are for the reconstruction of the fft from magnitude/phase and should give the 
// parameters used to construct the desired* parameters. The result is a filter kernel (or impulse response) that can 
// directly be used with the convolutiveFilter() method.
//=====================================================================================================================
double* SC_Synthesis::constructCustomFilter(double *desiredMagnitudeResponse, double *desiredPhaseResponse, int responseLength, int filterKernelLength, bool isFullLength, bool isDegrees, bool isLogarithmized) {
	double *desiredImpulseResponse, *filter, *xr, *xi;
	double oldTaperingLength = this->pTransform->getTaperingLength();
	unsigned int oldTaperingMode = this->pTransform->getTaperingMode();
	int i, M = filterKernelLength-1, T = (responseLength-1)*2;
	
	//switch off any unintentional signal processing stuff
	this->pTransform->setTaperingLength(0.0);
	this->pTransform->setTaperingMode(sclib::wndRectangle);
	
	//frequency-/phase-response => impulse-response
	MArray_1D(xr, responseLength, double, "SC_Synthesis.constructCustomFilter: xr");
	MArray_1D(xi, T, double, "SC_Synthesis.constructCustomFilter: xi");
	for (i = 0; i < T; i++) {
		if (i < responseLength) {
			xr[i] = desiredMagnitudeResponse[i];
		}
		xi[i] = desiredPhaseResponse[i];
	}
	this->pTransform->magnitudePhase2ft(xr, xi, isDegrees, isFullLength, isLogarithmized);
	desiredImpulseResponse = this->pTransform->ifft(xr, xi, this->pTransform->getFftLen()); 
	MFree_1D(xr);
	MFree_1D(xi);
		
	//shift (rotate) the signal M/2 points to the right
	sclib::rotateArray(desiredImpulseResponse, this->pTransform->getFftLen(), M/2);

	//truncate and window the signal, copy it to the output buffer
	MArray_1D(filter, filterKernelLength, double, "SC_Synthesis.constructCustomFilter: filter");
	for (i = 0; i < filterKernelLength; i++) {
		filter[i] = desiredImpulseResponse[i] * (0.54 - 0.46*cos(sclib::two_pi*(double)(i)/(double)(M))); //hamming-window
	}
	MFree_1D(desiredImpulseResponse);

	//restore class parameters	
	this->pTransform->setTaperingLength(oldTaperingLength);
	this->pTransform->setTaperingMode(oldTaperingMode);
	
	return filter;	
}

//=====================================================================================================================
// normalize the given signal to relativeMaxAmplitude% max. amplitude and write it to disc
//=====================================================================================================================
bool SC_Synthesis::writeSignal(double *dblSamples, unsigned long int sampleCount, unsigned int sampleRate, const char *fileName, double relativeMaxAmplitude) {
	unsigned long int s;
	bool res = true;
	short *samples;
	SC_Signal_WAVE *pSignal;
	double maxValue = 0.0;

	//find maximum amplitude value
	MArray_1D(samples, sampleCount, short, "SC_Synthesis.writeSignal: samples");
	for (s = 0; s < sampleCount; s++) {
		if (fabs(dblSamples[s]) > maxValue) {
			maxValue = fabs(dblSamples[s]); 
		}
	}

	//normalize to relativeMaxAmplitude% max. amplitude
	for (s = 0; s < sampleCount; s++) {
		samples[s] = (short)(sclib::round(relativeMaxAmplitude * std::numeric_limits<short>::max() * dblSamples[s]/maxValue)); //normalize & copy
	}

	//write the signal
	pSignal = new SC_Signal_WAVE();
	pSignal->setHeader(0, 1, 1, sampleRate, 16, 0);
	pSignal->setBuf_L(samples, sampleCount);
	if (pSignal->SaveSignal(fileName) <= 0) {
		res = false;
	}

	//clean up
	MFree_0D(pSignal); //also kills the samples buffer

	return res;
}

//=====================================================================================================================
// receives a spectral frame and adds pitch to it; this->pTransform->getFftLen() must return the correct spectrum size
// used to create this spectrum
//=====================================================================================================================
void SC_Synthesis::addPitch(double *spectrum, unsigned long int spectrumNr, int spectrumFrameSize, int spectrumFrameStep, SC_Conversion *pConverter, SV_Data *pPitch) {
	double deviation, relativeDeviation;
	
	if (pPitch != NULL) {
		//sclib::vectorOut("spect.txt", spectrum, this->pTransform->getFftLen()/2 +1, false, this->pTweak);
		for (unsigned long int f = pConverter->sample2audioFrame(pConverter->audioFrame2sample(spectrumNr, spectrumFrameSize, spectrumFrameStep), pPitch->Hdr.frameSize, pPitch->Hdr.frameStep); f <= sclib::min(pConverter->sample2audioFrame(pConverter->audioFrame2sample(spectrumNr, spectrumFrameSize, spectrumFrameStep, sclib::alignmentEnd), pPitch->Hdr.frameSize, pPitch->Hdr.frameStep), pPitch->Row-1); f++) { //for all pitch-frames in the range of this lpc-frame
			if (pPitch->Mat[f][0] > 0.0) { //only add pitch for voiced sounds
				for (unsigned long int d = 0; d < this->pTransform->getFftLen()/2 +1; d++) {
					if (d > 0) { //amplify/deamplify the spectral envelope up to 25% according to the distance of each bin to the next harmonic of f0; do nothing with the "dc coefficient"
						deviation = pConverter->index2frequency(d, this->pTransform->getFftLen()) / pPitch->Mat[f][0];
						relativeDeviation = fabs((double)(sclib::round(deviation)) - deviation); //this is now the relative deviation (between 0..0.5) of the current fft-bin's frequency from the next harmonic of f0
						spectrum[d] *= 1.25 - relativeDeviation; //this way, bins falling on a harmonic are amplified by 25%, and those bins maximally distant from any harmonic are deamplified by 25%
						//sclib::scalarOut("reldev.txt", 1.25 - relativeDeviation, this->pTweak, false, " ");
					}
				} //for all fft bins
				//sclib::stringOut("reldev.txt", "\n", this->pTweak, "");
				break; //so that different pitches do not overlap in one spectrum
			} //if voiced
		} //for all pitch frames in one lpc frame
		//sclib::vectorOut("pitch.txt", spectrum, this->pTransform->getFftLen()/2 +1, false, this->pTweak);
	} //if pitch is given

	return;
}

//=====================================================================================================================
// This methods receives an array of power-spectra (the number of spectra is automatically decuced from the 
// additionally given parameters frameSize, frameStep and originalSignalLength, which refers to the number of samples 
// the signal had from wich the spectral frames where computed and hence is also the number of returned samples). if 
// powerIsMagnitude is true, it also works with an array of magnitude-spectra. it is important that the same class 
// parameters (analysis window type and length, fft-length) for this class where used during specrum creation as now.
// the method works best when a hamming window is used due to it being non-null on all its coefficients. the method
// used to resynthesize the signal (i.e. reconstruct the pase) is a variant of overlap&add (OLA) called LSEE-MSTFTM
// (Least Square Error Estimation from Modified Short Time Fourier Transform Magnitude) due to Griffin, Lim, "Signal 
// Estimation from Modified Short-Time Fourier Transform", 1984
//=====================================================================================================================
double* SC_Synthesis::olaSynthesis(double **powerSpectra, unsigned int frameSize, unsigned int frameStep, unsigned long int originalSignalLength, bool powerIsMagnitude, double **originalPhase, unsigned int maxIterations, double errorAim) {
	unsigned int i, j, p, spectrumCount = sclib::getRowCount(originalSignalLength, frameSize, frameStep);
	unsigned long int reachedSamples = (spectrumCount-1)*frameStep + frameSize;
	double *initialPhaseSpec = NULL, *currentFrame = NULL, *newSignal = NULL, *currentSignal = NULL, *currentPhaseSpec = NULL, *currentMagnitudeSpec, *divisor = NULL, *window, *tmp;
	GN_WaveGen *pWaveGenerator = new GN_WaveGen();
	SC_MatrixFunctions matFunc;
	double error, wantedValue, averageMagnitude = 0.0;

	//initialize divisor for Griffin&Lim's equation 6 (sum of all the products of the analysis- and synthesis-window-coefficients for each sample)
	divisor = matFunc.zeros(originalSignalLength);
	window = matFunc.ones(frameSize);
	this->pTransform->tapering(window, frameSize, this->pTransform->getTaperingMode(), this->pTransform->getTaperingLength(), false, false); //this way, we get the window coefficients and regard all parameters, even the tapering-length, assuming that they have been used for the creation of the input spectra
	for (p = 0; p < spectrumCount; p++) {
		for (i = 0; i < frameSize; i++) {
			divisor[p*frameStep + i] += window[i] * window[i]; //product of coefficients of analysis- and synthesis windows, summed over all frames that include this sample
		}
	}

	//initialization: generate white noise and compute its phase spectrum
	tmp = pWaveGenerator->NoiseUni_Wave(frameSize);
	currentMagnitudeSpec = this->pTransform->magnitudeSpectrum(tmp, frameSize, initialPhaseSpec, false, false);
	MFree_0D(pWaveGenerator); //this kills the buffer pointed to by tmp, too
	MFree_1D(currentMagnitudeSpec);

	MArray_1D(currentSignal, originalSignalLength, double, "SC_Synthesis.olaSynthesis: currentSignal");
	MArray_1D(newSignal, originalSignalLength, double, "SC_Synthesis.olaSynthesis: newSignal");

	i = 0;
	do {
		//initialize the newSignal Buffer with zeros
		for (j = 0; j < originalSignalLength; j++) {
			newSignal[j] = 0.0;
		}
		error = 0.0;

		for (p = 0; p < spectrumCount; p++) {
			if (i > 0) { //in the 0th run, we start with the spectra and resynthesize the 0th currentSignal instead of computing spectra from (i-1)th currentSignal first...
				//fill current frame
				MArray_1D(currentFrame, frameSize, double, "SC_Synthesis.olaSynthesis: currentFrame");
				for (j = 0; j < frameSize; j++) {
					currentFrame[j] = currentSignal[p*frameStep + j];
				}

				//compute current power- and phase spectra
				currentMagnitudeSpec = this->pTransform->magnitudeSpectrum(currentFrame, frameSize, currentPhaseSpec, false, false);
				MFree_1D(currentFrame);

				//replace current power spectrum with the desired one
				for (j = 0; j < this->pTransform->getFftLen()/2 +1; j++) {
					wantedValue = (powerIsMagnitude==false) ? sqrt(powerSpectra[p][j]*(double)(frameSize)) : powerSpectra[p][j];
					error += fabs(wantedValue - currentMagnitudeSpec[j]); //this is the sum of all error magnitude
					currentMagnitudeSpec[j] = wantedValue;
				}
			} else { //...but instead we have to copy the initial phase guess to the right buffer on th 0th run; to save time, we combine it with the loop to copy the power
				error = std::numeric_limits<double>::max(); //just to not leave the loop too early
				MArray_1D(currentPhaseSpec, this->pTransform->getFftLen(), double, "SC_Synthesis.olaSynthesis: currentPhaseSpec"); //on the 0th run, it doesn't exist before
				MArray_1D(currentMagnitudeSpec, this->pTransform->getFftLen(), double, "SC_Synthesis.olaSynthesis: currentMagnitudeSpec"); //dito
				for (j = 0; j < this->pTransform->getFftLen(); j++) {
					currentPhaseSpec[j] = (originalPhase==NULL) ? initialPhaseSpec[j] : originalPhase[p][j];
					if (j < this->pTransform->getFftLen()/2 +1) {
						currentMagnitudeSpec[j] = (powerIsMagnitude==false) ? sqrt(powerSpectra[p][j]*(double)(frameSize)) : powerSpectra[p][j];
						averageMagnitude += currentMagnitudeSpec[j] / (double)(spectrumCount);
					}
				}
			}
	
			//convert power/phase spectrum to real/imaginary part (xr/xi) of ft
			this->pTransform->magnitudePhase2ft(currentMagnitudeSpec, currentPhaseSpec, false, false, false);
					
			//convert xr/xi to signal
			currentFrame = this->pTransform->ifft(currentMagnitudeSpec, currentPhaseSpec, this->pTransform->getFftLen());
			MFree_0D(currentMagnitudeSpec);
			MFree_0D(currentPhaseSpec);

			//window the resynthesized frame with the synthesis window (same as analysis window here according to Verhelst, 
			//"Overlap-Add Methods for Time-Scaling of Speech", 1999), then pack the just resynthesized frame into the newSignal 
			//buffer (to not overwrite the currentSignal which is in parts still needed by the next frames) by using the overlap&add 
			//strategy and formula 6 from Griffin, Lim, "Signal Estimation from Modified Short-Time Fourier Transform", 1984
			for (j = 0; j < frameSize; j++) {
				newSignal[p*frameStep + j] += window[j] * currentFrame[j];
			}
			MFree_1D(currentFrame);
		} //for all powerSpectra

		error = ((error/(double)(spectrumCount)) * 100.0) / averageMagnitude;
		if (this->verbose == true) {
			if (i > 0) {
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%5d: %#8.4f%%", i, error); //this error figure is the average error-magnitude expressed in percent of the average overall spectral magnitude
			} else {
				printf("\nOLA average magnitude error in round     0:       - %%");
			}
		}
		
		//kill the initial phase spectrum after the 0th round
		if (i == 0) {
			MFree_1D(initialPhaseSpec);
		}

		//finish the resynthesized ith signal by dividing all samples by the previously initialized (constant) divisor
		for (j = 0; j < reachedSamples; j++) {
			if (divisor[j] > 0.0) { //avoid div. by zero; it doesn't do any harm because when the divisor is zero, this means, one of the window-coefficients was zero, too, so the numerator will also be zero!
				newSignal[j] /= divisor[j];
			}
		}
		
		//switch the buffers used for the current and new signal in the next round
		tmp = currentSignal;
		currentSignal = newSignal; //currentSignal holds the current estimate of the reconstructed signal now
		newSignal = tmp;
	} while (i++ < maxIterations && error >= errorAim);
	
	if (this->verbose == true) {
		printf("\n");
	}
	
	MFree_1D(window);
	MFree_1D(currentFrame);
	MFree_1D(divisor);
	MFree_1D(newSignal);

	//add zeros for those samples that didn't fit into a last frame; we can't do any better!
	for (j = reachedSamples; j < originalSignalLength; j++) {
		currentSignal[j] = 0.0;
	}

	return currentSignal;
}

//====================================================================================================================
// creates a .wav-file of the given name containing a "hummed" melody following the pitch given in the first column of
// pPitch; if addHarmonics!=NULL, the (i+1)th harmonics (i=0 being the fundamental) is added to the hum if the array
// has the value 'true' at index-position i.
//====================================================================================================================
void SC_Synthesis::generateHumFromPitch(char *fileName, SV_Data *pPitch, bool *addHarmonics, unsigned int numberOfHarmonics, bool addVoicing) {
	short *samples = NULL;
	unsigned long int sampleCount;
	long int s, start, end, count = 0, lastT = 0;
	int halfStep = pPitch->Hdr.frameStep/2; 
	int halfSize = pPitch->Hdr.frameSize/2; 
	double yValue = 0.0, xValue = 0.0, sampleRate, arcsin = 0.0, lastFrequency = 0.0, maxValue = 0.0, maxAmplitude = 0.0;
	bool finishLastPeriod = false, *harmonics;
	unsigned int harmonicCount, h;
	int t, filterSize = 15;
	double *buffer, gain = 17393.5, lpcFilter[] =  {1.0, -1.96331, 1.63752, -1.15305, 0.942088, -1.01071, 1.75079, -2.047, 1.77946, -1.06579, 0.374125, -0.388501, 0.630771, -0.537089, 0.271788}; //lpc-coefficients for thilo's vocal tract in vowel-"a"-configuration
														//gain: 5074.53 lpc: {1.0, 1.55142, 0.824962, -1.0087, 1.37821, -1.38058, 1.21578, -0.690426, 0.787107, -1.0046, 1.35558, -1.16207, 0.450235, -0.701326, 0.566267}; //lpc-coefficients for thilo's vocal tract in vowel-"u"-configuration
	                                             
	if (pPitch != NULL && pPitch->Row > 0) {
		if (addHarmonics == NULL) {
			MArray_1D(harmonics, 1, bool, "SC_Synthesis.genaretaHumFromPitch: harmonics");
			harmonics[0] = true;
			harmonicCount = 1;
			maxAmplitude = 1.0;
		} else {
			harmonics = addHarmonics;
			harmonicCount = numberOfHarmonics;
			for (h = 0; h < harmonicCount; h++) {
				if (harmonics[h] == true) {
					maxAmplitude += 1.0;
				}
			}
		}

		sampleRate = (double)(pPitch->Hdr.sampleRate);
		
		sampleCount = ((pPitch->Row-1)*pPitch->Hdr.frameStep) + pPitch->Hdr.frameSize;
		MArray_1D(buffer, sampleCount, double, "SC_Synthesis.generateHumFromPitch: buffer");

		for (t = 0; t < pPitch->Row; t++) {
			start = (t > 0) ? halfSize-halfStep : 0;
			end = (t < pPitch->Row-1) ? halfSize+halfStep : pPitch->Hdr.frameSize;
			if (pPitch->Mat[t][0] > 0.0) {
				//find the correct x-value so that after possibly changing the frequency here the waveform is still continuous (without jumps)
				yValue = 0.0;
				for (h = 0; h < harmonicCount; h++) { //avoid aliasing by cancelling harmonics with frequency above the nyquist frequency
					if (harmonics[h] == true && (h+1)*pPitch->Mat[t][0] < sampleRate/2.0) {
						yValue += (h+1)*pPitch->Mat[t][0];
					}
				}
				xValue = (sampleRate * arcsin) / (yValue * sclib::two_pi) + 1;

				for (s = start; s < end; s++) {
					yValue = 0.0;
					arcsin = 0.0;
					for (h = 0; h < harmonicCount; h++) {
						if (harmonics[h] == true && (h+1)*pPitch->Mat[t][0] < sampleRate/2.0) { //avoid aliasing by cancelling harmonics with frequency above the nyquist frequency
							yValue += sin(sclib::two_pi * (h+1)*pPitch->Mat[t][0] * xValue/sampleRate);
							arcsin += sclib::two_pi * (h+1)*pPitch->Mat[t][0] * xValue/sampleRate; //needed to find correct xValue for the next frrequency later on
						}
					}
					buffer[count++] = 1.0/maxAmplitude * std::numeric_limits<short>::max() * yValue;
					xValue += 1.0;
				}
				finishLastPeriod = true;
				lastFrequency = pPitch->Mat[t][0];
				lastT = count;
			} else {
				for (s = start; s < end; s++) {
					if (finishLastPeriod == true) { //finnish the last cycle of a non-null pitch till it reaches zero (or changes sign), then continue with a 0-signal
						yValue = 0.0;
						arcsin = 0.0;
						for (h = 0; h < harmonicCount; h++) {
							if (harmonics[h] == true && (h+1)*lastFrequency < sampleRate/2.0) { //avoid aliasing by cancelling harmonics with frequency above the nyquist frequency
								yValue += sin(sclib::two_pi * (h+1)*lastFrequency * xValue/sampleRate);
								arcsin += sclib::two_pi * (h+1)*lastFrequency * xValue/sampleRate; //needed to find correct xValue for the next frrequency later on
							}
						}
						buffer[count++] = 1.0/maxAmplitude * std::numeric_limits<short>::max() * yValue;
						xValue += 1.0;
						if (sclib::sg(buffer[count-2]) != sclib::sg(buffer[count-1]) || buffer[count-1] == 0) {
							buffer[count-1] = 0.0;
							finishLastPeriod = false;
						}
					} else {
						buffer[count++] = 0.0;
						arcsin = 0.0; //needed to find correct xValue for the next frrequency later on
					}
				}
				yValue = 0.0;
			}
		}	

		//shape with lpc (vocal tract) filter to sound like a voice
		if (addVoicing == true) {
			buffer = sclib::directOutput(recursiveFilter(buffer, sampleCount, &gain, 1, lpcFilter, filterSize, maxValue, false, false), buffer);
		}

		//normalize & save signal
		if (writeSignal(buffer, count, pPitch->Hdr.sampleRate, fileName, 0.3) == false) {
			REPORT_ERROR(SVLIB_FileErr, "Couldn't write humming to given fileName");
		}	
		
		//clean up
		MFree_1D(buffer);
		if (harmonics != addHarmonics) {
			MFree_1D(harmonics);
		}
	}

	return;
}

//====================================================================================================================
// take the samples of the internal buffer, smooth spectrum by replacing the amplitude for each frequency by the mean 
// in its critical band and resynthesize the result as a wav-file
//====================================================================================================================
void SC_Synthesis::generateVoiceFromCB(char *fileName, short *samples, unsigned long int sampleCount, unsigned int sampleRate, unsigned int frameSize, unsigned int frameStep) {
	SC_Signal_WAVE *pVoice;
	short *newSamples;
  double *frame = NULL, *newFrame = NULL;
  double *powerSpectrum, *phaseSpectrum = NULL, *newPowerSpectrum = NULL;
  unsigned long int f, s, frameCount, b;
	double centerFrequency, leftFrequency, rightFrequency, bandwidth, mean, denominator;
	SC_Conversion converter(sampleRate);
	int bb, **bandBorders = NULL;
	double **spectra = NULL, *tmp;

	frameCount = sclib::getRowCount(sampleCount, frameSize, frameStep);
	MArray_1D(spectra, frameCount, double*, "SC_Synthesis.generateVoiceFromCB: spectra");

	//prepare a table that holds, for each frequency bin, the left and right indices of bins belonging to its band according to the critical bandwidth
	MArray_2D(bandBorders, (int)(this->pTransform->getFftLen()/2 +1), 2, int, "SC_Synthesis.generateVoiceFromCB: bandBorders");
	for (b = 0; b < (unsigned long int)(this->pTransform->getFftLen()/2 +1); b++) {
		centerFrequency = converter.index2frequency(b, this->pTransform->getFftLen());
		bandwidth = SC_Feature_FbE::getCriticalBandwidth(centerFrequency, centerFrequency, leftFrequency, rightFrequency, false);
		bandBorders[b][0] = converter.frequency2index(sclib::max(0, leftFrequency), this->pTransform->getFftLen());
		bandBorders[b][1] = converter.frequency2index(sclib::min(sampleRate/2, rightFrequency-0.0001), this->pTransform->getFftLen());
	}

	MArray_1D(frame, frameSize, double, "SC_Synthesis.generateVoiceFromCB: frame");
	for (f = 0; f < frameCount; f++) {
		for (s = 0; s < (unsigned long)(frameSize); s++) {
			frame[s] = samples[f*frameStep + s];
		}

		//compute the power spectrum for the current frame
		powerSpectrum = this->pTransform->powerSpectrum(frame, frameSize, false, false);

		//replace the energy value in each bin by the mean of all energies in the corresponding critical band
		MArray_1D(newPowerSpectrum, this->pTransform->getFftLen()/2 +1, double, "SC_Synthesis.generateVoiceFromCB: newPowerSpectrum");
		for (b = 0; b < (unsigned long int)(this->pTransform->getFftLen()/2 +1); b++) {
			mean = 0.0;
			denominator = (double)(bandBorders[b][1] - bandBorders[b][0] + 1.0);
			for (bb = bandBorders[b][0]; bb <= bandBorders[b][1]; bb++) {
				mean += powerSpectrum[bb] / denominator;
			}
			newPowerSpectrum[b] = mean;
		}
		MFree_1D(powerSpectrum);

		spectra[f] = newPowerSpectrum;
	}
	MFree_1D(frame);
	MFree_2D(bandBorders);

	//resynthesize the signal with phase-reconstruction by OLA method
	MArray_1D(newSamples, sampleCount, short, "SC_Synthesis.generateVoiceFromCB: newSamples");
	tmp = olaSynthesis(spectra, frameSize, frameStep, sampleCount, false, NULL, this->olaMaxIterations, this->olaErrorTarget);
	for (s = 0; s < sampleCount; s++) {
		newSamples[s] = (short)(sclib::round(tmp[s]));
	}
	MFree_1D(tmp);

	//free the spectral frames
	MFree_2Dex(spectra, frameCount);

	//write the signal
	pVoice = new SC_Signal_WAVE();
	pVoice->setHeader(0, 1, 1, sampleRate, 16, 0);
	pVoice->setBuf_L(newSamples, sampleCount);
	if (pVoice->SaveSignal(fileName) <= 0) {
		REPORT_ERROR(SVLIB_FileErr, "Couldn't write humming to given fileName");
	}
	MFree_0D(pVoice); //also kills the newSamples buffer

	return;
}

//====================================================================================================================
// uses olaSynthesis() to convert a series of (power) spectra back to a signal stored under the given file-name as a 
// wav file; if smooth==true, the spectra are smoothed before concatenating (i.e. 20% highest frequencies in the 
// spectral envelope are removed); if deLog==true, the spectra are assumed to be log-spectra and exponentialized 
// before resynthesis; if magnitudeIsPower==true, the spectra are assumed to be power-spectra instead of magnitude 
// spectra
//====================================================================================================================
void SC_Synthesis::spectrum2wav(char *fileName, SV_Data* pSpectrum, bool smooth, bool deLog, bool magnitudeIsPower, bool keepLimitations) {
	unsigned long int sampleCount = (pSpectrum->Row-1)*pSpectrum->Hdr.frameStep + pSpectrum->Hdr.frameSize;
	double **spectra = NULL, *dblSamples, **phase = NULL;
	int x, y, M = (unsigned int)floor(this->pTransform->getFftLen() - this->pTransform->getFftLen()/20.0); //remove the 20% highest frequencies in the spectral envelope
	bool hasPhase = (pSpectrum->Next!=NULL && pSpectrum->Next->Hdr.ID==sclib::featureSpectrum && pSpectrum->Next->Hdr.Signature[1]==1) ? true : false;

	//copy the data to a double matrix needed by olaSynthesis()
	MArray_1D(spectra, pSpectrum->Row, double*, "SC_Synthesis.spectrum2wav: spectra");
	if (hasPhase == true) {
		MArray_1D(phase, pSpectrum->Row, double*, "SC_Synthesis.spectrum2wav: phase");
	}
	for (y = 0; y < pSpectrum->Row; y++) {
		MArray_1D(spectra[y], pSpectrum->Col, double, "SC_Synthesis.spectrum2wav: spectra[y]");
		if (hasPhase == true) {
			MArray_1D(phase[y], pSpectrum->Next->Col, double, "SC_Synthesis.spectrum2wav: phase[y]");
		}
		for (x = 0; x < ((hasPhase==true)?pSpectrum->Next->Col:pSpectrum->Col); x++) {
			if (x < pSpectrum->Col) {
				spectra[y][x] = pSpectrum->Mat[y][x];
			}
			if (hasPhase == true) {
				phase[y][x] = pSpectrum->Next->Mat[y][x];
			}
		}
		if (smooth == true) {
			this->pTransform->smoothSpectrum(spectra[y], false, M, true);
		}
		if (deLog == true) {
			for (x = 0; x < pSpectrum->Col; x++) {
				spectra[y][x] = exp(spectra[y][x]);
			}
		}
	}

	//resynthesize the signal with phase-reconstruction by OLA method
	dblSamples = olaSynthesis(spectra, pSpectrum->Hdr.frameSize, pSpectrum->Hdr.frameStep, sampleCount, !magnitudeIsPower, phase, this->olaMaxIterations, this->olaErrorTarget);
	if (writeSignal(dblSamples, sampleCount, pSpectrum->Hdr.sampleRate, fileName, 0.3) == false) {
		REPORT_ERROR(SVLIB_FileErr, "Couldn't write resynthesized power spectra to given fileName");
	}

	//free the spectral frames and dbl-signal
	MFree_1D(dblSamples);
	MFree_2Dex(spectra, pSpectrum->Row);
	MFree_2Dex(phase, pSpectrum->Row);

	return;
}

//====================================================================================================================
// gets a set of mfcc vectors (and maybe corresponding pitch) and converts them back to spectra; the result-set is not 
// a matrix but an array pointing to #pMFCC->Row single spectrum arrays of length usedFftLen/2; the given parameters 
// mainly tell how the mfccs where created in order to correctly invert the process; references for ideas used in this 
// method are: Milner, Shao, "Clean Speech Reconstruction from MFCC Vectors and Fundamental Frequency using an 
// Integrated Front-End", 2005; Milner, Shao, "Speech Reconstruction from Mel-Frequency Cepstral Coefficients using a 
// Source-Filter Model", 2002
//====================================================================================================================
double** SC_Synthesis::mfcc2spectrum(SV_Data* pMFCC, double usedPreEmphasisFactor, unsigned int usedWindow, unsigned int usedFftLen, int usedFrequencyScale, double usedFilterbankLowEnd, double usedFilterbankHighEnd, unsigned int usedFilterbankSize, bool firstMfccIsDEnergy, SV_Data *pPitch) {
	double **spectra = NULL, *spectrum = NULL, *cepstrum = NULL, *filterOutput, *divisor = NULL;
	double **filterParam = SC_Feature_FbE::createSpacing(usedFrequencyScale, usedFilterbankLowEnd, (usedFilterbankHighEnd>0.0)?usedFilterbankHighEnd:pMFCC->Hdr.sampleRate/2.0, usedFilterbankSize);
	double lastMFCC0 = 0.0, filterValue;
	double preEmphasisImpulseResponse[] = {1.0, -1.0*usedPreEmphasisFactor}, *preEmphasisFrequencyResponse, *preEmphasisCepstrum, *filterbankAreaCepstrum;
	int t;
	unsigned int f, d, **idx = NULL, oldTaperingMode = this->pTransform->getTaperingMode(), oldFftLen = this->pTransform->getFftLen();
	SC_Conversion converter(pMFCC->Hdr.sampleRate);
	SC_MatrixFunctions matFunc;

	//get the frequency-response of the used preemphasis-filter in the cepstral domain to later on remove its influence on the mfccs (see Milner&Shao)
	this->pTransform->setFftLen(usedFftLen);
	this->pTransform->setTaperingMode(usedWindow);
	preEmphasisFrequencyResponse = this->pTransform->magnitudeSpectrum(preEmphasisImpulseResponse, 2, false, false);
	MArray_1D(preEmphasisCepstrum, usedFilterbankSize, double, "SC_Synthesis.mfcc2spectrum: preEmphasisCepstrum");
	for (f = 0; f < usedFilterbankSize; f++) {
		preEmphasisCepstrum[f] = SC_Feature_FbE::evalFilter(preEmphasisFrequencyResponse, filterParam[f][0], filterParam[f][1], filterParam[f][2], pMFCC->Hdr.sampleRate, usedFftLen);
    preEmphasisCepstrum[f] = (preEmphasisCepstrum[f] <= 0.0) ? -87.3 : log(preEmphasisCepstrum[f]); //this (as below) relies on the knowledge that the same thing is done when creating mfccs, not 20*log() or something alike...
	}
	preEmphasisCepstrum = sclib::directOutput(this->pTransform->dct(preEmphasisCepstrum, usedFilterbankSize, usedFilterbankSize), preEmphasisCepstrum, true); //finish by converting to the cepstral domain; use all cepstral coefficients
	MFree_1D(preEmphasisFrequencyResponse);	

	//get the area vector of the mfcc filters in the cepstral domain to kill its effect of spectral tilt towards high frequencies later on (see Milner&Shao)
	MArray_1D(filterbankAreaCepstrum, usedFilterbankSize, double, "SC_Synthesis.mfcc2spectrum: filterbankAreaCepstrum");
	for (f = 0; f < usedFilterbankSize; f++) {
		filterbankAreaCepstrum[f] = (filterParam[f][2]-filterParam[f][0])*1.0 / 2.0; //the height of each filter is 1.0, the form is triangular, so this is the formula for the area
		filterbankAreaCepstrum[f] = (filterbankAreaCepstrum[f] <= 0.0) ? -87.3 : log(filterbankAreaCepstrum[f]); //the log transforms this "filterbank-energy" vector into log-energy, as a first attempt to produce a cepstral vector
	}
	filterbankAreaCepstrum = sclib::directOutput(this->pTransform->dct(filterbankAreaCepstrum, usedFilterbankSize, usedFilterbankSize), filterbankAreaCepstrum, true); //finish producing an area-vector in the cepstral domain; use all cepstral coefficients

	//initialize the divisor for the overlap/add of the filter output to reconstruct the smoothed spectrum
	MArray_2D(idx, (int)(usedFilterbankSize), 2, unsigned int, "SC_Synthesis.mfcc2spectrum: idx");
	divisor = matFunc.zeros(this->pTransform->getFftLen()/2 +1);
	for (f = 0; f < usedFilterbankSize; f++) {
		idx[f][0] = converter.frequency2index(filterParam[f][0], this->pTransform->getFftLen()); //precompute indices so that we can save the function calls to the converter in the loop later on
		idx[f][1] = converter.frequency2index(filterParam[f][2], this->pTransform->getFftLen());
		for (d = idx[f][0]; d <= idx[f][1]; d++) {
			filterValue = SC_Feature_FbE::evalFilter(d, filterParam[f][0], filterParam[f][1], filterParam[f][2], pMFCC->Hdr.sampleRate, this->pTransform->getFftLen());
			divisor[d] += filterValue; // * filterValue;
		}
	}

	MArray_1D(spectra, pMFCC->Row, double*, "SC_Synthesis.mfcc2spectrum: cepstrum");
	MArray_1D(cepstrum, usedFilterbankSize, double, "SC_Synthesis.mfcc2spectrum: cepstrum");

	for (t = 0; t < pMFCC->Row; t++) {
		//fill up the original cepstrum and remove the effect of the preemphasis filter and the spectral tilt introduced by the spacing of the mel filterbank
		//unwrapping dEnergy only works if the chain of MFCCs is unbroken/unaltered an the 0th row is there, otherwise result is unpredictable (and ofton leads to a signal with strange loudness behaviour: silence and one short loud burst of energy)
		cepstrum[0] = (firstMfccIsDEnergy == true) ? lastMFCC0+pMFCC->Mat[t][0] : pMFCC->Mat[t][0]; //unwrap delta-energy //formerly (in human studie): 10.0; //TODO!!!
		for (d = 1; d < usedFilterbankSize; d++) {
			if (d < (unsigned int)(pMFCC->Col)) {
				cepstrum[d] = pMFCC->Mat[t][d];
			} else {
				cepstrum[d] = 0.0; //add zeros to fill up the possibly truncated cepstrum (see Milner&Shao)
			}
			cepstrum[d] -= filterbankAreaCepstrum[d] + preEmphasisCepstrum[d]; //remember: convolution in the time domain is multiplication in the frequency domain is addition in the log-cepstral domain, thats why inverse filtering is so easy here (see Milner&Shao)
		}
		lastMFCC0 = cepstrum[0];

		//convert cepstrum to log-filterbank output
		filterOutput = this->pTransform->idct(cepstrum, usedFilterbankSize);
    
		//invert log
		for (f = 0; f < usedFilterbankSize; f++) {
			filterOutput[f] = exp(filterOutput[f]); //this relies on the knowledge that log() is used while creating mfccs, not 20*log() or something alike...
		}

		//convert the linear, mel-scaled filterbank output to a magnitude-spectrum via overlap & add
		spectrum = matFunc.zeros(this->pTransform->getFftLen()/2 +1);
		for (f = 0; f < usedFilterbankSize; f++) {
			for (d = idx[f][0]; d <= idx[f][1]; d++) {
				spectrum[d] += filterOutput[f] * SC_Feature_FbE::evalFilter(d, filterParam[f][0], filterParam[f][1], filterParam[f][2], pMFCC->Hdr.sampleRate, usedFftLen);
			}
		}
    MFree_1D(filterOutput);

		//finish OLA
		for (d = 0; d < this->pTransform->getFftLen()/2 +1; d++) {
			if (divisor[d] > 0.0) {
				spectrum[d] /= divisor[d];
			}
		}

		//include pitch information, if wanted
		addPitch(spectrum, t, pMFCC->Hdr.frameSize, pMFCC->Hdr.frameStep, &converter, pPitch);

		//save the final power spectrum
		spectra[t] = spectrum;
	} //for each frame t

	this->pTransform->setFftLen(oldFftLen);
	this->pTransform->setTaperingMode(oldTaperingMode);

	MFree_2D(idx);
	MFree_1D(filterbankAreaCepstrum);
	MFree_1D(preEmphasisCepstrum);
	MFree_1D(cepstrum);
	MFree_2D(filterParam);
	MFree_1D(divisor);

	return spectra;
}

//====================================================================================================================
// uses mfcc2powerSpec() and olaSynthesis() to convert a series of mfccs back to a signal stored under the given 
// file-name as a wav file
//====================================================================================================================
void SC_Synthesis::mfcc2wav(char *fileName, SV_Data* pMFCC, double usedPreEmphasisFactor, unsigned int usedWindow, unsigned int usedFftLen, int usedFrequencyScale, double usedFilterbankLowEnd, double usedFilterbankHighEnd, unsigned int usedFilterbankSize, bool firstMfccIsDEnergy, unsigned int intermediateFrameCount, double steepness, SV_Data *pPitch, bool keepLimitations) {
	unsigned long int sampleCount;
	double **spectra = mfcc2spectrum(pMFCC, ((keepLimitations==false)?usedPreEmphasisFactor:0.0), usedWindow, usedFftLen, usedFrequencyScale, usedFilterbankLowEnd, usedFilterbankHighEnd, usedFilterbankSize, firstMfccIsDEnergy, pPitch), *dblSamples, **blended; //only remedy pre-emphasis if not wished to hear its effects
	unsigned int spectrumCount = pMFCC->Row;
	SC_FeatureHandler handler(NULL, false);

	//do blending, if wished
	if (intermediateFrameCount > 0) {
		blended = handler.blend(spectra, pMFCC->Row, usedFftLen/2 +1, intermediateFrameCount, steepness, 1, spectrumCount);
		MFree_2Dex(spectra, pMFCC->Row);
		spectra = blended;
	}

	//resynthesize the signal with phase-reconstruction by OLA method, normalize amplitude
	sampleCount = (spectrumCount-1)*pMFCC->Hdr.frameStep + pMFCC->Hdr.frameSize;
	dblSamples = olaSynthesis(spectra, pMFCC->Hdr.frameSize, pMFCC->Hdr.frameStep, sampleCount, true, NULL, this->olaMaxIterations, this->olaErrorTarget);
	if (writeSignal(dblSamples, sampleCount, pMFCC->Hdr.sampleRate, fileName, 0.3) == false) {
		REPORT_ERROR(SVLIB_FileErr, "Couldn't write resynthesized MFCCs to given fileName");
	}	

	//free the spectral frames and sbl-samples
	MFree_1D(dblSamples);
	MFree_2Dex(spectra, (int)(spectrumCount));

	return;
}

//====================================================================================================================
//	converts all the LPC coefficients in the feature set to smoothed log magnitude spectra as seen in many textbooks
//====================================================================================================================
SV_Data* SC_Synthesis::lpc2logSpectrum(SV_Data *pLPCs) {
	SV_Data *pSpectra = NULL;
	double *spectrum;
	int t, d, spectrumSize = this->pTransform->getFftLen()/2 +1;
		
	if (pLPCs != NULL && pLPCs->Row > 0 && pLPCs->Col <= (int)(this->pTransform->getFftLen())) {
		pSpectra = new SV_Data(pLPCs->Row, spectrumSize);
		pSpectra->Hdr = pLPCs->Hdr;
		pSpectra->Hdr.ID = sclib::featureSpectrum;
	
		for (t = 0; t < pLPCs->Row; t++) {
			spectrum = lpcSpectrum(pLPCs->Mat[t], pLPCs->Col - (pLPCs->Hdr.Signature[1]==1?1:0), pLPCs->Hdr.frameSize, pLPCs->Hdr.Signature[1]==1?pLPCs->Mat[t][pLPCs->Col-1]:1.0, true, false); //use the gain term if appended to the lpc feature set
			for (d = 0; d < spectrumSize; d++) {
				pSpectra->Mat[t][d] = (float)(spectrum[d]);
			}
			MFree_1D(spectrum);
		}
	}
	
	return pSpectra;	
}

//====================================================================================================================
// converts all the LPC coefficients in the feature set to smoothed spectra; sounds not as good as possible (as it
// does e.g. in mobile phone communication) because the excitation is not considered, just the vocal tract features 
// are made audible
//====================================================================================================================
double** SC_Synthesis::lpc2spectrum(SV_Data *pLPCs, double usedPreEmphasisFactor, unsigned int usedWindow, SV_Data *pPitch) {
	int t;
	double **spectra = NULL;
	unsigned int d, oldTaperingMode = this->pTransform->getTaperingMode();
	double preEmphasisImpuleResponse[] = {1.0, -1.0*usedPreEmphasisFactor}, *preEmphasisFrequencyResponse;
	SC_Conversion converter(pLPCs->Hdr.sampleRate);

	if (pLPCs != NULL && pLPCs->Row > 0 && pLPCs->Col <= (int)(this->pTransform->getFftLen())) {
		MArray_1D(spectra, pLPCs->Row, double*, "SC_Synthesis.lpc2spectrum: spectra");

		//get the frequency-response of the used preemphasis-filter in the spectral domain to later on remove its influence on the lpcs (see Milner&Shao)
		this->pTransform->setTaperingMode(usedWindow);
		preEmphasisFrequencyResponse = this->pTransform->magnitudeSpectrum(preEmphasisImpuleResponse, 2, false, false);

		for (t = 0; t < pLPCs->Row; t++) {
			spectra[t] = lpcSpectrum(pLPCs->Mat[t], pLPCs->Col - (pLPCs->Hdr.Signature[1]==1?1:0), pLPCs->Hdr.frameSize, pLPCs->Hdr.Signature[1]==1?pLPCs->Mat[t][pLPCs->Col-1]:1.0, false, true); //use the gain term if appended to the lpc feature set

			//remove preemphasis-influence
			for (d = 0; d < this->pTransform->getFftLen()/2 +1; d++) {
				if (preEmphasisFrequencyResponse[d] > 0.0) {
					spectra[t][d] /= preEmphasisFrequencyResponse[d];
				}
			}

			//include pitch information, if wanted
			addPitch(spectra[t], t, pLPCs->Hdr.frameSize, pLPCs->Hdr.frameStep, &converter, pPitch);
		} //for t

		MFree_1D(preEmphasisFrequencyResponse);
		this->pTransform->setTaperingMode(oldTaperingMode);

	} //if something to do
	
	return spectra;	
}

//====================================================================================================================
// uses lpc2powerSpec() and olaSynthesis() to convert a series of lpcs back to a signal stored under the given 
// file-name as a wav file
//====================================================================================================================
void SC_Synthesis::lpc2wav(char *fileName, SV_Data* pLPC, double usedPreEmphasisFactor, unsigned int usedWindow, unsigned int intermediateFrameCount, double steepness, SV_Data *pPitch, bool keepLimitations) {
	unsigned long int sampleCount;
	double **spectra, **blended, *dblSamples;
	unsigned int spectrumCount = pLPC->Row, oldFftLen = this->pTransform->getFftLen();
	SC_FeatureHandler handler(NULL, false);

	//convert lpc to spectra, adjust fft-length, if needed
	while (this->pTransform->getFftLen() < pLPC->Hdr.frameSize) {
		this->pTransform->setFftLen(this->pTransform->getFftLen()*2);
	}
	spectra = lpc2spectrum(pLPC, ((keepLimitations==false)?usedPreEmphasisFactor:0.0), usedWindow, pPitch); //only remedy pre-emphasis if not wished to hear ist effects

	//do blending, if wished
	if (intermediateFrameCount > 0) {
		blended = handler.blend(spectra, pLPC->Row, this->pTransform->getFftLen()/2 +1, intermediateFrameCount, steepness, 1, spectrumCount);
		MFree_2Dex(spectra, pLPC->Row);
		spectra = blended;
	}

	//resynthesize the signal with phase-reconstruction by OLA method, normalize amplitude
	sampleCount = (spectrumCount-1)*pLPC->Hdr.frameStep + pLPC->Hdr.frameSize;
	dblSamples = olaSynthesis(spectra, pLPC->Hdr.frameSize, pLPC->Hdr.frameStep, sampleCount, true, NULL, this->olaMaxIterations, this->olaErrorTarget);
	if (writeSignal(dblSamples, sampleCount, pLPC->Hdr.sampleRate, fileName, 0.3) == false) {
		REPORT_ERROR(SVLIB_FileErr, "Couldn't write resynthesized LPCs to given fileName");
	}	
	
	//free the spectral frames and dbl-samples
	MFree_1D(dblSamples);
	MFree_2Dex(spectra, (int)(spectrumCount));
	
	this->pTransform->setFftLen(oldFftLen);

	return;
}

//====================================================================================================================
// takes framed samples and reconnects them to a full signal, saving it to the given filename as a wav file
//====================================================================================================================
void SC_Synthesis::samples2wav(char *fileName, SV_Data *pSamples) {
	unsigned long int sampleCount, s, p, i;
	double *divisor, *dblSamples;
	SC_FeatureHandler handler(NULL, false);
	SC_MatrixFunctions matFunc;

	//resynthesize the signal by OLA method, determine max. amplitude
	sampleCount = (pSamples->Row-1)*pSamples->Hdr.frameStep + pSamples->Hdr.frameSize;
	divisor = matFunc.zeros(sampleCount);
	dblSamples = matFunc.zeros(sampleCount);
	for (p = 0; p < (unsigned long int)(pSamples->Row); p++) {
		for (i = 0; i < (unsigned long int)(pSamples->Col); i++) {
			divisor[p*pSamples->Hdr.frameStep + i] += 1.0;
			dblSamples[p*pSamples->Hdr.frameStep + i] += pSamples->Mat[p][i];
		}
	}
	for (s = 0; s < sampleCount; s++) {
		if (divisor[s] > 0.0) { //avoid div. by zero; it doesn't do any harm because when the divisor is zero, this means, one of the coefficients was zero, too, so the numerator will also be zero!
			dblSamples[s] /= divisor[s];
		}
	}
	MFree_1D(divisor);

	//write the signal to disc, free occupied memory
	if (writeSignal(dblSamples, sampleCount, pSamples->Hdr.sampleRate, fileName, 0.3) == false) {
		REPORT_ERROR(SVLIB_FileErr, "Couldn't write reunited samples to given fileName");
	}	
	MFree_1D(dblSamples);
	
	return;
}

//====================================================================================================================
// deduces from the featureType which specialized method is to be used and takes the needed parameters from the 
// tweakable parameters object; returns false if something went wrong, e.g. feature type unsupported
// if keepLimitations==true, some things (pre-emphasis at the moment) are not conversed in order to hear their effect
//====================================================================================================================
bool SC_Synthesis::feature2wav(char *fileName, SV_Data *pFeature, SV_Data *pPitch, bool keepLimitations) {
	bool result = true;

	if (pFeature->Hdr.ID == sclib::featureSamples) {
		samples2wav(fileName, pFeature);
	} else if (pFeature->Hdr.ID == sclib::featureLPC) {
		this->setTaperingMode(this->pTweak->featureLpc.window);
		lpc2wav(fileName, pFeature, this->pTweak->featureLpc.preEmphasizeFactor, this->pTweak->featureLpc.window, 0, 1.0, pPitch, keepLimitations);
	} else if (pFeature->Hdr.ID == sclib::featureSpectrum) {
		this->setFftLength(this->pTweak->featureSpectrum.FFTsize);
		this->setTaperingMode(this->pTweak->featureSpectrum.window);
		spectrum2wav(fileName, pFeature, false, this->pTweak->featureSpectrum.logarithmize, true, keepLimitations);
	} else if (pFeature->Hdr.ID == sclib::featureMFCC) {
		//this->setFftLength(this->pTweak->featureMfcc.fftSize); //is done inside mfcc2spectrum(), this is called inside mfcc2wav()
		//this->setTaperingMode(this->pTweak->featureMfcc.window);
		mfcc2wav(fileName, pFeature, this->pTweak->featureMfcc.preEmphasizeFactor, this->pTweak->featureMfcc.window, this->pTweak->featureMfcc.fftSize, this->pTweak->featureMfcc.sclib_frequencyScale, this->pTweak->featureMfcc.sclib_minFilterBankFrequency, this->pTweak->featureMfcc.sclib_maxFilterBankFrequency, this->pTweak->featureMfcc.filterBankSize, this->pTweak->featureMfcc.dEnergy, 0, 1.0, pPitch, keepLimitations);
	} else if (pFeature->Hdr.ID == sclib::featurePitch) {
		if (keepLimitations == false) {
			bool harmonics[] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
			generateHumFromPitch(fileName, pFeature, harmonics, 64, true);
		} else {
			generateHumFromPitch(fileName, pFeature, NULL, 0, false);
		}
	}	else {
		result = false;
	}
	
	return result;
}
