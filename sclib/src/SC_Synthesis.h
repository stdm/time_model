/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class that unites and encapsulates many algorithms for sound    */
/*        synthesis                                                       */
/*																																				*/
/*    Author  : Thilo Stadelmann															            */
/*    Date    : 12.06.2008																								*/
/**************************************************************************/

#ifndef __SC_Synthesis_H__
#define __SC_Synthesis_H__

#include "SC_Aux.h"
#include "SC_Transform.h"
#include "SC_Conversion.h"
#include "SC_TweakableParameters.h"
#include <SV_Feature.h>

class SC_Synthesis {

	private :

  protected :
	
		//=====================================================================================================================
		// creates and returns the smoothed log spectrum from the lpc coefficients of a signal-window; this is a template class 
		// to allow for float (as originating from an SV_Data feature set) and double (used in all other occasions) input.
		// recipe take from "Fundamentals of Speech Recognition", Rabiner/Juang, exercise/solution 3.4(c)
		//=====================================================================================================================
		template<class T> double* lpcSpectrum(T *lpcCoefficients, unsigned int coefficientCount, unsigned int frameSize, double gain = 1.0, bool logarithmize = true, bool givePower = false) {
			double **coeffSpec,*spectrum, *sig, real, imag, sumOfSquares;
			unsigned int i, oldTaperingMode = this->pTransform->getTaperingMode();
				
			MArray_1D(sig, coefficientCount+1, double, "SC_Synthesis.lpcSpectrum: sig");
			sig[0] = 1.0; //add 0th coefficient that is omitted in lpc generation because it always equals "1"
			for (i = 0; i < coefficientCount; i++) {
				sig[i+1] = (double)(lpcCoefficients[i]); //append coefficients
			}
			
			this->pTransform->setTaperingMode(sclib::wndRectangle); //switch windowing off for now
			coeffSpec = this->pTransform->fft(sig, coefficientCount+1);
			MFree_1D(sig);	
			this->pTransform->setTaperingMode(oldTaperingMode); //switch back to original state

			MArray_1D(spectrum, this->pTransform->getFftLen()/2 +1, double, "SC_Synthesis.lpcSpectrum: spectrum");
			for (i = 0; i < this->pTransform->getFftLen()/2 +1; i++) {
				//to get the frequency response (i.e. the complex spectrum) of the lpc-filter represented by its coefficients, 
				//we need to divide the spectrum of an impulse (g,0,0,0,...) (g is the gain here) by the spectrum of the signal just generated above; 
				//the impulse has a spectrum with its real part all equal to "g" and its imaginary part all equal to "0". hence:
				//divide complex number (g+i*0.0) by complex number coeffSpec[][i]; see Stoecker, "Taschenbuch Mathematischer Formeln und Moderner 
				//Verfahren", 3rd revision, p.582
				sumOfSquares = coeffSpec[i][0]*coeffSpec[i][0] + coeffSpec[i][1]*coeffSpec[i][1];
				real = (gain * coeffSpec[i][0]) / sumOfSquares;
				imag = (-1.0 * gain * coeffSpec[i][1]) / sumOfSquares;
				
				spectrum[i] = (givePower == false) ? sqrt(real*real + imag*imag) : (real*real + imag*imag)/(double)(frameSize);
				if (logarithmize == true) {
					spectrum[i] = log(spectrum[i]);
				}
			}
			MFree_2D(coeffSpec);
			
			return spectrum;
		}	

		//=====================================================================================================================
		// (recursive) filter (i.e. IIR) the given signal with the specified filter nominator (b) and denominator (a) 
		// coefficients as implemented in matlab's filter()-method:
		//   y(n) = (b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
		//                     - a(2)*y(n-1) - ... - a(na+1)*y(n-na)) / a(1)
		// a(1) is meant to equal 1;
		// if overwrite==true, the original sample values in the signal are changed and the poiter to the original signal is 
		// returned
		//=====================================================================================================================
		template<class T> double* recursiveFilter(double *signal, unsigned int signalLen, T *b, unsigned int bLen, T *a, unsigned int aLen, double &maxValue, bool overwrite = false, bool invert = false) {
			T *nominator, *denominator;
			double *newSig = NULL;
			unsigned int n, f, nomLen, denLen;
			
			if (invert == false) {
				nominator = b;
				nomLen = bLen;
				denominator = a;
				denLen = aLen;
			} else {
				nominator = a;
				nomLen = aLen;
				denominator = b;
				denLen = bLen;
			}
			
			maxValue = 0.0;
			if (signalLen > 0) {
				MArray_1D(newSig, signalLen, double, "SC_Synthesis.recursiveFilter: newSig");
				for (n = 0; n < signalLen; n++) {
					newSig[n] = 0.0;
					for (f = 0; f < sclib::min(n+1, nomLen); f++) {
						newSig[n] += (double)(nominator[f]) * signal[n-f];
					}
					for (f = 1; f < sclib::min(n+1, denLen); f++) {
						newSig[n] -= (double)(denominator[f]) * newSig[n-f];
					}
					newSig[n] /= (double)(denominator[0]);
					
					if (fabs(newSig[n]) > maxValue) {
						maxValue = newSig[n];
					}
				}
				
				if (overwrite == true) {
					for (n = 0; n < signalLen; n++) {
						signal[n] = newSig[n];
					}
					MFree_1D(newSig); //the buffer was needed anyway...
					newSig = signal;
				}
			}
			
			return newSig;
		}

		//=====================================================================================================================
		// (convolutive) filter (i.e. FIR) the given signal with the specified filter kernel (i.e. impulse response) as 
		// described  in Steven W. Smith, "Digital Signal Processing - A Practical Guide for Engineers and Scientists", chapter 
		// 6, Table 6-2 (AKA dspguide):
		//   y[i] = SUM_{j=0..M-1} h[j]*x[i-j]
		// y is the filtered signal, x is the input signal, and h is the filter kernel; if x has length N and h has length M, 
		// y has length N+M-2. This "naive" convolution has a runtime exponential in the filter kernel length and should be 
		// replaced with faster fft convolution (multiply the kernel's frequency response with the spectrum of the signal, 
		// synthesized via OLA) if it is longer than, say, 60 (see Smith, chapter 18).
		//=====================================================================================================================
		template<class T> double* convolutiveFilter(double *signal, unsigned int signalLen, T *h, unsigned int hLen) {
			double *newSig = NULL;
			unsigned int i, f, newLen = signalLen+hLen-2;
						
			if (signalLen > 0) {
				MArray_1D(newSig, newLen, double, "SC_Synthesis.convolutiveFilter: newSig");
				for (i = 0; i < newLen; i++) {
					newSig[i] = 0.0;
					for (f = 0; f < hLen; f++) {
						if (i-f >= 0 && i-f < signalLen) { 
							newSig[i] += (double)(h[f]) * signal[i-f];
						}
					}
				}
			}
			
			return newSig;
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
		double* fftFilter(double *signalMagnitude, double *signalPhase, unsigned int originalSignalLen, double *filterMagnitude, double *filterPhase, unsigned int originalFilterLen, unsigned int usedFftLen, unsigned int &resultLen, bool spectraArePolar = true);
		
		//=====================================================================================================================
		// interface to the fftFilter above that uses the signal and impulse response as inputs and does the necessary 
		// conversions; all comments above also apply here.
		//=====================================================================================================================
		template<class T> double* fftFilter(double *signal, unsigned int signalLen, T *h, unsigned int hLen, unsigned int &resultLen) {
			double *newSignal, *signalMag, *signalPhase, *filterMag, *filterPhase;
			unsigned long int oldFftLen = this->pTransform->getFftLen();
			unsigned int oldTaperingMode = this->pTransform->getTaperingMode();
			
			while (this->pTransform->getFftLen() < signalLen+hLen-2) { //avoid cricular convolution
				this->pTransform->setFftLen(this->pTransform->getFftLen()*2);
			}
			this->pTransform->setTaperingMode(sclib::wndRectangle); //switch off automatic windowing for now
			
			filterMag = this->pTransform->fft(h, hLen, filterPhase); //in fact, the ral and imaginary parts of the fourier transform are returned (i.e. rectangular form) instead of magnitude and phase as the variable names suggest (polar form)
			signalMag = this->pTransform->fft(signal, signalLen, signalPhase);
			
			newSignal = fftFilter(signalMag, signalPhase, signalLen, filterMag, filterPhase, hLen, this->pTransform->getFftLen(), resultLen, false);
						
			MFree_1D(filterMag);
			MFree_1D(filterPhase);
			MFree_1D(signalMag);
			MFree_1D(signalPhase);
						
			this->pTransform->setFftLen(oldFftLen); //restore previous state to avoid side effects
			this->pTransform->setTaperingMode(oldTaperingMode);
			
			return newSignal;
		}
		
		//=====================================================================================================================
		// returns the filter coefficients for a filter with the desired frequency/phase response; see Steven W. Smith,
		// "Digital Signal Processing - A Practical Guide for Engineers and Scientists", pp. 297-300. 
		// The last 3 parameters are for the reconstruction of the fft from magnitude/phase and should give the 
		// parameters used to construct the desired* parameters. The result is a filter kernel (or impulse response) that can 
		// directly be used with the convolutiveFilter() method.
		//=====================================================================================================================
		double* constructCustomFilter(double *desiredMagnitudeResponse, double *desiredPhaseResponse, int responseLength, int filterKernelLength, bool isFullLength, bool isDegrees, bool isLogarithmized);

		//=====================================================================================================================
		// normalize the given signal to relativeMaxAmplitude% max. amplitude and write it to disc
		//=====================================================================================================================
		bool writeSignal(double *dblSamples, unsigned long int sampleCount, unsigned int sampleRate, const char *fileName, double relativeMaxAmplitude = 0.95);

		//=====================================================================================================================
		// receives a spectral frame and adds pitch to it; this->pTransform->getFftLen() must return the correct spectrum size
		// used to create this spectrum
		//=====================================================================================================================
		void addPitch(double *spectrum, unsigned long int spectrumNr, int spectrumFrameSize, int spectrumFrameStep, SC_Conversion *pConverter, SV_Data *pPitch);

		//=====================================================================================================================
		// members
		//=====================================================================================================================
		unsigned int olaMaxIterations;
		double olaErrorTarget;
		bool verbose;
		SC_Transform *pTransform;
		SC_TweakableParameters *pTweak;
  		
	public :

    //====================================================================================================================
    // constructor/destructor
    //====================================================================================================================
		SC_Synthesis(SC_TweakableParameters *pTweak, unsigned int fftLength = 512, double taperingLength = 1.0, unsigned int taperingMode = sclib::wndHamming, unsigned int olaMaxIterations = 100, double olaErrorTarget = 4.0,  bool verbose = true);
		virtual ~SC_Synthesis(void);

    //====================================================================================================================
    // getter/setter
    //====================================================================================================================
		void setFftLength(unsigned int newLength) {this->pTransform->setFftLen(newLength); return;}
		void setTaperingLength(double newLength) {this->pTransform->setTaperingLength(newLength); return;}
		void setTaperingMode(unsigned int newMode) {this->pTransform->setTaperingMode(newMode); return;}
		void setOlaMaxIterations(unsigned int newCount) {this->olaMaxIterations = newCount; return;}
		void setOlaErrorTarget(double newTarget) {this->olaErrorTarget = newTarget; return;}

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
		double* olaSynthesis(double **powerSpectra, unsigned int frameSize, unsigned int frameStep, unsigned long int originalSignalLength, bool powerIsMagnitude = false, double **originalPhase = NULL, unsigned int maxIterations = 200, double errorAim = 4.0);
				
		//====================================================================================================================
		// creates a .wav-file of the given name containing a "hummed" melody following the pitch given in the first column of
		// pPitch; if addHarminics!=NULL, the (i+1)th harmonics (i=0 being the fundamental) is added to the hum if the array
		// has the value 'true' at index-position i.
		//====================================================================================================================
		void generateHumFromPitch(char *fileName, SV_Data *pPitch, bool *addHarmonics = NULL, unsigned int numberOfHarmonics = 0, bool addVoicing = false);

    //====================================================================================================================
		// take the samples of the internal buffer, smooth spectrum by replacing the amplitude for each frequency by the mean 
		// in its critical band and resynthesize the result as a wav-file
    //====================================================================================================================
		void generateVoiceFromCB(char *fileName, short *samples, unsigned long int sampleCount, unsigned int sampleRate, unsigned int frameSize, unsigned int frameStep);

		//====================================================================================================================
		// uses olaSynthesis() to convert a series of (power) spectra back to a signal stored under the given file-name as a 
		// wav file; if smooth==true, the spectra are smoothed before concatenating (i.e. 20% highest frequencies in the 
		// spectral envelope are removed); if deLog==true, the spectra are assumed to be log-spectra and exponentialized 
		// before resynthesis; if magnitudeIsPower==true, the spectra are assumed to be power-spectra instead of magnitude 
		// spectra
		//====================================================================================================================
		void spectrum2wav(char *fileName, SV_Data* pSpectrum, bool smooth = false, bool deLog = false, bool magnitudeIsPower = false, bool keepLimitations = false);

		//====================================================================================================================
		// gets a set of mfcc vectors (and maybe corresponding pitch) and converts them back to spectra; the result-set is not 
		// a matrix but an array pointing to #pMFCC->Row single spectrum arrays of length usedFftLen/2; the given parameters 
		// mainly tell how the mfccs where created in order to correctly invert the process; references for ideas used in this 
		// method are: Milner, Shao, "Clean Speech Reconstruction from MFCC Vectors and Fundamental Frequency using an 
		// Integrated Front-End", 2005; Milner, Shao, "Speech Reconstruction from Mel-Frequency Cepstral Coefficients using a 
		// Source-Filter Model", 2002
		//====================================================================================================================
		double** mfcc2spectrum(SV_Data* pMFCC, double usedPreEmphasisFactor, unsigned int usedWindow, unsigned int usedFftLen, int usedFrequencyScale, double usedFilterbankLowEnd, double usedFilterbankHighEnd, unsigned int usedFilterbankSize, bool firstMfccIsDEnergy, SV_Data *pPitch = NULL);

		//====================================================================================================================
		// uses mfcc2powerSpec() and olaSynthesis() to convert a series of mfccs back to a signal stored under the given 
		// file-name as a wav file
		//====================================================================================================================
		void mfcc2wav(char *fileName, SV_Data* pMFCC, double usedPreEmphasisFactor, unsigned int usedWindow, unsigned int usedFftLen, int usedFrequencyScale, double usedFilterbankLowEnd, double usedFilterbankHighEnd, unsigned int usedFilterbankSize, bool firstMfccIsDEnergy, unsigned int intermediateFrameCount = 0, double steepness = 1.0, SV_Data *pPitch = NULL, bool keepLimitations = false);
		
		//====================================================================================================================
		// converts all the LPC coefficients in the feature set to smoothed log magnitude spectra as seen in many textbooks
		//====================================================================================================================
		SV_Data* lpc2logSpectrum(SV_Data *pLPCs);

		//====================================================================================================================
		// converts all the LPC coefficients in the feature set to smoothed spectra; sounds not as good as possible (as it
		// does e.g. in mobile phone communication) because the excitation is not considered, just the vocal tract features 
		// are made audible
		//====================================================================================================================
		double** lpc2spectrum(SV_Data *pLPCs, double usedPreEmphasisFactor, unsigned int usedWindow, SV_Data *pPitch = NULL);

		//====================================================================================================================
		// uses lpc2powerSpec() and olaSynthesis() to convert a series of lpcs back to a signal stored under the given 
		// file-name as a wav file
		//====================================================================================================================
		void lpc2wav(char *fileName, SV_Data* pLPC, double usedPreEmphasisFactor, unsigned int usedWindow, unsigned int intermediateFrameCount = 0, double steepness = 1.0, SV_Data *pPitch = NULL, bool keepLimitations = false);

		//====================================================================================================================
		// takes framed samples and reconnects them to a full signal, saving it to the given filename as a wav file
		//====================================================================================================================
		void samples2wav(char *fileName, SV_Data *pSamples);

		//====================================================================================================================
		// deduces from the featureType which specialized method is to be used and takes the needed parameters from the 
		// tweakable parameters object; returns false if something went wrong, e.g. feature type unsupported
		// if keepLimitations==true, some things (pre-emphasis at the moment) are not conversed in order to hear their effect
		//====================================================================================================================
		bool feature2wav(char *fileName, SV_Data *pFeature, SV_Data *pPitch = NULL, bool keepLimitations = false);

		//====================================================================================================================
		// gets a signal (samples) along with extracted LPC and formant features on it; by inverse filtering, the effect of 
		// the first 'removeUntil' formants is removed from the spectrum; the result is a new signal that mainly consists of 
		// the voice rather than the speech (so it is hoped)
		// the bandwidth of the found formants can be stretched with setting bandwidthEnhancementFactor>1.0
		//====================================================================================================================
		template<class T> void removeFormantInfluence(T *samples, unsigned long int sampleCount, unsigned long int sampleRate, SV_Data *pFormants, unsigned int removeUntil = 2, double bandwidthEnhancementFactor = 1.0) {
			int t, filterLen = 512, i, frameStart, frameLen, lastFrameEnd = -1;
			double *filterMagnitude, *signalPhase, *signalMagnitude, lastPercentage = 0.0, *newSamples, *filteredFrame, *frame;
			SC_Conversion converter(sampleRate);
			unsigned int frameStep = converter.ms2sample(converter.sample2ms(pFormants->Hdr.frameStep, pFormants->Hdr.sampleRate), sclib::alignmentStart, sampleRate); //formant tracking maybe involves downsampling, be prepared for that! 
			unsigned int frameSize = converter.ms2sample(converter.sample2ms(pFormants->Hdr.frameSize, pFormants->Hdr.sampleRate), sclib::alignmentStart, sampleRate);	
			unsigned long int oldFftLen = this->pTransform->getFftLen();
			unsigned int oldTaperingMode = this->pTransform->getTaperingMode();
			unsigned int filteredLen, p, fMin, fMax, lastF = 0;

			if (this->verbose == true) {
				printf("\nFormant influence removal: ");
			}
			this->pTransform->setTaperingMode(sclib::wndRectangle);
			while (this->pTransform->getFftLen() < frameSize+filterLen-2) { //the size of a filtered signal is sigSize+filterSize-2 (it is indeed, and it is necessary to avoid circular convolution if the fft is used to implement a convolutive filter), see dspguide book pp.183-184 and p.314
				this->pTransform->setFftLen(this->pTransform->getFftLen()*2);
			}
			MArray_1D(newSamples, sampleCount+filterLen-2, double, "SC_Synthesis.removeFormantInfluence: newSamples");
			MArray_1D(frame, frameSize, double, "SC_Synthesis.removeFormantInfluence: frame");
			MArray_1D(filterMagnitude, this->pTransform->getFftLen()/2+1, double, "SC_Synthesis.removeFormantInfluence: filterMagnitude");
				
			for (t = 0; t < pFormants->Row; t++) {
				//filter construction
				//filterMagnitude = lpcSpectrum(pLPCs->Mat[t], pLPCs->Col-(pLPCs->Hdr.Signature[1]==1?1:0), pLPCs->Hdr.Signature[1]==1?pLPCs->Mat[t][pLPCs->Col-1]:1.0, false, false); //use the gain term if appended to the lpc feature set
				lastF = 0;
				for (p = 0; p < removeUntil; p++) {
					fMin = converter.frequency2index(sclib::max(pFormants->Mat[t][p] - pFormants->Mat[t][pFormants->Col/2 +p]*bandwidthEnhancementFactor, 0.0), this->pTransform->getFftLen());
					fMax = converter.frequency2index(sclib::min(pFormants->Mat[t][p] + pFormants->Mat[t][pFormants->Col/2 +p]*bandwidthEnhancementFactor, sampleRate/2.0), this->pTransform->getFftLen());
					for (i = (int)(lastF); i < (int)(fMin); i++) {
						filterMagnitude[i] = 1.0;
					}
					for (i = sclib::max(fMin, lastF); i <= (int)(fMax); i++) {
						filterMagnitude[i] = 0.0; //(filterMagnitude[i] > 0) ? 1.0/filterMagnitude[i] : 0.0;
					}
					lastF = sclib::max(lastF, fMax+1);
				}
				for (i = lastF; i < (int)(this->pTransform->getFftLen()/2+1); i++) {
					filterMagnitude[i] = 1.0;
				}

				//filtering
				frameStart = t*frameStep;
				frameLen = sclib::min(frameSize, sampleCount-frameStart);
				for (i = 0; i < frameLen; i++) {
					frame[i] = (double)(samples[frameStart+i]);
				}
				signalMagnitude = this->pTransform->magnitudeSpectrum(frame, frameLen, signalPhase, false, false);
				filteredFrame = fftFilter(signalMagnitude, signalPhase, frameLen, filterMagnitude, NULL, filterLen, this->pTransform->getFftLen(), filteredLen, true);
				MFree_1D(signalMagnitude);
				MFree_1D(signalPhase);
				
				//OLA synthesis of the overlapping frames (they overlap at least because of the filtering, which enlarges the resulting frame, see dspguide, pp. 314)
				for (i = frameStart; i < (int)(frameStart+filteredLen); i++) {
					newSamples[i] = (i > lastFrameEnd) ? filteredFrame[i-frameStart] : newSamples[i]+filteredFrame[i-frameStart]; //only add if there is overlap
				}
				lastFrameEnd = frameStart + filteredLen - 1;
				MFree_1D(filteredFrame);
				
				if (this->verbose == true) {
					lastPercentage = sclib::printPercentage(pFormants->Row, t, lastPercentage, 5.0, t==0);
				}
			}
			MFree_1D(frame);
			MFree_1D(filterMagnitude);
			
			//copy back to the real buffer; omit the last filterLen-2 samples generated by the filtering process
			for (t = 0; t < (int)(sclib::min(sampleCount, lastFrameEnd)); t++) {
				samples[t] = (T)(newSamples[t]);
			}
			MFree_1D(newSamples);
			
			if (this->verbose == true) {
				sclib::printPercentage(pFormants->Row, pFormants->Row, 0.0, 5.0, false);
			}	
			this->pTransform->setFftLen(oldFftLen);
			this->pTransform->setTaperingMode(oldTaperingMode);
			
			return;
		}
};

#endif
