/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for time-domain <-> frequency-domain transformations      */
/*																																				*/
/*		Background:																													*/
/*			- The transform-(fft)-algorithm is due to Prof.Dr.Oskar Hofmann,	*/
/*				FH Giessen-Friedberg, WS2003/2004 (Biosignalverarbeitung) 			*/
/*																																				*/
/*    Author  : Thilo Stadelmann															            */
/*    Date    : 26.10.2004																								*/
/**************************************************************************/

#ifndef __SC_Transform_H__
#define __SC_Transform_H__

#include <stdlib.h>
#include "SC_Aux.h"

class SC_Transform {

	private :

  protected :

		//====================================================================================================================
    // fills the signal in x with zeros till it reaches the length of this->fftLen
    //====================================================================================================================
    void zeroPadding(double *x, unsigned long int zeroStart, unsigned long int len);

    //====================================================================================================================
    // some auxiliary functions needed by transform()
    //====================================================================================================================
    long invert(int k, int m);
    void sortCoeff(double* xr, double* xi, unsigned int N, unsigned int M);
    void copySignal(double* to, unsigned long int len, const double* from = NULL);

    //====================================================================================================================
    // interface function that can either call the below written fastFourierTransform() (see there for parameter 
		// explanations or call the svlib version after conversion of parameters, if possible
    //====================================================================================================================
		bool ffTransform(const double* in_r, const double* in_i, double*& out_r, double*& out_i, unsigned long int originalLen, unsigned long int fftLen, unsigned short mode, unsigned int taperingMode, double taperingLength);

    //====================================================================================================================
    // The FFT algorithm itself with the follwing variants:
    //   mode = sclib::modeFourierCoefficients: fourier-coefficients
    //   mode = sclib::modeFourierTransform: fourier-transform
    //   mode = sclib::modeInverseFourierTransform: inverse fourier-transform
    // Input (real/imaginary) comes in in_* (and isn't altered)
    // Output is returned in out_* (the space is allocated by this function)
    // If computation fails, the return value is false, otherwise it's true
    //====================================================================================================================
    bool fastFourierTransform(const double* in_r, const double* in_i, double*& out_r, double*& out_i, unsigned long int originalLen, unsigned long int fftLen, unsigned short mode, unsigned int taperingMode, double taperingLength);
   
    //====================================================================================================================
    // low-pass filter for the dwt; the order of the ised filterbank-coefficients can be choosen; the result has half the
		// lenth of the original signal
    //====================================================================================================================
    double* dwtLowpass(double* signal,unsigned long int len, unsigned int daubechiesOrder = 2);

		//====================================================================================================================
		// high-pass filter for the dwt; the order of the ised filterbank-coefficients can be choosen; the result has half the
		// lenth of the original signal
    //====================================================================================================================
	  double* dwtHighpass(double* signal,unsigned long int len, unsigned int daubechiesOrder = 2);

		//====================================================================================================================
		// lowpass FIR filter coefficients (highpass: reverse the tap order and multiply by the sequence 1,-1, 1,-1,..) for 
		// Daubechies wavelets 1-38 (taken from: http://www.musicdsp.org/showone.php?id=20)
		// the n*2 coefficients for "Daubechies n" start at index 2*sclib::gaussSum(n-1), e.g. Daub4 resides in 
		// daubechiesCoefficients[12]..daubechiesCoefficients[19] 
		//====================================================================================================================
		static const double daubechiesCoefficients[1556];

		//====================================================================================================================
    // members
    //====================================================================================================================
    unsigned int taperingMode;
    unsigned long int fftLen;
    double taperingLength; //how many percent of the frame shall be affected by the window function (between 0 and 1)?
		char fftOrigin; //can be set to sclib::modeSClib or sclib::modeSVlib to choose an fft implementation
		
	public :

    //====================================================================================================================
    // constructor/destructor
    //====================================================================================================================
	  SC_Transform(unsigned long int fftLen = 512, unsigned int taperingMode = sclib::wndHamming, double taperingLength = 1.0, char fftOrigin = sclib::modeSClib);
		virtual ~SC_Transform(void);

    //====================================================================================================================
    // getter/setter
    //====================================================================================================================
		inline unsigned int getTaperingMode(void) {return this->taperingMode;}
		inline double getTaperingLength(void) {return this->taperingLength;}
		inline unsigned long int getFftLen(void) {return this->fftLen;}
		inline char getFftOrigin(void) {return this->fftOrigin;}
		inline void setTaperingMode(unsigned int newMode) {this->taperingMode = newMode; return;}
		inline void setTaperingLength(double newLength) {this->taperingLength = newLength; return;}
		inline void setFftLen(unsigned long int newLength) {this->fftLen = newLength; return;}
		inline void setFftOrigin(char newOrigin) {this->fftOrigin = newOrigin; return;}
		
    //====================================================================================================================
    // retruns the coefficients of the specified window of length len
    //====================================================================================================================
		double* window(unsigned int len, unsigned int mode);

    //====================================================================================================================
    // some window-functions, which return the coefficients of the desired window-type in a newly created array of lenght
    // len
    //====================================================================================================================
    double* hamming(unsigned int len);
    double* hanning(unsigned int len);
    double* bartlett(unsigned int len);
		double* rectangular(unsigned int len);

    //====================================================================================================================
		// Applies a window-function to the signal in x according to windowType
		// Tapering is normally done by the below transformation functions, so you don't need to call it yourself.
		// It's public just that you can taper a signal which you don't want to transform.
		//
		// taperingLength (between 0 and 1) tells how many percent of the signal shall be windowed (100% -> the whole signal
		// is subject to the window-function, 20% -> only the first and last 10% are smoothly drawn to zero)
		// If smoothSpectrum is true, only the second half of the window is used to smoothly draw the signal to zero
		// on a length of N samples
		// if invert is true, a previous tapering is undone by multiplying the signal with the inverted window function (1/w)
    //====================================================================================================================
    void tapering(double *x, unsigned int len, unsigned int windowType = sclib::wndHamming, double taperingLength = 0.2, bool smoothSpectrum = false, bool invert = false);

    //====================================================================================================================
    // some functions which convert the fourier-transform (xr/xi) into other quantities (magnitude/power/phase spectrum)
    // and back; the original data is overwritten with the result
    //====================================================================================================================
    void  ft2magnitudePhase(double*& xr, double*& xi, bool degrees = false, bool fullLength = true, bool logarithmize = false);
    void  ft2powerPhase(double*& xr, double*& xi, bool degrees = false, unsigned long int originalLength = 0, bool fullLength = true, bool logarithmize = false);
    void  ft2power(double*& xr, double*& xi, unsigned long int originalLength = 0, bool fullLength = true, bool logarithmize = false);
    void  magnitudePhase2ft(double*& magnitude, double*& phase, bool degrees = false, bool fullLength = true, bool logarithmized = false);
    void  powerPhase2ft(double*& power, double*& phase, bool degrees = false, unsigned long int originalLength = 0, bool fullLength = true, bool logarithmized = false);
    void  ft2magnitude(double*& xr, double*& xi, bool fullLength = true, bool logarithmize = false);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's fourier-coefficients in an 2-dimensional array (dim0:ai, dim1:bi)
    //====================================================================================================================
    double** fourierCoeff(double* signal, unsigned long int len);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's fourier-transform in an 2-dimensional array (col0:re, col1:im)
    //====================================================================================================================
    double** fft(double* signal, unsigned long int len);
		double* fft(double *signal, unsigned int len, double* &xi);
		
    //====================================================================================================================
    // Gets the fourier-transformed in an 2-dimensional array (dim0:re, dim1:im) of length len=fftLen/2, 
    // returns the signal
    //====================================================================================================================
    double* ifft(double** fourierTransform, unsigned long int len);
    double* ifft(double* ft_r, double* ft_i, unsigned long int len);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's power-spectrum of length fftLen/2 +1 or fftLen according to the 
    // parameter fullLength
    //====================================================================================================================
    double* powerSpectrum(double* signal, unsigned long int len, bool fullLength = false, bool logarithmize = false);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's power-spectrum and the corresponding phase-angle
    // space for the return-values is allocated by this function
    //====================================================================================================================
    double* powerSpectrum(double* signal, unsigned long int len, double* &phaseAngle, bool fullLength = false, bool logarithmize = false);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's magnitude-spectrum of length fftLen/2 +1
    //====================================================================================================================
    double* magnitudeSpectrum(double* signal, unsigned long int len, bool fullLength = false, bool logarithmize = false);

		//====================================================================================================================
		// Gets the signal of length "len", returns it's magnitude-spectrum of length fftLen/2 +1 and the corresponding phase-
		// angle of length fftLen; space for the return-values is allocated by this function
		//====================================================================================================================
		double* magnitudeSpectrum(double* signal, unsigned long int len, double* &phaseAngle, bool fullLength = false, bool logarithmize = false);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's smoothed power-spectrum
    // M is the count of coefficients to hold after fft of the raw spectrum and before ifft back to the smoothed spectrum
    //====================================================================================================================
    double* smoothedPowerSpectrum(double* signal, unsigned long int len, unsigned int M = 0);

		//====================================================================================================================
		// Gets the signal of length "len", returns it's smoothed magnitude-spectrum
		// M is the count of coefficients to hold after fft of the raw spectrum and before ifft back to the smoothed spectrum
		//====================================================================================================================
		double* smoothedMagnitudeSpectrum(double* signal, unsigned long int len, unsigned int M = 0);

    //====================================================================================================================
		// Gets the (power) spectrum of either full length (fftSize) or normal (fftSize/2 +1), returns it's smoothed version
		// M is the count of coefficients to hold after fft of the raw spectrum and before ifft back to the smoothed spectrum
		//====================================================================================================================
    double* smoothSpectrum(double* &powerSpec, bool fullLength = false, unsigned int M = 0, bool replace = false);

    //====================================================================================================================
    // Gets the signal of length "len", returns it's discrete cosine transform in an 1-dimensional array 
    // After: http://en.wikipedia.org/wiki/Discrete_cosine_transform
    // If coeffCount > 0, only so much coefficients are computed and returned (max. len nr. of coefficients)
    //====================================================================================================================
    double* dct(double* signal, unsigned long int len, unsigned long coeffCount = 0);

    //====================================================================================================================
    // Gets the cosine-transformed in an 1-dimensional array of length "len", returns the signal
    // After: http://en.wikipedia.org/wiki/Discrete_cosine_transform
    //====================================================================================================================
    double* idct(double* cosineTransform, unsigned long int len);

    //====================================================================================================================
    //	generates normalized correlation between two adjacent frames (means: they are non-overlapping and stored in con-
    //  catenated form in the twoFrames-array, with twoFrames[frameLegth] being the begining of the second frame) of 
    //  length frameLength; they are regarded two be the previous and current frame;
    //  the result has frameLength-1 elements, because no correlation needs to be computed for 0-shift
    //====================================================================================================================
    double* normalizedCorrelation(double* twoFrames, unsigned long int frameLength);

		//====================================================================================================================
		// Cepstrum: Real part of the inverse of  the IFFT of the log-power spectrum; the result has length fftlen
		// This version equals the rceps() function in matlab with to modifications: 
		//   First, here tapering is applied according to the caller's wishes
		//   Second, power-spectrum instead of magnitude spectrum is used before log() as in the paper "Cepstrum-based pitch 
		//   detection using a new statistical v/uv classification algorithm", ahmadi, spanias 1999
		//=====================================================================================================================
		double* cepstrum(double * signal, unsigned long int len);
		
    //====================================================================================================================
    // the discrete wavelet transform (dwt); it subsequqntly uses high- low-pass operations; it is implemented with a 
		// Daubechies filterbank of specifiyable order; the level determines the number of subsequent subband-generations
    //====================================================================================================================
	  double** dwt(double* signal, unsigned long int len, unsigned int level, unsigned int daubechiesOrder = 2);	

		//=====================================================================================================================
		// takes a power spectrum of size fftLen/2 +1 and creaters the full length version using the conjugate complex 
		// relationship between the samples; the newly created and returned power spectrum is of size fftLen; if 
		// replaceOriginal is true, the original pointer is freed and replaced by the new buffer
		//=====================================================================================================================
		double* powerSpec2fullPowerSpec(double* &powerSpec, bool replaceOriginal = false);
};

#endif
