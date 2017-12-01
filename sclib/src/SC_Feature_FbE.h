/**************************************************************************/
/*    Copied from:																												*/
/*      - This is a altered version of SV_Feature_MFCC	by Jialong He			*/
/*																																				*/
/*    Responsibility:																											*/
/*      - Class for extracting Filterbank Energys      										*/
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

#ifndef __SC_Feature_FbE_H__
#define __SC_Feature_FbE_H__

#include <SV_Feature.h>

class SC_Feature_FbE : public SV_Feature {

	private :

	protected :

    //====================================================================================================================
    // members which holds the desired scale for filterbank-computation (linear/mel/mod_mel/expolog) and the desired
    // result-type
    //====================================================================================================================
    unsigned char scale;
    unsigned char result;
    unsigned int  smoothSpectrum;

    //====================================================================================================================
    // members which define the lower and upper frequency-bound to be covered by the filterbank
    // standard is from 0.0 to half the sample-rate
    //====================================================================================================================
    double lowFrequency;
    double highFrequency;

    //====================================================================================================================
    // returns a pointer to the frequency-convert-function Hz->mode or vice versa
    //====================================================================================================================
    static double (*getFrequencyMappingHz2X(unsigned char mode))(double);
    static double (*getFrequencyMappingX2Hz(unsigned char mode))(double);

		//====================================================================================================================
    // this affects how much percent (0..1) of a frame a affected by the window: 100% (1) -> all the frame, 10% (0.1) ->
		// only the first 10th of the beginning and the last 10th of the end of the frame are drwan smoothly to zero
    //====================================================================================================================
		double taperingLength;
    
	public :

    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_FbE(int sampleRate, int frameSize, int frameStep, int filterBankSize, int FFTsize, int window, double preEmphasizeFactor, double lowestFrequency, double highestFrequency, int MFCCorder, bool dEnergy, unsigned char scale = sclib::scaleLinear, unsigned char result = sclib::resultLinear, unsigned int smoothing = sclib::smoothNone, double taperingLength = 1.0);
		virtual ~SC_Feature_FbE();

    //====================================================================================================================
    // some setters
    //====================================================================================================================
    void setLowFrequency(double fHz) {assert(fHz >= 0.0 && fHz < this->Para.SRate/2.0); this->lowFrequency = fHz; return;};
    void setHighFrequency(double fHz) {assert(fHz > 0.0 && fHz <= this->Para.SRate/2.0); this->highFrequency = fHz; return;};

    //====================================================================================================================
    // converts linear frequency (Hz) to linear frequeny (Hz)  (the identity function: silly, but needed :-)
    //====================================================================================================================
    static double Hz2Hz(double f);

    //====================================================================================================================
    // converts linear frequency (Hz) to mel
    //====================================================================================================================
    static double Hz2mel(double f);

    //====================================================================================================================
    // converts mel to linear frequency (Hz)
    //====================================================================================================================
    static double mel2Hz(double mel);

    //====================================================================================================================
    // converts linear frequency (Hz) to modified mel
    //====================================================================================================================
    static double Hz2mMel(double f);

    //====================================================================================================================
    // converts modified mel to linear frequency (Hz)
    //====================================================================================================================
    static double mMel2Hz(double mMel);

    //====================================================================================================================
    // converts linear frequency (Hz) to  expoLog
    //====================================================================================================================
    static double Hz2expoLog(double f);

    //====================================================================================================================
    // converts expoLog to linear frequency (Hz)
    //====================================================================================================================
    static double expoLog2Hz(double expoLog);

    //====================================================================================================================
    // converts linear frequency (Hz) to bark according to H. Traunmüller (1990) "Analytical expressions for the tonotopic 
		// sensory scale" J. Acoust. Soc. Am. 88: 97-100. 
    //====================================================================================================================
    static double Hz2bark(double f);

    //====================================================================================================================
    // converts bark to linear frequency (Hz)
		// see e.g. http://www.ling.su.se/staff/hartmut/bark.htm
    //====================================================================================================================
    static double bark2Hz(double z);

    //====================================================================================================================
		// returns, for a frequency in [Hz], its critical bandwidth in [Hz] according to the Bark scale; if usePublishedScale
		// is true, not the bandwidth of exact the given frequency f is returned, but the bandwidth of the bark band this 
		// frequency falls into as originally published is returned; additionally, the centerFrequancy xorresponding with this 
		// band is returned.
    //====================================================================================================================
    static double getCriticalBandwidth(double f, double &centerFrequency, double &leftEdge, double &rightEdge, bool usePublishedScale = false);

    //====================================================================================================================
    // returns, for a frequency in [Hz], its ERB in [Hz] according to B.C.J. Moore and B.R. Glasberg (1983) "Suggested 
		// formulae for calculating auditory-filter bandwidths and excitation patterns" J. Acoust. Soc. Am. 74: 750-753. 
		// the formula is valid for frequencies in the range 0.1-6.5kHz
    //====================================================================================================================
		static double getEquivalentRectangularBandwidth(double f);

    //====================================================================================================================
    // returns an array of the filter-parameters for filterbank spaced according to "mode" and the lowest and highest 
    // frequency to cover (the upper and lower end are the center-frequencies of the neighboring filters):
    // the first component holds the lower end frequency (Hz)
    //     second                    center-frequency (Hz)
    //     third                     upper end frequency (Hz)
    //====================================================================================================================
    static double** createSpacing(unsigned char mode, double fMinHz, double fMaxHz, int filterCount);

    //====================================================================================================================
    // evaluates the specified filter for the given spectrum. filterCenter is the filter's center frequency and 
    // filterMin/-Max or the left and right point where the filter reaches zero, all as returned by createSpacing()
    //
    // this is taken from m.slaneys auditory toolbox for matlab (calculation of mfccFilterWeights) with some changes: 
    // there, frequency = (idx * sampleRate) / fftSize, not frequency = (idx * sampleRate) / (2 * fftSize)
    //====================================================================================================================
    static double evalFilter(double* powerSpec, double filterMin, double filterCenter, double filterMax, double sampleRate, unsigned long int fftSize);

		//====================================================================================================================
		// returns the filter coefficient for the given triangular filter and the proposed idx into a (power-)spectrum. just
		// the same routines as in the other evalFilter() method.
		//====================================================================================================================
		static double evalFilter(unsigned int idx, double filterMin, double filterCenter, double filterMax, double sampleRate, unsigned long int fftSize);

    //====================================================================================================================
    // prints the current/desired filterbank-spacing to an ASCII-file
    //====================================================================================================================
    void spacingOut(void);
    void spacingOut(unsigned char mode, double fminHz, double fmaxHz);

    //====================================================================================================================
		// override base class method, return LogFilterbankEnergy vector sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
};

#endif
