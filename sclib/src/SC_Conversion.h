/**************************************************************************/
/*	This implements conversions between different measures (time-units,   */
/*  frame-sizes, ...)                                                     */
/*  Most are moved from SC_GroundTruth* here because they are more        */
/*  generally applicabel than within that old scope.                      */
/*																																				*/
/*  Attention: this projections are not(!) bidirectional identical under 	*/
/*            all circumstances!!!																		    */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.03.2007																								*/
/**************************************************************************/

#ifndef __SC_Conversion_H__
#define __SC_Conversion_H__

#include "SC_Api.h"
#include "SC_Aux.h"

class SCLIB_API SC_Conversion {

	private :

	protected:

		unsigned long int audioSampleRate;
		double videoFrameSize; //length of a videoframe in samples

	public :

	  SC_Conversion(unsigned long int audioSampleRate = 0, double videoFrameSize = 0.0);
		virtual ~SC_Conversion();

		//====================================================================================================================
		// Getter and setter
		//====================================================================================================================
		void setAudioSampleRate(unsigned long int newSampleRate);
		void setVideoFrameSize(double newFrameSize) {this->videoFrameSize = newFrameSize; return;}

		//====================================================================================================================
		// Methods to convert between the different time measurement units
		//====================================================================================================================
		unsigned long int audioFrame2ms(unsigned long int frame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment = sclib::alignmentStart);
		unsigned long int ms2audioFrame(unsigned long int ms, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment = sclib::alignmentStart);
		unsigned long int audioFrame2sample(unsigned long int frame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment = sclib::alignmentStart);
		unsigned long int sample2audioFrame(unsigned long int sample, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment = sclib::alignmentStart);
    unsigned long int sample2ms(unsigned long int sample, unsigned long int sampleRate = 0);
		unsigned long int	videoFrame2ms(unsigned long int frame, unsigned short alignment = sclib::alignmentStart);
		unsigned long int	ms2videoFrame(unsigned long int ms, unsigned short alignment = sclib::alignmentStart);
    unsigned long int videoFrame2sample(unsigned long int frame, unsigned short int alignment = sclib::alignmentStart);
    unsigned long int sample2videoFrame(unsigned long int sample, unsigned short int alignment = sclib::alignmentStart, unsigned long int sampleRate = 0);
    unsigned long int audioFrame2videoFrame(unsigned long int audioFrame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment = sclib::alignmentStart);
    unsigned long int videoFrame2audioFrame(unsigned long int videoFrame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment = sclib::alignmentStart);

		unsigned long int ms2sample(unsigned long int ms, unsigned short int alignment = sclib::alignmentStart, unsigned long int sampleRate = 0);
		unsigned long int ms2sample_f(float ms, unsigned short int alignment = sclib::alignmentStart, unsigned long int sampleRate = 0);
    
 		//====================================================================================================================
		// converts the index into a powerspectrum into the corresponding discrete frequency (in [Hz]) and vice versa
		// here it doesn't matter if the spectrum is normal sized (fftLength/2 +1) or fullLength (fftLength), because only the 
		// normal-sized part has non-redundant information;
		// it holds: 0 <= index <= fftLength/2, 0 <= frequency <= sampleRate/2
		//====================================================================================================================
    unsigned int frequency2index(double frequency, unsigned int fftLength);
    double index2frequency(unsigned int index, unsigned int fftLength);

 		//====================================================================================================================
		// converts the index into a cepstrum into the corresponding discrete quefrency (in [ms]) and vice versa
		// a fullength spectrum has fftLength entries, but the 2nd half is redundant
		//====================================================================================================================
    unsigned int quefrency2index(double quefrency, unsigned int fftLength);
    double index2quefrency(unsigned int index, unsigned int fftLength);

 		//====================================================================================================================
		// converts quefrency (in [ms]) to frequency (in [Hz]) and vice versa
		//====================================================================================================================
		double quefrency2frequency(double quefrency);
		double frequency2quefrency(double frequency);
};

#endif
