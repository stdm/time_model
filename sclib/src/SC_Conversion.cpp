/**************************************************************************/
/*	This implements conversions between different measures (time-units,   */
/*  frame-sizes, ...)                                                     */
/*  Most are moved from SC_Conversion* here because they are more         */
/*  generally applicabel than within that old scope.                      */
/*																																				*/
/*  Attention: this projections are not(!) bidirectional identical under 	*/
/*            all circumstances!!!																		    */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 19.03.2007																								*/
/**************************************************************************/

#include "SC_Conversion.h"
#include "SC_Aux.h" //for sclib::round()
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Conversion::SC_Conversion(unsigned long int audioSampleRate, double videoFrameSize) {
	this->audioSampleRate = audioSampleRate;
	this->videoFrameSize = videoFrameSize;

	if (this->videoFrameSize == 0.0) {
		this->videoFrameSize = this->audioSampleRate / 1000.0; //set a standard video frame size of 1 ms (in samples)
	}
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Conversion::~SC_Conversion() {

}

//====================================================================================================================
// Set the audio sample rate; if the video frame size is 0.0, also set it to be 1ms (in samples)
//====================================================================================================================
void SC_Conversion::setAudioSampleRate(unsigned long int newSampleRate) {
	this->audioSampleRate = newSampleRate;

	if (this->videoFrameSize == 0.0) {
		this->videoFrameSize = this->audioSampleRate / 1000.0; //set a standard video frame size of 1 ms (in samples)
	}

	return;
}

//====================================================================================================================
//	gives the time in ms of the starting (alignment = sclib::alignmentStart) or ending (alignment = 
//  sclib::alignmentEnd) of the frame; first frame is 0
//  frame is just the frame-number; audioFrameSize and audioFrameStep are given in samples
//====================================================================================================================
unsigned long int SC_Conversion::audioFrame2ms(unsigned long int frame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment) {
  unsigned long int sample = audioFrame2sample(frame, audioFrameSize, audioFrameStep, alignment);
  unsigned long int ms = sample2ms(sample);

  return ms;
}

//====================================================================================================================
//	gives the frame-nr (according to audioFrameSize and audioFrameStep given in samples) "nearest" to the specified
//  time in ms. "near" is specified by the alignment-parameter: 
//   - sclib::alignmentStart: gives the audioFrame who's beginning is nearest to the ms
//   - sclib::alignmentEnd: gives the audioFrame who's end is nearest to the ms
//  first frame is 0
//  frame is just the frame-number; audioFrameSize and audioFrameStep are given in samples
//====================================================================================================================
unsigned long int SC_Conversion::ms2audioFrame(unsigned long int ms, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment) {
  unsigned long int sample = ms2sample(ms, alignment);
  unsigned long int frame = sample2audioFrame(sample, audioFrameSize, audioFrameStep, alignment);

  return frame;
}

//====================================================================================================================
//	gives the first (if alignment = sclib::alignmentStart) or last (if alignment = sclib::alignmentEnd) sample in the 
//  frame (specified by audioFrameSize and audioFrameStep in samples); first frame is 0
//====================================================================================================================
unsigned long int	SC_Conversion::audioFrame2sample(unsigned long int frame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment) {
  switch (alignment) {
    case sclib::alignmentStart: {
      return frame * audioFrameStep; //first sample of the given frame
      break;
    }
    case sclib::alignmentEnd: {
      return (frame * audioFrameStep) + audioFrameSize - 1; //last sample of the given frame
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Type of alignment unknown");
      return frame * audioFrameStep; //first sample of the given frame, to return a somehow meaningful value
      break;
    }
  }
}

//====================================================================================================================
//	gives the frame-nr (according to audioFrameSize and audioFrameStep given in samples) "nearest" to the specified
//  sample. "near" is specified by the alignment-parameter: 
//   - sclib::alignmentStart: gives the audioFrame who's beginning is nearest to the sample
//   - sclib::alignmentEnd: gives the audioFrame who's end is nearest ti the sample
//  first sample is 0
//  frame is just the frame-number; audioFrameSize and audioFrameStep are given in samples
//  ATTENTION: the returned frame contains the given sample, so it's likely that the frame goes beyond or starts 
//  before the given sample; if you want the frame just starting/ending really before/after the sample, add/substract 
//  1 to/from the result!!!
//  if a sampleRate is given explicitly, it is used instead of the internal sampleRate
//====================================================================================================================
unsigned long int	SC_Conversion::sample2audioFrame(unsigned long int sample, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment) {
  switch (alignment) {
    case sclib::alignmentStart: {
      return sample / audioFrameStep; //most near frame-start to the left of a frame containing this sample; integer-division is correct for the first frame has number 0
      break;
    }
    case sclib::alignmentEnd: {
			return (sample < audioFrameSize) ? 0 : ((sample - audioFrameSize) / audioFrameStep) + 1; //most near frame-end to the right of a frame containing this sample
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Type of alignment unknown");
      return sample / audioFrameStep; //to return a somehow meaningful value...
      break;
    }
  }
}

//====================================================================================================================
//	converts ms to samples; first ms is 0; according to alignment, the first or last sample in the ms-period is 
//  returned; if a sampleRate is given explicitly, it is used instead of the internal sampleRate
//====================================================================================================================
unsigned long int SC_Conversion::ms2sample(unsigned long int ms, unsigned short int alignment, unsigned long int sampleRate) {
	unsigned int samplesPerMs = ((sampleRate > 0) ? sampleRate : this->audioSampleRate) / 1000;

  switch (alignment) {
    case sclib::alignmentStart: {
			return ms * samplesPerMs; //first sample of the given ms
      break;
    }
    case sclib::alignmentEnd: {
     	return (ms * samplesPerMs) + samplesPerMs - 1; //last sample of the given ms
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Type of alignment unknown");
 			return ms * samplesPerMs; //first sample of the given ms, to return a somehow meaningful value
      break;
    }
  }
}

//====================================================================================================================
//	converts ms to samples; first ms is 0; according to alignment, the first or last sample in the ms-period is 
//  returned; if a sampleRate is given explicitly, it is used instead of the internal sampleRate
//  in this version, ms can be a float (ms are partable, only samples are not!)
//====================================================================================================================
unsigned long int SC_Conversion::ms2sample_f(float ms, unsigned short int alignment, unsigned long int sampleRate) {
	unsigned int samplesPerMs = ((sampleRate > 0) ? sampleRate : this->audioSampleRate) / 1000;

  switch (alignment) {
    case sclib::alignmentStart: {
			return (unsigned long int)(ms*samplesPerMs); //first sample of the given ms (downrounded)
      break;
    }
    case sclib::alignmentEnd: {
			return (unsigned long int)(ms*samplesPerMs + samplesPerMs-1 + 0.5); //last sample of the given ms (uprounded)
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Type of alignment unknown");
			return (unsigned long int)(ms*samplesPerMs); //first sample of the given ms, to return a somehow meaningful value
      break;
    }
  }
}

//====================================================================================================================
//	converts samples to ms; first sample is 0; if a sampleRate is given explicitly, it is used instead of the internal 
//  sampleRate
//====================================================================================================================
unsigned long int SC_Conversion::sample2ms(unsigned long int sample, unsigned long int sampleRate) {
  unsigned int samplesPerMs = ((sampleRate > 0) ? sampleRate : this->audioSampleRate) / 1000;

  return sample / samplesPerMs; //yes, integer-division is correct here because first ms is also 0
}

//====================================================================================================================
//	gives the time in ms of the starting/ending (according to alignment) of the frame; first frame is 0
//====================================================================================================================
unsigned long int SC_Conversion::videoFrame2ms(unsigned long int frame, unsigned short alignment) {
  unsigned long int sample = videoFrame2sample(frame, alignment);
  unsigned long int ms = sample2ms(sample);

  return ms;
}

//====================================================================================================================
//	just a shortcut to ms2videoFrame(audioFrame2ms(audioFrame))
//====================================================================================================================
unsigned long int SC_Conversion::audioFrame2videoFrame(unsigned long int audioFrame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment) {
  return ms2videoFrame(audioFrame2ms(audioFrame, audioFrameSize, audioFrameStep, alignment));
}

//====================================================================================================================
//	just a shortcut to ms2audioFrame(videoFrame2ms9(videoFrame))
//====================================================================================================================
unsigned long int SC_Conversion::videoFrame2audioFrame(unsigned long int videoFrame, unsigned int audioFrameSize, unsigned int audioFrameStep, unsigned short alignment) {
  return ms2audioFrame(videoFrame2ms(videoFrame), audioFrameSize, audioFrameStep);
}

//====================================================================================================================
//	first ms is 0
//	gives the frame within which this ms lies
//====================================================================================================================
unsigned long int SC_Conversion::ms2videoFrame(unsigned long int ms, unsigned short alignment) {
  unsigned long int sample = ms2sample(ms, alignment);
  unsigned long int frame = sample2videoFrame(sample, alignment);

  return frame;
}

//====================================================================================================================
//	convert video-frame to sample; first frame and sample is 0
//====================================================================================================================
unsigned long int SC_Conversion::videoFrame2sample(unsigned long int frame, unsigned short int alignment) {
  switch (alignment) {
    case sclib::alignmentStart: {
      return sclib::round(frame * this->videoFrameSize); //first sample of the given frame
      break;
    }
    case sclib::alignmentEnd: {
      return sclib::round((frame * this->videoFrameSize) + this->videoFrameSize - 1); //last sample of the given frame
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Type of alignment unknown");
      return sclib::round(frame * this->videoFrameSize); //first sample of the given frame, to return a somehow meaningful value
      break;
    }
  }
}

//====================================================================================================================
//	convert sample to video-frame; first frame and sample is 0
//  if sampleRate!=0, this sampleRate is used instead of the internal one
//====================================================================================================================
unsigned long int SC_Conversion::sample2videoFrame(unsigned long int sample, unsigned short int alignment, unsigned long int sampleRate) {
	double vfs = (sampleRate > 0) ? this->videoFrameSize * ((double)(this->audioSampleRate) / (double)(sampleRate)) : this->videoFrameSize;

  switch (alignment) {
    case sclib::alignmentStart: {
      return sclib::round(sample / vfs); //most near frame-start to the left of a frame containing this sample; integer-division is correct for the first frame has number 0
      break;
    }
    case sclib::alignmentEnd: {
			return (sample < vfs) ? 0 : sclib::round(((sample - vfs) / vfs)) + 1; //most near frame-end to the right of a frame containing this sample
      break;
    }
    default: {
      REPORT_ERROR(SVLIB_BadArg, "Type of alignment unknown");
      return sclib::round(sample / vfs); //to return a somehow meaningful value...
      break;
    }
  }
}

//====================================================================================================================
// converts a frequency (in [Hz]) into an index into a powerspectrum
// it must hold: 0 <= frequency <= sampleRate/2
//====================================================================================================================
unsigned int SC_Conversion::frequency2index(double frequency, unsigned int fftLength) {
  double f_max, f_delta, idx, length, f;
	//here it doesn't matter if the spectrum is mormal sized (fftLength/2 +1) or fullLength (fftLength), because only the normal-sized part has non-redundant information

	if (frequency < 0.0 || frequency > this->audioSampleRate/2.0) {
		REPORT_ERROR(SVLIB_BadArg, "it must hold: 0 <= frequency <= sampleRate/2");
		f = sclib::getBetween(0.0, frequency, this->audioSampleRate/2.0);
	} else {
		f = frequency;
	}

	length = (double)(fftLength) / 2.0; //usable part of the spectrum is (fftLength/2)+1 long (idx 0 has frequency 0, idx ffLength/2 has the nyquist-frequency of sampleRate/2)
  f_max = (double)(this->audioSampleRate) / 2.0; //nyquist frequency
  f_delta = f_max / length; //increment from each index to the next

  idx = f / f_delta;
  
  return (unsigned long)sclib::round(idx);
}

//====================================================================================================================
// converts the index into a powerspectrum into the corresponding discrete frequency
// it must hold: 0 <= index <= fftLength/2
//====================================================================================================================
double SC_Conversion::index2frequency(unsigned int index, unsigned int fftLength) {
  double f_max, f_delta, freq, length;
	unsigned int i;
	//here it doesn't matter if the spectrum is normal sized (fftLength/2 +1) or fullLength (fftLength), because only the normal-sized part has non-redundant information

	if (index < 0 || index > fftLength/2) {
		REPORT_ERROR(SVLIB_BadArg, "it must hold: 0 <= index <= fftLength/2");
		i = (unsigned int)(sclib::getBetween(0.0, index, fftLength/2.0));
	} else {
		i = index;
	}	

	length = (double)(fftLength) / 2.0; //usable part of the spectrum is (fftLength/2)+1 long (idx 0 has frequency 0, idx ffLength/2 has the nyquist-frequency of sampleRate/2)
  f_max = (double)(this->audioSampleRate) / 2.0; //nyquist frequency
  f_delta = f_max / length; //increment from each index to the next

  freq = (double)(i) * f_delta;
  
  return freq;
}

//====================================================================================================================
// converts a quefrency (in [ms]) into an index into a cepstrum
// it must hold: 0 <= quefrency <= 1000*((fftLength/2)/sampleRate)
//====================================================================================================================
unsigned int SC_Conversion::quefrency2index(double quefrency, unsigned int fftLength) {
  double idx, q;

	if (quefrency < 0 || quefrency > (500.0*fftLength)/(double)(this->audioSampleRate)) {
		REPORT_ERROR(SVLIB_BadArg, "0 <= quefrency <= 1000*((fftLength/2)/sampleRate)");
		q = sclib::getBetween(0.0, quefrency, (500.0*fftLength)/(double)(this->audioSampleRate));
	} else {
		q = quefrency;
	}	

	idx = (q / 1000.0) * this->audioSampleRate;
 
  return (unsigned long)sclib::round(idx);
}

//====================================================================================================================
// converts the index into a cepstrum into the corresponding discrete quefrency (in [ms])
// it must hold: 0 <= index <= fftLength/2
//====================================================================================================================
double SC_Conversion::index2quefrency(unsigned int index, unsigned int fftLength) {
  double q;
	unsigned int i;

	if (index < 0 || index > fftLength/2) {
		REPORT_ERROR(SVLIB_BadArg, "it must hold: 0 <= index <= fftLength/2");
		i = (unsigned int)(sclib::getBetween(0.0, index, fftLength/2.0));
	} else {
		i = index;
	}	

  q = 1000.0 * ((double)(i) / (double)(this->audioSampleRate));
  
  return q;
}

//====================================================================================================================
// converts quefrency (in [ms]) to frequency (in [Hz])
//====================================================================================================================
double SC_Conversion::quefrency2frequency(double quefrency) {
	return (quefrency > 0.0) ? 1000.0 / quefrency : 0.0;
}

//====================================================================================================================
// converts frequency (in [Hz]) to quefrency (in [ms])
//====================================================================================================================
double SC_Conversion::frequency2quefrency(double frequency) {
	return (frequency > 0.0) ? 1000.0 / frequency : 0.0;
}
