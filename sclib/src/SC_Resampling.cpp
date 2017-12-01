/**************************************************************************/
/*    This class is an object-orientated wrapper around libsamplerate     */
/*    0.1.2                                                               */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 31.07.2006																								*/
/**************************************************************************/

#include <stdlib.h> //for NULL macro
#include "SC_Resampling.h"
#include "SC_Api.h"
#include "SC_Aux.h" //for sclib::bufferSize define
#include <SV_Error.h>
#include <GN_Filter.h>

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Resampling::SC_Resampling(SC_TweakableParameters *pTweak) {
#ifdef SC_USE_LIBSAMPLERATE
	int err;
	char errString[sclib::bufferSize];

	this->pTweak = pTweak;

	//medium quality sinc resampling is a lot faster, but only of neglectable lesser quality than the best quality sinc resampler
	//only one channel is supported 'cause otherwise the filtering using the GN_Filter class will fail
	this->pSRCstate = src_new(((this->pTweak != NULL && this->pTweak->resampling.fastConversion == true) ? SRC_SINC_MEDIUM_QUALITY : SRC_SINC_BEST_QUALITY), 1, &err);
	if (this->pSRCstate == NULL) {
		sprintf(errString, "%s\0", src_strerror(err));
		REPORT_ERROR(SVLIB_Fail, errString);
	}
#else
	REPORT_ERROR(SVLIB_Fail, "libsamplerate not available!");
#endif
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_Resampling::~SC_Resampling() {
#ifdef SC_USE_LIBSAMPLERATE
	this->pSRCstate = src_delete(this->pSRCstate);
#endif
}

//====================================================================================================================
//	Resample the given samples at the given sample-rate and do filtering to prevent from aliasing (only needed for not 
//  yet implementd simpler resamplers, so at the moment it's not needed), if wished. The new sample-rate has to be 
//  given in Hz (e.g. 4000 for a sample-rate of 4KHz). The memory for the new buffer is allocated within this method 
//  if outBuffer==NULL, otherwise it is used (and it is assumed than that it is big enough); in all cases the pointer 
//  to the result is returned; if inputFinished==false, the output is delayed (i.e. not as much samples are returned 
//  as would be expected but are cached for the next call to avoid discontinueities if maybe short frames are 
//  resampled consecutively) and completely returned when a last call with inputFinished==true is done.
//  TODO: Do gain normalization after complete signal has been resampled (not just framewise)
//====================================================================================================================
short* SC_Resampling::resample(short *inBuffer, unsigned long int length, double sampleRate, double newSampleRate, unsigned long int &newLength, short *outBuffer, bool inputFinished, bool doFiltering) {
#ifdef SC_USE_LIBSAMPLERATE
  short *newBuffer = outBuffer, tmp;
  float *fBufIn, *fBufOut;
  int res = 0;
  char errStr[sclib::bufferSize];
  SRC_DATA srcData;
  double ratio, cutOffFreq, gradient;
  GN_Filter FIR;

  ratio = newSampleRate / sampleRate;
  newLength = (unsigned long int)(ceil(length * ratio) + 1000.0); //make longer than necessary for internal caching if inputFinished==false
  
  if (src_is_valid_ratio(ratio)) {
    MArray_1D(fBufIn, length, float, "SC_Resampling.resample: fBufIn");
    MArray_1D(fBufOut, newLength, float, "SC_Resampling.resample: fBufOut");
    src_short_to_float_array(inBuffer, fBufIn, length);

    if (doFiltering == true) { //only needed for non-sinc resamplers
      cutOffFreq = (((newSampleRate/2.0)-1.0) / sampleRate) * 2.0; //cut off 1Hz below the Nyquist Frequency of the new samplerate, normalized to 1Hz
      if (cutOffFreq > 0.0 && cutOffFreq < 1.0) {
        FIR.LowPass(cutOffFreq);
        FIR.FIR_Filter(fBufIn, length);
      }
    }

    srcData.data_in = fBufIn;
    srcData.input_frames = length;
    srcData.data_out = fBufOut;
    srcData.output_frames = newLength;
    srcData.src_ratio = ratio;
		srcData.end_of_input = (inputFinished == true) ? 1 : 0; //output will be dalyed till the next call of this function at the benefit of no discontinueities between resampled frames

		res = src_process(this->pSRCstate, &srcData);
    if (res > 0) {
      sprintf(errStr, "%s\0", src_strerror(res));
      REPORT_ERROR(SVLIB_Fail, errStr);
			newLength = 0;
		} else {
			MFree_1D(fBufIn);
			newLength = (inputFinished == true) ? srcData.output_frames_gen + 1 : srcData.output_frames_gen; //somehow, there is always one sample too less for finished data... we correct that later
			if (newBuffer == NULL) { //othwerwise we use the memory already allocated by the caller and assume it's big enough
				MArray_1D(newBuffer, newLength, short, "SC_Resampling.resample: newBuffer");
			}
			src_float_to_short_array(fBufOut, newBuffer, newLength);
			MFree_1D(fBufOut);

			if (inputFinished == true) {
				//interpolate last sample, place it on second-last position because it is a linear interpolation between the last 2 new samples
				gradient = (double)(newBuffer[newLength-2] - newBuffer[newLength-3]);
				tmp = newBuffer[newLength-2];
				newBuffer[newLength-2] = (short)(floor(newBuffer[newLength-3] + (gradient * ratio))); //(newBuffer[newLength-3] + tmp) / 2;
				newBuffer[newLength-1] = tmp;
			}
		}
  } else {
    REPORT_ERROR(SVLIB_BadArg, "Resampling not possible due to invalid sample rate ratio");
    newLength = 0;
  }

	//reset internal cache if the input has come to an end
	if (inputFinished == true) {
		res = src_reset(this->pSRCstate);
    if (res > 0) {
      sprintf(errStr, "%s\0", src_strerror(res));
      REPORT_ERROR(SVLIB_Fail, errStr);
		}
	}

  return newBuffer;
#else
	return NULL;
#endif
}

//====================================================================================================================
//	The channels parameter must be set to the number of channels in the input signal; the resampler is not capable of 
//  doing channel number conversions; according to mode, all available channels are averaged, or one of the channels
//  is selected to become the mono signal; the outBuffer must be allocated outside this function; sampleCoutn refers 
//  to the number of samples at different instants of time, so it is the same for the multi-channel and the mono
//  signal; returns false on error (wrong selected channel)
//====================================================================================================================
bool SC_Resampling::toMono(short *inBuffer, short *outBuffer, unsigned long int sampleCount, unsigned int originalChannelCount, unsigned int mode) {
	unsigned long int x, z = 0;
	unsigned int c;
	bool res = true;

	if (mode == sclib::resampleAverage) {
		for (x = 0; x < sampleCount; x++) {
			outBuffer[x] = 0;
			for (c = 0; c < originalChannelCount; c++) {
				outBuffer[x] += inBuffer[z++];
			}
			outBuffer[x] /= originalChannelCount;
		}
	} else {
		if (originalChannelCount <= mode) {
			c = mode--; //convert wanted channel to an offset into the array
			for (x = 0; x < sampleCount; x++) {
				outBuffer[x] = inBuffer[x+c];
			}
		} else {
			res = false;
		}
	}

	return res;
}
