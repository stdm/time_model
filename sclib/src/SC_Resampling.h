/**************************************************************************/
/*    This class is an object-orientated wrapper around libsamplerate     */
/*    0.1.2                                                               */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 31.07.2006																								*/
/**************************************************************************/

#ifndef __SC_Resampling_H__
#define __SC_Resampling_H__

#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
#include "SC_Api.h"
#ifdef SC_USE_LIBSAMPLERATE
	#include <samplerate.h>
#endif
#include <SV_Error.h>

class SCLIB_API SC_Resampling {
  private:

  protected:
		
#ifdef SC_USE_LIBSAMPLERATE
		SRC_STATE	*pSRCstate;
#endif
    SC_TweakableParameters *pTweak;

  public:

    SC_Resampling(SC_TweakableParameters *pTweak = NULL);
    ~SC_Resampling();

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
    short* resample(short *inBuffer, unsigned long int length, double sampleRate, double newSampleRate, unsigned long int &newLength, short *outBuffer = NULL, bool inputFinished = true, bool doFiltering = false);

    //====================================================================================================================
    //	The channels parameter must be set to the number of channels in the input signal; the resampler is not capable of 
		//  doing channel number conversions; according to mode, all available channels are averaged, or one of the channels
		//  is selected to become the mono signal; the outBuffer must be allocated outside this function; sampleCount refers 
		//  to the number of samples at different instants of time, so it is the same for the multi-channel and the mono
		//  signal; returns false on error (wrong selected channel)
    //====================================================================================================================
		bool toMono(short *inBuffer, short *outBuffer, unsigned long int sampleCount, unsigned int originalChannelCount, unsigned int mode = sclib::resampleAverage);
};

#endif
