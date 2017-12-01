/**************************************************************************/
/*    Responsibiloity:																										*/
/*      - return an SC_Signal object able to read/write the specified     */
/*        audio-type                                                      */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2005																								*/
/**************************************************************************/

#ifndef __SC_SignalHandler_H__
#define __SC_SignalHandler_H__

#include "SC_Api.h"
#include "SC_Signal.h"
#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Signal_MPEG.h"
#include <SV_Data.h>

class SCLIB_API SC_SignalHandler {
	private :

	protected :

    unsigned short int signalType;
    SC_TweakableParameters *pTweak;

    //====================================================================================================================
		//	Hold a permanent copy of a SC_Signal_MPEG object as a hack for the non-working seeking in ffmpeg: If consecutive 
		//  segments should be loaded, the current state of this permanent object can be exploted to avoid decoding from the
		//  very beginning
		//====================================================================================================================
		SC_Signal_MPEG *pMPEG; 

    //====================================================================================================================
		//	guess the signal-type from the filename's extension or the files content,if the extension is ambigous
		//====================================================================================================================
    unsigned short int guessSignalType(const char *fileName);

		//====================================================================================================================
		// Meta-loader of SC_Signal_MPEG based signals: Instead of just opening the signal with a new (fresh-state) object,
		// this method uses a permanent SC_Signal_MPEG object that loads all content and then links it to the given pSignal or
		// sets it NULL in case of error; this way, the seeking problem of ffmpeg can be "workarounded" for the prominent case
		// that the nect consecutive call of the load-method just wants the segment that directly succeeds the previously
		// opened one
		//====================================================================================================================
		void loadMPEG(unsigned long int start, unsigned long int end, SC_Signal_MPEG* &pSignal);

	public :
				
		SC_SignalHandler(SC_TweakableParameters *pTweak, unsigned short int signalType = sclib::stGuess);
    virtual ~SC_SignalHandler();

    void setSignalType(unsigned short int newType) {this->signalType = newType; return;}

    //====================================================================================================================
    //	just opens the signal so that it's parameters can be extracted; if forceSampleRate>0, the signal (here, at least, 
		//  the parameters in the header) will be resampled to the given samplerate regardless of it's original sampleRate
    //====================================================================================================================
    SC_Signal* openSignal(const char* fileName, int forceSampleRate = 0, void* jniEnvironment = NULL, void* jStreamObject = NULL);

    //====================================================================================================================
    //	load the speech-signal between the given segment-boundarys (in samples)
    //  install also new debugprefix, if wished in the tweakable parameters, according to the fileName
		//  if forceSampleRate>0, the signal will be resampled to the given samplerate regardless of it's original sampleRate
    //====================================================================================================================
    SC_Signal* loadSignal(const char* fileName, unsigned long int segmentStart = 0, unsigned long int segmentEnd = 0, int forceSampleRate = 0);

    //====================================================================================================================
    // write an audio-file only incorporating samples between segmentStart and segmentEnd of specific type
    //         -"-         containing the samples in the given column of a linked list of sv_data-objects
    //         -"-         containing the samples in the given row (all rows for row==-1) of a sv_data-object
    //         -"-         containing the samples in the signal-buffer of given length, destructing signal-pointer!
    //====================================================================================================================
    long storeSignal(const char* fileName, unsigned long int	segmentStart, unsigned long int	segmentEnd, SC_Signal *pSigPar, SC_GroundTruth *pGT, unsigned long int type = sclib::noType, unsigned long int typesNot = sclib::noType, int origin = sclib::modeHypothesized, unsigned long int markGTboundary = sclib::noType, bool uniteTypes = false, bool uniteTypesNot = false);
    long storeSignal(const char* fileName, SV_Data* pSamples, unsigned short int col, SC_Signal *pSigPar);
    long storeSignal(const char* fileName, SV_Data* pSamples, SC_Signal *pSigPar, int row = -1);
    long storeSignal(const char* fileName, short* &signal, unsigned long int signalLength, SC_Signal *pSigPar);

    //====================================================================================================================
    //	just wrappers around the corresponding functions in SC_Signal, so that this class needs not to be exported
    //  (gives compiler warnings because the base-class from SV_Lib is not exported... ugly, but works ;-)
    //====================================================================================================================
    long saveSignal(const char* fileName, SC_Signal *pSignal, bool useDebugDir = true);
    bool implantSamples(SC_Signal* pSignal, unsigned long startSample, short* samples, unsigned long length);
		void exchangeSignalBuffer(SC_Signal *pSignal, short *newSignak, unsigned long int newSignalLength);

    //====================================================================================================================
    //	return a new signal object with the given signal resampled to the newSampleRate, with possibly done
    //  lowpass-filtering to avoid aliasing; if no resampling is necessary (because current and wished sampleRate already 
		//  are the same, a copy of pSignal is returned)
    //====================================================================================================================
    SC_Signal* resample(SC_Signal *pSignal, double newSampleRate, bool doFiltering = true);
};

#endif
