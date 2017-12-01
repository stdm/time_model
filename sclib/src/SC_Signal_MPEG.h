/**************************************************************************/
/*    Responisibility:																									  */
/*      - provide access to MPEG (video/MP2 audio) files and extract the  */
/*        audio-samples out of it                                         */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.04.2005																								*/
/**************************************************************************/

#ifndef __SC_Signal_MPEG_H__
#define __SC_Signal_MPEG_H__

#include "SC_Aux.h"
#include "SC_Signal.h"
#include "SC_Resampling.h"

#include "SC_Api.h"
#ifdef SC_USE_FFMPEG
extern "C" { //see http://lists.mplayerhq.hu/pipermail/ffmpeg-user/2007-June/009255.html for explanation
	#include <avcodec.h>
	#include <avformat.h>
}
#endif

class SC_Signal_MPEG : public SC_Signal {
	private :

	protected :

		typedef struct { //header information (filled by readHeader())
			unsigned short int originalChannelCount; //channels in the MPEG file
			unsigned short int originalSampleRate; //sample rate in the MPEG file
			unsigned short int channelCount; //channels as returned by LoadSignal(); may be different from originalChannelCount due to forceChannels parameter
			unsigned short int sampleRate; //sample rate as returned by LoadSignal(); may be different form originalSampleRate due to forceSampleRate parameter
			enum CodecID codecID; //avlib-internal codec id for the needed decoder (mp2, mp3, ...)
			int streamIndex; //index of the (first found) audio stream we want to decode in the MPEG file (that maybe contains more than one stream, audio and video for example)
			double videoFrameRate; //frame rate of the first found video-stream in the MPEG file (if any)
			int bitsPerSample; //bits per audiosample
			long int startTimeSTB; //start time (in [streamTimeBase] units) of the audio-stream
			double streamTimeBase; //time-base (in [s]) of the audio-stream
			double overallTimeBase; //time-base (in [s]) of the AVLib
		} SC_HeaderTypeMPEG;

#ifdef SC_USE_FFMPEG
    typedef struct { //parameters needed to store all information of the resampler
      ReSampleContext *pResampler;
			SC_Resampling *pHQResampler; //to do resampling (but not channel conversion) with higher quality than with the one provided by in ffmpeg
      int inChannels;
      int outChannels;
      int inSampleRate;
      int outSampleRate;
    } ResampleContext;
#endif

		static double lastPktStartTimeStamp; //timestamp (in [s]) of the beginning of the wanted (1st) audio-stream in the last decoded packet
		static double lastPktDuration; //duration (in[s]) of the audio-frame in the last decoded audio-packet
		static unsigned long int lastPktFirstSample; //first sample number (counting the first overall sample in the file as number 0) in the last decoded packet
		static unsigned long int lastPktNumber; //overall number of the last decoded packet (starting from 0)
		static char lastFileName[sclib::bufferSize]; //the file to which the lastPkt* static members refer

	  bool verbose; //make more debug output on the console if true
    bool avLibInitialized; //is the avLib already initialized? (used by initAVlib())
		bool hqResampling; //should a higher quality resampler be used?
		bool fastSeeking; //enables/diables seeking for the LoadSignal(..., start, end) method (doesn't work well at the moment)
		int forceSampleRate; //to convert whatever is included in the MPEG file to this samplerate and channelcount; 0 if nothing should be done
		int forceChannels; //must always be 1!

		//in these 3 variables samples are stored that are already decoded/resampled, but are no further needed by LoadSignal(..); when re-entering this function later and trying to read the directly following segment, they can be used so that seeking is not needed
		short *overreadBuffer; //buffer for this read-too-much samples
	  unsigned int overreadLength; //number of samples in that buffer
		unsigned long overreadBufferStartPos; //sample-number (in the whole file) of the first sample in the buffer, in output-samplerate units

#ifdef SC_USE_FFMPEG
    AVCodecContext *pDecoder; //pointer to the needed audio decoder (filled by initDecoder())
		AVFormatContext *pFormat; //pointer to the format-context of the loaded file (filled by readHeader())
    SC_Signal_MPEG::ResampleContext *pResampler; //pointer to a resampler & it's parameters (filled by initResampler())
#endif
		SC_Signal_MPEG::SC_HeaderTypeMPEG header; //file header information needed by LoadSignal() (filled by readHeader())

    //====================================================================================================================
    // initializes the AV-Lib codecs, if not already done
    //====================================================================================================================
#ifdef SC_USE_FFMPEG
    void initAVlib(void);
#endif

    //====================================================================================================================
    // initialize the decoder, if not already done
    //====================================================================================================================
#ifdef SC_USE_FFMPEG
    bool initDecoder(void);
#endif

    //====================================================================================================================
    // initialize the resampler, if not already done
    //====================================================================================================================
#ifdef SC_USE_FFMPEG
    bool initResampler(int inChannels, int outChannels, int inSampleRate, int outSampleRate, bool forceReload = false);
#endif

		//====================================================================================================================
		//  Resample the signal to another sampling rate and/or channel-count; returns the length of the output-buffer
		//  the resampler needs to be initialized beforehand by calling initResamper(); outSignal is to be allocated big 
		//  enough beforehand, the number of actually written samples is returned; sampleCount means the overall number of 
		//  samples at different time instants (invariant to more or less channels). It is assmumed that inSignal is of size 
		//  SCLIB_MPEG_INBUF_SIZE or a little greater for non-hq resampling. Set framePos==SCLIB_MODE_LASTFRAME when no more 
		//  audio is to be resampled so the cache is emptied and all samples are returned; set it to SCLIB_MODE_FIRSTFRAME for
		//  the first call on a new separate signal-part to intialize the cache of the non-hq resampler; set it to
		//  SCLIB_MODE_INTERMEDIATEFRAME in all other cases
		//====================================================================================================================
#ifdef SC_USE_FFMPEG
		int resample(short* inSignal, int sampleCount, short* outSignal, int framePos);
#endif

		//====================================================================================================================
		//  Decodes a single (the next) audio frame; readHeader() and initDecoder() must be called before
		//  Based on the code by M. Boehme to be found at http://www.inb.uni-luebeck.de/~boehme/using_libavcodec.html
		//  Sets lastFrame=true, when the end of the file is reached
		//  If onlyFreePacket==true, nopthing is done except freeing the static AVPacket variable
		//====================================================================================================================
#ifdef SC_USE_FFMPEG
		bool getNextAudioFrame(short *samples, int &len, bool &lastFrame, bool reset = false, bool onlyFreePacket = false);
#endif

		//====================================================================================================================
		//  Seeks to the last known sample/timestamp/packet correspondence (the last decoded packet by this class), if the
		//  wanted startSample is greater or equal as the last known startSample; this should be better as scanning the file 
		//  from the beginning. Returns true on a successefully seek,otherwise (failure or no seek) returns false. Seeking can 
		//  only be performed if the high quality resampler is used 'cause the other one produces unpredictable sample-numbers
		//====================================================================================================================
#ifdef SC_USE_FFMPEG
		bool seekNearer2sample(unsigned long int startSample);
#endif

	public :
   
		SC_Signal_MPEG(unsigned short int signalType, bool verbose = true, int forceSampleRate = 16000, int forceChannels = 1, bool hqResampling = true, bool fastSeeking = false);

    //====================================================================================================================
    // this is a copy-constructor, which copys the header-info of a given class
    //====================================================================================================================
    SC_Signal_MPEG(SC_Signal_MPEG* oldSignal);

    //====================================================================================================================
    // this constructor analyses the header without loading samples into memory
    //====================================================================================================================
    SC_Signal_MPEG(const char* fileName, unsigned short int signalType, bool verbose = true, int forceSampleRate = 16000, int forceChannels = 1, bool hqResampling = true, bool fastSeeking = false);
		
    //====================================================================================================================
    //  the destructor
    //====================================================================================================================
 		virtual ~SC_Signal_MPEG();
		
    //====================================================================================================================
    //  get & set
    //====================================================================================================================
    void setSignalType(unsigned short int newType) {this->signalType = newType; return;}
		SC_Signal_MPEG::SC_HeaderTypeMPEG* getHeader(void) {return &(this->header);}

		//====================================================================================================================
		//  extract info out of the given media-file and fill class members: fileName, channelCount, sampleRate, sampleCount,
		//  bitsPerSample, videFrameRate
		//====================================================================================================================
    bool readHeader(const char *fileName);

    //====================================================================================================================
    //  Load comlete signal into memory              
    //====================================================================================================================
 		virtual long LoadSignal(const char *FName);

    //====================================================================================================================
		//  Load signal from start-sample to end-sample into memory, assuming that the header is already read
		//  The given sample-boundaries 'start' and 'end' are expected to be given in terms of the output channel-copunt/
		//  sample-rate, not in terms of those originally found in the MPEG file (before possibly resampling)
    //====================================================================================================================
    virtual long LoadSignal(const char *fileName, unsigned long int start, unsigned long int end);  

    //====================================================================================================================
    //  Generate a new wav (!!!)-file with the header of this (adpted to the new data) and the samples stored in Buf_L
    //  because we don't want to encode a new MPEG file here, we just create a new WAV file instead
    //====================================================================================================================
    virtual long SaveSignal(const char *FName);
};

#endif
