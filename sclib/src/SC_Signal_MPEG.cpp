/**************************************************************************/
/*    Responisibility:																									  */
/*      - provide access to MPEG (video/MP2 audio) files and extract the  */
/*        audio-samples out of it                                         */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.04.2005																								*/
/**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <fstream>
#include "SC_Aux.h"
#include "SC_Signal_MPEG.h"
#include "SC_Signal_WAVE.h"

#define SCLIB_MPEG_INBUF_SIZE					4096

#define SCLIB_MODE_FIRSTFRAME					0
#define SCLIB_MODE_INTERMEDIATEFRAME	1
#define SCLIB_MODE_LASTFRAME					2

//====================================================================================================================
// definition of static members
//====================================================================================================================
double SC_Signal_MPEG::lastPktStartTimeStamp = 0.0;
double SC_Signal_MPEG::lastPktDuration = 0.0;
unsigned long int SC_Signal_MPEG::lastPktFirstSample = 0;
unsigned long int SC_Signal_MPEG::lastPktNumber = 0;
char SC_Signal_MPEG::lastFileName[sclib::bufferSize] = "\0";

//====================================================================================================================
// default constructor
//====================================================================================================================
SC_Signal_MPEG::SC_Signal_MPEG(unsigned short int signalType, bool verbose, int forceSampleRate, int forceChannels, bool hqResampling, bool fastSeeking) : SC_Signal() {
#ifdef SC_USE_FFMPEG
	this->pResampler = NULL;
  this->pDecoder = NULL;
	this->pFormat = NULL;
  this->avLibInitialized = false;
  this->signalType = signalType;
  this->verbose = verbose;
	this->forceSampleRate = forceSampleRate;
	this->forceChannels = forceChannels;
	this->hqResampling = hqResampling;
	this->fastSeeking = fastSeeking;

	this->header.bitsPerSample = 0;
	this->header.channelCount = 0;
	this->header.codecID = CODEC_ID_NONE;
	this->header.originalChannelCount = 0;
	this->header.originalSampleRate = 0;
	this->header.sampleRate = 0;
	this->header.streamIndex = 0;
	this->header.videoFrameRate = 0;

	this->overreadBuffer = NULL;
	this->overreadLength = 0;
	this->overreadBufferStartPos = 0;

	if (this->forceChannels != 1) {
		printf("Number of output channels must be 1");
		this->forceChannels = 1;
		this->failure = true;
	}
#else
	REPORT_ERROR(SVLIB_Fail, "ffmpeg not available");
#endif
}

//====================================================================================================================
// this constructor analyses the header without loading samples into memory
//====================================================================================================================
SC_Signal_MPEG::SC_Signal_MPEG(const char* fileName, unsigned short int signalType, bool verbose, int forceSampleRate, int forceChannels, bool hqResampling, bool fastSeeking) : SC_Signal() {
#ifdef SC_USE_FFMPEG
	this->pResampler = NULL;
  this->pDecoder = NULL;
	this->pFormat = NULL;
  this->avLibInitialized = false;
  this->signalType = signalType;
  this->verbose = verbose;
	this->forceSampleRate = forceSampleRate;
	this->forceChannels = forceChannels;
	this->hqResampling = hqResampling;
	this->fastSeeking = fastSeeking;

	this->header.bitsPerSample = 0;
	this->header.channelCount = 0;
	this->header.codecID = CODEC_ID_NONE;
	this->header.originalChannelCount = 0;
	this->header.originalSampleRate = 0;
	this->header.sampleRate = 0;
	this->header.streamIndex = 0;
	this->header.videoFrameRate = 0;

	this->overreadBuffer = NULL;
	this->overreadLength = 0;
	this->overreadBufferStartPos = 0;

	if (this->forceChannels != 1) {
		printf("Number of output channels must be 1");
		this->forceChannels = 1;
		this->failure = true;
	}

	if (readHeader(fileName) == false) {
		this->failure = true;
	}
#else
	REPORT_ERROR(SVLIB_Fail, "ffmpeg not available");
#endif
}

//====================================================================================================================
// this is a copy-constructor, which copys the header-info of a given class
//====================================================================================================================
SC_Signal_MPEG::SC_Signal_MPEG(SC_Signal_MPEG* oldSignal) : SC_Signal(oldSignal) {
#ifdef SC_USE_FFMPEG
	this->pResampler = NULL;
  this->pDecoder = NULL;
	this->pFormat = NULL; //can't be copied so easily... so, below we call readHeader() to accomplish this
  this->avLibInitialized = false;

	this->verbose = oldSignal->verbose;
	this->forceChannels = oldSignal->forceChannels;
	this->forceSampleRate = oldSignal->forceSampleRate;
	this->hqResampling = oldSignal->hqResampling;
	this->fastSeeking = oldSignal->fastSeeking;

	this->header = oldSignal->header;

	if (oldSignal->overreadBuffer != NULL && oldSignal->overreadLength > 0) {
		MArray_1D(this->overreadBuffer, oldSignal->overreadLength, short, "SC_Signal_MPEG: overreadBuffer");
	} else {
		this->overreadBuffer = NULL;
	}
	this->overreadLength = oldSignal->overreadLength;
	this->overreadBufferStartPos = oldSignal->overreadBufferStartPos;

	readHeader(this->fileName);
#else
	REPORT_ERROR(SVLIB_Fail, "ffmpeg not available");
#endif
}

//====================================================================================================================
// destructor 
//====================================================================================================================
SC_Signal_MPEG::~SC_Signal_MPEG() {
#ifdef SC_USE_FFMPEG
  if (this->pResampler != NULL) {
		if (this->pResampler->pResampler != NULL) {
			audio_resample_close(this->pResampler->pResampler);
		}
		MFree_0D(this->pResampler->pHQResampler);
    MFree_0D(this->pResampler);
  }
  if (this->pDecoder != NULL) {
    avcodec_close(this->pDecoder);
  }
	if (this->pFormat != NULL) {
		av_close_input_file(this->pFormat);
	}
	MFree_1D(this->overreadBuffer);
#endif
}

//====================================================================================================================
// initializes the AV-Lib codecs
//====================================================================================================================
#ifdef SC_USE_FFMPEG
void SC_Signal_MPEG::initAVlib(void) {
  if (this->avLibInitialized == false) {
    av_register_all();
    this->avLibInitialized = true;
  }

  return;
}
#endif

//====================================================================================================================
// initialize the decoder, if not already done; assumes that the header is already filled by readHeader()
//====================================================================================================================
#ifdef SC_USE_FFMPEG
bool SC_Signal_MPEG::initDecoder(void) {
  bool res = true;
  AVCodec *codec = NULL;

  //close the current decoder if it is for a different format
  if (this->pDecoder != NULL && this->pDecoder->codec_id != this->header.codecID) {
    avcodec_close(this->pDecoder);
    this->pDecoder = NULL;
  }

  if (this->pDecoder == NULL) {
    initAVlib();

    codec = avcodec_find_decoder(this->header.codecID);

    if (codec != NULL) {
      this->pDecoder = avcodec_alloc_context();
			if(codec->capabilities & CODEC_CAP_TRUNCATED) { //Inform the codec that we can handle truncated bitstreams -- i.e., bitstreams where frame boundaries can fall in the middle of packets
				this->pDecoder->flags|=CODEC_FLAG_TRUNCATED;
			}
      if (avcodec_open(this->pDecoder, codec) < 0) {
        printf(" <avcodec_open() returned error> ");
        MFree_0D(this->pDecoder);
      }
    } else {
      printf(" <avcodec_find_decoder() returned NULL> ");
    }

    res = (this->pDecoder != NULL) ? true : false;
  }

  return res;
}
#endif

//====================================================================================================================
// initialize the resampler, if not already done
//====================================================================================================================
#ifdef SC_USE_FFMPEG
bool SC_Signal_MPEG::initResampler(int inChannels, int outChannels, int inSampleRate, int outSampleRate, bool forceReload) {  
  bool res = true;

  initAVlib();

	//force re-initialization of the resampler to get rid of it's cache
  if (forceReload == true && this->pResampler != NULL) {
		if (this->pResampler->pResampler != NULL) {
			audio_resample_close(this->pResampler->pResampler);
		}
		MFree_0D(this->pResampler->pHQResampler);
    MFree_0D(this->pResampler);
  }

  //create space for the resampler & initialize it's parameters
  if (this->pResampler == NULL) {
    this->pResampler = new SC_Signal_MPEG::ResampleContext;
    this->pResampler->pResampler = NULL;
		this->pResampler->pHQResampler = (this->hqResampling == true) ? new SC_Resampling() : NULL;
  } 

	//set the parameters that are needed by both resamplers
  this->pResampler->inChannels = inChannels;
  this->pResampler->outChannels = outChannels;
  this->pResampler->inSampleRate = inSampleRate;
  this->pResampler->outSampleRate = outSampleRate;

  //create a new ffmpeg-resampler if it is needed and the previous parameters don't fit or there is no previous one
  if (this->hqResampling == false && (this->pResampler->pResampler == NULL || this->pResampler->inChannels != inChannels || this->pResampler->outChannels != outChannels || this->pResampler->inSampleRate != inSampleRate || this->pResampler->outSampleRate != outSampleRate)) {
		if (this->pResampler->pResampler != NULL) { //close old resampler, if existant
			audio_resample_close(this->pResampler->pResampler);
		}

		this->pResampler->pResampler = audio_resample_init(outChannels, inChannels, ((this->hqResampling == true && this->header.channelCount == 1) ? inSampleRate : outSampleRate), inSampleRate); //use this resampler only for channel conversion if hq resampling is wanted; don't use hq-resampler if more than 1 channel is wanted
    if (this->pResampler->pResampler == NULL) {
      MFree_0D(this->pResampler->pHQResampler);
			MFree_0D(this->pResampler);
      res = false;
    }
  }

  return res;
}
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
int SC_Signal_MPEG::resample(short* inSignal, int sampleCount, short* outSignal, int framePos) {
  int res = SVLIB_Fail, channelConversedLength;
	unsigned long int resampledLength;
	short *resampledSignal = NULL, *channelConversedSignal;
	double tmp;
	static int lag = 0;
  
  if (this->pResampler != NULL) {
		if (this->hqResampling == false || this->header.channelCount > 1) { //low quality resampling
			//audio_resample() is buggy (or implements a cache as libsamplerate but does not advertise it): sometimes, but always at first call, returns 7-8 samples too less
			//so, here comes a hack to overcomes the problem by artificially enlarging the input-signal with padding zeros for the last frame as suggested by Michel Bardiaux on the ffmpeg user mailing list
			if (framePos == SCLIB_MODE_FIRSTFRAME) {
				lag = 0;
			}
			if (framePos == SCLIB_MODE_FIRSTFRAME || framePos == SCLIB_MODE_INTERMEDIATEFRAME) { //just do normal...
		    res = audio_resample(this->pResampler->pResampler, outSignal, inSignal, sampleCount);
				tmp = (double)(sampleCount) * ((double)(this->header.sampleRate) / (double)(this->header.originalSampleRate));
				lag += sclib::round(tmp) - res; //count missing samples
			} else { //last frame => pad with zeros and apply hack
				res = sampleCount*this->header.originalChannelCount; //nr of entries filled with real data in the input buffer
				if (res < SCLIB_MPEG_INBUF_SIZE) { //don't write in someone else's memory...
					memset(inSignal + res, 0, (SCLIB_MPEG_INBUF_SIZE-res)*sizeof(short)); //set everything in the input-buffer after the real data to 0
				}
				res = audio_resample(this->pResampler->pResampler, outSignal, inSignal, sclib::max(sampleCount, SCLIB_MPEG_INBUF_SIZE/this->header.originalChannelCount));
				tmp = (double)(sampleCount) * ((double)(this->header.sampleRate) / (double)(this->header.originalSampleRate));
				res = sclib::round(tmp) + lag;
			}
		} else { //high quality resampling
			if (this->pResampler->inChannels != this->pResampler->outChannels) { //this is a simple all-channels-to-one-by-averaging channel down-converter
				pResampler->pHQResampler->toMono(inSignal, outSignal, sampleCount, this->header.originalChannelCount);
				channelConversedLength = sampleCount;
				channelConversedSignal = outSignal;
			} else {
				channelConversedLength = sampleCount;
				channelConversedSignal = inSignal;
			}

			if (this->pResampler->inSampleRate != this->pResampler->outSampleRate) { //use pResampler->HQResampler for sample rate conversion afterwards
				resampledSignal = this->pResampler->pHQResampler->resample(channelConversedSignal, channelConversedLength, this->pResampler->inSampleRate, this->pResampler->outSampleRate, resampledLength, (channelConversedSignal == inSignal) ? outSignal : NULL, (framePos == SCLIB_MODE_LASTFRAME) ? true : false); //write directly to outSignal if no channel-conversion was done
			} else {
				resampledLength = channelConversedLength;
				resampledSignal = channelConversedSignal;
			}

			if (outSignal != resampledSignal) { //write the final result if it has been changed
				for (int x = 0; x < (int)(resampledLength); x++) {
					outSignal[x] = resampledSignal[x];
				}
				MFree_1D(resampledSignal);
			}

			res = resampledLength;
		}
  }
  
  return res;
}
#endif

//====================================================================================================================
//  extract info out of the given media-file and fill class members: fileName, channelCount, sampleRate, sampleCount,
//  bitsPerSample, videFrameRate (all in the header), pFormat
//====================================================================================================================
bool SC_Signal_MPEG::readHeader(const char *fileName) {
#ifdef SC_USE_FFMPEG
	AVCodecContext *enc = NULL;
  AVFormatContext *ic = NULL;
	AVInputFormat *file_iformat = NULL;
  AVFormatParameters *ap = NULL;
  int err, i, ret;
	bool foundAudio = false, foundVideo = false, fileChange = false;

	this->failure = false;

	if (this->fileName == NULL || strcmp(this->fileName, fileName) != 0) {
		MFree_1D(this->fileName);
		MArray_1D(this->fileName, strlen(fileName) + 1, char, "SC_Signal_MPEG.readHeader: fileName");
		sprintf(this->fileName, "%s\0", fileName);
	}
	
	if (strcmp(this->fileName, SC_Signal_MPEG::lastFileName) != 0) {
		//empty cache on file change
		SC_Signal_MPEG::lastPktFirstSample = 0;
		SC_Signal_MPEG::lastPktNumber = 0;
		SC_Signal_MPEG::lastPktStartTimeStamp = 0.0;
		SC_Signal_MPEG::lastPktDuration = 0.0;
		sprintf(SC_Signal_MPEG::lastFileName, "%s", this->fileName);
		fileChange = true;

		MFree_1D(this->overreadBuffer);
		this->overreadLength = 0;
		this->overreadBufferStartPos = 0;
	}

	if (this->pFormat != NULL) {
		av_close_input_file(this->pFormat);
	}

  initAVlib();

  // open the input file with generic libav function 
  err = av_open_input_file(&ic, fileName, file_iformat, 0, ap);
	if (err < 0) { //error codes may be looked up in errno.h
		printf("Open MPEG file failed: %d", err);
		this->failure = true;
		return false;
  }
	this->pFormat = ic;

  // If not enough info to get the stream parameters, we decode the
  // first frames to get it. (used in mpeg case for example) 
  ret = av_find_stream_info(ic);
  if (ret < 0) {
		printf("Open MPEG file: Could not find codec parameters");
		this->failure = true;
    return false;
  }

	if (this->verbose == true && fileChange == true) { //only print ffmpeg header info for previously unseen files
		printf("\n");
		dump_format(ic, 0, fileName, false);
		printf("\n");
	}

	//go trough the stream and pick the necessary information from the first found streams
	for (i = 0; i < (int)(ic->nb_streams); i++) {
    enc = ic->streams[i]->codec;
		if (enc->codec_type == CODEC_TYPE_AUDIO && foundAudio == false) {
			this->header.channelCount = (this->forceChannels > 0) ? this->forceChannels : enc->channels;
      this->header.sampleRate = (this->forceSampleRate > 0) ? this->forceSampleRate : enc->sample_rate;
			this->header.originalChannelCount = enc->channels;
			this->header.originalSampleRate = enc->sample_rate;
			this->header.codecID = enc->codec_id;
			this->header.streamIndex = i;
			this->header.streamTimeBase = av_q2d(this->pFormat->streams[i]->time_base);
			this->header.overallTimeBase = 1.0 / (double)(AV_TIME_BASE);
			this->header.startTimeSTB = (long int)(this->pFormat->streams[i]->start_time); //in contrast to what is stated in ffmpeg's rare comments, this seems to be in streamTimeBase units rather than AV_TIME_BASE units

			//api header
			this->SigPar.NChannel = this->header.channelCount;
			this->SigPar.SRate = this->header.sampleRate;
			this->SigPar.Encode = 1;
			this->SigPar.StByte = 0;

			switch (enc->sample_fmt) {
				case SAMPLE_FMT_S16: 
					this->header.bitsPerSample = 16;
					break;
				case SAMPLE_FMT_S24: 
					this->header.bitsPerSample = 24;
					break;
				case SAMPLE_FMT_S32: 
					this->header.bitsPerSample = 32;
					break;
				case SAMPLE_FMT_U8: 
					this->header.bitsPerSample = 8;
					break;
				default: 
					this->header.bitsPerSample = 0;
					break;
			}
			foundAudio = true; //use the first found audio-stream
		} else if (foundVideo == false) {
			this->header.videoFrameRate = (ic->streams[i]->r_frame_rate.den && ic->streams[i]->r_frame_rate.num) ? av_q2d(ic->streams[i]->r_frame_rate) : 1.0/av_q2d(enc->time_base);
			foundVideo = true; //use the first found video-stream
		}
  }
  
  this->sampleCount = (unsigned long)((double)(ic->duration) * this->header.overallTimeBase * this->header.sampleRate); //first pair gives seconds

	return true;
#else
	return false;
#endif
}

//====================================================================================================================
//  Decodes a single (the next) audio frame; readHeader() and initDecoder() must be called before
//  Based on the code by M. Boehme to be found at http://www.inb.uni-luebeck.de/~boehme/using_libavcodec.html
//  Sets lastFrame=true, when the end of the file is reached
//  If onlyFreePacket==true, nothing is done except freeing the static AVPacket variable
//====================================================================================================================
#ifdef SC_USE_FFMPEG
bool SC_Signal_MPEG::getNextAudioFrame(short *samples, int &len, bool &lastFrame, bool reset, bool onlyFreePacket) {
	static AVPacket packet;
	static int bytesRemaining = 0;
	static uint8_t *rawData;
	static bool firstTime = true;
	int bytesDecoded;
	bool endOfFileReached = false;
	double dts, duration; 

	//reset if wanted (used if first called after recreation of the SC_Signal_MPEG object, because the static variables [most urgently bytesRemaining] don't get deleted if an object was created earlier)
	if (reset == true || onlyFreePacket == true) {
		bytesRemaining = 0;
		if (packet.data != NULL) {
			av_free_packet(&packet);	
		}
		firstTime = true;
	}

	if (onlyFreePacket == true) {
		return true;
	}

	if (this->pDecoder == NULL || this->pFormat == NULL) { //this must be initilaized by calling readHeader() and initDecoder() before
		return false;
	}

	if (firstTime) { 	//First time we're called, set packet.data to NULL to indicate it doesn't have to be freed
		firstTime = false;
		packet.data = NULL;
	}

	lastFrame = false;
	while (true) { //Decode packets until we have decoded a complete frame
		while (bytesRemaining > 0) { //Work on the current packet until we have decoded all of it
			//bytesDecoded = avcodec_decode_video(this->pDecoder, pFrame, &frameFinished, rawData, bytesRemaining);
			bytesDecoded = avcodec_decode_audio(this->pDecoder, samples, &len, rawData, bytesRemaining); //Decode the next chunk of data
			//bytesDecoded = avcodec_decode_audio2(this->pDecoder, samples, &thirdParam, rawData, bytesRemaining); //Decode the next chunk of data
			len /= sizeof(short); //convert #bytes -> #values

			if (bytesDecoded < 0)	{ //Was there an error?
				printf("Error while decoding audio frame\n");
				this->failure = true;
				return false;
			}

			bytesRemaining -= bytesDecoded;
			rawData += bytesDecoded;

			if (len > 0) { //Did we finish the current frame? Then we can return
				return true;
			}
		}

		do { //Read the next packet, skipping all packets that aren't for this stream
			if (packet.data!=NULL) { //Free old packet
				av_free_packet(&packet);
			}

			if (av_read_frame(this->pFormat, &packet) < 0) { //Read new packet
				lastFrame = true;
				endOfFileReached = true;
				break;
			}
		} while (packet.stream_index != this->header.streamIndex);

		//cache the timestamp/samplecount of the last decoded packet for faster random access to samples if successively called
		//packet's dts, pts and duration are (in contrast to what is stated in the ffmpeg's rare comments, where it tells it's in SV_TIME_BASE units)
		//in streamTimeBase units; pts is sometimes unsafe, so use dts
		dts = (double)(packet.dts - this->header.startTimeSTB);
		duration = (double)(packet.duration);
		//TODO: actually we don't need all of them... just the duration, rest is for debugging reasons
		SC_Signal_MPEG::lastPktNumber =  sclib::round(dts / duration);
		SC_Signal_MPEG::lastPktFirstSample = (SC_Signal_MPEG::lastPktNumber) * sclib::round(duration * this->header.streamTimeBase * this->header.sampleRate); //this doesnt't correspond fully with reality due to caching in the resampler... this value here would be correct if the resampler would not cache the data
		SC_Signal_MPEG::lastPktStartTimeStamp = (dts + this->header.startTimeSTB) * this->header.streamTimeBase;
		SC_Signal_MPEG::lastPktDuration = duration * this->header.streamTimeBase;

		if (endOfFileReached == false) {
			bytesRemaining = packet.size;
			rawData = packet.data;
		} else {
			break;
		}
	}
	
	//control reaches here only on EOF or error while reading packets
	//bytesDecoded=avcodec_decode_video(this->pDecoder, pFrame, &frameFinished, rawData, bytesRemaining);
	bytesDecoded = avcodec_decode_audio(this->pDecoder, samples, &len, rawData, bytesRemaining); //Decode the rest of the last frame
	//bytesDecoded = avcodec_decode_audio2(this->pDecoder, samples, &len, rawData, bytesRemaining); //Decode the rest of the last frame
	len /= sizeof(short); //convert #bytes -> #values

	if (packet.data != NULL) { //Free last packet
		av_free_packet(&packet);
	}

	return (len != 0);
}
#endif

//====================================================================================================================
//  Seeks to the last known sample/timestamp/packet correspondence (the last decoded packet by this class), if the
//  wanted startSample is greater or equal as the last known startSample; this should be better as scanning the file 
//  from the beginning. Returns true on a successefully seek,otherwise (failure or no seek) returns false. Seeking can 
//  only be performed if the high quality resampler is used 'cause the other one produces unpredictable sample-numbers
//====================================================================================================================
#ifdef SC_USE_FFMPEG
bool SC_Signal_MPEG::seekNearer2sample(unsigned long int startSample) {
	int res = SVLIB_Fail;
	//bool tmp; //only needed to give getNextAudioFrame a 3rd parameter...
	int64_t timeStamp = (int64_t)sclib::max(0, ((SC_Signal_MPEG::lastPktStartTimeStamp - SC_Signal_MPEG::lastPktDuration) / this->header.streamTimeBase) - (0.01 / this->header.streamTimeBase)); //start time of the packet before the last decoded packet in the audio-stream's timebase units (decode some previous frames to load the cache of the resampler)
	int64_t wantedTimeStamp = (int64_t)(this->header.startTimeSTB + (startSample / (double)(this->header.sampleRate)) / this->header.streamTimeBase); //timestamp of the startStample in the audio-stream's timebase units

	if (startSample > 0 && wantedTimeStamp >= timeStamp && strcmp(this->fileName, SC_Signal_MPEG::lastFileName) == 0) {
		//Here's an excerpt from a post by danielp on 2006-12-04 in the ffmpeg windows forum at
    //http://arrozcru.no-ip.org/ffmpeg_forum/viewtopic.php?t=227:
		//"From what you say, if I understand it correctly, you go back reading the same 
		// packet by issuing an av_seek_frame and the result you get back has the right 
		// time stamps but incorrect payload. 
		// In browsing the code for mpegps_read_packet, I saw that it increments the underlying 
		// buffer pointer. Since you're reading the same packet structure, I suppose that you 
		// need to reset this pointer. It seems that because you're reading the same packet 
		// twice, some data structure remember the offset to read from and continue at the 
		// point. (Maybe using av_free_packet or some function to reinit the stream will do the trick)."

		//so, here we go: => it does not help!!!
		//getNextAudioFrame(NULL, res, tmp, false, true); //free the internal AV_Packet variable as suggested 

		//now do the seek
		res = av_seek_frame(this->pFormat, this->header.streamIndex, timeStamp, AVSEEK_FLAG_BACKWARD);
	}

	return (res >= 0);
}
#endif

//====================================================================================================================
//  Load comlete signal into memory              
//====================================================================================================================
long SC_Signal_MPEG::LoadSignal(const char *FName) {
#ifdef SC_USE_FFMPEG
	short *signal, frameSignal[SCLIB_MPEG_INBUF_SIZE + FF_INPUT_BUFFER_PADDING_SIZE], resampledSignal[SCLIB_MPEG_INBUF_SIZE + FF_INPUT_BUFFER_PADDING_SIZE];
	int *length, frameLength, resampledLength, actualLength, framePos;
	unsigned long int overallSampleCount = 0, predictedSampleCount;
	bool firstTime = true, lastFrame;
	double tmp;

	this->failure = false;

	//reset maybe previously loaded data
	MFree_1D(this->Buf_L);
	this->Len = 0;
	MFree_1D(this->overreadBuffer);
	this->overreadLength = 0;
	this->overreadBufferStartPos = 0;

  //read numerous header data
	if (readHeader(FName) == false) {
		return SVLIB_Fail;
	}

  //instantiate the decoder
  if (initDecoder() == false) {
    return SVLIB_Fail;
  }

	//instantiate the resamper, if needed, and set new sampleRate/channelCount in the header
	if (this->header.originalChannelCount != this->header.channelCount || this->header.originalSampleRate != this->header.sampleRate) {
		if (initResampler(this->header.originalChannelCount, this->header.channelCount, this->header.originalSampleRate, this->header.sampleRate) == false) {
			return SVLIB_Fail;
		}
	}
	
	//set end of buffer to 0 (this ensures that no overreading happens for damaged mpeg streams)
  memset(frameSignal + SCLIB_MPEG_INBUF_SIZE, 0, FF_INPUT_BUFFER_PADDING_SIZE*sizeof(short)); 

	//which information shall be copied to the final buffer - the oiginal or the resampled one?
	length = (this->pResampler == NULL) ? &frameLength : &resampledLength;
	signal = (this->pResampler == NULL) ? frameSignal : resampledSignal;

	//allocate memory to hold the signal in the external available buffer
	MArray_1D(this->Buf_L, this->sampleCount, short, "SC_Signal_MPEG.LoadSignal: Buf_L");
	this->Len = this->sampleCount;

	//decode and resample all audio in the file
	while(getNextAudioFrame(frameSignal, frameLength, lastFrame, firstTime) == true)	{ //true implies frameLength > 0
		if (this->pResampler != NULL) { //we only have a resampler here (and then we have one for shure) if we need it
			tmp = (frameLength / (double)(this->header.originalChannelCount)) * ((double)(this->header.sampleRate) / (double)(this->header.originalSampleRate));
			predictedSampleCount = sclib::round(tmp);
			if (firstTime == true) { //determine where we stand (necessary for the resampler to handle caching)
				framePos = SCLIB_MODE_FIRSTFRAME;
			} else if (lastFrame == true || overallSampleCount + predictedSampleCount >= this->sampleCount) {
				framePos = SCLIB_MODE_LASTFRAME;
			} else {
				framePos = SCLIB_MODE_INTERMEDIATEFRAME;
			}
			resampledLength = resample(frameSignal, frameLength / this->header.originalChannelCount, resampledSignal, framePos);
			if (resampledLength <= 0) {
				printf("Error while resampling an audio frame");
				MFree_1D(this->Buf_L);
				this->Len = 0;
				this->failure = true;
				return SVLIB_Fail;
			}
		}
		
		//don't decode more samples than was stated earlier while analyzing the file-duration in readHeader()
		//strangely, there are ofton more samples to decode than was pre-computed, but it's typically less than 1 second, so omit these samples
		actualLength = (overallSampleCount + *length <= this->sampleCount) ? *length : this->sampleCount - overallSampleCount; 
		
		for (int x = 0; x < actualLength; x++) { //fill the final buffer
			this->Buf_L[overallSampleCount + x] = signal[x];
		}

		firstTime = false;
		overallSampleCount += actualLength;
		if (overallSampleCount >= this->sampleCount) {
			break; //omit the last few samples if the file is not read to the end but there are already as much samples read as can be concluded from the sample-rate and the file-duration in the header
		}
	}

  return overallSampleCount;
#else
	return SVLIB_Fail;
#endif
}

//====================================================================================================================
//  Load signal from start-sample to end-sample into memory, assuming that the header is already read
//  The given sample-boundaries 'start' and 'end' are expected to be given in terms of the output channel-count/
//  sample-rate, not in terms of those originally found in the MPEG file (before possibly resampling)
//====================================================================================================================
long SC_Signal_MPEG::LoadSignal(const char *fileName, unsigned long int start, unsigned long int end) {
#ifdef SC_USE_FFMPEG
	short *signal, frameSignal[SCLIB_MPEG_INBUF_SIZE + FF_INPUT_BUFFER_PADDING_SIZE], resampledSignal[SCLIB_MPEG_INBUF_SIZE + FF_INPUT_BUFFER_PADDING_SIZE];
	int *length, frameLength, resampledLength, actualLength, signalOffset, framePos, x;
	unsigned long int overallSampleCount = 0, actualSampleCount = 0, predictedSampleCount;
	bool firstTime = true, lastFrame, resetResampler = false;
	double tmp;

	this->failure = false;

	//reset maybe previously loaded data
	MFree_1D(this->Buf_L);
	this->Len = 0;

	if (start > this->sampleCount || end > this->sampleCount || start > end) {
		printf("Given sample-numbers are invalid or exceed the number of samples in the given file");
		this->failure = true;
		return SVLIB_Fail;
	}

  //is the header already read? (don't read it here 'cause it would destruct the pFromat meber and with this the current read position in the stream that can be exploited to overcome the seeking problem)
	if (this->header.codecID == CODEC_ID_NONE) {
		printf(" <no codec> ");
		return SVLIB_Fail;
	}

  //instantiate the decoder
  if (initDecoder() == false) {
    printf(" <no decoder> ");
    return SVLIB_Fail;
  }
	
	//set end of buffer to 0 (this ensures that no overreading happens for damaged mpeg streams)
  memset(frameSignal + SCLIB_MPEG_INBUF_SIZE, 0, FF_INPUT_BUFFER_PADDING_SIZE*sizeof(short)); 

	//allocate memory to hold the signal in the external available buffer
	MArray_1D(this->Buf_L, end-start+1, short, "SC_Signal_MPEG.LoadSignal: Buf_L");
	this->Len = end-start+1;

	//try to seek to the last known position so we can omit decoding the file from the beginning
	if (start == SC_Signal_MPEG::lastPktFirstSample + SC_Signal_MPEG::lastPktDuration*this->header.sampleRate) {
		//where are in perfectly the right state and just need to load the next frame 'cause the last one was consumed entirely and our segment starts directly at the frame to open
		overallSampleCount = start;
	} else if (start == this->overreadBufferStartPos && this->overreadLength > 0 && this->overreadBuffer != NULL) {
		//simple seek strategy: if we try to load a segment which directly comes after the last loaded sample, we don't need seeking and we don't need to decode from start
		for (x = 0; x < (int)(this->overreadLength); x++) {
			this->Buf_L[x] = this->overreadBuffer[x];
		}
		overallSampleCount = this->overreadBufferStartPos + this->overreadLength;
		actualSampleCount = this->overreadLength;
	} else if (this->fastSeeking == true &&	seekNearer2sample(start) == true) {
		//TODO: seeking somehow breaks up frame content... we think we read from the packet's start according to packet-dts 
		//      at the next getNextAudioFrame() call, but we are actually 2 frames ahead (in the simple example i tested) 
		//      without any chance of noticing it from code... update: seems to only mess up the timing in cooperation with resampling
		overallSampleCount = SC_Signal_MPEG::lastPktFirstSample;
		resetResampler = true;
	} else {
		//all seeking failed... read the stream from the beginning by re-opening it
		//if (av_seek_frame(this->pFormat, this->header.streamIndex, 0, AVSEEK_FLAG_BACKWARD) < 0) {
		if (readHeader(fileName) == false) {
			printf(" <reading header failed> ");
			return SVLIB_Fail;
		}
		resetResampler = true;
	}

	//instantiate the resamper, if needed, and set new sampleRate/channelCount in the header
	if (this->header.originalChannelCount != this->header.channelCount || this->header.originalSampleRate != this->header.sampleRate) {
		if (initResampler(this->header.originalChannelCount, this->header.channelCount, this->header.originalSampleRate, this->header.sampleRate, resetResampler) == false) {
			printf(" <no resampler> ");
			return SVLIB_Fail;
		}
	}

	//which information shall be copied to the final buffer - the oiginal or the resampled one?
	length = (this->pResampler == NULL) ? &frameLength : &resampledLength;
	signal = (this->pResampler == NULL) ? frameSignal : resampledSignal;

	//here, we already used the information of the last call, if it was available and we could use it. so free it and recreate it for this call only if it's necessary later
	MFree_1D(this->overreadBuffer);
	this->overreadLength = 0;
	this->overreadBufferStartPos = 0;

	//decode and resample all audio in the file
	while(getNextAudioFrame(frameSignal, frameLength, lastFrame, firstTime) == true)	{ //true implies frameLength > 0
		if (firstTime == true) { //because we maybe didn't seek far enough 
			//overallSampleCount = SC_Signal_MPEG::lastPktFirstSample; //this caused writing more samples than fitted into the Buf_L!!!
		}
		if (this->pResampler != NULL) { //we only have a resampler here (and then we have one for shure) if we need it
			tmp = (frameLength / (double)(this->header.originalChannelCount)) * ((double)(this->header.sampleRate) / (double)(this->header.originalSampleRate));
			predictedSampleCount = sclib::round(tmp);
			if (firstTime == true) { //determine where we stand
				framePos = SCLIB_MODE_FIRSTFRAME;
			} else if (lastFrame == true || overallSampleCount + predictedSampleCount >= end) {
				framePos = SCLIB_MODE_LASTFRAME;
			} else {
				framePos = SCLIB_MODE_INTERMEDIATEFRAME;
			}
			resampledLength = resample(frameSignal, frameLength / this->header.originalChannelCount, resampledSignal, framePos);
			if (resampledLength <= 0) {
				printf("Error while resampling an audio frame");
				MFree_1D(this->Buf_L);
				this->Len = 0;
				this->failure = true;
				return SVLIB_Fail;
			}
		} //resampling
		
		//only copy the data that lies between start and end, omit the rest
		//(we have to decode and resample it anyway to find the correct sample to start)
		if (overallSampleCount + *length >= start) {
			if (overallSampleCount > start) { //find the concrete borders of this frame's useful data
				signalOffset = 0; //the relevant data for this frame starts at it's beginning
				if (overallSampleCount + *length > end) {
					actualLength = end - overallSampleCount + 1;
				} else {
					actualLength = *length;
				}
			} else {
				signalOffset = start - overallSampleCount; //the relevant data for this frame maybe starts somewhere inside the frame, not directly at the beginning
				if (overallSampleCount + *length > end) {
					actualLength = end - start + 1;
				} else {
					actualLength = (overallSampleCount + *length) - start;
				}
			}

			if ((int)(actualSampleCount + actualLength) > this->Len) {
				REPORT_ERROR(SVLIB_BadData, "Read more samples from an MPEG signal than fitted into the buffer!");
			}

			for (x = 0; x < actualLength; x++) { //fill the final buffer
				this->Buf_L[actualSampleCount + x] = signal[signalOffset + x];
			}
			actualSampleCount += actualLength;
		} //loaded if samples reach into the area to load

		firstTime = false;
		overallSampleCount += *length;
		if (overallSampleCount >= end) {
			//fill "cache" of now unused samples for the next call, if it tries to load the really next part beginning with these unused samples
			this->overreadLength = *length - signalOffset - actualLength;
			
				MArray_1D(this->overreadBuffer, this->overreadLength+1, short, "SC_Signal_MPEG.LoadSignal: overreadBuffer");
				for (x = 0; x < (int)(this->overreadLength); x++) {
					this->overreadBuffer[x] = signal[signalOffset + actualLength + x];
				}
				this->overreadBufferStartPos = overallSampleCount - *length + signalOffset + actualLength;
			
			break; //this call is now finished
		}
	} //while decoding

  return overallSampleCount;
#else
	return SVLIB_Fail;
#endif
}

//====================================================================================================================
//  Generate a new wav (!!!) file with samples stored in Buf_L
//  because we don't want to encode a new MPEG file here, we just create a new WAV file instead
//====================================================================================================================
long SC_Signal_MPEG::SaveSignal(const char *FName) {
  long int res = 0;
  SC_Signal_WAVE *pWaveWriter = new SC_Signal_WAVE();

  //create a WAV-file header only with the necessary information to create a new WAV file
  pWaveWriter->setHeader(0, 1, this->SigPar.NChannel, this->SigPar.SRate, 16, 0); //use the SigPar.* parameters here (instead of header.*) that could have been modified from outside (together with resampling)

	SC_Signal_WAVE *pWaveHook = pWaveWriter;
	SC_Signal_MPEG *pHook = this;
	while (pHook != NULL) {
		pWaveHook->setBuf_L(pHook->GetBuf_L(), pHook->GetLen());
		pHook = (SC_Signal_MPEG*)(pHook->Next);
		if (pHook != NULL) {
			pWaveHook->Next = new SC_Signal_WAVE();
			pWaveHook = (SC_Signal_WAVE*)(pWaveHook->Next);
		}
	}	

  //pWaveWriter->setBuf_L(this->Buf_L, this->Len);
  res = pWaveWriter->SaveSignal(FName);
  //pWaveWriter->forgetBuf_L();

	pWaveHook = pWaveWriter;
	while (pWaveHook != NULL) {
		pWaveHook->forgetBuf_L();
		pWaveHook = (SC_Signal_WAVE*)(pWaveHook->Next);
	}

  return res;
}
