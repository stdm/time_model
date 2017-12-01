/**************************************************************************/
/*    Derived from:																												*/
/*      - SV_Signal to load data with RIFF WAVE header (Microsoft Windows */
/*        Audio Waveform format, PCM)																			*/
/*																																				*/
/*    Responsibility:																											*/
/*      - analyses MS-Windows WAV filetype (RIFF WAVE) header and loads		*/
/*        samples according to it																					*/
/*																																				*/
/*    Limitations:																												*/
/*      - only PCM format is supported																		*/
/*      - only 16 bit bitrate, mono channel is supported									*/
/*			- no byte-order-swapping at the moment, so only use on intel			*/
/*				machines!																												*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.02.2004																								*/
/**************************************************************************/

#ifndef __SC_Signal_WAVE_H__
#define __SC_Signal_WAVE_H__

#include "SC_Signal.h"
#include <SV_Error.h>

#define SCLIB_WAVE_HEADER_LENGTH 44

class SC_Signal_WAVE : public SC_Signal {
	private :
		
		typedef SV_Int32	DWORD;
		typedef SV_Int16	WORD;
		typedef SV_Int8		BYTE;
		typedef struct {
			 DWORD  fileSize;
			 WORD	  format;         /* 1 -> sclib::encodingPCM													 */
			 WORD   channelCount;   /* 1 (mono) or 2 (stereo, not supported)						 */
			 DWORD  sampleRate;  		/* 44100,. 22050, 11025, 8000, ...									 */
			 WORD   bitsPerSample;  /* 8 (not supported) or 16													 */
			 DWORD  dataLength;     /* length of DATA-chunk														   */
		} SC_HeaderTypeWAVE;

	protected :
		
    SC_Signal_WAVE::SC_HeaderTypeWAVE header;
		SV_DataIO::SV_DatatypeSizes sizeInFile, sizeInCode;

		void createWavFileSizes(SV_DataIO::SV_DatatypeSizes &sizes);
    int fillHeaderInfo(const char *FName);
		int Load_WAVE(const char *FName, short *Speech, unsigned long int start, unsigned long int end);
		
		long storeHeader(const char* fileName, int numSamples = 0);
		long storeSamples(const char* filename, short* pSamples, int seekTo = SCLIB_WAVE_HEADER_LENGTH);
		void makeLittleEndian(char* input, char* output, int length);

	public :
		
    //====================================================================================================================
    //  fill the header; fileSize and dataSize are only set if != 0
    //====================================================================================================================
    void setHeader(unsigned long fileSize = 0, unsigned short format = 1, unsigned short channelCount = 1, unsigned long sampleRate = 16000, unsigned short bitsPerSample = 16, unsigned long dataLength = 0);

		SC_Signal_WAVE();
    SC_Signal_WAVE(SC_Signal_WAVE* oldSignal);
    SC_Signal_WAVE(const char* fileName);
		virtual ~SC_Signal_WAVE();
		
 		long LoadSignal(const char *FName);
    long LoadSignal(const char *fileName, unsigned long int start, unsigned long int end);  
    long SaveSignal(const char *FName);
};

#endif
