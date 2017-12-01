/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Signal_jWAVE to load wav files via Java-streams							*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.01.2009																								*/
/**************************************************************************/

#ifndef __SC_Signal_jWAVE_H__
#define __SC_Signal_jWAVE_H__

#include "SC_Signal.h"
#include "SC_JStreamReader.h"
#include <SV_Error.h>

class SC_Signal_jWAVE : public SC_Signal {
	private :
		const unsigned int headerLength;
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
		
#ifdef SC_USE_JNI
		SC_JStreamReader jio;
#endif
    SC_Signal_jWAVE::SC_HeaderTypeWAVE header;
		SV_DataIO::SV_DatatypeSizes sizeInFile, sizeInCode;

		void createWavFileSizes(SV_DataIO::SV_DatatypeSizes &sizes);
    int fillHeaderInfo();
		int Load_WAVE(short *Speech, unsigned long int start, unsigned long int end);
		
	public :
		
    //====================================================================================================================
    //  fill the header; fileSize and dataSize are only set if != 0
    //====================================================================================================================
    void setHeader(unsigned long fileSize = 0, unsigned short format = 1, unsigned short channelCount = 1, unsigned long sampleRate = 16000, unsigned short bitsPerSample = 16, unsigned long dataLength = 0);

#ifdef SC_USE_JNI
		//====================================================================================================================
		// sets the stream object from which to read
		//====================================================================================================================
		void setStream(JNIEnv *env, jobject jStreamObject);
#endif

		SC_Signal_jWAVE();
    SC_Signal_jWAVE(SC_Signal_jWAVE* oldSignal);
#ifdef SC_USE_JNI
    SC_Signal_jWAVE(JNIEnv *env, jobject jStreamObject);
#endif
		virtual ~SC_Signal_jWAVE();
		
 		long LoadSignal(const char *FName);
    long LoadSignal(const char *fileName, unsigned long int start, unsigned long int end);  
    long SaveSignal(const char *FName);
};

#endif
