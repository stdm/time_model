/**************************************************************************/
/*    Derived from:																												*/
/*      - SV_Signal to permit some enhancements:                          */
/*        - convert short buffer to float                                 */
/*        - give acces to that buffer so that the signal doesn't have be  */
/*          loaded separately by each feature extraction class            */
/*        - write the buffer back to a file                               */
/*        - load specified part of a file                                 */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.04.2005																								*/
/**************************************************************************/

#ifndef __SC_Signal_H__
#define __SC_Signal_H__

#include "SC_TweakableParameters.h" //for SCLIB_ST_* constants
#include <SV_Signal.h>
#include <SV_Error.h>
#include <SV_Data.h>

class SC_Signal : public SV_Signal {
	private :

	protected :

		bool  failure; //true if a failure occured during construction
    char* fileName;
    unsigned long int sampleCount; //overall samplecount in the whole file (as opposed to the amount of loaded samples in "Len")
		unsigned short int signalType;

	public :
		SC_Signal *Next; //to form a linked list
				
		SC_Signal();

    //====================================================================================================================
    //  This is a copy-constructor, which copys the header-info of a given class
    //====================================================================================================================
    SC_Signal(SC_Signal* oldSignal);

    //====================================================================================================================
    //  This constructor analyses the header without loading samples into memory
    //====================================================================================================================
    SC_Signal(const char* fileName);
		
    char* getFileName(void) {return this->fileName;}
    unsigned long int getSampleCount(void) {return this->sampleCount;}
    void setSampleCount(unsigned long int newSampleCount) {this->sampleCount = newSampleCount; return;}

    virtual ~SC_Signal();
		
    //====================================================================================================================
    //  Load comlete signal into memory              
    //====================================================================================================================
 		virtual long LoadSignal(const char *FName) = 0;

    //====================================================================================================================
    //  Load signal from start-sample to end-sample into memory 
    //  this method expects the header already to be filled by the constructor!              
    //====================================================================================================================
    virtual long LoadSignal(const char *fileName, unsigned long int start, unsigned long int end) = 0;  

    //====================================================================================================================
    //  Generate a new audio-file with the header of this (adapted to the new data) and the samples stored in Buf_L
    //====================================================================================================================
    virtual long SaveSignal(const char *FName) = 0;

    //====================================================================================================================
    //  sets Buf_L to the given pointer
    //====================================================================================================================
		void setBuf_L(short* pSamples, unsigned long int length);

    //====================================================================================================================
    //  converts Buf_L to float and returns a SV_Data object containing it
    //====================================================================================================================
		SV_Data* convertBuf_L(void);

    //====================================================================================================================
    //  returns true, if the currently loaded signal can be replaced from start to start+sampleLength with the given 
    //  samples
    //====================================================================================================================
    bool implantSamples(unsigned long int start, short* pSamples, unsigned long int sampleLength);

    //====================================================================================================================
    //  Returns the type of signal this class represents
    //====================================================================================================================
		unsigned short int getSignalType(void) {return this->signalType;}

    //====================================================================================================================
    //  If the old signal shouldn't be freed when calling setBuf_L or destructing the class (maybe because it's just 
		//  linked from elsewhere), this can be called first to remove the link
    //====================================================================================================================
    void forgetBuf_L(void) {this->Buf_L = NULL; this->Len = 0; return;};

    //====================================================================================================================
    //  True if everything is ok (can be called after the header-reading constructor call to see if it worked and handle 
		//  errors outside)
    //====================================================================================================================
    bool isOk(void) {return !(this->failure);};
};

#endif
