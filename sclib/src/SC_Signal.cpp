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

#include <stdio.h>
#include <string.h>
#include <fstream>
#include "SC_Signal.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
// default constructor
//====================================================================================================================
SC_Signal::SC_Signal() : SV_Signal() {
  this->fileName = NULL;
  this->sampleCount = 0;
	this->failure = false;
	this->signalType = sclib::stGuess;
	this->Next = NULL;
}

//====================================================================================================================
//  This is a copy-constructor, which copys the header-info of a given class
//====================================================================================================================
SC_Signal::SC_Signal(SC_Signal* oldSignal) {
	//SV_Signal parameters
	this->SigPar = oldSignal->SigPar;

	//SC_Signal parameters
	if (oldSignal->fileName != NULL) {
		MArray_1D(this->fileName, sclib::bufferSize, char, "SC_Signal.SC_Signal: this->fileName");
		sprintf(this->fileName, "%s", oldSignal->fileName);
	} else {
		this->fileName = NULL;
	}
  this->sampleCount = oldSignal->sampleCount;
	this->failure = oldSignal->failure;
	this->signalType = oldSignal->signalType;
	this->Next = oldSignal->Next;
}

//====================================================================================================================
// destructor 
//====================================================================================================================
SC_Signal::~SC_Signal() {
	MFree_1D(this->fileName);
	this->Next = NULL;
}

//====================================================================================================================
//  sets Buf_L to the given pointer
//====================================================================================================================
void SC_Signal::setBuf_L(short* pSamples, unsigned long int length) {
  MFree_1D(this->Buf_L);
  this->Buf_L = pSamples; 
  this->Len = length; 
  return;
}

//====================================================================================================================
//  converts Buf_L to float and returns a SV_Data object containing it
//====================================================================================================================
SV_Data* SC_Signal::convertBuf_L(void) {
  SV_Data* pBuf_L = new SV_Data(this->Len, 1);

  for (int i = 0; i < this->Len; i++) {
    pBuf_L->Mat[i][0] = (float) this->Buf_L[i];
  }

  return pBuf_L;
}

//====================================================================================================================
//  returns true, if the currently loaded signal can be replaced from start to start+sampleLength with the given 
//  samples
//====================================================================================================================
bool SC_Signal::implantSamples(unsigned long int start, short* pSamples, unsigned long int sampleLength) {
  if (start+sampleLength > (unsigned long)this->Len) {
    return false;
  }

  for (unsigned long int i = start; i < start+sampleLength; i++) {
    this->Buf_L[i] = pSamples[i-start];
  }
  
  return true;
}
