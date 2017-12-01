/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Signal to load SPHERE headered(NIST) files with some comfort */
/*                                                                        */
/*    Also based heavily on Jialong He's SV_Signal_NIST class             */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 07.04.2006																								*/
/**************************************************************************/

#ifndef __SC_Signal_NIST_H__
#define __SC_Signal_NIST_H__

#include "SC_Signal.h"

class SC_Signal_NIST : public SC_Signal {

  private :

  protected :

    //====================================================================================================================
    // Important Header Info
    //====================================================================================================================
    typedef struct {
      short channel_count;
      short sample_n_bytes;
      long  sample_count;
      long  sample_rate;
      short SwapByte; //1: machine's byte order is different from data.
      short sample_coding;
    } SC_HeaderTypeSPHERE;

    SC_Signal_NIST::SC_HeaderTypeSPHERE header;

    //====================================================================================================================
    //  This procedure read a few inportant fields in SPHERE header. if success, return (1), otherwise return(0) 
    //====================================================================================================================
    int loadHeader(const char *FName, SC_Signal_NIST::SC_HeaderTypeSPHERE &HInfo);

    //====================================================================================================================
    //  This procedure can be called to load speech samples from NIST's SPHERE format data file.
    //  if sucess, return the number of samples in SigBuf, otherwise, return(-1).
    //  Speech buf allocated from outside
    //====================================================================================================================
    int LoadNIST(const char *FName, unsigned long int start, unsigned long int end, SC_Signal_NIST::SC_HeaderTypeSPHERE &HInfo, short* &buffer, int &length);

  public :

    //====================================================================================================================
	  //  Constructor
    //====================================================================================================================
	  SC_Signal_NIST();

    //====================================================================================================================
    //  This is a copy-constructor, which copys the header-info of a given class
    //====================================================================================================================
    SC_Signal_NIST(SC_Signal_NIST* oldSignal);

    //====================================================================================================================
    //  This constructor analyses the header without loading samples into memory
    //====================================================================================================================
    SC_Signal_NIST(const char* fileName);
    
    //====================================================================================================================
	  //  Destructor
    //====================================================================================================================
    virtual ~SC_Signal_NIST();

    //====================================================================================================================
    //  Load comlete signal into memory              
    //====================================================================================================================
	  virtual long LoadSignal(const char *FName);

    //====================================================================================================================
    //  Load signal from start-sample to end-sample into memory 
    //  this method expects the header already to be filled by the constructor!              
    //====================================================================================================================
    virtual long LoadSignal(const char *fileName, unsigned long int start, unsigned long int end);

    //====================================================================================================================
    //  Generate a new audio-file with the header of this (adapted to the new data) and the samples stored in Buf_L
    //====================================================================================================================
    virtual long SaveSignal(const char *FName);
};


#endif
