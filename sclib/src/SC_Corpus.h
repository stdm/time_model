/**************************************************************************/
/*    Responsibility:																											*/
/*      - base class for objects to encapsulates all algorithms and       */
/*        methods to handle a specific corpus of data and it's specific   */
/*        tasks. because the corpora are highly different and common      */
/*        taks and methods are seldom, there are only few methods in the  */
/*        base class                                                      */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.04.2006																								*/
/**************************************************************************/

#ifndef __SC_Corpus_H__
#define __SC_Corpus_H__

#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Signal.h"
#include "SC_SignalHandler.h"

class SCLIB_API SC_Corpus {
	
  private:

  protected:
	
    SC_TweakableParameters *pTweak;
    SC_GroundTruth *pGT;
		SC_SignalHandler *pSignalHandler; //we need it ias a member so it doesn't get destructed during consecutive load-operations, so that the MPEG-cache works...
 
  public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Corpus(SC_TweakableParameters* pTweak);
		virtual ~SC_Corpus();

    //====================================================================================================================
		// Provide Access to the GroundTruth
		//====================================================================================================================
    SC_GroundTruth* getGT(void) {return this->pGT;}

    //====================================================================================================================
		// Maybe the signals for a corpus can't just be loaded via the methods provided the by SC_SignalHandler class; so, 
    // for conformity reasons, this method is placed here and can be overwritten by the children, even to just call
    // the signal-handler
    // the given borders are references because they might get cahnged a little sclib::bit by this method (see SC_Corpus_TIMIT
    // for an example)
		//====================================================================================================================
    virtual SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries = false) = 0;

		//====================================================================================================================
    // write the given feature-set to a file
		//====================================================================================================================
    virtual void featureOut(char *fileName, char *featureName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data *pData, bool printHeader = false);

		//====================================================================================================================
		// outputs the wanted low-level-features of all voiced speech frames in the given segment, preceded by each frame's 
		// gt-speaker-id
		//====================================================================================================================
		virtual void speakerFeaturesOut(char *fileName, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pData, unsigned long int selectionMask = 0xFFFFFFFF);

		//====================================================================================================================
    // Return the fileName of the audio file containg the given startSample
		//====================================================================================================================
		virtual char* getCurrentAudioFileName(unsigned long int startSample);
};

#endif
