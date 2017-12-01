/**************************************************************************/
/*    Responsibility:																											*/
/*		  - derived from SC_GroundTruth to provide a groundtruth-object for */
/*        MIX2MAX's singer-recognition data                                */
/*      - because here is no real groundtruth, this class is lacking a    */
/*        framelist and only provides the time-conversion functions       */
/*      - the only knowledge of this class is the original sample-rate    */
/*        of the unseen signal on which the loaded feature are based      */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.04.2006																								*/
/**************************************************************************/

#ifndef __SC_GroundTruth_Wesley_H__
#define __SC_GroundTruth_Wesley_H__

#include "SC_Api.h"
#include "SC_GroundTruth.h"

class SCLIB_API SC_GroundTruth_Wesley : public SC_GroundTruth {

  private:

	protected :

		//====================================================================================================================
		// Initialize the frameList with it's new 4. column
		//====================================================================================================================
		virtual void initFrameList(void);

    //====================================================================================================================
		// Read Groundtruth from file(s) and initialize internal data structures (frameList etc.)
		//====================================================================================================================
    virtual bool readGroundTruth(void);

	public :

    //====================================================================================================================
		// Constructor, Destructor
		//====================================================================================================================
    SC_GroundTruth_Wesley(SC_TweakableParameters *pTweak, unsigned long int sampleRate);
		virtual ~SC_GroundTruth_Wesley();

    //====================================================================================================================
		// change the sampleRate
		//====================================================================================================================
    void setSampleRate(unsigned long int sampleRate) {this->audioSampleRate = sampleRate; return;}
};

#endif
