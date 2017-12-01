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

#include "SC_GroundTruth_Wesley.h"
#include <SV_DataIO.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_GroundTruth_Wesley::SC_GroundTruth_Wesley(SC_TweakableParameters *pTweak, unsigned long int sampleRate) : SC_GroundTruth(pTweak, NULL) {
  //some init-stuff
  this->gtType = sclib::gtWesley;
  this->audioFileName = NULL;
  this->pSignalPrototype = NULL;
  this->audioSampleCount = 0;
  this->internalFrameCount = 0;
  this->internalFrameSize = 0;
  this->pSpeakerMapping = NULL;
  for (int i = 0; i < sclib::maxSpeakers; i++) {
    this->speakerNames[i] = NULL;
  }
	this->uncertaintyRegion = 0;

  this->audioSampleRate = sampleRate;
	this->pConverter->setAudioSampleRate(this->audioSampleRate);
  
  initFrameList();
  //readGroundTruth();
}

//====================================================================================================================
//	destructor 
//====================================================================================================================
SC_GroundTruth_Wesley::~SC_GroundTruth_Wesley() {

}

//====================================================================================================================
//	initialize the framelist
//====================================================================================================================
void SC_GroundTruth_Wesley::initFrameList(void) {
  this->frameList = NULL;
	return;
}

//====================================================================================================================
// Read Groundtruth
//====================================================================================================================
bool SC_GroundTruth_Wesley::readGroundTruth(void) {
  return true;
}

