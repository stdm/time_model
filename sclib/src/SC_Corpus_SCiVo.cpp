/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle video-files for SCiVo         */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.04.2006																								*/
/**************************************************************************/

#include "SC_Corpus_SCiVo.h"
#include "SC_SignalHandler.h"
#include "SC_GroundTruth_SCiVo.h"
#include "SC_Aux.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Corpus_SCiVo::SC_Corpus_SCiVo(SC_TweakableParameters* pTweak, double videoFrameRate, const char *audioFileName, const char *sceneFileName, const char *segmentFileName) : SC_Corpus(pTweak) {
  this->pGT = new SC_GroundTruth_SCiVo(this->pTweak, videoFrameRate, audioFileName, sceneFileName, segmentFileName);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Corpus_SCiVo::~SC_Corpus_SCiVo() {

}

//====================================================================================================================
// provides the possibility to load the desired samples for conformity reasons (it is a straight-forward call of
// the corresponding method in SC_SegmentationHandler)
// the given borders reamin unchanged in this method
//====================================================================================================================
SC_Signal* SC_Corpus_SCiVo::loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries) {
  SC_Signal *pSignal = this->pSignalHandler->loadSignal(this->pGT->getAudioFileName(), segmentStart, sclib::min(segmentEnd, this->pGT->getAudioSampleCount()-1), this->pTweak->signalHandler.forceSampleRate);

  return pSignal;
}

