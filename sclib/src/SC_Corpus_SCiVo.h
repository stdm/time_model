/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle video-files for SCiVo         */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.04.2006																								*/
/**************************************************************************/

#ifndef __SC_Corpus_SCiVo_H__
#define __SC_Corpus_SCiVo_H__

#include "SC_Corpus.h"

class SCLIB_API SC_Corpus_SCiVo : public SC_Corpus {
	
  private:

  protected:
	
  public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Corpus_SCiVo(SC_TweakableParameters* pTweak, double videoFrameRate, const char *audioFileName, const char *sceneFileName, const char *segmentFileName);
		virtual ~SC_Corpus_SCiVo();

    //====================================================================================================================
		// provides the possibility to load the desired samples for conformity reasons (it is a straight-forward call of
    // the corresponding method in SC_SegmentationHandler)
    // the given borders reamin unchanged in this method
		//====================================================================================================================
    SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries = false);
};

#endif

