/**************************************************************************/
/*    Responsibility:																											*/
/*      - Provides the signal samples in a SV_Data container and chopped  */
/*        into frames                                                     */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 29.01.2009																								*/
/**************************************************************************/

#ifndef __SC_Feature_Samples_H__
#define __SC_Feature_Samples_H__

#include <SV_Feature.h>

class SC_Feature_Samples : public SV_Feature {

	private :

	protected :

		bool createPhase;
		bool logarithmize;

	public :

    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_Samples(int sampleRate, int frameLength, int frameStep);
		virtual ~SC_Feature_Samples();

    //====================================================================================================================
		// overrides base class method, returns framed samples
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
};

#endif
