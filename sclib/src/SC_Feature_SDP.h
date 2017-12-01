/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for extracting Symmetrized Dot Patterns per frame         */
/*																																				*/
/*    Author  : Bing Shi                           												*/
/*    Date    : 11.04.2006																								*/
/**************************************************************************/

#ifndef __SC_Feature_SDP_H__
#define __SC_Feature_SDP_H__

#include <SV_Feature.h>

class SC_Feature_SDP : public SV_Feature {

	private :

	protected :

    //====================================================================================================================
    // parameters for creating SDPs
    //====================================================================================================================
		int m;
		int lag;
		int color;
		int n;
		int pictureSize;
		int tau;

	public :
		
    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_SDP(int sampleRate, int frameLength, int frameStep, double preemphasize, int m, int lag, int color, int n, int pictureSize, int tau);
		virtual ~SC_Feature_SDP();

    //====================================================================================================================
		// override base class method, return cepstral peaks vector sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
};

#endif
