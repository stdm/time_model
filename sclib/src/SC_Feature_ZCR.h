/**************************************************************************/
/*    Class to extract Zero Crossing Rate Feature                         */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 12.02.2006																								*/
/**************************************************************************/

#ifndef __SC_Feature_ZCR_H__
#define __SC_Feature_ZCR_H__

#include <SV_Feature.h>

class SC_Feature_ZCR : public SV_Feature {

	private:

	protected:

    //====================================================================================================================
    //	A hardcoded implementation of the Chebyshev highpass filter as specified in "Cepstrum-Based Pitch Detection 
		//  Using a New Statistical V/UV Classification Algorithm" (Ahmadi, Spanias 1999): 9th order, 8kHz sample-rate, pass-
		//  band 100-4000Hz; the filter was generated using http://www-users.cs.york.ac.uk/~fisher/mkfilter 
    //====================================================================================================================
		void chebyshevHighPass(float* signal, long int length);

		bool useChebyshev;
		bool scaleResults;

	public:

		SC_Feature_ZCR(int sampleRate, int frameSize, int frameStep, double preemphasis, bool useChebyshev, bool scaleResults);
		virtual ~SC_Feature_ZCR();

		//====================================================================================================================
		//	extract ZCR per frame, normalized by 2*(frameSize-1)
		//====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
};

#endif
