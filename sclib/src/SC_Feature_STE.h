/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for extracting short time Energy								          */
/*			- Copy from SC_Feature_EnergyZCR					                        */
/*                                                                        */
/*    Author  : Bing Shi																					        */
/*    Date    : 11.04.2006																								*/
/**************************************************************************/

#ifndef __SC_Feature_STE_H__
#define __SC_Feature_STE_H__

#include <SV_Feature.h>

class SC_Feature_STE: public SV_Feature {

	private :

	protected :

    //====================================================================================================================
    //	A hardcoded implementation of the Butterworth bandpass filter as specified in "Cepstrum-Based Pitch Detection 
		//  Using a New Statistical V/UV Classification Algorithm" (Ahmadi, Spanias 1999): 9th order, 8kHz sample-rate, pass-
		//  band 200-3400Hz; the filter was generated using http://www-users.cs.york.ac.uk/~fisher/mkfilter 
    //====================================================================================================================
		void butterworthBandPass(float* signal, long int length);

		bool useButterworth;
		bool scaleResult; //per-frame => per-ms

	public :

		// constructor/destructor
		SC_Feature_STE(int sampleRate, int frameSize, int frameStep, double preEmphasize, bool useButterworth, bool scaleResult);
		virtual ~SC_Feature_STE();

		// override base class method, return ZCR vector sequence
		virtual SV_Data *ExtractFeature(void);
};

#endif
