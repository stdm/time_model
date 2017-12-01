/**************************************************************************/
/*    Copied from:																												*/
/*      - This is a altered version of SC_Feature_FbE               			*/
/*																																				*/
/*    Responsibility:																											*/
/*      - Class for extracting Powerspectrum features as used in the  		*/
/*				SC_Enhancement class.                                           */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 30.03.2005																								*/
/**************************************************************************/

#ifndef __SC_Feature_Spectrum_H__
#define __SC_Feature_Spectrum_H__

#include <SV_Feature.h>

class SC_Feature_Spectrum : public SV_Feature {

	private :

	protected :

		bool createPhase;
		bool logarithmize;

	public :

    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_Spectrum(int sampleRate, int frameLength, int frameStep, double preemphasize, int FFTsize, unsigned short window, bool logarithmize, bool createPhase);
		virtual ~SC_Feature_Spectrum();

    //====================================================================================================================
		// override base class method, return log-spectral vector sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
    
    //====================================================================================================================
		// the signal is already framed (each row in pFrames = one frame)... compute a new dataset containing their spectra
    //====================================================================================================================
    SV_Data *ExtractFeature(SV_Data *pFrames);
};

#endif
