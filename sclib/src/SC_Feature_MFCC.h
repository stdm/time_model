/**************************************************************************/
/*    Responsibility:																											*/
/*      - Class for extracting MFCC's via different algorithms: svlib's		*/
/*        original implementation or sclib's new one (via SC_Feature_FbE) */
/*																																				*/
/*    Author  : Thilo Stadelmann															            */
/*    Date    : 27.02.2004																								*/
/**************************************************************************/

#ifndef __SC_Feature_MFCC_H__
#define __SC_Feature_MFCC_H__

#include <SV_Feature.h>

class SC_Feature_MFCC : public SV_Feature {

	private :

	protected :

		int method;
		unsigned char coeffSelection[8];
		bool CMN;
		bool addDeltas; 
		bool addDeltaDeltas;
		int svlib_removeSilence; 
		unsigned char sclib_frequencyScale; 
		unsigned char sclib_smoothing;
		double sclib_minFilterBankFrequency;
		double sclib_maxFilterBankFrequency;

	public :

    //====================================================================================================================
		// constructor/destructor
    //====================================================================================================================
		SC_Feature_MFCC(int sampleRate, int frameSize, int frameStep, int MFCCorder, int filterBankSize, int fftSize, unsigned short int window, double preemphasize, unsigned char coeffSelection[8], int dEnergy, bool CMN, bool addDeltas, bool addDeltaDeltas, int method, unsigned char sclib_frequencyScale, unsigned char sclib_smoothing, double sclib_minFilterBankFrequency, double sclib_maxFilterBankFrequency);
		virtual ~SC_Feature_MFCC();

    //====================================================================================================================
		// override base class method, return MFCC vector sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);

    //====================================================================================================================
		// method to select some specific coefficients out of the feature-vector
		// selectionString is a bit-flag-map, where each setted bit corresponds to a col of pData which shouldn't be discarded
		// the first bit in selectionString[0] corresponds with the first column, the last bit in selectionstring[7] with the
		// 64th column.
		//====================================================================================================================
		void selectCoefficients(SV_Data* &pData, unsigned char selectionString[8]);
};

#endif
