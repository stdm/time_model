//************************************************************************
//    Calculate LPC residual.
//
//    Author  : Jun Zhou
//    Date    : April 10, 2006
//************************************************************************

#ifndef __SC_Feature_LPCresidual_H__
#define __SC_Feature_LPCresidual_H__

#include <SV_Feature.h>
#include "SC_Signal.h"

class SC_Feature_LPCresidual : public SV_Feature {

	private:

	protected:

	public:

		SC_Feature_LPCresidual(int sampleRate = 16000, int frameSize = 512, int frameStep = 256, double preemphasis = 0.97, int order = 4);
		virtual ~SC_Feature_LPCresidual();

		//====================================================================================================================
		//	extract LPC residual from Signal per frame, 
		//====================================================================================================================
		virtual SV_Data *ExtractFeature(void);

		//====================================================================================================================
		//	constructs a residual signal (that e.g. can be saved and listened to) from the extracted framed (overlapping) 
		//  residual samples
		//  by thilo
		//====================================================================================================================
	  SC_Signal* getResidualSignal(SV_Data *pResidual);

		//====================================================================================================================
		//	construct a residual signal (see above; stored in a one column SV_Data container with 1 sample per row) from the 
		//  given residual frames of voiced speech that conatins only those samples within a short window around the glottal
		//  closure instants; see Dhananjaya, Yegnanarayana, "Speaker Change Detection in Casual Conversations Using 
		//  Excitation Source Features", 2008; Naylor, Kounoudes, Gudnason, Brookes, "Estimation of Glottal Closure Instants 
		//  in Voiced Speech Using the DYPSA Algorithm", 2007
		//  phaseSlopeWindowSize is the size of the window used to compute the average energy weighted phase slope function in 
		//  [ms].
    //  windowSizeAroundGCI is the size of the window around the found GCIs from which the samples are taken for the final 
		//  result (in [ms]).
		//  by thilo
		//====================================================================================================================
		SV_Data* getResidualSignalAroundGCIs(SV_Data *pVoicedSpeechResidual, float phaseSlopeWindowSize = 3.0, float windowSizeAroundGCI = 5.0);

		//====================================================================================================================
		//	converts the result from getResidualSignalAroundGCIs() back to a frame-based representation with 1 frame per row
		//  in the result-set; 5ms frameSize should be used
		//  by thilo
		//====================================================================================================================
		SV_Data* residualSignal2frames(SV_Data *pResidualSignal, int frameSize, int frameStep = 1, bool doNormalization = true);

		//====================================================================================================================
		//	converts the result from getResidualSignal() back to a frame-based representation with 1 frame per row in the 
		//  result-set; 5ms frameSize should be used
		//  by thilo
		//====================================================================================================================
		SV_Data* residualSignal2frames(SC_Signal *pResidualSignal, int frameSize, int frameStep = 1, bool doNormalization = true);

		//====================================================================================================================
		//	(re)-set the clas members
		//  by thilo
		//====================================================================================================================
		void setParameters(int sampleRate, int frameSize, int frameStep, double preemphasis, int order);
};

#endif
