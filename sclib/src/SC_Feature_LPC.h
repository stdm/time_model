//************************************************************************
//    Calculate LPC coefficients.
//
//    Author  : Jun Zhou
//    Date    : March 22, 2006
//************************************************************************

#ifndef __SC_Feature_LPC_H__
#define __SC_Feature_LPC_H__

#include "SC_Aux.h"
#include <SV_Feature.h>

class SC_Feature_LPC : public SV_Feature {

private:

protected:

	bool computeGain; //by thilo: if true, the gain term is appended to each feature vector (in a last, additional column)

public:

	SC_Feature_LPC(int sampleRate, int frameSize, int frameStep, unsigned int window, double preemphasis, unsigned short LPCorder, bool computeGain);
	virtual ~SC_Feature_LPC();

  //====================================================================================================================
  //	extract LPC per frame
  //====================================================================================================================
	virtual SV_Data *ExtractFeature(void);
};

#endif
