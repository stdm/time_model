#ifndef __SC_Feature_BANDPERIODICITY_H__
#define __SC_Feature_BANDPERIODICITY_H__

#include <SV_Feature.h>

class SC_Feature_BandPeriodicity : public SV_Feature {

private:

protected:

  float *filterSignal(float* signal, unsigned long int len, double lowCut, double highCut);

  /*
  float* filterSignal(double low, double high);
  float calculateBP(double low, double high);
  float* generateCorrFunc(float* previousFrame, float* currentFrame);
  float* generateCorrFunc(int currentFrameOff, float *sig, int len);
  float* extractLocalMaxima(float *sig, int len);
  */

public:

	SC_Feature_BandPeriodicity(int sampleRate, int frameSize, int frameStep, double preemphasis, unsigned short int FFTsize, unsigned short int window);
	virtual ~SC_Feature_BandPeriodicity();

  //====================================================================================================================
  //	extract Band Periodicity feature per frame: the maximum peak in the normalized correlation fuction between the 
  //  current and the previous frame per frequency sub-band (0-1/8, 1/8-1/, 1/4-1/2, 1/2-1 normalized frequency)
  //====================================================================================================================
	virtual SV_Data *ExtractFeature(void);
};

#endif
