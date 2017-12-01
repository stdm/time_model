#ifndef __SC_Feature_BrightnessBandwidth_H__
#define __SC_Feature_BrightnessBandwidth_H__

#include <SV_Feature.h>

class SC_Feature_BrightnessBandwidth : public SV_Feature {

private:

protected:

public:

	SC_Feature_BrightnessBandwidth(int sampleRate, int frameSize, int frameStep, int fftLength, unsigned short window, double preemphasis);
	virtual ~SC_Feature_BrightnessBandwidth();

  //====================================================================================================================
  //	extract brightness & bandwidth per frame
  //====================================================================================================================
	virtual SV_Data *ExtractFeature(void);
};

#endif
