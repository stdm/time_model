#ifndef __SC_Feature_NFR_H__
#define __SC_Feature_NFR_H__

#include <SV_Feature.h>

class SC_Feature_NFR : public SV_Feature {

private:

protected:

  double nfrThreshold;

public:

	SC_Feature_NFR(int sampleRate, int frameSize, double preemphasis, double nfrThreshold, int FFTsize, unsigned short window);
	virtual ~SC_Feature_NFR();

  //====================================================================================================================
  //	extract the Noise Frame Ratio feature per frame; this means: it is not really a ratio, just a measure if a frame 
  //  is considered as noise (1) or not (0); the ratio has to computed later for a given sub-clip
  //====================================================================================================================
	virtual SV_Data *ExtractFeature(void);
};

#endif
