#ifndef __SC_FEATURE_SPECTRUMFLUX_H__
#define __SC_FEATURE_SPECTRUMFLUX_H__

#include <SV_Feature.h>

class SC_Feature_SpectrumFlux : public SV_Feature {

private:

protected:
	
public:

	SC_Feature_SpectrumFlux(int sampleRate, int frameSize, int frameStep, int fftLength, unsigned short window, double preemphasis);
	virtual ~SC_Feature_SpectrumFlux();

	//====================================================================================================================
  //	exptracts Spectrum Flux per frame: the avergage difference per frequency-bin in magnitude spectra of two succesive 
  //  frames
  //====================================================================================================================
  virtual SV_Data *ExtractFeature(void);
};

#endif
