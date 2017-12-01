#ifndef __SC_Feature_SubBandPower_H__
#define __SC_Feature_SubBandPower_H__

#include <SV_Feature.h>

class SC_Feature_SubBandPower : public SV_Feature {

  private:

  protected:

    //====================================================================================================================
    //  sum up the entries in the array from "from" to "to"
    //====================================================================================================================
    double sumArray(double *values, int from, int to);

  public:

	  SC_Feature_SubBandPower(int sampleRate, int frameSize, int frameStep, int fftLength, unsigned short window, double preemphasis);
	  virtual ~SC_Feature_SubBandPower();

    //====================================================================================================================
    //	extracts per frame: 1. Short Time Energy (derived from log of integral of powerspectrum)
    //                      2. Ratio of sub-band 1 (integral of 0   - 1/8 powerspectrum) to exp(STE)
    //                      3. Ratio of sub-band 2 (integral of 1/8 - 1/4 powerspectrum) to exp(STE)
    //                      4. Ratio of sub-band 3 (integral of 1/4 - 1/2 powerspectrum) to exp(STE)
    //                      5. Ratio of sub-band 4 (integral of 1/2 - 1   powerspectrum) to exp(STE)
    //====================================================================================================================
	  virtual SV_Data *ExtractFeature(void);
};

#endif
