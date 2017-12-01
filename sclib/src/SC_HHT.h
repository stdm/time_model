//simple method for Hilbert-Transform
//Author: Jun,Zhou
//Date:   March, 2007

#ifndef __SC_HHT_H__
#define __SC_HHT_H__

#include "../../svlib/src/SV_General.h"
#include "../../svlib/src/SV_Data.h"
#include "SC_TweakableParameters.h"

class SCLIB_API SC_HHT{

	private:
		SC_TweakableParameters *pTweak;

  protected:
	long int up_j;
	long int down_j;
	//double yp1;
	//double ypn;

	public :
		SC_HHT(SC_TweakableParameters *pTweak);
		virtual ~SC_HHT();

		//decomposition: Signal->IMFs
		SV_Data *emd(short *Signal, long int len);
		
		//calculate instantaneous frequency(based on Hz) pro sample
		SV_Data *insten_frequency(SV_Data *imf, long int row, long int col,double SampleRate,SV_Data *F, SV_Data *A);

		//another method to calculate instantaneous frequency(based on Hz) pro sample(only a similation for instantaneous frequency run quickly but not exactly und  )
		SV_Data *insten_frequency_op(SV_Data *imf, long int row, long int col,double SampleRate,SV_Data *F);
		
		//Hilbert-spectrum in matrix, value of element are the Amplitude in instantaneous frequency und time point
		SV_Data *hilbert_spectrum(SV_Data *frequency, SV_Data *H,SV_Data *imf, long int row, long int col);

		//get all the positions of maxima of signal, and store to the array, up_xs
		//get all the positions of minima of signal, and store to the array, down_xs
		//get all the values of maxima of signal, and store to the array, up_ys
		//get all the values of minima of signal, and store to the array, down_ys
		void max_min(double *h,long int *up_xs,long int *down_xs,double *up_ys,double *down_ys,long int Len);
		
		//interpolate the cubic-spline
		void interpolation(long int *xs,double *ys,long int length,long int Len,double *spline);

		//to test, whether the signal nearly momenton is
		bool tooSmall_simple(double *res,long int length,double alfa,double Thredhold);
};

#endif
