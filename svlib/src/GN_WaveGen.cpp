//************************************************************************
//    Implementation for WaveForm Generator 
//
//
//    Author  : Jialong HE
//    Date    : March 24, 1999
//************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SV_Error.h"
#include "GN_WaveGen.h"
#include "GN_Rand.h"

#define Pi 3.1415926535
static char SV_LibID[] = "Copyright (c) by Jialong He";
//=================================================================
//  Constructor 
//=================================================================
GN_WaveGen::GN_WaveGen() {
	
	Wave      = NULL;
	Length    = 0;

}

//=================================================================
//  Destructor 
//=================================================================
GN_WaveGen::~GN_WaveGen() {
	
	if (Wave != NULL) {
		MFree_1D(Wave);
		Length = 0;
	}

}

//=================================================================
//  Allocate memory for Wave buffer 
//=================================================================
void GN_WaveGen::AlloWaveBuf(int Len) {

	if (Length != Len) {

		if (Wave != NULL) {
			MFree_1D(Wave);
		}

		MArray_1D(Wave, Len, double, "WaveGen")
		Length = Len;
	}

}


//=================================================================
//  Sine Wave 
//=================================================================
double *GN_WaveGen::Tone_Wave(double Cycle, double Phase, int Len) {

	AlloWaveBuf(Len);
	
	for (int Cnt=0; Cnt<Length; Cnt++) {
		Wave[Cnt] = sin(Cycle*2.0*Pi*Cnt/Length + Phase);
	}


	return(Wave);
}

//=================================================================
//  Square Wave 
//=================================================================
double *GN_WaveGen::Square_Wave(int HighDur, int LowDur, int Len) {

	int Ind, Cnt;

	AlloWaveBuf(Len);
	//---------------------------------------
	// Fill Wave buffer with square wave
	//---------------------------------------
	Ind = 0;
	while (Ind < Length) {

		//--------------------------------
		// fill high part
		//--------------------------------
		Cnt = 0;
		while (Ind < Length && Cnt < HighDur) {
			Wave[Ind] = 1.0;
			Ind++;
			Cnt++;
		}

		//--------------------------------
		// fill low part
		//--------------------------------
		Cnt = 0;
		while (Ind < Length && Cnt < LowDur) {
			Wave[Ind] = 0.0;
			Ind++;
			Cnt++;
		}

	}   // while (Ind<Length)

	return(Wave);
}


//=================================================================
//  Sawtooth Wave 
//=================================================================
double *GN_WaveGen::Sawtooth_Wave(int RampDur, int LowDur, int Len) {

	int Ind, Cnt;

	AlloWaveBuf(Len);

	//---------------------------------------
	// Fill Wave buffer with square wave
	//---------------------------------------
	Ind = 0;
	while (Ind < Length) {

		//--------------------------------
		// fill ramp part
		//--------------------------------
		Cnt = 0;
		while (Ind < Length && Cnt < RampDur) {
			Wave[Ind] = double(Cnt) / RampDur;
			Ind++;
			Cnt++;
		}

		//--------------------------------
		// fill low part
		//--------------------------------
		Cnt = 0;
		while (Ind < Length && Cnt < LowDur) {
			Wave[Ind] = 0.0;
			Ind++;
			Cnt++;
		}

	}   // while (Ind<Length)

	return(Wave);
}


//=================================================================
//  Tone modulation Wave 
//=================================================================
double *GN_WaveGen::ModTone_Wave(double ModCyl, double CarCyl, int Len) {


	AlloWaveBuf(Len);

	//-----------------------------------
	// Generate Carrier wave
	//-----------------------------------
	for (int Cnt=0; Cnt<Length; Cnt++) {
		Wave[Cnt] =  (0.5 + 0.5*sin(ModCyl*2.0*Pi*Cnt/Length - Pi/2)) * // Modulator 
					 sin(CarCyl*2.0*Pi*Cnt/Length);             // Carrier
	}

	return(Wave);

}

//=================================================================
//  Uniform distribution noise (white noise)
//  0-1
//=================================================================
double *GN_WaveGen::NoiseUni_Wave(int Len) {

	GN_Rand RandEng;

	AlloWaveBuf(Len);

	for(int Cnt=0; Cnt<Length; Cnt++) {

		Wave[Cnt] = (double)RandEng.random() / RandEng.getmax();

	}

	return(Wave);
}


//=================================================================
//  Gaussian distribution noise 
//=================================================================
double *GN_WaveGen::NoiseGaus_Wave(double Mean, double Std, int Len) {

	GN_Rand RandEng;

	AlloWaveBuf(Len);
	for(int Cnt=0; Cnt<Length; Cnt++) {

		Wave[Cnt] = RandEng.rand_gaus(Mean, Std); // gaus distribution
	}

	return(Wave);
}



