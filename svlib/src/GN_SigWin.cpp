//************************************************************************
//    Implemention for signal window class. 
//
//
//    Author  : Jialong HE
//    Date    : March 24, 1999
//************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SV_Error.h"
#include "GN_SigWin.h"

#define Pi 3.1415926535
static char SV_LibID[] = "Copyright (c) by Jialong He";

//===========================================================
//  Constructor
//===========================================================
GN_SigWin::GN_SigWin (){

	SigWin	= NULL;
	Length	= 0;
	WinType	= HAMMING; 
}

//===========================================================
//  Destructor
//===========================================================
GN_SigWin::~GN_SigWin (){
	MFree_1D(SigWin);
}

//=================================================================
//  Allocate memory for Signal Window Buffer 
//=================================================================
void GN_SigWin::Alloc(int Len) {

	if (Length != Len) {

		if (SigWin != NULL) {
			MFree_1D(SigWin);
		}

		MArray_1D(SigWin, Len, double, "AlloSigWinBuf");
		Length = Len;
	}

}



//===========================================================
//  Hamming Window
//===========================================================
double *GN_SigWin::Hamming(int Len){

	double Factor;

	Alloc(Len);
	//------------------------------
	// Fill window buffer
	//------------------------------
    Factor = 2.0 * Pi / (Length-1.0);
	for (int Cnt = 0; Cnt<Length; Cnt++){
		SigWin[Cnt] = 0.54 - 0.46*cos(Factor * Cnt);
	}

	WinType	= HAMMING; 
	return (SigWin);
}


//===========================================================
//  Hanning Window
//===========================================================
double *GN_SigWin::Hanning(int Len){

	double Factor;

	Alloc(Len);
	//------------------------------
	// Fill window buffer
	//------------------------------
    Factor = 2.0 * Pi / (Length-1.0);
	for (int Cnt = 0; Cnt<Length; Cnt++){
		SigWin[Cnt] = 0.5 - 0.5*cos(Factor * Cnt);
	}

	WinType	= HANNING; 
	return (SigWin);
}

//===========================================================
//  Harris Window
//===========================================================
double *GN_SigWin::Harris(int Len){

	double Factor;

	Alloc(Len);
	//------------------------------
	// Fill window buffer
	//------------------------------
    Factor = 2.0 * Pi / Length;
	for (int Cnt = 0; Cnt<Length; Cnt++){
		SigWin[Cnt] = 0.35875 - 0.48829*cos(Factor*Cnt) + 
					0.14128*cos(2.0*Factor*Cnt) - 0.01168*cos(3.0*Factor*Cnt);
	}

	WinType	= HARRIS; 
	return (SigWin);
}

//===========================================================
//  Blackman Window
//===========================================================
double *GN_SigWin::Blackman(int Len){

	double Factor;

	Alloc(Len);
	//------------------------------
	// Fill window buffer
	//------------------------------
    Factor = 2.0 * Pi / (Length-1.0);
	for (int Cnt = 0; Cnt<Length; Cnt++){
		SigWin[Cnt] = 0.42 - 0.5*cos(Factor*Cnt) + 0.08*cos(2*Factor*Cnt);
	}

	WinType	= BLACKMAN; 
	return (SigWin);
}


//===========================================================
//  Triangle Window
//===========================================================
double *GN_SigWin::Triangle(int Len){

	int Cnt;	
	Alloc(Len);
	//------------------------------
	// Fill window buffer
	//------------------------------
	for (Cnt = 0; Cnt<Length/2; Cnt++){
		SigWin[Cnt] = double(Cnt) / (Length-1) * 2.0;
	}

	for (Cnt = Length/2; Cnt<Length; Cnt++){
		SigWin[Cnt] = 2.0 - double(Cnt) / (Length-1) * 2.0;
	}

	WinType	= TRIANGLE; 
	return (SigWin);
}


//===========================================================
//  Multiplying Signal by spcified window function
//  return Window function
//===========================================================
double *GN_SigWin::ApplyWindow(double *Sig, int Len, int NewType){

	if (Sig==NULL || Len<=0) {
		REPORT_ERROR(SVLIB_BadArg, "ApplyWindow");
	}

	//---------------------------------------
	// If current window function is different
	//---------------------------------------
	if (Length != Len || WinType != NewType) {
		switch (NewType) {
			case HAMMING: 
				Hamming(Len);
				break;

			case HANNING:
				Hanning(Len);
				break;

			case TRIANGLE:
				Triangle(Len);
				break;

			case BLACKMAN:
				Blackman(Len);
				break;

			case HARRIS:
				Harris(Len);
				break;
		}
	}

	//---------------------------------------
	// Multiplying signal by window
	//---------------------------------------
	for (int Cnt=0; Cnt<Len; Cnt++) {
		Sig[Cnt] *= SigWin[Cnt];
	}

	return (SigWin);
}


