//########################################################################
//  
// A C++ class library for automatic speech recognition and 
// speaker recognition (identification and verification). 
//
// This class library is provided "as is" without any express 
// or implied warranty of any kind with respect to this software. 
// In particular the author shall not be liable for any direct, 
// indirect, special, incidental or consequential damages arising 
// in any way from use of the software.
//
//
// Author   : Jialong He,  Copyright (C), all rights reserved. 
// Date     : May, 1999. 
//
// Contact  : Jialong_He@bigfoot.com, Jialong_He@homemail.com
// Web Page : www.bigfoot.com/~Jialong_He
//########################################################################

//************************************************************************
//    This file implements FFT, IFFT, DFT, IDFT and 
//    FFT for real signal.
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SV_General.h"
#include "SV_Error.h"
#include "GN_FFT.h"

#define Pi 3.1415926535

static char SV_LibID[] = "Copyright (c) by Jialong He";

//===========================================================
//  Constructor
//===========================================================
GN_FFT::GN_FFT (){



}


//===========================================================
//  Destructor
//===========================================================
GN_FFT::~GN_FFT (){



}

/**************************************************************************

log2 - base 2 logarithm

Returns base 2 log such that i = 2**ans where ans = log2(i).
if log2(i) is between two values, the larger is returned.

int Log2(int x)

*************************************************************************/
int Log2(int x) {

    unsigned int mask,i;

    if(x == 0) return(-1);      /* zero is an error, return -1 */

    x--;                        /* get the max index, x-1 */

    for(mask = 1 , i = 0 ; ; mask *= 2 , i++) {
        if(x == 0) return(i);   /* return log2 if all zero */
        x = x & (~mask);        /* AND off a bit */
    }
}

//===========================================================
//  In-place radix 2 decimation in time FFT
//  The result will be placed in Src buffer. The Len must be
//  2**m, such as 1024, 4096.
//===========================================================
void GN_FFT::fft(COMPLEX *x, int Len) {
	

    static COMPLEX *w=NULL;      /* used to store the w complex array */
    static int mstore = 0;       /* stores m for future reference */
    static int n = 1;            /* length of fft stored for future */
    COMPLEX u,temp,tm;
    COMPLEX *xi,*xip,*xj,*wptr;

    int m, i,j,k,l,le,windex;

    double arg,w_real,w_imag,wrecur_real,wrecur_imag,wtemp_real;


	//-------------------------------------------
	// Check if signal length is radix-2
	//-------------------------------------------
	m = Log2(Len);
	int Radix2 = 1 << m;  // find nearest radix-2 length
	if (Radix2 != Len) {
		REPORT_ERROR(SVLIB_BadArg, "Length is not radix-2");
	}

    //-----------------------------------------------

    if(m != mstore) {

    /* free previously allocated storage and set new m */

	if(w != NULL) delete(w);
	mstore = m;
	if(m == 0) return;       /* if m=0 then done */

    /* n = 2**m = fft length */

        n = 1 << m;
        le = n/2;

    /* allocate the storage for w */

	w = new COMPLEX[le-1];
	if(w==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
    }

	/* calculate the w values recursively */

        arg = 4.0*atan(1.0)/le;         /* PI/le calculation */
        wrecur_real = w_real = cos(arg);
        wrecur_imag = w_imag = -sin(arg);
        xj = w;
        for (j = 1 ; j < le ; j++) {
            xj->real = wrecur_real;
            xj->imag = wrecur_imag;
            xj++;
            wtemp_real = wrecur_real*w_real - wrecur_imag*w_imag;
            wrecur_imag = wrecur_real*w_imag + wrecur_imag*w_real;
            wrecur_real = wtemp_real;
        }
    }

	/* start fft */

    le = n;
    windex = 1;
    for (l = 0 ; l < m ; l++) {
        le = le/2;

	/* first iteration with no multiplies */

        for(i = 0 ; i < n ; i = i + 2*le) {
            xi = x + i;
            xip = xi + le;
            temp.real = xi->real + xip->real;
            temp.imag = xi->imag + xip->imag;
            xip->real = xi->real - xip->real;
            xip->imag = xi->imag - xip->imag;
            *xi = temp;
        }

	/* remaining iterations use stored w */

        wptr = w + windex - 1;
        for (j = 1 ; j < le ; j++) {
            u = *wptr;
            for (i = j ; i < n ; i = i + 2*le) {
                xi = x + i;
                xip = xi + le;
                temp.real = xi->real + xip->real;
                temp.imag = xi->imag + xip->imag;
                tm.real = xi->real - xip->real;
                tm.imag = xi->imag - xip->imag;             
                xip->real = tm.real*u.real - tm.imag*u.imag;
                xip->imag = tm.real*u.imag + tm.imag*u.real;
                *xi = temp;
            }
            wptr = wptr + windex;
        }
        windex = 2*windex;
    }            

	/* rearrange data by bit reversing */

    j = 0;
    for (i = 1 ; i < (n-1) ; i++) {
        k = n/2;
        while(k <= j) {
            j = j - k;
            k = k/2;
        }
        j = j + k;
        if (i < j) {
            xi = x + i;
            xj = x + j;
            temp = *xj;
            *xj = *xi;
            *xi = temp;
        }
    }
}



//===========================================================
//  In-place radix 2 decimation in time IFFT
//  The result will be placed in Src buffer. The Len must be
//  2**m, such as 1024, 4096.
//===========================================================
void GN_FFT::ifft(COMPLEX *x, int Len) {


    static COMPLEX *w=NULL;      /* used to store the w complex array */
    static int mstore = 0;       /* stores m for future reference */
    static int n = 1;            /* length of ifft stored for future */

    COMPLEX u,temp,tm;
    COMPLEX *xi,*xip,*xj,*wptr;

    int m, i,j,k,l,le,windex;

    double arg,w_real,w_imag,wrecur_real,wrecur_imag,wtemp_real, scale;

	//-------------------------------------------
	// Check if signal length is radix-2
	//-------------------------------------------
	m = Log2(Len);
	int Radix2 = 1 << m;  // find nearest radix-2 length
	if (Radix2 != Len) {
		REPORT_ERROR(SVLIB_BadArg, "Length is not radix-2");
	}

    //-----------------------------------------------
    if(m != mstore) {

		/* free previously allocated storage and set new m */

	if(w != NULL) delete(w);
	mstore = m;
	if(m == 0) return;       /* if m=0 then done */

	/* n = 2**m = inverse fft length */

        n = 1 << m;
        le = n/2;

	/* allocate the storage for w */

	w = new COMPLEX[le-1];
	if(w==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
	}

	/* calculate the w values recursively */

        arg = 4.0*atan(1.0)/le;         /* PI/le calculation */
        wrecur_real = w_real = cos(arg);
        wrecur_imag = w_imag = sin(arg);  /* opposite sign from fft */
        xj = w;
        for (j = 1 ; j < le ; j++) {
            xj->real = wrecur_real;
            xj->imag = wrecur_imag;
            xj++;
            wtemp_real = wrecur_real*w_real - wrecur_imag*w_imag;
            wrecur_imag = wrecur_real*w_imag + wrecur_imag*w_real;
            wrecur_real = wtemp_real;
        }
    }

	/* start inverse fft */

    le = n;
    windex = 1;
    for (l = 0 ; l < m ; l++) {
        le = le/2;

	/* first iteration with no multiplies */

        for(i = 0 ; i < n ; i = i + 2*le) {
            xi = x + i;
            xip = xi + le;
            temp.real = xi->real + xip->real;
            temp.imag = xi->imag + xip->imag;
            xip->real = xi->real - xip->real;
            xip->imag = xi->imag - xip->imag;
            *xi = temp;
        }

	/* remaining iterations use stored w */

        wptr = w + windex - 1;
        for (j = 1 ; j < le ; j++) {
            u = *wptr;
            for (i = j ; i < n ; i = i + 2*le) {
                xi = x + i;
                xip = xi + le;
                temp.real = xi->real + xip->real;
                temp.imag = xi->imag + xip->imag;
                tm.real = xi->real - xip->real;
                tm.imag = xi->imag - xip->imag;             
                xip->real = tm.real*u.real - tm.imag*u.imag;
                xip->imag = tm.real*u.imag + tm.imag*u.real;
                *xi = temp;
            }
            wptr = wptr + windex;
        }
        windex = 2*windex;
    }            

	/* rearrange data by bit reversing */

    j = 0;
    for (i = 1 ; i < (n-1) ; i++) {
        k = n/2;
        while(k <= j) {
            j = j - k;
            k = k/2;
        }
        j = j + k;
        if (i < j) {
            xi = x + i;
            xj = x + j;
            temp = *xj;
            *xj = *xi;
            *xi = temp;
        }
    }

	/* scale all results by 1/n */
    scale = (1.0/n);
    for(i = 0 ; i < n ; i++) {
        x->real = scale*x->real;
        x->imag = scale*x->imag;
        x++;
    }
}


/************************************************************

rfft - trig recombination real input FFT

Requires real array pointed to by x, pointer to complex
output array, y and the size of real FFT in power of
2 notation, m (size of input array and FFT, N = 2**m).
On completion, the COMPLEX array pointed to by y 
contains the lower N/2 + 1 elements of the spectrum.

void rfft(double *x, COMPLEX *y, int m)

***************************************************************/
void GN_FFT::rfft(double *x, COMPLEX *y, int Len) { 

    static    COMPLEX  *cf=NULL;
    static    int      mstore = 0;
    int       m,p,num,k;
    double     Realsum, Realdif, Imagsum, Imagdif;
    double    factor, arg;
    COMPLEX   *ck, *xk, *xnk, *cx;


	//-------------------------------------------
	// Check if signal length is radix-2
	//-------------------------------------------
	m = Log2(Len);
	int Radix2 = 1 << m;  // find nearest radix-2 length
	if (Radix2 != Len) {
		REPORT_ERROR(SVLIB_BadArg, "Length is not radix-2");
	}

    //-----------------------------------------------

    /* First call the fft routine using the x array but with
       half the size of the real fft */

    p = m - 1;
    cx = (COMPLEX *) x;
    fft(cx,	Len/2);

	/* Next create the coefficients for recombination, if required */

    num = 1 << p;    /* num is half the real sequence length.  */

    if (m!=mstore){
      if (cf != NULL) delete(cf);
      cf = new COMPLEX[num-1];

      if(cf==NULL){
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
      }

      factor = 4.0*atan(1.0)/num;
      for (k = 1; k < num; k++){
        arg = factor*k;
        cf[k-1].real = cos(arg);
        cf[k-1].imag = sin(arg);
      }
    }  

	/* DC component, no multiplies */
    y[0].real = cx[0].real + cx[0].imag;
    y[0].imag = 0.0;

	/* other frequencies by trig recombination */
    ck = cf;
    xk = cx + 1;
    xnk = cx + num - 1;
    for (k = 1; k < num; k++){
      Realsum = ( xk->real + xnk->real ) / 2;
      Imagsum = ( xk->imag + xnk->imag ) / 2;
      Realdif = ( xk->real - xnk->real ) / 2;
      Imagdif = ( xk->imag - xnk->imag ) / 2;

      y[k].real = Realsum + ck->real * Imagsum
                          - ck->imag * Realdif ;

      y[k].imag = Imagdif - ck->imag * Imagsum
                          - ck->real * Realdif ;
      ck++;
      xk++;
      xnk--;
    }
}

/***********************************************************************

dft - Discrete Fourier Transform
  
This function performs a straight DFT of N points on an array of
complex numbers whose first member is pointed to by Datain.  The
output is placed in an array pointed to by Dataout.

void dft(COMPLEX *Datain, COMPLEX *Dataout, int N)

*************************************************************************/
void GN_FFT::dft(COMPLEX *Datain, COMPLEX *Dataout, int N) {  // Len can be any

    int i,k,n,p;
    static int nstore = 0;      /* store N for future use */
    static COMPLEX *cf=NULL;    /* coefficient storage */
    COMPLEX *cfptr,*Dinptr;
    double arg;

	/* Create the coefficients if N has changed */

    if(N != nstore) {
	if(cf != NULL) delete(cf);    /* free previous */

	cf    = new COMPLEX[N];
	if (cf==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
	}	

        arg = 8.0*atan(1.0)/N;
        for (i=0 ; i<N ; i++) {
            cf[i].real = cos(arg*i);
            cf[i].imag = -sin(arg*i);
        }
    }

	/* Perform the DFT calculation */

    for (k=0 ; k<N ; k++) {

        Dinptr = Datain;
        Dataout->real = Dinptr->real;
        Dataout->imag = Dinptr->imag;
        Dinptr++;
        for (n=1; n<N; n++) {

        p = (int)((long)n*k % N);
            cfptr = cf + p;         /* pointer to cf modulo N */

            Dataout->real += Dinptr->real * cfptr->real
                             - Dinptr->imag * cfptr->imag;

            Dataout->imag += Dinptr->real * cfptr->imag
                             + Dinptr->imag * cfptr->real;
            Dinptr++;
        }
	Dataout++;          /* next output */
    }
}

/***********************************************************************

idft - Inverse Discrete Fourier Transform
  
This function performs an inverse DFT of N points on an array of
complex numbers whose first member is pointed to by Datain.  The
output is placed in an array pointed to by Dataout.
It returns nothing.

void idft(COMPLEX *Datain, COMPLEX *Dataout, int N)

*************************************************************************/
void GN_FFT::idft(COMPLEX *Datain, COMPLEX *Dataout, int N) {  // Len can be any


    int i,k,n,p;
    static int nstore = 0;      /* store N for future use */
    static COMPLEX *cf = NULL;         /* coefficient storage */
    COMPLEX *cfptr,*Dinptr;
    double arg;

	/* Create the coefficients if N has changed */

    if(N != nstore) {
	if(cf != NULL) delete (cf);    /* free previous */

	cf = new COMPLEX[N];
	if (cf == NULL) {
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
	}


	/* scale stored values by 1/N */
        arg = 8.0*atan(1.0)/N;
        for (i=0 ; i<N ; i++) {
            cf[i].real = (cos(arg*i)/(double)N);
            cf[i].imag = (sin(arg*i)/(double)N);
        }
    }

	/* Perform the DFT calculation */

    for (k=0 ; k<N ; k++) {

        Dinptr = Datain;
        Dataout->real = Dinptr->real * cf[0].real;
        Dataout->imag = Dinptr->imag * cf[0].real;
        Dinptr++;
        for (n=1; n<N; n++) {

        p = (int)((long)n*k % N);
            cfptr = cf + p;         /* pointer to cf modulo N */

            Dataout->real += Dinptr->real * cfptr->real
                             - Dinptr->imag * cfptr->imag;

            Dataout->imag += Dinptr->real * cfptr->imag
                             + Dinptr->imag * cfptr->real;
            Dinptr++;
        }
	Dataout++;          /* next output */
    }
}

//===========================================================
//  Convert Real/Imag representation to Mag/Phase. Put Mag
//  in In_Out.real and Phase in In_Out.imag.
//===========================================================
void GN_FFT::Mag_Angle(COMPLEX *In_Out, int Len){

	double Mag, Angle;

	for (int Cnt=0; Cnt<Len; Cnt++) {
	
		Mag   = In_Out[Cnt].real * In_Out[Cnt].real + In_Out[Cnt].imag * In_Out[Cnt].imag;
		Mag   = sqrt(Mag);
		Angle = 0;
		if (Mag != 0) {
  		   Angle = acos(In_Out[Cnt].real / Mag);
        } 

		//----------------------------------------
		// acos gives 0-pi, turn it to [-pi, pi]
		//----------------------------------------
		if (In_Out[Cnt].imag < 0) {
			Angle = -Angle;
		}

		//----------------------------------------
		// put Mag in real, and Angle in imag
		//----------------------------------------
		In_Out[Cnt].real = Mag;
		In_Out[Cnt].imag = Angle;

	}  // for (Cnt)
}

    
//===========================================================
//  Convert Mag/Angle representation to Real/Imag.
//  Mag is in In_Out[].real, Angle in .imag
//===========================================================
void GN_FFT::Real_Imag(COMPLEX *In_Out, int Len){

	double Real, Imag;

	for (int Cnt=0; Cnt<Len; Cnt++) {

		Real = In_Out[Cnt].real * cos( In_Out[Cnt].imag );
		Imag = In_Out[Cnt].real * sin( In_Out[Cnt].imag );
        In_Out[Cnt].real =  Real; 
        In_Out[Cnt].imag =  Imag; 
    }
}


//===========================================================
//  In-place COSINE transform for real sequence. 
//  The sequence length must be 2**m (because of using FFT).
//  (code ported from Matlab's dct.m) 
//===========================================================
void GN_FFT::dct(double *In_Out, int Len) {
	
	COMPLEX *Sig;
	int Cnt;

	//--------------------------------
	// Allocate Temp buffer
	//--------------------------------
	Sig = new COMPLEX[Len];
	if(Sig==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
    }


	//------------------------------------------------------
	// re-order sequence 
	// (0, 1, 2, 3, 4, 5, 6, 7) => (0, 2, 4, 6, 7, 5, 3, 1)
	//------------------------------------------------------
	for (Cnt=0; Cnt<Len/2; Cnt++) {
		Sig[Cnt].real = In_Out[Cnt*2];	
		Sig[Cnt].imag = 0.0;	
	}

	for (Cnt=Len/2; Cnt<Len; Cnt++) {
		Sig[Cnt].real = In_Out[Len*2 - Cnt*2 - 1];	
		Sig[Cnt].imag = 0.0;	
	}


    fft(Sig, Len);   // doing FFT

	//------------------------------------------------------
	// Multiplying FFT results with complex weights 
	//------------------------------------------------------
	for (Cnt=0; Cnt<Len; Cnt++) {

		In_Out[Cnt] = (cos(Cnt*Pi/(2*Len)) * Sig[Cnt].real + 
			                   sin(Cnt*Pi/(2*Len)) * Sig[Cnt].imag) 
					            / sqrt(2.0*Len) * 2;
	}

    In_Out[0] /= sqrt(2.0);  // First component scale by 

	delete(Sig);
};	


//===========================================================
//  In-place inverse COSINE transform for real sequence. 
//  The sequence length must be 2**m (because of using IFFT).
//  (code ported from Matlab's idct.m) 
//===========================================================
void GN_FFT::idct(double *In_Out, int Len) {

	COMPLEX *Sig;
	int Cnt;

	//--------------------------------
	// Allocate Temp buffer
	//--------------------------------
	Sig = new COMPLEX[Len];
	if(Sig==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "Not Enought Memory");
    }

	//------------------------------------------------------
	// Multiplying input by weight 
	//------------------------------------------------------
	for (Cnt=0; Cnt<Len; Cnt++) {
		Sig[Cnt].real =  cos(Cnt*Pi/(2*Len)) * In_Out[Cnt] * sqrt(2.0*Len);
		Sig[Cnt].imag =  sin(Cnt*Pi/(2*Len)) * In_Out[Cnt] * sqrt(2.0*Len);
	}

	Sig[0].real /= sqrt(2.0);
	Sig[0].imag /= sqrt(2.0);

	ifft(Sig, Len);

	//------------------------------------------------------
	// re-order sequence 
	// (0, 1, 2, 3, 4, 5, 6, 7) <= (0, 2, 4, 6, 7, 5, 3, 1)
	//------------------------------------------------------
	for (Cnt=0; Cnt<Len/2; Cnt++) {
		In_Out[Cnt*2] = Sig[Cnt].real;	
	}

	for (Cnt=Len/2; Cnt<Len; Cnt++) {
		In_Out[Len*2 - Cnt*2 - 1] = Sig[Cnt].real;	
	}

	delete(Sig);


};	

