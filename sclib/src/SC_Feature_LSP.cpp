/**************************************************************************/
/*    This class implements the Line Spectra Pairs feature, that is       */
/*    derived from LPCs (or the raw signal in this case, too). It         */
/*    includes the same information as LPCs but inherits more robustness. */
/*    The implementation is based on the corresponding speex code         */
/*    (http://www.speex.org/) or ephone MELP Proposed Federal Standard    */
/*    speech coder code and not much more than a wrapper around both.     */
/*    (see copyright etc. in the corresponding .cpp file)                 */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 12.02.2006																								*/
/**************************************************************************/

#include <math.h>
#include "SC_Aux.h"
#include "SC_Feature_LSP.h"
#include "SC_Feature_LPC.h"
#include <SV_Error.h>

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_LSP::SC_Feature_LSP(int sampleRate, int frameLength, int frameStep, unsigned int window, double preemphasize, int LPCorder, int method, double delta, int bisections, double minSeparation, int maxLoops) : SV_Feature(){
  this->Para.WinSz = frameLength;
	this->Para.StpSz = frameStep;
  this->Para.Alpha = preemphasize;
	this->Para.LPC_Order = LPCorder;
	this->Para.SRate = sampleRate;
	this->Para.HammingWin = window;

	this->method = method;
	this->delta = delta;
	this->bisections = bisections;
	this->minSeparation = minSeparation;
	this->maxLoops = maxLoops;
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_LSP::~SC_Feature_LSP() {
	
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data* SC_Feature_LSP::ExtractFeature(void) {
	SV_Data *pLPCs = NULL, *pLSPs = NULL;
	SC_Feature_LPC *pLPCmaker = new SC_Feature_LPC(this->Para.SRate, this->Para.WinSz, this->Para.StpSz, this->Para.HammingWin, this->Para.Alpha, this->Para.LPC_Order, false);

	//get LPCs
	pLPCmaker->setSignal(this->Sig, this->Len, false);
	pLPCs = pLPCmaker->ExtractFeature();
	pLPCmaker->setSignal(NULL, 0, false);
	MFree_0D(pLPCmaker);

	//convert them to LSPs
	if (pLPCs != NULL) {
		pLSPs = ExtractFeature(pLPCs);
		MFree_0D(pLPCs);
	}

  return pLSPs;
}

//====================================================================================================================
// the real wrapper for the speex-functions; turn LPCs into LSPs; return NULL on error
// a good explanation on what is done and how to set inout and parameters can be found here:
// http://www.mathworks.com/access/helpdesk/help/toolbox/dspblks/index.html?/access/helpdesk/help/toolbox/dspblks/lpctolsflspconversion.html
//====================================================================================================================
SV_Data* SC_Feature_LSP::ExtractFeature(SV_Data *pLPCs) {
	SV_Data *pData = NULL;
	int roots = 0, x, y;
	float *lpcs = NULL, *lsps = NULL;
	char *tmp = new char[2*pLPCs->Col*sizeof(float)]; //this is surely big enough
	bool needsOne = (pLPCs->Mat[0][0] != 1.0) ? true : false; 
  
	if (pLPCs != NULL) {
		pData = new SV_Data;
		if (pData==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		}

		pData->Row = pLPCs->Row;
		pData->Col = pLPCs->Col;
		pData->Alloc();
		pData->Hdr.frameSize = this->Para.WinSz;
		pData->Hdr.frameStep = this->Para.StpSz;
		pData->Hdr.sampleRate = this->Para.SRate;
		pData->Hdr.ID = sclib::featureLSP;
		pData->Hdr.Signature[1] = (char)(this->method); //encode which extractor was used

		if (needsOne == true) {
			MArray_1D(lpcs, pLPCs->Col+1, float, "SC_Feature_LSP.ExtractFeature: lpcs");
		}
		if (this->method == sclib::modeMELP) {
			MArray_1D(lsps, pLPCs->Col+1, float, "SC_Feature_LSP.ExtractFeature: lsps");
			for (x = 0; x < pLPCs->Col+1; x++) {
				lsps[x] = (float)(0.0); //initialize values to zero 'cause maybe not all roots will be found and we have no way here to detect which ones where missed... after the first run, not the zeros but the previous values of the missing coefficients will be sued which is an even better approximate value
			}
		}

		for (y = 0; y < pLPCs->Row; y++) {
			if (needsOne == true) { //the analysis algorithms below need the full LPC set including the first coeeficient, that always equals 1.0 and therefore is ofton ommitted in LPC parameter sets
				lpcs[0] = 1.0;
				for (x = 0; x < pLPCs->Col; x++) {
					lpcs[x+1] = pLPCs->Mat[y][x];
				}
			} else {
				lpcs = pLPCs->Mat[y];
			}
			if (this->method == sclib::modeSpeex) {
				//LPC to LSPs (x-domain) transform
#ifdef SC_USE_SPEEXLSP
				roots = lpc_to_lsp(lpcs, pLPCs->Col, pData->Mat[y], this->bisections, (float)(this->delta), tmp);
				for (x = roots; x < pLPCs->Col; x++) { //there are maybe fewer LSPs than LPCs, so set those uninitialized coefficients =-1 (minimal value of the roots, as it seems... but i'm unsure (TODO))
					pData->Mat[y][x] = -1.0;
				}
#else
				REPORT_ERROR(SVLIB_Fail, "Speex code not available");
				MFree_0D(pData);
#endif
			} else {
#ifdef SC_USE_MELPLSP
				roots = lpc_pred2lsp(lpcs, lsps, pLPCs->Col, (float)(this->minSeparation), (float)(this->delta), this->bisections, this->maxLoops); //roots is 1 on error (missed root) here, 0 on success
				for (x = 0; x < pLPCs->Col; x++) {
					pData->Mat[y][x] = lsps[x+1]; //the MELP code doesn't fill the first LSP in and begins with array-index [1]
				}
#else
				REPORT_ERROR(SVLIB_Fail, "MELP code not available");
				MFree_0D(pData);
#endif
			}
		}

		MFree_1D(lsps);
		if (needsOne == true) {
			MFree_1D(lpcs);
		}
	}
	MFree_1D(tmp);

  return pData;
}

//====================================================================================================================
// Below is the code (4 methods) directly (without modifications) taken from the speex project
//====================================================================================================================

#ifdef SC_USE_SPEEXLSP
/*---------------------------------------------------------------------------*\
Original copyright
	FILE........: AKSLSPD.C
	TYPE........: Turbo C
	COMPANY.....: Voicetronix
	AUTHOR......: David Rowe
	DATE CREATED: 24/2/93

Modified by Jean-Marc Valin

   This file contains functions for converting Linear Prediction
   Coefficients (LPC) to Line Spectral Pair (LSP) and back. Note that the
   LSP coefficients are not in radians format but in the x domain of the
   unit circle.

   Speex License:

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   
   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//Aligns the stack to a 'size' boundary; by thilo: changed typecast of stack from int to intptr_t for 64bit compliance
#define ALIGN(stack, size) ((stack) += ((size) - (intptr_t)(stack)) & ((size) - 1))

//Allocates 'size' elements of type 'type' on the stack
#define PUSH(stack, size, type) (ALIGN((stack),sizeof(type)),(stack)+=((size)*sizeof(type)),(type*)((stack)-((size)*sizeof(type))))

//Allocates a struct stack
#define PUSHS(stack, type) (ALIGN((stack),sizeof(long)),(stack)+=(sizeof(type)),(type*)((stack)-(sizeof(type))))

/*---------------------------------------------------------------------------*\

	FUNCTION....: lpc_to_lsp()

	AUTHOR......: David Rowe
	DATE CREATED: 24/2/93

    This function converts LPC coefficients to LSP
    coefficients.

\*---------------------------------------------------------------------------*/
int SC_Feature_LSP::lpc_to_lsp (float *a,int lpcrdr,float *freq,int nb,float delta, char *stack)
/*  float *a                    lpc coefficients                        */
/*  int lpcrdr                  order of LPC coefficients (10)          */
/*  float *freq                 LSP frequencies in the x domain         */
/*  int nb                      number of sub-intervals (4)             */
/*  float delta                 grid spacing interval (0.02)            */


{

    float psuml,psumr,psumm,temp_xr,xl,xr,xm=0;
    float temp_psumr/*,temp_qsumr*/;
    int i,j,m,flag,k;
    float *Q;                   /* ptrs for memory allocation           */
    float *P;
    float *px;                  /* ptrs of respective P'(z) & Q'(z)     */
    float *qx;
    float *p;
    float *q;
    float *pt;                  /* ptr used for cheb_poly_eval()
                                whether P' or Q'                        */
    int roots=0;                /* DR 8/2/94: number of roots found     */
    flag = 1;                   /*  program is searching for a root when,
                                1 else has found one                    */
    m = lpcrdr/2;               /* order of P'(z) & Q'(z) polynomials   */


    /* Allocate memory space for polynomials */
    Q = PUSH(stack, (m+1), float);
    P = PUSH(stack, (m+1), float);

    /* determine P'(z)'s and Q'(z)'s coefficients where
      P'(z) = P(z)/(1 + z^(-1)) and Q'(z) = Q(z)/(1-z^(-1)) */

    px = P;                      /* initialise ptrs                     */
    qx = Q;
    p = px;
    q = qx;
    *px++ = 1.0;
    *qx++ = 1.0;
    for(i=1;i<=m;i++){
        *px++ = a[i]+a[lpcrdr+1-i]-*p++;
        *qx++ = a[i]-a[lpcrdr+1-i]+*q++;
    }
    px = P;
    qx = Q;
    for(i=0;i<m;i++){
        *px = 2**px;
        *qx = 2**qx;
         px++;
         qx++;
    }
    px = P;                     /* re-initialise ptrs                   */
    qx = Q;

    /* Search for a zero in P'(z) polynomial first and then alternate to Q'(z).
    Keep alternating between the two polynomials as each zero is found  */

    xr = 0;                     /* initialise xr to zero                */
    xl = 1.0;                   /* start at point xl = 1                */


    for(j=0;j<lpcrdr;j++){
        if(j%2)                 /* determines whether P' or Q' is eval. */
            pt = qx;
        else
            pt = px;

        psuml = cheb_poly_eva(pt,xl,lpcrdr,stack);      /* evals poly. at xl    */
        flag = 1;
        while(flag && (xr >= -1.0)){
				//while(flag){
           float dd;
           /* Modified by JMV to provide smaller steps around x=+-1 */
           dd=(float)(delta*(1-.9*xl*xl)); //by thilo: added explicit cast to avoid warning
           if (fabs(psuml)<.2)
              dd *= .5;

           xr = xl - dd;                                /* interval spacing     */
            psumr = cheb_poly_eva(pt,xr,lpcrdr,stack);/* poly(xl-delta_x)       */
            temp_psumr = psumr;
            temp_xr = xr;

    /* if no sign change increment xr and re-evaluate poly(xr). Repeat til
    sign change.
    if a sign change has occurred the interval is bisected and then
    checked again for a sign change which determines in which
    interval the zero lies in.
    If there is no sign change between poly(xm) and poly(xl) set interval
    between xm and xr else set interval between xl and xr and repeat till
    root is located within the specified limits                         */

            if((psumr*psuml)<0.0){
                roots++;

                psumm=psuml;
                for(k=0;k<=nb;k++){
                    xm = (xl+xr)/2;             /* bisect the interval  */
                    psumm=cheb_poly_eva(pt,xm,lpcrdr,stack);
                    if(psumm*psuml>0.){
                        psuml=psumm;
                        xl=xm;
                    }
                    else{
                        psumr=psumm;
                        xr=xm;
                    }
                }

               /* once zero is found, reset initial interval to xr      */
               freq[j] = (xm);
               xl = xm;
               flag = 0;                /* reset flag for next search   */
            }
            else{
                psuml=temp_psumr;
                xl=temp_xr;
            }
        }
    }
    return(roots);
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: cheb_poly_eva()

	AUTHOR......: David Rowe
	DATE CREATED: 24/2/93

    This function evaluates a series of Chebyshev polynomials

\*---------------------------------------------------------------------------*/
float SC_Feature_LSP::cheb_poly_eva(float *coef,float x,int m,char *stack) //by thilo: removed 'static' specifier
/*  float coef[]  	coefficients of the polynomial to be evaluated 	*/
/*  float x   		the point where polynomial is to be evaluated 	*/
/*  int m 		order of the polynomial 			*/
{
    int i;
    float *T,sum;
    int m2=m>>1;

    /* Allocate memory for Chebyshev series formulation */
    T=PUSH(stack, m2+1, float);

    /* Initialise values */
    T[0]=1;
    T[1]=x;

    /* Evaluate Chebyshev series formulation using iterative approach  */
    /* Evaluate polynomial and return value also free memory space */
    sum = coef[m2] + coef[m2-1]*x;
    x *= 2;
    for(i=2;i<=m2;i++)
    {
       T[i] = x*T[i-1] - T[i-2];
       sum += coef[m2-i] * T[i];
    }
    
    return sum;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: lsp_to_lpc()

	AUTHOR......: David Rowe
	DATE CREATED: 24/2/93

    lsp_to_lpc: This function converts LSP coefficients to LPC
    coefficients.

\*---------------------------------------------------------------------------*/
void SC_Feature_LSP::lsp_to_lpc(float *freq,float *ak,int lpcrdr, char *stack)
/*  float *freq 	array of LSP frequencies in the x domain	*/
/*  float *ak 		array of LPC coefficients 			*/
/*  int lpcrdr  	order of LPC coefficients 			*/


{
    int i,j;
    float xout1,xout2,xin1,xin2;
    float *Wp;
    float *pw,*n1,*n2,*n3,*n4=NULL;
    int m = lpcrdr/2;

    Wp = PUSH(stack, 4*m+2, float);
    pw = Wp;

    /* initialise contents of array */

    for(i=0;i<=4*m+1;i++){       	/* set contents of buffer to 0 */
	*pw++ = 0.0;
    }

    /* Set pointers up */

    pw = Wp;
    xin1 = 1.0;
    xin2 = 1.0;

    /* reconstruct P(z) and Q(z) by  cascading second order
      polynomials in form 1 - 2xz(-1) +z(-2), where x is the
      LSP coefficient */

    for(j=0;j<=lpcrdr;j++){
       int i2=0;
	for(i=0;i<m;i++,i2+=2){
	    n1 = pw+(i*4);
	    n2 = n1 + 1;
	    n3 = n2 + 1;
	    n4 = n3 + 1;
	    xout1 = xin1 - 2*(freq[i2]) * *n1 + *n2;
	    xout2 = xin2 - 2*(freq[i2+1]) * *n3 + *n4;
	    *n2 = *n1;
	    *n4 = *n3;
	    *n1 = xin1;
	    *n3 = xin2;
	    xin1 = xout1;
	    xin2 = xout2;
	}
	xout1 = xin1 + *(n4+1);
	xout2 = xin2 - *(n4+2);
	ak[j] = (float)((xout1 + xout2)*0.5); //by thilo: added explicit cast to avoid warning
	*(n4+1) = xin1;
	*(n4+2) = xin2;

	xin1 = 0.0;
	xin2 = 0.0;
    }

}

/*Added by JMV
  Makes sure the LSPs are stable*/
void SC_Feature_LSP::lsp_enforce_margin(float *lsp, int len, float margin)
{
   int i;
   if (lsp[0]<margin)
      lsp[0]=margin;
   if (lsp[len-1]>sclib::pi-margin)
      lsp[len-1]=(float)(sclib::pi-margin); //by thilo: added explicit cast to avoid warning
   for (i=1;i<len-1;i++)
   {
      if (lsp[i]<lsp[i-1]+margin)
         lsp[i]=lsp[i-1]+margin;

      if (lsp[i]>lsp[i+1]-margin)
         lsp[i]= (float)(.5* (lsp[i] + lsp[i+1]-margin)); //by thilo: added explicit cast to avoid warning
   }
}
#endif 

//====================================================================================================================
// Below is the code (4 methods) taken from the TI's ephone / "MELP Proposed Federal Standard speech coder" project
//====================================================================================================================

#ifdef SC_USE_MELPLSP
/*
2.4 kbps MELP Proposed Federal Standard speech coder

version 1.2

Copyright (c) 1996, Texas Instruments, Inc.  

Texas Instruments has intellectual property rights on the MELP
algorithm.  The Texas Instruments contact for licensing issues for
commercial and non-government use is William Gordon, Director,
Government Contracts, Texas Instruments Incorporated, Semiconductor
Group (phone 972 480 7442).
*/

/* LPC_PRED2LSP
      get LSP coeffs from the predictor coeffs
      Input:
         a- the predictor coefficients
         p- the predictor order
      Output:
         w- the lsp coefficients
   Reference:  Kabal and Ramachandran

*/
int SC_Feature_LSP::lpc_pred2lsp(float *a,float *w,int p,float lsp_delta,float root_delta,int root_bisections,int clmp_max_loops) //by thilo: added lsp_delta (0.0), root_delta (0.00781250), root_bisections (4) and clmp_max_loops (10) parameter
{
    int i,p2;
    float **c;

    p2 = p/2;

    //MEM_2ALLOC(MALLOC,c,2,p2+1,float);
		MArray_2D(c, 2, p2+1, float, "SC_Feature_LSP.lpc_pred2lsp: c"); //by thilo
    c[0][p2] = c[1][p2] = 1.0;

    for(i=1; i <= p2; i++)
    {
        c[0][p2-i] = (a[i] + a[p+1-i] - c[0][p2+1-i]);
        c[1][p2-i] = c[1][p2+1-i] + a[i] - a[p+1-i];
    }
    c[0][0] /= 2.0;
    c[1][0] /= 2.0;

    i = lsp_roots(w,c,p2,root_delta,root_bisections); //by thilo: added root_* parameters

    if (i)
    {
        for(i=1; i <= p; i++)
            (void)fprintf(stderr,"%11.7f ",a[i]);
        (void)fprintf(stderr,"\n");
    }

    /* ensure minimum separation and sort */
    (void)lpc_clmp(w,lsp_delta,p,clmp_max_loops); //by thilo: used un-aliased function name; added max_loops parameter

    //MEM_2FREE(FREE,c);
		MFree_2D(c); //by thilo
    return(i);
} /* LPC_PRED2LSP */

/* LSP_ROOTS
        - find the roots of the two polynomials G_1(x) and G_2(x)
          the first root corresponds to G_1(x)
          compute the inverse cos (and these are the LSFs) */
int SC_Feature_LSP::lsp_roots(float *w,float **c,int p2,float delta,int bisections) //by thilo: added delta (0.00781250) and bisections parameter (4); removed static specifier before the function
{
    int i,k;
    float x,x0,x1,y,*ptr,g0,g1;

		//by thilo: replaced usage of macros BISECTIONS and DELTA with corresponding parameters; added explicit typecasts to float to avoid warnings

    w[0] = 0.0;

    ptr = c[0];
    x = 1.0;
    g0 = lsp_g(x,ptr,p2);

    for(k=1,x = (float)(1.0)-delta; x > -delta-(float)(1.0); x -= delta)
    {
        /* Search for a zero crossing */
        if (g0*(g1 = lsp_g(x,ptr,p2)) <= 0.0)
        {
            /* Search Incrementally using bisection */
            x0 = x+delta;
            x1 = x;

            for(i=0; i < bisections; i++)
            {
                x = (x0+x1)/(float)(2.0);
                y = lsp_g(x,ptr,p2);

                if(y*g0 < 0.0)
                {
                    x1 = x;
                    g1 = y;
                }
                else
                {
                    x0 = x;
                    g0 = y;
                }
            }
            /* Linear interpolate */
            x = (g1*x0-g0*x1)/(g1-g0);

            /* Evaluate the LSF */
            w[k] = (float)acos((double)x)/(float)(sclib::pi);

            ptr = c[k % 2];
            k++;
            if (k > 2*p2)
                return(0);
            g1 = lsp_g(x,ptr,p2);
        }
        g0 = g1;
    }
    (void)fprintf(stderr,"\n Error(lsp_roots): LSPs Not All Found\n");
    return(1);
} /* LSP_ROOTS */

/* G - compute the value of the Chebychev series
                sum c_k T_k(x) = x b_1(x) - b_2(x) + c_0
                b_k(x) = 2x b_{k+1}(x) - b_{k+2}(x) + c_k */
float SC_Feature_LSP::lsp_g(float x,float *c,int p2) //by thilo: removed static specifier before the function
{
    int i;
    float b[3];

    b[1] = b[2] = 0.0;

    for(i=p2; i > 0; i--)
    {
        b[0] = (float)(2.0)*x*b[1] - b[2] + c[i];
        b[2] = b[1];
        b[1] = b[0];
    }
    b[0] = x*b[1]-b[2]+c[0];
    return(b[0]);
} /* G */

/*
    Name: lpc_clmp- Sort and ensure minimum separation in LSPs.
    Aliases: lpc_clamp
    Description:
        Ensure that all LSPs are ordered and separated
        by at least delta.  The algorithm isn't guarenteed
        to work, so it prints an error message when it fails
        to sort the LSPs properly.
    Inputs:
        w- lsp vector (order p, w[1..p])
        delta- the clamping factor
        p- order of lpc filter
    Outputs:
        w- the sorted and clamped lsps
    Returns: NULL
    See_Also:
    Includes:
        spbstd.h
        lpc.h
    Bugs: 
        Currently only supports 10 loops, which is too
        complex and perhaps unneccesary.

    Systems and Info. Science Lab
    Copyright (c) 1995 by Texas Instruments, Inc.  All rights reserved.
*
*/
int SC_Feature_LSP::lpc_clmp(float *w, float delta, int p, int max_loops) //by thilo: added max_loops (10) parameter
{
    int i,j,unsorted;
    float tmp,d,step1,step2;

		//by thilo: replaced macro MAX_LOOPS with correspondig parameter; added explicit typecasts to float to avoid warnings

    /* sort the LSPs- for 10 loops, complexity is approximately 150 p */  
    for (j=0,unsorted=TRUE; unsorted && (j < max_loops); j++)
    {
        for(i=1,unsorted=FALSE; i < p; i++)
            if (w[i] > w[i+1])
            {
                tmp = w[i+1];
		w[i+1] = w[i];
		w[i] = tmp;
                unsorted = TRUE;
	    }
    }

    /* ensure minimum separation */
    if (!unsorted) 
    {

        for(j=0; j < max_loops; j++)
        {
            for(i=1; i < p; i++)
            {
                if ((d = w[i+1]-w[i]) < delta)
                {
                    step1 = step2 = (delta-d)/(float)(2.0);
                    if (i==1 && (w[i] < delta))
                    {
                        step1 = w[i]/(float)(2.0);
	            }
                    else if (i > 1)
                    {
                        if ((tmp = w[i] - w[i-1]) < delta)
                            step1 = 0;
                        else if (tmp < 2*delta)
                            step1 = (tmp-delta)/(float)(2.0);
		    }
                    if (i==(p-1) && (w[i+1] > ((float)(1.0)-delta)))
                    {
                        step2 = (1-w[i+1])/(float)(2.0);
		    }
                    else if (i < (p-1))
                    {
                        if ((tmp = w[i+2] - w[i+1]) < delta)
                            step2 = 0;
                        else if (tmp < 2*delta)
                            step2 = (tmp-delta)/(float)(2.0);
		    }
                    w[i] -= step1;
		    w[i+1] += step2;
	        }
	    }
        }
    }

    /* Debug: check if the minimum separation rule was met */
    for(j=1; j < p; j++)
      if ((w[j+1]-w[j]) < 0.99*delta)
          (void)fprintf(stderr,"%s: LSPs not separated by enough (line %d)\n",
              __FILE__,__LINE__);

    if (unsorted)
        (void)fprintf(stderr,"%s: LSPs still unsorted (line %d)\n",
		      __FILE__,__LINE__);

    return(0);
}
#endif
