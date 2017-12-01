//************************************************************************
//    Calculte LPC coefficients and LPC based cepstral coefficients.
//
//
//    Author  : Jialong HE
//    Date    : March 15, 1999
//************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //by thilo, for sqrt() in CalcLPC()
#include "GN_LPC.h"
#include "SV_Error.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//=================================================================
//  Default Constructor
//=================================================================
GN_LPC::GN_LPC() {
	LPC_Order = 0;
}

//=================================================================
//  Default Destructor 
//=================================================================
GN_LPC::~GN_LPC() {


}

//=================================================================
//  Calculate Auto-Correlation coefficients 
//=================================================================
double GN_LPC::AutoCorrelation(float *Signal, int Samples, double *CorCoeff, int Order) { //by thilo: changed return type from void to double to return original CorCoeff[0] so to be able to reverse the normalization applied here
	int  Cnt, LoopCnt;
	double originalZerothCorrCoeff; //by thilo

	for (Cnt = 0; Cnt < Order; Cnt++) {
		CorCoeff[Cnt] = 0.0;
		for (LoopCnt=0; LoopCnt<(Samples-Cnt); LoopCnt++) {
			CorCoeff[Cnt] += Signal[LoopCnt] * Signal[LoopCnt + Cnt];
    }
  };

	originalZerothCorrCoeff = CorCoeff[0]; //by thilo

	/*---------------------------------------*/
	/* Now normalizing the autocorrelation.  */
	/* all divided by first element          */
	/*---------------------------------------*/
  if (CorCoeff[0] != 0) {
     for (Cnt = Order-1; Cnt >= 0; Cnt--) CorCoeff[Cnt] /= CorCoeff[0];
  }
    
  return originalZerothCorrCoeff; //by thilo
}

//=================================================================
//  Using  Durbin's recursive method calculating LPC coefficients
//  from AutoCorrelation coefficients.
//=================================================================
void GN_LPC::DUrbin(double *CorrCoeff, double *LPC_Coef, double *Reflect) {

   int       Cnt, j;
   double    Error, alpha;
   double    *E, *LPCArray;

   //-----------------------------------
   // Tempory buffer for prediction error
   //-----------------------------------
   MArray_1D(E, (LPC_Order+2), double, "Durbin_E");
   MArray_1D(LPCArray, (LPC_Order+2), double, "Durbin_LPC");
   
   //---------------------------------------
   // Recursive calculated LPC coefficients
   //---------------------------------------
   LPCArray[1] = 1.0;
   E[0]   = CorrCoeff[0];
   if ( E[0] != 0)  Reflect[1] = -CorrCoeff[1] / E[0];
   alpha = E[0] * (1.0 - Reflect[1]*Reflect[1]);
   LPCArray[2] = Reflect[1];

   for (Cnt=2; Cnt<=LPC_Order; Cnt++) {
      Error = 0.0;
      for (j = 1; j <= Cnt; j++) {
				Error += LPCArray[j] * CorrCoeff[Cnt + 1 - j];
			}
	  
			if (alpha == 0) {
				//REPORT_ERROR(SVLIB_DivBy0, "Divided by 0!");
				Reflect[Cnt] = -Error/0.00001; //by thilo to avoid error
			}
			else {Reflect[Cnt] = -Error/alpha;}

			alpha *= (1.0 - Reflect[Cnt]*Reflect[Cnt]);
			LPCArray[Cnt + 1] = Reflect[Cnt];

			for (j = 2; j <= Cnt; j++)
				E[j] = LPCArray[j] + Reflect[Cnt]*LPCArray[Cnt + 2 - j];

			for (j = 2; j <= Cnt; j++)
				LPCArray[j] = E[j];
	 };

   //---------------------------------------------
   // Re-organize LPC coefficients, index from 0
   //---------------------------------------------
   for (Cnt = 0; Cnt < LPC_Order; Cnt++)
      LPC_Coef[Cnt] = LPCArray[Cnt + 2];

   MFree_1D(E);
   MFree_1D(LPCArray);
}

//=================================================================
//  Derive LPC coefficients for a segment of signal
//=================================================================
double GN_LPC::CalcLPC (float* Sig, int Num, double* LPC_Coef, int Order) { //by thilo: changed return value from void to double to return the LPC filter's gain
	double *AutoCorr, *Reflect;
	double originalZerothCorrCoeff, gain = 0.0; //by thilo

  //-----------------------------------------
	// Allocate Temp memory
	//-----------------------------------------
  MArray_1D(Reflect, (Order+1), double, "Reflect"); 
  MArray_1D(AutoCorr, (Order+1), double, "AutoCorr"); 

	LPC_Order = Order;
	//by thilo: catch new return value of unnormalized zeroth autocorr coefficient that is needed later on for gain computation
  originalZerothCorrCoeff = AutoCorrelation(Sig, Num, AutoCorr, Order+1);  // p+1 AutoCoeff -> p LPC

	//----------------------------------------------------
	// ***NOTE*** LPC coefficients should be (1, a1, a2, a3, ...)
	// but the first "1" is not returned in LPC_Coef.
	//----------------------------------------------------
	DUrbin(AutoCorr, LPC_Coef, Reflect); 
	
	//block by thilo: compute gain; heres a description how to do that, submitted by cranen on Mon, 2006-12-11 18:04 to http://www.isca-students.org/lpc_gain_term:
	// "An elaborate description of how the Gain (G) for the LPC model can be computed can be found in 
	//  Rabiner, L.R. and Schafer,R.W. (1978) Digital processing of speech signals, Prentice Hall, Englewood Cliffs, p. 404-406.
  //  The result in formula:
  //  G^2=R_n(0)-Sum_{k=1}^{p} \a_{k} R_{n}(k)
  //  with 
  //  G gain factor
  //  R_{n}(k) k-th autocorrelation coefficient at time instant (i.e. sample index) n
  //  a_k LPC coefficients"
	for (int i = 0; i < LPC_Order; i++) {
		gain += LPC_Coef[i] * AutoCorr[i+1];
	}
	gain = sqrt((1.0 - gain) * originalZerothCorrCoeff); //1.0 is the normalized zeroth autocorr coefficient here
	//end block by thilo
    
	MFree_1D(Reflect);
	MFree_1D(AutoCorr);
	
	return gain; //by thilo
}

//=================================================================
//  Convert LPC coefficients to LPC based cepstral coefficients
//=================================================================
void GN_LPC::Lpc2Cep (double* LPCArray, int LpcOrder, double* Cep, int CepOrder) {

  int    CepCnt, Cnt;
  double sum;
 
  /*-----------------------------------------------------*/
  /* Calculate Complex cepstrum  with recursive relation */
  /* See. Rabiner's book                                  */
  /*-----------------------------------------------------*/
  for (CepCnt=0; CepCnt<CepOrder; CepCnt++) {

    //-------------------------------
	// First term in the formula    
	//-------------------------------
	if (CepCnt < LpcOrder) {sum = LPCArray[CepCnt];}
	else sum = 0.0;               /* CepOrder can be larger than LpcOrder */

    //-------------------------------
	// Second term in the formula    
	//-------------------------------
    for (Cnt=0; Cnt < CepCnt; Cnt++) {
     if (CepCnt-Cnt < LpcOrder)
	   sum += Cep[Cnt] * LPCArray[CepCnt-Cnt] * (Cnt+1) / (CepCnt+1);
    }

    Cep[CepCnt] = sum;
  }  /* CepCnt */

}


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "SV_Lib.h"


#define DELTA  0.00781250
#define BISECTIONS 4

#ifdef _MSC_VER
  //by thilo: 
  //M_PI is only defined in math.h (for ms-vc++) if _USE_MATH_DEFINES is defined beforehand (what is here not the case)
  //is that different under gcc?
  #define M_PI       3.141592653589793 
#endif

//==============================================================
// G - compute the value of the Chebychev series
//                sum c_k T_k(x) = x b_1(x) - b_2(x) + c_0
//                b_k(x) = 2x b_{k+1}(x) - b_{k+2}(x) + c_k 
//==============================================================
static double lsp_g(double x, double *c,int p2) {

	int i;
	double b[3];

    b[1] = b[2] = 0.0;

    for(i=p2; i > 0; i--) {

        b[0] = 2.0*x*b[1] - b[2] + c[i];
        b[2] = b[1];
        b[1] = b[0];
    }

    b[0] = x*b[1]-b[2]+c[0];
    return(b[0]);
} 

//==============================================================
// LSP_ROOTS
//        - find the roots of the two polynomials G_1(x) and G_2(x)
//          the first root corresponds to G_1(x)
//          compute the inverse cos (and these are the LSFs) 
//==============================================================
int lsp_roots(double *w, double **c, int p2) {

    int i,k;
    double x,x0,x1,y,*ptr,g0,g1;

    w[0] = 0.0;

    ptr = c[0];
    x = 1.0;
    g0 = lsp_g(x,ptr,p2);

    for(k=1,x = 1.0-DELTA; x > -DELTA-1.0; x -= DELTA) {

        /* Search for a zero crossing */
        if (g0*(g1 = lsp_g(x,ptr,p2)) <= 0.0) {

            /* Search Incrementally using bisection */
            x0 = x+DELTA;
            x1 = x;

            for(i=0; i < BISECTIONS; i++)
            {
                x = (x0+x1)/2.0;
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
            w[k] = (double)acos((double)x)/M_PI;

            ptr = c[k % 2];
            k++;
            if (k > 2*p2)
                return(0);

            g1 = lsp_g(x,ptr,p2);
        }
        g0 = g1;
    }

	REPORT_ERROR(SVLIB_Fail, "LSPs Not All Found");
    return(1);
} 

//==============================================================
//  Convert LPC coefficients to Line Spectrum Frequency (LSF)
//  Input    : LPC coefficients,  
//  Outnput  : LSF angle range is from 0 - 1, (divided by PI) 
//  both array has p+1 elements, w[0..p], a[0, p]; 
//
//  Reference:  Kabal and Ramachandran
//==============================================================
void GN_LPC::Lpc2Lsf(double *LpcCoef, double *LsfCoef, int LpcOrder) {
//void GN_LPC::Lpc2Lsf(double *a, double *w, int p) {


    int i, p, p2;
    double **c, *a, *w;

	//--------------------------------------------
	// append first element "1" before LpcCoef 
	//--------------------------------------------
	p = LpcOrder + 1;
	MArray_1D(a, p+5, double, "Lpc2Lsf");
	MArray_1D(w, p+5, double, "Lpc2Lsf");
	a[0] = 1.0;
	w[0] = 0.0;
	for (i=1; i<=p; i++) {
		a[i] = LpcCoef[i-1];
		w[i] = 0.0;
	}

    p2 = p/2;

    MArray_2D(c, 2, (p2+1), double, "LPC2LSF");

	c[0][p2] = c[1][p2] = 1.0;

    for(i=1; i <= p2; i++) {
        c[0][p2-i] = (a[i] + a[p+1-i] - c[0][p2+1-i]);
        c[1][p2-i] = c[1][p2+1-i] + a[i] - a[p+1-i];
    }

    c[0][0] /= 2.0;
    c[1][0] /= 2.0;

    i = lsp_roots(w, c, p2);

	//--------------------------------------------
	// discard first element "0" in w 
	//--------------------------------------------
	for (i=1; i<=p; i++) {
		LsfCoef[i-1] = w[i];
	}

	MFree_2D(c);
	MFree_1D(a);
	MFree_1D(w);

} /* LPC_PRED2LSP */



//==============================================================
//  Convert Line Spectrum Frequency to LPC coefficients
//  Input  : LSF angle (divided by PI, so the range is from 0 - 1) 
//  Output : LPC coefficients,  
//  both array has p+1 elements, w[0..p], a[0, p]; 
//
//  Reference:  Kabal and Ramachandran
//==============================================================
void GN_LPC::Lsf2Lpc(double *LsfCoef, double *LpcCoef, int LpcOrder) {
//void GN_LPC::Lsf2Lpc(double *w, double *a, int p) {

    int i,j,k,p,p2;
    double **f,c[2], *a, *w;

	//--------------------------------------------
	// append first element "0" at the begining  
	//--------------------------------------------
	p = LpcOrder + 1;
	MArray_1D(a, p+5, double, "Lsf2Lpc");
	MArray_1D(w, p+5, double, "Lsf2Lpc");
	w[0] = 0.0;
	for (i=1; i<=p; i++) {
		w[i] = LsfCoef[i-1];
	}


	a[0] = 1.0;
    p2 = p/2;
	MArray_2D(f, 2, (p2+1), double, "LSF2LPC");

    f[0][0] = f[1][0] = 1.0;
    f[0][1] = (double)-2.0*cos((double)w[1]*M_PI);
    f[1][1] = (double)-2.0*cos((double)w[2]*M_PI);

    k = 3;

    for(i=2; i <= p2; i++)
    {
        c[0] = (double)-2.0*cos((double)w[k++]*M_PI);
        c[1] = (double)-2.0*cos((double)w[k++]*M_PI);
        f[0][i] = f[0][i-2];
        f[1][i] = f[1][i-2];

        for(j=i; j >= 2; j--)
        {
            f[0][j] += c[0]*f[0][j-1]+f[0][j-2];
            f[1][j] += c[1]*f[1][j-1]+f[1][j-2];
        }
        f[0][1] += c[0]*f[0][0];
        f[1][1] += c[1]*f[1][0];
    }

    for(i=p2; i > 0; i--)
    {
        f[0][i] += f[0][i-1];
        f[1][i] -= f[1][i-1];

        a[i] = 0.50*(f[0][i]+f[1][i]);
        a[p+1-i] = 0.50*(f[0][i]-f[1][i]);
    }

	//--------------------------------------------
	// discard first element a(0) 
	//--------------------------------------------
	for (i=1; i<=p; i++) {
		LpcCoef[i-1] = a[i];
	}

    MFree_2D(f);
	MFree_1D(a);
	MFree_1D(w);
}

#ifdef TEST_LSF2LPC
//================================================
// test LPC <-> LSF
//================================================
void main (void ) {

	GN_LPC Eng;

	int Cnt;
    const int Order = 10;

	double LpcCoef[Order] = {-0.0161, -0.0435, -0.0371, -0.0336, -0.0612,
					-0.0158, -0.0485, 0.0316, -0.0920, -0.0260};

	double LsfCoef[Order], LPC_new[Order];

    Eng.Lpc2Lsf (LpcCoef, LsfCoef, Order);
	Eng.Lsf2Lpc (LsfCoef, LPC_new, Order);

	for (Cnt=0; Cnt<Order; Cnt++) {
	
		cout << LsfCoef[Cnt] << " " << LPC_new[Cnt] << endl;
	}    

}
#endif 
// TEST_LSF2LPC

//========================================================================
//  Calculate prediction residual of a signal with given LPC coefficients
//  R(n) = S(n) + a(0) S(n-1) + a(1)S(n-2) + ...+ a(k)S(n-k-1)
//  The return signal (Rst) should have same length as Sig and allocate
//  memory outside this procedure.
//========================================================================
void GN_LPC::Residual(float* Sig, float* Rst, int Num, double* LpcCoef, int LpcOrder) {


	int Ind;

	for (int SCnt=0; SCnt<Num; SCnt++) {
		Rst[SCnt] = Sig[SCnt];
		Ind = 0;
		while ( (SCnt>=Ind+1) & (Ind<LpcOrder)) {
			Rst[SCnt] += (float)LpcCoef[Ind] * Sig[SCnt-Ind-1];
			Ind++;
		}
	}
}
