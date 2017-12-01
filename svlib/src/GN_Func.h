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
//    This class provides sevearl algorithms to deal with 
//    single variable function,   y = f(x).
//
//    Author  : Jialong HE
//    Date    : May 17, 1999
//
//************************************************************************
#ifndef __GN_Func_H__
#define __GN_Func_H__

#include "SV_General.h"

//===================================================================
//  This class provides some method to math functions
//===================================================================
class GN_Func {

private :

protected:


public :

	GN_Func();
	virtual ~GN_Func();
  
	//--------------------------------------------------------------
	//  Obtain the machine's EPSILON
	//  i.e. the smallest positive number which, been added to 1., 
	//  yields the result other than 1.
	//  
	//  1.0 + EPSILON != 1.0
	//  
	//  EPSILON = 2.22e-16
	//--------------------------------------------------------------
	double Epsilon(void);

	//--------------------------------------------------------------
	//  Find a x value, so that f(x) = 0 within range (x1, x2).
	//  tol is toleranc value, default is 0.001
	//  Example:  Root = Fzero(sin, -1, 1);
	//--------------------------------------------------------------
	double Fzero(double (*func)(double x), double x1, double x2, double tol = 0.001);

	//--------------------------------------------------------------
	//  Find x value within the range (x1, x2) so that f(x) is minimized.
	//  tol is toleranc value, default is 0.001
	//  Example:  Xmin = Fmin(sin, 2, 6);
	//--------------------------------------------------------------
	double Fmin (double (*func)(double x), double x1, double x2, double tol = 0.001);       	        /* An estimate to the min location*/ 

	//--------------------------------------------------------------
	//  Numerically integrate f(x) from a to b
	//  
	//  Example:  Val = Finte(sin, 0, 3.14);
	//--------------------------------------------------------------
	double Finte(double (*func)(double x) , double x1, double x2);
	
	//--------------------------------------------------------------
	//   Calculate error function (Gaussian integration)
	//
	//                           x                            
	//                           -                            
	//                 2         | |          2              
	//   erf(x)  =  --------     |    exp( - t  ) dt.        
	//              sqrt(pi)   | |                           
	//                          -                            
	//	                         0                           	
	//                                                       
	//  erfc(x) = 1 - erf(x)
	//  inverf(y) find x that y = erf(x)
	//  
	//  NOTE : erf(x) calcuated from poly expansion
	//--------------------------------------------------------------
	double erf( double x);
	double erfc( double x);
	double inverf( double y);


	//===========================================================
	// Function fitting
	//  
	// Given a set of data points X[1..ndata], Y[1..ndata]  
	// with individual standard deviations W[1..ndata].
	// find fitting coefficients Coef for the approximation
	// Y = f(X) with
	//
	//  Y= Sum( Coef(i) * afunc(i, X) )
	//
	//  Funcs(x, afunc, NCoef), user provided functions.
	//  return NCoef base fucntion values in "afunc"
	//   
	//===========================================================
	double FitFuncs (double *X, double *Y, double *W, int ndata,
							double *Coef, int NCoef,
 						    void (*Funcs)(double a, double* afunc, int ma));

	//===========================================================
	//  Find all roots of a polynomial
	//
	//  Polynomial is given sum (C(i) * x^i), i = 0..M;
	//  M roots in Roots[0..M-1].
	//
	//  NOTE *C has M+1 elements from 0..M   
	//===========================================================
	void PolyRoots (COMPLEX *C, COMPLEX *Roots, int M);


};   // class GN_Func


#endif

