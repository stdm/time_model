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
//    Several common matrix operations.
//
//
//    Author  : Jialong HE
//    Date    : March 25, 1999
//
//************************************************************************
#ifndef __GN_Matrix_H__
#define __GN_Matrix_H__

#include "SV_General.h"

//===================================================================
//  Some commonly used Matrix and Vector operations
//===================================================================
class GN_Matrix {

private :

protected:


public :

	GN_Matrix();
	virtual ~GN_Matrix();

	//-----------------------------------
	//  Return determinant of SrcMat
	//-----------------------------------
	double Det (double **SrcMat, int Dim);

	//-----------------------------------
	// Replace SrcMat with its inverse
	// return 0 if success, otherwise -1 
	//-----------------------------------
	int Inv (double **SrcMat, int Dim);

	//------------------------------------------------------
	// Find Symmetric Matrix Eigen value and Eigen Vectors
	// column eigen vectors put in **EigVec
	// eigen values in *EigVal
	//------------------------------------------------------
	void Eigen(double **SymMat, double **EigVec, double *EigVal, int Dim);
	void eigens( double A[], double RR[], double E[], int N ); //by thilo: made eigens() a class-member method

	//=========================================================
	//  Solve linear equation (Ax=B) by Gauss-Jordan elimination.
	//  A[0..n-1][0..n-1] is the input matrix. 
	//  B[0..n-1][0..m-1] is input containing the m right-hand side vectors 
	//  On output, A is replaced by its inverse matrix, 
	//  and B is replaced by the corresponding set of solution vectors.
	//=========================================================
	void lsolver(double **A, double **B, int n, int m);

	//-------------------------------------------------
	// Singular value decomposition (SVD)
	// 
	// Decompose any (m-by-n) matrix: Mat = U * S * V' 
	// 
	// U  is a m-by-m symmetric matrix (not returned)
	//
	// US is a m-by-n matrix of U*S
	// S  is a m-by-n diagonal matrix
	// V  is a n-by-n matrix
	//
	//  NOTE: program return U*S rather than U
	//  if (m != n), decompose might be not unique.
	//-------------------------------------------------
	void SVD (double **Mat, double **US, double **S, double **V, int m, int n);



};   // class GN_Matrix


#endif

