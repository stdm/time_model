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
//    Train/Test Orthogonal Gaussian Mixture Model (OGMM).
//
//
//    Author  : Jialong HE
//    Date    : May 11, 1999
//************************************************************************
#ifndef __SV_Model_GMM_H__
#define __SV_Model_GMM_H__

#include <iostream>
#include "SV_General.h"
#include "SV_Model.h"
#include "SV_Data.h"

//===================================================================
//  Defination of OGMM model
//===================================================================
class SV_Model_GMM : public SV_Model {

private :

protected :
	//------------------------
	// Model's parameters
	//------------------------
	int   VecDim;                         /* current dim after mapping */
	int   MixNum;                         /* mixture number in GC */
	int   OrthT;                          /* if apply orthogonal transform */
	double **MeanMat;                     /* Mean vector of Gaus comp. */
	double **VariMat;                     /* Variance vector of Gaus Comp. */
	double *WgtVec;                       /* Mixture weight vector */
	double **OMat;                        /* Orthogonal tran. Matrix */


	//--------------------------------
	// Other info from training data
	//--------------------------------
	double **CMat; 	                      // overall covariance matrix 
	double *MVec; 	                      // overall mean of training data 

	void FreeModel(void);      // free dynamic allocate memory
	void AllocModel(void);	   // allocate memory for OMat, WgtVec, MeanMat, VariMat 
	int  Split(int CurrMixNum);
	void reestimate(float **DataMatrix, int DataNum, int CurrCode);
	double Log_Likelihood(float *TestVec, int Dim);

public :

	//------------------------------- 
	// Used for setting training par.
	//------------------------------- 
	int Verbose;
	int Mixtures;
	int WithOrth;
	int MaxIter;
	int SplitMethod;
	int RandSeed;

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Model_GMM();
	virtual ~SV_Model_GMM();

	//----------------------------------------- 
	// provide access to protected members
	//----------------------------------------- 
	double **GetMeanMat() {return(MeanMat);} 
	double **GetVariMat() {return(VariMat);} 
	double *GetWgtVec() {return(WgtVec);}
	int   GetVecDim() {return(VecDim);}
	int   GetMixNum(){return(MixNum);} 

	//------------------------------- 
	// public methods
	//-------------------------------
	virtual int SaveModel(void);
	virtual SV_Model *LoadModel(void);
	virtual int TrainModel(SV_Data *TrainData); //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	virtual SV_Data *TestModel(SV_Data *TestData);

	//------------------------------- 
	// dump model's parameters in ASCII
	//-------------------------------
	friend ostream& operator<< (ostream& os, SV_Model_GMM& Data);

};

#endif
