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
//    Train/Test uni-modal multi-variate Gauss model.
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//
//************************************************************************
#ifndef __SV_Model_Gaus_H__
#define __SV_Model_Gaus_H__

#include <iostream>
#include "SV_General.h"
#include "SV_Model.h"
#include "SV_Data.h"

using namespace std;

//===================================================================
//  Defination of Gaus model
//===================================================================
class SV_Model_Gaus : public SV_Model {

private :

protected :
	double **ICov;
	double *MVec;
	double DetV;
	int	   Dim;

	double AvgL(double **Cx, double *Mx, double **Cy, double *My, int Dim);

public :

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Model_Gaus();
	virtual ~SV_Model_Gaus();

	//------------------------------- 
	// Model IO methods
	//-------------------------------
	virtual int SaveModel(void);
	virtual SV_Model *LoadModel(void);

	//------------------------------- 
	// Train and Test This model
	//-------------------------------
	virtual int TrainModel(SV_Data *TrainData); //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	virtual SV_Data *TestModel(SV_Data *TestData);

	//------------------------------- 
	// dump model's parameter for Debug
	//-------------------------------
	friend ostream& operator<< (ostream& os, SV_Model_Gaus& Data);

	int getDim(void) {return this->Dim;}; //by thilo
	double** getInvertedCov(void) {return this->ICov;}; //by thilo
	double* getMean(void) {return this->MVec;}; //by thilo
};

#endif
