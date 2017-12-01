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
//    Train/Test VQ codebook (LBG algorithm).
//
//    Author  : Jialong He
//    Date    : May 7, 1999
//
//************************************************************************
#ifndef __SV_Model_VQ_H__
#define __SV_Model_VQ_H__

#include <iostream>
#include "SV_General.h"
#include "SV_Model.h"
#include "SV_Data.h"

using namespace std;

//===================================================================
//  Defination of VQ model
//===================================================================
class SV_Model_VQ : public SV_Model {

private :

protected :
	float **CBook;
	int	   CNum, Dim;

	float LBG_Codebook(float **codebook, int numcodes, float **data, int numdata, int dimension);
	int winning_cell(float **CodeBook, float *OneVec, int numcv, int dimension, float *Dist);
	void centroid (float **data, float *centvec, int *LabelVec, int numdata, int dimension, int Label);
	void Split_Codebook(float **OneBook, float **OneClass, int CodeNum, int OneNum, int DataDim);

public :
	//------------------------------- 
	// Training parameters
	//------------------------------- 
	int Verbose;
	int MaxIter;
	int SplitMethod;
	int CBSize;
	int RandSeed;

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Model_VQ();
	virtual ~SV_Model_VQ();

	//---------------------------------------- 
	// provide access to model's parameters 
	//----------------------------------------
	int GetCNum(void ) {return (CNum);}
	int GetCDim(void ) {return (Dim);}
	float **GetCBook(void) {return (CBook);}

	//------------------------------- 
	// public methods
	//-------------------------------
	virtual int SaveModel(void);
	virtual SV_Model *LoadModel(void);
	virtual int TrainModel(SV_Data *TrainData); //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	virtual SV_Data *TestModel(SV_Data *TestData);

	//----------------------------------------- 
	// dump model's parameter in ASCII
	//-----------------------------------------
	friend ostream& operator<< (ostream& os, SV_Model_VQ& Data);

};

#endif
