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
//    Train/Test continuous Density HMM.
//
//
//    Author  : Jialong HE
//    Date    : May 11, 1999
//    Modified: 29.04.2008 by thilo: added initial state distribution
//************************************************************************
#ifndef __SV_Model_CHMM_H__
#define __SV_Model_CHMM_H__

#include <iostream>
#include "SV_General.h"
#include "SV_Model.h"
#include "SV_Data.h"

/*----------------------------------------------*/
/* define a Gaussian Mixture structure          */
/*----------------------------------------------*/
typedef struct {
   double wgt;           /* component weight */
   double *mean;         /* mean vector of Gaussian function */
   double *vari;         /* variance of Gaussian function */
} MixType;


//===================================================================
//  CHMM model class
//===================================================================
class SV_Model_CHMM : public SV_Model {

private :

protected :

  //------------------------
	// Model's parameters
	//------------------------
	int   StaNum;                         /* States in the HMM */
	int   VecDim;                         /* vector dimension for Mean */
	int   MixNum;                         /* mixture number in one state */
	int   OrthT;                          /* if apply orthogonal transform */
	double **OMat;                        /* Orthogonal tran. Matrix */
	double **Tran;                        /* Hold transition prob. */
	MixType **Emit;                       /* Hold emision prob. */
  double *initialStateProbability;      //by thilo: initial state distribution (equation (9) in rabiner's paper)

	//--------------------------------
	// Other info from training data
	//--------------------------------
	double **CMat; 	                      // overall covariance matrix 
	double *MVec; 	                      // overall mean of training data 
	int VecNum;                           /* number of training vectors  */
	int PatNum;                           /* number of training patterns  */
	int MaxLen;                           /* Max length of training strings */ 
	double LMean;                         /* Likelihood mean and variance */
	double LVari;                         /* used for verification */

	void FreeModel(void);		            	// free dynamic allocate memory
	void AllocModel(void);			          // allocate memory for OMat, WgtVec, MeanMat, VariMat 
	void InitModel(bool isLeftRightModel = true); //by thilo: added parameter to drive initialization of initial state probabilities
	double viterbi(double **Tran, double **MixOut, int* Path, int Len);
	void emission(double **MixOut, MixType **Emit, float **Data, int SeqLen);
	double forward(double **Alpha, double *Scale, double **Tran, double **MixOut, int Len);
	void backward(double **Beta, double *Scale, double **Tran, double **MixOut, int Len);

public :

	//------------------------------------------
	// used for setting training parameters
	//------------------------------------------
	int Verbose;
	int Mixtures;
	int WithOrth;
	int MaxIter;
	int RandSeed;
	int NState;
	int **ConfMat;					// hold HMM structure (for Init Tran)
	bool isLeftRightModel;  // by thilo: if true, the model is considered to be a left-2 right model (i.e. the initial state prob. of state 0 is 1, all others are 0; all other constraints have to be communicated via the ConfMat structure, i.e. max. jump size to avoid large changes in state indices, see. rabiner, equation (45-48b))

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Model_CHMM();
	virtual ~SV_Model_CHMM();

	//------------------------------- 
	// access protected members
	//-------------------------------
	int   GetStaNum() {return(StaNum);}
	int   GetVecDim() {return(VecDim);}
	int   GetMixNum() {return(MixNum);}
	double **GetTran() {return(Tran);}
	MixType **GetEmit() {return(Emit);}
	double *getInitialStateProbability() {return this->initialStateProbability;} //by thilo
	double doViterbi(double **mixtureOutput, int *bestPath, int pathLength) {return viterbi(this->Tran, mixtureOutput, bestPath, pathLength);} //by thilo
	void getMixtureEmissions(double **mixtureOutput, float **features, int featureCount) {return emission(mixtureOutput, this->Emit, features, featureCount);} //by thilo

	//------------------------------- 
	// public methods
	//-------------------------------
	virtual int SaveModel(void);
	virtual SV_Model *LoadModel(void);
	virtual int TrainModel(SV_Data *TrainData); //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	virtual SV_Data *TestModel(SV_Data *TestData);

	//------------------------------- 
	// dump model's parameter for Debug
	//-------------------------------
	friend ostream& operator<< (ostream& os, SV_Model_CHMM& Data);

};

#endif
