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
//    Train/Test discrete Density HMM.
//
//
//    Author  : Jialong HE
//    Date    : May 13, 1999
//************************************************************************
#ifndef __SV_Model_DHMM_H__
#define __SV_Model_DHMM_H__

#include <iostream>
#include "SV_General.h"
#include "SV_Model_VQ.h"
#include "SV_Data.h"

using namespace std;

const double MIN_PROB = 0.000001;	// never let alpha or beta equal zero

//=====================================
// State Class
//=====================================
class state {

public:
  double* recur_out; //  output prob. for each symbol when self-loop
  double* next_out;  //  output prob. for each symbol when jump to next state
  double recur_trans;
  double next_trans;
  int SymNum;

  state(int num_symbols);
 ~state();
  double set_recur_out(int symbol, double prob=-1.0);
  double set_next_out(int symbol, double prob=-1.0);
  double set_recur_trans(double prob=-1.0);
  double set_next_trans(double prob=-1.0);

}; 


//===================================================================
//  DHMM model, the base model is VQ
//===================================================================
class SV_Model_DHMM : public SV_Model_VQ {

private :
	                
protected :

    //------------------------
	// Model's parameters
	//------------------------
	state** States;     // column vector of states (one col matrix)
	int   StaNum;       /* States in the HMM */

	//------------------------------------------------
	// Convert vector sequences to discrete symbols
	// with current Codebook.
	//------------------------------------------------
	int trellis_width;
	double** alpha;                   // matrix for alpha trellis
	double** beta;                    // matrix for beta trellis
	double** gamma_next;              // matrix of gamma transition probs
	double** gamma_recur;             // matrix of gamma recurrence probs
	double*  scaling_factors;         // array of alpha and beta scaling factors
	double* a_numer_sum_recur;        // array of numerators for a_ij's
	double* a_numer_sum_next;         // array of numerators for a_ij's
	double* a_denom_sum_recur;        // array of denomonators for a_ij's
	double* a_denom_sum_next;         // array of denomonators for a_ij's
	double** b_numer_sum_recur;       // array of numerators for b_ij's
	double** b_numer_sum_next;        // array of numerators for b_ij's

	short **strings;
	int *string_len;
	int num_strings;
	int SymNum;
	int GenSymbol(SV_Data* pData);
	
	void FreeModel(void);			// free dynamic allocate memory
	void AllocModel(void);			// allocate memory for state
	void InitModel(void);			// allocate memory for state
	double BatchTrain(void);
	double set_cumulative_ab_counts(void);
	void rescale_alphas(int col);
	void rescale_betas(int col);
	double alpha_F(short* symbol_array, int symbol_count );
	double beta_I(short* symbol_array, int symbol_count);
	void compute_gamma(short* symbol_array, int symbol_count );
	double a_numer(int i, int j, int symbol_count);
	double a_denom(int i, int j, int symbol_count );
	double b_numer(int i, int j, int sym, short *symbol_array, int symbol_count );
	double reestimate(void);
	double test(short* string, int symbol_count );
	
public :

	int NState;
	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	SV_Model_DHMM();
	virtual ~SV_Model_DHMM();

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
	friend ostream& operator<< (ostream& os, SV_Model_DHMM& Data);

};

#endif

