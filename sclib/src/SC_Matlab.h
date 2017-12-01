/**************************************************************************/
/*	This class encapsulates communication with Matlab compiled libraries  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 11.08.2009																								*/
/**************************************************************************/

#ifndef __SC_Matlab_H__
#define __SC_Matlab_H__

#include "SC_Api.h"
#include <SV_Data.h>

typedef void* MatlabRealMatrix;

class SC_Matlab {

	private :

	protected:
		static bool matlabIsInitialized; //true after Matlab has been initialized for the first time

	public :

    //====================================================================================================================
		// constructor/destructor
    //====================================================================================================================
	  SC_Matlab();
		virtual ~SC_Matlab();

    //====================================================================================================================
		// initializes Matlab, if not already done
    //====================================================================================================================
		static bool initialize(void);

    //====================================================================================================================
		// terminates Matlab (only allowed once at runtime, better do never...)
    //====================================================================================================================
		static bool terminate(void);

    //====================================================================================================================
		// allocates memory for the matlab data structures (mxArray) and copies the values; 
		// results must be casted to mxArray*
    //====================================================================================================================
		MatlabRealMatrix scalar2matlab(double scalar);
		MatlabRealMatrix vector2matlab(double *vector, unsigned long int dim, bool rowVector = true);
		MatlabRealMatrix vector2matlab(float *vector, unsigned long int dim, bool rowVector = true);
		MatlabRealMatrix matrix2matlab(double **matrix, unsigned long int len, unsigned long int dim);
		MatlabRealMatrix matrix2matlab(float **matrix, unsigned long int len, unsigned long int dim);

		//====================================================================================================================
		// receives a Matlab matrix and returns its values in a newly created c++ object/matrix/array/scalar
		// parameter must be a mxArray* casted to void*
		//====================================================================================================================
		SV_Data* matlab2svdata(MatlabRealMatrix m, bool pivot = false);
		double** matlab2matrix(MatlabRealMatrix m, bool pivot = false);
		double* matlab2array(MatlabRealMatrix m);
		double matlab2scalar(MatlabRealMatrix m);

		//====================================================================================================================
		// releases the memory behind the Matlab matrix m
		//====================================================================================================================
		void freeMatrix(MatlabRealMatrix m);

		//====================================================================================================================
		// wrappers around functions available in libMatlab
		//====================================================================================================================
		MatlabRealMatrix sc_mlfL1qc_logbarrier(MatlabRealMatrix x0, MatlabRealMatrix A, MatlabRealMatrix b, MatlabRealMatrix epsilon, MatlabRealMatrix lbtol, MatlabRealMatrix mu);
		MatlabRealMatrix sc_mlfMfcc(MatlabRealMatrix input, MatlabRealMatrix samplingRate, MatlabRealMatrix frameRate);
		//bool mlfMelfcc(mxArray** cepstra, mxArray** aspectrum, mxArray** pspectrum, mxArray* samples, mxArray* sr, mxArray* varargin);
		MatlabRealMatrix sc_mlfL1eq_pd(MatlabRealMatrix x0, MatlabRealMatrix A, MatlabRealMatrix b);
};


#endif
