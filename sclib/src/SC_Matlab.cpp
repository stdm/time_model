/**************************************************************************/
/*	This class encapsulates communication with Matlab compiled libraries  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 11.08.2009																								*/
/**************************************************************************/

#include "SC_Matlab.h"
#include <SV_Error.h>
#ifdef SC_USE_MATLAB
	#include <libMatlab.h>
#endif

bool SC_Matlab::matlabIsInitialized = false;

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Matlab::SC_Matlab() {

}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Matlab::~SC_Matlab() {

}

//====================================================================================================================
// initializes Matlab, if not already done
//====================================================================================================================
bool SC_Matlab::initialize(void) {
#ifdef SC_USE_MATLAB
	if (SC_Matlab::matlabIsInitialized == false) {
		char *pStrings[]={"-nojvm","-nojit"};
		if (!mclInitializeApplication(const_cast<const char**>(pStrings), 2)) {
			REPORT_ERROR(SVLIB_Fail, "Could not initialize MCR for the application.\n");
			return false;
		}

		if(!libMatlabInitialize()) {
			REPORT_ERROR(SVLIB_Fail, "Could not initialize the Matlab-derived library.\n");
			return false;
		}

		SC_Matlab::matlabIsInitialized = true;
	}

	return true;
#else
	return false;
#endif
}

//====================================================================================================================
// terminates Matlab (only allowed once at runtime, better do never...)
//====================================================================================================================
bool SC_Matlab::terminate(void) {
#ifdef SC_USE_MATLAB
	if (SC_Matlab::matlabIsInitialized == true) {
		libMatlabTerminate();

		if (!mclTerminateApplication()) {
			return false;
		}

		SC_Matlab::matlabIsInitialized = false;
	}
	
	return true;
#else
	return false;
#endif
}

//====================================================================================================================
// allocates memory for the matlab data structure (1*1 mxArray) and copies the values
// result must be casted to mxArray*
//====================================================================================================================
MatlabRealMatrix SC_Matlab::scalar2matlab(double scalar) {
#ifdef SC_USE_MATLAB
	mxArray *m = mxCreateDoubleMatrix(1, 1, mxREAL);
	memcpy(mxGetPr(m), &scalar, 1*sizeof(double));

	return (MatlabRealMatrix)(m);
#else
	return NULL;
#endif
}



//====================================================================================================================
// allocates memory for the matlab data structure (1*dim mxArray in case of rowVector==true) and copies the values
// result must be casted to mxArray*
//====================================================================================================================
MatlabRealMatrix SC_Matlab::vector2matlab(double *vector, unsigned long int dim, bool rowVector) {
#ifdef SC_USE_MATLAB
	mxArray *m;
	if (rowVector == true) {
		m = mxCreateDoubleMatrix(1, dim, mxREAL);
	} else {
		m = mxCreateDoubleMatrix(dim, 1, mxREAL);
	}
	memcpy(mxGetPr(m), vector, dim*sizeof(double));

	return (MatlabRealMatrix)(m);
#else
	return NULL;
#endif
}

//====================================================================================================================
// allocates memory for the matlab data structure (1*dim mxArray in case of rowVector==true) and copies the values
// result must be casted to mxArray*
//====================================================================================================================
MatlabRealMatrix SC_Matlab::vector2matlab(float *vector, unsigned long int dim, bool rowVector) {
#ifdef SC_USE_MATLAB
	mxArray *m;
	if (rowVector == true) {
		m = mxCreateDoubleMatrix(1, dim, mxREAL);
	} else {
		m = mxCreateDoubleMatrix(dim, 1, mxREAL);
	}
	
	double *values = mxGetPr(m);
	for (unsigned long int x = 0; x < dim; x++) {
		values[x] = vector[x];
	}	

	return (MatlabRealMatrix)(m);
#else
	return NULL;
#endif
}

//====================================================================================================================
// allocates memory for the matlab data structure (len*dim mxArray) and copies the values
// result must be casted to mxArray*
//====================================================================================================================
MatlabRealMatrix SC_Matlab::matrix2matlab(double **matrix, unsigned long int len, unsigned long int dim) {
#ifdef SC_USE_MATLAB
	mxArray *m = mxCreateDoubleMatrix(len, dim, mxREAL);
	double *values = mxGetPr(m);
	
	for (unsigned long int y = 0; y < len; y++) {
		for (unsigned long int x = 0; x < dim; x++) {
			values[len*x + y] = matrix[y][x];
		}
	}

	return (MatlabRealMatrix)(m);
#else
	return NULL;
#endif
}

//====================================================================================================================
// allocates memory for the matlab data structure (len*dim mxArray) and copies the values
// result must be casted to mxArray*
//====================================================================================================================
MatlabRealMatrix SC_Matlab::matrix2matlab(float **matrix, unsigned long int len, unsigned long int dim) {
#ifdef SC_USE_MATLAB
	mxArray *m = mxCreateDoubleMatrix(len, dim, mxREAL);
	double *values = mxGetPr(m);
	
	for (unsigned long int y = 0; y < len; y++) {
		for (unsigned long int x = 0; x < dim; x++) {
			values[len*x + y] = (double)(matrix[y][x]);
		}
	}

	return (MatlabRealMatrix)(m);
#else
	return NULL;
#endif
}

//====================================================================================================================
// receives a Matlab matrix and returns its values in a newly created SV_Data object
// parameter must be a mxArray* casted to void*
//====================================================================================================================
SV_Data* SC_Matlab::matlab2svdata(MatlabRealMatrix m, bool pivot) {
#ifdef SC_USE_MATLAB
	SV_Data *pRes = NULL;
	mxArray *input = (mxArray*)(m);

	if (input != NULL) {
		double *values = mxGetPr(input);
		int len, dim;

		if (pivot == false) {
			len = (int)(mxGetM(input));
			dim = (int)(mxGetN(input)); //mxGetN() return numbers of columns
		} else {
			len = (int)(mxGetN(input));
			dim = (int)(mxGetM(input));
		}

		if (len>0 && dim>0) {
			pRes = new SV_Data(len, dim);
			for(int y = 0; y < len; y++) {
				for(int x = 0; x < dim; x++) {
					if (pivot == false) {
						pRes->Mat[y][x] = (float)(values[len*x + y]); //experimentally proven that this yields the correct pivot
					} else {
						pRes->Mat[y][x] = (float)(values[dim*y + x]); //experimentally proven that this yields the correct pivot
					}
				}
			}
		}
	}

	return pRes;
#else
	return NULL;
#endif
}

//====================================================================================================================
// receives a Matlab matrix and returns its values in a newly created double matrix
// parameter must be a mxArray* casted to void*
//====================================================================================================================
double** SC_Matlab::matlab2matrix(MatlabRealMatrix m, bool pivot) {
#ifdef SC_USE_MATLAB
	mxArray *input = (mxArray*)(m);
	double **res = NULL;

	if (input != NULL) {
		double *values = mxGetPr(input);
		int len, dim;

		if (pivot == false) {
			len = (int)(mxGetM(input));
			dim = (int)(mxGetN(input)); //mxGetN() return numbers of columns
		} else {
			len = (int)(mxGetN(input));
			dim = (int)(mxGetM(input));
		}

		if (len>0 && dim>0) {
			MArray_2D(res, len, dim, double, "SC_Matlab.matlab2matrix: res");
			for(int y = 0; y < len; y++) {
				for(int x = 0; x < dim; x++) {
					if (pivot == false) {
						res[y][x] = (float)(values[len*x + y]); //experimentally proven that this yields the correct pivot
					} else {
						res[y][x] = (float)(values[dim*y + x]); //experimentally proven that this yields the correct pivot
					}
				}
			}
		}
	}

	return res;
#else
	return NULL;
#endif
}

//====================================================================================================================
// receives a Matlab matrix and returns its values in a newly created double array
// parameter must be a mxArray* casted to void*
//====================================================================================================================
double* SC_Matlab::matlab2array(MatlabRealMatrix m) {
#ifdef SC_USE_MATLAB
	mxArray *input = (mxArray*)(m);
	double *res = NULL;

	if (input != NULL) {
		double *values = mxGetPr(input);
		int len = (int)(mxGetM(input)), dim = (int)(mxGetN(input));
		
		if (len == 1) {
			MArray_1D(res, dim, double, "SC_Matlab::matlab2array: res");
			for (int i = 0; i < dim; i++) {
				res[i] = values[i];
			}
		} else if (dim == 1) {
			MArray_1D(res, len, double, "SC_Matlab::matlab2array: res");
			for (int i = 0; i < len; i++) {
				res[i] = values[i];
			}
		}
	}

	return res;
#else
	return NULL;
#endif
}

//====================================================================================================================
// receives a Matlab matrix and returns its values in a newly created double scalar
// parameter must be a mxArray* casted to void*
//====================================================================================================================
double SC_Matlab::matlab2scalar(MatlabRealMatrix m) {
#ifdef SC_USE_MATLAB
	mxArray *input = (mxArray*)(m);
	double *values = mxGetPr(input);

	return (m!=NULL) ? values[0] : -1.0;
#else
	return -1.0;
#endif
}

//====================================================================================================================
// releases the memory behind the Matlab matrix m
//====================================================================================================================
void SC_Matlab::freeMatrix(MatlabRealMatrix m) {
#ifdef SC_USE_MATLAB
	mxDestroyArray((mxArray*)(m));
#endif

	return;
}

//====================================================================================================================
// calls Slaney's MFCC extraction routine with only the cepstra as return values
//====================================================================================================================
MatlabRealMatrix SC_Matlab::sc_mlfMfcc(MatlabRealMatrix input, MatlabRealMatrix samplingRate, MatlabRealMatrix frameRate) {
#ifdef SC_USE_MATLAB
	mxArray *ceps = NULL;
	bool res = mlfMfcc(1, &ceps, NULL, NULL, NULL, NULL, (mxArray*)(input), (mxArray*)(samplingRate), (mxArray*)(frameRate)); 

	if (res == false) {
		mxDestroyArray(ceps);
		ceps = NULL;
	}

	return (MatlabRealMatrix)(ceps);
#else
	return NULL;
#endif
}

//====================================================================================================================
// calls Romberg's l1-Magic qc-logbarrier solver and returns the estimated sparse vector
//====================================================================================================================
MatlabRealMatrix SC_Matlab::sc_mlfL1qc_logbarrier(MatlabRealMatrix x0, MatlabRealMatrix A, MatlabRealMatrix b, MatlabRealMatrix epsilon, MatlabRealMatrix lbtol, MatlabRealMatrix mu) {
#ifdef SC_USE_MATLAB
	mxArray *xp = NULL;
	mxArray *tmp1 = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxArray *tmp2 = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxArray *tmp3 = mxCreateDoubleMatrix(1, 1, mxREAL);

	bool res = mlfL1qc_logbarrier(1, &xp, (mxArray*)(x0), (mxArray*)(A), tmp1, (mxArray*)(b), (mxArray*)(epsilon), (mxArray*)(lbtol), (mxArray*)(mu), tmp2, tmp3);

	mxDestroyArray(tmp1);
	mxDestroyArray(tmp2);
	mxDestroyArray(tmp3);

	if (res == false) {
		mxDestroyArray(xp);
		xp = NULL;
	}

	return (MatlabRealMatrix)(xp);
#else
	return NULL;
#endif
}

//====================================================================================================================
// calls Romberg's l1-Magic eq-pd solver and returns the estimated sparse vector
//====================================================================================================================
MatlabRealMatrix SC_Matlab::sc_mlfL1eq_pd(MatlabRealMatrix x0, MatlabRealMatrix A, MatlabRealMatrix b) {
#ifdef SC_USE_MATLAB
	mxArray *xp = NULL;
	mxArray *tmp1 = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxArray *pdtol = (mxArray*)(scalar2matlab(1e-3));
	mxArray *pdmaxiter = (mxArray*)(scalar2matlab(50));
	mxArray *cgtol = (mxArray*)(scalar2matlab(1e-8));
	mxArray *cgmaxiter = (mxArray*)(scalar2matlab(200));

	bool res = mlfL1eq_pd(1, &xp, (mxArray*)(x0), (mxArray*)(A), tmp1, (mxArray*)(b), pdtol, pdmaxiter, cgtol, cgmaxiter);

	mxDestroyArray(tmp1);
	mxDestroyArray(pdtol);
	mxDestroyArray(pdmaxiter);
	mxDestroyArray(cgtol);
	mxDestroyArray(cgmaxiter);

	if (res == false) {
		mxDestroyArray(xp);
		xp = NULL;
	}

	return (MatlabRealMatrix)(xp);
#else
	return NULL;
#endif
}
