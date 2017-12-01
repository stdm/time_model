/*
 * MATLAB Compiler: 4.9 (R2008b)
 * Date: Tue Aug 11 11:03:40 2009
 * Arguments: "-B" "macro_default" "-W" "lib:libMatlab" "-d"
 * "D:\Data\Matlab\libMatlab\src" "-T" "link:lib" "-v"
 * "D:\Data\Matlab\l1magic-1.1\Optimization\l1qc_logbarrier.m"
 * "D:\Data\Matlab\AuditoryToolbox\mfcc.m" "D:\Data\Matlab\rastamat\melfcc.m"
 * "D:\Data\Matlab\l1magic-1.1\Optimization\l1eq_pd.m" 
 */

#ifndef __libMatlab_h
#define __libMatlab_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libMatlab
#define PUBLIC_libMatlab_C_API __global
#else
#define PUBLIC_libMatlab_C_API /* No import statement needed. */
#endif

#define LIB_libMatlab_C_API PUBLIC_libMatlab_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libMatlab
#define PUBLIC_libMatlab_C_API __declspec(dllexport)
#else
#define PUBLIC_libMatlab_C_API __declspec(dllimport)
#endif

#define LIB_libMatlab_C_API PUBLIC_libMatlab_C_API


#else

#define LIB_libMatlab_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libMatlab_C_API 
#define LIB_libMatlab_C_API /* No special import/export declaration */
#endif

extern LIB_libMatlab_C_API 
bool MW_CALL_CONV libMatlabInitializeWithHandlers(mclOutputHandlerFcn error_handler,
                                                  mclOutputHandlerFcn print_handler);

extern LIB_libMatlab_C_API 
bool MW_CALL_CONV libMatlabInitialize(void);

extern LIB_libMatlab_C_API 
void MW_CALL_CONV libMatlabTerminate(void);



extern LIB_libMatlab_C_API 
void MW_CALL_CONV libMatlabPrintStackTrace(void);


extern LIB_libMatlab_C_API 
bool MW_CALL_CONV mlxL1qc_logbarrier(int nlhs, mxArray *plhs[],
                                     int nrhs, mxArray *prhs[]);

extern LIB_libMatlab_C_API 
bool MW_CALL_CONV mlxMfcc(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_libMatlab_C_API 
bool MW_CALL_CONV mlxMelfcc(int nlhs, mxArray *plhs[],
                            int nrhs, mxArray *prhs[]);

extern LIB_libMatlab_C_API 
bool MW_CALL_CONV mlxL1eq_pd(int nlhs, mxArray *plhs[],
                             int nrhs, mxArray *prhs[]);

extern LIB_libMatlab_C_API 
long MW_CALL_CONV libMatlabGetMcrID() ;



extern LIB_libMatlab_C_API bool MW_CALL_CONV mlfL1qc_logbarrier(int nargout
                                                                , mxArray** xp
                                                                , mxArray* x0
                                                                , mxArray* A
                                                                , mxArray* At
                                                                , mxArray* b
                                                                , mxArray* epsilon
                                                                , mxArray* lbtol
                                                                , mxArray* mu
                                                                , mxArray* cgtol
                                                                , mxArray* cgmaxiter);

extern LIB_libMatlab_C_API bool MW_CALL_CONV mlfMfcc(int nargout, mxArray** ceps
                                                     , mxArray** freqresp
                                                     , mxArray** fb
                                                     , mxArray** fbrecon
                                                     , mxArray** freqrecon
                                                     , mxArray* input
                                                     , mxArray* samplingRate
                                                     , mxArray* frameRate);

extern LIB_libMatlab_C_API bool MW_CALL_CONV mlfMelfcc(int nargout
                                                       , mxArray** cepstra
                                                       , mxArray** aspectrum
                                                       , mxArray** pspectrum
                                                       , mxArray* samples
                                                       , mxArray* sr
                                                       , mxArray* varargin);

extern LIB_libMatlab_C_API bool MW_CALL_CONV mlfL1eq_pd(int nargout
                                                        , mxArray** xp
                                                        , mxArray* x0
                                                        , mxArray* A
                                                        , mxArray* At
                                                        , mxArray* b
                                                        , mxArray* pdtol
                                                        , mxArray* pdmaxiter
                                                        , mxArray* cgtol
                                                        , mxArray* cgmaxiter);

#ifdef __cplusplus
}
#endif

#endif
