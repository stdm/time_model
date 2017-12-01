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
//    Interface to Error Handling Functions
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_Error_H__
#define __SV_Error_H__

#define REPORT_ERROR(ErrorCode, ErrorMsg) {(*FunPoint)(ErrorCode, ErrorMsg,  __FILE__, __LINE__);}
void SV_SetErrorHandler (void (*UserHandler)(int ErrorCode, const char* ErrorMsg, const char* FName, int LName) ); //by thilo: added "const" specifier to avoid warnings
void SV_DefaultHandler(int ErrorCode, const char* ErrorMsg, const char* FName, int LName); //by thilo: added "const" specifier to avoid warnings
extern void (*FunPoint)(int ErrorCode, const char* ErrorMsg, const char* FName, int LName); //by thilo: added "const" specifier to avoid warnings

//=========================================================
//  Define ErrorCode
//=========================================================
#define SVLIB_Ok      0
#define SVLIB_Fail    -1
#define SVLIB_NoMem   -2
#define SVLIB_BadArg  -3
#define SVLIB_DivBy0  -4
#define SVLIB_NoPara  -5
#define SVLIB_BadData -6
#define SVLIB_FileErr -7
#define SVLIB_NoInv   -8



//------------------------------------------------------
//  Allocate/Release 1D and 2D Array (Matrix)
//
//   MArray_1D(Target, Dim, type, Msg);
//   MArray_2D(Target, Dim1, Dim2, type, Msg);
//   MFree_1D(Target)
//   MFree_2D(Target)
//
//   Example:
//     float **matrix;
//     MArray_2D(matrix, 10, 10, float, "memory for matrix");
//       .... (access matrix[][]) .....
//     MFree_2D(matrix);
//
//    NOTE: unlike my C version which use calloc()
//          this version use new operator, array is not 
//          clear to ZERO
//------------------------------------------------------
//   Author : Jialong He
//   Date   : June 20, 1994
//   Modify : April 13, 1995  add Macro version
//------------------------------------------------------

//------------//
// 1-D Array  //
//------------//
#define MArray_1D(Target, Dim, Type, Msg) { \
	    if (((Dim) <= 0) || (Target = new Type[Dim]) == NULL) \
		{(*FunPoint)(SVLIB_NoMem, Msg, __FILE__, __LINE__);} }

//----------------------------------------------
// for an array of Object, must use delete []
//----------------------------------------------
#define MFree_1D(Target) {delete [] (Target); Target = NULL;}


//------------//
// 2-D Array  //
//------------//
#define MArray_2D(Target, Dim1, Dim2, Type, Msg) {  \
	    Type *_Pnt; int _Cnt;  \
		if ( (Dim1)<=0 || (Dim2) <= 0||(_Pnt = new Type[(Dim1)*(Dim2)])==NULL) {(*FunPoint)(SVLIB_NoMem, Msg, __FILE__, __LINE__);} \
		if ( (Target = new Type* [Dim1]) == NULL) {(*FunPoint)(SVLIB_NoMem, Msg, __FILE__, __LINE__);} \
        for(_Cnt = 0; _Cnt < (Dim1) ; _Cnt++) {Target[_Cnt] = _Pnt + _Cnt * (Dim2);}  \
	 }

#define MFree_2D(Target) { if (Target != NULL) {delete [] (*Target);  delete [] (Target); Target = NULL;}}

#endif

