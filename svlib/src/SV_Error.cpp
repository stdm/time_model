//************************************************************************
//    Implement default error handling procedure.
//
//
//    Author  : Jialong HE
//    Date    : March 14, 1999
//************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include "SV_Error.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==================================================
// This is the default Error Handler Function
//==================================================
void SV_DefaultHandler(int ErrorCode, const char* ErrorMsg, const char* FName, int LNum) { //by thilo: added "const" specifier to avoid warnings

   fprintf (stderr, "%s (%d) [%s:%d]\n", ErrorMsg, ErrorCode, FName, LNum);
   exit(ErrorCode);

}

//=============================================================
// Point to Default Error Handler Function
//
//  int (*Func)(float) : Func is a pointer to function varible.
//  
//  The pointed function takes a float varible and return an int
//
//============================================================
void (*FunPoint)(int ErrorCode, const char* ErrorMsg, const char* FName, int LName) = (FunPoint == NULL) ? SV_DefaultHandler : FunPoint; //by thilo: added "const" specifier to avoid warnings

//==================================================
// User selected error handler
//==================================================
void SV_SetErrorHandler (void (*UserHandler)(int ErrorCode, const char* ErrorMsg, const char* FName, int LName) )  { //by thilo: added "const" specifier to avoid warnings

  FunPoint = UserHandler;

}


