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
//    Put common function here. This file is included
//    by other header files. 
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#ifndef __SV_General_H__
#define __SV_General_H__

#include <stdlib.h>

//===================================================================
//  Some constant definitions
//===================================================================
#define LARGE_D	1.0e+30
#define LARGE	float(1.0e+30)
#define TINY	2.22e-16
#define	TRUE		1
#define FALSE		0
#define REAL   double

//-----------------------------------
// COMPLEX structure
//-----------------------------------
typedef struct {
    double real, imag;
} COMPLEX;


//===================================================================
//  Bit operations
//===================================================================
#define BitSet(arg,posn) ((arg) | (1L << (posn)))
#define BitClr(arg,posn) ((arg) & ~(1L << (posn)))
#define BitFlp(arg,posn) ((arg) ^ (1L << (posn)))



int getopt(int argc, char **argv, char *optstring);


#endif

