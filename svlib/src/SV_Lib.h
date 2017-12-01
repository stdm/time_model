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
//    This is a master include file.
//
//    Author  : Jialong HE
//    Date    : March 14, 1999
//************************************************************************
#ifndef __SVLIB_H__
#define __SVLIB_H__

#include "SV_General.h"
#include "SV_Error.h"
#include "SV_Signal.h"
#include "SV_Signal_NIST.h"
#include "SV_Model.h"
#include "SV_Model_DTW.h"
#include "SV_Model_VQ.h"
#include "SV_Model_Gaus.h"
#include "SV_Model_GMM.h"
#include "SV_Model_CHMM.h"
#include "SV_Model_DHMM.h"
#include "SV_Feature.h"
#include "SV_Feature_LPCC.h"
#include "SV_Feature_MFCC.h"
#include "SV_Feature_Pitch.h"
#include "SV_Data.h"
#include "SV_DataIO.h"
#include "GN_Func.h"
#include "GN_FFT.h"
#include "GN_Filter.h"
#include "GN_LPC.h"
#include "GN_Matrix.h"
#include "GN_Rand.h"
#include "GN_SigWin.h"
#include "GN_WaveGen.h"


#endif
