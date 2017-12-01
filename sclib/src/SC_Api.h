/**************************************************************************/
/*    Some definitions for library-export/-import as well as some         */
/*    definitions that control the usage of external libraries            */
/*                                                                        */
/*    Comment out the defines if a library is unavailable - this will     */
/*    surely destroy functionality but the sclib will compile anyway;     */
/*    only do this if you know that the unavailable parts will not be     */
/*    used by the caller of this truncated sclib! For the compilation to  */
/*    work correctly, the additional libraries have to be removed from    */
/*    the linker call as well.                                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 26.09.2006																								*/
/**************************************************************************/

#ifndef __SC_Api_H__
#define __SC_Api_H__

//====================================================================================================================
// definition of the correct .so/.dll-import/-export macros
//====================================================================================================================
#ifdef _WIN32
  #ifdef SCLIB_EXPORTS
    #define SCLIB_API __declspec(dllexport)
    #define SCLIB_API_C extern "C" SCLIB_API
  #else
    #define SCLIB_API __declspec(dllimport)
    #define SCLIB_API_C extern "C" SCLIB_API
  #endif
#else
  #define SCLIB_API
  #define SCLIB_API_C extern "C"
#endif

//====================================================================================================================
// the ffmpeg libraries (libavcodec.so/avcodec-51.dll, libavformat.so/avformat-51-dll, libavutil.so/avutil-49.dll)
// are used to open all file formats different from Microsoft RIFF Wave files and NIST SPHERE files for loading the 
// signal
//====================================================================================================================
//#define SC_USE_FFMPEG

//====================================================================================================================
// libsamplerate.so/.dll is used for resampling loaded signals to a specific samplerate needed by feature-extractors 
// or algorithms as the LZL audiotype classificator
//====================================================================================================================
//#define SC_USE_LIBSAMPLERATE

//====================================================================================================================
// libeasybmp.so/easybmp.dll is used for writing bitmap files, e.g. for similarity matrixes, PDE plots and SDP 
// features visualization
//====================================================================================================================
//#define SC_USE_EASYBMP

//====================================================================================================================
// libfann.so/.dll is used for modelling/classification via neural networks
//====================================================================================================================
//#define SC_USE_FANN

//====================================================================================================================
// below is third party code that is directly incorporated into the sclib code but can be switched on/off just as with 
// the external libraries before
//====================================================================================================================

//--------------------------------------------------------------------------------------------------------------------
// Yossi Rubner's reference implementation of the Earth Mover's Distance
//--------------------------------------------------------------------------------------------------------------------
#define SC_USE_EARTHMOVERSDISTANCE

//--------------------------------------------------------------------------------------------------------------------
// For convenience of compilation, the single source file of the libsvm was directly copied herein; it implements a
// Support Vector Machine library
//--------------------------------------------------------------------------------------------------------------------
#define SC_USE_LIBSVM

//--------------------------------------------------------------------------------------------------------------------
// The ESPS/Talkin get_f0() pitch tracker as implemented in SC_Feature_Pitch.cpp
//--------------------------------------------------------------------------------------------------------------------
#define SC_USE_ESPSPITCH

//--------------------------------------------------------------------------------------------------------------------
// The ESPS/Talkin formant tracker as implemented in SC_Feature_Formant.cpp
//--------------------------------------------------------------------------------------------------------------------
#define SC_USE_ESPSFORMANT

//--------------------------------------------------------------------------------------------------------------------
// The Speex implementation of LSP feature extraction
//--------------------------------------------------------------------------------------------------------------------
#define SC_USE_SPEEXLSP

//--------------------------------------------------------------------------------------------------------------------
// The ePhone/MELP implementation of LSP feature extraction
//--------------------------------------------------------------------------------------------------------------------
#define SC_USE_MELPLSP

//--------------------------------------------------------------------------------------------------------------------
// Uses the Java Native Interface for direct coupling of this lib with Java
//--------------------------------------------------------------------------------------------------------------------
//#define SC_USE_JNI

//--------------------------------------------------------------------------------------------------------------------
// Several routines written in Matlab and deployed via the Matlab compiler and deploytool (Windows only at the 
// moment): 
//  - The MFCC feature extraction routine form Malcolm Slaneys (Auditory Toolbox)
//  - The MFCCs by Dan Ellis (Rastamat)
//  - The l1-minimization algorithms from Justin Romberg (l1-Magic)
//--------------------------------------------------------------------------------------------------------------------
//#define SC_USE_MATLAB

#endif
