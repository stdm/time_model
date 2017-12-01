/*************************************************
*                                                *
*  EasyBMP Cross-Platform Windows Bitmap Library * 
*                                                *
*  Author: Paul Macklin                          *
*   email: macklin01@users.sourceforge.net       *
* support: http://easybmp.sourceforge.net        *
*                                                *
*          file: EasyBMP_VariousBMPutilities.h   *
*    date added: 05-02-2005                      *
* date modified: 12-01-2006                      *
*       version: 1.06                            *
*                                                *
*   License: BSD (revised/modified)              *
* Copyright: 2005-6 by the EasyBMP Project       * 
*                                                *
* description: Various utilities.                *
*                                                *
*************************************************/

#ifndef _EasyBMP_VariousBMPutilities_h_
#define _EasyBMP_VariousBMPutilities_h_

#include "EasyBMP_Api.h" //by thilo

EASYBMP_API BMFH GetBMFH( const char* szFileNameIn );
EASYBMP_API BMIH GetBMIH( const char* szFileNameIn );
EASYBMP_API void DisplayBitmapInfo( const char* szFileNameIn );
EASYBMP_API int GetBitmapColorDepth( const char* szFileNameIn );
EASYBMP_API void PixelToPixelCopy( BMP& From, int FromX, int FromY,  
                       BMP& To, int ToX, int ToY);
EASYBMP_API void PixelToPixelCopyTransparent( BMP& From, int FromX, int FromY,  
                                  BMP& To, int ToX, int ToY,
                                  RGBApixel& Transparent );
EASYBMP_API void RangedPixelToPixelCopy( BMP& From, int FromL , int FromR, int FromB, int FromT, 
                             BMP& To, int ToX, int ToY );
EASYBMP_API void RangedPixelToPixelCopyTransparent( 
     BMP& From, int FromL , int FromR, int FromB, int FromT, 
     BMP& To, int ToX, int ToY ,
     RGBApixel& Transparent );
EASYBMP_API bool CreateGrayscaleColorTable( BMP& InputImage );

EASYBMP_API bool Rescale( BMP& InputImage , char mode, int NewDimension );

#endif
