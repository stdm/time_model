/*************************************************
*                                                *
*  EasyBMP Cross-Platform Windows Bitmap Library * 
*                                                *
*  Authors: Paul Macklin                         *
*                                                *
*   email: pmacklin@math.uci.edu                 *
*                                                *
*    file: AdditionalSource.cpp                  *
*    date: 9-27-2006                             *
* version: 1.05.00                               *
*                                                *
*   License: BSD (revised)                       *
* Copyright: 2006 by the EasyBMP Project         * 
*                                                *
* description: Test of the EasyBMP Geometry and  *
*              Font extensions.                  *
*                                                *
*************************************************/

#include "AdditionalSource.h"

bool DoSomethingFonty( void )
{
 BMP Temp;
 Temp.SetSize(100,100);
 RGBApixel Black;
 Black.Red = 0;
 Black.Green = 0;
 Black.Blue = 0;
 PrintLetter( Temp,'A',2,2,96,Black );
 return true;
}
