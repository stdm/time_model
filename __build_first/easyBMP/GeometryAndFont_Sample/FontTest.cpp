/*************************************************
*                                                *
*  EasyBMP Cross-Platform Windows Bitmap Library * 
*                                                *
*  Authors: Paul Macklin                         *
*                                                *
*   email: pmacklin@math.uci.edu                 *
*                                                *
*    file: FontTest.cpp                          *
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

#include "EasyBMP.h"
#include "EasyBMP_Geometry.h"
#include "EasyBMP_Font.h"
#include "AdditionalSource.h"

using namespace std;

int main( int argc, char* argv[] )
{
 RGBApixel LightGray;
 LightGray.Red = 192;
 LightGray.Green = 192;
 LightGray.Blue = 192;

 // create a baseline image and clear it to a background color
 
 BMP Window;
 Window.SetSize( 640 , 480 );
 for( int j=0 ; j < Window.TellHeight() ; j++ )
 {
  for( int i=0 ; i < Window.TellWidth() ; i++ )
  {
   *Window(i,j) = LightGray;
  }
 }

 // Create a font color, set a "cursor position", create a text string, 
 // and get ready to print!
 
 int CursorI = 3;
 int CursorJ = 3;
 int FontHeight = 4;
 
 RGBApixel FontColor; 
 FontColor.Red = (ebmpBYTE) ebmpRound( (double) CursorJ * 255.0 / (double) Window.TellHeight() );
 FontColor.Green = 0;
 FontColor.Blue = 255 - FontColor.Red; 

 char TextString1 [1024];
 char TextString2 [1024];
 
 strcpy( TextString1 , " 4: abcdefghijklmnopqrstuvwxyz`1234567890-=[]\\;',./");
 strcpy( TextString2 , "   ABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&* ( ) _+{}|:\"<>?");
 TextString2[1] = COPYRIGHT_SYMBOL;
 
 // start printing the message in different colors and sizes
 
 while( CursorJ+2*FontHeight+5 < Window.TellHeight() )
 {
  FontColor.Red = (ebmpBYTE) ebmpRound( (double) CursorJ * 255.0 / (double) Window.TellHeight() );
  FontColor.Green = 0;
  FontColor.Blue = 255 - FontColor.Red; 
  PrintString( Window, TextString1, CursorI, CursorJ, FontHeight, FontColor );
  CursorJ += (FontHeight+5);

  FontColor.Red = (ebmpBYTE) ebmpRound( (double) CursorJ * 255.0 / (double) Window.TellHeight() );
  FontColor.Green = 0;
  FontColor.Blue = 255 - FontColor.Red; 
  PrintString( Window, TextString2, CursorI, CursorJ, FontHeight, FontColor );
  CursorJ += (FontHeight+5);
  
  FontHeight += 1;
  char FontHeightText [3];
  sprintf( FontHeightText, "%2i", FontHeight );
  TextString1[0] = FontHeightText[0];
  TextString1[1] = FontHeightText[1];
  
  DoSomethingFonty();
 }

 Window.WriteToFile( "FontTest.bmp" );

 return 0;
}
