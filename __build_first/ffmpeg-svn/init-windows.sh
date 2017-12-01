#!/bin/bash

#ffmpeg on windows init (SVN rev. 9375 proven to work)
#pre-todo:  read http://arrozcru.no-ip.org/ffmpeg_wiki/tiki-index.php?page=Shared
#           install mingw and msys
#           load ms-visual-c++ environment variables into msys by using the following first line in msys-bat (change path' respectively): call "C:\Program Files\Microsoft Visual Studio .NET 2003\Vc7\bin\vcvars32.bat"
#           start msys-shell (msys.bat)
#           go to the ffmpeg directory
#post-todo: make
#           make install
#           move the following file from PREFIX/bin to PREFIX: avcodec-51.dll avcodec-51.lib avcodec-51.dll avcodec-51.lib avutil-49.dll avutil-49.lib
#           delete PREFIX/pkgconfig, PREFIX/lib, PREFIX/bin
#           get windows-versions of inttypes.h and stdint.h and copy them to PREFIX

if [ -e "Makefile" ]; then 
  make distclean
fi

./configure --prefix=/d/Data/Code/Sprecherklassifikation/dlls --libdir=/d/Data/Code/Sprecherklassifikation/dlls --incdir=/d/Data/Code/Sprecherklassifikation/dlls/ffmpeg --enable-memalign-hack --enable-shared --disable-static --disable-debug --disable-ffserver --disable-ffplay
