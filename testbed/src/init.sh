#!/bin/bash

#testbed init
#pre-todo:  "svn up ." on complete "Sprecherklassifikation" directory
#           build svlib
#post-todo: make
#           make install

if [ -e "Makefile" ]; then 
	make distclean
fi

if [ ! -e "configure" ]; then 
	aclocal
	autoheader
	autoconf
	
	touch NEWS
	touch README
	touch AUTHORS
	touch ChangeLog
	
	libtoolize --automake
	
	automake -a
fi

if [ ! -e "$HOME/private" ]; then 
	mkdir "$HOME/private"
fi

export CXXFLAGS="$CXXFLAGS -I$HOME/private/include"
export LDFLAGS="$LDFLAGS -L$HOME/private/lib"

./configure --prefix="$HOME/private"

