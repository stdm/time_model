#!/bin/bash

#libsamplerate 0.1.2 init
#pre-todo:  
#post-todo: make
#           make install

if [ -e "Makefile" ]; then 
	make distclean
fi

if [ ! -e "$HOME/private" ]; then 
	mkdir "$HOME/private"
fi

./configure --prefix=$HOME/private
