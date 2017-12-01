#!/bin/bash

#easyBMP init
#pre-todo:  
#post-todo: make
#           make install

if [ -e "Makefile" ]; then 
	make clean
fi

if [ ! -e "$HOME/private" ]; then 
	mkdir "$HOME/private"
fi
