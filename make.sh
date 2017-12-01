#/bin/bash

#update from svn
svn up .

#make & install svlib
cd svlib/src
make
make install

#make & install sclib
cd ../../sclib/src
make
make install

#make testbed
cd ../../testbed/src
make
