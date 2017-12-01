CC = gcc
CFLAGS = -O3 -pipe -fomit-frame-pointer -funroll-all-loops -s
PREFIX = $(HOME)/private

all: libeasybmp.a libeasybmp.so

static: libeasybmp.a

shared: libeasybmp.so

libeasybmp.so: EasyBMP.cpp EasyBMP_Font.cpp EasyBMP_Geometry.cpp EasyBMP.h EasyBMP_BMP.h EasyBMP_DataStructures.h EasyBMP_VariousBMPutilities.h EasyBMP_Api.h EasyBMP_Font.h EasyBMP_Geometry.h
	$(CC) $(CFLAS) -fPIC -c EasyBMP_Font.cpp -o EasyBMP_Font.o
	$(CC) $(CFLAS) -fPIC -c EasyBMP_Geometry.cpp -o EasyBMP_Geometry.o
	$(CC) $(CFLAS) -fPIC -c EasyBMP.cpp -o EasyBMP.o
	$(CC) $(CFLAS) -shared -o libeasybmp.so EasyBMP.o EasyBMP_Font.o EasyBMP_Geometry.o
	
libeasybmp.a: EasyBMP.cpp EasyBMP_Font.cpp EasyBMP_Geometry.cpp EasyBMP.h EasyBMP_BMP.h EasyBMP_DataStructures.h EasyBMP_VariousBMPutilities.h EasyBMP_Api.h EasyBMP_Font.h EasyBMP_Geometry.h
	$(CC) $(CFLAS) -fPIC -c EasyBMP_Font.cpp -o EasyBMP_Font.o
	$(CC) $(CFLAS) -fPIC -c EasyBMP_Geometry.cpp -o EasyBMP_Geometry.o
	$(CC) $(CFLAGS) -c EasyBMP.cpp -o EasyBMP.o
	ar rc libeasybmp.a EasyBMP.o EasyBMP_Font.o EasyBMP_Geometry.o
	ranlib libeasybmp.a

clean:
	rm -f *.o
	rm -f *.a
	rm -f *.so
	
install: libeasybmp.so EasyBMP.h EasyBMP_BMP.h EasyBMP_DataStructures.h EasyBMP_VariousBMPutilities.h
	cp libeasybmp.so $(PREFIX)/lib
	cp libeasybmp.a $(PREFIX)/lib
	cp EasyBMP.h $(PREFIX)/include
	cp EasyBMP_BMP.h $(PREFIX)/include
	cp EasyBMP_DataStructures.h $(PREFIX)/include
	cp EasyBMP_VariousBMPutilities.h $(PREFIX)/include
	cp EasyBMP_Font.h $(PREFIX)/include
	cp EasyBMP_Geometry.h $(PREFIX)/include
	cp EasyBMP_Api.h $(PREFIX)/include

