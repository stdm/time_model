CC = g++
CFLAGS = -O3 -s

OBJECTS = FontTest.o AdditionalSource.o EasyBMP.o EasyBMP_Geometry.o EasyBMP_Font.o 

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o FontTest $(OBJECTS)
	
FontTest.o: FontTest.cpp
	cp ../EasyBMP.h .
	cp ../EasyBMP.cpp .
	cp ../EasyBMP_DataStructures.h .
	cp ../EasyBMP_BMP.h .
	cp ../EasyBMP_VariousBMPutilities.h .
	cp ../EasyBMP_Geometry.* .
	cp ../EasyBMP_Font.* .
	$(CC) $(CFLAGS) -c FontTest.cpp
	
EasyBMP.o: EasyBMP.cpp
	$(CC) $(CFLAGS) -c EasyBMP.cpp
	
EasyBMP_Font.o: EasyBMP_Font.cpp
	$(CC) $(CFLAGS) -c EasyBMP_Font.cpp
	
EasyBMP_Geometry.o: EasyBMP_Geometry.cpp
	$(CC) $(CFLAGS) -c EasyBMP_Geometry.cpp
	
AdditionalSource.o: AdditionalSource.cpp
	$(CC) $(CFLAGS) -c AdditionalSource.cpp
	
clean:
	rm EasyBMP.h
	rm EasyBMP.cpp
	rm EasyBMP_DataStructures.h
	rm EasyBMP_BMP.h
	rm EasyBMP_VariousBMPutilities.h
	rm EasyBMP_Geometry.*
	rm EasyBMP_Font.*
	rm *.bmp
	rm *.o
