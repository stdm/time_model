//########################################################################
//  
// A C++ class library for automatic speech recognition and 
// speaker recognition (identification and verification). 
//
// This class library is provided "as is" without any express 
// or implied warranty of any kind with respect to this software. 
// In particular the author shall not be liable for any direct, 
// indirect, special, incidental or consequential damages arising 
// in any way from use of the software.
//
//
// Author   : Jialong He,  Copyright (C), all rights reserved. 
// Date     : May, 1999. 
//
// Contact  : Jialong_He@bigfoot.com, Jialong_He@homemail.com
// Web Page : www.bigfoot.com/~Jialong_He
//########################################################################

//************************************************************************
//    Disk input/output for SV_Data.
//
//
//    Author  : Jialong HE
//    Date    : April 27, 1999
//************************************************************************
#ifndef __SV_DataIO_H__
#define __SV_DataIO_H__

#include <iostream> //by thilo
#include <fstream>
#include <typeinfo>
#include "SV_General.h"
#include "SV_Data.h"
#include "SV_Error.h" //by thilo: for SVLIB_Fail macro

//block by thilo
#ifndef _MSC_VER
	#include "sys/types.h"
	typedef __int8_t SV_Int8;
	typedef __int16_t SV_Int16;
	typedef __int32_t SV_Int32;
	typedef __int64_t SV_Int64;
#else
	typedef __int8 SV_Int8;
	typedef __int16 SV_Int16;
	typedef __int32 SV_Int32;
	typedef __int64 SV_Int64;
#endif
//end by thilo

//--------------------------------
// Open mode for data stream 
//--------------------------------
#define  WRITE_REC	1   // append
#define  READ_REC	2	// read

#define SVLIB_LITTLE_ENDIAN 0 //by thilo
#define SVLIB_BIG_ENDIAN 1 //by thilo

//===================================================================
//  This is the base class of loading speech file into memory.
//  It can be used to derive a class to load signal in specific format 
//  such as NIST format.
//===================================================================
class SV_DataIO {

private :
	fstream		DFile;

protected :

public :

	//------------------------------- 
	// Constructor / Destructor
	//-------------------------------
	SV_DataIO();
	virtual ~SV_DataIO();

	//------------------------------- 
	// Open / Close a data file
	//-------------------------------
	void CloseFile(void);
	void OpenFile(const char *FName, int Mode); //by thilo: added const

	//------------------------------- 
	// Write a data record to current
	// opened file.
	//-------------------------------
	int  PutDataRec(SV_Data &Data);

	//------------------------------- 
	// Get data from opened file.
	//-------------------------------
	SV_Data *GetDataRec(void);					// get next record
	SV_Data *GetDataRec(int RecID);				// get next record with RecID
	SV_Data *GetAllRec(void);					// get all record
	SV_Data *GetAllRec(int RecID);				// get all with ID
	SV_Data *oldGetDataRec(void);					//by thilo: get next record
	SV_Data *oldGetDataRec(int RecID);				//by thilo: get next record with RecID
	SV_Data *oldGetAllRec(void);					//by thilo get all record
	SV_Data *oldGetAllRec(int RecID);				//by thilo: get all with ID
  long int getAllRecRowCount(int &cols);   //by thilo: returns the overall number of rows (and the number of cols, if equal in all recs, otherwise 0) of all the datarecs in the file

	//====================================================================================================================
	//	By thilo: What follows is code that helps making *ALL* binary file i/o (not only of SV_Data objects, but also for
	//            models, signals, ...) read-/writable on systems with different endianness and/or 32bit/64bit
	//            It has to be put here (in the SV_Lib) 'cause some binary i/o already happens here...
	//  TODO: needs speed optimizations (simplest form: just use old methods when no interchange is necessary, as in 
	//        readArray())!!!
	//====================================================================================================================
	typedef struct {
		unsigned char boolSize;
		unsigned char charSize;
		unsigned char shortSize;
		unsigned char intSize;
		unsigned char longSize;
		unsigned char longLongSize;
		unsigned char floatSize;
		unsigned char doubleSize;
		unsigned char longDoubleSize;
		unsigned char endianness;
	} SV_DatatypeSizes; //a struct containing sizes (in bytes) for all relevant datatypes as read/written form/to files

	SV_DataIO::SV_DatatypeSizes& getCurrentDatatypeSizes(SV_DataIO::SV_DatatypeSizes &sizes); //return the datatype-sizes of this machine at runtime
	void invertByteOrder(char *value, unsigned char byteCount); //inverts the byteorder of the given value (e.g. to correct for different endianness)
	unsigned char getEndianness(void); //returns endianness of current machine
	bool swapBytes(char* buffer, int swapLength, unsigned long int numberOfBytes); //regard swapLength bytes in the buffer as a unit and invert their order, for all those units in the buffer; return false if numberOfBytes/swapLength is not an integer
	int getNotUsedHeaderSize(SV_DataIO::SV_DatatypeSizes sizes); //returns the size of the NotUsed-member of the header-struct based on the knowledge of the number and type of members and their sizes as given in "sizes"
	int readMachineHeader(istream *file, SV_DataIO::SV_DatatypeSizes &sizes, bool restorePosition = false); //each binary file gets a header containing information on the datatype-sizes (and endianness) of the machine that wrote it; this is extracted and returned; if restorePosition==true, the original read-cursor is restored before returning
  int writeMachineHeader(ostream *file, SV_DatatypeSizes sizes); //the header is saved at the beginning of the file (which is assumed to be empty)
	int consumeBytes(istream *file, unsigned int bytes); //read and discard #bytes bytes

	//====================================================================================================================
	// returns the size of the type of the given scalar according to the given pSizes-structure; also fills the two bools 
	// telling if the giventype is a floating-point- or integer-type
	//====================================================================================================================
	template<typename T> unsigned char getDatatypeSize(T scalar, SV_DataIO::SV_DatatypeSizes sizes, bool &isInt, bool &isFloat) {
		if (typeid(scalar) == typeid(bool)) {
			isInt = false;
			isFloat = false;
			return sizes.boolSize;
		} else if (typeid(scalar) == typeid(char) || typeid(scalar) == typeid(unsigned char)) {
			isInt = false;
			isFloat = false;
			return sizes.charSize;
		} else if (typeid(scalar) == typeid(short int) || typeid(scalar) == typeid(unsigned short int)) {
			isInt = true;
			isFloat = false;
			return sizes.shortSize;
		} else if (typeid(scalar) == typeid(int) || typeid(scalar) == typeid(unsigned int)) {
			isInt = true;
			isFloat = false;
			return sizes.intSize;
		} else if (typeid(scalar) == typeid(long int) || typeid(scalar) == typeid(unsigned long int)) {
			isInt = true;
			isFloat = false;
			return sizes.longSize;
		} else if (typeid(scalar) == typeid(long long int) || typeid(scalar) == typeid(unsigned long long int)) {
			isInt = true;
			isFloat = false;
			return sizes.longLongSize;
		} else if (typeid(scalar) == typeid(float)) {
			isInt = false;
			isFloat = true;
			return sizes.floatSize;
		} else if (typeid(scalar) == typeid(double)) {
			isInt = false;
			isFloat = true;
			return sizes.doubleSize;
		} else if (typeid(scalar) == typeid(long double)) {
			isInt = false;
			isFloat = true;
			return sizes.longDoubleSize;
		}

		return 0; //none of the previous compares matched, so the given type is unknown... indicate that with a deduced length of 0 bytes
	}

	//====================================================================================================================
	// returns the size of the type of the given scalar according to the given pSizes-structure (same as above without the 
	// isInt/isFloat information)
	//====================================================================================================================
	template<typename T> unsigned char getDatatypeSize(T scalar, SV_DataIO::SV_DatatypeSizes sizes) {
		if (typeid(scalar) == typeid(bool)) {
			return sizes.boolSize;
		} else if (typeid(scalar) == typeid(char) || typeid(scalar) == typeid(unsigned char)) {
			return sizes.charSize;
		} else if (typeid(scalar) == typeid(short int) || typeid(scalar) == typeid(unsigned short int)) {
			return sizes.shortSize;
		} else if (typeid(scalar) == typeid(int) || typeid(scalar) == typeid(unsigned int)) {
			return sizes.intSize;
		} else if (typeid(scalar) == typeid(long int) || typeid(scalar) == typeid(unsigned long int)) {
			return sizes.longSize;
		} else if (typeid(scalar) == typeid(long long int) || typeid(scalar) == typeid(unsigned long long int)) {
			return sizes.longLongSize;
		} else if (typeid(scalar) == typeid(float)) {
			return sizes.floatSize;
		} else if (typeid(scalar) == typeid(double)) {
			return sizes.doubleSize;
		} else if (typeid(scalar) == typeid(long double)) {
			return sizes.longDoubleSize;
		}

		return 0; //none of the previous compares matched, so the given type is unknown... indicate that with a deduced length of 0 bytes
	}

	//====================================================================================================================
	// write a scalar to the current position of the file, move the write-cursor; returns nr. of bytes processed or 
	// SVLIB_Fail
	//====================================================================================================================
	template<typename T> int writeScalar(ostream *file, T scalar) {
		file->write((char*)(&scalar), sizeof(scalar));
		return (file->good() == TRUE) ? sizeof(scalar) : SVLIB_Fail;
	}

	//====================================================================================================================
	// write a scalar with explicitly given output-format to the current position of the file, move the write-cursor; 
	// returns nr. of bytes written or SVLIB_Fail
	//====================================================================================================================
	template<typename T> int writeScalar(ostream *file, T scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		bool isInt, isFloat, success = true;
		unsigned char sizeInFile = getDatatypeSize(scalar, fileSizes, isInt, isFloat), sizeInCode = getDatatypeSize(scalar, currentSizes);;

		//care for different datatype-sizes 
		if (sizeInCode != sizeInFile) { //sizes differ in file and the current program
			if (isInt == true) {
				switch (sizeInFile) {
					case 1: {
						SV_Int8 i8 = (SV_Int8)(scalar);
						file->write((char*)(&i8), sizeInFile);
						break;
					}
					case 2: {
						SV_Int16 i16 = (SV_Int16)(scalar);
						if (fileSizes.endianness != currentSizes.endianness) {
							invertByteOrder((char*)&i16, sizeInFile);
						}
						file->write((char*)(&i16), sizeInFile);
						break;
					}
					case 4: {
						SV_Int32 i32 = (SV_Int32)(scalar);
						if (fileSizes.endianness != currentSizes.endianness) {
							invertByteOrder((char*)&i32, sizeInFile);
						}
						file->write((char*)(&i32), sizeInFile);
						break;
					}
					case 8: {
						SV_Int64 i64 = (SV_Int64)(scalar);
						if (fileSizes.endianness != currentSizes.endianness) {
							invertByteOrder((char*)&i64, sizeInFile);
						}
						file->write((char*)(&i64), sizeInFile);
						break;
					}
					default:
						REPORT_ERROR(SVLIB_BadData, "Unhandled data size to write to file");
						success = false;
				}
			} else if (isFloat == true) {
				REPORT_ERROR(SVLIB_BadData, "Unhandled float-size differences while writing to file");
				success = false;
			} else { //bool or char... should never hit!
				REPORT_ERROR(SVLIB_BadData, "Unhandled data-type/-size (not int, not float) to write to file");
				success = false;
			}
		} else { //sizes are identical, so just write it out
			T tmp = scalar;
			if (fileSizes.endianness != currentSizes.endianness) {
				invertByteOrder((char*)&tmp, sizeInFile);
			}
			file->write((char*)(&tmp), sizeInFile);
		}

		return (file->good() == TRUE && success == true) ? sizeInFile : SVLIB_Fail;
	}

	//====================================================================================================================
	// write a boolean scalar with explicitly given output-format to the current position of the file, move the 
	// write-cursor; returns nr. of bytes written or SVLIB_Fail
	//====================================================================================================================
	template<typename T> int writeScalar(ostream *file, bool scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		if (currentSizes.boolSize != 1 || fileSizes.boolSize != 1) {
			REPORT_ERROR(SVLIB_BadData, "Can't handle boolean values not 1 byte long");
		}

		file->write((char*)(&scalar), 1);

		return (file->good() == TRUE) ? 1 : SVLIB_Fail;
	}

	//====================================================================================================================
	// write a character scalar with explicitly given output-format to the current position of the file, move the 
	// write-cursor; returns nr. of bytes written or SVLIB_Fail
	//====================================================================================================================
	template<typename T> int writeScalar(ostream *file, char scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		if (currentSizes.charSize != 1 || fileSizes.charSize != 1) {
			REPORT_ERROR(SVLIB_BadData, "Can't handle charcters not 1 byte long");
		}
		
		file->write((char*)(&scalar), 1);
		
		return (file->good() == TRUE) ? 1 : SVLIB_Fail;
	}

	//====================================================================================================================
	// read a non-boolean/char scalar from the current position of the file, move the read-cursor; returns nr of bytes 
	// actually read or SVLIB_Fail
	//====================================================================================================================
	template<typename T> int readScalar(istream *file, T &scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		char *value;
		bool isFloat, isInt, success = true;
		unsigned char sizeInFile = getDatatypeSize(scalar, fileSizes, isInt, isFloat), sizeInCode = getDatatypeSize(scalar, currentSizes);

		//file->read((char*)(&scalar), sizeof(T)); //old way of doing it => machine-dependant!!!

		if (sizeInFile == 0) {
			REPORT_ERROR(SVLIB_BadData, "Unhandled data type to read from file");
			return SVLIB_Fail; //the size of the given type as stored in the file couldn't be deduced because it's unknown
		}

		//read as many bytes from the file as where stored for this variable
		value = new char[sizeInFile];
		file->read(value, sizeInFile);

		//care for endianness-issues
		if (currentSizes.endianness != fileSizes.endianness) {
			invertByteOrder(value, sizeInFile);
		}

		//care for different datatype-sizes 
		if (sizeInCode != sizeInFile) { //sizes differ in file and the current program
			if (isInt == true) {
				switch (sizeInFile) {
					case 1:
						scalar = (T)(*(SV_Int8*)value);
						break;
					case 2:
						scalar = (T)(*(SV_Int16*)value);
						break;
					case 4:
						scalar = (T)(*(SV_Int32*)value);
						break;
					case 8:
						scalar = (T)(*(SV_Int64*)value);
						break;
					default:
						REPORT_ERROR(SVLIB_BadData, "Unhandled data size to read from file");
						success = false;
				}
			} else if (isFloat == true) {
				REPORT_ERROR(SVLIB_BadData, "Unhandled float-size differences while reading from file");
				success = false;
			} else { //bool or char... should never hit!
				REPORT_ERROR(SVLIB_BadData, "Unhandled data-type/-size (not int, not float) to read from file");
				success = false;
			}
		} else { //sizes are identical, so just doing a typecast will do
			scalar = *(T*)(value);
		}

		delete [] value;

		return (file->good() == TRUE && success == true) ? sizeInFile : SVLIB_Fail;
	}

	//====================================================================================================================
	// read a boolean scalar from the current position of the file, move the read-cursor; returns nr of bytes actually read 
	// or SVLIB_Fail
	//====================================================================================================================
	int readScalar(istream *file, bool &scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		if (currentSizes.boolSize != 1 || fileSizes.boolSize != 1) {
			REPORT_ERROR(SVLIB_BadData, "Can't handle boolean values not 1 byte long");
		}

		file->read((char*)(&scalar), 1);

		return (file->good() == TRUE) ? 1 : SVLIB_Fail;
	}

	//====================================================================================================================
	// read a character scalar from the current position of the file, move the read-cursor; returns nr of bytes actually 
	// read or SVLIB_Fail
	//====================================================================================================================
	int readScalar(istream *file, char &scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		if (currentSizes.charSize != 1 || fileSizes.charSize != 1) {
			REPORT_ERROR(SVLIB_BadData, "Can't handle characters not 1 byte long");
		}

		file->read((char*)(&scalar), 1);

		return (file->good() == TRUE) ? 1 : SVLIB_Fail;
	}

	//====================================================================================================================
	// just read as much bytes from the file as the desired type would need
	//====================================================================================================================
	template<typename T> int consumeScalar(istream *file, T scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		char *value;
		unsigned char sizeInFile = getDatatypeSize(scalar, fileSizes);

		if (sizeInFile == 0) {
			REPORT_ERROR(SVLIB_BadData, "Unhandled data type to read from file");
			return SVLIB_Fail; //the size of the given type as stored in the file couldn't be deduced because it's unknown
		}

		//read as many bytes from the file as where stored for this variable
		value = new char[sizeInFile];
		file->read(value, sizeInFile);
		
		delete [] value;

		return (file->good() == TRUE) ? sizeInFile : SVLIB_Fail;
	}

	//====================================================================================================================
	// write a dim-array to the current position of the file, move the write-cursor
	//====================================================================================================================
	template<typename T> int writeArray(ostream *file, T* array, unsigned long int dim) {
		int success, res = 0;
		bool ok = true;
		
		//file->write((char*)(array), dim * sizeof(T));

		for (unsigned long int d = 0; d < dim; d++) {
			success = writeScalar(file, array[d]);
			if (success != SVLIB_Fail) { //write as much as wanted from the file so that the write-cursor is where it is supposed to be after the function-call, regardless of overall success
				res += success;
			} else {
				ok = false;
			}
		}

		return (ok==true) ? res : SVLIB_Fail;
	}

	//====================================================================================================================
	// write a dim-array to the current position of the file, move the write-cursor, using explicit type sizes
	//====================================================================================================================
	template<typename T> int writeArray(ostream *file, T* array, unsigned long int dim, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		int success, res = 0;
		bool ok = true;
		unsigned char codeSize = getDatatypeSize(array[0], currentSizes);

		if (codeSize == getDatatypeSize(array[0], fileSizes) && (currentSizes.endianness == fileSizes.endianness || codeSize == 1)) {
			file->write((char*)(array), dim * sizeof(T)); //old way of doing it => machine-dependant, but ok here because of the above checks and MUCH faster!!!
			res = dim * sizeof(T);
		} else {
			for (unsigned long int d = 0; d < dim; d++) {
				success = writeScalar(file, array[d], currentSizes, fileSizes);
				if (success != SVLIB_Fail) { //write as much as wanted from the file so that the write-cursor is where it is supposed to be after the function-call, regardless of overall success
					res += success;
				} else {
					ok = false;
				}
			}
		}

		return (ok==true) ? res : SVLIB_Fail;
	}

	//====================================================================================================================
	// read a dim-array from the current position of the file, move the read-cursor
	//====================================================================================================================
	template<typename T> int readArray(istream *file, T* array, unsigned long int dim, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) { 
		int success, res = 0;
		unsigned char codeSize = getDatatypeSize(array[0], currentSizes);
		
		if (codeSize == getDatatypeSize(array[0], fileSizes) && (currentSizes.endianness == fileSizes.endianness || codeSize == 1)) {
			file->read((char*)(array), dim * sizeof(T)); //old way of doing it => machine-dependant, but ok here because of the above checks and MUCH faster!!!
			res = dim * sizeof(T);
		} else {
			for (unsigned long int d = 0; d < dim; d++) {
				success = readScalar(file, array[d], currentSizes, fileSizes);
				if (success != SVLIB_Fail) { //read as much as wanted from the file so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
					res += success;
				}
			}
		}

		return res;
	}

	//====================================================================================================================
	// consume as many bytes from the file as a dim-array would need
	//====================================================================================================================
	template<typename T> int consumeArray(istream *file, T* array, unsigned long int dim, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) { 
		int success, res = 0;
		T tmp = 0;

		//speed up consuming compared to the version below: once and for all, get the number of bytes of a single element to consume and then over-seek as much bytes as there are elements in the array
		success = consumeScalar(file, tmp, currentSizes, fileSizes);
		if (success != SVLIB_Fail) {
			file->seekg((dim-1)*success, ios::cur);
			res = dim * success;
		} else {
			res = SVLIB_Fail;
		}

		/*
		for (unsigned long int d = 0; d < dim; d++) {
			success = consumeScalar(file, tmp, currentSizes, fileSizes);
			if (success != SVLIB_Fail) { //read as much as wanted from the file so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
				res += success;
			}
		}
		*/

		return res;
	}

	//====================================================================================================================
	// write a row*col-matrix to the current position of the file, move the write-cursor
	//====================================================================================================================
	template<typename T> int writeMatrix(ostream *file, T** matrix, unsigned long int rows, unsigned long int cols) {
		int success, res = 0;
		
		//file->write((char*)(matrix[0]), cols * rows * sizeof(T));

		for (unsigned long int y = 0; y < rows; y++) {
			for (unsigned long int x = 0; x < cols; x++) {
				success = writeScalar(file, matrix[y][x]);
				if (success != SVLIB_Fail) { //write as much as wanted from the file so that the write-cursor is where it is supposed to be after the function-call, regardless of overall success
					res += success;
				}
			}
		}

		return res;
	}

	//====================================================================================================================
	// write a row*col-matrix to the current position of the file, move the write-cursor, using explicit sizes
	//====================================================================================================================
	template<typename T> int writeMatrix(ostream *file, T** matrix, unsigned long int rows, unsigned long int cols, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		int success, res = 0;
		
		//file->write((char*)(matrix[0]), cols * rows * sizeof(T));

		for (unsigned long int y = 0; y < rows; y++) {
			for (unsigned long int x = 0; x < cols; x++) {
				success = writeScalar(file, matrix[y][x], currentSizes, fileSizes);
				if (success != SVLIB_Fail) { //write as much as wanted from the file so that the write-cursor is where it is supposed to be after the function-call, regardless of overall success
					res += success;
				}
			}
		}

		return res;
	}

	//====================================================================================================================
	// read a row*col-matrix from the current position of the file, move the read-cursor
	//====================================================================================================================
	template<typename T> int readMatrix(istream *file, T** matrix, unsigned long int rows, unsigned long int cols, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		int success, res = 0;
		unsigned char codeSize = getDatatypeSize(matrix[0][0], currentSizes);

		//file->read((char*)(matrix[0]), cols * rows * sizeof(T)); //old way of doing it => machine-dependant!!!

		if (codeSize == getDatatypeSize(matrix[0][0], fileSizes) && (currentSizes.endianness == fileSizes.endianness || codeSize == 1)) {
			file->read((char*)(matrix[0]), cols * rows * sizeof(T)); //old way of doing it => machine-dependant, but ok here because of the above checks and MUCH faster!!!
			res = cols * rows * sizeof(T);
		} else {
			for (unsigned long int y = 0; y < rows; y++) {
				for (unsigned long int x = 0; x < cols; x++) {
					success = readScalar(file, matrix[y][x], currentSizes, fileSizes);
					if (success != SVLIB_Fail) { //read as much as wanted from the file so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
						res += success;
					}
				}
			}
		}

		return res;
	}

	//====================================================================================================================
	// consume as many bytes from the file as a row*col-matrix would need (and discard the values, matrix is just to 
	// provide a type)
	//====================================================================================================================
	template<typename T> int consumeMatrix(istream *file, T** matrix, unsigned long int rows, unsigned long int cols, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes fileSizes) {
		int success, res = 0;
		T tmp = 0;

		//speed up consuming compared to the version below: once and for all, get the number of bytes of a single element to consume and then over-seek as much bytes as there are elements in the matrix
		success = consumeScalar(file, tmp, currentSizes, fileSizes);
		if (success != SVLIB_Fail) {
			file->seekg(rows*cols*success - success, ios::cur);
			res = rows * cols * success;
		} else {
			res = SVLIB_Fail;
		}

		/*
		for (unsigned long int y = 0; y < rows; y++) {
			for (unsigned long int x = 0; x < cols; x++) {
				success = consumeScalar(file, tmp, currentSizes, fileSizes);
				if (success != SVLIB_Fail) { //read as much as wanted from the file so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
					res += success;
				}
			}
		}
		*/

		return res;
	}
	//====================================================================================================================
	//	End by thilo
	//====================================================================================================================

};   // class SV_DataIO


#endif

