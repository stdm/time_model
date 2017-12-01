/**************************************************************************/
/*	Encapsulates methods to read from a Java sream via JNI; most code is  */
/*  taken from SV_DataIO																									*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.12.2009																								*/
/**************************************************************************/

#ifndef __SC_JStreamReader_H__
#define __SC_JStreamReader_H__

#include "SC_Api.h"

#ifdef SC_USE_JNI

#include <jni.h>
#include <SV_DataIO.h>

class SC_JStreamReader {

	private :

		JNIEnv *env; //the java envirnment (pointer to)
		jobject jStreamObject; //the stream object itself (pointer to)
		jclass jStreamClass; //class of the stream object
		jmethodID midRead; //method-id for read()
		jmethodID midReadByteArray; //method-id for readByteArray()
		SV_DataIO io;

	protected:

	public :

	  SC_JStreamReader(JNIEnv *env = NULL, jobject jStreamObject = NULL);
		virtual ~SC_JStreamReader();

		//====================================================================================================================
		// sets the JNI environment and corresponding members
		//====================================================================================================================
		void setEnvironment(JNIEnv *env, jobject jStreamObject);
		JNIEnv* getEnvironment(void) {return this->env;}
		jobject getStream(void) {return this->jStreamObject;}

		//====================================================================================================================
		// returns true if the JNI environment and stream is set properly
		//====================================================================================================================
		bool isInitialized(void) {return (this->env!=NULL && this->jStreamObject!=NULL);}

		//====================================================================================================================
		// read and discard #bytes bytes
		//====================================================================================================================
		int consumeBytes(unsigned int bytes);

		//====================================================================================================================
		// read a non-boolean/char scalar from the current position of the stream, move the read-cursor; returns nr of bytes 
		// actually read or SVLIB_Fail
		//====================================================================================================================
		template<typename T> int readScalar(T &scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) {
			bool isFloat, isInt, success = true;
			unsigned char sizeInStream = io.getDatatypeSize(scalar, streamSizes, isInt, isFloat), sizeInCode = io.getDatatypeSize(scalar, currentSizes);

			if (sizeInStream == 0) {
				REPORT_ERROR(SVLIB_BadData, "Unhandled data type to read from jStream");
				return SVLIB_Fail; //the size of the given type as stored in the jStream couldn't be deduced because it's unknown
			}

			//read as many bytes from the jStream as where stored for this variable
			jbyteArray jTmp = this->env->NewByteArray(sizeInStream);
			jint length = this->env->CallIntMethod(this->jStreamObject, this->midReadByteArray, jTmp, 0, sizeInStream);
			jbyte *value = this->env->GetByteArrayElements(jTmp, NULL);

			//care for endianness-issues
			if (currentSizes.endianness != streamSizes.endianness) {
				io.invertByteOrder((char*)(value), sizeInStream);
			}

			//care for different datatype-sizes 
			if (sizeInCode != sizeInStream) { //sizes differ in stream and the current program
				if (isInt == true) {
					switch (sizeInStream) {
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
							REPORT_ERROR(SVLIB_BadData, "Unhandled data size to read from jStream");
							success = false;
					}
				} else if (isFloat == true) {
					REPORT_ERROR(SVLIB_BadData, "Unhandled float-size differences while reading from jStream");
					success = false;
				} else { //bool or char... should never hit!
					REPORT_ERROR(SVLIB_BadData, "Unhandled data-type/-size (not int, not float) to read from jStream");
					success = false;
				}
			} else { //sizes are identical, so just doing a typecast will do
				scalar = *(T*)(value);
			}

			this->env->ReleaseByteArrayElements(jTmp, value, 0);
			this->env->DeleteLocalRef(jTmp);

			return (length==sizeInStream && success==true) ? sizeInStream : SVLIB_Fail;
		}

		//====================================================================================================================
		// read a boolean scalar from the current position of the stream, move the read-cursor; returns nr of bytes actually read 
		// or SVLIB_Fail
		//====================================================================================================================
		int readScalar(bool &scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) {
			if (currentSizes.boolSize != 1 || streamSizes.boolSize != 1) {
				REPORT_ERROR(SVLIB_BadData, "Can't handle boolean values not 1 byte long");
			}

			jint byte = this->env->CallIntMethod(this->jStreamObject, this->midRead);
			scalar = (byte > 0);

			return (byte != -1) ? 1 : SVLIB_Fail;
		}

		//====================================================================================================================
		// read a character scalar from the current position of the stream, move the read-cursor; returns nr of bytes actually 
		// read or SVLIB_Fail
		//====================================================================================================================
		int readScalar(char &scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) {
			if (currentSizes.charSize != 1 || streamSizes.charSize != 1) {
				REPORT_ERROR(SVLIB_BadData, "Can't handle characters not 1 byte long");
			}

			jint byte = this->env->CallIntMethod(this->jStreamObject, this->midRead);
			scalar = (char)(byte);

			return (byte != -1) ? 1 : SVLIB_Fail;
		}

		//====================================================================================================================
		// just read as much bytes from the stream as the desired type would need
		//====================================================================================================================
		template<typename T> int consumeScalar(T scalar, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) {
			char *value;
			unsigned char sizeInStream = io.getDatatypeSize(scalar, streamSizes);

			if (sizeInStream == 0) {
				REPORT_ERROR(SVLIB_BadData, "Unhandled data type to read from jStream");
				return SVLIB_Fail; //the size of the given type as stored in the stream couldn't be deduced because it's unknown
			}

			//read as many bytes from the stream as where stored for this variable
			jbyteArray jTmp = this->env->NewByteArray(sizeInStream);
			jint length = this->env->CallIntMethod(this->jStreamObject, this->midReadByteArray, jTmp, 0, sizeInStream);
			this->env->DeleteLocalRef(jTmp);

			return (length == sizeInStream) ? sizeInStream : SVLIB_Fail;
		}

		//====================================================================================================================
		// read a dim-array from the current position of the stream, move the read-cursor
		//====================================================================================================================
		template<typename T> int readArray(T* array, unsigned long int dim, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) { 
			int success, res = 0;
			unsigned char sizeInCode = io.getDatatypeSize(array[0], currentSizes), sizeInStream = io.getDatatypeSize(array[0], streamSizes);

			if (sizeInCode == sizeInStream && (currentSizes.endianness == streamSizes.endianness || sizeInCode == 1)) {
				jbyteArray jTmp = this->env->NewByteArray(dim*sizeInStream);
				jint length = this->env->CallIntMethod(this->jStreamObject, this->midReadByteArray, jTmp, 0, dim*sizeInStream);
				this->env->GetByteArrayRegion(jTmp, 0, length, (jbyte*)(array));
				this->env->DeleteLocalRef(jTmp);
				res = length;
			} else {
				for (unsigned long int d = 0; d < dim; d++) {
					success = readScalar(array[d], currentSizes, streamSizes);
					if (success != SVLIB_Fail) { //read as much as wanted from the stream so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
						res += success;
					}
				}
			}

			return res;
		}

		//====================================================================================================================
		// consume as many bytes from the stream as a dim-array would need
		//====================================================================================================================
		template<typename T> int consumeArray(T* array, unsigned long int dim, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) { 
			int success, res = 0;
			T tmp = 0;
			unsigned char sizeInCode = io.getDatatypeSize(array[0], currentSizes), sizeInStream = io.getDatatypeSize(array[0], streamSizes);

			if (sizeInCode == sizeInStream && (currentSizes.endianness == streamSizes.endianness || sizeInCode == 1)) {
				jbyteArray jTmp = this->env->NewByteArray(dim*sizeInStream);
				jint length = this->env->CallIntMethod(this->jStreamObject, this->midReadByteArray, jTmp, 0, dim*sizeInStream);
				this->env->DeleteLocalRef(jTmp);
				res = length;
			} else {
				for (unsigned long int d = 0; d < dim; d++) {
					success = consumeScalar(tmp, currentSizes, streamSizes);
					if (success != SVLIB_Fail) { //read as much as wanted from the stream so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
						res += success;
					}
				}
			}

			return res;
		}

		//====================================================================================================================
		// read a row*col-matrix from the current position of the stream, move the read-cursor
		//====================================================================================================================
		template<typename T> int readMatrix(T** matrix, unsigned long int rows, unsigned long int cols, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) {
			int success, res = 0;

			for (unsigned long int y = 0; y < rows; y++) {
				for (unsigned long int x = 0; x < cols; x++) {
					success = readScalar(matrix[y][x], currentSizes, streamSizes);
					if (success != SVLIB_Fail) { //read as much as wanted from the stream so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
						res += success;
					}
				}
			}

			return res;
		}

		//====================================================================================================================
		// consume as many bytes from the file as a row*col-matrix would need (and discard the values, matrix is just to 
		// provide a type)
		//====================================================================================================================
		template<typename T> int consumeMatrix(T** matrix, unsigned long int rows, unsigned long int cols, SV_DataIO::SV_DatatypeSizes currentSizes, SV_DataIO::SV_DatatypeSizes streamSizes) {
			int success, res = 0;
			T tmp = 0;

			for (unsigned long int y = 0; y < rows; y++) {
				for (unsigned long int x = 0; x < cols; x++) {
					success = consumeScalar(tmp, currentSizes, streamSizes);
					if (success != SVLIB_Fail) { //read as much as wanted from the stream so that the read-cursor is where it is supposed to be after the function-call, regardless of overall success
						res += success;
					}
				}
			}

			return res;
		}
};

#endif
#endif
