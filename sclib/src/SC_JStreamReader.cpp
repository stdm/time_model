/**************************************************************************/
/*	Encapsulates methods to read from a Java sream via JNI; most code is  */
/*  taken from SV_DataIO																									*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 15.12.2009																								*/
/**************************************************************************/

#include "SC_JStreamReader.h"
#include <SV_Error.h>

#ifdef SC_USE_JNI
//====================================================================================================================
//	constructor; ptweak is only needed for storing the results in the specified debug-dir
//====================================================================================================================
SC_JStreamReader::SC_JStreamReader(JNIEnv *env, jobject jStreamObject) {
	setEnvironment(env, jStreamObject);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_JStreamReader::~SC_JStreamReader() {

}

//====================================================================================================================
// sets the JNI environment and corresponding members
//====================================================================================================================
void SC_JStreamReader::setEnvironment(JNIEnv *env, jobject jStreamObject) {
	this->env = env;
	this->jStreamObject = jStreamObject;

	if (this->env!=NULL && this->jStreamObject!=NULL) {
		this->jStreamClass = this->env->GetObjectClass(this->jStreamObject); 
		this->midRead = this->env->GetMethodID(this->jStreamClass, "read", "()I");
		this->midReadByteArray = this->env->GetMethodID(this->jStreamClass, "read", "([BII)I");
	} else {
		this->jStreamClass = NULL;
		this->midRead = NULL;
		this->midReadByteArray = NULL;
	}

	return;
}

//====================================================================================================================
//  Read and discard 'bytes' bytes
//====================================================================================================================
int SC_JStreamReader::consumeBytes(unsigned int bytes) {
	jbyteArray value = this->env->NewByteArray(bytes);
	jint length = this->env->CallIntMethod(this->jStreamObject, this->midReadByteArray, value, 0, bytes);
	this->env->DeleteLocalRef(value);

	return (length == bytes) ? (int)(bytes) : SVLIB_Fail;
}
#endif
