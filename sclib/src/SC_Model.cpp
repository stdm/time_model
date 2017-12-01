/**************************************************************************/
/*    Derived from:																												*/
/*      - SV_Model as a base-class for new speaker models									*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.02.2006																								*/
/**************************************************************************/

#include <math.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <iomanip>
#include <string.h>
#include "SC_Model.h"
#include "SC_TweakableParameters.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model::SC_Model(SC_TweakableParameters *pTweak) : SV_Model() {
  this->pTweak = pTweak;
	this->trainingDataCount = 0;
	sprintf(this->lastUsedFileName, "");
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model::SC_Model(const SC_Model& pParent) : SV_Model() {
	this->Hdr = pParent.Hdr;
	this->Next = pParent.Next;
  this->pTweak = pParent.pTweak;
	this->trainingDataCount = pParent.trainingDataCount;
	sprintf(this->lastUsedFileName, "%s", pParent.lastUsedFileName);
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model::~SC_Model() {

}

//====================================================================================================================
// assignment operator
//====================================================================================================================
SC_Model& SC_Model::operator=(const SC_Model& pParent) {
	if (this != &pParent) {
		this->Hdr = pParent.Hdr;
		this->Next = pParent.Next;
		this->pTweak = pParent.pTweak;
		this->trainingDataCount = pParent.trainingDataCount;
		sprintf(this->lastUsedFileName, "%s", pParent.lastUsedFileName);
	}

	return *this;
}

//====================================================================================================================
// opens the file just as the svlib version, but manages to remember the given filename in a new class member
//====================================================================================================================
void SC_Model::OpenFile(const char *FName, int Mode) {

	if (Mode == WRITE_MODEL) {
		this->DFile.open(FName, ios::out|ios::binary);  // truncate
	}	else if (Mode == APPEND_MODEL) {
		this->DFile.open(FName, ios::app|ios::binary);  // append
	}	else if (Mode == READ_MODEL) {
		this->DFile.open(FName, ios::in|ios::binary);   // read
	}	else {
		REPORT_ERROR(SVLIB_FileErr, "OpenFile"); 
	}
	
	if (DFile.fail()) {
		REPORT_ERROR(SVLIB_FileErr, "OpenFile");
	} else {
		sprintf(this->lastUsedFileName, "%s", FName);
	}
	
	return;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& operator<< (ostream& OutS, SC_Model& Data) {
	return Data.modelOut(OutS);
}