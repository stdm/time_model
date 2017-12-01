//************************************************************************
//    Implementation of SV_Data class
//
//
//    Author  : Jialong HE
//    Date    : April 27, 1999
//************************************************************************
#include <iostream>
#include <string.h>
#include "SV_Data.h"
#include "SV_Error.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Data::SV_Data() {

/*
	ClassSig = 196273;
	//-------------------------------
	// Data default header contents
	//-------------------------------
	Hdr.Signature[0]	= 'J';
	Hdr.Signature[1]	= 'H';
	Hdr.Signature[2]	= 62;
	Hdr.Signature[3]	= 7;
	Hdr.Signature[4]	= 3;
	Hdr.Signature[5]	= 64;
	Hdr.Signature[6]	= 1;
	Hdr.Signature[7]	= 21;

	Hdr.ByteOrder		= 0x01020304;
	Hdr.Version			= float(0.1);	
	Hdr.ID				= -1;

	strcpy(Hdr.Name, "NoName");   // maximum 16 chars
	
	Row = 0;
	Col = 0;
	Mat  = NULL;
	Next = NULL;
*/ //by thilo: instead:
	init(this);
}

//by thilo, initializes the class as the default-constructor did before
void SV_Data::init(SV_Data* pData) {
	pData->ClassSig = 196273;
	//-------------------------------
	// Data default header contents
	//-------------------------------
	pData->Hdr.Signature[0]	= 0; //by thilo; was 'J';
	pData->Hdr.Signature[1]	= 0; //by thilo; was 'H';
	pData->Hdr.Signature[2]	= 0; //by thilo; was 62;
	pData->Hdr.Signature[3]	= 0; //by thilo; was 7;
	pData->Hdr.Signature[4]	= 0; //by thilo; was 3;
	pData->Hdr.Signature[5]	= 0; //by thilo; was 64;
	pData->Hdr.Signature[6]	= 0; //by thilo; was 1;
	pData->Hdr.Signature[7]	= 0; //by thilo; was 21;

	pData->Hdr.ByteOrder		= 0x01020304;
	pData->Hdr.Version			= float(0.1);	
	pData->Hdr.ID						= -1;
  
  pData->Hdr.frameSize = 0; //by thilo
  pData->Hdr.frameStep = 0; //by thilo
  pData->Hdr.sampleRate = 0; //by thilo

	strcpy(pData->Hdr.Name, "NoName");   // maximum 16 chars
	
	pData->Row = 0;
	pData->Col = 0;
	pData->Mat  = NULL;
	pData->Next = NULL;

  pData->justLinked = false;

  return;
}

//==========================================
// Constructor that allocate memory
//==========================================
SV_Data::SV_Data(int M, int N) {

	//SV_Data(); //by thilo: instead:
	init(this);

	Row = M;
	Col = N;
  Alloc();

}

/*
//==========================================
// Copy Constructor
//==========================================
SV_Data::SV_Data(const SV_Data &Src) {

	int I, J;


	ClassSig = 196273;
	Hdr = Src.Hdr;

	Row = Src.Row;
	Col = Src.Col;
  Next= Src.Next;	
	Mat = NULL;    // not constructor called, must set to NULL 

  Alloc();

	//-------------------------------
	// copy data matrix
	//-------------------------------
	for (I=0; I<Row; I++) {
		for (J=0; J<Col; J++) {
			Mat[I][J] = Src.Mat[I][J];
		}
	}

}
*/

//==========================================
//by thilo: copy constructor, that is able 
//to just link to the data-matrix of the 
//parent-object instead of copying it
//==========================================
SV_Data::SV_Data(const SV_Data& Src, bool linkOnly) {
	int I, J;

  ClassSig = 196273;
	Hdr = Src.Hdr;

	Row = Src.Row;
	Col = Src.Col;
  Next= Src.Next;	

  if (linkOnly == false) {
    Mat = NULL;    // not constructor called, must set to NULL 
    Alloc();

	  //-------------------------------
	  // copy data matrix
	  //-------------------------------
	  for (I=0; I<Row; I++) {
		  for (J=0; J<Col; J++) {
			  Mat[I][J] = Src.Mat[I][J];
		  }
	  }
    this->justLinked = false;
  } else {
    this->Mat = Src.Mat;
    this->justLinked = true;
  }
}

//==========================================
//by tilo: copy constructor to copy a 
//submatrix; min is >=, max is <
//==========================================
SV_Data::SV_Data(const SV_Data& Src, int minLen, int maxLen, int minDim, int maxDim) {
	int I, J;

  ClassSig = 196273;
	Hdr = Src.Hdr;

	Row = maxLen - minLen;
	Col = maxDim - minDim;
  Next= Src.Next;	

  Mat = NULL;    // not constructor called, must set to NULL 
  Alloc();

	for (I = 0; I < this->Row; I++) {
		for (J = 0; J < this->Col; J++) {
			Mat[I][J] = Src.Mat[minLen+I][minDim+J];
		}
	}
  this->justLinked = false;
}

//==========================================
// destructor
//==========================================
SV_Data::~SV_Data() {
	if (Mat != NULL && this->justLinked == false) {
		MFree_2D(Mat);
	} else { //block by thilo
		this->Mat = NULL;
	}

	ClassSig = 0;
}

//==========================================
// Allocate memory for Mat
//==========================================
void SV_Data::Alloc(void) {

	if (Mat != NULL && this->justLinked == false) {
		MFree_2D(Mat);
	}

  this->justLinked = false;
	MArray_2D(Mat, Row, Col, float, "Data_Mat");
}

//==========================================
// Test Class's Signature
// Return 1: valid, 0: invalid (deleted)
//==========================================
int SV_Data::Valid(void) {

	if (ClassSig == 196273) {return(1);}
	else {return(0);}
	
}

//=============================================
// Dump a data record in ASCII ot cout
//=============================================
ostream& operator<< (ostream& OutS, SV_Data& Data) {
	int Row, Col;

	if (Data.Valid()) {

		OutS << Data.Row <<" "<<Data.Hdr.ID <<" "<< Data.Col << endl;
		for (Row=0; Row<Data.Row; Row++) {
			for (Col=0; Col<Data.Col; Col++)
				OutS<< Data.Mat[Row][Col]<< " ";
			OutS<<endl;
		}
	}
	return(OutS);
}

//==========================================
// overload assignment operator
//==========================================
SV_Data &SV_Data::operator= (const SV_Data &Src) {
	int I, J;

	if (this->justLinked == false) { //by thilo
		MFree_2D(Mat); //by thilo: maybe there already is a matrix which we need to delete first!
	} //by thilo

	Hdr = Src.Hdr;
	Row = Src.Row;
	Col = Src.Col;
  Alloc();
  Next = Src.Next;	

  //-------------------------------
	// copy data matrix
	//-------------------------------
	for (I=0; I<Row; I++) {
		for (J=0; J<Col; J++) {
			Mat[I][J] = Src.Mat[I][J];
		}
	}

	return (*this);
}

//====================================================================================================================
// by thilo:
//
// Merge linked list data into one great feature-container;
// Take not the whole linked list till it's end, but only maxSegments segments
// Ignore maxSegments, if maxSegments == 0
//
// if 0>=dim>Cols, only the specified column is merged, resulting in a one-column feature set
//
// this was formerly placed in the SV_Model/SC_Model classes, but makes more sense here...
//====================================================================================================================
SV_Data* SV_Data::MergeData(unsigned long int maxSegments, int dim) {
	SV_Data *pData = this, *pTData, *DataCurr = NULL;
	int TotalVec = 0;
	int Row, Col, Dim;
	unsigned long int count = 0;

  if (pData == NULL) {REPORT_ERROR(SVLIB_BadArg, "\nNo Data in list to merge\n");}
	
  //break here if nothing is to do
  //if (maxSegments == 1) {return pData;}
  //DON'T DO THIS because it is assummed throughout the programm that the result is a new dataset (which must be deleted)!

	// Count total vectors
	DataCurr = pData;
	Dim = pData->Col;
	while ((DataCurr != NULL) && ((count < maxSegments) || (maxSegments == 0))) {
    if (Dim != DataCurr->Col) {REPORT_ERROR(SVLIB_BadArg, "\nFeature-dimension differs in linked list to merge\n");}
		count++;
		TotalVec += DataCurr->Row; 
		DataCurr = DataCurr->Next;
	}  // while

	// Allocate memory for New SV_Data
	pTData = new SV_Data;
	pTData->Row = TotalVec;
	pTData->Col = (dim>=0 && dim<Dim) ? 1 : Dim;
	pTData->Alloc();
	pTData->Hdr = this->Hdr; //by thilo

	//block commented by thilo, see above
  //pTData->Hdr.frameSize = pData->Hdr.frameSize;
  //pTData->Hdr.frameStep = pData->Hdr.frameStep;
  //pTData->Hdr.ID = pData->Hdr.ID;

	// Accumulate patterns
	DataCurr = pData;
	TotalVec = 0;
	count = 0;
	while ((DataCurr != NULL) && ((count < maxSegments) || (maxSegments == 0))) {
		for (Row=0; Row<DataCurr->Row; Row++) { //accumulate vectors in this pattern
			if (dim>=0 && dim<Dim) {
				pTData->Mat[TotalVec+Row][0] = DataCurr->Mat[Row][dim];
			} else {
				for (Col=0; Col<Dim; Col++) {
					pTData->Mat[TotalVec+Row][Col] = DataCurr->Mat[Row][Col];
				}
			}
		}
		TotalVec += DataCurr->Row; 
		DataCurr = DataCurr->Next;
		count++;
	}  // while

	return(pTData);
}

