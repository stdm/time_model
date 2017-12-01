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
//    a class to hold a data pattern (matrix). 
//
//
//    Author  : Jialong HE
//    Date    : April 27, 1999
//************************************************************************
#ifndef __SV_Data_H__
#define __SV_Data_H__
#include <iostream>

#include "SV_General.h"

using namespace std;

#define DHLen 64         // data pattern header length
#define NotUsedDHLen (DHLen - 24*sizeof(char) - 2*sizeof(long) - sizeof(float) - 3*sizeof(unsigned int)) //by thilo to achieve a 64byte header on almost all machines...; attention: this knowledge is also used inside SV_DataIO, so also check there if something is changed here!
//--------------------------------------------------------- 
// Each record has 64 bytes RecHdr, followed by a matrix
//--------------------------------------------------------- 
typedef struct {
	char  Signature[8];		  // identifying a data record
	                        // by thilo: convention for usage
	                        //          [0] flag specifying for which algorithm (i.e. with which parameter-set) the data was extracted, 0 for non-specific/global data
	                        //          [1] flag specifying used extraction method, if the user can select between different ones
													//          [2] flag specifying if this is a trajectory-set (>0; 2 means 2nd. order trajectory set...) or not (=0)
													//          [3] framesPerTrajectory, if this feature set is a trajectory set
													//          [4] trajectoryStep, if this feature set is a trajectory set
	long  ByteOrder;			  // init to 0x1, 0x02, 0x03, 0x04
	float	Version;			    // e.g., 1.3
	long	ID;					      // record's class ID
	char	Name[16];			    // record's Label (Name)
	//char	NotUsed[DHLen-36];	// for extension
  unsigned int frameSize; //by thilo: in samples
  unsigned int frameStep; //by thilo: in samples
  unsigned int sampleRate; //by thilo: e.g. 16000 for 16kHz
	char  NotUsed[NotUsedDHLen]; //by thilo
} REC_HDR;

//--------------------------------------------------------- 
// Data Record class
//--------------------------------------------------------- 
class SV_Data {

private :

	void init(SV_Data* pData);		//by thilo, initializes the class as the default-constructor did before
	long ClassSig;              // Class Signature
  bool justLinked;              //by thilo, to remember if the data matrix was just linked with the copy constructor, so it doesn't get destructed if it is borrowed

public :
	REC_HDR	Hdr;				// header info
	float	**Mat;				// data matrix
	int		Row;				// Row of matrix
	int		Col;				// Col of matrix
	SV_Data *Next;				// point to next DataRec

	//----------------------------------
	// Constructor/Destructor
	//----------------------------------
	SV_Data();					// constructor
	SV_Data(int M, int N);		// automatic allocate Matrix
	//SV_Data(const SV_Data& Src); // COPY constructor
  SV_Data(const SV_Data& Src, bool linkOnly = false); //by thilo: copy constructor, that is able to just link to the data-matrix of the parent-object instead of copying it
	SV_Data(const SV_Data& src, int minLen, int maxLen, int minDim, int maxDim); //by tilo: copy constructor to copy a submatrix; min is >=, max is <

	virtual ~SV_Data();			// destructor
	void Alloc(void);			// allocate memory for Mat
	int  Valid(void);			// 1: yes, 0: no. Test Class's Signature

	//----------------------------------
	// overload <<, dump Data in ASCII
	//----------------------------------
	friend ostream& operator<< (ostream& os, SV_Data& Data);

	//----------------------------------
	// overload assignment operator =
	//----------------------------------
	SV_Data &operator= (const SV_Data &Src);

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
  SV_Data* MergeData(unsigned long int maxSegments = 0, int dim = -1);

	void setJustLinked(bool newStatus) {this->justLinked = newStatus; return;} //by thilo: to allow external control of this status (if the datamat needs to be exchanged, for example)
	bool getJustLinked(void) {return this->justLinked;} //by thilo
};

#endif

