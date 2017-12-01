//************************************************************************
//    The implementation of the base class of models
//
//    Author  : Jialong HE
//    Date    : May 3, 1999
//************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SV_Model.h"
#include "SV_Error.h"
#include "GN_Rand.h"
#include "vector.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model::SV_Model() {

	//-------------------------------
	// Default model's header contents
	//-------------------------------
	Hdr.Signature[0]	= 'L';
	Hdr.Signature[1]	= 'H';
	Hdr.Signature[2]	= 19;
	Hdr.Signature[3]	= 95;
	Hdr.Signature[4]	= 1;
	Hdr.Signature[5]	= 30;
	Hdr.Signature[6]	= 0;
	Hdr.Signature[7]	= 0;

	Hdr.ByteOrder		= 0x01020304;
	Hdr.Version			= float(0.1);	
	Hdr.ID				= -1;
	Hdr.ModelType	    = MT_GAUS;

	strcpy(Hdr.Name, "NoName");   // maximum 16 chars
	
	Next = NULL;

	//----------------------
	// Multivariate Gaussian
	//----------------------
	FirstCall = 1; 
	LogTab = NULL;
	
	ClassSig = 1964121;

}

//==========================================
// default destructor
//==========================================
SV_Model::~SV_Model() {
	if (LogTab != NULL) {
		MFree_1D(LogTab); 
	};

	ClassSig = 0;
}

//=======================================================*/
//  Test if current class is valid        
//  return 1: yes, it is valid             
//  return 0: no, it is not valid, maybe deleted.
//=======================================================*/
int SV_Model::Valid (void) {
	if (ClassSig == 1964121) {return (1);}
	else {return(0);}
}


//==========================================
// Close DFile, used for refopen 
//==========================================
void SV_Model::CloseFile(void) {
	DFile.close();	
}

//===============================================
// Open stream DFile for READ_REC or WRITE_REC 
//===============================================
void SV_Model::OpenFile(const char *FName, int Mode) { //by thilo: added const

	if (Mode == WRITE_MODEL) {
		DFile.open(FName, ios::out|ios::binary);  // truncate
	}
	else if (Mode == APPEND_MODEL) {
		DFile.open(FName, ios::app|ios::binary);  // append
	}
	else if (Mode == READ_MODEL) {
		DFile.open(FName, ios::in|ios::binary);   // read
	}
	else {REPORT_ERROR(SVLIB_FileErr, "OpenFile"); }
	
	if (DFile.fail()) {
		REPORT_ERROR(SVLIB_FileErr, "OpenFile");
	}

}

//===============================================
// Save Model's header to current model file 
// if success, return SVLIB_Ok, otherwise return
// SVLIB_Fail
//===============================================
int SV_Model::SaveHdr(void) {

	//------------------------------------
	// Write Model Hdr into file
	//------------------------------------
	//DFile.write((char*)(&Hdr), MHLen);
	//by thilo:
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes sizes;
	io.getCurrentDatatypeSizes(sizes);
	int bytes = io.writeMachineHeader(&(this->DFile), sizes);
	bytes += io.writeArray(&(this->DFile), this->Hdr.Signature, 8);
	bytes += io.writeScalar(&(this->DFile), this->Hdr.ByteOrder);
	bytes += io.writeScalar(&(this->DFile), this->Hdr.Version);
	bytes += io.writeScalar(&(this->DFile), this->Hdr.ID);
	bytes += io.writeScalar(&(this->DFile), this->Hdr.ModelType);
	bytes += io.writeArray(&(this->DFile), this->Hdr.Name, 16);
	bytes += io.writeArray(&(this->DFile), this->Hdr.NotUsed, NotUsedMHLen);
	//end by thilo
	if (DFile.good() != TRUE) {
		return (SVLIB_Fail);
	}
	else {
		return (SVLIB_Ok);
	}
}

//===============================================
// LoadModel's header to current model file 
// if success, return SVLIB_Ok, otherwise return
// SVLIB_Fail
//===============================================
int SV_Model::LoadHdr(SV_DataIO::SV_DatatypeSizes &sizes) { //by thilo: returns also machine-dependant header

	//------------------------------------
	// Load Model Hdr into file
	//------------------------------------
	//DFile.read((char*)(&Hdr), MHLen);
	//by thilo:
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes;
	io.getCurrentDatatypeSizes(codeSizes);
	int bytes = io.readMachineHeader(&(this->DFile), sizes, true);
	if (bytes > 0) {
		io.consumeBytes(&(this->DFile), bytes);
	} else {
		bytes = 0;
	}
	bytes += io.readArray(&(this->DFile), this->Hdr.Signature, 8, codeSizes, sizes);
	bytes += io.readScalar(&(this->DFile), this->Hdr.ByteOrder, codeSizes, sizes);
	bytes += io.readScalar(&(this->DFile), this->Hdr.Version, codeSizes, sizes);
	bytes += io.readScalar(&(this->DFile), this->Hdr.ID, codeSizes, sizes);
	bytes += io.readScalar(&(this->DFile), this->Hdr.ModelType, codeSizes, sizes);
	bytes += io.readArray(&(this->DFile), this->Hdr.Name, 16, codeSizes, sizes);
	bytes += io.consumeArray(&(this->DFile), this->Hdr.NotUsed, getNotUsedHeaderSize(sizes), codeSizes, sizes);
	//end by thilo
	if (DFile.good() != TRUE) {
		return (SVLIB_Fail);
	}
	else {
		return (SVLIB_Ok);
	}
}

//===============================================
// by thilo: LoadModel's header to current model 
// file; if success, return SVLIB_Ok, otherwise 
// return SVLIB_Fail
//===============================================
int SV_Model::oldLoadHdr(SV_DataIO::SV_DatatypeSizes &sizes) {

	//------------------------------------
	// Load Model Hdr into file
	//------------------------------------
	SV_DataIO io;
	io.readMachineHeader(&(this->DFile), sizes, true); //if there is no header (what we assume), a standard one is returned
	DFile.read((char*)(&Hdr), MHLen);
	if (DFile.good() != TRUE) {
		return (SVLIB_Fail);
	}
	else {
		return (SVLIB_Ok);
	}
}

//====================================================================================================================
//  by thilo: Returns the size of the NotUsed-member of the header-struct based on the knowledge of the number and 
//  type of members and their sizes as given in "sizes"
//====================================================================================================================
int SV_Model::getNotUsedHeaderSize(SV_DataIO::SV_DatatypeSizes sizes) {
	int bytes = MHLen - 8*sizes.charSize - 3*sizes.longSize - sizes.floatSize -  16*sizes.charSize; //attention: thi suses internal knowledge of how the header is constructed - only co-change!

	return bytes;
}

//===============================================
// Merge a linked list into one SV_Data with Mat
// contains concatenate of all Mats in pData
//===============================================
SV_Data *SV_Model::MergeData(SV_Data *pData) {

	SV_Data *pTData, *DataCurr = NULL;

	int TotalVec = 0;
	int Row, Col, Dim;

	//------------------------------------
	// for safety, check input variables
	//------------------------------------
	if (pData == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "MergeData"); 	
	}

	DataCurr = pData;
	Dim		 = pData->Col;
	//------------------------------------
	// Count total vectors
	//------------------------------------
	while (DataCurr != NULL) {
		if (Dim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "MergeData"); 	
		}

		TotalVec += DataCurr->Row; 
		DataCurr = DataCurr->Next;
	}  // while

	//------------------------------------
	// Allocate memory for New SV_Data
	//------------------------------------
	pTData = new SV_Data;
	pTData->Row = TotalVec;
	pTData->Col = Dim;
  pTData->Alloc();
  pTData->Hdr.frameSize = pData->Hdr.frameSize; //by thilo
  pTData->Hdr.frameStep = pData->Hdr.frameStep; //by thilo
  pTData->Hdr.ID = pData->Hdr.ID; //by thilo

	DataCurr = pData;
	TotalVec = 0;
	//------------------------------------
	// Accumulate patterns
	//------------------------------------
	while (DataCurr != NULL) {

		//------------------------------------
		// accumulate vectors in this pattern
		//------------------------------------
		for (Row=0; Row<DataCurr->Row; Row++) {
			for (Col=0; Col<Dim; Col++) {
				pTData->Mat[TotalVec+Row][Col] = DataCurr->Mat[Row][Col];
			}
		}
		TotalVec += DataCurr->Row; 
		DataCurr = DataCurr->Next;
	}  // while

	return(pTData);
}



//==========================================
// calculate mean vector of linked list
//==========================================
double *SV_Model::MeanVec(SV_Data *pData) {

	SV_Data *DataCurr = NULL;
	int TotalVec = 0;
	int Row, Col, Dim;
	double *MVec;

	//------------------------------------
	// for safety, check input variables
	//------------------------------------
	if (pData == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "MeanVec"); 	
	}

	DataCurr = pData;
	Dim		 = pData->Col;
	MArray_1D(MVec, Dim, double, "MeanVec: MVec");

	//------------------------------------
	// clear MVec to zero
	//------------------------------------
	for (Col=0; Col<Dim; Col++) {
		MVec[Col] = 0.0;
	}

	//------------------------------------
	// Accumulate patterns
	//------------------------------------
	while (DataCurr != NULL) {

		TotalVec += DataCurr->Row; 
		if (Dim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "Different Dim"); 	
		}

		//------------------------------------
		// accumulate vectors in this pattern
		//------------------------------------
		for (Row=0; Row<DataCurr->Row; Row++) {
			for (Col=0; Col<Dim; Col++) {
				MVec[Col] += DataCurr->Mat[Row][Col];
			}
		}
		DataCurr = DataCurr->Next;
	}  // while

	//------------------------------------
	// normalized by total vectors
	//------------------------------------
	for (Col=0; Col<Dim; Col++) {
		MVec[Col] /= TotalVec;
	}

	return (MVec);
};


//==========================================
// calculate covariance matrix of linked list 
//==========================================
double **SV_Model::CovMatrix(SV_Data *pData) {

	SV_Data *DataCurr = NULL;
	double **CMat, *Diff, *MVec;
	int TotalVec = 0;
	int Row, Col, Cnt, Dim;

	//------------------------------------
	// for safety, check input variables
	//------------------------------------
	if (pData == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "CovMatrix"); 	
	}

	DataCurr = pData;
	Dim		 = pData->Col;
	MArray_2D(CMat, Dim, Dim, double, "CovMatrix: CMat");
	MArray_1D(Diff, Dim, double, "CovMatrix: Diff");

	//------------------------------------
	// clear CMat to zero
	//------------------------------------
	for (Row=0; Row<Dim; Row++) {
		for (Col=0; Col<Dim; Col++) {
			CMat[Row][Col] = 0.0;
		}
	}

	MVec = MeanVec(pData);   // obtain mean vector
	//------------------------------------
	// Accumulate patterns
	//------------------------------------
	while (DataCurr != NULL) {

		TotalVec += DataCurr->Row; 
		if (Dim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "Dim difference"); 	
		}

		/*------------------------------------*/
		/* calculate the covariance matrix    */
		/*------------------------------------*/
		for (Row=0; Row<DataCurr->Row; Row++) {
			for (Col=0; Col<Dim; Col++) {
				Diff[Col] = DataCurr->Mat[Row][Col] - MVec[Col];
			}

			for (Col=0; Col<Dim; Col++) {
				for (Cnt=0; Cnt<Dim; Cnt++) {
					CMat[Col][Cnt] +=  Diff[Col] * Diff[Cnt];
				}
			}
		}
		DataCurr = DataCurr->Next;
	}  // while

	//------------------------------------
	// normalize
	//------------------------------------
	for (Row=0; Row<Dim; Row++) {
		for (Col=0; Col<Dim; Col++) {
			CMat[Row][Col] /= TotalVec;
		}
	}

    MFree_1D(MVec);
    MFree_1D(Diff);
	return(CMat);
};


//============================================================
//  evaluate Multivariate Gaussian function 
//============================================================
double SV_Model::Gaussian (double *Mean, double *Vari, float *Vector, int Dim) {


   int DimCnt, Ind, Cnt;
   double X, Dif, Prob;
   double Sum = 0.0, Det = 0.0;

   //--------------------------------------------
   // for speed reason, using table looking 
   // rather than calculate log each time
   //--------------------------------------------
   if (FirstCall) {      
     FirstCall = 0;
     MArray_1D(LogTab, ITEMS, double, "LogTab");
     X = STEPSIZE;
     for (Cnt=0; Cnt<ITEMS; Cnt++) {
       LogTab[Cnt] = log(X);
       X += STEPSIZE;
     }
   }

   for (DimCnt=0; DimCnt<Dim; DimCnt++) {
     Dif  = Mean[DimCnt] - Vector[DimCnt];
     Sum += Dif * Dif / Vari[DimCnt];
     if (Vari[DimCnt] >= RANGE) Det += log(Vari[DimCnt]);
     else {
       Ind = (int) (Vari[DimCnt] / STEPSIZE);  /* Table Looking */
       Det += LogTab[Ind];
     }
   }

   Prob = (double)Dim * log(8.0 * atan(1.0)) + Det + Sum;
   if (Prob > 600) Prob = 600;
   else if (Prob < -600) Prob = -600;

   return (exp(-Prob/2.0));
}



//=================================================
// Orthogonal transform data by TranMat
//  
//=================================================
void SV_Model::OrthTrans(SV_Data *pData, double **TransMat) {

	SV_Data *DataCurr = NULL;
	int Row, Col, DimCnt, Dim;
	double *OneVec;

	//------------------------------------
	// for safety, check input variables
	//------------------------------------
	if (pData == NULL || TransMat == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "MeanVec"); 	
	}

	DataCurr = pData;
	Dim		 = pData->Col;
	
	MArray_1D(OneVec, Dim, double, "OrthTrans:: OneVec");

	//------------------------------------
	// Accumulate patterns
	//------------------------------------
	while (DataCurr != NULL) {

		if (Dim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "Different Dim"); 	
		}

		//------------------------------------
		// orthogonal transform this vector
		//------------------------------------
		for (Row=0; Row<DataCurr->Row; Row++) {
			for (Col=0; Col<Dim; Col++) {
				OneVec[Col] = 0.0;
				for (DimCnt=0; DimCnt<Dim; DimCnt++) {
					 OneVec[Col] +=  DataCurr->Mat[Row][DimCnt] * TransMat[DimCnt][Col];
				}
			}

			/*----------------------------------*/
			/* Copy back transformed vector     */
			/*----------------------------------*/
			for (Col=0; Col<Dim; Col++) {
				DataCurr->Mat[Row][Col] = (float)OneVec[Col];
			}
		}  // for (Row)


		DataCurr = DataCurr->Next;
	}  // while

	MFree_1D(OneVec);
};



/*========================================================*/
/*  This procedure find the nearest code vector           */
/*  and return the index of this code vector              */
/*========================================================*/
int SV_Model::nearest_code (float **CodeBook, float *testvector, int CodeNum, int Dim) {

	float MinDist = LARGE, CurrDist;
	int   CodeCnt, MinIndex, CurrIndex;

	CurrIndex = 0;
	for (CodeCnt = 0; CodeCnt<CodeNum; CodeCnt++) {

		DIST_VEC(CodeBook[CodeCnt], testvector, CurrDist, Dim, float, float);
		if (CurrDist < MinDist) {
			MinDist  = CurrDist;
			MinIndex = CurrIndex;
		}  /* if */
		CurrIndex++;
	}  /* for */

	return(MinIndex);
}

//==========================================
//  Count patterns in the linked list 
//==========================================
int SV_Model::CountPat(SV_Data *pData) {

	SV_Data *DataCurr = pData;
	int PatNum = 0;
	
	//-----------------------------
	// Count number of patterns
	//-----------------------------
	while (DataCurr != NULL) {
		PatNum++;
		DataCurr = DataCurr->Next;
	}

	return (PatNum);
}




/*================================================================*/
/* This procedure first ramp a array from 0 to MaxNum-1           */
/* and then random the sequence of this array. This function      */
/* will be used to random the order of a sequence                 */
/* a int array will allocated which can be released by free()     */
/*                                                                */
/*                                                                */
/* Author : Jialong                                               */
/* Date   : July 29, 1994                                         */
/*================================================================*/
int *SV_Model::RandSeq(int MaxNum, int Seed) {

  int  MaxIndex, Cnt;
  int  RandPnt, FindPnt, RealPnt;
  int *RndVec, *IndexArry;

  GN_Rand RandEng;
  RandEng.srandom(Seed);
  
  MArray_1D(IndexArry, MaxNum, int, "IndexArray");
  MArray_1D(RndVec, MaxNum, int, "RndVec");
  for (Cnt=0; Cnt<MaxNum; Cnt++) IndexArry[Cnt] = Cnt;

  MaxIndex = MaxNum;

  for (Cnt=0; Cnt<MaxNum; Cnt++) {

     RandPnt = RandEng.random() % MaxIndex;
     FindPnt = 0;
     RealPnt = 0;

     while ((FindPnt<=RandPnt) && (RealPnt<MaxNum)) {
	  if (IndexArry[RealPnt] != -1) FindPnt++;
	  RealPnt++;
     };     /* While */

  RealPnt = RealPnt - 1;         /* eliminating the effect of last RealRnt++*/
  RndVec[Cnt] = RealPnt;
  IndexArry[RealPnt] = -1;       /* Exclude this vector from next search */
  MaxIndex = MaxIndex - 1;

  } /* for Cnt */

  MFree_1D(IndexArry);
  return (RndVec);

}


//==========================================
// base function  
//==========================================
int SV_Model::TrainModel(SV_Data *pData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	
	cerr << "TrainModel must be overloaded!";

  return -1;
}

//==========================================
// base function  
//==========================================
SV_Data* SV_Model::TestModel(SV_Data *pData) {
	
	cerr << "TestModel must be overloaded!";
	return(NULL);
}

//==========================================
// base function  
//==========================================
int SV_Model::SaveModel(void) {
	
	cerr << "SaveModel must be overloaded!";
	return(0);
}

//==========================================
// base function  
//==========================================
SV_Model *SV_Model::LoadModel(void) {
	
	cerr << "LoadModel must be overloaded!";
	return (NULL);
}

