//************************************************************************
//    Implementation of VQ model
//
//
//    Author  : Jialong HE
//    Date    : May 11, 1999
//************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include "SV_Model_VQ.h"
#include "SV_Error.h"
#include "GN_Rand.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model_VQ::SV_Model_VQ() {
	CBook = NULL;
	CNum  = 0;
	Dim   = 0;
	Verbose = 1;
	SplitMethod = 0;
	MaxIter = 100;
	CBSize  = 8;
	RandSeed = 1999;
}

//==========================================
// default destructor
//==========================================
SV_Model_VQ::~SV_Model_VQ() {
	if (CBook != NULL) {
		MFree_2D(CBook);
	}
}

//===============================================
// Save the model to current model file
// if success, return total bytes writed,
// otherwise, REPORT_ERROR
// by thilo: replaced all calls to the write()-
//           function
//===============================================
int SV_Model_VQ::SaveModel(void) {

	int RtCode;
	//int TotalByte; //by thilo
	int bytes; //by thilo
	SV_DataIO io; //by thilo

	Hdr.ModelType = MT_VQ;
    RtCode = SaveHdr();
	if (RtCode == SVLIB_Fail) {
		REPORT_ERROR(SVLIB_Fail, "Save VQ model Failed!");
	}
	
	//--------------------------
	// VQ model's parameters
	//--------------------------
	//DFile.write((char*)(&Dim), sizeof(int));
	bytes = io.writeScalar(&(this->DFile), this->Dim);
	//DFile.write((char*)(&CNum), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->CNum);
	//DFile.write((char*)(CBook[0]), CNum*Dim*sizeof(float));
	bytes += io.writeMatrix(&(this->DFile), this->CBook, this->CNum, this->Dim);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save VQ Model Failed!");
	}

	//TotalByte = MHLen + 2*sizeof(int) +  CNum*Dim*sizeof(float);
	//return(TotalByte);
	return bytes + MHLen; //MHLen may be incorrect...
}

//===========================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
// by thilo: replaced all read()-calls
//===========================================================
SV_Model * SV_Model_VQ::LoadModel(void) {

	int RtCode, NewDim, NewSize;
	int bytes; //by thilo
	SV_DataIO io; //by thilo
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes; //by thilo
	io.getCurrentDatatypeSizes(codeSizes); //by thilo

    RtCode = LoadHdr(fileSizes); //by thilo: catch machine-dependant header
	if (RtCode == SVLIB_Fail) {
		return(NULL);
	}

	//--------------------------
	// check if valid header
	//--------------------------
	if (Hdr.ModelType != MT_VQ) {
		return(NULL);
	}

	//--------------------------
	// VQ model's parameters
	//--------------------------
	//DFile.read((char*)(&NewDim), sizeof(int));
	bytes = io.readScalar(&(this->DFile), NewDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewDim == 0) {
		return(NULL);
	}

	
	//DFile.read((char*)(&NewSize), sizeof(int));
	bytes += io.readScalar(&(this->DFile), NewSize, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewSize == 0) {
		return(NULL);
	}

	//--------------------------------------------
	// Allocate CBook
	//--------------------------------------------
	if (NewDim != Dim || NewSize != CNum) {
		MFree_2D(CBook);
		CNum = NewSize;
		Dim  = NewDim;
		MArray_2D(CBook, CNum, Dim, float, "LoadVQ");

	}; 
	
	//DFile.read((char*)(CBook[0]), CNum*Dim*sizeof(float));
	bytes += io.readMatrix(&(this->DFile), this->CBook, this->CNum, this->Dim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {
		return(NULL);
	}

	return(this);
}


//===============================================
// Generate Codebook with the LBG algorithm
//===============================================
int SV_Model_VQ::TrainModel(SV_Data *pData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	
	SV_Data *pTData;

	if (pData == NULL) {
		//REPORT_ERROR(SVLIB_BadArg, "VQ TrainModel"); //by thilo
		return SVLIB_Fail; //by thilo
	}

	//-------------------------------
	// allocate memory for codebook
	//-------------------------------
	if (CBook != NULL) {MFree_1D(CBook);}
	Dim = pData->Col;
	CNum = CBSize;      // num of codes
	MArray_2D(CBook, CNum, Dim, float, "VQ: TrainModel");
		
	//-------------------------------------------
	// merge link list pData into one 
	//-------------------------------------------
	pTData = MergeData(pData);  // remember release pTData

	if (SplitMethod) {
        Split_Codebook(CBook, pTData->Mat,  CNum, pTData->Row, Dim);
	}
	else {
		GN_Rand REng;
		double *MVec;
		int Row, Col;
		double Purt;

		REng.srandom(RandSeed);
		//-------------------------------------
		// Init codebook = mean + rand
		//-------------------------------------
		MVec = MeanVec(pTData);
		for (Row=0; Row<CNum; Row++) {
			for (Col=0; Col<Dim; Col++) {
				Purt = double (REng.random()) / double(REng.getmax() - 0.5);
				CBook[Row][Col] = float(MVec[Col] + Purt*pTData->Mat[Row][Col]);
			}
		}

        LBG_Codebook(CBook, CNum, pTData->Mat, pTData->Row, Dim);
		MFree_1D(MVec);
	}

	delete(pTData);

  return 0;
}

//===============================================
// Test codebook
//===============================================
SV_Data* SV_Model_VQ::TestModel(SV_Data *pData) {

	SV_Data *Score, *DataCurr;
	int PatNum, NInd, PatCnt, VecCnt;
	float AvgDist, Dist;

	if (pData->Col != Dim) {
		REPORT_ERROR(SVLIB_BadData, "Model and Data have different Dim!");
	}

	PatNum = CountPat(pData);
	//----------------------------
	// allocate score space
	//----------------------------
	Score = new SV_Data;
	Score->Row = PatNum;	
	Score->Col = 1;	         // one score for each pattern
	Score->Alloc();
	
	//----------------------------
	// test score for each pattern
	//----------------------------
	DataCurr = pData;
	for (PatCnt=0; PatCnt<PatNum; PatCnt++) {
		AvgDist = 0.0;
		for (VecCnt=0; VecCnt<DataCurr->Row; VecCnt++) {
			NInd = winning_cell(CBook, DataCurr->Mat[VecCnt], CNum, Dim, &Dist);
			AvgDist += float(sqrt(Dist));	// sum of Eucledian distance.		
		}
		
		Score->Mat[PatCnt][0] = AvgDist / DataCurr->Row;
		DataCurr = DataCurr->Next;		   // point to next
	}

	return (Score);
}

//===============================================
// Generate Codebook with spliting method
//===============================================
void SV_Model_VQ::Split_Codebook(float **CodeBook, float **DataMat,  int CodeNum, int DataRow, int DataDim) {

  int   DataCnt, DimCnt;        /* loop counter variables      */
  int   CdCnt, CurrCode;
  float Orig, Sum;

  /*----------------------------------------------*/
  /* First code vector is centroid of all vectors */
  /*----------------------------------------------*/
  for (DimCnt=0; DimCnt<DataDim; DimCnt++) {
    Sum = 0.0;
    for (DataCnt=0; DataCnt<DataRow; DataCnt++)
       Sum += DataMat[DataCnt][DimCnt];
    CodeBook[0][DimCnt] = Sum / (float) DataRow;
  }

  CurrCode = 1;  /* this one code vector */
  /*------------------------------------*/
  /* loop until splitting is sufficient */
  /*------------------------------------*/
  while (CurrCode * 2 <= CodeNum) {

    /*-----------------------------*/
    /* split codes by perturbation */
    /*-----------------------------*/
    for (CdCnt=0; CdCnt<CurrCode; CdCnt++)
      for (DimCnt=0; DimCnt<DataDim; DimCnt++) {
	 Orig = CodeBook[CdCnt][DimCnt];
	 CodeBook[CdCnt + CurrCode][DimCnt] = Orig + float(0.00001) * Orig;
	 CodeBook[CdCnt][DimCnt] = Orig - float(0.00001) * Orig;
      }

    CurrCode *= 2;    /* double code book size */

    LBG_Codebook(CodeBook, CurrCode, DataMat, DataRow, DataDim);
  } /* while */
}

//==========================================
// This the Engine of LBG training algorithm
//==========================================
float SV_Model_VQ::LBG_Codebook(float **codebook, int numcodes, float **data, int numdata, int dimension) {

  int   CodeCnt, VecCnt;        /* loop counter variables      */
  int   count = 0;
  int   *Labels;
  float ThisError, LastError=LARGE;
  float Dist;


  MArray_1D(Labels, numdata, int, "labels");

  /*------------------------------------*/
  /* loop until splitting is sufficient */
  /*------------------------------------*/
  while (count<MaxIter) {

    ThisError = 0.0;
    /*------------------------*/
    /* computer winning cells */
    /*------------------------*/
    for (VecCnt=0; VecCnt<numdata; VecCnt++) {
      Labels[VecCnt] = winning_cell(codebook, data[VecCnt], numcodes, dimension, &Dist);
      ThisError += Dist * Dist;
    }

    ThisError /= ((float) numdata  * (float) dimension);
    if (LastError - ThisError<TINY) break;

    for (CodeCnt=0; CodeCnt<numcodes ; CodeCnt++) {
      centroid(data, codebook[CodeCnt], Labels, numdata, dimension, CodeCnt);
	}

    count++;
    LastError = ThisError;
    if (Verbose) {
		cerr << ThisError << endl;
	}
  
  }

  MFree_1D(Labels);
  return(ThisError);
}

//---------------------------------------------
//  find centroid of data belong to spcified code
//---------------------------------------------
void SV_Model_VQ::centroid (float **data, float *centvec, int *LabelVec, int numdata, int dimension, int Label) {

  int VecCnt, DimCnt;
  int numlabel = 0;
	
  for (DimCnt=0; DimCnt< dimension; DimCnt++)
	  centvec[DimCnt] = 0.0;

  /*------------------------------------------------*/
  /* Accumulate vectors belong to this code vector  */
  /*------------------------------------------------*/

  for (VecCnt=0; VecCnt<numdata; VecCnt++) {
		if (LabelVec[VecCnt] == Label) {
			for (DimCnt=0; DimCnt<dimension; DimCnt++) {
				centvec[DimCnt] += data[VecCnt][DimCnt];
			}
			numlabel++;          /* Vector counter increase */
		}
	}

	if (numlabel == 0) {
		numlabel = 1;
	}

	for (DimCnt=0; DimCnt<dimension; DimCnt++) {
      centvec[DimCnt] /= numlabel;
	}
}


/*
 *==================================================
 *  Find out nearest code vector for given OneVec
 *  Return the index of this code vector
 *==================================================
 */
int SV_Model_VQ::winning_cell(float **CodeBook, float *OneVec, int numcv, int dimension, float *Dist) {

  int   CodeCnt, DimCnt, Nearest=0;
  float EucliDist,  MinError=LARGE;

  for (CodeCnt = 0; CodeCnt < numcv; CodeCnt++) {
    EucliDist = 0.0;
    for (DimCnt = 0; DimCnt < dimension; ++DimCnt)
       EucliDist += (CodeBook[CodeCnt][DimCnt] - OneVec[DimCnt]) *
	       (CodeBook[CodeCnt][DimCnt] - OneVec[DimCnt]);
    if (EucliDist < MinError) {
       MinError = EucliDist;
       Nearest  = CodeCnt;
    }
  }    /* CodeCnt */

  (*Dist) = MinError;
  return (Nearest);

}

//=============================================
// Dump model's parameter in ASCII 
//=============================================
ostream& operator<< (ostream& OutS, SV_Model_VQ& Data) {

	int Row, Col;

	OutS << Data.CNum <<" "<<Data.Dim << endl;
	for (Row=0; Row<Data.CNum; Row++) {
		for (Col=0; Col<Data.Dim; Col++)
			OutS<<Data.CBook[Row][Col]<< " ";
		OutS<<endl;
	}
	
	return(OutS);
}



