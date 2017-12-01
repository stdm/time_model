//************************************************************************
//    Implementation of uni-modal multi-variate Gaus model
//                                                                
//                -p/2        -1/2             T  -1              
//     G(x) = (2pi)    * Det(C)     * exp((x-u)  C  (x-u) / -2)   
//
//       where C is the covariance matrix, u is the mean vector.                
//                                                                
//
//    Author  : Jialong HE
//    Date    : May 11, 1999
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "SV_Model_Gaus.h"
#include "SV_Error.h"
#include "GN_Matrix.h"
#include "vector.h"


static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model_Gaus::SV_Model_Gaus() {
	
	ICov = NULL;
	MVec = NULL;
	Dim  = 0;
}

//==========================================
// default destructor
//==========================================
SV_Model_Gaus::~SV_Model_Gaus() {
	if (ICov != NULL) {
		MFree_2D(ICov);
	}

	if (MVec != NULL) {
		MFree_1D(MVec);
	}

}


//===============================================
// Generate unimodal, multvariate Gaussian model
// a gaussian model is characterised by 
// (1) mean vector, (2) inverse covariance matrix
// (3) det of covariance matrix 
//===============================================
int SV_Model_Gaus::TrainModel(SV_Data *pData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)

	double **MCov, *pMVec;
	GN_Matrix MEng;
	int RtCode;

	Dim  = pData->Col;			  // vector dimension
	//-----------------------------------
	// Obtain mean vectors
	//-----------------------------------
	pMVec= MeanVec(pData);		  // Mean vector
	if (MVec != NULL) {MFree_1D(MVec);}
	MVec = pMVec;

	//-----------------------------------
	// Obtain inverse of covariance mat
	//-----------------------------------
	MCov = CovMatrix(pData);          // get covariance matrix	
	DetV = MEng.Det (MCov, Dim);      // Det of MCov

	RtCode = MEng.Inv (MCov, Dim);    // Inv of MCov

	if (RtCode == 0)   { // success
		if (ICov != NULL) {MFree_2D(ICov);}
		ICov = MCov;  // MCov allocated in CovMatrix;	
	}
	else {
		REPORT_ERROR(SVLIB_NoInv, "Train GausModel");
	}

  return 0;
}

//===============================================
// Save the model to current model file
// if success, return total bytes writed,
// otherwise, REPORT_ERROR
// by thilo: changed all write()-calls
//===============================================
int SV_Model_Gaus::SaveModel(void) {

	int RtCode;
	//int TotalByte; //by thilo
	int bytes; //by thilo
	SV_DataIO io; //by thilo

	Hdr.ModelType = MT_GAUS;
    RtCode = SaveHdr();
	if (RtCode == SVLIB_Fail) {
		REPORT_ERROR(SVLIB_Fail, "Save Gaus Model:");
	}
	
	//--------------------------
	// Gaus model's parameters
	//--------------------------
	//DFile.write((char*)(&Dim), sizeof(int));
	bytes = io.writeScalar(&(this->DFile), this->Dim);
	//DFile.write((char*)(&DetV), sizeof(double));
	bytes += io.writeScalar(&(this->DFile), this->DetV);
	//DFile.write((char*)(MVec), Dim*sizeof(double));
	bytes += io.writeArray(&(this->DFile), this->MVec, this->Dim);
	//DFile.write((char*)(ICov[0]), Dim*Dim*sizeof(double));
	bytes += io.writeMatrix(&(this->DFile), this->ICov, this->Dim, this->Dim);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save Gaus Model:");
	}

	//TotalByte = MHLen + sizeof(int) + sizeof(double)*(Dim+1) +  Dim*Dim*sizeof(double);
	//return(TotalByte);
	return bytes + MHLen; //MHLen maybe an incorrect size...
}

//==============================================================
// Read model parameters from current opened model file
// if success, return (this) pointer,  if fail, return (NULL)
// by thilo: changed all read()-calls
//==============================================================
SV_Model *SV_Model_Gaus::LoadModel(void) {

	int RtCode;
	int NewDim;
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
	if (Hdr.ModelType != MT_GAUS) {
		return(NULL);
	}

	//--------------------------
	// Gaus model's parameters
	//--------------------------
	//DFile.read((char*)(&NewDim), sizeof(int));
	bytes = io.readScalar(&(this->DFile), NewDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewDim == 0) {
		return(NULL);
	}

	//DFile.read((char*)(&DetV), sizeof(double));
	bytes += io.readScalar(&(this->DFile), DetV, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {
		return(NULL);
	}

	//--------------------------------------------
	// Load Mean Vec and Inv Covariance matrix
	//--------------------------------------------
	if (NewDim != Dim) {
		MFree_2D(ICov); 
		MFree_1D(MVec);

		Dim = NewDim;
		MArray_1D(MVec, Dim, double, "Gaus LoadModel");
		MArray_2D(ICov, Dim, Dim, double, "Gaus LoadModel");
	}

	//DFile.read((char*)(MVec), Dim*sizeof(double));
	bytes += io.readArray(&(this->DFile), this->MVec, this->Dim, codeSizes, fileSizes);
	//DFile.read((char*)(ICov[0]), Dim*Dim*sizeof(double));
	bytes += io.readMatrix(&(this->DFile), this->ICov, this->Dim, this->Dim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {
		return(NULL);
	}

	return(this);
}

//===============================================
// Test unimodal, multvariate Gaussian model
//===============================================
SV_Data *SV_Model_Gaus::TestModel(SV_Data *pData) {

	SV_Data *DataCurr, *pOnePat;
	SV_Data *Score;

	int PatNum = 0, PatCnt;

	double *DataMVec, **DataCMat, LL;

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

		pOnePat  = DataCurr; 
		//---------------------------------------------------
		// since CovMatrix, MeanVec calculate linked list
		// Here I just want one pattern, set Next = NULL
		//---------------------------------------------------
		DataCurr = DataCurr->Next;		   // point to next
		pOnePat ->Next = NULL;			  // 
		DataCMat = CovMatrix(pOnePat);    // allocate CMat 
		DataMVec = MeanVec(pOnePat);      // allocate MVec
     
		//------------------------------
		// Log likelihood
		//------------------------------
		LL = AvgL(ICov, MVec, DataCMat, DataMVec, Dim);
    LL = -0.5 * (LL + log(DetV) + Dim*log(6.28));
		Score->Mat[PatCnt][0] = float(LL);
	
		MFree_1D(DataMVec);
		MFree_2D(DataCMat);
	}

	return(Score);
}
	

//--------------------------------------------------------*/
// Calculate Score L = tr(Cy Cx) + (My-Mx)' Cx (My - Mx)   */
//--------------------------------------------------------*/
double SV_Model_Gaus::AvgL(double **Cx, double *Mx, double **Cy, double *My, int Dim) {

    double *Diff;
    double  Right, AccuLi = 0.0;
    int Row, Col;
	int ScoreType = 0;   // 0: use all, 1: just mean, 2: just covariance

    MArray_1D(Diff, Dim, double, "Diff");
    SUB_VEC(My, Mx, Diff, Dim, double, double, double);

    /*-------------------------*/
    /* Right =  Cx * (My-Mx)   */
    /*-------------------------*/
    if (ScoreType == 0 || ScoreType == 1)
      for (Row=0; Row<Dim; Row++) {
		DOTP_VEC(Diff, Cx[Row], Right, Dim, double, double);
        AccuLi += Right * Diff[Row];
      }

    /*-------------------------*/
    /* Trace(Cy * Cx)          */
    /*-------------------------*/
    if (ScoreType == 0 || ScoreType == 2)
		for (Row=0; Row<Dim; Row++) {
			for (Col=0; Col<Dim; Col++) {
				AccuLi += Cy[Row][Col] * Cx[Col][Row];
			}
		}

    MFree_1D(Diff);
    return (AccuLi);
};


//=============================================
// Dump a data record in ASCII ot cout
//=============================================
ostream& operator<< (ostream& OutS, SV_Model_Gaus& Data) {

	int Row, Col;

	OutS << "Dim: "<<Data.Dim << endl;
	OutS << "Det: "<<Data.DetV << endl;

	OutS << "--- Mean Vector---"<< endl;
	for (Col=0; Col<Data.Dim; Col++) {
		OutS<<Data.MVec[Col]<< " ";
	}

	OutS << endl << "--- Inv Cov --"<< endl;
	for (Row=0; Row<Data.Dim; Row++) {
		for (Col=0; Col<Data.Dim; Col++)
			OutS<<Data.ICov[Row][Col]<< " ";
		OutS<<endl;
	}
	
	return(OutS);
}
