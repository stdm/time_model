//************************************************************************
//    Implementation of OGMM
//
//
//    Author  : Jialong HE
//    Date    : May 11, 1999
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include "SV_Model_GMM.h"
#include "SV_Error.h"
#include "GN_Matrix.h"
#include "GN_Rand.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model_GMM::SV_Model_GMM() {

	OrthT = 1;
	OMat  = NULL; 

	MixNum  = 0; 
	VecDim  = 0; 
	WgtVec  = NULL; 
	MeanMat = NULL; 
	VariMat = NULL;


	CMat    = NULL;
	MVec    = NULL;

	//--------------------------------------
	// Default paras for generate model
	//--------------------------------------
	Verbose  = 1;  // display messages
	Mixtures = 8;  // number of mixtures
	WithOrth = 1;  // with orthogonal transform
	MaxIter  = 100;  // with orthogonal transform
	SplitMethod = 0;
	RandSeed = 1999;

}

//==========================================
// default destructor
//==========================================
SV_Model_GMM::~SV_Model_GMM() {
	
	FreeModel();

}

//==============================================
// Depends on MixNum, VecDim,  allocate memory
// for OMat, WgtVec, MeanMat, VariMat
//==============================================
void SV_Model_GMM::AllocModel(void) {

	if (MixNum != 0 && VecDim != 0) {
		MArray_2D(OMat, VecDim, VecDim, double, "OMat");	
		MArray_2D(MeanMat, MixNum, VecDim, double, "MeanMat");	
		MArray_2D(VariMat, MixNum, VecDim, double, "MeanMat");	
		MArray_1D(WgtVec, MixNum, double, "WgtVec");	
	}
	else {
		REPORT_ERROR(SVLIB_BadArg, "GMM:: AllocModel");
	}	
};	 

//==============================================
// Release all dynamicly allocate memory
//==============================================
void SV_Model_GMM::FreeModel(void) {

	MFree_1D(WgtVec);
	MFree_2D(OMat);
	MFree_2D(MeanMat);
	MFree_2D(VariMat);
	MixNum = 0;
	VecDim = 0;
};

//===============================================
// Generate OGMM
//===============================================
int SV_Model_GMM::TrainModel(SV_Data *pData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	
	SV_Data *pTData;

	int Row, Col, MixCnt, DimCnt, DataDim, CurrMix;

	if (pData == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "GMM TrainModel");
	}

	DataDim = pData->Col;
  OrthT   = WithOrth;
	//-------------------------------------
	// re-allocate memory if necessary
	//-------------------------------------
	if (VecDim != DataDim || MixNum != Mixtures) {
		FreeModel();
		VecDim = DataDim;
		MixNum = Mixtures;
		AllocModel();
	}

	//-------------------------------------------
	// merge linked list pData into one pattern 
	//-------------------------------------------
	pTData = MergeData(pData);  // remember release pTData
	CMat   = CovMatrix(pTData);   // overall covariance matrix

	//--------------------------------
	// Orthogonal transform data
	//--------------------------------
	if (OrthT) {
		GN_Matrix MEng;
		double *EigVal;

		MArray_1D(EigVal, VecDim, double, "TrainModel::EigVal");
		OrthT = 1;
		//----------------------------------------
		// find orthogonal trans mat OMat and 
		// transform pTData
		//----------------------------------------
		MEng.Eigen(CMat, OMat, EigVal, VecDim);  
		OrthTrans(pTData, OMat);
		
		//---------------------------------
		// New Covariance is a diagonal one
		//---------------------------------
		for (Row=0; Row<VecDim; Row++) {
			for (Col=0; Col<VecDim; Col++) {
				CMat[Row][Col] = 0.0;		
			}
			CMat[Row][Row] = EigVal[Row];		
		}

		MFree_1D(EigVal);
  }   // if (OrthT)

	MVec = MeanVec(pTData);     // overall mean vector

	if (SplitMethod) {

	  CurrMix = 1;
		//-------------------------------
		// Init first mixture components
		//-------------------------------
		for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
			MeanMat[0][DimCnt] = MVec[DimCnt];
			VariMat[0][DimCnt] = CMat[DimCnt][DimCnt];
		}
		WgtVec[0] = 1;

		//-------------------------------
		// generate GMM by spliting
		//-------------------------------
		while (CurrMix * 2 <= MixNum) {
			CurrMix = Split(CurrMix);
			if (Verbose) {
				cerr << "Current Mixture -->:" << CurrMix << endl;
			}
			reestimate(pTData->Mat, pTData->Row, CurrMix);
		}  // while
	} 
	else {
		GN_Rand REng;
		double Purt;

		REng.srandom(RandSeed);
		//-------------------------------
		// Init all mixture components
		//-------------------------------
		for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
			for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
				Purt = double (REng.random()) / double(REng.getmax() - 0.5);
				MeanMat[MixCnt][DimCnt] = MVec[DimCnt] + Purt*pTData->Mat[MixCnt][DimCnt];
				VariMat[MixCnt][DimCnt] = CMat[DimCnt][DimCnt];
			}
			WgtVec[MixCnt] = 1.0/MixNum;
		}
		reestimate(pTData->Mat, pTData->Row, MixNum);
	}
	
	MFree_2D(CMat);
	MFree_1D(MVec);
	delete(pTData);

  return 0;
}

//=====================================================
//  double mixture numbers by spliting
//  memory must large enough to split
//=====================================================
int SV_Model_GMM::Split(int CdNum) {

  int   CdCnt, DimCnt;        /* loop counter variables      */
  double Orig;

  for (CdCnt=0; CdCnt<CdNum; CdCnt++) {
     for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
		Orig = MeanMat[CdCnt][DimCnt];
		MeanMat[CdCnt + CdNum][DimCnt] = Orig + 0.01 * VariMat[CdCnt][DimCnt];
		MeanMat[CdCnt][DimCnt] = Orig - 0.01 * VariMat[CdCnt][DimCnt];
		VariMat[CdCnt + CdNum][DimCnt] = VariMat[CdCnt][DimCnt]; /* Dup */
     }

     /*----------------------*/
     /* weight as half       */
     /*----------------------*/
     WgtVec[CdCnt] /= 2.0;
     WgtVec[CdCnt + CdNum]= WgtVec[CdCnt];
  }

  return (CdNum * 2);    /* double size */
}


//============================================================
//  EM algorithm estimate MeanMat, VariMat, and WgtVec
//============================================================
void SV_Model_GMM::reestimate(float **DataMatrix, int DataNum, int CurrCode) {

   double ThisProb = -LARGE / 2.0, LastProb = -LARGE;

   double **NewBook, **NewVari, *NewWgt;
   double ProbSum, *CurrProb, TmpFloat;
   int DataCnt, CodeCnt, DimCnt;
   int Iter = 0, ThresNum = 0;

   MArray_2D(NewBook, CurrCode, VecDim, double, "NewBook");
   MArray_2D(NewVari, CurrCode, VecDim, double, "NewVari");
   MArray_1D(NewWgt,  CurrCode, double, "NewWgt");
   MArray_1D(CurrProb, CurrCode, double, "CurrProb");

   /*-----------------------------------*/
   /* repeat until no prob. increasing  */
   /*-----------------------------------*/
   while ( ThresNum <= 3 && Iter <= MaxIter ) {

     Iter++;
     /*---------------------------------------------------------*/
     /* count the iteration number coutinuously lower the thres */
     /*---------------------------------------------------------*/
     if ((ThisProb - LastProb) / fabs(ThisProb) < 0.001) ThresNum++;
     else ThresNum = 0;

    /*------------------------------------------*/
    /*  clear accumulating buffer               */
    /*------------------------------------------*/
    for (CodeCnt = 0; CodeCnt<CurrCode; CodeCnt++) {
      for (DimCnt = 0; DimCnt<VecDim; DimCnt++) {
		    NewBook[CodeCnt][DimCnt] = 0.0;
		    NewVari[CodeCnt][DimCnt] = 0.0;
      }
      NewWgt[CodeCnt] = 0.0;
    }

	  LastProb = ThisProb; ThisProb = 0.0;
    //------------------------------------------*/
    //  for each data vector 
    //------------------------------------------*/
    for (DataCnt = 0; DataCnt<DataNum; DataCnt++) {  /* one loop */

      /*--------------------------------------*/
      /*  calculate p(i|x_t, lamda)           */
      /*--------------------------------------*/
      ProbSum = 0.0;
      for (CodeCnt = 0; CodeCnt<CurrCode; CodeCnt++) {
		    CurrProb[CodeCnt] = WgtVec[CodeCnt] * Gaussian (MeanMat[CodeCnt], VariMat[CodeCnt], DataMatrix[DataCnt], VecDim);
		    ProbSum += CurrProb[CodeCnt];
      }
		
	  //-----------------------------------
	  // normalize by sum
	  //-----------------------------------
      for (CodeCnt = 0; CodeCnt<CurrCode; CodeCnt++) {
		    CurrProb[CodeCnt] /= ProbSum;
	    }
	  
      ThisProb += log(ProbSum);

      /*--------------------------------------*/
      /*  calculate p_i, u_i, delta_i         */
      /*--------------------------------------*/
      for (CodeCnt = 0; CodeCnt<CurrCode; CodeCnt++) {
		    for (DimCnt = 0; DimCnt<VecDim; DimCnt++) {
			    TmpFloat = CurrProb[CodeCnt] * DataMatrix[DataCnt][DimCnt];
			    NewBook[CodeCnt][DimCnt] += TmpFloat;
			    NewVari[CodeCnt][DimCnt] += TmpFloat * DataMatrix[DataCnt][DimCnt];
		    }
		    NewWgt[CodeCnt] += CurrProb[CodeCnt];    /* last column, p_i */
      }
    }   /* for (DataCnt) */

    if (Verbose) {
		cerr << ThisProb/DataNum << endl;
	}

    /*------------------------------------------*/
    /*  generate new mean and variance book     */
    /*------------------------------------------*/
    for (CodeCnt = 0; CodeCnt<CurrCode; CodeCnt++) {
      for (DimCnt = 0; DimCnt<VecDim; DimCnt++) {
	      TmpFloat = NewBook[CodeCnt][DimCnt] / NewWgt[CodeCnt];
	      MeanMat[CodeCnt][DimCnt] = TmpFloat;
	      VariMat[CodeCnt][DimCnt] = NewVari[CodeCnt][DimCnt] / NewWgt[CodeCnt] - TmpFloat * TmpFloat;
	      //------------------------------------------
	      // Variance floor should depends on overall 
	      // variance, not fixed value.
	      // 0.01 * CMat[DimCnt][DimCnt]
	      // 
	      //  F. Bimbot
	      //  http://www.PTT-Telecom.nl/cave
	      //------------------------------------------
	      if (VariMat[CodeCnt][DimCnt] < 0.01*CMat[DimCnt][DimCnt]) {
		      VariMat[CodeCnt][DimCnt] = 0.01*CMat[DimCnt][DimCnt];
        }
      }
      WgtVec[CodeCnt] = NewWgt[CodeCnt] / (double) DataNum;
    }

   }    /* while */

   MFree_2D(NewBook);
   MFree_2D(NewVari);
   MFree_1D(NewWgt);
   MFree_1D(CurrProb);
}


//--------------------------------------------------*/
//  Evaluated log-likelihood for a given vector     */
//--------------------------------------------------*/
double SV_Model_GMM::Log_Likelihood(float *TestVec, int Dim) {

	float *OneVec;
	double ProbSum;
	int Col, DimCnt, MixCnt;

	MArray_1D(OneVec, VecDim, float, "OneVec");
	/*--------------------------------------*/
	/* Transform or copy data to OneVec     */
	/*--------------------------------------*/
	if (OrthT)  {
		for (Col=0; Col<VecDim; Col++) {
			OneVec[Col] = 0.0;
			for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
				OneVec[Col] += float(TestVec[DimCnt] * OMat[DimCnt][Col]);
			}
		}
	}
	else {  
		for (Col=0; Col<VecDim; Col++) {
			OneVec[Col] = TestVec[Col];
		}	
	}

	
  /*--------------------------------------*/
  /*  calculate p(i|x_t, lamda)           */
  /*--------------------------------------*/
  ProbSum = 0.0;
  for (MixCnt = 0; MixCnt<MixNum; MixCnt++)
     ProbSum += WgtVec[MixCnt] * Gaussian(MeanMat[MixCnt], VariMat[MixCnt], OneVec, VecDim);


  MFree_1D(OneVec);
  return( log(ProbSum) );

};



//================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
// by thilo: changed all calls to write()
//================================================================
int SV_Model_GMM::SaveModel(void) {

	int RtCode;
	//int TotalByte; //by thilo
	int bytes; //by thilo
	SV_DataIO io; //by thilo

	Hdr.ModelType = MT_GMM;
    RtCode = SaveHdr();
	if (RtCode == SVLIB_Fail) {
		REPORT_ERROR(SVLIB_Fail, "Save GMM model Failed!");
	}
	
	//--------------------------
	// 
	//--------------------------
	//DFile.write((char*)(&VecDim), sizeof(int));
	bytes = io.writeScalar(&(this->DFile), this->VecDim);
	//DFile.write((char*)(&MixNum), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->MixNum);
	//DFile.write((char*)(&OrthT), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->OrthT);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save GMM Model Failed!");
	}
	//TotalByte = MHLen + 3*sizeof(int);


	//-------------------------------------
	// write mean and vari matrices
	//-------------------------------------
	//DFile.write((char*)(MeanMat[0]), MixNum*VecDim*sizeof(double));
	bytes += io.writeMatrix(&(this->DFile), this->MeanMat, this->MixNum, this->VecDim);
	//DFile.write((char*)(VariMat[0]), MixNum*VecDim*sizeof(double));
	bytes += io.writeMatrix(&(this->DFile), this->VariMat, this->MixNum, this->VecDim);
	//TotalByte += 2*MixNum*VecDim*sizeof(double);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save GMM Model Failed!");
	}


	//------------------------------------------------
	// write weight and orthogonal transform matrix
	//------------------------------------------------
	//DFile.write((char*)WgtVec, MixNum * sizeof(double));
	bytes += io.writeArray(&(this->DFile), this->WgtVec, this->MixNum);
	//DFile.write((char*)OMat[0], VecDim * VecDim *sizeof(double));
	bytes += io.writeMatrix(&(this->DFile), this->OMat, this->VecDim, this->VecDim);
	//TotalByte += MixNum * sizeof(double) + VecDim * VecDim *sizeof(double);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save GMM Model Failed!");
	}

	//return(TotalByte);
	return bytes + MHLen; //MHLen is maybe not correct...
}

//===========================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
// by thilo: changed all calls to read()
//===========================================================
SV_Model * SV_Model_GMM::LoadModel(void) {

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
	if (Hdr.ModelType != MT_GMM) {
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

	//DFile.read((char*)(&NewSize), sizeof(int));
	bytes += io.readScalar(&(this->DFile), NewSize, codeSizes, fileSizes);
	//DFile.read((char*)(&OrthT), sizeof(int));
	bytes += io.readScalar(&(this->DFile), this->OrthT, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewSize == 0) {
		return(NULL);
	}

	//--------------------------------------------
	// Load Mean Vec and Inv Covariance matrix
	//--------------------------------------------
	if (NewDim != VecDim || NewSize != MixNum) {

		MFree_2D(MeanMat);
		MFree_2D(VariMat);
		MFree_2D(OMat);
		MFree_1D(WgtVec);

		VecDim = NewDim;
		MixNum = NewSize;
		AllocModel();

	}; 

	//-------------------------------------
	// read mean and vari matrices
	//-------------------------------------
	//DFile.read((char*)(MeanMat[0]), MixNum*VecDim*sizeof(double));
	bytes += io.readMatrix(&(this->DFile), this->MeanMat, this->MixNum, this->VecDim, codeSizes, fileSizes);
	//DFile.read((char*)(VariMat[0]), MixNum*VecDim*sizeof(double));
	bytes += io.readMatrix(&(this->DFile), this->VariMat, this->MixNum, this->VecDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Load GMM Model Failed!");
	}

	//------------------------------------------------
	// read weight and orthogonal transform matrix
	//------------------------------------------------
	//DFile.read((char*)WgtVec, MixNum * sizeof(double));
	bytes += io.readArray(&(this->DFile), this->WgtVec, this->MixNum, codeSizes, fileSizes);
	//DFile.read((char*)OMat[0], VecDim * VecDim *sizeof(double));
	bytes += io.readMatrix(&(this->DFile), this->OMat, this->VecDim, this->VecDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Load GMM Model Failed!");
	}

	return(this);
}

//===============================================
// Test GMM, give average log likelihood
//===============================================
SV_Data* SV_Model_GMM::TestModel(SV_Data *pData) {

	SV_Data *Score, *DataCurr;
	int PatNum, PatCnt, VecCnt;
	double LL;

	if (pData->Col != VecDim) {
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
		
		LL = 0.0;
		for (VecCnt=0; VecCnt<DataCurr->Row; VecCnt++) {
			LL += Log_Likelihood(DataCurr->Mat[VecCnt], VecDim);
		}
		
		Score->Mat[PatCnt][0] = float(LL) / DataCurr->Row;
		DataCurr = DataCurr->Next;		   // point to next
	}

	return (Score);
}

//=============================================
// Dump model's parameter in ASCII 
//=============================================
ostream& operator<< (ostream& OutS, SV_Model_GMM& Data) {

	int Row, Col;

	OutS << Data.MixNum <<" "<<Data.VecDim << endl;
	for (Row=0; Row<Data.MixNum; Row++) {

		OutS << Data.WgtVec[Row] << endl;  // weight
		//-----------------------------
		// Mean vector
		//-----------------------------
		for (Col=0; Col<Data.VecDim; Col++)
			OutS<<Data.MeanMat[Row][Col]<< " ";
		OutS<<endl;

		//-----------------------------
		// Variance vector
		//-----------------------------
		for (Col=0; Col<Data.VecDim; Col++)
			OutS<<Data.VariMat[Row][Col]<< " ";
		OutS<<endl;

	}
	
	//--------------------------
	// Orthogonal transform
	//--------------------------
	if (Data.OrthT) {
		for (Row=0; Row<Data.VecDim; Row++) {
			for (Col=0; Col<Data.VecDim; Col++) {
				OutS<<Data.OMat[Row][Col]<< " ";		
			}
			OutS << "\n";
		}
	}


	return(OutS);
}



