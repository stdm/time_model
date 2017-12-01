//************************************************************************
//    Implementation of CHMM
//
//
//    Author  : Jialong HE
//    Date    : May 11, 1999
//    Modified: 29.04.2008 by thilo: added initial state distribution
//
//    TODO : when length of training vector sequence < StaNum, Problem
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "SV_Model_CHMM.h"
#include "SV_Error.h"
#include "GN_Matrix.h"
#include "GN_Rand.h"
#include "SV_DataIO.h" //by thilo: to achieve machine-independant file io

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model_CHMM::SV_Model_CHMM() {
	OrthT = 1;
	OMat  = NULL; 

	StaNum  = 0;
	MixNum  = 0; 
	VecDim  = 0; 
	Tran	= NULL; 
	Emit	= NULL; 

	VecNum  = 0; 
	LMean   = 0; 
	LVari   = 0; 
	PatNum  = 0;                         /* number of training patterns  */
	MaxLen  = 0;                         /* Max length of training strings */ 

	//--------------------------------------
	// Default paras for generate model
	//--------------------------------------
	Verbose  = 1;  // display messages
	Mixtures = 8;  // number of mixtures
	WithOrth = 1;  // with orthogonal transform
	MaxIter  = 100;  // with orthogonal transform
	RandSeed = 1999;

	CMat = NULL; 	                      // overall covariance matrix 
	MVec = NULL; 	                      // overall mean of training data 
	ConfMat = NULL;

	this->initialStateProbability = NULL; //by thilo
	this->isLeftRightModel = true; //by thilo
}

//==========================================
// default destructor
//==========================================
SV_Model_CHMM::~SV_Model_CHMM() {
	FreeModel();
}


//==============================================
// Allocate memory for OMat, Tran, Emit
// Depends on StaNum, MixNum, VecDim, 
//==============================================
void SV_Model_CHMM::AllocModel(void) {
	int StaCnt, MixCnt;

	MArray_2D(OMat, VecDim, VecDim, double, "OMat");		// orth trans
	MArray_2D(Tran, StaNum, StaNum, double, "CHMM Trans");  // transion mat
	MArray_1D(this->initialStateProbability, StaNum, double, "SV_Model_CHMM.AllocModel: initialStateProbability");  //by thilo

	//----------------------------------
	// Output probabilities
	//----------------------------------
	MArray_2D(Emit, StaNum, MixNum, MixType, "Emit");
	for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
		for (MixCnt=0; MixCnt<MixNum; MixCnt++)  {
			MArray_1D(Emit[StaCnt][MixCnt].mean, VecDim, double, "mean");
			MArray_1D(Emit[StaCnt][MixCnt].vari, VecDim, double, "vari");
		}
	}
	
	return; //by thilo: i like returns at the end of procedures... :-)
}	 

//==============================================
// Release all dynamicly allocate memory
//==============================================
void SV_Model_CHMM::FreeModel(void) {
	int StaCnt, MixCnt;

	if (OMat != NULL) {MFree_2D(OMat);}
	if (Tran != NULL) {MFree_2D(Tran);}
	if (CMat != NULL) {MFree_2D(CMat);}
	if (MVec != NULL) {MFree_1D(MVec);}
	if (this->initialStateProbability != NULL) {MFree_1D(this->initialStateProbability);} //by thilo
	if (Emit != NULL) { 
		for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
			for (MixCnt=0; MixCnt<MixNum; MixCnt++)  {
				MFree_1D(Emit[StaCnt][MixCnt].mean);
				MFree_1D(Emit[StaCnt][MixCnt].vari);
			}
		}
		MFree_2D(Emit);
	}

	StaNum = 0;
	MixNum = 0;
	VecDim = 0;

	return; //by thilo: i like returns at the end of procedures... :-)
}

//------------------------------------------
// Init Tran and Emit. 
// CMat and MVec must be ready
//------------------------------------------
void SV_Model_CHMM::InitModel(bool isLeftRightModel) { //by thilo: added parameter to drive initialization of initial state probabilities) {
	int Row, Col;
	int DimCnt, StaCnt, MixCnt;
	double Sum, Purt;

	GN_Rand REng;

	REng.srandom(RandSeed);
	//---------------------------------------*/
	// Init transition prob. to random value */
	// keep sum in each row to be 1.0        */
	//---------------------------------------*/
	for (Row=0; Row<StaNum; Row++) {
		//block by thilo
		if (isLeftRightModel == true) {
			this->initialStateProbability[Row] = (Row == 0) ? 1.0 : 0.0; //equation (46)
		} else {
			this->initialStateProbability[Row] = 1.0 / (double)(this->StaNum); //uniform start-state probability as suggested by rabiner's paper, bottom-right corner of page 273 in section V.C
		}
		//end by thilo
		Sum = 0.0;
		for (Col=0; Col<StaNum; Col++) {
			if (ConfMat[Row][Col] != 0) {  /* only none ZERO element */
				Tran[Row][Col] = (double) REng.random() + 1.0;
				Sum += Tran[Row][Col];
			}
			else {
				Tran[Row][Col] = 0.0;
			}
		}

		//----------------------------------------------------*/
		// Normalize sum to 1, to meet statistical costraint  */
		//----------------------------------------------------*/
		for (Col=0; Col<StaNum; Col++) {
				Tran[Row][Col] /= Sum;
		}
	}  // for (Row)

	//--------------------------------------------------------*/
	// Initalized mixture components in each state            */
	//--------------------------------------------------------*/
	for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
		for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
			Emit[StaCnt][MixCnt].wgt = 1.0 / (float)MixNum;
			for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
				Purt = double (REng.random()) / double(REng.getmax() - 0.5);
				Emit[StaCnt][MixCnt].vari[DimCnt] = CMat[DimCnt][DimCnt];
				Emit[StaCnt][MixCnt].mean[DimCnt] = MVec[DimCnt] + 0.1 * Purt * CMat[DimCnt][DimCnt];
			}
		}  /* for (MixCnt) */
	}

	return; //by thilo: i like returns at the end of procedures... :-)
}

/*------------------------------------------------------------------*/
/* This is the kernal program which does EM estimation for          */
/* parameters in CDHMM. Model structure can be any type, but all    */
/* Gaussian components have diagnal covariance matrices.            */
/* Formula number refer to Paper by Rabiner (IEEE proceeding)       */
/*                                                                  */
/* Note: some external variables are needed. They are defined in    */
/* main program.                                                    */
/*                                                                  */
/* Author : Jialong He                                              */
/* Date   : 23-01-96                                                */
/* Modify : 03-01-97; add orthogonal transform                      */
/*          29.04.2008 by thilo; added initial state probability    */
/*------------------------------------------------------------------*/

/*-----------------------------------------------------------------*/
/* Viterbi algorithm: formula (32a - 35, 104-105c)                 */
/* Note: the content of MixOut will be replaced by its log value   */
/*-----------------------------------------------------------------*/
double SV_Model_CHMM::viterbi(double **Tran, double **MixOut, int* Path, int Len) {
  int **Back;                 /* Back tracking matrix */
  double **LTran, **Delta;    /* Log version of Tran, Delta is Alpha */
  double Prob, MaxValue;
  int LenCnt, StaCnt, MaxInd, Ind;

  MArray_2D(Back,  StaNum, MaxLen,   int, "Back");
  MArray_2D(LTran, StaNum, StaNum, double, "LTran");
  MArray_2D(Delta, StaNum, MaxLen, double, "Delta");

  /*--------------------------------------------------------------*/
  /* Preprocessing: turn Tran, Emit to their log version          */
  /*--------------------------------------------------------------*/
  for (StaCnt=0; StaCnt<StaNum; StaCnt++)  {

    /*--------------------------------------------*/
    /* Log version Transition probability matrix  */
    /*--------------------------------------------*/
    for (Ind=0; Ind<StaNum; Ind++)
      if (Tran[StaCnt][Ind] <= 0.0 ) LTran[StaCnt][Ind] = -LARGE;
      else  LTran[StaCnt][Ind]  = log(Tran[StaCnt][Ind]);

    /*--------------------------------------------*/
    /* Log version emission probability            */
    /*--------------------------------------------*/
    for (LenCnt=0; LenCnt<Len; LenCnt++)
      if (MixOut[StaCnt][LenCnt] <= 0.0 ) MixOut[StaCnt][LenCnt] = -LARGE;
      else  MixOut[StaCnt][LenCnt] = log(MixOut[StaCnt][LenCnt]);
	}     /* for (StaCnt) */

  /*--------------------------------------------------------------*/
  /* Init Delta, formula (105a)                                   */
  /*--------------------------------------------------------------*/
  //for (Ind=0; Ind<StaNum; Ind++)  Delta[Ind][0] = -LARGE; //commented out by thilo
  //Delta[0][0] = MixOut[0][0]; //commented out by thilo
	//block by thilo:
	for (Ind = 0; Ind < this->StaNum; Ind++) {
		if (this->initialStateProbability[Ind] > 0.0) {
			Delta[Ind][0] = log(this->initialStateProbability[Ind]) + MixOut[Ind][0]; //equation (105a) (now including the initial state probabilities)
		} else {
			Delta[Ind][0] = -LARGE + MixOut[Ind][0];
		}
		Back[Ind][0] = 0; //equation (32b)
	}
	//end by thilo

  /*--------------------------------------------------------------*/
  /* Induction to get Delta for t=1..Len-1, formula (105b)        */
  /*--------------------------------------------------------------*/
  for (LenCnt=1; LenCnt<Len; LenCnt++) {
    for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
      MaxValue = -LARGE;
      MaxInd   = 0;
      /*----------------------------------*/
      /* find Max prob. for this state    */
      /*----------------------------------*/
      for (Ind=0; Ind<StaNum; Ind++) {
		  	Prob = Delta[Ind][LenCnt-1] + LTran[Ind][StaCnt];
				if (Prob>MaxValue) {
					MaxValue = Prob;
					MaxInd   = Ind;
				}
      }
      Delta[StaCnt][LenCnt] = MaxValue + MixOut[StaCnt][LenCnt];
      Back[StaCnt][LenCnt]  = MaxInd;
    }  /* for (StaCnt) */
  }  /* for (LenCnt) */

  /*----------------------------------------------------*/
  /* Find out the last max prob and Ind, formula (105c) */
  /*----------------------------------------------------*/
  MaxValue = Delta[0][Len-1];
  MaxInd   = 0;
  for (StaCnt=1; StaCnt<StaNum; StaCnt++)
    if  (Delta[StaCnt][Len-1] > MaxValue) {
      MaxValue = Delta[StaCnt][Len-1];
      MaxInd   = StaCnt;
    }

  /*-----------------------------------------------*/
  /* Backtracking the maximum path,  formula (35)  */
  /*-----------------------------------------------*/
  //Path[0] = 0;  /* left-to-right model, start from first state */ //by thilo: commented out, see below
  for (LenCnt=LenCnt-1; LenCnt>=0; LenCnt--) { //by thilo: was LenCnt>0, but we initilaized Back[i][0] properly and can use it now because we don't necessarily have a left-to-right model here
    Path[LenCnt] = MaxInd;
    MaxInd = Back[MaxInd][LenCnt];
  }

  MFree_2D(LTran);
  MFree_2D(Back);
  MFree_2D(Delta);

  return(MaxValue);
}

/*-----------------------------------------------------------------*/
/* Given vector sequence, calculate output from Gaussian mixture   */
/* (emission probability) for each state at different frame (t)    */
/*-----------------------------------------------------------------*/
void SV_Model_CHMM::emission(double **MixOut, MixType **Emit, float **Data, int SeqLen) {
  int SeqCnt, StaCnt, MixCnt;
  double Wgt, *Mean, *Vari;

  for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
    for (SeqCnt=0; SeqCnt<SeqLen; SeqCnt++) {
			MixOut[StaCnt][SeqCnt] = 0.0;    /* clear for mixture acc.*/
			for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
				Wgt  = Emit[StaCnt][MixCnt].wgt;   /* simplified notation */
				Mean = Emit[StaCnt][MixCnt].mean;
				Vari = Emit[StaCnt][MixCnt].vari;
				MixOut[StaCnt][SeqCnt] += Wgt * Gaussian(Mean, Vari, Data[SeqCnt], VecDim);
			}
		}  // for (SeqCnt)
	}  // for (StaCnt)

	return; //by thilo: i like returns at the end of procedures... :-)
}

/*-----------------------------------------------------------------*/
/* Forward variable (Alpha) and scaling factor (Scale) from given  */
/* transition prob. (Tran) and output prob. (MixOut). formula 19-21*/
/*-----------------------------------------------------------------*/
double SV_Model_CHMM::forward(double **Alpha, double *Scale, double **Tran, double **MixOut, int Len) {
  int LenCnt, StaCnt, Ind;
  double Prob, Sum;

  /*--------------------------------------------------------------*/
  /* Init Alpha and Scaling factor, here is left-to-right model   */
  /* may init to some other values. formula (19)                  */
  /*--------------------------------------------------------------*/
  //for (Ind=0; Ind<StaNum; Ind++) { //commented out by thilo
	//  Alpha[Ind][0] = 0.0; //commented out by thilo
  //} //commented out by thilo
  //Alpha[0][0] = 1;//commented out by thilo
	//block by thilo
	Scale[0] = 0.0;
	for (Ind = 0; Ind < this->StaNum; Ind++) {
		Alpha[Ind][0] = this->initialStateProbability[Ind] * MixOut[Ind][0]; //equation (19), this time using the initial state probabilities because we can't be sure we have a left-to-right model here!
		Scale[0] += Alpha[Ind][0]; //equation (91)
	}
	//end by thilo

  //Scale[0]    = MixOut[0][0]; //commented out by thilo; has been treated above in the context of not necessarily having a left-to-right model here
  if (Scale[0]<= 0.0) {Prob = -LARGE;}
  else {Prob = log(Scale[0]);}

  /*--------------------------------------------------------------*/
  /* Induction to get Alpha for t=1..Len-1    formula (20)        */
  /*--------------------------------------------------------------*/
  for (LenCnt=1; LenCnt<Len; LenCnt++) {
    Scale[LenCnt] = 0.0;
    for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
      Sum = 0.0;
      for (Ind=0; Ind<StaNum; Ind++) {
				Sum += Alpha[Ind][LenCnt-1] * Tran[Ind][StaCnt];
		  }
      Alpha[StaCnt][LenCnt] = Sum * MixOut[StaCnt][LenCnt];
      Scale[LenCnt] += Alpha[StaCnt][LenCnt];  /* scaling factor */
    }

    /*-----------------------------------------------*/
    /* Scale Alpha for this (t=LenCnt), formula (91) */
    /*-----------------------------------------------*/
    if (Scale[LenCnt] == 0) {
      cerr <<"Error: Scaling factor is ZERO\n";
		  REPORT_ERROR(SVLIB_DivBy0, "Scale");
    }

    for (Ind=0; Ind<StaNum; Ind++) {
      Alpha[Ind][LenCnt] /= Scale[LenCnt];
		}

    if (Scale[LenCnt] <= 0) {Prob += -LARGE;}
    else {Prob += log(Scale[LenCnt]);}  /* Accu. Prob. for this vector */
  }  /* for (LenCnt) */

  return(Prob);
}

/*-----------------------------------------------------------------*/
/* Calculate scaled backward varible (Beta) from  scaling (Scale)  */
/* transison prob (Tran) and output prob. (MixOut)                 */
/*-----------------------------------------------------------------*/
void SV_Model_CHMM::backward(double **Beta, double *Scale, double **Tran, double **MixOut, int Len) {
  int LenCnt, StaCnt, Ind;

  /*--------------------------------------------------------------*/
  /* Init Beta, left-to-right model, formula (24)                 */
  /*--------------------------------------------------------------*/
	//by thilo: we don't assume left-to-right model architecture here any more, but it seems we don't need to change anything...
  for (Ind=0; Ind<StaNum; Ind++) {
	  Beta[Ind][Len-1] = 1.0 / Scale[Len-1];
  }

  /*--------------------------------------------------------------*/
  /* Induction to get Beta for t=Len-2..0, formula (25)           */ //by thilo: comment was wrong, stated t=Len-2..1 previously
  /*--------------------------------------------------------------*/
  for (LenCnt=Len-1; LenCnt>0; LenCnt--)
    for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
      Beta[StaCnt][LenCnt-1] = 0.0;
      for (Ind=0; Ind<StaNum; Ind++)
				Beta[StaCnt][LenCnt-1] += Tran[StaCnt][Ind] * MixOut[Ind][LenCnt] * Beta[Ind][LenCnt];
      Beta[StaCnt][LenCnt-1] /= Scale[LenCnt-1];  /* scale Beta */
    }  /* for (StaCnt) */

	return; //by thilo: i like returns at the end of procedures... :-)
}

/*-------------------------------------------------*/
/* reestimation procedure for CDHMM parameters     */
/* using formula (52-54), (111),                   */
/*-------------------------------------------------*/
int  SV_Model_CHMM::TrainModel(SV_Data *pData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)
  MixType **NEmit;
  SV_Data *DataCurr;
  double **MixOut, **NTran, **Alpha, **Beta, *Scale; 
  double ThisProb = -LARGE/2.0, LastProb = -LARGE;
  double Sum, Ratio, RatioD, Gama, Tmp;
	double *newPi, newPiNorm; //by thilo: new initial state probability and normalization factor

  int StaCnt, MixCnt, DimCnt;
  int FrmNum, FrmCnt, Ind, StrNum, Iter=0, ThresNum = 0;
  double Wgt, *Mean, *Vari;

	//--------------------------------
	// Check valid of pData
	//--------------------------------
	if (pData == NULL) {
		REPORT_ERROR(SVLIB_BadData, "SV_Model_CHMM:: TrainModel"); 	
	}

	//---------------------------------------------
	// allocate memory for OMat, Tran, Emit
	// need StaNum, MixNum, VecDim
	//---------------------------------------------
	FreeModel();

	StaNum = NState;
	VecDim = pData->Col;
	MixNum = Mixtures;     /* mixture number in one state */
  OrthT  = WithOrth;
	if (StaNum<=0 || VecDim <=0 || MixNum <= 0) {
		REPORT_ERROR(SVLIB_BadArg, "Wrong StaNum or VecDim or MixNum!"); 	
	}

	AllocModel();
	CMat = CovMatrix(pData);     // overall covariance

	//--------------------------------
	// Orthogonal transform data
	//--------------------------------
	if (OrthT) {
		GN_Matrix MEng;
		double *EigVal;

		MArray_1D(EigVal, VecDim, double, "TrainModel::EigVal");
		//----------------------------------------
		// find orthogonal trans mat OMat and 
		// transform pTData
		//----------------------------------------
		MEng.Eigen(CMat, OMat, EigVal, VecDim);  
		OrthTrans(pData, OMat);
		
		//---------------------------------
		// New Covariance is a diagonal one
		//---------------------------------
		for (int Row=0; Row<VecDim; Row++) {
			for (int Col=0; Col<VecDim; Col++) {
				CMat[Row][Col] = 0.0;		
			}
			CMat[Row][Row] = EigVal[Row];		
		}

		MFree_1D(EigVal);
  }   // if (OrthT)

	MVec = MeanVec(pData);
	InitModel(this->isLeftRightModel); //by thilo: added parameter

	//--------------------------------------------------------------
	// find out max length (MaxLen) and number of pattern (PatNum)
	//--------------------------------------------------------------
	DataCurr = pData;
	MaxLen = -1;
	PatNum = 0;
	while (DataCurr != NULL) {
		if (VecDim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "SV_Model_CHMM:: TrainModel"); 	
		}

		if (DataCurr->Row < this->StaNum) { //by thilo: warn if the problem described by jialong he in the trailing header-TODO is detected
			REPORT_ERROR(SVLIB_BadData, "C-HMM needs not less training vectors than states!"); //by thilo
		} //by thilo
		
		if (DataCurr->Row > MaxLen) {
			MaxLen = DataCurr->Row;
		}

		PatNum++;
		DataCurr = DataCurr->Next;
	}  // while

	//------------------------------------------------------------------*/
	// Allocate some buffers:                                           */
	//   Alpha : forward variable                                       */
	//   Beta  : backward variable                                      */
	//   MixOut: emission prob of one sequence                          */
	//   Scale : scaling factor                                         */
	//   NTran : new estimated transison prob.                         */
	//   NEmit : new Gaussian mixture model                            */
	//------------------------------------------------------------------*/
	MArray_2D(MixOut, StaNum, MaxLen, double, "MixOut");
	MArray_2D(Alpha,  StaNum, MaxLen, double, "Alpha");
	MArray_2D(Beta,   StaNum, MaxLen, double, "Beta");
	MArray_2D(NTran,  StaNum, StaNum, double, "NTran");
	MArray_2D(NEmit,  StaNum, MixNum, MixType, "NEmit");
	MArray_1D(Scale,  MaxLen, double, "Scale");
	MArray_1D(newPi,  this->StaNum, double, "SV_Model_CHMM.TrainModel: newPi"); //by thilo: new initial state probability

	for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
		for (MixCnt=0; MixCnt<MixNum; MixCnt++)  {
			MArray_1D(NEmit[StaCnt][MixCnt].mean, VecDim, double, "mean");
			MArray_1D(NEmit[StaCnt][MixCnt].vari, VecDim, double, "vari");
		}
	}

	//---------------------------------------------------------*/
	// repeat until no prob. increasing, or enough iteration   */
	//---------------------------------------------------------*/
	while ( ThresNum <= 3 && Iter <= MaxIter ) {

		/*---------------------------------------------------------*/
		/* count the iteration number coutinuously lower the thres */
		/*---------------------------------------------------------*/
		if ((ThisProb - LastProb) / fabs(ThisProb) < 0.001) {
			ThresNum++;
		}
		else {ThresNum = 0;}

		StrNum   = 0;
		LastProb = ThisProb;
		ThisProb = 0.0;
		Iter++;

		/*------------------------------------------------------*/
		/* Clear NTran, NEmit buffer for accumulating           */
		/*------------------------------------------------------*/
		for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
			for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
				NEmit[StaCnt][MixCnt].wgt = 0.0;
				for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
					NEmit[StaCnt][MixCnt].mean[DimCnt] = 0.0;
					NEmit[StaCnt][MixCnt].vari[DimCnt] = 0.0;
				}
			}
			for (Ind=0; Ind<StaNum; Ind++) {
				NTran[StaCnt][Ind] = 0.0;
			}
			newPi[StaCnt] = 0.0; //by thilo
		}
		newPiNorm = 0.0; //by thilo

		/*--------------------------------------------------------------*/
		/* One Epoch, accumulating for all strings belong to this class */
		/*--------------------------------------------------------------*/
		DataCurr = pData;
		while (DataCurr != NULL) {

			StrNum++;     /* count string of this class */
			FrmNum = DataCurr->Row;
			emission(MixOut, Emit, DataCurr->Mat, FrmNum);
			ThisProb += forward(Alpha, Scale, Tran, MixOut, FrmNum);
			backward(Beta, Scale, Tran, MixOut, FrmNum);

			//block by thilo: reestimate the initial state probability pi by equation (40a) and (27) 
			for (Ind = 0; Ind < this->StaNum; Ind++) {
				if (this->isLeftRightModel == false) {
					newPi[Ind] += Alpha[Ind][0] * Beta[Ind][0]; //this is gamma_t(j) from equation (27), which equals pi_i (equation (40a)), just without the normalization which follows later on
					newPiNorm += newPi[Ind];
				}
			}
			//end by thilo

			/*-----------------------------------------*/
			/* reestimation of Transition prob. NTran  */
			/*-----------------------------------------*/
			for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
				for (Ind=0; Ind<StaNum; Ind++) {
					if (Tran[StaCnt][Ind] == 0.0) continue;  /* estimate nonzero element*/
					Sum = 0.0;
					for (FrmCnt=0; FrmCnt<FrmNum-1; FrmCnt++) {
						Sum += Alpha[StaCnt][FrmCnt] * MixOut[Ind][FrmCnt+1] * Beta[Ind][FrmCnt+1];
					}
					NTran[StaCnt][Ind] += Sum * Tran[StaCnt][Ind];
				}
			}

			/*------------------------------------------*/
			/* reestimation Gaussian mixture parameters */
			/*------------------------------------------*/
			for (FrmCnt=0; FrmCnt<FrmNum; FrmCnt++) {
				/*--------------------------------------------*/
				/*  Denominator of Ratio by forward/backward  */
				/*--------------------------------------------*/
				RatioD = 0.0;
				for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
					RatioD += Alpha[StaCnt][FrmCnt] * Beta[StaCnt][FrmCnt];
				}

				if (RatioD == 0.0)  {
					REPORT_ERROR(SVLIB_DivBy0, "RatioD");
				}
				
				for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
					/*--------------------------------------------*/
					/*  Ratio by forward/backward variables       */
					/*--------------------------------------------*/
					Ratio = Alpha[StaCnt][FrmCnt] * Beta[StaCnt][FrmCnt] / RatioD;
					for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
						Wgt  = Emit[StaCnt][MixCnt].wgt;
						Mean = Emit[StaCnt][MixCnt].mean;
						Vari = Emit[StaCnt][MixCnt].vari;
						Gama = Gaussian (Mean, Vari, DataCurr->Mat[FrmCnt], VecDim) * Wgt * Ratio / MixOut[StaCnt][FrmCnt];
						
						/*---------------------------------------------*/
						/* accumulating Wgt, Mean, Vari for this frame */
						/* Numerator of formula (52, 53, 54)           */
						/*             2      2            2           */
						/* using E(x-u)  = E(X) -2uE(X) + u            */
						/*---------------------------------------------*/
						NEmit[StaCnt][MixCnt].wgt += Gama;   /* weight */
						for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
							Tmp = Gama * DataCurr->Mat[FrmCnt][DimCnt];
							NEmit[StaCnt][MixCnt].mean[DimCnt] += Tmp;
							NEmit[StaCnt][MixCnt].vari[DimCnt] += Tmp * DataCurr->Mat[FrmCnt][DimCnt];
						} /* for (DimCnt) */
					}   /* for (MixCnt) */
				}     /* for (StaCnt) */
			}       /* for (FrmCnt) */

			DataCurr = DataCurr->Next;
		}         /* while (DataCurr!=NULL) */
		
		/*---------------------------------------------------------*/
		/* Update model's transition prob. Tran, from NTran.       */
		/* Each column normalized its SUM to 1.0, formula (111)    */
		/*---------------------------------------------------------*/
		for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
			Sum = 0.0;
			for (Ind=0; Ind<StaNum; Ind++) {
				Sum += NTran[StaCnt][Ind];
			}
			for (Ind=0; Ind<StaNum; Ind++) {
				Tran[StaCnt][Ind] = NTran[StaCnt][Ind] / Sum;
			}
			if (this->isLeftRightModel == false) { //by thilo
				initialStateProbability[StaCnt] = newPi[StaCnt] / newPiNorm; //by thilo
			} //by thilo
		}
		
		/*---------------------------------------------------------*/
		/* Update model's Gaussian mixture parameters Emit from    */
		/* NEmit with proper normalization  (52-54)                */
		/*---------------------------------------------------------*/
		for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
			/*------------------------------*/
			/* Normalize weight sum to 1.0  */
			/*------------------------------*/
			Sum = 0.0;
			for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
				Sum += NEmit[StaCnt][MixCnt].wgt;
			}

			for (MixCnt=0; MixCnt<MixNum; MixCnt++) {
				for (DimCnt=0; DimCnt<VecDim; DimCnt++) {
					Tmp = NEmit[StaCnt][MixCnt].mean[DimCnt] / NEmit[StaCnt][MixCnt].wgt;
					Emit[StaCnt][MixCnt].mean[DimCnt] = Tmp;
					Emit[StaCnt][MixCnt].vari[DimCnt] = NEmit[StaCnt][MixCnt].vari[DimCnt] / NEmit[StaCnt][MixCnt].wgt - Tmp * Tmp;

					/*-----------------------------*/
					/* prevent too small variance  */
					/*-----------------------------*/
					if (Emit[StaCnt][MixCnt].vari[DimCnt] <= 0.0001)
						Emit[StaCnt][MixCnt].vari[DimCnt] = 0.0001;
				}
				Emit[StaCnt][MixCnt].wgt = NEmit[StaCnt][MixCnt].wgt / Sum;
				
				/*-----------------------------------------------------*/
				/* prevent outlier component causing numerical problem */
				/*-----------------------------------------------------*/
				if (Emit[StaCnt][MixCnt].wgt<0.0001) {
					Emit[StaCnt][MixCnt].wgt = 0.0001;
				}
			}   /* for (MixCnt) */
		}     /* for (StaCnt) */

		ThisProb /= (double)StrNum;
		if (Verbose) {
			cerr << ThisProb << endl;
		}

	}   /* while */

	/*------------------------------*/
	/* free allocated memory        */
	/*------------------------------*/
	for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
		for (MixCnt=0; MixCnt<MixNum; MixCnt++)  {
			MFree_1D(NEmit[StaCnt][MixCnt].mean);
			MFree_1D(NEmit[StaCnt][MixCnt].vari);
		}
	}

	MFree_2D(NEmit);
	MFree_2D(NTran);
	MFree_2D(MixOut);
	MFree_2D(Alpha);
	MFree_2D(Beta);
	MFree_1D(Scale);
	MFree_1D(newPi); //by thilo

  return 0;
}

//===============================================
// Test CHMM, give average log likelihood
//===============================================
SV_Data* SV_Model_CHMM::TestModel(SV_Data *pData) {
	double **MixOut, **Alpha,  *Scale; 
	SV_Data *Score, *DataCurr, *pOnePat;
	int PatNum, PatCnt, MaxLen;
	int Row, Col;
	double LL;


	if (pData->Col != VecDim) {
		REPORT_ERROR(SVLIB_BadData, "Model and Data have different Dim!");
	}

	//--------------------------------------------------------------
	// find out max length (MaxLen) and number of pattern (PatNum)
	//--------------------------------------------------------------
	DataCurr = pData;
	MaxLen = -1;
	PatNum = 0;
	while (DataCurr != NULL) {
		if (VecDim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "SV_Model_CHMM:: TestModel"); 	
		}
		
		if (DataCurr->Row > MaxLen) {
			MaxLen = DataCurr->Row;
		}

		PatNum++;
		DataCurr = DataCurr->Next;
	}  // while

	//----------------------------
	// allocate score space
	//----------------------------
	Score = new SV_Data;
	Score->Row = PatNum;	
	Score->Col = 1;	         // one score for each pattern
	Score->Alloc();


	MArray_2D(MixOut,  StaNum, MaxLen, double, "MixOut");
	MArray_2D(Alpha,   StaNum, MaxLen, double, "Alpha");
	MArray_1D(Scale,   MaxLen, double, "Scale");
	
	//--------------------------------------------------------
	// For orthogonal transform, keep original data unchanged
	//--------------------------------------------------------
	pOnePat = new SV_Data;
	pOnePat->Row = MaxLen;    // long enough to hold any pattern
	pOnePat->Col = VecDim;
	pOnePat->Alloc();

	//----------------------------
	// test score for each pattern
	//----------------------------
	DataCurr = pData;
	for (PatCnt=0; PatCnt<PatNum; PatCnt++) {
		
		//------------------------------------------
		// make a copy for orthogonal transfom 
		//------------------------------------------
		for (Row=0; Row<DataCurr->Row; Row++) {
			for (Col=0; Col<DataCurr->Col; Col++) {
				pOnePat->Mat[Row][Col] = DataCurr->Mat[Row][Col];
			}
		}
		pOnePat->Next = NULL;  // just one pattern

		if (OrthT) {
			OrthTrans(pOnePat, OMat); 
		}
	
		emission(MixOut, Emit, pOnePat->Mat, DataCurr->Row);
		LL = forward(Alpha, Scale, Tran, MixOut, DataCurr->Row);	
//		LL = viterbi(Tran, MixOut, CurrPath, DataCurr->Row);

		Score->Mat[PatCnt][0] = float(LL) / DataCurr->Row;
		DataCurr = DataCurr->Next;		   // point to next
	}

	MFree_2D(MixOut);
	MFree_2D(Alpha);
	MFree_1D(Scale);
	delete (pOnePat);

	return (Score);
}

//================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
// by thilo: changed all .write-calls
//================================================================
int SV_Model_CHMM::SaveModel(void) {
	int RtCode, StaCnt, MixCnt;
	//int TotalByte; //by thilo
	int bytes; //by thilo
	SV_DataIO io; //by thilo

	Hdr.ModelType = MT_CHMM;
    RtCode = SaveHdr();
	if (RtCode == SVLIB_Fail) {
		REPORT_ERROR(SVLIB_Fail, "Save C-HMM model Failed!");
	}
	
	//--------------------------
	//  save simple variables
	//--------------------------
	//DFile.write((char*)(&StaNum), sizeof(int));
	bytes = io.writeScalar(&(this->DFile), this->StaNum);
	//DFile.write((char*)(&VecDim), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->VecDim);
	//DFile.write((char*)(&MixNum), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->MixNum);
	//DFile.write((char*)(&OrthT), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->OrthT);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save C-HMM Model Failed!");
	}
	//TotalByte = MHLen + 4*sizeof(int);

	//------------------------------------------------
	// write tran. and orthogonal transform matrix
	//------------------------------------------------
	bytes += io.writeArray(&(this->DFile), this->initialStateProbability, this->StaNum); //new by thilo
	//DFile.write((char*)Tran[0], StaNum * StaNum * sizeof(double));
	bytes += io.writeMatrix(&(this->DFile), this->Tran, this->StaNum, this->StaNum);
	//DFile.write((char*)OMat[0], VecDim * VecDim *sizeof(double));
	bytes += io.writeMatrix(&(this->DFile), this->OMat, this->VecDim, this->VecDim);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save C-HMM Model Failed!");
	}
	//TotalByte += (MixNum * MixNum + VecDim * VecDim) * sizeof(double);

	//------------------------------------------------
	// Save emission prob.
	//------------------------------------------------
	for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
		for (MixCnt=0; MixCnt<MixNum; MixCnt++)  {
			//DFile.write((char*)&(Emit[StaCnt][MixCnt].wgt), sizeof(double));
			bytes += io.writeScalar(&(this->DFile), Emit[StaCnt][MixCnt].wgt);
			//DFile.write((char*)(Emit[StaCnt][MixCnt].mean), VecDim * sizeof(double));
			bytes += io.writeArray(&(this->DFile), Emit[StaCnt][MixCnt].mean, VecDim);
			//DFile.write((char*)(Emit[StaCnt][MixCnt].vari), VecDim * sizeof(double));
			bytes += io.writeArray(&(this->DFile), Emit[StaCnt][MixCnt].vari, VecDim);
		}
	}

	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save C-HMM Model Failed!");
	}
	//TotalByte += (1 + VecDim * 2) * MixNum * StaNum * sizeof(double);

	//return(TotalByte);
	return bytes + MHLen; //maybe not correct size because of MHLen instead of real written bytes for header
}

//===========================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
// by thilo: changed all read()-calls
//===========================================================
SV_Model * SV_Model_CHMM::LoadModel(void) {
	int RtCode, NewDim, NewSize, NewSta;
	int StaCnt, MixCnt;
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
	if (Hdr.ModelType != MT_CHMM) {
		return(NULL);
	}

	//--------------------------
	// CHMM model's parameters
	//--------------------------
	//DFile.read((char*)(&NewSta), sizeof(int));
	bytes = io.readScalar(&(this->DFile), NewSta, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewSta <= 0) {
		return(NULL);
	}

	//DFile.read((char*)(&NewDim), sizeof(int));
	bytes = io.readScalar(&(this->DFile), NewDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewDim <= 0) {
		return(NULL);
	}

	//DFile.read((char*)(&NewSize), sizeof(int));
	bytes = io.readScalar(&(this->DFile), NewSize, codeSizes, fileSizes);
	//DFile.read((char*)(&OrthT), sizeof(int));
	bytes = io.readScalar(&(this->DFile), this->OrthT, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewSize <= 0) {
		return(NULL);
	}



	//--------------------------------------------
	// Allocate model's space 
	//--------------------------------------------
	if (NewSta != StaNum || NewDim != VecDim || NewSize != MixNum) {

		FreeModel();
		StaNum = NewSta;
		VecDim = NewDim;
		MixNum = NewSize;
		AllocModel();
	}

	//------------------------------------------------
	// read tran. and orthogonal transform matrix
	//------------------------------------------------
	bytes = io.readArray(&(this->DFile), this->initialStateProbability, this->StaNum, codeSizes, fileSizes); //new by thilo
	//DFile.read((char*)Tran[0], StaNum * StaNum * sizeof(double));
	bytes = io.readMatrix(&(this->DFile), this->Tran, this->StaNum, this->StaNum, codeSizes, fileSizes);
	//DFile.read((char*)OMat[0], VecDim * VecDim *sizeof(double));
	bytes = io.readMatrix(&(this->DFile), this->OMat, this->VecDim, this->VecDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Load C-HMM Model Failed!");
	}

	//------------------------------------------------
	// Load emission prob.
	//------------------------------------------------
	for (StaCnt=0; StaCnt<StaNum; StaCnt++) {
		for (MixCnt=0; MixCnt<MixNum; MixCnt++)  {
			//DFile.read((char*)&(Emit[StaCnt][MixCnt].wgt), sizeof(double));
			bytes = io.readScalar(&(this->DFile), this->Emit[StaCnt][MixCnt].wgt, codeSizes, fileSizes);
			//DFile.read((char*)(Emit[StaCnt][MixCnt].mean), VecDim * sizeof(double));
			bytes = io.readArray(&(this->DFile), this->Emit[StaCnt][MixCnt].mean, VecDim, codeSizes, fileSizes);
			//DFile.read((char*)(Emit[StaCnt][MixCnt].vari), VecDim * sizeof(double));
			bytes = io.readArray(&(this->DFile), this->Emit[StaCnt][MixCnt].vari, VecDim, codeSizes, fileSizes);
		}
	}

	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Load C-HMM Model Failed!");
	}

	return(this);
}

//=============================================
// Dump model's parameter in ASCII 
//=============================================
ostream& operator<< (ostream& OutS, SV_Model_CHMM& Data) {

	int StaCnt, MixCnt, DimCnt, Row, Col;

	OutS << "#(States) #(Mixtures per state) #(Dimension)" << endl; //by thilo
	OutS << Data.StaNum << " " <<Data.MixNum <<" "<<Data.VecDim << endl;
	for (Row=0; Row<Data.StaNum; Row++) {
		for (Col=0; Col<Data.StaNum; Col++) {
			OutS<<Data.Tran[Row][Col]<< " ";
		}
		OutS<<endl;
	}


	for (StaCnt=0; StaCnt<Data.StaNum; StaCnt++) {
		OutS << "State: " << StaCnt << ", initial probability: " << Data.initialStateProbability[StaCnt] << endl; //by thilo
		for (MixCnt=0; MixCnt<Data.MixNum; MixCnt++) {
			OutS << " Mixture: " << MixCnt << ", weight: "; //by thilo
			OutS << Data.Emit[StaCnt][MixCnt].wgt << endl;
			OutS << " Mean: "; //by thilo
			for (DimCnt=0; DimCnt<Data.VecDim; DimCnt++) {
				OutS << Data.Emit[StaCnt][MixCnt].mean[DimCnt] << " ";
			}
			OutS << endl;
			OutS << " Variance: "; //by thilo
			for (DimCnt=0; DimCnt<Data.VecDim; DimCnt++) {
				OutS << Data.Emit[StaCnt][MixCnt].vari[DimCnt] << " ";
			}
			OutS << "\n\n";
		}

	}
	
	return(OutS);
}
