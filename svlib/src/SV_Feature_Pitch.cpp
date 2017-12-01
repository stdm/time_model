//************************************************************************
//    
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "SV_Feature_Pitch.h"
#include "SV_Error.h"
#include "GN_Filter.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
#define BlankRate      0.4
#define Detectors	    8         /* pulse detector number */

float **PulseTrain;
float PitchMatrix[6][Detectors];
int NPulse;
float PitchAv;	

//----------------------------
// function prototype
//----------------------------
int Impulse_Trains(float* Data, int WinSz);
void Fill_PitchMatrix(void);		
float SelectMaxCoin (void);  // return pitch
void median_smooth(float *PitchArray, int Len);

//==========================================
// default constructor
//==========================================
SV_Feature_Pitch::SV_Feature_Pitch() : SV_Feature(){

	Para.WinSz		= 512;			// default 256 samples
	Para.StpSz		= 256;			// default 128 samples
	Para.RmvSilence = 0;			// 1: remove silence part, 0: intact
	Para.Smooth     = 1;			// 1: smooth pitch contour with 3-point median filter
	PulseTrain		= NULL;
	PitchAv			= float(8.0);   // initial guess for pitch (ms) 

}

//==========================================
// default desstructor
//==========================================
SV_Feature_Pitch::~SV_Feature_Pitch() {

	if (PulseTrain != NULL) {
		MFree_2D(PulseTrain);
	}

}

//-------------------------------------------------
// Calculate pitch contour for a given signal
//-------------------------------------------------
SV_Data *SV_Feature_Pitch::ExtractFeature(void) {

	SV_Data *DataSet = NULL;
	float	*SigBuf, *pSig;
	long	SigLen;


	//--------------------------------
	// Check if Signal is loaded
	//--------------------------------
	if (IsSigLoaded()) { 
		SigBuf = GetSig();
		SigLen = GetLen();
	}
	else return (NULL);

	if (Para.RmvSilence) {
		RemoveSilence(0.3, 32);
		SigLen = GetLen();  // new SigLen
	}

	//--------------------------------
	// Allocate memory for DataSet
	//--------------------------------
	DataSet = new SV_Data;
	if (DataSet==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}
	DataSet->Row = (SigLen - Para.WinSz + Para.StpSz) / Para.StpSz;
	DataSet->Col = 1;  // only feature 
	DataSet->Alloc();

	//------------------------------------
	// Allocate memory for PulseTrain
	//------------------------------------
	if(PulseTrain != NULL) {
		MFree_2D(PulseTrain);
	}
	MArray_2D(PulseTrain, 10, Para.WinSz, float, "PulseTrain");

	//-------------------------------------------------
	// bandpass filtering signal (100 - 900 Hz)
	// normalized to sampling rate to 2 Hz (Nyquest 1 Hz
	//--------------------------------------------------
	GN_Filter FIR;
	//------------------------------------------
	// WARNING: integer operation 100/Para.SRate*2 is ZERO (0)
	//------------------------------------------
	FIR.BandPass(100.0/Para.SRate*2.0, 900.0/Para.SRate*2.0);
	FIR.FIR_Filter(SigBuf, SigLen);  // filtering data

	//--------------------------------------------------
	// Extract pitch period from each frame  
	//--------------------------------------------------
	pSig = SigBuf;
	for (int FrmCnt=0; FrmCnt<DataSet->Row; FrmCnt++) {
		NPulse = Impulse_Trains(pSig, Para.WinSz);
		Fill_PitchMatrix();
		DataSet->Mat[FrmCnt][0] = SelectMaxCoin ();
		pSig += Para.StpSz;
	}

	//--------------------------------------------------
	// Smooth pitch contour  
	//--------------------------------------------------
	if (Para.Smooth) {
		median_smooth(DataSet->Mat[0], DataSet->Row);
	}

	return(DataSet);
}

/*-------------------------------------------------------------*/
/*  I make use of 8 pulse trains instead of 6                  */
/*  (Gold & Rabiner's paper                                    */
/*  M4 is a train of local maximum if M1, M8 is the local      */
/*  minimum train of M5 M1, M2, M3, M5, M6 and M7              */
/*  corresponding to original impluse trains                   */
/*-------------------------------------------------------------*/
int Impulse_Trains(float* Data, int WinSz) {

	int SampleCnt, PulseCnt, PositiveCnt=0, NegtiveCnt=0;
	float CurrentSample, TmpVal;

	/*---------------------------------------------------*/
	/* find positive and negtive pulse train 1 - 8       */
	/*---------------------------------------------------*/
	for(SampleCnt=1; SampleCnt < (WinSz-1); SampleCnt++) {

		CurrentSample = Data[ SampleCnt ];
		/*-----------------------------------------------------*/
		/* for a local maximum, generate positive pulse (0-3)  */
		/*-----------------------------------------------------*/
		if ((Data[ SampleCnt - 1 ] < CurrentSample) &&
			(Data[ SampleCnt + 1 ] <= CurrentSample )) {

			/*---------- Positive train 1 -----------------------*/
			/* for each local maximum, there is a pluse          */
			/*---------------------------------------------------*/
			PulseTrain[0][ PositiveCnt ] = CurrentSample;
			PulseTrain[8][ PositiveCnt ] = float(SampleCnt);

			/*---------- Positive train 2 -----------------------*/
			/*  difference this peak with last peak, if<0, set 0 */
			/*---------------------------------------------------*/
			if (PositiveCnt==0) {TmpVal = CurrentSample;}
			else {TmpVal = CurrentSample - PulseTrain[0][PositiveCnt-1];}

			if (TmpVal>0) {PulseTrain[1][ PositiveCnt ] = TmpVal;}
			else  {PulseTrain[1][ PositiveCnt ] = 0;}

			/*---------- Positive train 3 -----------------------*/
			/* if left and right peak less than this peak,       */
			/*---------------------------------------------------*/
			if (PositiveCnt >= 2) {
				TmpVal = PulseTrain[0][PositiveCnt - 1];
				if ((PulseTrain[0][PositiveCnt-2]<=TmpVal)&& (CurrentSample<=TmpVal)) {
					PulseTrain[2][ PositiveCnt - 1 ] = TmpVal;
				} else  {PulseTrain[2][ PositiveCnt - 1 ] = 0;}
			}

			/*---------- Positive train 4 -----------------------*/
			/* if left and right peak large than this peak,      */
			/*---------------------------------------------------*/
			if (PositiveCnt >= 2) {
				TmpVal = PulseTrain[0][PositiveCnt - 1];
				if ((PulseTrain[0][PositiveCnt-2]>=TmpVal) && (CurrentSample>=TmpVal)) {
					PulseTrain[3][ PositiveCnt - 1 ] = TmpVal;
				}  else {PulseTrain[3][ PositiveCnt - 1 ] = 0;}

			}
			PositiveCnt++;
		}   /* if positive */
//===============================================================
		/*-----------------------------------------------*/
		/* a local minimum, generate pulse train  5 - 8  */
		/*-----------------------------------------------*/
		if ((Data[ SampleCnt - 1 ] > CurrentSample) &&
			(Data[ SampleCnt + 1 ] >= CurrentSample)) {

			/*---------- Negtive train 5 -----------------------*/
			/* each minimum, there are is pulse                 */
			/*--------------------------------------------------*/
			PulseTrain[4][ NegtiveCnt ] = CurrentSample;
			PulseTrain[9][ NegtiveCnt ] = float(SampleCnt);  /* Pulse position */

			/*---------- Negtive train 6 -----------------------*/
			/* Difference of this and last pulses               */
			/*--------------------------------------------------*/
			if (NegtiveCnt==0) {TmpVal = CurrentSample;}
			else {TmpVal = CurrentSample - PulseTrain[4][NegtiveCnt-1];}

			if (TmpVal<0) {PulseTrain[6][ NegtiveCnt ] = TmpVal;}
			else  {PulseTrain[6][ NegtiveCnt ] = 0;}

			/*---------- Positive train 7 ----------------------*/
			/* both side large                                  */
			/*--------------------------------------------------*/
			if (NegtiveCnt >= 2) {
				TmpVal = PulseTrain[4][NegtiveCnt - 1];
				if ((PulseTrain[4][NegtiveCnt-2]>=TmpVal)&& (CurrentSample>=TmpVal)){
					PulseTrain[7][ NegtiveCnt - 1 ] = TmpVal;
				}  else {PulseTrain[7][ NegtiveCnt - 1 ] = 0;}
			}

			/*---------- Negtive train 8 -----------------------*/
			/* both side small                                  */
			/*--------------------------------------------------*/
			if (NegtiveCnt >= 2) {
				TmpVal = PulseTrain[4][NegtiveCnt - 1];
				if ((PulseTrain[4][NegtiveCnt-2]<=TmpVal)&&(CurrentSample<=TmpVal)) {
					PulseTrain[5][ NegtiveCnt - 1 ] = TmpVal;
				} else {PulseTrain[5][ NegtiveCnt - 1 ] = 0;}
			}
			NegtiveCnt++;
		}  /* if negtive */

	}  /* for SampleCnt*/


   /*-----------------------------------------------*/
   /* select minimum pulse number of +/-            */
   /*-----------------------------------------------*/
   if (PositiveCnt<=NegtiveCnt) PulseCnt = PositiveCnt;
   else  PulseCnt = NegtiveCnt;
	
   return (PulseCnt);
}

/*----------------------------------------------------------------*/
/* filling 6x8 pitch matrix                                       */
/* Line0 of PitchMatrix is the pitch estimation from current (t)  */
/* frame of speech data. Line 1, 2 are (t-1) and (t-2) pitch data */
/* Line3 = Line0 + Line1; Line4=Line1+Line2                       */
/* Line5 = Line0 + Line1 + Line2                                  */
/*----------------------------------------------------------------*/
void SV_Feature_Pitch::Fill_PitchMatrix(void) {


   for (int Cnt=0; Cnt<Detectors; Cnt++) {
      PitchMatrix[2][Cnt] = PitchMatrix[1][Cnt];
      PitchMatrix[1][Cnt] = PitchMatrix[0][Cnt];
      if (Cnt<4) {      /* positive trains */
		PitchMatrix[0][Cnt] = Pitch_Detector(PulseTrain[Cnt], PulseTrain[8], NPulse);
      }
		else {            /* negtive trains */
       PitchMatrix[0][Cnt] = Pitch_Detector(PulseTrain[Cnt],PulseTrain[9], NPulse);
		
	  }
      PitchMatrix[3][Cnt] = PitchMatrix[0][Cnt]+PitchMatrix[1][Cnt];
      PitchMatrix[4][Cnt] = PitchMatrix[1][Cnt]+PitchMatrix[2][Cnt];
      PitchMatrix[5][Cnt] = PitchMatrix[0][Cnt]+PitchMatrix[1][Cnt]+PitchMatrix[2][Cnt];
   }

}

/*-----------------------------------------------------*/
/*  Calculate average pitch from a pulse of train      */
/* It return the averaged period in samples per (ms)   */
/*-----------------------------------------------------*/
float SV_Feature_Pitch::Pitch_Detector (float *PulseVal, float *PulsePos, int PulseNum) {

	int Cnt;
	float GoingIndex;
	float Coeff = -0.003f;
	float LastPulsePos, LastPulseVal;
	float GoingVal, ConstantVal, BlankLen;


	LastPulseVal = (float)fabs(PulseVal[0]);
	LastPulsePos = PulsePos[0];
	BlankLen = float(BlankRate * PitchAv * Para.SRate/1000.0);
	
	for (Cnt=1; Cnt<PulseNum; Cnt++) {

		GoingIndex = PulsePos[Cnt] - LastPulsePos;
		if (GoingIndex > BlankLen) {     /* exceed blanking time */
			ConstantVal = PulsePos[Cnt] - LastPulsePos - BlankLen * PitchAv * Coeff;
			GoingVal = float(LastPulseVal * exp(ConstantVal));
			if (GoingVal<fabs(PulseVal[Cnt])) {      /* trigged next pulse */
				PitchAv = float((PitchAv + (PulsePos[Cnt]-LastPulsePos) / (Para.SRate/1000.0)) / 2.0);

				/*-------------------------------*/
				/* prevent extrame values        */
				/*-------------------------------*/
				if (PitchAv  > 15 ) PitchAv = 15.0;
				if (PitchAv  < 4)   PitchAv = 4.0;
				
				LastPulseVal = float(fabs(PulseVal[Cnt]));
				LastPulsePos = PulsePos[Cnt];
				BlankLen = float(BlankRate * PitchAv * Para.SRate/1000.0);

			} /* if GoingVal*/
		}  /* if GoingIndex */
	}    /* for Cnt */

   return (PitchAv);
}

/*-----------------------------------------------------*/
/*  Select Max coin from PitchMatrix and return Pitch  */
/*  and the sum of conincedences of all pitch detectors */
/*-----------------------------------------------------*/
float SelectMaxCoin (void) {

  int DetectorCnt, MaxValue, SumCoin;
  int Cnt, Row, Col, BiasCnt;
  int Bias[4] = {2, 3, 6, 8};
  int ThisCoin[4]={-2, -3, -6, -8};      /* set negtive bias */
  int MaxCoin[Detectors]={0, 0, 0, 0, 0, 0, 0, 0};

  float  Threhold[4];
  float  CurrPitch, CurrDiff;
  int SameMaxCnt;
  float SameMaxAvg;         /* average of pitch which have same Coin */

  MaxValue  = -100;         /* set a very samll coincidence */
  for (DetectorCnt=0; DetectorCnt<Detectors; DetectorCnt++) {
     CurrPitch = PitchMatrix[0][DetectorCnt];
     /*--------------------------------------------------*/
     /*  Set Threhold array in (ms) depend current pitch */
     /*--------------------------------------------------*/
     for (Cnt=0; Cnt<4; Cnt++)
       Threhold[Cnt]  = float((Cnt+1.0) / 10.0 * ceil (CurrPitch / 3.1));

     /*-------------------------*/
     /* Calculate Coincidence   */
     /*-------------------------*/
     for(Col=0; Col<Detectors; Col++)   /* each detector */
       for(Row=0; Row<6; Row++) {         /* each line of PitchMatrix */
	 CurrDiff = (float)fabs(CurrPitch - PitchMatrix[Row][Col]);
	 for (BiasCnt=0; BiasCnt<4; BiasCnt++)
	   if (CurrDiff < Threhold[BiasCnt] ) ThisCoin[BiasCnt]++;
       }

     /*---------------------------------------------------------*/
     /* select Max Coincidence from 4 bias and clear ThisCoin[] */
     /*---------------------------------------------------------*/
     MaxCoin[DetectorCnt] = ThisCoin[ 0 ];
     for(Cnt=0; Cnt<4; Cnt++) {
      if (ThisCoin[Cnt] > MaxCoin[DetectorCnt])
	       MaxCoin[DetectorCnt] = ThisCoin[Cnt];
      ThisCoin[Cnt] = -Bias[Cnt];
     };
     /*--------------------------------*/
     /* find maximum coincident values */
     /*--------------------------------*/
     if (MaxValue<MaxCoin[DetectorCnt])
	   MaxValue = MaxCoin[DetectorCnt];

  }  /* for DetectorCnt */

  /*----------------------------------------------------------*/
  /* Select pitch with Max Coincidence and sum all Coin       */
  /*----------------------------------------------------------*/
  SameMaxCnt = 0;       /* Cnt of detectors which have same Max output */
  SameMaxAvg = 0.0;
  SumCoin    = 0;
  for (Cnt=0; Cnt<Detectors; Cnt++){
     SumCoin += MaxCoin[Cnt];
     /*----------------------------------------------------*/
     /*  find all pitch which all have Max Coincidences    */
     /*----------------------------------------------------*/
     if (MaxValue==MaxCoin[Cnt]) {
	SameMaxAvg += PitchMatrix[0][Cnt];
	SameMaxCnt++;
     }
  }  /* for Cnt */

   return(SameMaxAvg/SameMaxCnt);
   
}   /* SelectMaxCoin */

/*-------------------------------------------------*/
/* smooth pitch contour by nonlinear median filter */
/* consider only three frames                      */
/*-------------------------------------------------*/
void median_smooth(float *PitchArray, int Len) {

	int  Cnt;
	float *Orig, AvgPitch, Ratio;
	
	MArray_1D(Orig, Len+2, float, "Orig");

	//---------------------------------
	// average pitch
	//---------------------------------
	AvgPitch = 0.0;
	for (Cnt=0; Cnt<Len; Cnt++) {
		AvgPitch += PitchArray[Cnt];
		Orig[Cnt+1] = PitchArray[Cnt];   /**/
	}
	AvgPitch /= Len;      /* Average pitch */

	/*-----------------------------------------*/
	/* try to eliminate pitch double and half  */
	/*-----------------------------------------*/
	for (Cnt=1; Cnt<=Len; Cnt++) {
		Ratio = Orig[Cnt+1] / AvgPitch;
		if (Ratio>=1.8) Orig[Cnt+1] /= 2;
		if (Ratio<=0.4) Orig[Cnt+1] *= 2;
	}

	for (Cnt=0; Cnt<Len; Cnt++) {
		if (Orig[Cnt] >= Orig[Cnt+1] && Orig[Cnt] >= Orig[Cnt+2])
			PitchArray[Cnt] = Orig[Cnt+1]>Orig[Cnt+2]?Orig[Cnt+1]:Orig[Cnt+2];

		if (Orig[Cnt+1] >= Orig[Cnt] && Orig[Cnt+1] >= Orig[Cnt+2])
			PitchArray[Cnt] = Orig[Cnt]>Orig[Cnt+2]?Orig[Cnt]:Orig[Cnt+2];

		if (Orig[Cnt+2] >= Orig[Cnt] && Orig[Cnt+2] >= Orig[Cnt+1])
			PitchArray[Cnt] = Orig[Cnt]>Orig[Cnt+1]?Orig[Cnt]:Orig[Cnt+1];

	}
	
	MFree_1D(Orig);
};


