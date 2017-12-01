//    Implementation for DTW comparision
//
//
//    Author  : Jialong HE
//    Date    : May 21, 1999
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "SV_Model_DTW.h"
#include "SV_Error.h"
#include "vector.h"


static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model_DTW::SV_Model_DTW() {
	RefData = NULL;
	NumRef  = 0;
}

//==========================================
// default destructor
//==========================================
SV_Model_DTW::~SV_Model_DTW() {


}


//=======================================================
// DTW is a non-parametric model, it do not
// need training, just store training data
//
// NOTE: not copy, just store a point to training data 
//======================================================
int SV_Model_DTW::TrainModel(SV_Data *TrainData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	
	NumRef  = CountPat(TrainData);
	RefData = TrainData;

  return 0;
};

/*==============================================================*/
/* This procedure implements the dynamic time wraping algorithm */
/* to compare two patten sequences. The maximum slope is 2 and  */
/* the minimum slope is 1/2. Local continue constraints are:    */
/*                                                              */
/*     |  D(x-2, y-1) + [d(x-1, y) + d(x,y)] / 2		*/
/* min <  D(x-1, y-1) + d(x,y)                                  */
/*     |  D(x-1, y-2) + [d(x, y-1) + d(x,y)] / 2                */
/*                                                              */
/*                                                              */
/* Ref: L.R. Rabiner, B.H.Juang, "Fundamentals of speech        */
/*      recognition," New Jersey: Prentice Hall, 1993           */
/*--------------------------------------------------------------*/
/*  Author : Jialong He                                         */
/*  Date   : April 19, 1995                                     */
/*                                                              */
/*==============================================================*/
/*-------------------------------------------------------------*/
/* TestFrame[][], TestNum, test vector buffer and number       */
/* RefFrame[][], RefNum, reference vector buffer and number    */
/* Dim: vector dimension                                       */
/*-------------------------------------------------------------*/

double SV_Model_DTW::dtw_dist(float **TestFrame, int TestNum, float **RefFrame, int RefNum, int Dim) {

  float **GDist, *LDist, *Weight;
  float TopPath, MidPath, BotPath;
  float LocalDist, LocalDist_1;
  float Ratio;
	
  double DTW_Dist;

  int TestCnt, RefCnt;

  /*------------------------------------------*/
  /* The maximum (minimum) slope is 2 (1/2)   */
  /*------------------------------------------*/
  if ((( (float)TestNum / (float)RefNum) > 2.0) ||
     (( (float)TestNum / (float)RefNum) < 0.5))
    return (LARGE);   /* can't wraping, return a large distance */

  MArray_2D(GDist, TestNum, RefNum, float, "GDist");
  MArray_1D(LDist, RefNum, float, "LDist");
  MArray_1D(Weight, TestNum, float, "Weight");

  /*------------------------------------------*/
  /* Initialize all element in GDist matrix   */
  /* NOTE: C++ new NOT init to ZERO           */	
  /*------------------------------------------*/
  for(RefCnt = 0; RefCnt < RefNum; RefCnt++) {
	LDist[RefCnt] = float(0.0);
  }

  for(TestCnt = 0; TestCnt < TestNum; TestCnt++) {
	  for(RefCnt = 0; RefCnt < RefNum; RefCnt++) {
		GDist[TestCnt][RefCnt] = LARGE;
      }
    Weight[TestCnt] = (float)exp(-(double)TestCnt/TestNum);
  }
  /*------------------------------------------*/
  /* possible start point                     */
  /*------------------------------------------*/
  for(TestCnt = 0; TestCnt < 5; TestCnt++) {
     DIST_VEC(TestFrame[TestCnt], RefFrame[0], GDist[TestCnt][0], Dim, float, float);
     DIST_VEC(TestFrame[TestCnt], RefFrame[1], GDist[TestCnt][1], Dim, float, float);
  }

  for(RefCnt = 0; RefCnt < 5; RefCnt++) {
     DIST_VEC(TestFrame[0], RefFrame[RefCnt], GDist[0][RefCnt], Dim, float, float);
     DIST_VEC(TestFrame[1], RefFrame[RefCnt], GDist[1][RefCnt], Dim, float, float);
  }

  /*------------------------------------------------------*/
  /* check all test frames against all reference frames   */
  /*------------------------------------------------------*/
  for(TestCnt = 2; TestCnt < TestNum; TestCnt++) {

   DIST_VEC(TestFrame[TestCnt], RefFrame[1], LocalDist_1, Dim, float, float);
   for (RefCnt = 2; RefCnt < RefNum; RefCnt++) {

    /*-------------------------------------------------*/
    /* apply global path constraints, within two line  */
    /* y = k*x +/- T0 *RefNum  , k=0.3                 */
    /*-------------------------------------------------*/
    Ratio = (float)fabs (RefCnt - (float)RefNum / (float) TestNum * TestCnt);

    if (Ratio < 0.3 * RefNum)  {

      DIST_VEC(TestFrame[TestCnt], RefFrame[RefCnt], LocalDist, Dim, float, float);
/*      LocalDist *= Weight [TestCnt];   */
      TopPath = GDist[TestCnt-2][RefCnt-1] + (LocalDist + LDist[RefCnt]) / float(2.0);
      MidPath = GDist[TestCnt-1][RefCnt-1] + LocalDist;
      BotPath = GDist[TestCnt-1][RefCnt-2] + (LocalDist + LocalDist_1) / float(2.0);
      LocalDist_1 = LocalDist;
      LDist[RefCnt] = LocalDist;

      /*------------------------------------------*/
      /* find minimum among three path            */
      /*------------------------------------------*/
      if (TopPath < MidPath && TopPath < BotPath)
	GDist[TestCnt][RefCnt] = TopPath;
      else if (MidPath < BotPath)
	GDist[TestCnt][RefCnt] = MidPath;
      else
	GDist[TestCnt][RefCnt] = BotPath;

    }  /* if (Ratio) */
   }   /* for (RefCnt) */
  }    /* for (TestCnt) */

  DTW_Dist = GDist[TestNum - 1 ][RefNum - 1] / TestNum; 	

  MFree_1D(Weight);
  MFree_1D(LDist);
  MFree_2D(GDist);
  return(DTW_Dist);  /* Normalized distance */
}


//===============================================
// return distance between two vector sequence.
// The distance is obtained by DTW 
//===============================================
double SV_Model_DTW::DTW_Comp(SV_Data &First, SV_Data& Second) {

	double DTW_dist;	
	//-----------------------------
	// Simply check valid dim
	//-----------------------------
	if (First.Col == Second.Col) {
		DTW_dist = dtw_dist(First.Mat, First.Row, Second.Mat, Second.Row, First.Col);
		return(DTW_dist);
	}
	else {
		REPORT_ERROR(SVLIB_BadArg, "DTW_Comp");	
		return (LARGE_D);
	}

}






