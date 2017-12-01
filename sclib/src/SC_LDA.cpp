//************************************************************************
//  LDA-analysis   
//  calculate Transform Matrix and 
//  transform Feature vector to lower dimension
//
//
//    Author  : Jun Zhou
//    Date    : April 29, 2006
//************************************************************************


#include <stdlib.h>
#include <stdio.h>
#include "SC_LDA.h"
#include <SV_Error.h>

//=================================================================
//  Default Constructor
//=================================================================
SC_LDA::SC_LDA() {
	}

//=================================================================
//  Default Destructor 
//=================================================================
SC_LDA::~SC_LDA() {
}

//calculte Transformmatrix for 2 classes, and return the resulting Matrix
SV_Data* SC_LDA::getTransMatrix(SV_Data *class1, SV_Data *class2)
{
  float *x1,*x2,*w1,*w2;
	float *X,*T,*W;
	long int min;

	if( class1->Row >= class2->Row){
		min = class2->Row;
	}
	else {
		min = class1->Row;
	}

	MArray_1D(x1,min, float, "x1");
	
	MArray_1D(x2,min, float, "x2");
	
	MArray_1D(w1,min,float,"w1");
   	
	MArray_1D(w2,min,float,"w2");

	MArray_1D(X,min,float,"X");

	MArray_1D(T,min,float,"T");

	MArray_1D(W,min,float,"W");
	
	SV_Data	*result=new SV_Data(min,class1->Col);

	for(int i=0;i<min;i++){
	x1[i] = 0.0;
	x2[i] = 0.0;
	w1[i] = 0.0;
	w2[i] = 0.0;
	X[i] = 0.0;
	T[i] = 0.0;
	W[i] = 0.0;
	}
	
	if (result==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");}

		for(int FrmCnt=0;FrmCnt<min;FrmCnt++)
		{
     //calculate means of class1 and class2
			for(int n=0;n<247;n++)
			{
				x1[FrmCnt]+=class1->Mat[FrmCnt][n];
				x2[FrmCnt]+=class2->Mat[FrmCnt][n];
			}

			x1[FrmCnt]=x1[FrmCnt]/class1->Col;
			x2[FrmCnt]=x2[FrmCnt]/class2->Col;

		//calculate covariances of class1 and class2
			for(int n=0;n<class1->Col;n++)
			{
				w1[FrmCnt]+=(class1->Mat[FrmCnt][n]-x1[FrmCnt])*(class1->Mat[FrmCnt][n]-x1[FrmCnt]);
				w2[FrmCnt]+=(class2->Mat[FrmCnt][n]-x2[FrmCnt])*(class2->Mat[FrmCnt][n]-x2[FrmCnt]);
			}
     
		//calculate u3=0.5*u1+0.5*u2
			for(int n=0;n<class1->Col;n++)
			{
				X[FrmCnt]+=class1->Mat[FrmCnt][n];
			}
			for(int n=0;n<class2->Col;n++)
			{
				X[FrmCnt]+=class2->Mat[FrmCnt][n];
			}
			X[FrmCnt]=X[FrmCnt]/(class1->Col+class2->Col);

		//calculate within-class
			for(int n=0;n<class1->Col;n++)
			{
				T[FrmCnt]+=(class1->Mat[FrmCnt][n]-X[FrmCnt])*(class1->Mat[FrmCnt][n]-X[FrmCnt]);
			}
			for(int n=0;n<class2->Col;n++)
			{
				T[FrmCnt]+=(class2->Mat[FrmCnt][n]-X[FrmCnt])*(class2->Mat[FrmCnt][n]-X[FrmCnt]);
			}
        T[FrmCnt]=T[FrmCnt]/(class1->Col+class2->Col);

		//calculate between-class
			W[FrmCnt]=(w1[FrmCnt]+w2[FrmCnt])/(class1->Col+class2->Col);

			for(int j=0;j<result->Col;j++){
				result->Mat[FrmCnt][j]=0;
			}

			if(FrmCnt<=result->Col){
				result->Mat[FrmCnt][FrmCnt]=T[FrmCnt]/W[FrmCnt];
			}
		}//end of FrmCnt 
	MFree_1D(x1);
  MFree_1D(w1);
  MFree_1D(x2);
  MFree_1D(w2);
  MFree_1D(X);
  MFree_1D(T);
  MFree_1D(W);

  return result;
}


float* SC_LDA::transform(float *Vector, SV_Data *transMatrix, int dim) {
  float *newVector = NULL;
  int maxRows = min(transMatrix->Row, dim);

  MArray_1D(newVector, dim, float, "SC_LDA.transform: newVector");

  for(int m = 0; m < maxRows; m++) { 
    newVector[m] = 0.0;
		for(int n = 0; n < transMatrix->Col; n++) {
	    newVector[m] += transMatrix->Mat[m][n] * Vector[n];
	  }
  }

  return newVector;
}
