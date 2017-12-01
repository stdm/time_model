/**************************************************************************/
/*    emd.h                                                               */
/*                                                                        */
/*    Last update: 3/24/98                                                */
/*                                                                        */
/*    An implementation of the Earth Movers Distance.                     */
/*    Based of the solution for the Transportation problem as described   */
/*    in "Introduction to Mathematical Programming" by F. S. Hillier and  */
/*    G. J. Lieberman, McGraw-Hill, 1990.                                 */
/*                                                                        */
/*    Copyright (C) 1998 Yossi Rubner                                     */
/*    Computer Science Department, Stanford University                    */
/*    E-Mail: rubner@cs.stanford.edu                                      */
/*    URL: http://vision.stanford.edu/~rubner                             */
/*                                                                        */
/*    This class is based on the above-mentioned C code by Yossi Rubner.  */
/*    It is nothing more than an object-orientaded wrapper around it with */
/*    the following changes:																							*/
/*     - The fixed static arrays of MAX_SIG_SIZE size have been	changed		*/
/*       to be allocated dynamically to reduce the memory-sage for large  */
/*			 signatures                                                       */
/*     - SC_Signature and SC_Centroid are now classes (no more structs)		*/
/*       and therefore have a method for ground-distance computation, so	*/
/*       the need for a function pointer to the distance measure has			*/
/*       vanished																													*/
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 13.06.2005																								*/
/**************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SC_EMD.h"
#include "SC_Aux.h"

//====================================================================================================================
//	The constructor is new, but rather simple...
//====================================================================================================================
SC_EMD::SC_EMD(SC_TweakableParameters *pTweak) : infinity(1e20), epsilon(1e-6) {
	this->pTweak = pTweak;
#ifdef SC_USE_EARTHMOVERSDISTANCE
	this->_C = NULL;
  this->_X = NULL;
  this->_IsX = NULL;
  this->_RowsX = NULL;
  this->_ColsX = NULL;
#else
	REPORT_ERROR(SVLIB_Fail, "EMD code not available");
#endif
}

//====================================================================================================================
//	The destructor, too :)
//====================================================================================================================
SC_EMD::~SC_EMD() {
#ifdef SC_USE_EARTHMOVERSDISTANCE
  MFree_2D(this->_C); //normally this should already be destructed, but for safety, we do it again :-)
  MFree_1D(this->_X);
  MFree_2D(this->_IsX);
  MFree_2D(this->_RowsX);
  MFree_2D(this->_ColsX);
#endif
}

/******************************************************************************
float emd(SC_Signature *Signature1, SC_Signature *Signature2, SC_Flow *Flow, 
          int *FlowSize)
  
where

   Signature1, Signature2  Pointers to signatures that their distance we want
              to compute.
   Flow       (Optional) Pointer to a vector of SC_Flow (defined in emd.h) 
              where the resulting flow will be stored. Flow must have n1+n2-1
              elements, where n1 and n2 are the sizes of the two signatures
              respectively.
              If NULL, the flow is not returned.
   FlowSize   (Optional) Pointer to an integer where the number of elements in
              Flow will be stored
              
******************************************************************************/
float SC_EMD::emd(SC_Signature *Signature1, SC_Signature *Signature2, SC_Flow *Flow, int *FlowSize)
{
#ifdef SC_USE_EARTHMOVERSDISTANCE
  int itr;
  double totalCost;
  float w;
  node2_t *XP;
  SC_Flow *FlowP;
  //node1_t U[SCLIB_EMD_MAX_SIG_SIZE1], V[SCLIB_EMD_MAX_SIG_SIZE1];
  node1_t *U = NULL, *V = NULL;

  w = init(Signature1, Signature2);

  #if SCLIB_EMD_DEBUG_LEVEL > 1
    printf("\nINITIAL SOLUTION:\n");
    printSolution();
  #endif
 
  if (this->_n1 > 1 && this->_n2 > 1)  /* IF this->_n1 = 1 OR this->_n2 = 1 THEN WE ARE DONE */
  {

		MArray_1D(U, this->_n1+1, node1_t, "SC_EMD.emd: U");
		MArray_1D(V, this->_n2+1, node1_t, "SC_EMD.emd: V");

    for (itr = 1; itr < (int)(this->pTweak->emd.maxIterations); itr++)
    {
      /* FIND BASIC VARIABLES */
			if (findBasicVariables(U, V) != true) { //by thilo: catch error and return on failure; TODO: does that prohibit a crash?
				break; 
			}
  	  
      /* CHECK FOR OPTIMALITY */
      if (isOptimal(U, V))
        break;
  	  
      /* IMPROVE SOLUTION */
      newSol();
  	  
      #if SCLIB_EMD_DEBUG_LEVEL > 1
          printf("\nITERATION # %d \n", itr);
          printSolution();
      #endif
    }

    MFree_1D(U);
    MFree_1D(V);

    if (itr == this->pTweak->emd.maxIterations)
	    fprintf(stderr, "emd: Maximum number of iterations has been reached (%d)\n", this->pTweak->emd.maxIterations);
  }

  MFree_2D(this->_IsX);

  /* COMPUTE THE TOTAL FLOW */
  totalCost = 0;
  if (Flow != NULL)
    FlowP = Flow;
  for(XP=this->_X; XP < this->_EndX; XP++)
  {
		if (XP == this->_EnterX)  /* this->_EnterX IS THE EMPTY SLOT */
			continue;
    if (XP->i == Signature1->getN() || XP->j == Signature2->getN())  /* DUMMY FEATURE */
			continue;
    if (XP->val == 0)  /* ZERO FLOW */
			continue;

    totalCost += (double)XP->val * this->_C[XP->i][XP->j];
    if (Flow != NULL)
		{
		  FlowP->from = XP->i;
			FlowP->to = XP->j;
			FlowP->amount = (float)(XP->val);
			FlowP++;
		}
  }
  if (Flow != NULL)
    *FlowSize = (int)(FlowP-Flow);

	#if SCLIB_EMD_DEBUG_LEVEL > 0
		printf("\n*** OPTIMAL SOLUTION (%d ITERATIONS): %f ***\n", itr, totalCost);
	#endif

  MFree_2D(this->_C);
  MFree_1D(this->_RowsX);
  MFree_1D(this->_ColsX);
  MFree_1D(this->_X);

  /* RETURN THE NORMALIZED COST == EMD */
  return (float)(totalCost / w);
#else
	return (float)(SVLIB_Fail);
#endif
}

#ifdef SC_USE_EARTHMOVERSDISTANCE
/**********************
   init
**********************/
float SC_EMD::init(SC_Signature *Signature1, SC_Signature *Signature2)
{
  int i, j;
  double sSum, dSum, diff;
  SC_Centroid *P1, *P2;
  //double S[SCLIB_EMD_MAX_SIG_SIZE1], D[SCLIB_EMD_MAX_SIG_SIZE1];
  double *S, *D;
 
  this->_n1 = Signature1->getN();
  this->_n2 = Signature2->getN();

  if (this->_n1 > (int)(this->pTweak->emd.maxSigSize) || this->_n2 > (int)(this->pTweak->emd.maxSigSize))
    {
      //fprintf(stderr, "emd: Signature size is limited to %d\n", SCLIB_EMD_MAX_SIG_SIZE);
      //exit(1);
      char* errMsg = new char[80];
      sprintf(errMsg, "emd: Signature size is limited to %d", this->pTweak->emd.maxSigSize);
      REPORT_ERROR(SVLIB_BadArg, errMsg);
      delete errMsg;
    }

  MArray_2D(this->_C, this->_n1+1, this->_n2+1, float, "SC_EMD.emd: _C");

  /* COMPUTE THE DISTANCE MATRIX */
  this->_maxC = 0;
  for(i=0; i < this->_n1; i++, P1++)
    for(j=0; j < this->_n2; j++, P2++) 
      {
				P1=Signature1->getCentroid(i);
				P2=Signature2->getCentroid(j);
	      this->_C[i][j] = (float)(P1->getDistance(P2));
	      if (this->_C[i][j] > this->_maxC)
	        this->_maxC = this->_C[i][j];
      }

  MArray_1D(this->_RowsX, this->_n1+1, node2_t*, "SC_EMD.init: _RowsX");
  MArray_1D(this->_ColsX, this->_n2+1, node2_t*, "SC_EMD.init: _ColsX");
  MArray_1D(S, this->_n1+1, double, "SC_EMD.init: S");
  MArray_1D(D, this->_n2+1, double, "SC_EMD.init: D");

  /* SUM UP THE SUPPLY AND DEMAND */
  sSum = 0.0;
  for(i=0; i < this->_n1; i++)
    {
      S[i] = Signature1->getWeight(i);
      sSum += Signature1->getWeight(i);
      this->_RowsX[i] = NULL;
    }
  dSum = 0.0;
  for(j=0; j < this->_n2; j++)
    {
      D[j] = Signature2->getWeight(j);
      dSum += Signature2->getWeight(j);
      this->_ColsX[j] = NULL;
    }

  /* IF SUPPLY DIFFERENT THAN THE DEMAND, ADD A ZERO-COST DUMMY CLUSTER */
  diff = sSum - dSum;
  if (fabs(diff) >= this->epsilon * sSum)
    {
      if (diff < 0.0)
	{
	  for (j=0; j < this->_n2; j++)
	    this->_C[this->_n1][j] = 0;
	  S[this->_n1] = -diff;
	  this->_RowsX[this->_n1] = NULL;
	  this->_n1++;
	}
      else
	{
	  for (i=0; i < this->_n1; i++)
	    this->_C[i][this->_n2] = 0;
	  D[this->_n2] = diff;
	  this->_ColsX[this->_n2] = NULL;
	  this->_n2++;
	}
    }

  MArray_2D(this->_IsX, this->_n1+1, this->_n2+1, char, "SC_EMD.init: _IsX");
  MArray_1D(this->_X, this->_n1+this->_n2+2, node2_t, "SC_EMD.init: _X");

  /* INITIALIZE THE BASIC VARIABLE STRUCTURES */
  for (i=0; i < this->_n1; i++)
    for (j=0; j < this->_n2; j++)
	this->_IsX[i][j] = 0;
  this->_EndX = this->_X;
   
  this->_maxW = sSum > dSum ? sSum : dSum;

  /* FIND INITIAL SOLUTION */
  russel(S, D);

  MFree_1D(S);
  MFree_1D(D);

  this->_EnterX = this->_EndX++;  /* AN EMPTY SLOT (ONLY this->_n1+this->_n2-1 BASIC VARIABLES) */

  return sSum > dSum ? (float)dSum : (float)sSum;
}


/**********************
    findBasicVariables
 **********************/
bool SC_EMD::findBasicVariables(node1_t *U, node1_t *V) //by thilo: changed return value from void to bool to indicate error
{
  int i, j, found;
  int UfoundNum, VfoundNum;
  node1_t u0Head, u1Head, *CurU, *PrevU;
  node1_t v0Head, v1Head, *CurV, *PrevV;

  /* INITIALIZE THE ROWS LIST (U) AND THE COLUMNS LIST (V) */
  u0Head.Next = CurU = U;
  for (i=0; i < this->_n1; i++)
  {
    CurU->i = i;
    CurU->Next = CurU+1;
    CurU++;
  }
  (--CurU)->Next = NULL;
  u1Head.Next = NULL;

  CurV = V+1;
  v0Head.Next = this->_n2 > 1 ? V+1 : NULL;
  for (j=1; j < this->_n2; j++)
  {
    CurV->i = j;
    CurV->Next = CurV+1;
    CurV++;
  }
  (--CurV)->Next = NULL;
  v1Head.Next = NULL;

  /* THERE ARE this->_n1+this->_n2 VARIABLES BUT ONLY this->_n1+this->_n2-1 INDEPENDENT EQUATIONS,
     SO SET V[0]=0 */
  V[0].i = 0;
  V[0].val = 0;
  v1Head.Next = V;
  v1Head.Next->Next = NULL;

  /* LOOP UNTIL ALL VARIABLES ARE FOUND */
  UfoundNum=VfoundNum=0;
  while (UfoundNum < this->_n1 || VfoundNum < this->_n2)
  {

		#if SCLIB_EMD_DEBUG_LEVEL > 3
			printf("UfoundNum=%d/%d,VfoundNum=%d/%d\n",UfoundNum,this->_n1,VfoundNum,this->_n2);
			printf("U0=");
			for(CurU = u0Head.Next; CurU != NULL; CurU = CurU->Next)
				printf("[%d]",CurU-U);
			printf("\n");
			printf("U1=");
			for(CurU = u1Head.Next; CurU != NULL; CurU = CurU->Next)
				printf("[%d]",CurU-U);
			printf("\n");
			printf("V0=");
			for(CurV = v0Head.Next; CurV != NULL; CurV = CurV->Next)
				printf("[%d]",CurV-V);
			printf("\n");
			printf("V1=");
			for(CurV = v1Head.Next; CurV != NULL; CurV = CurV->Next)
				printf("[%d]",CurV-V);
				printf("\n\n");
		#endif
      
    found = 0;
    if (VfoundNum < this->_n2)
		{
			/* LOOP OVER ALL MARKED COLUMNS */
			PrevV = &v1Head;
			for (CurV=v1Head.Next; CurV != NULL; CurV=CurV->Next)
	    {
	      j = CurV->i;
	      /* FIND THE VARIABLES IN COLUMN j */
	      PrevU = &u0Head;
	      for (CurU=u0Head.Next; CurU != NULL; CurU=CurU->Next)
				{
					i = CurU->i;
					if (this->_IsX[i][j])
					{
						/* COMPUTE U[i] */
						CurU->val = this->_C[i][j] - CurV->val;
						/* ...AND ADD IT TO THE MARKED LIST */
						PrevU->Next = CurU->Next;
						CurU->Next = u1Head.Next != NULL ? u1Head.Next : NULL;
						u1Head.Next = CurU;
						CurU = PrevU;
					}
					else
						PrevU = CurU;
				}
	      PrevV->Next = CurV->Next;
	      VfoundNum++;
	      found = 1;
			}
		}
		if (UfoundNum < this->_n1)
		{
			/* LOOP OVER ALL MARKED ROWS */
			PrevU = &u1Head;
			for (CurU=u1Head.Next; CurU != NULL; CurU=CurU->Next)
	    {
	      i = CurU->i;
	      /* FIND THE VARIABLES IN ROWS i */
	      PrevV = &v0Head;
	      for (CurV=v0Head.Next; CurV != NULL; CurV=CurV->Next)
				{
					j = CurV->i;
					if (this->_IsX[i][j])
					{
						/* COMPUTE V[j] */
						CurV->val = this->_C[i][j] - CurU->val;
						/* ...AND ADD IT TO THE MARKED LIST */
						PrevV->Next = CurV->Next;
						CurV->Next = v1Head.Next != NULL ? v1Head.Next: NULL;
						v1Head.Next = CurV;
						CurV = PrevV;
					}
					else
						PrevV = CurV;
				}
	      PrevU->Next = CurU->Next;
	      UfoundNum++;
	      found = 1;
	    }
		}
    if (! found)
    {
			//fprintf(stderr, "emd: Unexpected error in findBasicVariables!\n");
			//fprintf(stderr, "This typically happens when the this->epsilon defined in\n");
			//fprintf(stderr, "emd.h is not right for the scale of the problem.\n");
			//exit(1);
			REPORT_ERROR(SVLIB_BadData, "emd: Unexpected error in findBasicVariables! (This typically happens when the this->epsilon defined in emd.h is not right for the scale of the problem)");
			return false;
		}
  }

	return true; //by thilo: return success
}




/**********************
    isOptimal
 **********************/
int SC_EMD::isOptimal(node1_t *U, node1_t *V)
{    
  double delta, deltaMin;
  int i, j, minI, minJ;

  /* FIND THE MINIMAL Cij-Ui-Vj OVER ALL i,j */
  deltaMin = this->infinity;
  for(i=0; i < this->_n1; i++)
    for(j=0; j < this->_n2; j++)
      if (! this->_IsX[i][j])
	{
	  delta = this->_C[i][j] - U[i].val - V[j].val;
	  if (deltaMin > delta)
	    {
              deltaMin = delta;
	      minI = i;
	      minJ = j;
	    }
	}

#if SCLIB_EMD_DEBUG_LEVEL > 3
  printf("deltaMin=%f\n", deltaMin);
#endif

   if (deltaMin == this->infinity)
     {
       //fprintf(stderr, "emd: Unexpected error in isOptimal.\n");
       //exit(0);
        REPORT_ERROR(SVLIB_BadData, "emd: Unexpected error in isOptimal.");
     }
   
   this->_EnterX->i = minI;
   this->_EnterX->j = minJ;
   
   /* IF NO NEGATIVE deltaMin, WE FOUND THE OPTIMAL SOLUTION */
   return deltaMin >= -this->epsilon * this->_maxC;

/*
   return deltaMin >= -this->epsilon;
 */
}


/**********************
    newSol
**********************/
void SC_EMD::newSol()
{
    int i, j, k;
    double xMin;
    int steps;
    //node2_t *Loop[2*SCLIB_EMD_MAX_SIG_SIZE1], *CurX, *LeaveX;
    node2_t **Loop, *CurX, *LeaveX;
 
#if SCLIB_EMD_DEBUG_LEVEL > 3
    printf("EnterX = (%d,%d)\n", this->_EnterX->i, this->_EnterX->j);
#endif

    /* ENTER THE NEW BASIC VARIABLE */
    i = this->_EnterX->i;
    j = this->_EnterX->j;
    this->_IsX[i][j] = 1;
    this->_EnterX->NextC = this->_RowsX[i];
    this->_EnterX->NextR = this->_ColsX[j];
    this->_EnterX->val = 0;
    this->_RowsX[i] = this->_EnterX;
    this->_ColsX[j] = this->_EnterX;

    MArray_1D(Loop, this->_n1+this->_n2+2, node2_t*, "SC_EMD.newSol: Loop");

    /* FIND A CHAIN REACTION */
    steps = findLoop(Loop);

    /* FIND THE LARGEST VALUE IN THE LOOP */
    xMin = this->infinity;
    for (k=1; k < steps; k+=2)
      {
	if (Loop[k]->val < xMin)
	  {
	    LeaveX = Loop[k];
	    xMin = Loop[k]->val;
	  }
      }

    /* UPDATE THE LOOP */
    for (k=0; k < steps; k+=2)
      {
	Loop[k]->val += xMin;
	Loop[k+1]->val -= xMin;
      }

#if SCLIB_EMD_DEBUG_LEVEL > 3
    printf("LeaveX = (%d,%d)\n", LeaveX->i, LeaveX->j);
#endif

    /* REMOVE THE LEAVING BASIC VARIABLE */
    i = LeaveX->i;
    j = LeaveX->j;
    this->_IsX[i][j] = 0;
    if (this->_RowsX[i] == LeaveX)
      this->_RowsX[i] = LeaveX->NextC;
    else
      for (CurX=this->_RowsX[i]; CurX != NULL; CurX = CurX->NextC)
	if (CurX->NextC == LeaveX)
	  {
	    CurX->NextC = CurX->NextC->NextC;
	    break;
	  }
    if (this->_ColsX[j] == LeaveX)
      this->_ColsX[j] = LeaveX->NextR;
    else
      for (CurX=this->_ColsX[j]; CurX != NULL; CurX = CurX->NextR)
	if (CurX->NextR == LeaveX)
	  {
	    CurX->NextR = CurX->NextR->NextR;
	    break;
	  }

    MFree_1D(Loop);

    /* SET this->_EnterX TO BE THE NEW EMPTY SLOT */
    this->_EnterX = LeaveX;
}



/**********************
    findLoop
**********************/
int SC_EMD::findLoop(node2_t **Loop)
{
  int i, steps;
  node2_t **CurX, *NewX;
  //char IsUsed[2*SCLIB_EMD_MAX_SIG_SIZE1]; 
  char *IsUsed; 
 
  MArray_1D(IsUsed, this->_n1+this->_n2+2, char, "SC_EMD.findLoop: IsUsed");

  for (i=0; i < this->_n1+this->_n2; i++)
    IsUsed[i] = 0;

  CurX = Loop;
  NewX = *CurX = this->_EnterX;
  IsUsed[this->_EnterX-this->_X] = 1;
  steps = 1;

  do
    {
      if (steps%2 == 1)
	      {
	        /* FIND AN UNUSED X IN THE ROW */
	        NewX = this->_RowsX[NewX->i];
	        while (NewX != NULL && IsUsed[NewX-this->_X])
	          NewX = NewX->NextC;
	      }
            else
	      {
	        /* FIND AN UNUSED X IN THE COLUMN, OR THE ENTERING X */
	        NewX = this->_ColsX[NewX->j];
	        while (NewX != NULL && IsUsed[NewX-this->_X] && NewX != this->_EnterX)
	          NewX = NewX->NextR;
	        if (NewX == this->_EnterX)
	          break;
 	      }

      if (NewX != NULL)  /* FOUND THE NEXT X */
        {
	        /* ADD X TO THE LOOP */
	        *++CurX = NewX;
	        IsUsed[NewX-this->_X] = 1;
	        steps++;
          #if SCLIB_EMD_DEBUG_LEVEL > 3
	          printf("steps=%d, NewX=(%d,%d)\n", steps, NewX->i, NewX->j);    
          #endif
        }
      else  /* DIDN'T FIND THE NEXT X */
        {
	        /* BACKTRACK */
	        do
	          {
	            NewX = *CurX;
	            do 
	              {
		              if (steps%2 == 1)
		                NewX = NewX->NextR;
		              else
		              NewX = NewX->NextC;
	              } 
              while (NewX != NULL && IsUsed[NewX-this->_X]);
        	     
	            if (NewX == NULL)
	              {
		              IsUsed[*CurX-this->_X] = 0;
		              CurX--;
		              steps--;
	              }
	          } 
          while (NewX == NULL && CurX >= Loop);
	 
          #if SCLIB_EMD_DEBUG_LEVEL > 3
	          printf("BACKTRACKING TO: steps=%d, NewX=(%d,%d)\n",
		          steps, NewX->i, NewX->j);    
          #endif
           IsUsed[*CurX-this->_X] = 0;
	        *CurX = NewX;
	        IsUsed[NewX-this->_X] = 1;
        }     
    } 
  while(CurX >= Loop);
  
  if (CurX == Loop)
    {
      //fprintf(stderr, "emd: Unexpected error in findLoop!\n");
      //exit(1);
      REPORT_ERROR(SVLIB_BadData, "emd: Unexpected error in findLoop!");
			return steps;
    }
  #if SCLIB_EMD_DEBUG_LEVEL > 3
    printf("FOUND LOOP:\n");
    for (i=0; i < steps; i++)
      printf("%d: (%d,%d)\n", i, Loop[i]->i, Loop[i]->j);
  #endif

  MFree_1D(IsUsed);

  return steps;
}



/**********************
    russel
**********************/
void SC_EMD::russel(double *S, double *D)
{
  int i, j, found, minI, minJ;
  double deltaMin, oldVal, diff;
  //double Delta[SCLIB_EMD_MAX_SIG_SIZE1][SCLIB_EMD_MAX_SIG_SIZE1];
  double **Delta;
  //node1_t Ur[SCLIB_EMD_MAX_SIG_SIZE1], Vr[SCLIB_EMD_MAX_SIG_SIZE1];
  node1_t *Ur, *Vr;
  node1_t uHead, *CurU, *PrevU;
  node1_t vHead, *CurV, *PrevV;
  node1_t *PrevUMinI, *PrevVMinJ, *Remember;

  MArray_1D(Ur, this->_n1+1, node1_t, "SC_EMD.russel: Ur");
  MArray_1D(Vr, this->_n2+1, node1_t, "SC_EMD.russel: Vr");

  /* INITIALIZE THE ROWS LIST (Ur), AND THE COLUMNS LIST (Vr) */
  uHead.Next = CurU = Ur;
  for (i=0; i < this->_n1; i++)
    {
      CurU->i = i;
      CurU->val = -this->infinity;
      CurU->Next = CurU+1;
      CurU++;
    }
  (--CurU)->Next = NULL;
  
  vHead.Next = CurV = Vr;
  for (j=0; j < this->_n2; j++)
    {
      CurV->i = j;
      CurV->val = -this->infinity;
      CurV->Next = CurV+1;
      CurV++;
    }
  (--CurV)->Next = NULL;
  
  /* FIND THE MAXIMUM ROW AND COLUMN VALUES (Ur[i] AND Vr[j]) */
  for(i=0; i < this->_n1 ; i++)
    for(j=0; j < this->_n2 ; j++)
      {
	float v;
	v = this->_C[i][j];
	if (Ur[i].val <= v)
	  Ur[i].val = v;
	if (Vr[j].val <= v)
	  Vr[j].val = v;
      }
  
  MArray_2D(Delta, this->_n1+1, this->_n2+1, double, "SC_EMD.russel: Delta");

  /* COMPUTE THE Delta MATRIX */
  for(i=0; i < this->_n1 ; i++)
    for(j=0; j < this->_n2 ; j++)
      Delta[i][j] = this->_C[i][j] - Ur[i].val - Vr[j].val;

  /* FIND THE BASIC VARIABLES */
  do
    {
#if SCLIB_EMD_DEBUG_LEVEL > 3
      printf("Ur=");
      for(CurU = uHead.Next; CurU != NULL; CurU = CurU->Next)
	printf("[%d]",CurU-Ur);
      printf("\n");
      printf("Vr=");
      for(CurV = vHead.Next; CurV != NULL; CurV = CurV->Next)
	printf("[%d]",CurV-Vr);
      printf("\n");
      printf("\n\n");
#endif
 
      /* FIND THE SMALLEST Delta[i][j] */
      found = 0; 
      deltaMin = this->infinity;      
      PrevU = &uHead;
      for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
	{
	  int i;
	  i = CurU->i;
	  PrevV = &vHead;
	  for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
	    {
	      int j;
	      j = CurV->i;
	      if (deltaMin > Delta[i][j])
		{
		  deltaMin = Delta[i][j];
		  minI = i;
		  minJ = j;
		  PrevUMinI = PrevU;
		  PrevVMinJ = PrevV;
		  found = 1;
		}
	      PrevV = CurV;
	    }
	  PrevU = CurU;
	}
      
      if (! found)
	break;

      /* ADD X[minI][minJ] TO THE BASIS, AND ADJUST SUPPLIES AND COST */
      Remember = PrevUMinI->Next;
      addBasicVariable(minI, minJ, S, D, PrevUMinI, PrevVMinJ, &uHead);

      /* UPDATE THE NECESSARY Delta[][] */
      if (Remember == PrevUMinI->Next)  /* LINE minI WAS DELETED */
	{
	  for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
	    {
	      int j;
	      j = CurV->i;
	      if (CurV->val == this->_C[minI][j])  /* COLUMN j NEEDS UPDATING */
		{
		  /* FIND THE NEW MAXIMUM VALUE IN THE COLUMN */
		  oldVal = CurV->val;
		  CurV->val = -this->infinity;
		  for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
		    {
		      int i;
		      i = CurU->i;
		      if (CurV->val <= this->_C[i][j])
			CurV->val = this->_C[i][j];
		    }
		  
		  /* IF NEEDED, ADJUST THE RELEVANT Delta[*][j] */
		  diff = oldVal - CurV->val;
		  if (fabs(diff) < this->epsilon * this->_maxC)
		    for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
		      Delta[CurU->i][j] += diff;
		}
	    }
	}
      else  /* COLUMN minJ WAS DELETED */
	{
	  for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
	    {
	      int i;
	      i = CurU->i;
	      if (CurU->val == this->_C[i][minJ])  /* ROW i NEEDS UPDATING */
		{
		  /* FIND THE NEW MAXIMUM VALUE IN THE ROW */
		  oldVal = CurU->val;
		  CurU->val = -this->infinity;
		  for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
		    {
		      int j;
		      j = CurV->i;
		      if(CurU->val <= this->_C[i][j])
			CurU->val = this->_C[i][j];
		    }
		  
		  /* If NEEDED, ADJUST THE RELEVANT Delta[i][*] */
		  diff = oldVal - CurU->val;
		  if (fabs(diff) < this->epsilon * this->_maxC)
		    for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
		      Delta[i][CurV->i] += diff;
		}
	    }
	}
    } while (uHead.Next != NULL || vHead.Next != NULL);

    MFree_2D(Delta);
    MFree_1D(Vr);
    MFree_1D(Ur);
}




/**********************
    addBasicVariable
**********************/
void SC_EMD::addBasicVariable(int minI, int minJ, double *S, double *D, 
			     node1_t *PrevUMinI, node1_t *PrevVMinJ,
			     node1_t *UHead)
{
  double T;
  
  if (fabs(S[minI]-D[minJ]) <= this->epsilon * this->_maxW)  /* DEGENERATE CASE */
    {
      T = S[minI];
      S[minI] = 0;
      D[minJ] -= T; 
    }
  else if (S[minI] < D[minJ])  /* SUPPLY EXHAUSTED */
    {
      T = S[minI];
      S[minI] = 0;
      D[minJ] -= T; 
    }
  else  /* DEMAND EXHAUSTED */
    {
      T = D[minJ];
      D[minJ] = 0; 
      S[minI] -= T; 
    }

  /* X(minI,minJ) IS A BASIC VARIABLE */
  this->_IsX[minI][minJ] = 1; 

  this->_EndX->val = T;
  this->_EndX->i = minI;
  this->_EndX->j = minJ;
  this->_EndX->NextC = this->_RowsX[minI];
  this->_EndX->NextR = this->_ColsX[minJ];
  this->_RowsX[minI] = this->_EndX;
  this->_ColsX[minJ] = this->_EndX;
  this->_EndX++;

  /* DELETE SUPPLY ROW ONLY IF THE EMPTY, AND IF NOT LAST ROW */
  if (S[minI] == 0 && UHead->Next->Next != NULL)
    PrevUMinI->Next = PrevUMinI->Next->Next;  /* REMOVE ROW FROM LIST */
  else
    PrevVMinJ->Next = PrevVMinJ->Next->Next;  /* REMOVE COLUMN FROM LIST */
}





/**********************
    printSolution
**********************/
void SC_EMD::printSolution()
{
  node2_t *P;
  double totalCost;

  totalCost = 0;

#if SCLIB_EMD_DEBUG_LEVEL > 2
  printf("SIG1\tSIG2\tFLOW\tCOST\n");
#endif
  for(P=this->_X; P < this->_EndX; P++)
    if (P != this->_EnterX && this->_IsX[P->i][P->j])
      {
#if SCLIB_EMD_DEBUG_LEVEL > 2
	printf("%d\t%d\t%f\t%f\n", P->i, P->j, P->val, this->_C[P->i][P->j]);
#endif
	totalCost += (double)P->val * this->_C[P->i][P->j];
      }

  printf("COST = %f\n", totalCost);
}
#endif
