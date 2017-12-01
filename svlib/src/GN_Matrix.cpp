//************************************************************************
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include "SV_Error.h"
#include "GN_Matrix.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//===========================================================
//  Constructor
//===========================================================
GN_Matrix::GN_Matrix (){



}


//===========================================================
//  Destructor
//===========================================================
GN_Matrix::~GN_Matrix (){



}


/*----------------------------------------------------------*/
/* This file provides two matrix functions: inverse and det */
/*                                                          */
/* Usages:                                                  */
/*  (1) rtcode = Inverse(double **Mat, int Size)            */
/*          rtcode = 0 sucess, result in Mat,               */
/*          otherwise rtcode = -1, singular Matrix          */
/*                                                          */
/*  (2) DetValue = Det(double **Mat, int Size)              */
/*          DetValue is double precision determinant        */
/*          contents of Mat will be destroied               */
/*                                                          */
/*                                                          */
/*                                                          */
/*  Author : Jialong He                                     */
/*  Data   : March 04, 1996                                 */
/*----------------------------------------------------------*/
//===========================================================
//  Matrix inverse
//===========================================================
int GN_Matrix::Inv(double **Mat, int Dim){
	double big, pivot_inverse, temp, abs_element;
	int *pivot_flag, *swap_col, *swap_row;
	int Cnt, row, col, swap, irow, icol;

	MArray_1D(pivot_flag, Dim, int, "pivot_flag");
	MArray_1D(swap_row, Dim, int, "swap_row");
	MArray_1D(swap_col, Dim, int, "swap_col");

	//-------------------------
	// Clear to Zero
	//-------------------------
	for(Cnt = 0 ; Cnt < Dim ; Cnt++) {
		pivot_flag[Cnt] = 0;      
		swap_row[Cnt] = 0;      
		swap_col[Cnt] = 0;      
   }

  for(Cnt = 0 ; Cnt < Dim ; Cnt++) {
	  /* find the biggest pivot element */
		big = 0.0;
		for(row = 0 ; row < Dim ; row++) {
			if(!pivot_flag[row]) {       /* only unused pivots */
				for(col = 0 ; col < Dim ; col++) {
					if(!pivot_flag[col]) {
						abs_element = fabs(Mat[row][col]);
						if(abs_element >= big) {
							big = abs_element;
							irow = row;
              icol = col;
            }
          }
        }
      }
    }
    pivot_flag[icol]++;    /* mark this pivot as used */

		/* swap rows to make this diagonal the biggest absolute pivot */
    if(irow != icol) {
	    for(col = 0 ; col < Dim ; col++) {
				temp = Mat[irow][col];
				Mat[irow][col] = Mat[icol][col];
				Mat[icol][col] = temp;
      }
    }

		/* store what we swaped */
    swap_row[Cnt] = irow;
    swap_col[Cnt] = icol;

		/* bad news if the pivot is zero */
    if(Mat[icol][icol] == 0.0) {
			//REPORT_ERROR(SVLIB_Fail, "ERROR: SINGULAR MATRIX"); //by thilo: instead return -1
			return -1;
    }

		/* divide the row by the pivot */
    pivot_inverse = 1.0/Mat[icol][icol];
    Mat[icol][icol] = 1.0;        /* pivot = 1 to avoid round off */
    for(col = 0 ; col < Dim ; col++)
			Mat[icol][col] = Mat[icol][col]*pivot_inverse;

		/* fix the other rows by subtracting */
    for(row = 0 ; row < Dim ; row++) {
			if(row != icol) {
				temp = Mat[row][icol];
	      Mat[row][icol] = 0.0;
	      for(col = 0 ; col < Dim ; col++)
					Mat[row][col] = Mat[row][col]-Mat[icol][col]*temp;
			}
		}
  }

	/* fix the affect of all the swaps for final answer */
  for(swap = Dim-1 ; swap >= 0 ; swap--) {
		if(swap_row[swap] != swap_col[swap]) {
	    for(row = 0 ; row < Dim ; row++) {
				temp = Mat[row][swap_row[swap]];
				Mat[row][swap_row[swap]] = Mat[row][swap_col[swap]];
				Mat[row][swap_col[swap]] = temp;
	    }
		}
  }

  MFree_1D(pivot_flag);
  MFree_1D(swap_row);
  MFree_1D(swap_col);
  return(0);
}



/*------------------------------------------*/
/*           Matrix determinant             */
/*                                          */
/* Please note: after swapping, Mat[0]      */
/* is no longer point to the start point of */
/* memory chunk. When free memory, there may*/
/* have problem, so please save this point  */
/* NOW: I make a copy from OrigMat          */
/*------------------------------------------*/
//===========================================================
//  Matrix Determinant
//===========================================================
double GN_Matrix::Det(double **OrigMat, int Dim){
  double *a_ptr, **Mat, *MatPnt;
  double det, big, pivot_inverse, temp, abs_element;
  int row, col, swap_row, pivot;

  /*-----------------------------------------------------*/
  /* make a copy and keep the original matrix unchanged */
  /*-----------------------------------------------------*/
  MArray_2D(Mat, Dim, Dim, double, "A copy of Mat");

  MatPnt = Mat[0];       /* keep a pointer to start position */
  for (row=0; row<Dim; row++)
    for (col=0; col<Dim; col++)
			Mat[row][col] = OrigMat[row][col];

  det = 1.0;
  for(pivot = 0; pivot < Dim-1; pivot++) {
		big = fabs(Mat[pivot][pivot]);
		swap_row = 0;
		for(row = pivot + 1; row < Dim ; row++) {
	    abs_element = fabs(Mat[row][pivot]);
	    if(abs_element > big) {
				swap_row = row;
				big = abs_element;
	    }
		}

		/* unless swap_row is still zero we must swap two rows */
    if(swap_row != 0) {
			a_ptr = Mat[pivot];
	    Mat[pivot] = Mat[swap_row];
	    Mat[swap_row] = a_ptr;
	    det = -det * Mat[pivot][pivot];
		}

		else  /* the determinant by the product of the pivots */
	    det = det * Mat[pivot][pivot];

		/* if almost singular matrix, give up now */
		if(fabs(det) < 1.0e-300) {
			//REPORT_ERROR(SVLIB_DivBy0, "Singular matrix, no Det!"); //by thilo: instead: return 0, the correct det for a singular matrix
			Mat[0] = MatPnt; //restore start position //by thilo: to kill memory leak in case of error
			MFree_2D(Mat);
			return 0;
		}

		pivot_inverse = 1.0/Mat[pivot][pivot];
		for(col = pivot + 1; col < Dim; col++)
	    Mat[pivot][col] = Mat[pivot][col] * pivot_inverse;

		for(row = pivot + 1; row < Dim; row++) {
			temp = Mat[row][pivot];
			for(col = pivot + 1; col < Dim; col++)
				Mat[row][col] = Mat[row][col] - Mat[pivot][col] * temp;
		}
  }

  det = det * Mat[Dim-1][Dim-1];

  Mat[0] = MatPnt;  /* restore start position */
  MFree_2D(Mat);
  return(det);
}

//===========================================================
//  Symmetric Matrix Eigen Values and EigenVectors
//===========================================================
/*                          eigens.c
 *
 *	Eigenvalues and eigenvectors of a real symmetric matrix
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double A[n*(n+1)/2], EV[n*n], E[n];
 * void eigens( A, EV, E, n );
 *
 *
 *
 * DESCRIPTION:
 *
 * The algorithm is due to J. vonNeumann.
 *                   -     -
 * A[] is a symmetric matrix stored in lower triangular form.
 * That is, A[ row, column ] = A[ (row*row+row)/2 + column ]
 * or equivalently with row and column interchanged.  The
 * indices row and column run from 0 through n-1.
 *
 * EV[] is the output matrix of eigenvectors stored columnwise.
 * That is, the elements of each eigenvector appear in sequential
 * memory order.  The jth element of the ith eigenvector is
 * EV[ n*i+j ] = EV[i][j].
 *
 * E[] is the output matrix of eigenvalues.  The ith element
 * of E corresponds to the ith eigenvector (the ith row of EV).
 *
 * On output, the matrix A will have been diagonalized and its
 * orginal contents are destroyed.
 *
 * ACCURACY:
 *
 * The error is controlled by an internal parameter called RANGE
 * which is set to 1e-10.  After diagonalization, the
 * off-diagonal elements of A will have been reduced by
 * this factor.
 *
 * ERROR MESSAGES:
 *
 * None.
 *
 */
/*
Copyright 1973, 1991 by Stephen L. Moshier
Copyleft version.
*/

void GN_Matrix::eigens( double A[], double RR[], double E[], int N ) { //by thilo: made eigens() a class-member method

	int IND, L, LL, LM, M, MM, MQ, I, J,IA, LQ;
int IQ, IM, IL, NLI, NMI;
double ANORM, ANORMX, AIA, THR, ALM, ALL, AMM, X, Y;
double SINX, SINX2, COSX, COSX2, SINCS, AIL, AIM;
double RLI, RMI;
static double RANGE = 1.0e-10; /*3.0517578e-5;*/


/* Initialize identity matrix in RR[] */
for( J=0; J<N*N; J++ )
	RR[J] = 0.0;
MM = 0;
for( J=0; J<N; J++ )
	{
	RR[MM + J] = 1.0;
	MM += N;
	}

ANORM=0.0;
for( I=0; I<N; I++ )
	{
	for( J=0; J<N; J++ )
		{
		if( I != J )
			{
			IA = I + (J*J+J)/2;
			AIA = A[IA];
			ANORM += AIA * AIA;
			}
		}
	}
if( ANORM <= 0.0 )
	goto done;
ANORM = sqrt( ANORM + ANORM );
ANORMX = ANORM * RANGE / N;
THR = ANORM;

while( THR > ANORMX )
{
THR=THR/N;

do
{ /* while IND != 0 */
IND = 0;

for( L=0; L<N-1; L++ )
	{

for( M=L+1; M<N; M++ )
	{
	MQ=(M*M+M)/2;
	LM=L+MQ;
	ALM=A[LM];
	if( fabs(ALM) < THR )
		continue;

	IND=1;
	LQ=(L*L+L)/2;
	LL=L+LQ;
	MM=M+MQ;
	ALL=A[LL];
	AMM=A[MM];
	X=(ALL-AMM)/2.0;
	Y=-ALM/sqrt(ALM*ALM+X*X);
	if(X < 0.0)
		Y=-Y;
	SINX = Y / sqrt( 2.0 * (1.0 + sqrt( 1.0-Y*Y)) );
	SINX2=SINX*SINX;
	COSX=sqrt(1.0-SINX2);
	COSX2=COSX*COSX;
	SINCS=SINX*COSX;

/*	   ROTATE L AND M COLUMNS */
for( I=0; I<N; I++ )
	{
	IQ=(I*I+I)/2;
	if( (I != M) && (I != L) )
		{
		if(I > M)
			IM=M+IQ;
		else
			IM=I+MQ;
		if(I >= L)
			IL=L+IQ;
		else
			IL=I+LQ;
		AIL=A[IL];
		AIM=A[IM];
		X=AIL*COSX-AIM*SINX;
		A[IM]=AIL*SINX+AIM*COSX;
		A[IL]=X;
		}
	NLI = N*L + I;
	NMI = N*M + I;
	RLI = RR[ NLI ];
	RMI = RR[ NMI ];
	RR[NLI]=RLI*COSX-RMI*SINX;
	RR[NMI]=RLI*SINX+RMI*COSX;
	}

	X=2.0*ALM*SINCS;
	A[LL]=ALL*COSX2+AMM*SINX2-X;
	A[MM]=ALL*SINX2+AMM*COSX2+X;
	A[LM]=(ALL-AMM)*SINCS+ALM*(COSX2-SINX2);
	} /* for M=L+1 to N-1 */
	} /* for L=0 to N-2 */

	}
while( IND != 0 );

} /* while THR > ANORMX */

done:	;

/* Extract eigenvalues from the reduced matrix */
L=0;
for( J=1; J<=N; J++ )
	{
	L=L+J;
	E[J-1]=A[L-1];
	}
}

/*------------------------------------------*/
/* functions that needed by eignes          */
/*------------------------------------------*/
/* Multiply r (rows) by c (columns) matrix A on the left
 * by column vector V of dimension c on the right
 * to produce a (column) vector Y output of dimension r.
 */
void mvmpy(int  r, int c, double *A, double *V, double *Y ) {

double s;
double *pA, *pV, *pY;
int i, j;

pA = A;
pY = Y;
for( i=0; i<r; i++ )
	{
	pV = V;
	s = 0.0;
	for( j=0; j<c; j++ )
		{
		s += *pA++ * *pV++;
		}
	*pY++ = s;
	}
}


/* Multiply an r (rows) by c (columns) matrix A on the left
 * by a c (rows) by r (columns) matrix B on the right
 * to produce an r by r matrix Y.
 */
void mmmpy( int r, int c, double *A, double *B, double *Y ) {

register double s;
double *pA, *pB, *pY, *pt;
int i, j, k;

pY = Y;
pB = B;
for( i=0; i<r; i++ )
	{
	pA = A;
	for( j=0; j<r; j++ )
		{
		pt = pB;
		s = 0.0;
		for( k=0; k<c; k++ )
			{
			s += *pA++ * *pt;
			pt += r; /* increment to next row underneath */
			}
		*pY++ = s;
		}
	pB += 1;
	}
}


/* Transpose the n by n square matrix A and put the result in T.
 * T may occupy the same storage as A.
 */
void mtransp( int n, double *A, double *T ) {

int i, j, np1;
double *pAc, *pAr, *pTc, *pTr, *pA0, *pT0;
double x;

np1 = n+1;
pA0 = A;
pT0 = T;
for( i=0; i<n-1; i++ ) /* row index */
	{
	pAc = pA0; /* next diagonal element of input */
	pAr = pAc + n; /* next row down underneath the diagonal element */
	pTc = pT0; /* next diagonal element of the output */
	pTr = pTc + n; /* next row underneath */
	*pTc++ = *pAc++; /* copy the diagonal element */
	for( j=i+1; j<n; j++ ) /* column index */
		{
		x = *pAr;
		*pTr = *pAc++;
		*pTc++ = x;
		pAr += n;
		pTr += n;
		}
	pA0 += np1; /* &A[n*i+i] for next i */
	pT0 += np1; /* &T[n*i+i] for next i */
	}
*pT0 = *pA0; /* copy the diagonal element */
}


/* Return maximum off-diagonal element of n by n square matrix A
 */
double maxoffd(int  n,  double *A ) {

	double e, x;
	int i, j, nm1;
	double *pA;

nm1 = n-1;
e = 0.0;
pA = A;
for( i=0; i<nm1; i++ )
	{
	++pA; /* skip over the diagonal element */
	for( j=0; j<n; j++ )
		{
		x = *pA++;
		if( x < 0 )
			x = -x;
		if( x > e )
			e = x;
		}
	}
return( e );
}




/* Unpack symmetric matrix T stored in lower triangular form
 * into a symmetric n by n square matrix S.
 */
void tritosquare(int  n,  double T[], double S[]) {

double *pT;
int i, j, ni, nj;

/* offset to (i,j) element is (j*j+j)/2 + i */
pT = T;
ni = 0;
for( i=0; i<n; i++ )
	{
	nj = 0;
	for( j=0; j<i; j++ )
		{
		S[ni+j] = *pT;
		S[nj+i] = *pT++;
		nj += n;
		}
	S[ni+i] = *pT++;
	ni += n;
	}
}


/*--------------------------------------------------*/
/* Another interface which is same as my eignesym.c */
/*--------------------------------------------------*/
void GN_Matrix::Eigen(double **SymMat, double **EigVec, double *EigVal, int Dim) {

   double  *A_Mat, TmpDouble;
   int     Len = Dim * (Dim+1) / 2;
   int     Row, Col, Ind = 0;

   MArray_1D(A_Mat, Len, double, "A_Mat");

   /*--------------------------------*/
   /* copy SymMat to A_Mat           */
   /*--------------------------------*/
   for (Row=0; Row<Dim; Row++)
     for (Col=0; Col<=Row; Col++)
	A_Mat[Ind++] = SymMat[Row][Col];

   eigens(A_Mat, EigVec[0], EigVal, Dim);

   /*--------------------------------------------------------------*/
   /* tranpose ROW eigen vector to COLUMN eigne vector so that     */
   /* they are same as my eigensym.c                               */
   /*--------------------------------------------------------------*/
   for (Row=0; Row<Dim; Row++)
     for (Col=0; Col<Row; Col++) {
	 TmpDouble = EigVec[Row][Col];
	 EigVec[Row][Col] = EigVec[Col][Row];
	 EigVec[Col][Row] = TmpDouble;
     }

   MFree_1D (A_Mat);
}


/*****************************************************************

Thanx for the help, the SVD route to the pseudoinverse works like a charm!!

I did find another SVD algorithm in JC Nash's "Compact Numerical Methods for
Computers" (I think it was published in 1990).  I was not able to get the
routine in NR to work so I couldn't get a timing comparison between the two
of them.  I've checked Nash's routine against Mathematica and it's been doing
fine so far with the matricies I've tested.  I've translated Nash's code from
Pascal to C, and I've made a few minor modifications; here's the function along 
with a driver program.  Nash's program is well commented but I've left these 
out because I can't type fast. 

A couple of things to note: A needs to have twice as much room allocated for
it (2*n + 2*m) since the W in the svd function requires this (part of a 
rotation algorithm).  After the routine has run W contains two maticies
of the decomposition  A = USV'.  The first nRow rows contain the product US
and the next nCol rows contain V (not V').  Z is equal to the vector of the
sqares of the diagonal elements of S.

compile file with : gcc -o svd svd.c -lm

Let me know if you have any difficulties,
Thanx again,
Bryant
****************************************************************/
void svd_engi(double **W, double *Z, int nRow, int nCol) {

  int i, j, k, EstColRank, RotCount, SweepCount, slimit;
  double eps, e2, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2;
  eps = 1.0e-6;         /* some small tolerance value */
  slimit = nCol/4;
  if (slimit < 6.0)  slimit = 6;
  SweepCount = 0;
  e2 = 10.0 * nRow * eps * eps;
  tol = eps * 0.1;
  EstColRank = nCol;

  for (i=0; i<nCol; i++)
    for (j=0; j<nCol; j++)
    {
			W[nRow+i][j] = 0.0;
			W[nRow+i][i] = 1.0;
    }
  

	RotCount = EstColRank * (EstColRank-1)/2;
  while (RotCount != 0 && SweepCount <= slimit)
  {
    RotCount = EstColRank*(EstColRank-1)/2;
    SweepCount++;
    for (j=0; j<EstColRank-1; j++)
		{
			for (k=j+1; k<EstColRank; k++)
	    {
	      p = q = r = 0.0;
	      for (i=0; i<nRow; i++)
				{
					x0 = W[i][j]; y0 = W[i][k];
					p += x0*y0; q += x0*x0; r += y0*y0;
				}
	      Z[j] = q; Z[k] = r;
	      if (q >= r)
				{
				  if (q<=e2*Z[0] || fabs(p)<=tol*q) RotCount--;
					else
					{
						p /= q; r = 1 - r/q; vt = sqrt(4*p*p+r*r);
						c0 = sqrt(.5*(1+r/vt)); s0 = p/(vt*c0);
						for (i=0; i<nRow+nCol; i++)
						{
							d1 = W[i][j]; d2 = W[i][k];
							W[i][j] = d1*c0+d2*s0; W[i][k] = -d1*s0+d2*c0;
						}
					}
				}
	      else
				{
					p /= r; q /= (r-1); vt = sqrt(4*p*p+q*q);
					s0 = sqrt(.5*(1-q/vt));
					if (p<0) s0 = -s0;
					c0 = p/(vt*s0);
					for (i=0; i<nRow+nCol; i++)
					{
						d1 = W[i][j]; d2 =  W[i][k];
						W[i][j] = d1*c0+d2*s0; W[i][k] = -d1*s0+d2*c0;
					}
				}
	    }
		}
    while (EstColRank>=3 && Z[EstColRank]<=Z[0]*tol+tol*tol)
			EstColRank--;
  }
	
  if (SweepCount > slimit) {
		REPORT_ERROR(SVLIB_Fail, "SweepCount > slimit");
  }

}

//===========================================================
//  Singular value decomposition
//===========================================================
void GN_Matrix::SVD (double **Mat, double **US, double **S, double **V, int m, int n) {
	
	int Row, Col;
	//--------------------------------------------
	// allocate extra space and make a copy 
	//--------------------------------------------
	double **W, *Z;
	MArray_2D(W, (m+n), (n), double, "SVD:: W");
	MArray_1D(Z, (m+n), double, "SVD:: Z");
	for (Row=0; Row<m; Row++) {
		for (Col=0; Col<n; Col++) {
			W[Row][Col] = Mat[Row][Col];
		}
	}


	//-----------------------------------------
	// Do real work
	//-----------------------------------------
	svd_engi(W, Z, m, n);

	//-----------------------------------------
	// Copy back decomposed components
	//-----------------------------------------
	for (Row=0; Row<m; Row++) {
		for (Col=0; Col<n; Col++) {
			S[Row][Col] = 0.0;
		}
		S[Row][Row] = sqrt(Z[Row]);
	}
	
	for (Row=0; Row<m; Row++) {
		for (Col=0; Col<n; Col++) {
			US[Row][Col] = W[Row][Col];
		}
	}

	for (Row=0; Row<n; Row++) {
		for (Col=0; Col<n; Col++) {
			V[Row][Col] = W[Row+m][Col];
		}
	}

	MFree_2D(W);
	MFree_1D(Z);
}

//=========================================================
//
//  Linear equation solution by Gauss-Jordan elimination.
//  a[0..n-1][0..n-1] is the input matrix. b[0..n-1][0..m-1] is input containing
//  the m right-hand side vectors (Ax=b). On output, a is replaced by its
//  inverse matrix, and b is replaced by the corresponding set of solution
//  vectors.
//=========================================================
void GN_Matrix::lsolver(double **a, double **b, int n, int m) {

	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv, temp;

	MArray_1D(indxc, n, int, "indxc");
	MArray_1D(indxr, n, int, "indxr");
	MArray_1D(ipiv, n, int, "ipiv");
	for (j=0;j<n;j++) ipiv[j]=0;

	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) {
				        REPORT_ERROR(SVLIB_Fail, "GAUSSJ: Singular Matrix-1");						
					};
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) {
				temp = a[irow][l];
				a[irow][l] = a[icol][l];
				a[icol][l] = temp;
			}
			for (l=0;l<m;l++) {
				temp = b[irow][l];
				b[irow][l] = b[icol][l];
				b[icol][l] = temp;
			} 
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) {
	        REPORT_ERROR(SVLIB_Fail, "GAUSSJ: Singular Matrix-2");						
		};
		
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++) {
				temp = a[k][indxr[l]];
				a[k][indxr[l]] = a[k][indxc[l]];
				a[k][indxc[l]] = temp;
			}
	}

	MFree_1D(ipiv);
	MFree_1D(indxr);
	MFree_1D(indxc);
}



