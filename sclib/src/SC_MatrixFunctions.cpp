/**************************************************************************/
/*    Derived from:																												*/
/*      - GN_Matrix																												*/
/*																																				*/
/*    Responsibility:																											*/
/*      - encapsulates some commonly used mathematical functions					*/
/*			- most of them are matrix-related, so this class is derived from	*/
/*				GN_Matrix to augment it's functionality													*/
/*			- there are also some statistical functions like covariance				*/
/*				and mean calculation																						*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.03.2004																								*/
/**************************************************************************/

#include <math.h>
#include "SC_MatrixFunctions.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_MatrixFunctions::SC_MatrixFunctions() {
	this->pSolver = new GN_Matrix();
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_MatrixFunctions::~SC_MatrixFunctions() {
	MFree_0D(this->pSolver);
}

//====================================================================================================================
// to hide the sclib namespace in this header
//====================================================================================================================
double SC_MatrixFunctions::getRandom(double mean, double variance) {
	return sclib::getRandomizer()->rand_gaus(mean, variance);
}

//====================================================================================================================
// Interface to more sophisticated methods provided by the GN_Matrix object
//====================================================================================================================
double SC_MatrixFunctions::det(double** squareMatrix, int dim) {
	return this->pSolver->Det(squareMatrix, dim);
}

int SC_MatrixFunctions::inv(double** squareMatrix, int dim) {
	return this->pSolver->Inv(squareMatrix, dim);
}

void SC_MatrixFunctions::eigenDecomposition(double **symmetricMatrix, int dim, double **eigenVectors, double *eigenValues) {
	this->pSolver->Eigen(symmetricMatrix, eigenVectors, eigenValues, dim);

	return;
}

//====================================================================================================================
// This version of eigenDecomposition overwrites the covariance-matrix to save the eigentvectors (thus saving some 
// memory); the code herein is merely copied from GN_Matrix.Eigen()
//====================================================================================================================
void SC_MatrixFunctions::eigenDecomposition(double **symmetricMatrix, int dim, double *eigenValues) {
	double  *A_Mat, TmpDouble;
  int Len = dim * (dim+1) / 2;
  int Row, Col, Ind = 0;

	MArray_1D(A_Mat, Len, double, "SC_MatrixFunctions::eigenDecomposition: A_Mat");

	//copy symmetricMatrix to A_Mat  (in the format needed by GN_Matrix' internal functions)
	for (Row=0; Row<dim; Row++) {
		for (Col=0; Col<=Row; Col++) {
			A_Mat[Ind++] = symmetricMatrix[Row][Col];
		}
	}

	//do the real work
	this->pSolver->eigens(A_Mat, symmetricMatrix[0], eigenValues, dim);
  MFree_1D(A_Mat);

  //tranpose ROW eigen vector to COLUMN eigne vector so that they are same as my eigensym.c                               */
	for (Row=0; Row<dim; Row++) {
		for (Col=0; Col<Row; Col++) {
			TmpDouble = symmetricMatrix[Row][Col];
			symmetricMatrix[Row][Col] = symmetricMatrix[Col][Row];
			symmetricMatrix[Col][Row] = TmpDouble;
		}
	}

	return;
}

////====================================================================================================================
////  Symmetrizes a square matrix by copying the upper triangle upon the lower one
////====================================================================================================================
//void SC_MatrixFunctions::symmetrize(double** matrix, unsigned long int dim) {
//  for (unsigned long int x = 1; x < dim; x++) {
//    for (unsigned long int y = 0; y < x; y++) {
//      matrix[y][x] = matrix[x][y];
//    }
//  }
//
//  return;  
//}

//====================================================================================================================
// Returns the log(det(squareMatrix)) via the Cholesky decomposition when squareMatrix is near singular and computing
// it the standard way would produce log(0); squareMatrix needs to be positive definite, like e.g. a covariance matrix
//====================================================================================================================
double SC_MatrixFunctions::logDet(double **squareMatrix, int dim, bool checkAssumptions) {
	double **R = choleskyDecomposition(squareMatrix, dim, checkAssumptions);
	double d = 0.0;

	//from http://www.mathworks.com/matlabcentral/fileexchange/22026:
	//The key idea of the implementation is based on the mathematical fact that the determinant of a triangular matrix equals 
	//the product of its diagonal elements. Hence, the matrix's log-determinant is equal to the sum of their logarithm values. 
	//Therefore, we can effectively tackle the problem by computing sum-of-log rather than computing log-of-product.
	//For positive definite matrices, cholesky factorization can be used.
	for (int i = 0; i < dim; i++) {
		d += sclib::sLog(R[i][i]);
	}

	MFree_2D(R);

	return 2.0 * d;
}

//====================================================================================================================
//	Checks whether the (real) matrix is positive definite; it is iff all eigenvalues of  it's symmetric part are 
//  positive
//====================================================================================================================
bool SC_MatrixFunctions::isPositiveDefinite(double **matrix, unsigned long int dim) {
  double **symm, **eigenVectors, *eigenValues, minimum;
  bool res = false;

  MArray_2D(eigenVectors, (long)dim, (long)dim, double, "SC_MatrixFunctions.isPositiveDefinite: eigenVectors");
  MArray_1D(eigenValues, dim, double, "SC_MatrixFunctions.isPositiveDefinite: eigenValues");

  //get symmetric part of matrix and it's eigenvalues
  symm = getSymmetricPart(matrix, dim);
  this->pSolver->Eigen(symm, eigenVectors, eigenValues, dim);
  MFree_2D(eigenVectors);
  MFree_2D(symm);

  //look for the smallest eigenvalue... is it positive?
  minimum = min(eigenValues, dim);
  res = (minimum >= 0.0) ? true : false;

  return res;
}

//====================================================================================================================
//  cholesky decomposition: A=L*L' returns the Cholesky factor L, a lower triangular matrix also known as A's 
//  "square-root". A must be symmetric and positive definite
//  
//  In general for i=0..dim-1 and j=i+1..dim-1 : 
//    l_ii = sqrt( a_ii - sum^{i-1}_{k=0}(l_ik)^2 )
//    l_ji = ( a_ji - sum^{i-1}_{k=0}l_jk*l_ik ) / l_ii
//  Because A is symmetric and positive definite, the expression under the square root is always positive, and all 
//  l_ij are real. 
//
//  The computation is based upon the ideas from the choldc-function from "numerical recipes in c"
//====================================================================================================================
double** SC_MatrixFunctions::choleskyDecomposition(double** matrix, unsigned long int dim, bool checkAssumptions) {
  double **L, sum;
  long int i, j, k;
  bool success = true;

  if (checkAssumptions == true) {
    if ((isPositiveDefinite(matrix, dim) == false) || (isSymmetric(matrix, dim) == false)) {
      return NULL;
    }
  }

  L = zeros(dim, dim);

  for (i=0; i < (long)dim; i++) {
    for (j=i; j < (long)dim; j++) {
      for (sum = matrix[i][j], k = i-1; k >= 0; k--) {sum -= L[i][k] * L[j][k];}
      if (i == j) {
        if (sum <= 0.0) {success = false;}
        L[i][i] = sqrt(sum);
      } else {
        L[j][i] = sum / L[i][i];
      }
    }
  }

  if (success == false) {
    MFree_2D(L);
    return NULL;
  }
  
  return L;
}

//====================================================================================================================
//  QR decomposition: A=QR where A is len*dim, Q (len*dim) is orthogonal and R (dim*dim) is upper triangular
//  taken from (but modifyed):
//
//  http://www.phys.au.dk/~bassler/numeric/ex6_2.c
//  Numerical Methods - Ugeseddel 6, Ex. 2&3
//  21.07.2004 Niels Bassler <bassler@phys.au.dk>
//  Implement QR/QL algorithm with shifts and deflation for real symmetric matrices.
//  this code is written so all matrices are of the form a[row][col]
//====================================================================================================================
void SC_MatrixFunctions::qrDecomposition(double** A, unsigned long int len, unsigned long int dim, double** &Q, double** &R) {
  unsigned long int x, y, xx;
  double norm, **tmpA = copy(A, len, dim); //A would get modified, so work on a copy instead

	MFree_2D(R);
  MFree_2D(Q);
  MArray_2D(R, (long)(dim), (long)(dim), double, "SC_MatrixFunctions.qrDecomposition: R");
  MArray_2D(Q, (long)(len), (long)(dim), double, "SC_MatrixFunctions.qrDecomposition: Q");

  for (x = 0; x < dim; x++) {
		norm = 0.0;
		for (y = 0; y < len; y++) {
			norm += tmpA[y][x] * tmpA[y][x];
		}
    R[x][x] = sqrt(norm);

    for (y = 0; y < len; y++) {
			Q[y][x] = tmpA[y][x] / R[x][x];
    }
        
    for (xx = x+1; xx < dim; xx++) {
			norm = 0.0;
			for (y = 0; y < len; y++) {
				norm += Q[y][x] * tmpA[y][xx];
			}
			if (!sclib::isFinite(norm)) {
				xx = xx;
			}
      R[x][xx] = norm; //the dot-product of Q's x'th column with A's xx'th column

      for (y = 0; y < len; y++) {
	      tmpA[y][xx] -= Q[y][x] * norm;
      }
    }
  }

  //make R upper triangular
  triu(R, dim);

	MFree_2D(tmpA);

  return;
}

//====================================================================================================================
//  This function provides an interface to the SV_Lib SVD-function, which also returns the len*dim-matrix U instead 
//  of just U*S. If U*S is not invertible, U is returned as NULL; calculating U can be switched off by setting needU
//  to false.
//====================================================================================================================
void SC_MatrixFunctions::SVD(double** M, double** &U, double** &US, double** &S, double** &V, unsigned long int len, unsigned long int dim, bool needU) {
  double **Sinv;
  int res;
  
  MFree_2D(U);
  MFree_2D(US);
  MFree_2D(S);
  MFree_2D(V);
  MArray_2D(US, (long)len, (long)dim, double, "SC_MatrixFunctions.SVD: US");
  MArray_2D(S, (long)len, (long)dim, double, "SC_MatrixFunctions.SVD: S");
  MArray_2D(V, (long)dim, (long)dim, double, "SC_MatrixFunctions.SVD: V");

  this->pSolver->SVD(M, US, S, V, len, dim);

  if (len == dim && needU == true) {
    //calc the inverse of S: S^(-1)
    Sinv = copy(S, dim, dim);
    res = this->pSolver->Inv(Sinv, dim);

    if (res == 0) {    
      //U = US * S^(-1)
      U = mult(US, Sinv, len, dim, dim, dim);
    } else {
      //S is singular => no inverse => no U => return U as NULL
      U = NULL;
    }
  }

  return;
}

//====================================================================================================================
//	return the standard deviation according to the parameters of the given inverse gaussian distribution
//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
//====================================================================================================================
double SC_MatrixFunctions::igStd(double mean, double lambda) {
	return sqrt((mean*mean*mean) / lambda);
}

//====================================================================================================================
//	return the standard deviation vector according to the parameters of the given inverse gaussian distribution
//  see Kotti, Benetos, Kotropoulos, "Compuationally Efficient and Robust BIC-Based Speaker Segmentation", 2008
//====================================================================================================================
double* SC_MatrixFunctions::igStd(double *mean, double *lambda, unsigned long int dim) {
	double *std = NULL;

	MArray_1D(std, dim, double, "SC_MatrixFunctions.igStd: std");

	for (unsigned long int x = 0; x < dim; x++) {
		std[x] = sqrt((mean[x]*mean[x]*mean[x]) / lambda[x]);
	}

	return std;
}
