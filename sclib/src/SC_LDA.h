//************************************************************************
//  LDA-analysis   
//  calculate Transform Matrix from 2 Class
//  and transform Feature-Vektor to lower dimension.
//
//    Author  : Jun Zhou
//    Date    : April 29, 2006
//************************************************************************

#ifndef __SC_LDA_H__
#define __SC_LDA_H__

#include <SV_General.h>
#include <SV_Data.h>

class SC_LDA {

private :

protected :

public :
	//-------------------------------
	// Constructor/Destructor
	//-------------------------------
	SC_LDA();
	virtual ~SC_LDA();

	//calculate Transform Matrix for 2 classes
	SV_Data* getTransMatrix(SV_Data *class1, SV_Data *class2);

  //transform Feature-vector to lower dimension "dim"
  float* transform(float *Vector, SV_Data *transMatrix, int dim);
};// class SC_LDA
#endif
