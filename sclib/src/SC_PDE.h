/**************************************************************************/
/*    Responsibility:																											*/
/*		  - For a given Data-Matrix, compute the non-parametric Pareto      */
/*        Density Estimation (PDE). The code is based on F.Moerchen's     */
/*        Java-implementation.                                            */
/*                                                                        */
/*      - Literature:                                                     */
/*          Ultsch, A.: Optimal Density Estimation in Data containing     */
/*          Clusters of unknown Structure, Technical Report No. 34, Dept. */
/*          of Mathematics and Computer Science, University of Marburg,   */
/*          Germany, (2003)                                               */
/*          Ultsch, A.: Pareto Density Estimation: Probablity Density     */
/*          Estimation for Knowledge Discovery, In Innovations in         */
/*          Classification, Data Science, and Information Systems,        */
/*          Proceedings 27th Annual Conference of the German              */
/*          Classification Society (GfKl 2003), (2004)                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.02.2006																								*/
/**************************************************************************/

#ifndef __SC_PDE_H__
#define __SC_PDE_H__

#include "SC_TweakableParameters.h"
#include "SC_DistanceMeasures.h"
#include "SC_MatrixFunctions.h"
#include <SV_Data.h>

class SC_PDE {
	private :

  protected :

    SC_TweakableParameters *pTweak;
    SC_DistanceMeasures *pDist;
    SC_MatrixFunctions *pMat;

    //====================================================================================================================
    //	calculate the PDE for the given radius by just operating on the pairwise-distance-matrix (stored in the array)
    //====================================================================================================================
    double* getDensities(double* distances, unsigned long int distanceLength, unsigned long int dataLength, double radius, bool normalize = false);

    //====================================================================================================================
    //	find the optimal pareto radius for this data-set (represented by it's pairwise distances only)
    //====================================================================================================================
    double findParetoRadius(double* distances, unsigned long int distanceLength, unsigned long int dataLength, unsigned long int dataDim);

  public :
		
    SC_PDE(SC_TweakableParameters *pTweak = NULL);
    virtual ~SC_PDE();

    //====================================================================================================================
    //	in case the complete estimated density function is not needed but just the optimal pareto radius, this function
    //  can be called
    //====================================================================================================================
    double estimateParetoRadius(SV_Data *pData);

    //====================================================================================================================
    //  this function estimates the optimal pareto radius (if given but <= 0.0; otherwise it uses the given one) and 
    //  returns a new	 SV_Data object with the last column containing the PDE (empirical density function), all previous
    //  columns containing the original data
    //====================================================================================================================
    SV_Data* estimateDensityFunction(SV_Data *pData, double &radius);

    //====================================================================================================================
    //	for a given dataset and it's Pareto-radius, this function returns the density value of the test-vector of 
    //  dimensionality dim
    //  pData can also be an PDE estimate (with the density in the last column), if dim = pData->Col-1
    //====================================================================================================================
    double getDensity(SV_Data *pData, double radius, float* test, unsigned short int dim);

    //====================================================================================================================
    //	Visualizes the given densities in a bitmap; treturn false if writing failed
    //====================================================================================================================
		bool drawDensities(SV_Data *pDensities, char* fileName, int xSize = 512, int ySize = 512, double paretoRadius = 0.0);
		bool drawDensities(SV_Data **pMarginalDensities, unsigned int numberOfMarginals, char* fileName, int xSize = 512, int ySize = 512, double* paretoRadii = NULL);
};

#endif
