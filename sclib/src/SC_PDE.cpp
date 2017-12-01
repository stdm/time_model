/**************************************************************************/
/*    Responsibility:																											*/
/*		  - For a given Data-Matrix, copute the non-parametric Pareto       */
/*        Density Estimation (PDE). The code is based on F. Moerchen's    */
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

#include "SC_PDE.h"
#include "SC_Aux.h"
#include "SC_MatrixFunctions.h"

#include "SC_Api.h"
#ifdef SC_USE_EASYBMP
#include <EasyBMP.h>
#include <EasyBMP_Font.h>
#endif

#define SCLIB_PARETO_SIZE 0.2013

//====================================================================================================================
//	The constructor is new, but rather simple...
//====================================================================================================================
SC_PDE::SC_PDE(SC_TweakableParameters *pTweak) {
  this->pTweak = pTweak;
  this->pDist = new SC_DistanceMeasures(this->pTweak);
  this->pMat = new SC_MatrixFunctions();
}

//====================================================================================================================
//	The destructor, too :)
//====================================================================================================================
SC_PDE::~SC_PDE() {
  MFree_0D(this->pDist);
  MFree_0D(this->pMat);
}

//====================================================================================================================
//	find the optimal pareto radius for this data-set (represented by it's pairwise distances only
//====================================================================================================================
double SC_PDE::findParetoRadius(double* distances, unsigned long int distanceLength, unsigned long int dataLength, unsigned long int dataDim) {
  double radius, *percentiles = NULL, *densities = NULL;
  unsigned int percentile = 18, lastPercentile = percentile;
	double diff = 0.0, lastDiff = 1.0;
	double medianDensity, upperPercentile = 50.0, lowerPercentile = 2.0;
	bool finished = false;
	SC_MatrixFunctions mFunc;
  
  if (dataDim == 1) { //for 1D-data, the optimal radius is just the 18th percentile of the pairwise distances

    radius = sclib::percentile(distances, distanceLength, percentile, false); //don't need to compute all percentiles here, so sclib::percentile() is enough
		if (radius <= 0.0) {
			radius = 0.0;
			radius = mFunc.min(distances, distanceLength, NULL, &radius); //if the 18th percentile is 0.0, pick the next greater value from distances
			if (radius <= 0.0) {
				radius = 1.0;
			}
		}

  } else { //for multi-D-data, we have to search along the distance-percentiles for the one leading (when used as the radius) to an average of 20.13% of the overall points in the hyperspheres
    
    percentiles = sclib::percentiles(distances, distanceLength, false);
  
    while (!finished && (percentile > lowerPercentile) && (percentile < upperPercentile)) {
      radius = percentiles[percentile];
      densities = getDensities(distances, distanceLength, dataLength, radius, false); //densities are unnormalized here: just the count of points in hypersphere
      
      medianDensity = sclib::median(densities, dataLength, true) / (double)(dataLength); //median density
      diff = medianDensity - SCLIB_PARETO_SIZE;
      finished = (fabs(diff) > fabs(lastDiff)) || (diff * lastDiff < 0.0); //finished if distance got larger or sign changed

      if (fabs(diff) > fabs(lastDiff)) { //if distance got larger go back one step
        percentile = lastPercentile;
      }

      if (!finished) {
        lastDiff = diff;
        if (diff > 0) { //adjust percentile towards optimum
          upperPercentile = percentile;
          percentile = (unsigned int)(sclib::round((lowerPercentile + percentile) / 2.0));
        } else {
          lowerPercentile = percentile;
          percentile = (unsigned int)(sclib::round((upperPercentile + percentile) / 2.0));
        }
      } //if (!finished) ...

      MFree_1D(densities);
    } //while (!finished) ...

    //adjust percentile with factor for unknown number of clusters
    percentile = percentile / 3; //TODO: correction factor for cluster-count; factor 1/3 is within the 95% confidence interval for 1-13 clusters
		radius = percentiles[percentile];
		if (radius <= 0.0) {
			radius = 0.0;
			radius = mFunc.min(percentiles, 101, NULL, &radius); //if the percentile is 0.0, pick the next greater value from percentiles
			if (radius <= 0.0) {
				radius = 1.0;
			}
		}

    MFree_1D(percentiles);
  }

  return radius;  
}

//====================================================================================================================
//	calculate the PDE for the given radius and the given dimension (column) of the data-matrix
//====================================================================================================================
double* SC_PDE::getDensities(double* distances, unsigned long int distanceLength, unsigned long int dataLength, double radius, bool normalize) {
  double *densities = NULL;
  long int sum;

  if (radius <= 0.0) {
    REPORT_ERROR(SVLIB_BadArg, "Pareto-Radius can't be zero");
  } else {
    MArray_1D(densities, dataLength, double, "SC_PDE.getDensities: densities");
    
    for (unsigned long int t = 0; t < dataLength; t++) { //reuse previously computed pairwise distances
      sum = -1; // do not count the point itself
      for (unsigned long int tt = 0; tt < dataLength; tt++) {
        if (this->pDist->get1DMatrixValue(distances, dataLength, t, tt) <= radius) {  //was: distanceLength
          sum++;
        }
      }
      densities[t] = (normalize == false) ? (double)(sum) : (double)(sum) / (double)(dataLength);
    }
  } //radius > 0.0

  return densities;
}

//====================================================================================================================
//	in case the complete estimated density function is not needed but just the optimal pareto radius, this function
//  can be called
//====================================================================================================================
double SC_PDE::estimateParetoRadius(SV_Data *pData) {
  double radius, *distances = NULL, *densities = NULL;
  unsigned long int distanceLength = pData->Row*(pData->Row-1)/2;
  
  if (pData->Col == 1) {
    distances = this->pDist->getPairwiseDistances(pData->Mat, pData->Row, pData->Col, 0); //1D distance computation
  } else {
    distances = this->pDist->getPairwiseDistances(pData->Mat, pData->Row, pData->Col); //multi-D distance computation
  }
  
  radius = findParetoRadius(distances, distanceLength, pData->Row, pData->Col);
  MFree_1D(distances);
 
  return radius;
}

//====================================================================================================================
//  this function estimates the optimal pareto radius (if given but <= 0.0; otherwise it uses the given on) and 
//  returns a new	 SV_Data object with the last column containing the PDE (empirical density function), all previous
//  columns containing the original data
//====================================================================================================================
SV_Data* SC_PDE::estimateDensityFunction(SV_Data *pData, double &radius) {
  double *distances = NULL, *densities = NULL;
  unsigned long int x, y, distanceLength = pData->Row*(pData->Row-1)/2;
  SV_Data *pPDE = NULL;
  
  if (pData->Col == 1) {
    distances = this->pDist->getPairwiseDistances(pData->Mat, pData->Row, pData->Col, 0); //1D distance computation
  } else {
    distances = this->pDist->getPairwiseDistances(pData->Mat, pData->Row, pData->Col); //multi-D distance computation
  }
 
  if (radius <= 0.0) {
    radius = findParetoRadius(distances, distanceLength, pData->Row, pData->Col);
  }

  densities = getDensities(distances, distanceLength, pData->Row, radius, true);
  MFree_1D(distances);
 
  //add the calculated density estimation as the new last column to a copy of the original data-set
  pPDE = new SV_Data(pData->Row, pData->Col + 1);
  for (x = 0; x < (unsigned long int)(pPDE->Row); x++) {
    for (y = 0; y < (unsigned long int)(pPDE->Col); y++) {
      if (y < (unsigned long int)(pPDE->Col - 1)) {
        pPDE->Mat[x][y] = pData->Mat[x][y];
      } else {
        pPDE->Mat[x][y] = (float)(densities[x]);
      }
    }
  }
  MFree_1D(densities);

  return pPDE;
}

//====================================================================================================================
//	for a given dataset and it's Pareto-radius, this function returns the density value of the test-vector of 
//  dimensionality dim
//  pData can also be an PDE estimate (with the density in the last column), if dim = pData->Col-1
//====================================================================================================================
double SC_PDE::getDensity(SV_Data *pData, double radius, float* test, unsigned short int dim) {
  double density = 0.0;
  unsigned long int sum = 0;

  if (dim == 1) { //only 1 data-column and 1 containing the density estimation
    for (unsigned long int t = 0; t < (unsigned long int)(pData->Row); t++) {
      if (this->pDist->euclid(pData->Mat[t][0], *test) <= radius) { //1D distance computation
        sum++;
      }
    }
  } else {
    for (unsigned long int t = 0; t < (unsigned long int)(pData->Row); t++) {
      if (this->pDist->euclid(pData->Mat[t], test, dim) <= radius) { //multi-D distance computation
        sum++;
      }
    }
  }

  density = (double)(sum) / (double)(pData->Row);

  return density;
}

//====================================================================================================================
//	Visualizes the given densities in a bitmap
//====================================================================================================================
bool SC_PDE::drawDensities(SV_Data *pDensities, char* fileName, int xSize, int ySize, double paretoRadius) {
#ifdef SC_USE_EASYBMP
	int x, y, z, xCoord, yCoord, barHeight, h, scaledRadius;
	double *yMax, *yMin, **preBitMap;
	char label[sclib::bufferSize];
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	BMP bitMap;
	RGBApixel color; 
	int highestOrderOfMagnitude, lowestOrderOfMagnitude;
	double mean, sd;

	preBitMap = pFunc->zeros(ySize, xSize);
	yMax = pFunc->max(pDensities->Mat, pDensities->Row, pDensities->Col);
	yMin = pFunc->min(pDensities->Mat, pDensities->Row, pDensities->Col);
	scaledRadius = sclib::round(sclib::scaleToInterval((float)paretoRadius, (float)(yMin[0]), (float)(yMax[0]), (float)(0.0), (float)(xSize-1)));

	//draw 2-dimensional distribuions jointly; draw only marginals for high-dimensional data
	if (pDensities->Col == 3) { //2d-data

		for (y = 0; y < pDensities->Row; y++) {
			xCoord = sclib::round(sclib::scaleToInterval(pDensities->Mat[y][0], (float)(yMin[0]), (float)(yMax[0]), (float)(0.0), (float)(xSize-1))); //1st col => x-dim
			yCoord = ySize - sclib::round(sclib::scaleToInterval(pDensities->Mat[y][1], (float)(yMin[1]), (float)(yMax[1]), (float)(0.0), (float)(ySize-1))); //2nd col => y-dim
			for (x = xCoord-scaledRadius; x <= xCoord+scaledRadius; x++) { //draw densities in the whole circular area of the pareto-radius, not just for the point
				for (z = yCoord-scaledRadius; z <= yCoord+scaledRadius; z++) {
					if (floor(sqrt((double)(x*x+z*z))) <= scaledRadius) {
						preBitMap[z][x] += pDensities->Mat[y][2]; //the density of the point is in the last column; accumulate densities because several distinct points in feature-space may be mapped to the same pixel in bitmap-space
					}
				}
			}
		}

	} else if (pDensities->Col == 2) { //1d-data (density is expressed here in the height (y) coordinate, not in color)

		for (y = 0; y < pDensities->Row; y++) {
			xCoord = sclib::round(sclib::scaleToInterval(pDensities->Mat[y][0], (float)(yMin[0]), (float)(yMax[0]), (float)(5.0), (float)(xSize-5))); //1st col => x-dim
			yCoord = ySize - sclib::round(sclib::scaleToInterval(pDensities->Mat[y][1], (float)(0.0), (float)(yMax[1]), (float)(5.0), (float)(ySize-5))); //2nd col => density
			//for (x = xCoord-scaledRadius; x <= xCoord+scaledRadius; x++) { //draw densities in the whole area of the pareto-radius, not just for the point
				preBitMap[yCoord][xCoord] += 1; 
			//}
		}

	} else { //>=3 dimensional data

		barHeight = ySize / (pDensities->Col - 1); //draw one horizontal bar to represent each marginal distribution
		for (y = 0; y < pDensities->Row; y++) {
			for (x = 0; x < pDensities->Col-1; x++) {
				xCoord = sclib::round(sclib::scaleToInterval(pDensities->Mat[y][x], (float)(yMin[x]), (float)(yMax[x]), (float)(0.0), (float)(xSize-1)));
				for (h = x*barHeight; h < (x+1)*barHeight; h++) {
					yCoord = h;
					for (z = xCoord-scaledRadius; z <= xCoord+scaledRadius; z++) { //draw densities in the whole area of the pareto-radius, not just for the point
						preBitMap[yCoord][z] += pDensities->Mat[y][pDensities->Col-1]; //the density of the point is in the last column; accumulate densities because several distinct points in feature-space may be mapped to the same pixelk in bitmap-space
					}
				}
			}
		}

	}

	bitMap.SetSize(xSize, ySize);
	bitMap.SetBitDepth(24);

	//convert pre-bitmap to final picture (pre-bitmap was necessry to accumulate densities of points that get the same pixel-location on the bitmap but have distinct locations (and densities) in original space
	for (y = 0; y < ySize; y++) {
		for (x = 0; x < xSize; x++) {
			if (preBitMap[y][x] == 0.0) {
				bitMap(x, y)->Red = 255;
				bitMap(x, y)->Green = 255;
				bitMap(x, y)->Blue = 255;
			} else {
				if (pDensities->Col != 2) {
					bitMap(x, y)->Red = 50;
					bitMap(x, y)->Green = 50;
					bitMap(x, y)->Blue = (ebmpBYTE)(sclib::round(sclib::scaleToInterval(preBitMap[y][x], yMin[pDensities->Col-1], yMax[pDensities->Col-1], 50.0, 255.0)));
				} else {
					bitMap(x, y)->Red = 0;
					bitMap(x, y)->Green = 0;
					bitMap(x, y)->Blue = (ebmpBYTE)(255); //density is encoded in the y coordinate in 1d, not in color
				}
			}
		}
	}
	/*//red 10x10-corner beginning at 0/0 to indicate this point
	for (y = 0; y < 10; y++) {
		for (x = 0; x < 10; x++) {
			bitMap(x, y)->Red = 255;
			bitMap(x, y)->Green = 0;
			bitMap(x, y)->Blue = 0;
		}
	}*/

	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;

	//label the axis
	if (pDensities->Col == 3) { //2d
		//x-axis
		sprintf(label, "d=0; %f", yMin[0]);
		PrintString(bitMap, label, 5, ySize/2-4, 8, color); //start 5 pixel right from the left side and in the middle between top & bottom of the current image

		sprintf(label, "%f", yMax[0]);
		PrintString(bitMap, label, xSize-5-8*(int)(strlen(label)), ySize/2-4, 8, color); //stop 5 pixel left from the right side and in the middle between top & bottom of the current image

		//y-axis
		sprintf(label, "d=1; %f", yMin[1]);
		PrintString(bitMap, label, xSize/2-4*(int)(strlen(label)), ySize-5-8, 8, color);

		sprintf(label, "%f", yMax[1]);
		PrintString(bitMap, label, xSize/2-4*(int)(strlen(label)), 5, 8, color);
	} else if (pDensities->Col == 2) { //1d
		//x-axis
		sprintf(label, "%f", yMin[0]);
		PrintString(bitMap, label, 5, ySize-5-8, 8, color);

		sprintf(label, "%f", yMax[0]);
		PrintString(bitMap, label, xSize-5-8*(int)(strlen(label)), ySize-5-8, 8, color);

		//print markers for prominent points (i.e. orders of magnitude that fall into the xmax-xmin-interval)
		highestOrderOfMagnitude = (int)(floor(log(yMax[0])/log(10.0)));
		lowestOrderOfMagnitude = (int)(ceil(sclib::max(0.0, log(yMin[0])/log(10.0))));
		for (x = lowestOrderOfMagnitude; x <= highestOrderOfMagnitude; x++) {
			xCoord = sclib::round(sclib::scaleToInterval((float)(pow(10.0, (double)(x))), (float)(yMin[0]), (float)(yMax[0]), (float)(5.0), (float)(xSize-5))); //10^x
			for (y = 0; y < ySize; y++) { //draw a red line for each order of magnitude
				bitMap(xCoord, y)->Red = 255;
				bitMap(xCoord, y)->Green = 0;
				bitMap(xCoord, y)->Blue = 0;				
			}
			sprintf(label, "%.1f", pow(10.0, double(x)));
			PrintString(bitMap, label, xCoord+2, ySize/2-4, 8, color);
		}

		//print green vertical line for mean+2sd
		mean = *(pFunc->mean(pDensities->Mat, pDensities->Row, pDensities->Col, 0, 0, 0, 1));
		sd = *(pFunc->std(pDensities->Mat, pDensities->Row, pDensities->Col, &mean, 0, 0, 0, 1));
		xCoord = sclib::round(sclib::scaleToInterval((float)(mean+2*sd), (float)(yMin[0]), (float)(yMax[0]), (float)(5.0), (float)(xSize-5))); //mean + 2*sd
		for (y = 0; y < ySize; y++) { //draw a green line for each order of magnitude
			bitMap(xCoord, y)->Red = 0;
			bitMap(xCoord, y)->Green = 255;
			bitMap(xCoord, y)->Blue = 0;				
		}
		sprintf(label, "%.1f", mean+2*sd);
		PrintString(bitMap, label, xCoord+2, ySize/2-4, 8, color);

		//y-axis
		sprintf(label, "max-density: %f", yMax[1]);
		PrintString(bitMap, label, 5, 5, 8, color);
	} else {
		for (x = 0; x < (int)(pDensities->Col-1); x++) { //>=3d
			sprintf(label, "d=%d; %f", x, yMin[x]);
			PrintString(bitMap, label, 5, x*barHeight+10, 8, color); //start 5 pixel right from the left side and 2 pixel under the upper border of each bar with the min-label

			sprintf(label, "%f", yMax[x]);
			PrintString(bitMap, label, xSize-5-8*(int)(strlen(label)), x*barHeight+10, 8, color); //stop 5 pixel left from the right side and 2 pixel under the upper border of each bar with the max-label (assuming 4 pixels character-width)
		}
	}

	bitMap.WriteToFile(fileName);

	MFree_2D(preBitMap);
	MFree_1D(yMin);
	MFree_1D(yMax);
	MFree_0D(pFunc);

	return true;
#else
	return false;
#endif
}

//====================================================================================================================
//	Visualizes the given marginal densities (each array entry in pMarginalDensities holds a datamatrix with one data 
//  column and the respective densitiy-column)
//====================================================================================================================
bool SC_PDE::drawDensities(SV_Data **pMarginalDensities, unsigned int numberOfMarginals, char* fileName, int xSize, int ySize, double *paretoRadii) {
#ifdef SC_USE_EASYBMP
	int x, y, z, xCoord, yCoord, barHeight, h, *scaledRadius;
	double *yMax, *yMin, **preBitMap, *minDensity, *maxDensity;
	char label[sclib::bufferSize];
	SC_MatrixFunctions *pFunc = new SC_MatrixFunctions();
	BMP bitMap;
	RGBApixel color; 

	preBitMap = pFunc->zeros(ySize, xSize);
	MArray_1D(minDensity, numberOfMarginals, double, "SC_PDE.drawDesities: minDensity");
	MArray_1D(maxDensity, numberOfMarginals, double, "SC_PDE.drawDesities: maxDensity");
	MArray_1D(yMax, numberOfMarginals, double, "SC_PDE.drawDesities: yMax");
	MArray_1D(yMin, numberOfMarginals, double, "SC_PDE.drawDesities: yMax");
	MArray_1D(scaledRadius, numberOfMarginals, int, "SC_PDE.drawDesities: scaledRadius");
	for (x = 0; x < (int)(numberOfMarginals); x++) {
		yMax[x] = pFunc->maxOfCol(pMarginalDensities[x]->Mat, pMarginalDensities[x]->Row, 0);
		yMin[x] = pFunc->minOfCol(pMarginalDensities[x]->Mat, pMarginalDensities[x]->Row, 0);
		maxDensity[x] = pFunc->maxOfCol(pMarginalDensities[x]->Mat, pMarginalDensities[x]->Row, 1);
		minDensity[x] = pFunc->minOfCol(pMarginalDensities[x]->Mat, pMarginalDensities[x]->Row, 1);
		scaledRadius[x] = sclib::round(sclib::scaleToInterval((float)((paretoRadii != NULL) ? paretoRadii[x] : 0.0), (float)(yMin[0]), (float)(yMax[0]), (float)(0.0), (float)(xSize-1)));
	}

	//draw 2-dimensional distributions jointly; draw only marginals for high-dimensional data
	barHeight = ySize / (numberOfMarginals); //draw one horizontal bar to represent each marginal distribution
	for (y = 0; y < pMarginalDensities[0]->Row; y++) { //assume all marginals have same number of rows (must be the case!)
		for (x = 0; x < (int)(numberOfMarginals); x++) {
			xCoord = sclib::round(sclib::scaleToInterval(pMarginalDensities[x]->Mat[y][0], (float)(yMin[x]), (float)(yMax[x]), (float)(0.0), (float)(xSize-1)));
			for (h = x*barHeight; h < (x+1)*barHeight; h++) {
				yCoord = h;
				for (z = xCoord-scaledRadius[x]; z <= xCoord+scaledRadius[x]; z++) { //draw densities in the whole area of the pareto-radius, not just for the point
					preBitMap[yCoord][z] += pMarginalDensities[x]->Mat[y][1]; //the density of the point is in the last column; accumulate densities because several distinct points in feature-space may be mapped to the same pixel in bitmap-space
				}
			}
		}
	}
	
	bitMap.SetSize(xSize, ySize);
	bitMap.SetBitDepth(24);

	//convert pre-bitmap to final picture (pre-bitmap was necessry to accumulate densities of points that get the same pixel-location on the bitmap but have distinct locations (and densities) in original space
	for (y = 0; y < ySize; y++) {
		for (x = 0; x < xSize; x++) {
			if (preBitMap[y][x] == 0.0) {
				bitMap(x, y)->Red = 255;
				bitMap(x, y)->Green = 255;
				bitMap(x, y)->Blue = 255;
			} else {
				bitMap(x, y)->Red = 50;
				bitMap(x, y)->Green = 50;
				bitMap(x, y)->Blue = (ebmpBYTE)(sclib::round(sclib::scaleToInterval(preBitMap[y][x], minDensity[x], maxDensity[x], 50.0, 255.0)));
			}
		}
	}

	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;

	//label the axis
	for (x = 0; x < (int)(numberOfMarginals); x++) {
		sprintf(label, "d=%d; %f", x, yMin[x]);
		PrintString(bitMap, label, 5, x*barHeight+10, 8, color); //start 5 pixel right from the left side and 2 pixel under the upper border of each bar with the min-label

		sprintf(label, "%f", yMax[x]);
		PrintString(bitMap, label, xSize-5-8*(int)(strlen(label)), x*barHeight+10, 8, color); //stop 5 pixel left from the right side and 2 pixel under the upper border of each bar with the max-label (assuming 4 pixels character-width)
	}

	bitMap.WriteToFile(fileName);

	MFree_2D(preBitMap);
	MFree_1D(yMin);
	MFree_1D(yMax);
	MFree_1D(maxDensity);
	MFree_1D(minDensity);
	MFree_1D(scaledRadius);
	MFree_0D(pFunc);

	return true;
#else 
	return false;
#endif
}
