/**************************************************************************/
/*    This class clusters the given SC_Signatures using K-Means           */
/*                                                                        */
/*    Author  : Bing Shi																									*/
/*    Date    : 31.07.2006																								*/
/**************************************************************************/
#ifndef __SC_Clusterer_H__
#define __SC_Clusterer_H__

#include <list>
#include <limits>
#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_Centroid_Point.h"
#include "SC_Centroid.h"
#include "SC_Signature.h"
#include "SC_DistanceMeasures.h"
#include <SV_Error.h>
#include <SV_Data.h>

class SCLIB_API SC_Clusterer {
  private:
				
  protected:
		 SC_TweakableParameters *pTweak;
		 bool verbose; 
		
		//====================================================================================================================
		// Means of PData using simple K-Means
		// get: centroids: instants of Signatures these belong in one same Clusterer
		//			simpleCountK: Number of Clusters to cluster these instants
		//			ST: Number of centroids
		//		  countSign: count of Signatures
		//      weight: weights of these instants
		// return: kMeans of centroids
		//====================================================================================================================
		SC_Signature * getMeanSignatureSimpleKMeans(unsigned int simpleCountK, unsigned long int ST, unsigned int countSign, SC_Centroid **centroids, double *weights, bool useRandom, int maxIteration);

		//====================================================================================================================
		// read Centroids and its weights from given Signatures into Array of Centroids and weights
		// get: pCentroids: target Array for centorids from given Signatures
		//			pWeight: target Array for weight from given Signatures
		//			ST: Number of centroids
		//		  start: start indicise of Signatures
		//      pSign: Source Signatures
		//====================================================================================================================
		void  readCoordAndWeight(SC_Centroid** &pCentroids, double* &pWeight, SC_Signature* pSign, unsigned long int &start);

		//====================================================================================================================
		//print the Coordinates of Centroids of a Signature into a file 
		//====================================================================================================================
		void printCoord(char* fileName, SC_Signature *pData);

		//====================================================================================================================
		//print the Coordinates of Centroids into a file 
		//====================================================================================================================
		void printCoord(char *fileName, SC_Centroid_Point **cent, long int length);

		//====================================================================================================================
		//print the weights of a Signature into a file 
		//====================================================================================================================
		void printWeight(char* fileName, SC_Signature *pData);
		
		//====================================================================================================================
		//create the KMeanlist with random wise
		//get: kMeanList: target Array to save the KMean-List. KMeansList[i] = K; i->(Nr. of Signatures), K->(Nr. of Clusterer)
		//		 N: target Array to save that how many Signatures in a Cluster. N[i] = n; i->(Nr. of Clusterer), n = Num of Signatures in this Clusterer
		//		 countK: Number of Clusterers
		//		 T: Num Of Signatures
		//====================================================================================================================
		void randomDistribute(unsigned int *kMeanList, unsigned int *N, unsigned int countK, long int T);

		//====================================================================================================================
		//Calculate MeanSignature of a Clusterer using PCA
		//get: vectorMatrix: vectorMatrix[0]-> a vector of Coodinates of Centroides of a Signature.
		//return: Mean of Signatures in the same Clusterer
		//====================================================================================================================
		SC_Signature* getMeanSignaturePCA(double **vectorMatrix, int lenghthOFVector, int numOfvector); 

		//====================================================================================================================
		//convert SV_Data to Vector
		//get: pData
		//return: 1D - double Array 
		//====================================================================================================================
		double* sv_data2vector(SV_Data *pData);
		
		//====================================================================================================================
		//convert Vector SV_Data to 
		//get: 1D - double Array 
		//return: SV_Data
		//====================================================================================================================
		SV_Data* vector2sv_data(double *vector, int lengthOfVector);
		
		//====================================================================================================================
		//Quick-Sort these elements and its indices of an Array 
		//get: a: to be sorted Array
		//	 : target Array to save the sorted indices of elments
		//   : l: left index of the Array
		//   : r: right index of the Array
		//====================================================================================================================
		void sortIndex(double *a, int *indices, int l, int r);

		//====================================================================================================================
		//sort shapeList with kriterium of num of points in the shape
		//get: shapeList
		//	 : numOfPoints: the num of point of a shape
		//====================================================================================================================
		void sortShapeList(int *shapeList, int l, int r ,int *numOfPoints);

		//====================================================================================================================
		// Means of Signatures in a same Clusterer using simple K-Means
		// get: shapeList: Nr. of Signatures these belong to one same Clusterer
		//			numOfShape: Num of Signatures these belong to one same Clusterer
		//			countK: Number of Clusters to cluster these instants
		//		  pSignatures: all Signatures these will be clustered
		// return: kMeans of centroids
		//====================================================================================================================
		SC_Signature* getMeanSignature(int* shapeList, int numOfShape, int countK, SC_Signature** pSignatures);

		//====================================================================================================================
		// Means of Signatures in a same Clusterer using simple K-Means
		// get: shapeList: Nr. of Signatures these belong to one same Clusterer
		//			numOfShape: Num of Signatures these belong to one same Clusterer
		//			countK: Number of Clusters to cluster these instants
		//		  pData: all vectors in this pData these will be clustered
		//			numOfPoints: number of Points (item's value > 0.0) in a vector of pData
		//			sumOfVector: sum of the vector of pData
		// return: kMeans of centroids
		//====================================================================================================================
		SC_Signature* getMeanSignature(int* shapeList, int numOfShape, int countK, SV_Data* pData, int* numOfPoints, float* sumOfVector);

		//====================================================================================================================
		// check, if the coodinates has been allready appeared in the Array MeanPionts.
		// get: pCurrentMeanCoodinates: to be checked element
		//			pMeanPoints: to be checked Array
		//			dim: dim of this Coodinate
		//		  length: length of pMeanPoints
		// return: true->appeared, false-> not appeared
		//====================================================================================================================
		bool isContain(double* pCurrentMeanCoodinates, SC_Centroid_Point **pMeanPoints, int dim, int length);
		
		//====================================================================================================================
		//Struct to describe a Pixel(or a SC_Centroid_Point)
		//====================================================================================================================
		struct Pixel{
			double                 dist;
			double*                coordinates;
			int										 indexInVector;
			int                    dim;
			double                 weight;
			
			void setDist(double newDist){dist = newDist;}
			void setDim(int newDim){dim = newDim;}
			void setCoordinates(double* newCoordinates){coordinates = newCoordinates;}
			void setIndexInVector(int newIndex){indexInVector = newIndex;}
			void setWeight(double newWeight){weight = newWeight;}

			double getDist(){return dist;}
			double* getCoordinates(){return coordinates;}
			int getIndexInVector(){return indexInVector;}
			int getDim(){return dim;}
			double getWeight(){return weight;}
    };

		//====================================================================================================================
		// Calculate the miniDist between a point and an other point in the givn Signature
		// return: this piont in Struct Pixel with its Coordinate, weight, and dalculated dist
		//====================================================================================================================
		SC_Clusterer::Pixel getMiniDistPixel(SC_Centroid_Point* pPoint, SC_Signature* pSignature);

		//====================================================================================================================
		// Calculate the miniDist between a point and an other point in the givn vector
		// numOfPionts: num of points (item's value > 0.0)in this vector
		// sumOfVector: sum of vector in order to cal. the weight for each point of vector
		// return: this piont in Struct Pixel with its Coordinate, weight, and dalculated dist
		//====================================================================================================================
		SC_Clusterer::Pixel getMiniDistPixel(float* pPoint, int dimOfPiont, float* vector, float sumOfVector, int dimOfVector);

		//====================================================================================================================
		// Initializes kMeans randomly, returning the found means
		//====================================================================================================================
		SV_Data* randomInitialization(unsigned int k, SV_Data *pData);

		//====================================================================================================================
		// Initializes kMeans deterministically employing Var-Part algorithm from "In Search of Deterministic Methods for 
		// Initializing k-Means and Gaussian Mixture Clustering", Su, Dy, 2006
		// The method returns the found means
		//====================================================================================================================
		SV_Data* deterministicInitialization(unsigned int k, SV_Data *pData);

  public:

    SC_Clusterer(SC_TweakableParameters *pTweak, bool verbose = true);
    ~SC_Clusterer();

		//====================================================================================================================
		//This method is clustering the gived SC_Signatures using K-Means
		// get the SC_Signature ** pData
		// return kMeans of pData
		// countK : Number of KMeans(Clusters) T: Number of Signatures
		// useRandom :assign each signature randomly or no randomly via the useRandom to a class mean, remember this mapping 
		//            in kMeansList. true->randomly  false->not randomly
		// maxIteration: the max iterations for running of this alg     
		// iterationToIniKmeanList: the max iterations to ini the kmeanlist
		//====================================================================================================================
		SC_Signature** kMeans(int k,  SC_Signature **pData, int T , double varianceLimit, bool useRandom, int maxIteration, int iterationToIniKMeanList);

		//====================================================================================================================
		// This method is clustering the given data using K-Means
		// get the SV_Data * pData
		// return k Means of pData
		// countK : Number of KMeans(Clusters) T: Number of Signatures
		// useRandom :assign each signature randomly or no randomly via the useRandom to a class mean, remember this mapping 
		//            in kMeansList. true->randomly  false->not randomly
		// maxIteration: the max iterations for running of this alg     
		// iterationToIniKmeanList: the max iterations to ini the kmeanlist
		//====================================================================================================================
		SC_Signature** kMeans(int k, SV_Data *pData, int T, double varianceLimit, bool useRandom, int maxIteration, int iterationToIniKMeanList);

		//====================================================================================================================
		// Takes the data in pData and returns the k Centroids determined using #iterations of k-means with Euclidean distance
		//====================================================================================================================
		SV_Data* kMeans(unsigned int k, SV_Data *pData, unsigned int iterations = 10, bool randomInit = false, unsigned int repeatedInitializations = 1, bool replaceWithMedoids = false);
		SV_Data* kMeans(unsigned int k, SV_Data *pData, unsigned short int* &clusterMapping, unsigned long int* &N, unsigned int iterations = 10, bool randomInit = false, unsigned int repeatedInitializations = 1, bool replaceWithMedoids = false);

		//====================================================================================================================
		// Gets a set of "templates" in pMeans and replaces each vector in pData with the nearest template vector according
		// to the Euclidean distance; better search strategy than brute force: use ANN (approximate nearest neighbor search)
		//====================================================================================================================
		void pickNearest(SV_Data *pData, SV_Data *pMeans);

		//====================================================================================================================
		// returns the index of the nearest template (brute force)
		//====================================================================================================================
		template<class T> unsigned int pickNearest(T *vector, SV_Data *pMeans) {
			double minDist = std::numeric_limits<double>::max(), dist;
			int minIdx, i;

			for (i = 0; i < pMeans->Row; i++) {
				dist = SC_DistanceMeasures::euclid(vector, pMeans->Mat[i], pMeans->Col);
				if (dist < minDist) {
					minDist = dist;
					minIdx = i;
				}
			}

			return (unsigned int)(minIdx);
		}
};
#endif
