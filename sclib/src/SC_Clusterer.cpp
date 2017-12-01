/**************************************************************************/
/*    This class clusters the given SC_Signatures using K-Means           */
/*                                                                        */
/*    Author  : Bing Shi																									*/
/*    Date    : 31.07.2006																								*/
/**************************************************************************/

//#include <stdlib.h> //for NULL macro
#include <math.h>
#include <vector>
#include "SC_Clusterer.h"
#include "SC_Aux.h" //for sclib::bufferSize define
#include "SC_Centroid_Point.h"
#include "SC_SDP.h"
#include "SC_Timer.h"
#include "SC_MatrixFunctions.h"
#include <SV_Error.h>
#include <GN_Filter.h>
#include <SV_Data.h>

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_Clusterer::SC_Clusterer(SC_TweakableParameters *pTweak, bool verbose) {
	this->pTweak = pTweak;
	this->verbose = verbose;
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_Clusterer::~SC_Clusterer() {

}

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
SC_Signature** SC_Clusterer::kMeans (int countK, SC_Signature **pData, int T , double varianceLimit, bool useRandom, int maxIteration, int iterationToIniKmeanList) {
	SC_Signature** means;
	SC_Signature** tempMeans;
	long int  t, minDist;
	countK = (countK <= T)? countK : T;
	unsigned int *kMeansList, *bestKMeanList; // KMeansList[i] = K; i->(Nr. of Signatures), K->(Nr. of Clusterer)
	unsigned int *N, *bestN; //N[i] = n; i->(Nr. of Clusterer), n = Num of Signatures in this Clusterer
	unsigned int k, r;
	double* dist; // distant between two Signatures
	unsigned long int *meanCN, *coordCount;//contains the count of signatures in the same cluster
	//unsigned long int maxCNOfClusters = 0; // the max number of Centroids over all cluster
	int *shapeList; // the Numbers of Signatures in one same Cluster, the max length is as length of K-MeansList;
	int i;
	//unsigned long int ** pBMP;
	SC_DistanceMeasures* distMeasures = new SC_DistanceMeasures(this->pTweak);
	SC_Timer* time = new SC_Timer(countK+1, this->pTweak);
	SC_MatrixFunctions *pFunctions = new SC_MatrixFunctions();
	SC_SDP *pSDP = new SC_SDP(this->pTweak->featureSdp.m, this->pTweak->featureSdp.lag, this->pTweak->featureSdp.color, this->pTweak->featureSdp.n, this->pTweak->featureSdp.pictureSize, this->pTweak->featureSdp.tau);
	
	MArray_1D(tempMeans, countK, SC_Signature*, "SC_Clusterer.KMeans:: means");
	MArray_1D(dist, countK, double, "SC_Clusterer.KMeans:: dist");
	MArray_1D(kMeansList, T, unsigned int, "SC_Clusterer.KMeans:: kMeansList");
	MArray_1D(N, countK, unsigned int, "SC_Clusterer.KMeans:: N");
	MArray_1D(meanCN, countK, unsigned long int, "SC_Clusterer.KMeans::meanCN");
	MArray_1D(coordCount, countK,unsigned long  int, "SC_Clusterer.KMeans:: CoordCount");
	
	//init N, meanCN and dist with zeros
	time->startTimer(1); // time0 is for the all Alg.
	for (k = 0; k < (unsigned int)countK; k++) {
		N[k] = 0;
		meanCN[k] = 0;
		dist[k] = 0.0;
		tempMeans[k] = NULL;
	}

	if (useRandom){
		printf("\n			Phase IV_2:\tIni-KeamList...");
		double *sumOfDist, *bestDist;
		MArray_1D(sumOfDist, countK, double, "SC_Clusterer::KMeans. sumOfDist");
		MArray_1D(bestDist, countK, double, "SC_Clusterer::KMeans. bestDist");
		MArray_1D(bestKMeanList, T, unsigned int, "SC_Clusterer::KMeans. bestKMeanList");
		MArray_1D(bestN, countK, unsigned int, "SC_Clusterer::KMeans. bestN");

		for(k = 0; k < (unsigned int)countK; k++){
			bestDist[k] = std::numeric_limits<double>::max();
		}
		
		randomDistribute(kMeansList, N, countK, T);

		for (k = 0; k < (unsigned int)countK; k++) {
				bestKMeanList[k] = kMeansList[k];
			}

		//serch the best KMeanList
		for(i = 0; i < iterationToIniKmeanList; i++){
			
			printf("\n				Ini-KeamList Iteration %d...", i);
			// ini the meanCN and coordCount;
			for (k = 0; k < (unsigned int)countK; k++) {
				coordCount[k] = 0;
			}
			for (t = 0; t < T; t++ ) {
				coordCount[kMeansList[t]] += pData[t]->getN();
			}
			for (k = 0; k < (unsigned int)countK; k++) {
				meanCN[k] = coordCount[k] / N[k];
			}

			//calculate the K-Means
			for (k = 0; k < (unsigned int)countK; k++) {
				MArray_1D(shapeList, N[k], int, "SC_Clusterer::KMeans.shapeList");
				int index = 0;
				for (t = 0; t < T; t++) {
					if (k == kMeansList[t]){
						shapeList[index++]= t;
					}
				}// t = 1...T
				MFree_0D(tempMeans[k]);
				tempMeans[k] = getMeanSignature(shapeList, N[k], meanCN[k], pData);
				MFree_1D(shapeList);
			}//k = 1...K

			//calculate the sum of dist
			pFunctions->fillWithValue(sumOfDist, countK, 0.0);
			for (t = 0; t < T; t++) {
				sumOfDist[kMeansList[t]] += distMeasures->beigi(tempMeans[kMeansList[t]], pData[t]);
			}//t = 0 .... T
			for (k = 0; k < (unsigned int)countK; k++ ){
				sumOfDist[k] /= N[k];
			}//k = 0 ..

			//updata the bestKmeanList
			int n = 0;
			for(k = 0; k < (unsigned int)countK; k++){
				if(sumOfDist[k] < bestDist[k])
					n++;
			}
			if(n >= countK * 0.66){
				for(k = 0; k < (unsigned int)countK; k++){
					bestDist[k] = sumOfDist[k];
				}
				for(t = 0; t < T; t++){
					bestKMeanList[t] = kMeansList[t];
				}
				for(k = 0; k < (unsigned int)countK; k++){
					bestN[k] = N[k];
				}
			}

			//freshen the kMeanList
			randomDistribute(kMeansList, N, countK, T);
			while(pFunctions->equals(bestKMeanList, countK, kMeansList)){
				randomDistribute(kMeansList, N, countK, T);
			}
			//sclib::vectorOut("kmeanList1.txt",kMeansList,T);
			//sclib::vectorOut("bestList.txt",bestKMeanList,T);
		}
		
		//save the bestKmeanList into K-MeanList
		for(t = 0; t < T; t++){
			kMeansList[t] = bestKMeanList[t];
		}
		for(k = 0; k < (unsigned int)countK; k++){
			N[k] = bestN[k];
		}
		MFree_1D(bestKMeanList);
		MFree_1D(bestN);
		MFree_1D(sumOfDist);
		MFree_1D(bestDist);
		for(k = 0; k < (unsigned int)countK; k++){
			MFree_0D(tempMeans[k]);
		}
		MFree_1D(tempMeans);

	}else{
		for (t = 0; t < T; t++) {
			r = (unsigned short int)(t % countK);
			N[r]++; //count the number of signatures used to estimate this mean of the cluster; 
			kMeansList[t] = r; //remember to which signature was assigned
		}
	}

	MArray_1D(means, countK, SC_Signature*, "SC_Clusterer.KMeans:: means");
	for (k = 0; k < (unsigned int)(countK); k++) {
		means[k] = NULL;
	}

	printf("\n			Phase IV_2:\tClustering...");

	for (i = 0; i < maxIteration; i++) {
		printf("\n				Clutering Iteration %d...", i);
		// ini the meanCN and coordCount;
		for (k = 0; k < (unsigned int)countK; k++) {
			coordCount[k] = 0;
		}

		for (t = 0; t < T; t++ ) {
			coordCount[kMeansList[t]] += pData[t]->getN();
		}
		for (k = 0; k < (unsigned int)countK; k++) {
			meanCN[k] = coordCount[k] / N[k];
		}

		//calculate the K-Means
		for (k = 0; k < (unsigned int)countK; k++) {
			printf("\n					Clusterer %d...", i);
			MArray_1D(shapeList, N[k], int, "SC_Clusterer::KMeans.shapeList");
			int index = 0;
			for (t = 0; t < T; t++) {
				if (k == kMeansList[t]){
					shapeList[index++]= t;
				}
			}// t = 1...T
			MFree_0D(means[k]);
			means[k] = getMeanSignature(shapeList, N[k], meanCN[k], pData);	
			MFree_1D(shapeList);
			
		}//k = 1...K

		//calculate the minDis		
		for (t = 0; t < T; t++) {
			minDist = 0;
			for (k = 0; k < (unsigned int)countK; k++ ){
				dist[k] = distMeasures->beigi(means[k], pData[t]);
				if (dist[k] <= dist[minDist]) {
					minDist = k;
				}
			}//k = 0 ... countK

			//update the KMeanList
			if(N[kMeansList[t]] > 1 && kMeansList[t]!= minDist){
				N[kMeansList[t]]--; //there is now one member less in this cluster
				N[minDist]++; //and this cluster has got one more
				kMeansList[t] = (unsigned short)minDist; //note that this Signature is now a member of the new Mean
			}//end if
		}//t = 0... T	
		sclib::vectorOut("kmeanList.txt",kMeansList,T);
	}

	time->stopTimer(1);
	time->logIt();
	MFree_1D(shapeList);
	MFree_1D(dist);
	MFree_1D(N);
	MFree_1D(kMeansList);
	MFree_1D(meanCN);
	MFree_1D(coordCount);
	MFree_0D(distMeasures);
	MFree_0D(time);
	MFree_0D(pFunctions);
	MFree_0D(pSDP);

	

	return means;
}
//====================================================================================================================
// Means of PData using simple K-Means
// get: centroids: instants of Signatures these belong in one same Clusterer
//			simpleCountK: Number of Clusters to cluster these instants
//			ST: Number of centroids
//		  countSign: count of Signatures
//      weight: weights of these instants
// return: kMeans of centroids
//====================================================================================================================
SC_Signature* SC_Clusterer::getMeanSignatureSimpleKMeans(unsigned int simpleCountK, unsigned long int ST, unsigned int countSign, SC_Centroid **centroids, double *weights, bool useRandom, int maxIteration){
	SC_Signature* pSignature = NULL;
	unsigned long int  t, k, miniDist;
	simpleCountK = (simpleCountK <= ST)? simpleCountK : ST;
	unsigned int *simpleKMeansList;
	double *simpleDist, *meanWeight, *weightNorFactor;// Factor of the normalize for every Clusterer
	unsigned int *simpleN, r;//contains the count of coordinates in the same cluster, r is the Number of SimpleKMeans;
	int i, n;
	double sumWeight = 0.0, w, wMean;
	SC_Timer* ts = new SC_Timer(5, this->pTweak);

	MArray_1D(simpleKMeansList, ST, unsigned  int, "SC_Clusterer.simpleKMeans: simpleKMeansList");
	MArray_1D(simpleDist, simpleCountK, double, "SC_Clusterer.simpleKMeans: simpleDist");
	MArray_1D(meanWeight, simpleCountK, double, "SC_Clusterer.simpleKMeans: meanWeight");
	MArray_1D(weightNorFactor, simpleCountK, double, "SC_Clusterer.simpleKMeans: weightNorFactor");
	MArray_1D(simpleN, simpleCountK, unsigned int, "SC_Clusterer.simpleKMeans: simpleN");

	switch (centroids[0]->getCentroidType()) {

		case sclib::centroidPoint: {
			ts->startTimer(1);// for init
			SC_Centroid_Point **centrPoint = (SC_Centroid_Point **) centroids;
			SC_Centroid_Point **simpleMean;
			double *point, *pointMean;
			int dim = centrPoint[0]->getDim();
			double tempPoint[2];// cood for the Centroid
			MArray_1D(simpleMean,simpleCountK, SC_Centroid_Point *, "SC_Clusterer.simpleKMeans: simpleMean");
			unsigned short int oldGd = this->pTweak->distanceMeasure.groundDistance;
			this->pTweak->distanceMeasure.groundDistance = sclib::dmEuclid;

			//printCoord("centroids.txt", centrPoint, 2);

			//init with zeros
			for (k = 0; k < simpleCountK; k++) {
				simpleN[k] = 0;
				for (i = 0; i < dim; i++) {
					tempPoint[i] = 0.0;
				}
				simpleMean[k] = new SC_Centroid_Point(this->pTweak, dim, tempPoint, false);
				simpleDist[k] = 0.0;
				meanWeight[k] = 0.0;
				weightNorFactor[k] = 0.0;
			}

			if (useRandom){
				randomDistribute(simpleKMeansList, simpleN, simpleCountK, ST); 			
			}else{
				for (t = 0; t < ST; t++) {
					r = (unsigned int)(t % simpleCountK);
					simpleN[r]++; //count the number of coord used to estimate this simpleMean of the cluster; 
					simpleKMeansList[t] = r; //remember to which Coord was assigned
				}
			}

			//normalize of weights, mean weight and weightNorFactor
			for (t = 0; t < ST; t++){
				weights[t] /= countSign; 	//normalize these weight
				weightNorFactor[simpleKMeansList[t]] += weights[t];
				meanWeight[simpleKMeansList[t]] += weights[t] ;
			}

			//mean of meanWeight and the Factor for the Normalize of weights
			for (k = 0; k < simpleCountK; k++) {
				meanWeight[k] /= simpleN[k];
			}

			//intial Mean Point
			for(t = 0; t < ST; t++){
				point = centrPoint[t]->getCoordinate(); // set the coord
				pointMean = simpleMean[simpleKMeansList[t]]->getCoordinate();
				for (int d = 0; d < dim; d++)	{	
					pointMean[d] +=  point[d] * weights[t] * (1 / weightNorFactor[simpleKMeansList[t]]);
				}
				simpleMean[simpleKMeansList[t]]->setCoordinate(pointMean);
			}

			ts->stopTimer(1);
			ts->startTimer(2);

			//simpel k-means to update means
			for (n = 0; n < maxIteration; n++) {
				for (t = 0; t < ST; t++) {
					ts->startTimer(3);
					miniDist = 0;
					point = centrPoint[t]->getCoordinate();
					w = weights[t];
					//find the closest centroid (mean) for the Coords, remember it's nr in "minDist"
					for (k = 0; k < simpleCountK; k++) {
						ts->startTimer(4);
						pointMean = simpleMean[k]->getCoordinate();
						simpleDist[k] = simpleMean[k]->getDistance(centrPoint[t]);
						ts->stopTimer(4);
						if(simpleDist[k] < simpleDist[miniDist]) {
							miniDist = k;
						}
					}//k = 0....simpleCountK
					
					//move the point form old Cluster to new Cluster
					ts->startTimer(5);
					if((simpleN[simpleKMeansList[t]] > 1) && (simpleKMeansList[t] != miniDist) ) {
						// remove this point
						pointMean = simpleMean[simpleKMeansList[t]]->getCoordinate();
						wMean = meanWeight[simpleKMeansList[t]];
						double factor = weightNorFactor[simpleKMeansList[t]];
						
						for (i = 0; i < dim; i++) {
							double p = pointMean[i];
							double temp = (p/(1/factor)) - (point[i]*w);
							double temp1 = (1/(factor-w));
							pointMean[i] = temp * temp1;
							if(pointMean[i] < 0.0){
								pointMean[i] = 0.0;
							}
						}
						simpleMean[simpleKMeansList[t]]->setCoordinate(pointMean);

						weightNorFactor[simpleKMeansList[t]] -= w;
						meanWeight[simpleKMeansList[t]] =  (1.0 / (simpleN[simpleKMeansList[t]] - 1.0)) * (simpleN[simpleKMeansList[t]] * wMean - weights[t]);
						
						//add it to new one
						pointMean = simpleMean[miniDist]->getCoordinate();
						wMean = meanWeight[miniDist] ;
						factor = weightNorFactor[miniDist];
						for (i = 0; i < dim; i++) {
							double p = pointMean[i];
							double temp = (p/(1/factor)) + (point[i]*w);
							double temp1 = (1/(factor+w));
							pointMean[i] = temp * temp1;
							//pointMean[i] = ((pointMean[i]/(1/factor)) + (point[i]*w)) * (1/(factor+w));
						}
						simpleMean[miniDist]->setCoordinate(pointMean);

						weightNorFactor[miniDist] += w;
						meanWeight[miniDist] =  (1.0 / (simpleN[miniDist] + 1.0)) * (simpleN[miniDist] *wMean + weights[t]);
						

						simpleN[simpleKMeansList[t]]--; //there is now one member less in this cluster
						simpleN[miniDist]++; //and this cluster has got one more
						simpleKMeansList[t] = miniDist; //note that this point is now a member of the new simpleMean
					}
					ts->stopTimer(5);
					ts->stopTimer(3);
				}// t = 0...ST
			}// n = 0...maxIteration

			ts->stopTimer(2);
		
			//calculate the sum of Mean Weights in order to normalize them
			double sumOfMWeight = 0.0;
			for(k = 0; k < simpleCountK; k++) {
				sumOfMWeight += meanWeight[k];
			}
			for(k = 0; k < simpleCountK; k++) {
				meanWeight[k] /= sumOfMWeight;
			}
			//printCoord("simpleMean.txt", simpleMean, simpleCountK);
			
			//simpleMean of  Signature
			pSignature = new SC_Signature((SC_Centroid**)simpleMean, meanWeight, simpleCountK, false, 1.0);
			MFree_1D(simpleKMeansList);
			MFree_1D(simpleDist);
			MFree_1D(simpleN);

			this->pTweak->distanceMeasure.groundDistance = oldGd;
			break;
															 }

		default: {
			pSignature = NULL;
			break;
						 }
	}

	MFree_1D(meanWeight);

	return pSignature;
}

//====================================================================================================================
// read Centroids and its weights from given Signatures into Array of Centroids and weights
// get: pCentroids: target Array for centorids from given Signatures
//			pWeight: target Array for weight from given Signatures
//			ST: Number of centroids
//		  start: start indicise of Signatures
//      pSign: Source Signatures
//====================================================================================================================
void  SC_Clusterer::readCoordAndWeight(SC_Centroid** &pCentroids, double* &pWeight, SC_Signature* pSign, unsigned long int &start) {
	long int i, length = pSign->getN();
	for(i = 0; i < length; i++){
		pCentroids[start+i] = pSign->getCentroid(i);
		pWeight[start+i] = pSign->getWeight(i);
	}
	start = start + length;
}

//====================================================================================================================
//print the Coordinates of Centroids of a Signature into a file 
//====================================================================================================================
void SC_Clusterer::printCoord(char* fileName, SC_Signature *pData){
	int t;
	for(t = 0; t < pData->getN(); t++){
	 SC_Centroid_Point** temp = 	(SC_Centroid_Point**)pData->getCentroids();
		sclib::vectorOut(fileName, temp[t]->getCoordinate(), 2);
	}
}

//====================================================================================================================
//print the weights of a Signature into a file 
//====================================================================================================================
void SC_Clusterer::printWeight(char* fileName, SC_Signature *pData){
		double* temp = pData->getWeight();
		sclib::vectorOut(fileName,temp, pData->getN(),true);
}

//====================================================================================================================
//print the Coordinates of Centroids into a file 
//====================================================================================================================
void SC_Clusterer::printCoord(char* fileName, SC_Centroid_Point **cent, long int length) {
	long int t;
	int *temp;
	double *temp1;
	MArray_1D(temp, 2, int, "");
	for(t = 0; t < length; t++){
		temp1  = cent[t]->getCoordinate();
		for(int d = 0; d < 2; d++){
			temp[d] = (int)temp1[d];
		}
		sclib::vectorOut(fileName, temp, 2);
	}
}

//====================================================================================================================
//create the KMeanlist with random wise
//get: kMeanList: target Array to save the KMean-List. KMeansList[i] = K; i->(Nr. of Signatures), K->(Nr. of Clusterer)
//		 N: target Array to save that how many Signatures in a Cluster. N[i] = n; i->(Nr. of Clusterer), n = Num of Signatures in this Clusterer
//		 countK: Number of Clusterers
//		 T: Num Of Signatures
//====================================================================================================================
void SC_Clusterer::randomDistribute(unsigned int *kMeansList, unsigned  int *N, unsigned int countK, long int T){
	unsigned int* array1;
	unsigned int r, temp, k;
	MArray_1D(array1, T,  unsigned int, "SC_Clusterer.KMeans:: temp");
	//ini this array1
	for ( int i = 0; i < T; i++)
	{
		array1[i] = i;
	}

	//permotation this contents of this array1 with random wise
	for (int i = 0;  i < T; i++)
	{
		r = (int)(floor(((float)(rand()) / (float)(RAND_MAX+1)) * T));
		temp = array1[i];
		array1[i] = array1[r];
		array1[r] = temp;
	}

	for(k = 0; k< countK; k++){
		N[k] = 0;
	}

	unsigned int avgCount = T / countK;
	long int n = 0;
	//give this contents to every clusterer
	for (k = 0; k < countK; k++){
		for( unsigned int i = 0; i < avgCount; i++){
			kMeansList[array1[n]] = (unsigned)k;
			N[k]++;
			n++;
		}
	}
	k = 0;
	for (n; n < T; n++){
		kMeansList[array1[n]] = k;
		N[k]++;
		k++;
	}
	MFree_1D(array1);
}

//====================================================================================================================
//Calculate MeanSignature of a Clusterer using PCA
//get: vectorMatrix: vectorMatrix[0]-> a vector of Coodinates of Centroides of a Signature.
//return: Mean of Signatures in the same Clusterer
//====================================================================================================================
SC_Signature* SC_Clusterer::getMeanSignaturePCA(double **vectors, int lengthOfVector, int numOfVector){
	double *meanOfvectors, *eigenValues;
	double **covAndEigen;
	double *meanPictureVector;
	SC_SDP *pSDP = new SC_SDP(this->pTweak->featureSdp.m, this->pTweak->featureSdp.lag, this->pTweak->featureSdp.color, this->pTweak->featureSdp.n, this->pTweak->featureSdp.pictureSize, this->pTweak->featureSdp.tau);
	int x, y;
	
	SC_MatrixFunctions *pFunctions = new SC_MatrixFunctions(); 
	
	//normalize all vectors
	meanOfvectors = pFunctions->mean(vectors, numOfVector, lengthOfVector);

	for(x = 0; x < numOfVector; x++){
		for(y = 0; y < lengthOfVector; y++){
			vectors[x][y] -= meanOfvectors[y];
		}
	}

	covAndEigen = pFunctions->cov(vectors, numOfVector, lengthOfVector, meanOfvectors);

	MArray_1D(eigenValues, lengthOfVector, double, "SC_Clusterer::meanPicturePCA.eigenValues");
	pFunctions->eigenDecomposition(covAndEigen, lengthOfVector, eigenValues);

	int* indices = NULL;
	MArray_1D(indices, lengthOfVector, int, "SC_Clusterer::meanPicturePCA.indices");
	for(x = 0; x < lengthOfVector; x++){
		indices[x] = x;
	}

	//select the eigenVector with max eigenValue as mean-picture
	int index = pFunctions->maxIdx(eigenValues, lengthOfVector);
	//meanPictureVector = pFunctions->getRow(covAndEigen, lengthOfVector, lengthOfVector, index);
	
	sortIndex(eigenValues, indices, 0, lengthOfVector-1);
	
	meanPictureVector = pFunctions->getRow(covAndEigen, lengthOfVector, lengthOfVector, indices[lengthOfVector-1]);
	
	//meanPicture to Signature
	SV_Data *pData = vector2sv_data(meanPictureVector, lengthOfVector);
	
	SC_Signature *pSignature = pSDP->sv_data2signature(pData);
	
	MFree_0D(pFunctions);
	MFree_0D(pSDP);
	MFree_1D(indices);
	MFree_1D(eigenValues);
	MFree_2D(covAndEigen);
	
	return pSignature;
}
//====================================================================================================================
//convert SV_Data to Vector
//get: pData
//return: 1D - double Array 
//====================================================================================================================
double* SC_Clusterer::sv_data2vector(SV_Data *pData){
	int x, y;
	double *vector;
	if(pData->Row != pData->Col){REPORT_ERROR(SVLIB_BadArg, "The Matrix of pData must be a symmetric matrix");}

	MArray_1D(vector, pData->Row * pData->Col, double, "SC_Clusterer::sv_data2vector.vector");
	int index = 0;
	for(x = 0; x < pData->Row; x++){
		for(y = 0; y < pData->Col; y++){
			vector[index++] = pData->Mat[x][y];
		}
	}
	return vector;	
}

//====================================================================================================================
//convert Vector SV_Data to 
//get: 1D - double Array 
//return: SV_Data
//====================================================================================================================
SV_Data* SC_Clusterer::vector2sv_data(double *vector, int lengthOfVector){
	SV_Data * pData = NULL;
	pData = new SV_Data();
	if (pData==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}

	pData->Row = (int)sqrt((double)lengthOfVector);
	pData->Col = (int)sqrt((double)lengthOfVector);
	if(pData->Row != pData->Col){REPORT_ERROR(SVLIB_BadArg, "The Matrix of pData must be a symmetric matrix");}

	pData->Alloc();
	pData->Hdr.frameSize = this->pTweak->featureSdp.frameSize;
	pData->Hdr.frameStep = this->pTweak->featureSdp.frameStep;
  unsigned long x, y;
	
	// copy the vector to pData
	int index = 0;
	for (x = 0; x < (unsigned long)pData->Row; x++){
		for(y = 0; y < (unsigned long)pData->Col; y++){
			if(vector[index] != 0.0)
				pData->Mat[x][y] = 1.0;
			else
				pData->Mat[x][y] = 0.0;
			index++;
		}
	}
 return pData;
}

//====================================================================================================================
//Quick-Sort these elements and its indices of an Array 
//get: a: to be sorted Array
//	 : target Array to save the sorted indices of elments
//   : l: left index of the Array
//   : r: right index of the Array
//====================================================================================================================
void SC_Clusterer::sortIndex(double *a, int *indices, int l, int r){
	int i, j;
	double	x;		//storage for comparison-element
	double	temp; //storage for temp-element for triangle-exchange
	int tempI;
	if ( l>=r ) { //nothing to do here
		return;
	} 

	i=l;
	j=r;
	x = a[(l+r)/2];

	//start of division-step
	while (i <= j) {
		while (a[i] < x) {i++;}
		while (a[j] > x) {j--;}
		if (i<=j)	{ //do the triangle-exchange
			temp = a[i];
			a[i] = a[j];
			a[j] = temp;
			tempI = indices[i];
			indices[i] = indices[j];
			indices[j] = tempI;
			i++; 
			j--;
		}
	}

	sortIndex (a,indices, l, j);
	sortIndex (a,indices, i, r);
	return;
}

//====================================================================================================================
//sort shapeList with kriterium of num of points in the shape
//get: shapeList
//	 : numOfPoints: the num of point of a shape
//====================================================================================================================
void SC_Clusterer::sortShapeList(int *shapeList, int l, int r, int *numOfPoints){
	int i, j;
	int	x;		//storage for comparison-element
	int	temp; //storage for temp-element for triangle-exchange
	if ( l>=r ) { //nothing to do here
		return;
	} 

	i=l;
	j=r;
	x = shapeList[(l+r)/2];

	//start of division-step
	while (i <= j) {
		while (numOfPoints[shapeList[i]] < numOfPoints[x]) {i++;}
		while (numOfPoints[shapeList[j]] > numOfPoints[x]) {j--;}
		if (i<=j)	{ //do the triangle-exchange
			temp = shapeList[i];
			shapeList[i] = shapeList[j];
			shapeList[j] = temp;
			i++; 
			j--;
		}
	}

	sortShapeList (shapeList, l, j, numOfPoints);
	sortShapeList (shapeList, i, r, numOfPoints);
	return;
}


//====================================================================================================================
// Means of Signatures in a same Clusterer using simple K-Means
// get: shapeList: Nr. of Signatures these belong to one same Clusterer
//			numOfShape: Num of Signatures these belong to one same Clusterer
//			countK: Number of Clusters to cluster these instants
//		  pSignatures: all Signatures these will be clustered
// return: kMeans of centroids
//====================================================================================================================
SC_Signature* SC_Clusterer::getMeanSignature(int* shapeList, int numOfShape, int countK, SC_Signature** pSignatures){
	SC_Signature *pSignature = NULL;
	
	if(numOfShape == 1){
		SC_Centroid_Point *pPoint = NULL;
		SC_Centroid_Point **pMeanPoints = NULL;
		double* pMeanWeigths;
		int dim, numOfCentroid = pSignatures[shapeList[0]]->getN();
		MArray_1D(pMeanPoints,countK, SC_Centroid_Point*, "SC_Clusterer::meanSignature.meanPoint");

		for(int k = 0; k < numOfCentroid; k++){
			pPoint  = (SC_Centroid_Point*)(pSignatures[shapeList[0]]->getCentroid(k));
			dim = pPoint->getDim();
			pMeanPoints[k] = new SC_Centroid_Point(this->pTweak, dim, pPoint->getCoordinate(), false);
		}
		pMeanWeigths=  pSignatures[shapeList[0]]->getWeight();
		pSignature = new SC_Signature((SC_Centroid**)pMeanPoints, pMeanWeigths, numOfCentroid, false,  pSignatures[shapeList[0]]->getSmallestUnnormalizedWeight());
		return pSignature;
	}else {
		SC_Clusterer::Pixel **pPixels;
		int x, y, k;
		SC_Signature *pMaxShape;
		int numOfPixel, temp, index, dim;
		double *sumOfminiDist, *pTempPoint, *pMeanWeigths, *pCoodinates, *pMeanCoodinates, sumOfWeights, factor = 0.0;
		int *indices;
		SC_Centroid_Point* pPoint;
		SC_Centroid_Point **pMeanPoints;
		SC_MatrixFunctions *pFunctions = new SC_MatrixFunctions();
		
		//sort the shapeList with num of pixels
		for(int i = 0; i < numOfShape; i++){
			for(int j = i; j < numOfShape; j++){
				if(pSignatures[shapeList[i]]->getN() > pSignatures[shapeList[j]]->getN()){
					temp = shapeList[i];
					shapeList[i] = shapeList[j];
					shapeList[j] = temp;
				}
			}
		}

		pMaxShape = pSignatures[shapeList[numOfShape-1]];
		numOfPixel = pMaxShape->getN();

		MArray_2D(pPixels, numOfPixel, numOfShape-1, SC_Clusterer::Pixel,"SC_Clusterer::getMeanSignature.pixels");
		sumOfminiDist = pFunctions->initVector(numOfPixel, 0.0);

		for(x = 0; x < numOfPixel; x++){
			pPoint = (SC_Centroid_Point*)pMaxShape->getCentroid(x);
			for(y = 0; y < numOfShape-1; y++){
				pPixels[x][y] = getMiniDistPixel(pPoint, pSignatures[shapeList[y]]);
				sumOfminiDist[x] += pPixels[x][y].getDist();
				//factor += pPixels[x][y].getWeight();
			}
		}

		MArray_1D(indices, numOfPixel, int, "SC_Clusterer::getMeanSignature.indices");
		for(x = 0;  x < numOfPixel; x++){
			indices[x] = x;
		}

		//sort the sumOfMinDist
		sortIndex(sumOfminiDist, indices, 0, numOfPixel-1);
		MFree_1D(sumOfminiDist);
		
		dim = pPoint->getDim();
		MArray_1D(pMeanPoints,countK, SC_Centroid_Point*, "SC_Clusterer::meanSignature.meanPoint");
		MArray_1D(pMeanWeigths,countK, double, "SC_Clusterer::meanSignature.pMeanWeigths");
		MArray_1D(pTempPoint, dim, double, "SC_Clusterer::meanSignature.pTempPoint");

		//init with zeros
		for (k = 0; k < countK; k++) {
			for (int d = 0; d < dim; d++) {
				pTempPoint[d] = 0.0;
			}
			pMeanPoints[k] = new SC_Centroid_Point(this->pTweak, dim, pTempPoint, false);
			pMeanWeigths[k] = 0.0; 
		}
		
		//select the frist countK indices of pixel and col the mean
		int j = 0; 
		double *pCurrentMeanCoodinates = pFunctions->initVector(dim, 0.0);
		double currentMeanWeights;
		double weight = 0.0; 
		for(int i = 0; i < countK; i++){
			do{
					index = indices[j];
					
					for(int d = 0; d < dim; d++){
						pCurrentMeanCoodinates[d] = 0.0;
					}
					currentMeanWeights = 0.0;
					
					sumOfWeights = pMaxShape->getWeight(index);
					for(y = 0; y < numOfShape - 1; y++){
						weight = pPixels[index][y].getWeight();
						sumOfWeights += weight;
					}

					pPoint = (SC_Centroid_Point*)pMaxShape->getCentroid(index);
					pCoodinates = pPoint->getCoordinate();
					weight = pMaxShape->getWeight(index);
					for(int d = 0; d < dim; d++){
						pCurrentMeanCoodinates[d] +=  pCoodinates[d] * weight * (1 / sumOfWeights);
					}

					currentMeanWeights += weight;

					for(y = 0; y < numOfShape - 1; y++){
						pCoodinates = pPixels[index][y].getCoordinates();
						weight = pPixels[index][y].getWeight();
						for(int d = 0; d < dim; d++){
							pCurrentMeanCoodinates[d] +=  pCoodinates[d] * weight * (1 / sumOfWeights);
						}
						currentMeanWeights += weight;
						currentMeanWeights /= numOfShape;
					}
					j++;
			}while((i > 0) && isContain(pCurrentMeanCoodinates, pMeanPoints, dim, i)&& ((countK - i) < numOfPixel - j));
			
			pMeanCoodinates = pMeanPoints[i]->getCoordinate();
			for(int d = 0; d < dim; d++){
				pMeanCoodinates[d] = pCurrentMeanCoodinates[d];
			}
			//pMeanPoints[i]->setCoordinate(pMeanCoodinates);
			pMeanWeigths[i] = currentMeanWeights;
		}

		//normalize the weights
		sumOfWeights = 0.0;
		for(int i = 0; i < countK; i++){
			sumOfWeights += pMeanWeigths[i];
		}

		for(int i = 0; i < countK; i++){
			pMeanWeigths[i] /= sumOfWeights ;
		}

		sumOfWeights = 0.0;
		for(int i = 0; i < countK; i++){
			sumOfWeights += pMeanWeigths[i];
		}

		pSignature = new SC_Signature((SC_Centroid**)pMeanPoints, pMeanWeigths, countK, true, 1.0);
		pSignature->setJustLink(false, false); //this way, the pMeanPoints get killed when the signature gets killed but we create no copy here

		/*for(x = 0; x < numOfPixel; x++){
			for(y = 0; y < numOfShape-1; y++){
				MFree_0D(*pPixels[x][y]);
			}
		}*/
		MFree_2D(pPixels);
		MFree_0D(pFunctions);
		MFree_1D(pTempPoint);
		MFree_1D(indices);
		MFree_1D(pCurrentMeanCoodinates);

		return pSignature;
	}
}

//====================================================================================================================
// Calculate the miniDist between a point and an other point in the givn Signature
// return: this piont in Struct Pixel with its Coordinate, weight, and dalculated dist
//====================================================================================================================
SC_Clusterer::Pixel SC_Clusterer::getMiniDistPixel(SC_Centroid_Point* pPoint, SC_Signature* pSignature){
	double miniDist = std::numeric_limits<double>::max(), dist = 0.0;
	int index = -1;
	SC_Clusterer::Pixel pPixle;
	SC_Centroid_Point* temp;

	unsigned short int oldGd = this->pTweak->distanceMeasure.groundDistance;
	this->pTweak->distanceMeasure.groundDistance = sclib::dmEuclid;

	for(int n = 0; n < pSignature->getN(); n++){
		dist = pPoint->getDistance((SC_Centroid_Point*)pSignature->getCentroid(n));
		if(dist < miniDist){
			miniDist = dist;
			index = n;
		}
	}
	this->pTweak->distanceMeasure.groundDistance = oldGd;

	temp = (SC_Centroid_Point*)(pSignature->getCentroid(index));
	//pPixle = * new SC_Clusterer::Pixel();
	pPixle.setCoordinates(temp->getCoordinate());
	pPixle.setDim(temp->getDim());
	pPixle.setDist(miniDist);
	pPixle.setWeight(pSignature->getWeight(index));
	return pPixle;
}

//====================================================================================================================
// Calculate the miniDist between a point and an other point in the givn vector
// numOfPionts: num of points (item's value > 0.0)in this vector
// sumOfVector: sum of vector in order to cal. the weight for each point of vector
// return: this piont in Struct Pixel with its Coordinate, weight, and dalculated dist
//====================================================================================================================
SC_Clusterer::Pixel SC_Clusterer::getMiniDistPixel(float* pPoint, int dimOfPiont, float* vector, float sumOfVector, int dimOfVector){
	double miniDist = std::numeric_limits<double>::max(), dist = 0.0;
	int squareSize = (int)(sqrt((double)dimOfVector));
	int index = -1;
	SC_Clusterer::Pixel pPixle;
	float coord[2];

	for(int n = 0; n < dimOfVector; n++){
		if(vector[n] > 0){
			coord[0] = (float)(n / squareSize);
			coord[1] = (float)(n % squareSize);
			dist = SC_DistanceMeasures::euclid(pPoint, coord, dimOfPiont);
			if(dist < miniDist){
				miniDist = dist;
				index = n;
			}
		}
	}
	pPixle.setIndexInVector(index);
	pPixle.setDim(dimOfPiont);
	pPixle.setDist(miniDist);
	pPixle.setWeight(vector[index] / sumOfVector);
	return pPixle;
}
//====================================================================================================================
// check, if the coodinates has been allready appeared in the Array MeanPionts.
// get: pCurrentMeanCoodinates: to be checked element
//			pMeanPoints: to be checked Array
//			dim: dim of this Coodinate
//		  length: length of pMeanPoints
// return: true->appeared, false-> not appeared
//====================================================================================================================
bool SC_Clusterer::isContain(double* pCurrentMeanCoodinates, SC_Centroid_Point **pMeanPoints, int dim, int length){
	int res, a, b;
	double* pTemp = NULL;

	for(int i = 0; i < length; i++){
		pTemp = pMeanPoints[i]->getCoordinate();
		res = 0;
		for(int d = 0; d < dim; d++){
			a = (int)pCurrentMeanCoodinates[d]; 
			b = (int)pTemp[d];
			if(a == b){
				res++;
			}
		}
		if(res == dim) return true;
	}
	return false;
}

//====================================================================================================================
//This method is clustering the gived SC_Signatures using K-Means
// get the SV_Data * pData
// return kMeans of pData
// countK : Number of KMeans(Clusters) T: Number of Signatures
// useRandom :assign each signature randomly or no randomly via the useRandom to a class mean, remember this mapping 
//            in kMeansList. true->randomly  false->not randomly
// maxIteration: the max iterations for running of this alg     
// iterationToIniKmeanList: the max iterations to ini the kmeanlist
//====================================================================================================================
SC_Signature** SC_Clusterer::kMeans(int countK, SV_Data *pData, int T, double varianceLimit, bool useRandom, int maxIteration, int iterationToIniKmeanList) {
	SC_Signature** means;
	SC_Signature** tempMeans;
	long int  t, minDist;
	countK = (countK <= T)? countK : T;
	unsigned int *kMeansList, *bestKMeanList, *tempKMeanList;// KMeansList[i] = K; i->(Nr. of Signatures), K->(Nr. of Clusterer)
	unsigned int *N, *bestN; //N[i] = n; i->(Nr. of Clusterer), n = Num of Signatures in this Clusterer
	unsigned int k, r;
	double *dist; // distant between two Signatures
	unsigned long int *meanCN, *coordCount;//contains the count of signatures in the same cluster
	//unsigned long int maxCNOfClusters = 0; // the max number of Centroids over all cluster
	int i, *shapeList; // the Numbers of Signatures in one same Cluster, the max length is as length of K-MeansList;
	float *sumOfVector; //sum of vector in order to normalize weights of every vector in the pData
	int *numOfPionts; // number of Pionts in a vector of the pData
	int dimOfvector = pData->Col;
	float *pVector; // point of a Vector in the pData;
	double *wS, *wV;
	
	SC_DistanceMeasures* distMeasures = new SC_DistanceMeasures(this->pTweak);
	SC_Timer* time = new SC_Timer(countK+1, this->pTweak);
	SC_MatrixFunctions *pFunctions = new SC_MatrixFunctions();
	SC_SDP *pSDP = new SC_SDP(this->pTweak->featureSdp.m, this->pTweak->featureSdp.lag, this->pTweak->featureSdp.color, this->pTweak->featureSdp.n, this->pTweak->featureSdp.pictureSize, this->pTweak->featureSdp.tau);

	MArray_1D(tempMeans, countK, SC_Signature*, "SC_Clusterer.KMeans:: means");
	MArray_1D(dist, countK, double, "SC_Clusterer.KMeans:: dist");
	MArray_1D(kMeansList, T, unsigned int, "SC_Clusterer.KMeans:: kMeansList");
	MArray_1D(tempKMeanList, T, unsigned int, "SC_Clusterer.KMeans:: tempKMeanList");
	MArray_1D(N, countK, unsigned int, "SC_Clusterer.KMeans:: N");
	MArray_1D(meanCN, countK, unsigned long int, "SC_Clusterer.KMeans::meanCN");
	MArray_1D(coordCount, countK,unsigned long  int, "SC_Clusterer.KMeans:: CoordCount");
	MArray_1D(sumOfVector, T, float, "SC_Clusterer.KMeans:: sumOfvector");
	MArray_1D(numOfPionts, T, int, "SC_Clusterer.KMeans:: numOfPixels");
	MArray_1D(wS, dimOfvector, double, "SC_Clusterer.KMeans:: wS");
	MArray_1D(wV, dimOfvector, double, "SC_Clusterer.KMeans:: wV");

	//init N, meanCN and dist with zeros
	time->startTimer(1); // time0 is for the all Alg.
	for (k = 0; k < (unsigned int)countK; k++) {
		N[k] = 0;
		meanCN[k] = 0;
		dist[k] = 0.0;
		tempMeans[k] = NULL;
	}
	
	//col the sumOfvector and numOfPixels
	for(t = 0; t < T; t++){
		pVector = pData->Mat[t];
		sumOfVector[t] = 0;
		numOfPionts[t] = 0;
		for(i = 0; i < dimOfvector; i++){
			if(pVector[i] > 0){
				sumOfVector[t] += pVector[i];
				numOfPionts[t] += 1;
			}
		}
	}

	if (useRandom){
		printf("\n\n	Phase IV_1:\tIni-KeamList...");
		double sumOfDist, bestDist = std::numeric_limits<double>::max();
		//MArray_1D(sumOfDist, countK, double, "SC_Clusterer::KMeans. sumOfDist");
		//MArray_1D(bestDist, countK, double, "SC_Clusterer::KMeans. bestDist");
		MArray_1D(bestKMeanList, T, unsigned int, "SC_Clusterer::KMeans. bestKMeanList");
		MArray_1D(bestN, countK, unsigned int, "SC_Clusterer::KMeans. bestN");

		randomDistribute(kMeansList, N, countK, T);

		//serch the best KMeanList
		for(i = 0; i < iterationToIniKmeanList; i++){
			
			printf("\n	--- Ini-KeamList Iteration %d...", i);
			// ini the meanCN and coordCount;
			for (k = 0; k < (unsigned int)countK; k++) {
				coordCount[k] = 0;
			}
			for (t = 0; t < T; t++ ) {
				coordCount[kMeansList[t]] += numOfPionts[t];
			}
			for (k = 0; k < (unsigned int)countK; k++) {
				meanCN[k] = coordCount[k] / N[k];
			}

			//calculate the K-Means
			for (k = 0; k < (unsigned int)countK; k++) {
				MArray_1D(shapeList, N[k], int, "SC_Clusterer::KMeans.shapeList");
				int index = 0;
				for (t = 0; t < T; t++) {
					if (k == kMeansList[t]){
						shapeList[index++]= t;
					}
				}// t = 1...T
				MFree_0D(tempMeans[k]);
				tempMeans[k] = getMeanSignature(shapeList, N[k], meanCN[k], pData, numOfPionts, sumOfVector);
				MFree_1D(shapeList);
			}//k = 1...K

			//calculate the sum of dist
			//pFunctions->fillWithValue(sumOfDist, countK, 0.0);
			sumOfDist = 0.0;
			for (t = 0; t < T; t++) {
				sumOfDist += distMeasures->beigi(tempMeans[kMeansList[t]], pData->Mat[t], dimOfvector, wS, wV);
			}//t = 0 .... T

			//updata the bestKmeanList
			if(sumOfDist < bestDist){
				bestDist = sumOfDist;
				for(t = 0; t < T; t++){
					bestKMeanList[t] = kMeansList[t];
				}
				for(k = 0; k < (unsigned int)countK; k++){
					bestN[k] = N[k];
				}
			}

			sclib::vectorOut("kmeanList1.txt",kMeansList,T);
			sclib::vectorOut("bestList.txt",bestKMeanList,T);

			//freshen the kMeanList
			randomDistribute(kMeansList, N, countK, T);
			while(pFunctions->equals(bestKMeanList, countK, kMeansList)){
				randomDistribute(kMeansList, N, countK, T);
				sclib::vectorOut("kmeanList2.txt",kMeansList,T);
			}
			
			printf("done!\n");
		}
		
		//save the bestKmeanList into K-MeanList
		for(t = 0; t < T; t++){
			kMeansList[t] = bestKMeanList[t];
		}
		for(k = 0; k < (unsigned int)countK; k++){
			N[k] = bestN[k];
		}
		MFree_1D(bestKMeanList);
		MFree_1D(bestN);

		for(k = 0; k < (unsigned int)countK; k++){
			MFree_0D(tempMeans[k]);
		}
		MFree_1D(tempMeans);
		printf("done!\n");

	}else{
		for (t = 0; t < T; t++) {
			r = (unsigned short int)(t % countK);
			N[r]++; //count the number of signatures used to estimate this mean of the cluster; 
			kMeansList[t] = r; //remember to which signature was assigned
		}
	}

	MArray_1D(means, countK, SC_Signature*, "SC_Clusterer.KMeans:: means");
	for (k = 0; k < (unsigned int)(countK); k++) {
		means[k] = NULL;
	}

	printf("\n	Phase IV_2:\tClustering...");

	sclib::vectorOut("kmeanList.txt",kMeansList,T);
	for (i = 0; i < maxIteration; i++) {
		printf("\n	--- Clutering Iteration %d...", i);
		// ini the meanCN and coordCount;
		for (k = 0; k < (unsigned int)countK; k++) {
			coordCount[k] = 0;
		}

		for (t = 0; t < T; t++ ) {
			coordCount[kMeansList[t]] += numOfPionts[t];
			tempKMeanList[t] = kMeansList[t];
		}
		for (k = 0; k < (unsigned int)countK; k++) {
			meanCN[k] = coordCount[k] / N[k];
		}

		//calculate the K-Means
		for (k = 0; k < (unsigned int)countK; k++) {
			printf("\n	--- Clusterer %d...", k);
			MArray_1D(shapeList, N[k], int, "SC_Clusterer::KMeans.shapeList");
			int index = 0;
			for (t = 0; t < T; t++) {
				if (k == kMeansList[t]){
					shapeList[index++]= t;
				}
			}// t = 1...T
			MFree_0D(means[k]);
			means[k] = getMeanSignature(shapeList, N[k], meanCN[k], pData, numOfPionts, sumOfVector);	
			MFree_1D(shapeList);
			printf("done!\n");
		}//k = 1...K

		//calculate the minDis	
		printf("\n	MiniDis ...");
		for (t = 0; t < T; t++) {
			minDist = 0;
			for (k = 0; k < (unsigned int)countK; k++ ){
				dist[k] = distMeasures->beigi(means[k], pData->Mat[t], dimOfvector, wS, wV);
				if (dist[k] < dist[minDist]) {
					minDist = k;
				}
			}//k = 0 ... countK

			//update the KMeanList
			if(N[kMeansList[t]] > 1 && kMeansList[t]!= minDist){
				N[kMeansList[t]]--; //there is now one member less in this cluster
				N[minDist]++; //and this cluster has got one more
				kMeansList[t] = (unsigned short)minDist; //note that this Signature is now a member of the new Mean
			}//end if
		}//t = 0... T	
		printf("done!\n");

		sclib::vectorOut("kmeanList.txt",kMeansList,T);
		
		//check, if the current KMeanList equlas to the old KMeanlist.
		if(pFunctions->equals(kMeansList, T, tempKMeanList)){
			time->stopTimer(1);
			time->logIt();
			MFree_1D(shapeList);
			MFree_1D(dist);
			MFree_1D(N);
			MFree_1D(kMeansList);
			MFree_1D(meanCN);
			MFree_1D(coordCount);
			MFree_0D(distMeasures);
			MFree_0D(time);
			MFree_0D(pFunctions);
			MFree_0D(pSDP);
			MFree_1D(numOfPionts);
			MFree_1D(sumOfVector);
			MFree_1D(tempKMeanList);
			MFree_1D(wS);
			MFree_1D(wV);
			printf("done!\n");
			return means;
		}
		printf("done!\n");
	}

	time->stopTimer(1);
	time->logIt();
	MFree_1D(shapeList);
	MFree_1D(dist);
	MFree_1D(N);
	MFree_1D(kMeansList);
	MFree_1D(meanCN);
	MFree_1D(coordCount);
	MFree_0D(distMeasures);
	MFree_0D(time);
	MFree_0D(pFunctions);
	MFree_0D(pSDP);
	MFree_1D(numOfPionts);
	MFree_1D(sumOfVector);
	MFree_1D(tempKMeanList);
	MFree_1D(wS);
	MFree_1D(wV);
	return means;
}


//====================================================================================================================
// Means of Signatures in a same Clusterer using simple K-Means
// get: shapeList: Nr. of Signatures these belong to one same Clusterer
//			numOfShape: Num of Signatures these belong to one same Clusterer
//			countK: Number of Clusters to cluster these instants
//		  pData: all vectors in this pData these will be clustered
//			sumOfPixels: number of pixels in a vector of pData
//			sumOfvector: sum of the vector of pData
// return: kMeans of centroids
//====================================================================================================================
SC_Signature* SC_Clusterer::getMeanSignature(int* shapeList, int numOfShape, int countK, SV_Data* pData, int* numOfPoints, float* sumOfvector){
	SC_Signature *pSignature = NULL;
	SC_SDP *pSDP = new SC_SDP(this->pTweak->featureSdp.m, this->pTweak->featureSdp.lag, this->pTweak->featureSdp.color, this->pTweak->featureSdp.n, this->pTweak->featureSdp.pictureSize, this->pTweak->featureSdp.tau);
	int dimOfVector = pData->Col;
	
	if (numOfShape == 1) {
		pSignature = pSDP->vector2signature(pData->Mat[shapeList[0]], dimOfVector);
		return pSignature;
	} else {
		SC_Clusterer::Pixel **pPixels, *pPixelsOfMaxShape;
		int x, y, k;
		float *pMaxShape;
		int maxPixel, index, dim;
		double *sumOfminiDist, *pTempPoint, *pMeanWeigths, *pMeanCoodinates, sumOfWeights, factor = 0.0;
		int *indices;
		float point[2];
		SC_Centroid_Point **pMeanPoints;
		SC_MatrixFunctions *pFunctions = new SC_MatrixFunctions();
		int squareSize = (int)(sqrt((double)dimOfVector));
		
		//sort the shapeList with num of points
		sortShapeList(shapeList, 0, numOfShape-1, numOfPoints);

		pMaxShape = pData->Mat[shapeList[numOfShape-1]];
		maxPixel = numOfPoints[shapeList[numOfShape-1]];

		MArray_1D(pPixelsOfMaxShape, maxPixel, SC_Clusterer::Pixel,"SC_Clusterer::getMeanSignature.pPixelsOfMaxShape");
		MArray_2D(pPixels, maxPixel, numOfShape-1, SC_Clusterer::Pixel,"SC_Clusterer::getMeanSignature.pixels");
		MArray_1D(indices, maxPixel, int, "SC_Clusterer::getMeanSignature.indices");
		sumOfminiDist = pFunctions->initVector(maxPixel, 0.0);


		//col these dis between every pixel in maxShap and each other
		int j = 0;
		for(x = 0; x < dimOfVector; x++){
			if(pMaxShape[x] > 0){
				pPixelsOfMaxShape[j].setIndexInVector(x);
				pPixelsOfMaxShape[j].setDim(2);
				pPixelsOfMaxShape[j].setWeight(pMaxShape[x] / sumOfvector[shapeList[numOfShape-1]]);
				point[0] = (float) (x / squareSize);
				point[1] = (float) (x % squareSize);
				for(y = 0; y < numOfShape-1; y++){
					float *vector = pData->Mat[shapeList[y]];
					pPixels[j][y] = getMiniDistPixel(point, 2, vector, sumOfvector[shapeList[y]], dimOfVector);
					sumOfminiDist[j] += pPixels[j][y].getDist();
				}
				j++;
			}
		}

		for(x = 0;  x < maxPixel; x++){
			indices[x] = x;
		}

		//sort the sumOfMinDist
		sortIndex(sumOfminiDist, indices, 0, maxPixel-1);
		MFree_1D(sumOfminiDist);
		
		dim = 2;
		MArray_1D(pMeanPoints,countK, SC_Centroid_Point*, "SC_Clusterer::meanSignature.meanPoint");
		MArray_1D(pMeanWeigths,countK, double, "SC_Clusterer::meanSignature.pMeanWeigths");
		MArray_1D(pTempPoint, dim, double, "SC_Clusterer::meanSignature.pTempPoint");

		//init with zeros
		for (k = 0; k < countK; k++) {
			for (int d = 0; d < dim; d++) {
				pTempPoint[d] = 0.0;
			}
			pMeanPoints[k] = new SC_Centroid_Point(this->pTweak, dim, pTempPoint, false);
			pMeanWeigths[k] = 0.0; 
		}
		
		//select the frist countK indices of pixel and col the mean
		j = 0; 
		double *pCurrentMeanCoodinates = pFunctions->initVector(dim, 0.0);
		double currentMeanWeights;
		double weight = 0.0; 
		for(int i = 0; i < countK; i++){
			do{
					index = indices[j];
					
					for(int d = 0; d < dim; d++){
						pCurrentMeanCoodinates[d] = 0.0;
					}
					currentMeanWeights = 0.0;
					
					sumOfWeights = pPixelsOfMaxShape[index].getWeight();
					for(y = 0; y < numOfShape - 1; y++){
						weight = pPixels[index][y].getWeight();
						sumOfWeights += weight;
					}
					
					point[0] = (float)(pPixelsOfMaxShape[index].getIndexInVector() / squareSize);
					point[1] = (float)(pPixelsOfMaxShape[index].getIndexInVector() % squareSize);
					weight = pPixelsOfMaxShape[index].getWeight();
					for(int d = 0; d < dim; d++){
						pCurrentMeanCoodinates[d] +=  point[d] * weight * (1 / sumOfWeights);
					}

					currentMeanWeights += weight;

					for(y = 0; y < numOfShape - 1; y++){
						int indexInVector = pPixels[index][y].getIndexInVector();
						point[0] = (float)(indexInVector / squareSize);
						point[1] = (float)(indexInVector % squareSize);
						weight = pPixels[index][y].getWeight();
						for(int d = 0; d < dim; d++){
							pCurrentMeanCoodinates[d] +=  point[d] * weight * (1 / sumOfWeights);
						}
						currentMeanWeights += weight;
						currentMeanWeights /= numOfShape;
					}
					j++;
			}while((i > 0) && isContain(pCurrentMeanCoodinates, pMeanPoints, dim, i)&& ((countK - i) < maxPixel - j));
			
			pMeanCoodinates = pMeanPoints[i]->getCoordinate();
			for(int d = 0; d < dim; d++){
				pMeanCoodinates[d] = pCurrentMeanCoodinates[d];
			}
			pMeanWeigths[i] = currentMeanWeights;
		}

		//normalize the weights
		sumOfWeights = 0.0;
		for(int i = 0; i < countK; i++){
			sumOfWeights += pMeanWeigths[i];
		}

		for(int i = 0; i < countK; i++){
			pMeanWeigths[i] /= sumOfWeights ;
		}
		pSignature = new SC_Signature((SC_Centroid**)pMeanPoints, pMeanWeigths, countK, true, 1.0);
		pSignature->setJustLink(false, false); //this way, the pMeanPoints get killed when the signature gets killed but we create no copy here

		MFree_2D(pPixels);
		MFree_1D(pPixelsOfMaxShape);
		MFree_0D(pFunctions);
		MFree_1D(pTempPoint);
		MFree_1D(indices);
		MFree_1D(pCurrentMeanCoodinates);

		return pSignature;
	}
}

//====================================================================================================================
// Initializes kMeans randomly, returning the found means
//====================================================================================================================
SV_Data* SC_Clusterer::randomInitialization(unsigned int k, SV_Data *pData) {
	unsigned int t, d, i, r, idx;
	SV_Data *pMeans;
	vector<long int> indexes;
	bool unique, allEqual;

	//init with zeros
	pMeans = new SV_Data(k, pData->Col);
	for (i = 0; i < k; i++) {
		for (d = 0; d < (unsigned int)(pData->Col); d++) {
			pMeans->Mat[i][d] = 0.0f;
		}
	}

	//assign each feature-vector randomly (or heuristically) to a class
	for (t = 0; t < (unsigned int)(pData->Row); t++) {
		indexes.push_back(t);
	}

	//find k distinct feature medoids
	for (i = 0; i < k; i++) {
		do {
			//draw the index of the next medoid (centroid) at random and delete it from the repository of available indexes
			unique = true;
			r = sclib::rand((unsigned int)(indexes.size())-1);
			idx = indexes[r];
			indexes.erase(indexes.begin() + r);

			//check if the corresponding vector is unique
			allEqual = false;
			for (t = 0; t < i; t++) {
				allEqual = true;
				for (d = 0; d < (unsigned int)(pData->Col); d++) {
					if (pMeans->Mat[t][d] != pData->Mat[idx][d]) {
						allEqual = false;
						break;
					}
				}
				if (allEqual == true) {
					break;
				}
			}
			if (allEqual == true) {
				unique = false;
			}

			if (unique==false && indexes.empty()==true) {
				break;
			}
		} while (unique == false);

		if (unique==false && indexes.empty()==true) { //proceed with less means
			printf("Error: Too big k in k-means for this dataset, can't find %d different clusters, using %d!", k, i);
			SV_Data *pNewMeans = new SV_Data(i, pMeans->Col);
			for (unsigned int x = 0; x < i; x++) {
				for (d = 0; d < (unsigned int)(pMeans->Col); d++) {
					pNewMeans->Mat[x][d] = pMeans->Mat[x][d];
				}
			}
			MFree_0D(pMeans);
			pMeans = pNewMeans;
			break;
		}

		//copy that point as the next centroid
		for (d = 0; d < (unsigned int)(pData->Col); d++) {
			pMeans->Mat[i][d] = pData->Mat[idx][d];
		}
	}

	return pMeans;
}

//====================================================================================================================
// Initializes kMeans deterministically employing Var-Part algorithm from "In Search of Deterministic Methods for 
// Initializing k-Means and Gaussian Mixture Clustering", Su, Dy, 2006
// The method returns the found means
//====================================================================================================================
SV_Data* SC_Clusterer::deterministicInitialization(unsigned int k, SV_Data *pData) {
	SC_MatrixFunctions matFunc;
	SV_Data *pMeans;
	double *variance, *mean, *SSE;
	unsigned int currentK=1, cj=0, dp, t, d, j;
	float **pCjIdx;
	unsigned int *mapping, *vectorCount;
	
	pMeans = new SV_Data(k, pData->Col);

	if (k == 1) {
		mean = matFunc.mean(pData->Mat, pData->Row, pData->Col);
		for (d = 0; d < (unsigned int)(pData->Col); d++) {
			pMeans->Mat[0][d] = (float)(mean[d]);
		}
		MFree_1D(mean);
	} else {
		mapping = matFunc.initVector(pData->Row, (unsigned int)(0));
		SSE = matFunc.initVector(pMeans->Row, 0.0);
		vectorCount = matFunc.initVector(pMeans->Row, (unsigned int)(0));
		vectorCount[cj] = pData->Row;
		
		while (currentK < k) { //iteratively split the cluster Cj that has maximum SSE (squared error) measure; split it perpendicular to the dimension of largest variance
			if (currentK > 1) {
				//build a memory-saving matrix representation of all vectors belonging to Cj
				MArray_1D(pCjIdx, vectorCount[cj], float*, "SC_Clusterer.varPart: pCjIdx");
				j = 0;
				for (t = 0; t < (unsigned int)(pData->Row); t++) {
					if (mapping[t] == cj) {
						pCjIdx[j++] = pData->Mat[t];
					}
				}
			} else {
				pCjIdx = pData->Mat;
			}

			//find dimension to split on
			mean = matFunc.mean(pCjIdx, vectorCount[cj], pData->Col); //mean in the choosen cluster to plit
			variance = matFunc.variance(pCjIdx, vectorCount[cj], pData->Col, mean); //variance in the choosen cluster to split
			dp = matFunc.maxIdx(variance, pData->Col); //dimension with largest variance
			MFree_1D(variance);
			if (pCjIdx != pData->Mat) {
				MFree_1D(pCjIdx);
			}

			//re-map vectors to the two new clusters and re-estimate means
			vectorCount[cj] = 0; //the "lower" part of the splitted cluster will keep the original cluster's id
			vectorCount[currentK] = 0; //the "upper" part will get the next free one
			for (d = 0; d < (unsigned int)(pMeans->Col); d++) {
				pMeans->Mat[cj][d] = 0.0f;
				pMeans->Mat[currentK][d] = 0.0f;
			}
			for (t = 0; t < (unsigned int)(pData->Row); t++) {
				if (mapping[t] == cj) {
					if (pData->Mat[t][dp] <= mean[dp]) {
						mapping[t] = cj;
						vectorCount[cj]++;
						for (d = 0; d < (unsigned int)(pMeans->Col); d++) {
							pMeans->Mat[cj][d] += pData->Mat[t][d];
						}
					} else {
						mapping[t] = currentK;
						vectorCount[currentK]++;
						for (d = 0; d < (unsigned int)(pMeans->Col); d++) {
							pMeans->Mat[currentK][d] += pData->Mat[t][d];
						}
					}
				}
			}
			MFree_1D(mean);
			for (d = 0; d < (unsigned int)(pMeans->Col); d++) { //renormalize means
				pMeans->Mat[cj][d] /= (float)(vectorCount[cj]);
				pMeans->Mat[currentK][d] /= (float)(vectorCount[currentK]);
			}

			//check if there are vectors in the new cluster otherwise complain and proceed with less means
			if (vectorCount[currentK] < 2) {
				printf("Error: Too big k in k-means for this dataset, can't find %d different clusters, using %d!", k, currentK);
				SV_Data *pNewMeans = new SV_Data(currentK, pMeans->Col);
				for (unsigned int x = 0; x < currentK; x++) {
					for (d = 0; d < (unsigned int)(pMeans->Col); d++) {
						pNewMeans->Mat[x][d] = pMeans->Mat[x][d];
					}
				}
				MFree_0D(pMeans);
				pMeans = pNewMeans;
				break;
			}

			//re-estimate error
			SSE[cj] = 0.0;
			SSE[currentK] = 0.0;
			for (t = 0; t < (unsigned int)(pData->Row); t++) {
				if (mapping[t] == cj) {
					SSE[cj] += SC_DistanceMeasures::squaredError(pData->Mat[t], pMeans->Mat[cj], pMeans->Col);
				} else if (mapping[t] == currentK) {
					SSE[currentK] += SC_DistanceMeasures::squaredError(pData->Mat[t], pMeans->Mat[currentK], pMeans->Col);
				}
			}

			//find next cluster to split
			cj = matFunc.maxIdx(SSE, pMeans->Row);

			currentK++;
		}

		MFree_1D(mapping);
		MFree_1D(SSE);
		MFree_1D(vectorCount);
	}

	return pMeans;
}

//====================================================================================================================
// Takes the data in pData and returns the k Centroids determined using #iterations of k-means with Euclidean distance
//====================================================================================================================
SV_Data* SC_Clusterer::kMeans(unsigned int k, SV_Data *pData, unsigned int iterations, bool randomInit, unsigned int repeatedInitializations, bool replaceWithMedoids) {
	unsigned short int *clusterMapping = NULL;
	unsigned long int *N = NULL;
	SV_Data *pMeans = kMeans(k, pData, clusterMapping, N, iterations, randomInit, repeatedInitializations);

	MFree_1D(clusterMapping);
	MFree_1D(N);

	return pMeans;
}

SV_Data* SC_Clusterer::kMeans(unsigned int k, SV_Data *pData, unsigned short int* &clusterMapping, unsigned long int* &N, unsigned int iterations, bool randomInit, unsigned int repeatedInitializations, bool replaceWithMedoids) {
  unsigned long int t, x, T = pData->Row, *bestN;
	unsigned short int i, d, *newMapping, D = pData->Col, *bestMapping;
	unsigned int y;
	double dist, minDist;
	SV_Data *pMeans, *pBestMeans; 
	double lastPercentage = 0.0, overallDist, bestOverallDist = std::numeric_limits<double>::max();
	bool changed;
	vector<long int> indexes;

	if (this->verbose == true) {
		printf("\n  Starting kMeans...");
	}

	MFree_0D(clusterMapping);
	MArray_1D(clusterMapping, T, unsigned short int, "SC_Clusterer.kMean: clusterMapping");
	MArray_1D(newMapping, T, unsigned short int, "SC_Clusterer.kMean: newMapping");
	MArray_1D(bestMapping, T, unsigned short int, "SC_Clusterer.kMean: bestMapping");
	MFree_1D(N);
	MArray_1D(N, k, unsigned long int, "SC_Clusterer.kMean: N");
	MArray_1D(bestN, k, unsigned long int, "SC_Clusterer.kMean: bestN");

	if (randomInit == false) {
		repeatedInitializations = 1;
	}

	//repeat the whole process several times when random init is wished, return best result in terms of minimum overall distance of points to their centroids
	for (y = 0; y < repeatedInitializations; y++) {
		if (this->verbose == true) {
			printf("\n  kMeans initialization %d/%d: ", y+1, repeatedInitializations);
		}
		
		//initialize
		for (i = 0; i < k; i++) {
			N[i] = 0;
		}
		pMeans = (randomInit==true) ? randomInitialization(k, pData) : deterministicInitialization(k, pData);
		k = pMeans->Row; //may have changed (decreased) during init. if this dataset doesn't justify that many clusters

		//#iterations of k-means to update means
		for (x = 0; x < iterations; x++) {
			if (this->verbose == true) {
				printf("\n    kMeans iteration %d/%d: ", x+1, iterations);
			}

			//find the closest centroid (mean) for each feature-vector, remember it's cluster-nr in the newMapping
			overallDist = 0.0;
			changed = false;
			for (i = 0; i < k; i++) {
				N[i] = 0;
			}
			for (t = 0; t < T; t++) {
				minDist = std::numeric_limits<double>::max();
				for (i = 0; i < k; i++) {
					dist = SC_DistanceMeasures::squaredError(pData->Mat[t], pMeans->Mat[i], D);
					if (dist < minDist) {
						newMapping[t] = i;
						minDist = dist;
					}
				}
				N[newMapping[t]]++;
				overallDist += minDist;
				if (newMapping[t]!=clusterMapping[t] || x==0) {
					changed = true;
				}
			}

			//stop current initialization's run if no change in mapping occured
			if (changed == false) {
				if (this->verbose == true) {
					printf(" (overall error: %f; no further changes)", overallDist);
				}
				break;
			}
	
			//reestimate centroids
			for (i = 0; i < k; i++) {
				for (d = 0; d < D; d++) {
					pMeans->Mat[i][d] = 0.0f;
				}
			}
			for (t = 0; t < T; t++) {
				for (d = 0; d < D; d++) {
					pMeans->Mat[newMapping[t]][d] += (float)((1.0/(double)(N[newMapping[t]])) * (double)(pData->Mat[t][d]));
				}
				clusterMapping[t] = newMapping[t];
			} 

			if (this->verbose == true) {
				printf(" (overall error: %f)", overallDist);
			}
		} //x=0..kMeansCount

		if (overallDist < bestOverallDist) {
			bestOverallDist = overallDist;
			for (i = 0; i < k; i++) {
				bestN[i] = N[i];
			}
			for (t = 0; t < T; t++) {
				bestMapping[t] = clusterMapping[t];
			}
			pBestMeans = new SV_Data(*pMeans, false);
		}
		MFree_0D(pMeans);
	} //y=0..repeatedInitializations

	//remove emtpy clusters from pMeans
	y = 0;
	for (i = 0; i < k; i++) {
		if (bestN[i] > 0) {
			y++;
		}
	}
	pMeans = new SV_Data(y, D);
	y = 0;
	for (i = 0; i < k; i++) {
		if (bestN[i] > 0) {
			for (d = 0; d < D; d++) {
				pMeans->Mat[y][d] = pBestMeans->Mat[i][d];
			}
			if (y < i) { //remap vectors to new cluster number
				for (t = 0; t < T; t++) {
					if (bestMapping[t] == i) {
						bestMapping[t] = y;
					}
				}
			}
			y++;
		}
	}
	MFree_1D(N);
	N = bestN;
	MFree_1D(clusterMapping);
	clusterMapping = bestMapping;

	//replace each centroid with the nearest medoid, if wished
	if (replaceWithMedoids = true) {
		for (i = 0; i < y; i++) {
			minDist = std::numeric_limits<double>::max();
			for (t = 0; t < T; t++) {
				dist = SC_DistanceMeasures::squaredError(pMeans->Mat[i], pData->Mat[t], D);
				if (dist < minDist) {
					minDist = dist;
					x = t;					
				}
			}
			for (d = 0; d < D; d++) {
				pMeans->Mat[i][d] = pData->Mat[x][d];
			}
		}
	}

	//clean up
	MFree_0D(pBestMeans);
	MFree_1D(newMapping);

	if (this->verbose == true) {
		printf("\n  kMeans (%d/%d found) finished. Best overall error: %f", y, k, bestOverallDist);
	}

	return pMeans;
}

//====================================================================================================================
// Gets a set of "templates" in pMeans and replaces each vector in pData with the nearest template vector according
// to the Euclidean distance; better search strategy than brute force: use ANN (approximate nearest neighbor search)
//====================================================================================================================
void SC_Clusterer::pickNearest(SV_Data *pData, SV_Data *pMeans) {
	double minDist, dist;
	int minIdx, t, i, d;

	if (pData->Col != pMeans->Col) {
		REPORT_ERROR(SVLIB_BadArg, "pickNearest needs similar vector dimensions for test- and template-data");
	} else {
		for (t = 0; t < pData->Row; t++) {
			minDist = std::numeric_limits<double>::max();
			for (i = 0; i < pMeans->Row; i++) {
				dist = SC_DistanceMeasures::euclid(pData->Mat[t], pMeans->Mat[i], pData->Col);
				if (dist < minDist) {
					minDist = dist;
					minIdx = i;
				} //if min dist
			} //for i
			for (d = 0; d < pData->Col; d++) {
				pData->Mat[t][d] = pMeans->Mat[minIdx][d];
			} //for d
		} //for t
	} //if #columns equals

	return;
}
