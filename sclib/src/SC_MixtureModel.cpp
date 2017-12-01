/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel as a base-class for new mixture models				  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.04.2004																								*/
/**************************************************************************/

#include <math.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <iomanip>
#include <string.h>
#include "SC_MixtureModel.h"
#include "SC_TweakableParameters.h"
#include "SC_Clusterer.h"

SC_Gauss SC_MixtureModel::gaussSolver(20.0, 0.001);

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel::SC_MixtureModel(SC_TweakableParameters *pTweak) : SC_Model(pTweak), dummyWeight(1.0), dummyMean(-1.0E10), dummyVariance(1.0E-10) {
  this->mean = NULL;
  this->variance = NULL;
  this->sd = NULL;
  this->weight = NULL;
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_MixtureModel::SC_MixtureModel(const SC_MixtureModel& pParent) : SC_Model(pParent), dummyWeight(1.0), dummyMean(-1.0E10), dummyVariance(1.0E-10) {
  this->mixtureCount = pParent.mixtureCount;
	this->dim = pParent.dim;
  this->maxEMsteps = pParent.maxEMsteps;
	
	this->weight = NULL;
	this->variance = NULL;
	this->sd = NULL;
	this->mean = NULL;

	if (this->mixtureCount > 0) {
		if (pParent.weight != NULL) {
			MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel.SC_MixtureModel: this->weight");
		}
		if (pParent.variance != NULL) {
			MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel.SC_MixtureModel: this->variance");
		}
		if (pParent.sd != NULL) {
			MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel.SC_MixtureModel: this->sd");
		}
		if (pParent.mean != NULL) {
			MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel.SC_MixtureModel: this->mean");
		}
		
		for (unsigned short int m = 0; m < this->mixtureCount; m++) {
			if (this->weight != NULL) {
				this->weight[m] = pParent.weight[m];
			}

			if (this->dim > 0) {
				for (unsigned short int d = 0; d < this->dim; d++) {
					if (this->variance != NULL) {
						this->variance[m][d] = pParent.variance[m][d];
					}
					if (this->sd != NULL) {
						this->sd[m][d] = pParent.sd[m][d];
					}
					if (this->mean != NULL) {
						this->mean[m][d] = pParent.mean[m][d];
					}
				}
			}
		}
	}
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_MixtureModel::~SC_MixtureModel() {
  MFree_2D(this->mean);
  MFree_2D(this->variance);
  MFree_2D(this->sd);
  MFree_1D(this->weight);
}

//====================================================================================================================
// assignment operator
//====================================================================================================================
SC_MixtureModel& SC_MixtureModel::operator=(const SC_MixtureModel& pParent) {
	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		this->mixtureCount = pParent.mixtureCount;
		this->dim = pParent.dim;
		this->maxEMsteps = pParent.maxEMsteps;
		
		MFree_1D(this->weight);
		MFree_2D(this->variance);
		MFree_2D(this->sd);
		MFree_2D(this->mean);

		if (this->mixtureCount > 0) {
			if (pParent.weight != NULL) {
				MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel.operator=: this->weight");
			}
			if (pParent.variance != NULL) {
				MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel.operator=: this->variance");
			}
			if (pParent.sd != NULL) {
				MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel.operator=: this->sd");
			}
			if (pParent.mean != NULL) {
				MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel.operator=: this->mean");
			}
			
			for (unsigned short int m = 0; m < this->mixtureCount; m++) {
				if (this->weight != NULL) {
					this->weight[m] = pParent.weight[m];
				}

				if (this->dim > 0) {
					for (unsigned short int d = 0; d < this->dim; d++) {
						if (this->variance != NULL) {
							this->variance[m][d] = pParent.variance[m][d];
						}
						if (this->sd != NULL) {
							this->sd[m][d] = pParent.sd[m][d];
						}
						if (this->mean != NULL) {
							this->mean[m][d] = pParent.mean[m][d];
						}
					}
				}
			}
		}
	}

	return *this;
}

//====================================================================================================================
// model-initialization prior to parameter-estimation is neccessary, but not understood well. as described by
// rose/reynolds, i here use random means, followed by kMeansCount iterations of k-means to cluster the data somehow
// meaningfull into the mixtures, which  should provide broad acoustic classes.
//
// this method also allocates memory for the mean- and variance-vectors
//====================================================================================================================
void SC_MixtureModel::initParameters(SV_Data *pData, double varianceLimit, unsigned int kMeansCount) {
	int t, i, d;
	unsigned short int *clusterMapping = NULL;
	unsigned long int *N = NULL;
	SC_Clusterer clusterer(this->pTweak, false);
	SV_Data *pMeans = clusterer.kMeans(this->mixtureCount, pData, clusterMapping, N, kMeansCount, false, 5, false);

	if (pMeans->Row < this->mixtureCount) {
		printf("Error: Mixture count probably too high, reducing from %d to %d!", this->mixtureCount, pMeans->Row);
		this->mixtureCount = pMeans->Row;
		MFree_2D(this->variance);
		MFree_2D(this->sd);
		MFree_2D(this->mean);
		MFree_1D(this->weight);
    MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: sd");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_GMM: mean");
    MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_GMM: weight");
	}

	//copy the means to internal storage, initialize variances (where has this been done before???)
	for (i = 0; i < this->mixtureCount; i++) {
		for (d = 0; d < this->dim; d++) {
			this->mean[i][d] = pMeans->Mat[i][d];
			this->variance[i][d] = 0.0;
		}
	}
	MFree_0D(pMeans);

	//1. step to estimate variances
	for (t = 0; t < pData->Row; t++) {
		for (d = 0; d < this->dim; d++) {
			this->variance[clusterMapping[t]][d] += (pData->Mat[t][d] - this->mean[clusterMapping[t]][d]) * (pData->Mat[t][d] - this->mean[clusterMapping[t]][d]);
		}
	}

	//estimate true weights, 2. step to estimate variances, variance-limiting
	for (i = 0; i < this->mixtureCount; i++) {
		this->weight[i] = (double)(N[i]) / (double)(pData->Row); //to fulfill the constraint SUM_i=1..M weight[i] = 1
		for (d = 0; d < this->dim; d++) {
			this->variance[i][d] /= (double)(N[i] - 1);
      if ((this->variance[i][d] < varianceLimit) || (!sclib::isFinite(this->variance[i][d]))) { //do variance-limiting
        this->variance[i][d] = varianceLimit;
      }
      this->sd[i][d] = sqrt(this->variance[i][d]);
		}
	}

	MFree_1D(clusterMapping);
	MFree_1D(N);
	
	sortParameters();

	return;
/*
	unsigned long int t, x, T = pData->Row;
	unsigned short int i, d, r, minDist, *kMeansList;
	double *dist;
	unsigned long int *N; //contains the count of feature-vectors which formed a mixture

	//srand((unsigned int)time(NULL));
	MArray_1D(dist, this->mixtureCount, double, "SC_MixtureModel.initParameters: dist");
	MArray_1D(kMeansList, pData->Row, unsigned short int, "SC_MixtureModel.initParameters: kMeansList");
	MArray_1D(N, this->mixtureCount, unsigned long int, "SC_MixtureModel.initParameters: N");

	//init with zeros
	for (i = 0; i < this->mixtureCount; i++) {
		N[i] = 0;
		for (d = 0; d < this->dim; d++) {
			this->mean[i][d] = 0.0;
			this->variance[i][d] = 0.0;
		}
	}

	//assign each feature-vector randomly to a class mean, remember this mapping in kMeansList
	for (t = 0; t < T; t++) {
		r = (unsigned short)(t % this->mixtureCount);
		//r = sclib::rand(this->mixtureCount-1); //(unsigned short int)(floor(((float)(rand()) / (float)(RAND_MAX+1)) * this->mixtureCount)); //generate a random number between 0..M-1, M=mixtureCount
		N[r]++; //count the number of vectors used to estimate this mixture; 
		kMeansList[t] = r; //remember to which mixture this vector was assigned
		for (d = 0; d < this->dim; d++) {
			this->mean[r][d] += pData->Mat[t][d];
		}
	}

	//normalize the estimated mean
	for (i = 0; i < this->mixtureCount; i++) {
		for (d = 0; d < this->dim; d++) {
			this->mean[i][d] /= N[i];
		}
	}

	//kMeansCount iterations of k-means to update means
	for (x = 0; x < kMeansCount; x++) {
		for (t = 0; t < T; t++) {
			
			//find the closest centroid (mean) for this feature-vector, remember it's mixture-nr in "minDist"
			minDist = 0;
			for (i = 0; i < this->mixtureCount; i++) {
				dist[i] = 0.0;
				for (d = 0; d < this->dim; d++) {
					dist[i] += (pData->Mat[t][d] - this->mean[i][d]) * (pData->Mat[t][d] - this->mean[i][d]); //dist[i] = SQRT(SUM_d=0..D( (x-µ)^2 ))
				}
				dist[i] = sqrt(dist[i]);
				if (dist[i] < dist[minDist]) {minDist = i;}
			}

			//move the point from old centroid to new centroid
			if (N[kMeansList[t]] > 2) {
				for (d = 0; d < this->dim; d++) {
					this->mean[kMeansList[t]][d] = (1.0/(N[kMeansList[t]]-1.0)) * ((N[kMeansList[t]]*this->mean[kMeansList[t]][d]) - pData->Mat[t][d]); //remove this point from old centroid
					this->mean[minDist][d] = (1.0/(N[minDist]+1.0)) * ((N[minDist]*this->mean[minDist][d]) + pData->Mat[t][d]); //add it to new one
				}
				N[kMeansList[t]]--; //there is now one member less in this mixture-mean
				N[minDist]++; //and this mixture has got one more
				kMeansList[t] = minDist; //note that this point is now a member of the new mixture-mean				
			}

		}  //t=0..T
	} //x=0..kMeansCount

	//1. step to estimate variances
	for (t = 0; t < T; t++) {
		for (d = 0; d < this->dim; d++) {
			this->variance[kMeansList[t]][d] += (pData->Mat[t][d] - this->mean[kMeansList[t]][d]) * (pData->Mat[t][d] - this->mean[kMeansList[t]][d]);
		}
	}

	//estimate true weights, 2. step to estimate variances, variance-limiting
	for (i = 0; i < this->mixtureCount; i++) {
		this->weight[i] = (double)(N[i]) / (double)(pData->Row); //to fulfill the constraint SUM_i=1..M weight[i] = 1
		for (d = 0; d < this->dim; d++) {
			this->variance[i][d] /= (double)(N[i] - 1);
      if ((this->variance[i][d] < varianceLimit) || (!sclib::isFinite(this->variance[i][d]))) { //do variance-limiting
        this->variance[i][d] = varianceLimit;
      }
      this->sd[i][d] = sqrt(this->variance[i][d]);
		}
	}
	
	sortParameters();

	MFree_1D(dist);
	MFree_1D(kMeansList);
	MFree_1D(N);

	return;
*/
}

//====================================================================================================================
// Sort the model-parameters according to the weights
//====================================================================================================================
void SC_MixtureModel::sortParameters() {
	unsigned short int i, d;
	SV_Data *temp;
	double *oldWeight, **oldMean, **oldVariance, **oldSd;
	
	temp = new SV_Data(this->mixtureCount, 2);
	MArray_1D(oldWeight, this->mixtureCount, double, "SC_MixtureModel.sortParameters: oldWeight");
	MArray_2D(oldMean, this->mixtureCount, this->dim, double, "SC_MixtureModel.sortParameters: oldMean");
	MArray_2D(oldVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel.sortParameters: oldVariance");
  MArray_2D(oldSd, this->mixtureCount, this->dim, double, "SC_MixtureModel.sortParameters: oldSd");

	for (i = 0; i < this->mixtureCount; i++) {
		temp->Mat[i][0] = (float)(this->weight[i]);
		temp->Mat[i][1] = (float)(i); //remember the old index in the 2nd column of the data-matrix
		oldWeight[i] = this->weight[i];
		for (d = 0; d < this->dim; d++) {
			oldMean[i][d] = this->mean[i][d];
			oldVariance[i][d] = this->variance[i][d];
      oldSd[i][d] = this->sd[i][d];
		}
	}

	//sort the list of weights; remember the old index in the 2nd column of the data-matrix
	//sclib::quickSort(tempt, 0, this->mixtureCount - 1, 2, 0);
  sclib::quickSort(temp->Mat, 0, this->mixtureCount-1, 2, 0);

	for (i = 0; i < this->mixtureCount; i++) {
		this->weight[i] = oldWeight[(int)(temp->Mat[i][1])];
		for (d = 0; d < this->dim; d++) {
			this->mean[i][d] = oldMean[(int)(temp->Mat[i][1])][d];
			this->variance[i][d] = oldVariance[(int)(temp->Mat[i][1])][d];
      this->sd[i][d] = oldSd[(int)(temp->Mat[i][1])][d];
		}
	}

	MFree_0D(temp);
	MFree_1D(oldWeight);
	MFree_2D(oldMean);
	MFree_2D(oldVariance);
  MFree_2D(oldSd);

	return;
}

//====================================================================================================================
// Remove the indexed mixture-component from the models mean/variance/weight-vectors/matrices
//====================================================================================================================
void SC_MixtureModel::killMixture(unsigned short int mixture) {
  unsigned short int i, i_idx = 0, d,  I = this->mixtureCount, D = this->dim;
  double *newWeight, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I-1, double, "SC_MixtureModel.killMixture: newWeight");
  MArray_2D(newVariance, I-1, D, double, "SC_MixtureModel.killMixture: newVariance");
  MArray_2D(newMean, I-1, D, double, "SC_MixtureModel.killMixture: newMean");
  MArray_2D(newSd, I-1, D, double, "SC_MixtureModel.killMixture: newSd");

  for (i = 0; i < I; i++) { //copy mixtures, leave the one to kill out
    if (i != mixture) {
      newWeight[i_idx] = this->weight[i] * (1.0  / (1.0 - this->weight[mixture])); //renormalize weight
      for (d = 0; d < D; d++) {
        newVariance[i_idx][d] = this->variance[i][d];
        newMean[i_idx][d] = this->mean[i][d];
        newSd[i_idx][d] = this->sd[i][d];
      }
			if (I == 2 && newWeight[i_idx] == 0.0) { //if components had no weight at all...
				newWeight[i_idx] = 1.0;
			}
      i_idx++;
    }
  }

  MFree_1D(this->weight);
  MFree_2D(this->variance);
  MFree_2D(this->mean);
  MFree_2D(this->sd);

  this->weight = newWeight;
  this->variance = newVariance;
  this->mean = newMean;
  this->sd = newSd;
  this->mixtureCount--;

  return;
}

//====================================================================================================================
// Add one new mixture component with all means=0, variances=varianceLimit, weight=1/mixtureCount
//====================================================================================================================
void SC_MixtureModel::addMixture(void) {
  unsigned short int i, d, I = this->mixtureCount, D = this->dim;
  double *newWeight, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I+1, double, "SC_MixtureModel.addMixture: newWeight");
  MArray_2D(newVariance, I+1, D, double, "SC_MixtureModel.addMixture: newVariance");
  MArray_2D(newMean, I+1, D, double, "SC_MixtureModel.addMixture: newMean");
  MArray_2D(newSd, I+1, D, double, "SC_MixtureModel.addMixture: newSd");

  //add new mixture
  newWeight[I] = 1.0 / (this->mixtureCount+1.0);
  for (d = 0; d < D; d++) {
    newVariance[I][d] = 1.0;
    newMean[I][d] = 0.0;
    newSd[I][d] = 1.0;
  }

  //recalc old mixtures
  for (i = 0; i < I; i++) { 
    newWeight[i] = this->weight[i] * (1.0 - newWeight[I]);
    for (d = 0; d < D; d++) {
      newVariance[i][d] = this->variance[i][d];
      newMean[i][d] = this->mean[i][d];
      newSd[i][d] = this->sd[i][d];
    }
  }

  MFree_1D(this->weight);
  MFree_2D(this->variance);
  MFree_2D(this->mean);
  MFree_2D(this->sd);

  this->weight = newWeight;
  this->variance = newVariance;
  this->mean = newMean;
  this->sd = newSd;
  this->mixtureCount++;

  return;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_MixtureModel::modelOut(ostream& os) {
	unsigned short int x, y;

	os.setf(ios::fixed|ios::basefield);
  os.precision(5);

	os << "Mixtures: " << this->mixtureCount << endl;
	os << "Dimension: " << this->dim << endl;
	
	//parameters in matlab syntax for simple visualization
	os << "w=[";
	for (x = 0; x < this->mixtureCount-1; x++) {
		os << this->weight[x] << " ";
	}
	os << this->weight[this->mixtureCount-1] << "];" << endl;
	os << "m=[";
	for (x = 0; x < this->mixtureCount; x++) {
		for (y = 0; y < this->dim-1; y++) {
			os << this->mean[x][y] << " ";
		}
		if (x < this->mixtureCount-1) {
			os << this->mean[x][this->dim-1] << ";";
		} else {
			os << this->mean[x][this->dim-1] << "];" << endl;
		}
	}
	os << "v=[";
	for (x = 0; x < this->mixtureCount; x++) {
		for (y = 0; y < this->dim-1; y++) {
			os << this->variance[x][y] << " ";
		}
		if (x < this->mixtureCount-1) {
			os << this->variance[x][this->dim-1] << ";";
		} else {
			os << this->variance[x][this->dim-1] << "];" << endl;
		}
	}

	//and again in human readable form
	os << "MixNr\tWeight\tMean" << endl;
	for (x = 0; x < this->mixtureCount; x++) {
		os << x << "\t" << setw(7) << this->weight[x] << "\t";
		for (y = 0; y < this->dim; y++) {
			os << setw(10) << this->mean[x][y] << " ";
		}
		os << endl;
	}
	os << endl << "MixNr\t\tVariances" << endl;
	for (x = 0; x < this->mixtureCount; x++) {
		os << x << "\t\t";
		for (y = 0; y < this->dim; y++) {
				os << setw(10) << this->variance[x][y] << " ";
		}	
		os << endl;
	}
	
	os << endl;

	return(os);
}

//====================================================================================================================
// Step trough a linked list of models (starting from this one) and return the greatest nr. of mixture-components,
// that one single model in the list holds 
//====================================================================================================================
unsigned int SC_MixtureModel::getMaxMixturesInList() {
	unsigned int	maxN	= 0;
	SC_MixtureModel*			pHook = this;

	while (pHook != NULL) {
		maxN	= (pHook->getMixtureCount() > maxN) ? pHook->getMixtureCount() : maxN;
		pHook	= (SC_MixtureModel*)(pHook->Next);
	}

	return maxN;
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes written, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_MixtureModel::SaveModel(void) {
 	int res, x, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel model Failed!");}

  //write mixtureCount and featureDim
	bytes = io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel Model Failed!");}

  //write weight-vector
	bytes += io.writeArray(&(this->DFile), this->weight, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel Model Failed!");}

	//write mean and variance matrices
	for (x = 0; x < this->mixtureCount; x++) {
		bytes += io.writeArray(&(this->DFile), this->mean[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->variance[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->sd[x], this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel Model Failed!");}
	}

  //write trainingDataCount
  bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel Model Failed!");}

  return bytes + MHLen; //MHLen is just an estimate...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model * SC_MixtureModel::LoadModel(void) {
	int res, x;
	unsigned short int newDim = 0, newMixtureCount = 0;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

	//read mixtureCount and featureDim
	io.readScalar(&(this->DFile), newDim, codeSizes, fileSizes);
	if ((DFile.good() != TRUE) || (newDim == 0)) {return(NULL);}

	io.readScalar(&(this->DFile), newMixtureCount, codeSizes, fileSizes);
	if ((DFile.good() != TRUE || (newMixtureCount == 0))) {return(NULL);}

	//prepare this model for changing it's parameters
	if ((newDim != this->dim) || (newMixtureCount != this->mixtureCount)) {
		MFree_1D(this->weight);
		MFree_2D(this->mean);
		MFree_2D(this->variance);
    MFree_2D(this->sd);

		this->dim = newDim;
		this->mixtureCount = newMixtureCount;

 		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel.loadModel: weight");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel.loadModel: mean");
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel.loadModel: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel.loadModel: sd");
	};

	//read weight, mean and variance
	io.readArray(&(this->DFile), this->weight, this->mixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel Failed!");}
	for (x = 0; x < this->mixtureCount; x++) {
		io.readArray(&(this->DFile), this->mean[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->variance[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->sd[x], this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel Failed!");}
	}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel Failed!");}

  this->maxEMsteps = 150;

	return(this);
}

//====================================================================================================================
// Load model's parameter from current opened model file (stored in the old-fashioned way)
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model * SC_MixtureModel::oldLoadModel(void) {
	int res, x, newDim = 0, newMixtureCount = 0;
	SV_DataIO::SV_DatatypeSizes sizes;

  //read header
  res = oldLoadHdr(sizes);
	if (res == SVLIB_Fail) {return(NULL);}

	//read mixtureCount and featureDim
	DFile.read((char*)(&newDim), sizeof(unsigned short int));
	if ((DFile.good() != TRUE) || (newDim == 0)) {return(NULL);}

  DFile.read((char*)(&newMixtureCount), sizeof(unsigned short int));
	if ((DFile.good() != TRUE || (newMixtureCount == 0))) {return(NULL);}

	//prepare this model for changing it's parameters
	if ((newDim != this->dim) || (newMixtureCount != this->mixtureCount)) {
		MFree_1D(this->weight);
		MFree_2D(this->mean);
		MFree_2D(this->variance);
    MFree_2D(this->sd);

		this->dim = newDim;
		this->mixtureCount = newMixtureCount;

 		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel.loadModel: weight");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel.loadModel: mean");
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel.loadModel: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel.loadModel: sd");
	};

	//read weight, mean and variance
	DFile.read((char*)(this->weight), this->mixtureCount * sizeof(double));
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel Failed!");}
	for (x = 0; x < this->mixtureCount; x++) {
		DFile.read((char*)(this->mean[x]), this->dim * sizeof(double));
		DFile.read((char*)(this->variance[x]), this->dim * sizeof(double));
		DFile.read((char*)(this->sd[x]), this->dim * sizeof(double));
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel Failed!");}
	}

  //read trainingDataCount
  DFile.read((char*)(&(this->trainingDataCount)), sizeof(unsigned long int));
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel Failed!");}

  this->maxEMsteps = 150;

	return(this);
}

//====================================================================================================================
// create a (linked, not copied) SC_Signature-view on this model (for distance computation )
//====================================================================================================================
SC_Signature* SC_MixtureModel::toSignature(void) {
  SC_Centroid **centroids = NULL;
  SC_Signature *pSignature;
  
  MArray_1D(centroids, this->mixtureCount, SC_Centroid*, "SC_MixtureModel.toSignature: centroids");

  for (unsigned short int i = 0; i < this->mixtureCount; i++)   {
		centroids[i] = new SC_Centroid_Gaussian(this->pTweak, this->dim, this->mean[i], this->variance[i]);
  }

	pSignature = new SC_Signature(centroids, this->weight, this->mixtureCount);

  return pSignature;
}

//====================================================================================================================
// destruct a signature created by this model's toSignature()-method
//====================================================================================================================
void SC_MixtureModel::killSignature(SC_Signature *pSignature) {
  MFree_0D(pSignature);

  return;
}

//====================================================================================================================
// To compute the sd (srqt(variance)) out of the given variance scalar
//====================================================================================================================
double SC_MixtureModel::variance2sd(double variance) {
  return sqrt(variance);
}

//====================================================================================================================
// To compute the sd (srqt(variance)) out of the given variance vector
//====================================================================================================================
double* SC_MixtureModel::variance2sd(double *variance, unsigned int dim) {
  double *sd = NULL;

  MArray_1D(sd, dim, double, "SC_MixtureModel.variance2sd: sd");

  for (unsigned int i = 0; i < dim; i++) {
    sd[i] = sqrt(variance[i]);
  }

  return sd;
}

//====================================================================================================================
// To compute the sd (srqt(variance)) out of the given variance matrix (1 variance vector of a component per row)
//====================================================================================================================
double** SC_MixtureModel::variance2sd(double **variance, unsigned int len, unsigned int dim) {
  double **sd = NULL;

  MArray_2D(sd, (long int)len, dim, double, "SC_MixtureModel.variance2sd: sd");

  for (unsigned int i = 0; i < len; i++) {
    for (unsigned int d = 0; d < dim; d++) {
      sd[i][d] = sqrt(variance[i][d]);
    }
  }

  return sd;
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
// this method follows the instructions in "A Monte-Carlo Method for Score Normalization in Automatic Speaker 
// Verification using Kullback-Leibler Distances" by Ben, Blouet, Bimbot on ICASSP 2002
//====================================================================================================================
SV_Data* SC_MixtureModel::drawSamplesFromDistribution(unsigned long int count) {
  unsigned long int t;
  unsigned short int i, d; 
  SV_Data *pSamples = NULL;

  if (this->mixtureCount > 0 && this->dim > 0) {
    pSamples = new SV_Data(count, this->dim);
    
    for (t = 0; t < count; t++) {
      //1. randomly draw a mixture-index respecting the mixture weights
			i = sclib::drawIndexFromDistribution(this->weight, this->mixtureCount);

      //2. generate a gaussian random variable according to the parameters of the i'th mixture selected above
      for (d = 0; d < this->dim; d++) {
				pSamples->Mat[t][d] = (float)(sclib::getRandomizer()->rand_gaus(this->mean[i][d], this->sd[i][d])); //this uses the Box-Muller algorithm inside as in the paper
      }
    }
  }

  return pSamples;
}
