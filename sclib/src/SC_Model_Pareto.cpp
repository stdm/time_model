/**************************************************************************/
/*    Intended to model a speaker's voice's probability density via the   */
/*    non-parametric Pareto Density Estimation technique by Alfred Ultsch */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.02.2006																								*/
/**************************************************************************/

#include <time.h>
#include "SC_Model_Pareto.h"
#include "SC_Centroid_Point.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Model_Pareto::SC_Model_Pareto(SC_TweakableParameters *pTweak, SC_Model_Pareto *pBackgound, bool useMarginalDistributions) : SC_Model(pTweak) {
  this->pBackground = pBackgound;
  this->pPDEsolver = new SC_PDE(this->pTweak);
	this->Hdr.ModelType = sclib::mtPareto;

  this->useMarginalDistributions = useMarginalDistributions;
  this->dim = 0;
  this->pDensities = NULL;
  this->paretoRadius = NULL;
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_Model_Pareto::SC_Model_Pareto(const SC_Model_Pareto& pParent) : SC_Model(pParent) {
  this->pBackground = pParent.pBackground;
  this->pPDEsolver = new SC_PDE(this->pTweak); //this class holds no parameters, so it can be created newly and doesn't need to be copied

  this->useMarginalDistributions = pParent.useMarginalDistributions;
  this->dim = pParent.dim;
  if (this->dim > 0) {
    MArray_1D(this->pDensities, this->dim, SV_Data*, "SC_Model_Pareto: pDensities");
    MArray_1D(this->paretoRadius, this->dim, double, "SC_Model_Pareto: paretoRadius");
    for (unsigned int i = 0; i < this->dim; i++) {
      this->pDensities[i] = new SV_Data((SV_Data&)*(pParent.pDensities[i]));
      this->paretoRadius[i] = pParent.paretoRadius[i];
    }
  } else {
    this->pDensities = NULL;
    this->paretoRadius = NULL;
  }
}

//====================================================================================================================
// destructor
//====================================================================================================================
SC_Model_Pareto::~SC_Model_Pareto() {
  for (unsigned int i = 0; i < this->dim; i++) {
    MFree_0D(this->pDensities[i]);
  }
  MFree_1D(this->pDensities);
  MFree_1D(this->paretoRadius);
  MFree_0D(this->pPDEsolver);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_Model_Pareto& SC_Model_Pareto::operator=(const SC_Model_Pareto& pParent) {
  unsigned int i;

	if (this != &pParent) {
		this->SC_Model::operator=(pParent);

		//Destruct old parameters
		for (i = 0; i < this->dim; i++) {
			MFree_0D(this->pDensities[i]);
		}
		MFree_1D(this->pDensities);
		MFree_1D(this->paretoRadius);

		this->pBackground = pParent.pBackground;

		this->useMarginalDistributions = pParent.useMarginalDistributions;
		this->dim = pParent.dim;
		if (this->dim > 0) {
			MArray_1D(this->pDensities, this->dim, SV_Data*, "SC_Model_Pareto: pDensities");
			MArray_1D(this->paretoRadius, this->dim, double, "SC_Model_Pareto: paretoRadius");
			for (i = 0; i < this->dim; i++) {
				this->pDensities[i] = new SV_Data((SV_Data&)*(pParent.pDensities[i]));
				this->paretoRadius[i] = pParent.paretoRadius[i];
			}
		}
	}

  return *this;
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the second one and return a new model
//====================================================================================================================
SC_Model* SC_Model_Pareto::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_Model_Pareto*)pSecond, false));
}
SC_Model_Pareto* SC_Model_Pareto::combineModels(SC_Model_Pareto* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_Model_Pareto::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_Model_Pareto*)pFirst, (SC_Model_Pareto*)pSecond, keepFirstsNext));
}
SC_Model_Pareto* SC_Model_Pareto::combineModels(SC_Model_Pareto* pFirst, SC_Model_Pareto* pSecond, bool keepFirstsNext) {
  SC_Model_Pareto *pNewModel = NULL;
  SV_Data *pData1 = NULL, *pData2 = NULL;

  if (pSecond == NULL || (pFirst->getDim() == 0 && pSecond->getDim() == 0)) {
    pNewModel = new SC_Model_Pareto(*(SC_Model_Pareto*)pFirst);
  } else { //there's nop good way to combine the parameters without reestimation...
    pNewModel = new SC_Model_Pareto(pFirst->getTweak(), NULL, pFirst->getMarginalDistributionUsage());
    pData1 = pFirst->getTrainingData();
    pData2 = pSecond->getTrainingData();
    pData1->Next = pData2;
    pNewModel->TrainModel(pData1, 2);
    MFree_0D(pData2);
    MFree_0D(pData1);
  }

  return pNewModel;
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_Model_Pareto::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_Model_Pareto*)pSecond, pSpeechFrames, segmentsToMerge, (SC_Model_Pareto*)pBackgroundModels));
}
SC_Model_Pareto* SC_Model_Pareto::combineModels(SC_Model_Pareto* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model_Pareto* pBackgroundModels) {
  SC_Model_Pareto *pNewModel = new SC_Model_Pareto(this->pTweak, pBackgroundModels, this->useMarginalDistributions);
	
	//bug: pSpeechFrames is already meant to contain training data for both models in a linked list with segmentsToMerge segments!
	//SV_Data *pThisTrainingData = getTrainingData();
  //pThisTrainingData->Next = pSpeechFrames;
  //pNewModel->TrainModel(pThisTrainingData, 2);

	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

  //MFree_0D(pThisTrainingData);

  return pNewModel;
}

//====================================================================================================================
// Convert the original Training data, which is still there in the pDensities-member, back into it's original form
// (1 single SV_Data object without densitiy-values in the last col)
//====================================================================================================================
SV_Data* SC_Model_Pareto::getTrainingData(void) {
  SV_Data *pTrainingData = NULL;
  unsigned int i;
  int x, y;

  if (this->dim > 0) {
    if (this->dim == 1) { //mutlivariate case
      pTrainingData = new SV_Data(this->pDensities[0]->Row, this->pDensities[0]->Col-1);
      for (y = 0; y < this->pDensities[0]->Row; y++) {
        for (x = 0; x < this->pDensities[0]->Col-1; x++) {
          pTrainingData->Mat[y][x] = this->pDensities[0]->Mat[y][x];
        }
      }
    } else { //marginal case
      pTrainingData = new SV_Data(this->pDensities[0]->Row, this->dim);
      for (y = 0; y < this->pDensities[0]->Row; y++) {
        for (i = 0; i < this->dim; i++) {
          pTrainingData->Mat[y][i] = this->pDensities[i]->Mat[y][0];
        }
      }
    }
  }

  return pTrainingData;
}

//====================================================================================================================
// Train a model; here, only the first segment of a linked list of feature-vectors is used
//====================================================================================================================
int SC_Model_Pareto::TrainModel(SV_Data *TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Train a model; here, one can specify how many segments of a linked list of feature-vectors should be merged
// (0 means all)
//====================================================================================================================
int	SC_Model_Pareto::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
	unsigned int i;
  long int j;
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge), *pSingleCol = NULL;

  //destruct old model, if any
  if (this->dim > 0) { 
    for (i = 0; i < this->dim; i++) {
      MFree_0D(this->pDensities[i]);
    }
    MFree_1D(this->pDensities);
    MFree_1D(this->paretoRadius);  
  }

  //create new one
  this->trainingDataCount = pCompleteData->Row;

  if (this->useMarginalDistributions == true) {
    this->dim = pCompleteData->Col;
  } else {
    this->dim = 1;
  }

  MArray_1D(this->pDensities, this->dim, SV_Data*, "SC_Model_Pareto.TrainModel: pDensities");
  MArray_1D(this->paretoRadius, this->dim, double, "SC_Model_Pareto.TrainModel: paretoRadius");

  //estimate optimal radius and corresponding density for each dim
  if (this->dim == 1) { //multivariate case
    this->pDensities[0] = this->pPDEsolver->estimateDensityFunction(pCompleteData, this->paretoRadius[0]);
  } else { //train each dimension indepentantly
    pSingleCol = new SV_Data(pCompleteData->Row, 1);
    for (i = 0; i < this->dim; i++) {
			this->paretoRadius[i] = 0.0;
      for (j = 0; j < pCompleteData->Row; j++) { //create a copy only of the current dimension/col
        pSingleCol->Mat[j][0] = pCompleteData->Mat[j][i];
      }
      this->pDensities[i] = this->pPDEsolver->estimateDensityFunction(pSingleCol, this->paretoRadius[i]);
    }  
    MFree_0D(pSingleCol);
  }

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  return SVLIB_Ok;
}

//====================================================================================================================
// Test the model while merging no segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_Pareto::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;
}

//====================================================================================================================
// Test the model with specified nr. of segments in the linked list of feature-vectors
//====================================================================================================================
SV_Data* SC_Model_Pareto::TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) {
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(TestData) == 1) ? TestData : TestData->MergeData(segmentsToMerge); //merge all feature-vectors in a linked list together
  SV_Data *pScore = new SV_Data(1, 1);
  double logLikelihood = 0.0, res;
  long int j;
  unsigned int i;

  if (this->dim > 0) {
    for (i = 0; i < this->dim; i++) {
      for (j = 0; j < pCompleteData->Row; j++) {
        res = this->pPDEsolver->getDensity(this->pDensities[i], this->paretoRadius[i], pCompleteData->Mat[j], this->pDensities[i]->Col-1);
        logLikelihood += log(res);
      }
    }
    //TODO: invert?!?
    pScore->Mat[0][0] = (float)(logLikelihood) / (float)(pCompleteData->Row); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...
  } else {
    REPORT_ERROR(SVLIB_Fail, "Can't test on an untrained model");
    MFree_0D(pScore);
  }
  
  if (TestData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  return pScore;
}

//====================================================================================================================
// load a model from file
//====================================================================================================================
SV_Model*	SC_Model_Pareto::LoadModel(void) {
	int res, newDim = 0;
  unsigned int x;
  SV_Data *pTempData = NULL;
	SV_DataIO io;
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
	io.getCurrentDatatypeSizes(codeSizes);

  //read header
  res = LoadHdr(fileSizes);
	if (res == SVLIB_Fail) {return(NULL);}

	//read dim and useMarginalDistribution
	io.readScalar(&(this->DFile), newDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {return(NULL);}

	io.readScalar(&(this->DFile), useMarginalDistributions, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {return(NULL);}

	//prepare this model for changing it's parameters
	if ((newDim != this->dim)) {
    for (unsigned int i = 0; i < this->dim; i++) {
      MFree_0D(this->pDensities[i]);
    }
    MFree_1D(this->pDensities);
    MFree_1D(this->paretoRadius);
  
		this->dim = newDim;
    if (this->dim > 0) {
      MArray_1D(this->pDensities, this->dim, SV_Data*, "SC_Model_Pareto.LoadModel: pDensities");
      MArray_1D(this->paretoRadius, this->dim, double, "SC_Model_Pareto.LoadModel: paretoRadius");
    }
  }

  //read radius-vector
	io.readArray(&(this->DFile), this->paretoRadius, this->dim, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Pareto Model Failed!");}

	//read densities
	for (x = 0; x < this->dim; x++) {
    pTempData = new SV_Data();

		//header
		io.readArray(&(this->DFile), pTempData->Hdr.Signature, 16, codeSizes, fileSizes);
		io.readScalar(&(this->DFile), pTempData->Hdr.ByteOrder, codeSizes, fileSizes);
		io.readScalar(&(this->DFile), pTempData->Hdr.Version, codeSizes, fileSizes);	
		io.readScalar(&(this->DFile), pTempData->Hdr.ID, codeSizes, fileSizes);
		io.readArray(&(this->DFile), pTempData->Hdr.Name, 16, codeSizes, fileSizes);
		io.readScalar(&(this->DFile), pTempData->Hdr.frameSize, codeSizes, fileSizes);
		io.readScalar(&(this->DFile), pTempData->Hdr.frameStep, codeSizes, fileSizes);
		io.readScalar(&(this->DFile), pTempData->Hdr.sampleRate, codeSizes, fileSizes);
		io.readArray(&(this->DFile), pTempData->Hdr.NotUsed, 16, codeSizes, fileSizes);

		//data, part 1
		io.readScalar(&(this->DFile), pTempData->Row, codeSizes, fileSizes);
		io.readScalar(&(this->DFile), pTempData->Col, codeSizes, fileSizes);

    if (DFile.good() != TRUE) {
		  MFree_0D(pTempData);
      REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Pareto Model Failed!");
      return NULL;
	  }

    pTempData->Alloc(); 
		io.readMatrix(&(this->DFile), pTempData->Mat, pTempData->Row, pTempData->Col, codeSizes, fileSizes);

    if (DFile.good() != TRUE) {
		  MFree_0D(pTempData);
      REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Pareto Model Failed!");
      return NULL;
	  }

    this->pDensities[x] = pTempData;
  }

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_Model_Pareto Failed!");}

	return(this);
}

//====================================================================================================================
// save a model to file
//====================================================================================================================
int	SC_Model_Pareto::SaveModel(void) {
	unsigned int x;
 	int res, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Pareto model Failed!");}

  //write dim and useMarginalDistribution
	bytes = io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->useMarginalDistributions);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Pareto Model Failed!");}

  //write radius-vector
	bytes += io.writeArray(&(this->DFile), this->paretoRadius, this->dim);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Pareto Model Failed!");}

	//write densities
  for (x = 0; x < this->dim; x++) {
		//header
		bytes += io.writeArray(&(this->DFile), this->pDensities[x]->Hdr.Signature, 16);
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Hdr.ByteOrder);
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Hdr.Version);	
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Hdr.ID);
		bytes += io.writeArray(&(this->DFile), this->pDensities[x]->Hdr.Name, 16);
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Hdr.frameSize);
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Hdr.frameStep);
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Hdr.sampleRate);
		bytes += io.writeArray(&(this->DFile), this->pDensities[x]->Hdr.NotUsed, 16);

		//data
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Row);
		bytes += io.writeScalar(&(this->DFile), this->pDensities[x]->Col);
		bytes += io.writeMatrix(&(this->DFile), this->pDensities[x]->Mat, this->pDensities[x]->Row, this->pDensities[x]->Col);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Pareto Model Failed!");}
	}

  //write trainingDataCount
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_Model_Pareto Model Failed!");}

  return bytes + MHLen; //MHLen may be incorrect...
}

//====================================================================================================================
// create a (copied, because the oiginal data is float, not double) SC_Signature-view on this model (for distance 
// computation ), and destruct it
//====================================================================================================================
SC_Signature* SC_Model_Pareto::toSignature(void) {
  SC_Centroid **centroids;
  SC_Signature *pSignature = NULL; 
  double *weights, *point;
  long int x, X, y, Y;
  
  if (this->dim > 0) {
    Y = this->pDensities[0]->Row;
		X = (this->dim == 1) ? this->pDensities[0]->Col - 1 : this->dim;

		MArray_1D(centroids, Y, SC_Centroid*, "SC_Model_Pareto.toSignature: centroids");
    MArray_1D(weights, Y, double, "SC_Model_Pareto.toSignature: weights");
   
		for (y = 0; y < Y; y++) {
			MArray_1D(point, X, double, "SC_Model_Pareto.toSignature: point");
			weights[y] = (this->dim == 1) ? 1.0 : this->pDensities[0]->Mat[y][this->pDensities[0]->Col-1];
			for (x = 0; x < X; x++) {
				point[x] = (this->dim == 1) ? this->pDensities[x]->Mat[y][0] : this->pDensities[0]->Mat[y][x];
				if (this->dim == 1) {
					weights[y] *= this->pDensities[x]->Mat[y][1]; 
				}
			}
			centroids[y] = new SC_Centroid_Point(this->pTweak, X, point);
		}

		pSignature = new SC_Signature(centroids, weights, Y, true);
  }

  return pSignature;
}

void SC_Model_Pareto::killSignature(SC_Signature *pSignature) {
  if (pSignature != NULL) {
		MFree_0D(pSignature);
  }
  
  return;
}

//====================================================================================================================
// Dump model's parameter in ASCII 
//====================================================================================================================
ostream& SC_Model_Pareto::modelOut(ostream& os) {
	unsigned long int x, y, z;
	char bmpFileName[sclib::bufferSize];

	os.setf(ios::fixed|ios::basefield);
  os.precision(5);

  os << "MarginalDistributionUsage: " << this->useMarginalDistributions << endl;
  os << "Dimensions: " << this->dim << endl << endl;

  if (this->dim > 0) {
    
    os << "Pareto-Radius per dim: " << endl;
    for (x = 0; x < this->dim; x++) {
      os << setw(10) << this->paretoRadius[x] << "\t\t\t";
    }
    os << endl;

    os << "Training-Data and Density per dim: " << endl;
    for (z = 0; z < (unsigned long)(this->pDensities[0]->Row); z++) {
      for (x = 0; x < this->dim; x++) {
        for (y = 0; y < (unsigned long)(this->pDensities[x]->Col); y++) {
          os << setw(10) << this->pDensities[x]->Mat[z][y] << "\t";
        }
        os << "\t";
      }
      os << endl;
    }

		//draw the distribution as an image (sadly, we can't get the filename corresponding with the open stream...)
		sprintf(bmpFileName, "%spdeModel_%sbmp", this->pTweak->debug.debugDir, tmpnam(NULL)+1);
		if (this->dim == 1) {
			this->pPDEsolver->drawDensities(this->pDensities[0], bmpFileName, 512, 512, this->paretoRadius[0]);
		} else {
			this->pPDEsolver->drawDensities(this->pDensities, this->dim, bmpFileName, 512, 512, this->paretoRadius);
		}
		os << "see " << bmpFileName << " for a visual representation of this PDE model.\n";

  } else {
    os << "The model is still untrained..." << endl;
  }

	os << endl;

	return(os);
}

//====================================================================================================================
// draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
//====================================================================================================================
SV_Data* SC_Model_Pareto::drawSamplesFromDistribution(unsigned long int count) {
  unsigned long int t;
  unsigned short int i, d, D; 
  double randomValue, densitySum;
  SV_Data *pSamples = NULL;

  if (this->pDensities != NULL) {
    D = (this->dim == 1) ? this->pDensities[0]->Col - 1 : this->dim;
    pSamples = new SV_Data(count, D);
    
    for (t = 0; t < count; t++) {
      if (this->dim == 1) { //multidimensional density is modelled here
				randomValue = sclib::rand(0.0, 1.0);
        densitySum = 0.0;
        i = 0;
        while (densitySum < randomValue) {
          densitySum += this->pDensities[0]->Mat[i][D-1];
          i++;
        }
        i--;
        for (d = 0; d < D; d++) {
          pSamples->Mat[t][d] = this->pDensities[0]->Mat[i][d]; //TODO: this way we only draw samples we already saw (during training), this should be more "smeared"... (but: I think we never will use this method, it's just here to fullfill the interface somehow ;-)
        }
      } else { //only the marginals are modelled
        for (d = 0; d < D; d++) {
					randomValue = sclib::rand(0.0, 1.0);
          densitySum = 0.0;
          i = 0;
          while (densitySum < randomValue) {
            densitySum += this->pDensities[d]->Mat[i][2];
            i++;
          }
          i--;
          pSamples->Mat[t][d] = this->pDensities[d]->Mat[i][1]; //TODO: same as above...
        }
      }
    }
  }

  return pSamples;
}
