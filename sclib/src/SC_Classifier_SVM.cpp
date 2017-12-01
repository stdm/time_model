/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates an SVM-classifier based on the libsvm-SVM-         */
/*        implementation v2.81 (with additions by Hsuan-Tien Lin to       */
/*        support weighted examples for c-svc and epsilon-svr)            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.03.2006																								*/
/**************************************************************************/

#include <vector>
#include <math.h>
#include "SC_Classifier_SVM.h"
#include "SC_DistanceMeasures.h"
#include <SV_DataIO.h>
#include <SV_Error.h>

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Classifier_SVM::SC_Classifier_SVM(SC_TweakableParameters* pTweak, bool doParameterSearch, bool doScaling, bool verbose) : SC_ClassifierWithWeights(pTweak, doScaling, verbose) {
	this->classifierType = sclib::ctSVM;
	this->pSVM = new SC_SVM(verbose);
  this->pClassifier = NULL;
  this->pTraingData = NULL;
  this->doParameterSearch = doParameterSearch;
	this->pClassWeightList = NULL;
	this->justLinked = false;
}

//====================================================================================================================
//	copy-constructor
//====================================================================================================================
SC_Classifier_SVM::SC_Classifier_SVM(const SC_Classifier_SVM& pParent, bool justLink) : SC_ClassifierWithWeights(pParent) {
	this->classifierType = sclib::ctSVM;
	this->pSVM = new SC_SVM(this->verbose);
	this->doParameterSearch = pParent.doParameterSearch;
	if (justLink == true) {
		this->pClassifier = pParent.pClassifier;
		this->pTraingData = pParent.pTraingData;
		this->justLinked = true;
	} else {
		this->pClassifier = copySVMmodel(pParent.pClassifier);
		this->pTraingData = copySVMproblem(pParent.pTraingData);
		this->justLinked = false;
	}
	this->pClassWeightList = sclib::copyLinkedList(pParent.pClassWeightList);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Classifier_SVM::~SC_Classifier_SVM() {
	if (this->justLinked == false) {
		killSVMproblem(this->pTraingData);
		if (this->pClassifier != NULL) {
			this->pSVM->svm_destroy_param(&this->pClassifier->param);
			if (this->pClassifier->free_sv == 42) { //dirty hack to circumvent the problem that we don't link to the training data when we copied the classifier, so we have to free it the way we created it
				for (int i = 0; i < this->pClassifier->l; i++) {
					free(this->pClassifier->SV[i]);
				}
				free(this->pClassifier->SV);
				this->pClassifier->SV = NULL;
				this->pClassifier->free_sv = 0; //svm_destroy_model() has to do nothing with the SVs now
			}
			this->pSVM->svm_destroy_model(this->pClassifier);
		}
	}
  MFree_0D(this->pSVM);
	sclib::destructLinkedList(this->pClassWeightList);
}

//====================================================================================================================
//	reduce a training-set by returning a new one with maxSize vectors drawn at random
//====================================================================================================================
SC_SVMproblem* SC_Classifier_SVM::reduceTrainingSet(SC_SVMproblem *pProblem, unsigned long int maxSize) {
  long int i, idx;
	double randomValue, weightsSoFar, weightsOverall = 0.0, weightSum = 0.0;
	SC_SVMproblem *pReduced;
	vector<long int> *idxVec;
 
	pReduced = new SC_SVMproblem;
	pReduced->l = maxSize; //select a subset randomly without recline 

  if (this->verbose == true) {
		printf("\n    Reducing training-set to %02.1f%% original size from %d to %d points ", (double)(pReduced->l)/(double)(pProblem->l)*100.0, pProblem->l, pReduced->l);
	}

  MArray_1D(pReduced->y, pReduced->l, double, "SC_Classifier_SVM.reduceTrainingSet: pReduced->y");
  MArray_1D(pReduced->x, pReduced->l, SC_SVMnode*, "SC_Classifier_SVM.reduceTrainingSet: pReduced->x");
	if (pProblem->W != NULL) {
		MArray_1D(pReduced->W, pReduced->l, double, "SC_Classifier_SVM.reduceTrainingSet: pReduced->W");
	} else {
		pReduced->W = NULL;
	}
  idxVec = new vector<long int>(pProblem->l);

  for (i = 0; i < pProblem->l; i++) { //create a list of all possibly selectable indices of the original problem
    idxVec->push_back(i);
		if (pProblem->W != NULL) {
			weightsOverall += pProblem->W[i]; //get overall sum of weights (should be 1, but who knows...)
		}
  }

  for (i = 0; i < pReduced->l; i++) {
		if (pReduced->W == NULL) { //randomly select an *index* of the still available indices of the original problem
			randomValue = sclib::rand(0.0, (double)(idxVec->size()-1));
			idx = sclib::round(randomValue);
		} else { //select an index regarding the distribution of the examples in the weights
			randomValue = sclib::rand(0.0, weightsOverall); //randomly select a *weight-sum*, that, when reached, determines the index to use
			weightsSoFar = 0.0;
			idx = 0;
			while (weightsSoFar < randomValue) {
				weightsSoFar += pProblem->W[(*idxVec)[idx]];
				idx++;
			}
			if (idx > 0) {
				idx--;
			}
			weightsOverall -= pProblem->W[(*idxVec)[idx]];
		}

		pReduced->y[i] = pProblem->y[(*idxVec)[idx]];
    pReduced->x[i] = copySVMnodes(pProblem->x[(*idxVec)[idx]]);
		if (pReduced->W != NULL) {
			pReduced->W[i] = pProblem->W[(*idxVec)[idx]];
			weightSum += pReduced->W[i];
		}

    idxVec->erase(idxVec->begin() + idx); //remove the just selected index from the list of available (still selectable) indices
  }
	if (pReduced->W != NULL) {
		for (i = 0; i< pReduced->l; i++) { //renormalize weights to sum up to 1
			pReduced->W[i] /= weightSum;
		}
	}

  MFree_0D(idxVec);

	return pReduced;
}

//====================================================================================================================
//	find suitable parameter settings accroding to the tweakable parameters
//====================================================================================================================
SC_SVMparameter* SC_Classifier_SVM::getParameters(SC_SVMproblem *pProblem) {
  long int i, j, count = 0, overallCount, dim = 2; //dim=2 so if no problem is given still a meaningful value (0.5) for gamma will be used
  int fold = this->pTweak->classifierSvm.cvFolds; //nr. of cross-validation folds
	double bestGamma, bestC, bestAccuray, erg, *target;
  double gammaStart, gammaEnd, gammaStep, CStart, CEnd, CStep;
  double finished, product, maxSquaredNorm = 0.0, avgNorm = 1.0; //avgNorm=1 so if no problem is given still a meaningful value (1) for C will be used
	double increment;
  SC_SVMparameter *pParam = new SC_SVMparameter;
  SC_SVMproblem *pCVproblem = pProblem; //in case of small enough a problem, use all examples for parameter-search

	//get dimensionality of data by looking at the greatest index in all rows (remember: a svm-problem is a sparse data matrix)
	//also get the maximum squared euclidean norm of the examples, which is used to compute an upper bound for C during parameter search as stated on the svmlight website ("The value of C should typically be less than 1000 times 1.0/(the maximum squared Euclidian length of your feature vectors)", http://svmlight.joachims.org/)
	if (pProblem != NULL) {
		avgNorm = 0.0;
		for (i = 0; i < pProblem->l; i++) {
			j = 0;
			product = 0.0;
			while (pProblem->x[i][j].index != -1) {
				product += pProblem->x[i][j].value * pProblem->x[i][j].value;
				j++;
			}
			avgNorm += sqrt(product) / (double)(pProblem->l);
			if (product > maxSquaredNorm) {
				maxSquaredNorm = product;
			}
			if (j > 0) {
				if (pProblem->x[i][j-1].index > dim) {
					dim = pProblem->x[i][j-1].index;
				}
			}
		}
	}

	//default values
	pParam->svm_type = this->pTweak->classifierSvm.svm_type;
	pParam->kernel_type = this->pTweak->classifierSvm.kernel_type;
	pParam->degree = this->pTweak->classifierSvm.degree;
	pParam->gamma = (this->pTweak->classifierSvm.gamma > 0) ? this->pTweak->classifierSvm.gamma : 1.0/(double)(dim); // 1/dim as a default if 0 is provided in the parameters
	pParam->coef0 = this->pTweak->classifierSvm.coef0;
	pParam->nu = this->pTweak->classifierSvm.nu;
	pParam->cache_size = this->pTweak->classifierSvm.cache_size;
	pParam->C = (this->pTweak->classifierSvm.C > 0) ? this->pTweak->classifierSvm.C : 1.0/avgNorm; //the default parameter of svmlight, which works well there according to ralph ewerth
	pParam->eps = (this->pTweak->classifierSvm.eps == 0.0) ? ((pParam->svm_type==SCLIB_SVM_TYPE_NUSVC)?0.00001:0.001) : this->pTweak->classifierSvm.eps; //if pTweak gives 0.0, use the type-dependant defaults of the libsvm-guys (http://www.cise.ufl.edu/~fu/Teach/Learn/SVM.txt)
	pParam->p = this->pTweak->classifierSvm.p;
	pParam->shrinking = this->pTweak->classifierSvm.shrinking;
	pParam->probability = this->pTweak->classifierSvm.probability;
	classWeights2parameters(pParam); //handle class-weights

  //if parameter-search is wished, do grid-search for best parameters using cross-validation:
  if (this->doParameterSearch == true && pProblem != NULL) {

    if (this->verbose == true) {
			printf("\nStarting SVM-best-parameter-search using cross-validation:");
		}

    //randomly select samples from the SVMproblem to create a smaller subset of the data
    if (this->pTweak->classifierSvm.cvMaxDatasetSize > 0 && pProblem->l > this->pTweak->classifierSvm.cvMaxDatasetSize) {
			pCVproblem = reduceTrainingSet(pProblem, this->pTweak->classifierSvm.cvMaxDatasetSize);
    }

    //do two iterations, first coarse, then finer near the previously found maximum
    MArray_1D(target, pCVproblem->l, double, "SC_Classifier_SVM.getParameters: target");
    for (int x = 0; x < 2; x++) {

      if (x == 0) { //parameters for coarse search
				gammaStart = (this->pTweak->classifierSvm.cvCoarseGammaMin == 0.0) ? sclib::log2(pParam->gamma) : this->pTweak->classifierSvm.cvCoarseGammaMin;
				gammaEnd = (this->pTweak->classifierSvm.cvCoarseGammaMax == 0.0) ? sclib::log2(pParam->gamma) : this->pTweak->classifierSvm.cvCoarseGammaMax;
        gammaStep = this->pTweak->classifierSvm.cvCoarseGammaStep;
        CStart = this->pTweak->classifierSvm.cvCoarseCMin;
        CEnd = sclib::min(this->pTweak->classifierSvm.cvCoarseCMax, sclib::max(CStart, sclib::log2(1000.0/maxSquaredNorm))); //1000/maxSquaredNorm is an upper bound for C as stated on the svmlite homepage: "The value of C should typically be less than 1000 times 1.0/(the maximum squared Euclidian length of your feature vectors)", http://svmlight.joachims.org/
        CStep = this->pTweak->classifierSvm.cvCoarseCStep;
        overallCount = (long int)ceil(((1 + fabs(gammaEnd - gammaStart) / gammaStep) * (1 + fabs(CEnd - CStart) / CStep)) + (((this->pTweak->classifierSvm.cvFineGammaRadius*2 + 1) / this->pTweak->classifierSvm.cvFineGammaStep) * ((this->pTweak->classifierSvm.cvFineCRadius*2 + 1) / this->pTweak->classifierSvm.cvFineCStep)));
        bestGamma = (pParam->gamma != 0.0) ? sclib::log2(pParam->gamma) : 0.0;
        bestC = (pParam->C != 0.0) ? sclib::log2(pParam->C) : 0.0;
        bestAccuray = -1.0;
      } else { //parameters for finer search
        overallCount = (long int)ceil((1 + fabs(gammaEnd - gammaStart) / gammaStep) * (1 + fabs(CEnd - CStart) / CStep));
        gammaStart = sclib::max(gammaStart, bestGamma-this->pTweak->classifierSvm.cvFineGammaRadius);
        gammaEnd = sclib::min(gammaEnd, bestGamma+this->pTweak->classifierSvm.cvFineGammaRadius);
        gammaStep = this->pTweak->classifierSvm.cvFineGammaStep;
        CStart = sclib::max(CStart, bestC-this->pTweak->classifierSvm.cvFineCRadius);
        CEnd = sclib::min(CEnd, bestC+this->pTweak->classifierSvm.cvFineCRadius);
        CStep = this->pTweak->classifierSvm.cvFineCStep;
        overallCount += (long int)ceil((1 + fabs(gammaEnd - gammaStart) / gammaStep) * (1 + fabs(CEnd - CStart) / CStep)); //maybe the overallCount decreased a little bit as compared to the calculation in the first loop because of the min/max-functions involved in finding the start/end, so recalculate it
      }
      
      //do grid search for best C and gamma in the defined range
	    for (double c = CStart; c <= CEnd; c += CStep) {
		    for (double g = gammaStart; g <= gammaEnd; g += gammaStep) {
			    pParam->C = pow(2.0, c);
			    pParam->gamma = pow(2.0, g);
    			
          if (checkParameters(pCVproblem, pParam) == true) {
			      this->pSVM->svm_cross_validation(pCVproblem, pParam, fold, target);

            //compute accuracy as the ratio of correctly predicted labels to the complete number of labels
						increment = 1.0 / (double)(pCVproblem->l);
            erg = 0.0;
            for (int i = 0; i < pCVproblem->l; i++) {
							if (pCVproblem->W != NULL) {
								increment = pCVproblem->W[i]; //if possible, calculate weighted accuracy
							}
              erg += (pCVproblem->y[i] == target[i]) ? increment : 0.0;
            }

            count++;
            finished = (double)(count)/(double)(overallCount) * 100.0;
            if (this->verbose == true) {
							printf("\nCV-rate: %02.5f,\tlog-C=%02.1f,\tlog-gamma=%02.1f (%02.1f%% finished)", erg, c, g, finished);
						}

			      if(erg > bestAccuray) {
				      bestAccuray = erg;
				      bestC = c; //save logarithms of these parameters, for better readability
				      bestGamma = g;
			      }
          } else {
            REPORT_ERROR(SVLIB_BadData, "Wrong parameter settings for SVM");
          }
 		    }
      }
      if (this->verbose == true && x == 0) {
				printf("\nCoarse CV search finished, doing fine search...\n");
			}
    }
    MFree_1D(target);
    
    if (this->verbose == true) {
			printf("\nBest parameters found on %02.1f%% of the Trainingsset: log-C=%02.1f, log-gamma=%02.1f, CV-rate=%02.5f", (double)(pCVproblem->l)/(double)(pProblem->l)*100.0, bestC, bestGamma, bestAccuray);
		}
    pParam->C = pow(2.0, bestC);
    pParam->gamma = pow(2.0, bestGamma);

    if (pCVproblem != pProblem) {
      killSVMproblem(pCVproblem);
    }
  } //if doParameterSearch

	if (this->verbose == true) {
		printf("\nparameters found for SVM training: gamma=%f, C=%f, eps=%f", pParam->gamma, pParam->C, pParam->eps);
	}

  return pParam;
}

//====================================================================================================================
//	check the given parameter-set for correctness/applicability
//====================================================================================================================
bool SC_Classifier_SVM::checkParameters(SC_SVMproblem *pProblem, SC_SVMparameter *pParam) {
  const char *report = this->pSVM->svm_check_parameter(pProblem, pParam); //are this parameters suitable?
  bool res = true;
  
  if (report != NULL) {
    char *errStr = new char[sclib::bufferSize];
    sprintf(errStr, "Wrong parameters for SVM: %s", report);
    REPORT_ERROR(SVLIB_BadArg, errStr);
    MFree_1D(errStr);
    res = false;
  } 

  return res;
}

//====================================================================================================================
//  return a gamma-parameter used in one-class SVM computation that is optimal with respect to the given data
//  this also creates a scale-matrix if scaling is switched on for this classifier
//====================================================================================================================
double SC_Classifier_SVM::getOneClassGammaParameter(SV_Data *pData) {
	SC_SVMproblem *pProblem;
	SC_SVMparameter *pParam = getParameters(); //get only standard-parameters from pTweak without standard grid serach
	double gamma;
	
	pParam->svm_type = SCLIB_SVM_TYPE_ONECLASS; //assure one-class nu-SVM
	pParam->kernel_type = SCLIB_SVM_KERNEL_RBF; //RBF (gaussian) kernel as in Desobry et al's paper on kernel change detection
	pParam->probability = 0; //probability-outputs are not supported yet for one-class SVMs

  if (this->doScaling == true) {
    this->pScale = findScalingParameters(pData);
  }
  pProblem = svData2svmProblem(pData); //convert to accepted data format
	gamma = getOneClassGammaParameter(pProblem, pParam);

  this->pSVM->svm_destroy_param(pParam);
	MFree_0D(pParam);
	killSVMproblem(pProblem);

	return gamma;
}

//====================================================================================================================
//  return a gamma-parameter used in one-class SVM computation that is optimal with respect to the given data
//====================================================================================================================
double SC_Classifier_SVM::getOneClassGammaParameter(SC_SVMproblem *pProblem, SC_SVMparameter *pParam) {
	double res, gamma, count;
	double avgDist = 0.0, minDist = std::numeric_limits<double>::max(), maxDist = std::numeric_limits<double>::max()*-1.0;
	int t, tt, r, iterations;
	double lastPercentage = 0.0;
	SC_SVMproblem *pCVproblem = pProblem;

	if (this->verbose == true) {
		printf("\n  Search for best gamma in one-class SVM: ");
	}

	//find min, max and average pairwise distance among feature vectors
	count = (double)(pCVproblem->l)*(double)(pCVproblem->l-1) / 2.0; //i.e. number of loops for parwise distances
	if (count < (double)(this->pTweak->classifierSvm.oneClassGammaSearchMaxIterations*this->pTweak->classifierSvm.oneClassGammaSearchRepeats)) {
		if (this->verbose == true) {
			printf("\n    Computing all %f pairwise distances: ", count);
		}
		for (t = 0; t < pCVproblem->l-1; t++) {
			for (tt = t+1; tt < pCVproblem->l; tt++) {
				res = SC_DistanceMeasures::euclid(pCVproblem->x[t], pCVproblem->x[tt]);
				avgDist += res / count; //divide here to avoid overflow...
				if (res<minDist && res>0.0) {
					minDist = res;
				} else if (res > maxDist) {
					maxDist = res;
				}
			}
			if (this->verbose == true) {
				lastPercentage = sclib::printPercentage(pCVproblem->l-1, t, lastPercentage, 1.0, t==0);
			}
		}
	} else {
		if (this->verbose == true) {
			printf("\n    Computing %d/%f pairwise distances at random, make %d repetitions: ", this->pTweak->classifierSvm.oneClassGammaSearchMaxIterations, count, this->pTweak->classifierSvm.oneClassGammaSearchRepeats);
		}
		//if the number of elements is so high that computing all pairwise distances is too burdensome, we use a heuristic:
		//randomly (with recline) select pairs of vectors and accumualte their distance for the mean distance
		count = (double)(this->pTweak->classifierSvm.oneClassGammaSearchMaxIterations * this->pTweak->classifierSvm.oneClassGammaSearchRepeats);
		for (r = 0; r < this->pTweak->classifierSvm.oneClassGammaSearchRepeats; r++) {
			if (this->verbose == true) {
				printf("r=%d ", r);
			}
			iterations = 0;
			do {
				t = sclib::rand(pCVproblem->l-1);
				do { //if we only choose tt>t, we get a wrong distribution of distances (the distances betwen the last 2 vetors will occur very ofton while distances between the first and any other vector are sampled only sparsely); this is due to rand() yielding uniformly distributed reandom numbers
					tt = sclib::rand(pCVproblem->l-1);
				} while (tt == t); //t must be different from tt

				res = SC_DistanceMeasures::euclid(pCVproblem->x[t], pCVproblem->x[tt]);
				avgDist += res / count;
				if (res<minDist && res>0.0) {
					minDist = res;
				} else if (res > maxDist) {
					maxDist = res;
				}
				if (this->verbose == true) {
					lastPercentage = sclib::printPercentage(this->pTweak->classifierSvm.oneClassGammaSearchMaxIterations, iterations, lastPercentage, 1.0, iterations==0);
				}
			} while (iterations++ < this->pTweak->classifierSvm.oneClassGammaSearchMaxIterations);
		}
	} //if count < maxIterations...
	if (this->verbose == true) {
		printf("\n    Found distances: min=%f, average=%f, max=%f", minDist, avgDist, maxDist);
	}
	
	double maxGamma = log(1.0 / (minDist*minDist)); //gamma = 1/s^2 
	double minGamma = log(1.0 / (maxDist*maxDist));
	double step = (this->pTweak->classifierSvm.cvCoarseGammaMax-this->pTweak->classifierSvm.cvCoarseGammaMin) / this->pTweak->classifierSvm.cvCoarseGammaStep;
	double positiveRate, bestPositiveRate, difference, smallestDifference = std::numeric_limits<double>::max();
	double *target;
	double svFraction;
	double lastDiff;
	int noImprovementCount;

	//TODO: hack: take minGamma as gamma
	//maxGamma = minGamma-1.0;
	//gamma = exp(minGamma+step);

	//if the training-set gets shrinked here, repeat the search to see if sampling is good
	for (int k = 0; k < this->pTweak->classifierSvm.cvFolds; k++) {
		//randomly select samples from the SVMproblem to create a smaller subset of the data
		if (this->pTweak->classifierSvm.cvMaxDatasetSize > 0 && pProblem->l > this->pTweak->classifierSvm.cvMaxDatasetSize) {
			pCVproblem = reduceTrainingSet(pProblem, this->pTweak->classifierSvm.cvMaxDatasetSize);
		}

		//tweak the gamma parameter (=1/(2*sigma^2) in the gaussian sense) as proposed by Desobry et al: 
		//sigma should be one order of magnitude smaller than the average distance between training vectors
		//gamma = 0.5 * pow(avgDist/10.0, 2.0);
		// => this method doesn't seem to work on our data!!!

		//New method from D. Tax "One-Class Classification - Concept-Learning in the Absence of Counter-Examples", 2001, chapter 2.3-2.4:
		//In this thesis, a technique called SVDD is introduced which is quite similar to one-class SVM (chapter 2.5 shows this). It
		//also uses a gaussian (rbf) kernel with a kernel parameter s (gamma=1/s^2), which can be tweaked between the bounds s_min=minDist,
		//s_max=maxDist (the former one giving a tight parzen window density estimator with small gaussian kernels around almost every training 
		//point as a SV, the latter giving a loose hypershphere with few SVs). We compute a minimum and a maximum gamma from this bounds and
		//step trough them in equal intervals (taken from the tweakable parameters if coarse gamma values are given there, otherwise we 
		//part in 25 parts) in search for the gamma that finds a number of positively classified training examples most near the wanted 1-nu.
		if (step == 0.0) {
			step = (maxGamma - minGamma) / 25;
		}
		if (this->verbose == true) {
			printf("\n    Starting %d/%d search for best gamma in the range log(gamma)=%f..%f, step=%f: ", k+1, ((pCVproblem==pProblem)?1:this->pTweak->classifierSvm.cvFolds), minGamma, maxGamma, step);
		}
		this->pSVM->setVerbose(false);

		//do grid search for best gamma
		noImprovementCount = 0;
		lastDiff = 1.0;
		double oldG = pParam->gamma;
		pParam->gamma = 0.0;
		MArray_1D(target, pCVproblem->l, double, "SC_Classifier_SVM.getOneClassGammaParameter: target");
		for (double g = minGamma; g <= maxGamma; g += step) {
			pParam->gamma = exp(g);
			
			if (checkParameters(pCVproblem, pParam) == true) {
				this->pSVM->svm_cross_validation(pCVproblem, pParam, this->pTweak->classifierSvm.cvFolds, target, svFraction);

				//compute ratio of positives
				positiveRate = 0.0;
				for (t = 0; t < pCVproblem->l; t++) {
					if (sclib::round(target[t]) == sclib::labelPositive) {
						positiveRate += 1.0;
					}
				}
				positiveRate /= (double)(pCVproblem->l);

				//compare it with nu
				difference = fabs(positiveRate - (1.0-this->pTweak->classifierSvm.nu));

				if (difference<smallestDifference && positiveRate>0.0) {
					gamma = pParam->gamma;
					bestPositiveRate = positiveRate;
					smallestDifference = difference;
				}
				if (difference-lastDiff >= 0.0) {
					noImprovementCount++;
				} else {
					noImprovementCount = 0;
				}

				if (this->verbose == true) {
					printf("\n      gamma=%f => pos=%f diff=%f sv=%f", pParam->gamma, positiveRate, difference, svFraction);
				}

				if (noImprovementCount > 3) {
					if (this->verbose == true) {
						printf(" (%dx no improvement => out)", noImprovementCount);
					}
					break;
				}
				lastDiff = difference;
			} else {
				REPORT_ERROR(SVLIB_BadData, "Wrong parameter settings for SVM");
			}
		}
		MFree_1D(target);
		pParam->gamma = oldG;
		this->pSVM->setVerbose(this->verbose);

		if (pCVproblem != pProblem) {
			killSVMproblem(pCVproblem);
		} else {
			break;
		}
	}

	if (this->verbose == true) {
		printf("\n  Search for best gamma finished: gamma=%f", gamma);
		sclib::scalarOut("gamma.txt", gamma, this->pTweak, true, "\n");
	}

	return gamma;
}

//====================================================================================================================
//	train a one-class nu-SVM (more of a model than a classifier)
//====================================================================================================================
int SC_Classifier_SVM::trainOneClass(SV_Data *pData) {
  SC_SVMproblem *pProblem = NULL;
	int res;

  if (this->doScaling == true) {
    this->pScale = findScalingParameters(pData);
  }
  pProblem = svData2svmProblem(pData); //convert to accepted data format
	res = trainOneClass(pProblem);

	if (this->pClassifier == NULL) { //don't kill the problem if a model was created, because it links to the problem's data
    killSVMproblem(pProblem);
	}

	return res;
}

//====================================================================================================================
//	train a one-class nu-SVM (more of a model than a classifier) using the SVM-Problem data structure directly as 
//  input
//====================================================================================================================
int SC_Classifier_SVM::trainOneClass(SC_SVMproblem *pProblem) {
  SC_SVMparameter *pParam = NULL;

  //destroy previously trained classifier (and corresponding data), if any
  MFree_0D(this->pScale);
  killSVMproblem(this->pTraingData);
  if (this->pClassifier != NULL) { 
		this->pSVM->svm_destroy_param(&this->pClassifier->param);
    this->pSVM->svm_destroy_model(this->pClassifier);
  }
	this->isTrained = false;
	this->classCount = 0;	

	pParam = getParameters(); //get only standard-parameters from pTweak without standard grid serach
	pParam->svm_type = SCLIB_SVM_TYPE_ONECLASS; //assure one-class nu-SVM
	pParam->kernel_type = SCLIB_SVM_KERNEL_RBF; //RBF (gaussian) kernel as in Desobry et al's paper on kernel change detection
	pParam->probability = 0; //probability-outputs are not supported yet for one-class SVMs

	if (this->doParameterSearch == true) {
		pParam->gamma = getOneClassGammaParameter(pProblem, pParam);
	}
  
  if (checkParameters(pProblem, pParam) == true) {
    this->pClassifier = this->pSVM->svm_train(pProblem, pParam);
  }
  this->pSVM->svm_destroy_param(pParam);
	MFree_0D(pParam);
	this->pClassifier->param.nr_weight = 0; //because the referenced parameter-set was killed above
	this->pClassifier->param.weight = NULL;
	this->pClassifier->param.weight_label = NULL;

  if (this->pClassifier != NULL) {
    this->pTraingData = pProblem;
    this->isTrained = true;
		this->classCount = 1;
  }

  return (this->pClassifier != NULL) ? SVLIB_Ok : SVLIB_Fail;
}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//====================================================================================================================
int SC_Classifier_SVM::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative) {
	return trainTwoClass(pPositive, pNegative, NULL, NULL);
}

//====================================================================================================================
//	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
//  in the two-class-case, the constants SCIB_LABEL_POSITIVE/sclib::labelNegative should be used to indicate the 
//  classes, and they normaly should evaluate to +1/-1, respectively
//  the weights in the two arrays together (corresponding to rows in the SV_Data objects) are meant to sum up to 1
//====================================================================================================================
int SC_Classifier_SVM::trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double *positiveWeights, double *negativeWeights) {
  SC_SVMproblem *pProblem = new SC_SVMproblem;
  SC_SVMparameter *pParam = NULL;

	if (positiveWeights != NULL || negativeWeights != NULL) {
		REPORT_ERROR(SVLIB_BadArg, "providing weights per example for SVM training seems to be corrupt (bugs in this specific libsvm extension)!");
	}

  //destroy previously trained classifier (and corresponding data), if any
  MFree_0D(this->pScale);
  killSVMproblem(this->pTraingData);
  if (this->pClassifier != NULL) {
		this->pSVM->svm_destroy_param(&this->pClassifier->param);
    this->pSVM->svm_destroy_model(this->pClassifier);
  }
  this->isTrained = false;
	this->classCount = 0;

  if (this->doScaling == true) {
    SV_Data *pHook = pPositive->Next;
    SV_Data *pCompleteData;
    pPositive->Next = pNegative;
    pCompleteData = pPositive->MergeData(2);
    this->pScale = findScalingParameters(pCompleteData);
    pPositive->Next = pHook;
  }
  pProblem = svData2svmProblem(pPositive, pNegative, positiveWeights, negativeWeights); //convert to accepted data format
	pParam = getParameters(pProblem); //get best parameter-settings according to the tweakable parameters
  
  if (checkParameters(pProblem, pParam) == true) {
    this->pClassifier = this->pSVM->svm_train(pProblem, pParam);
  }
  this->pSVM->svm_destroy_param(pParam);
	MFree_0D(pParam);
	this->pClassifier->param.nr_weight = 0; //because the referenced parameter-set was killed above
	this->pClassifier->param.weight = NULL;
	this->pClassifier->param.weight_label = NULL;

  if (this->pClassifier == NULL) { //don't kill the problem if a model was created, because it links to the problem's data
    killSVMproblem(pProblem);
    MFree_0D(this->pScale);
    return SVLIB_Fail;
  } else {
    this->pTraingData = pProblem;
    this->isTrained = true;
		this->classCount = 2;
    return SVLIB_Ok;
  }
}

//====================================================================================================================
//	train a classifier for distinguishing between several classes
//  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
//  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
//  respective row of pData.
//====================================================================================================================
int SC_Classifier_SVM::trainMultiClass(SV_Data *pData, int *classes) {
	return trainMultiClass(pData, classes, NULL);
}

//====================================================================================================================
//	train a classifier for distinguishing between several classes
//  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
//  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
//  respective row of pData.
//  the weights (corresponding to rows in the SV_Data object) are meant to sum up to 1
//====================================================================================================================
int SC_Classifier_SVM::trainMultiClass(SV_Data *pData, int *classes, double *weights) {
  SC_SVMproblem *pProblem = new SC_SVMproblem;
  SC_SVMparameter *pParam = NULL;

	if (weights != NULL) {
		REPORT_ERROR(SVLIB_BadArg, "providing weights per example for SVM training seems to be corrupt (bugs in this specific libsvm extension)!");
	}

  //destroy previously trained classifier (and corresponding data), if any
  MFree_0D(this->pScale);
  killSVMproblem(this->pTraingData);
  if (this->pClassifier != NULL) { 
		this->pSVM->svm_destroy_param(&this->pClassifier->param);
    this->pSVM->svm_destroy_model(this->pClassifier);
  }
	this->classCount = 0;

  if (this->doScaling == true) {
    this->pScale = findScalingParameters(pData);
  }
  pProblem = svData2svmProblem(pData, classes, weights); //convert to accepted data format
  pParam = getParameters(pProblem); //get best parameter-settings according to the tweakable parameters
  
  if (checkParameters(pProblem, pParam) == true) {
    this->pClassifier = this->pSVM->svm_train(pProblem, pParam);
  }
  this->pSVM->svm_destroy_param(pParam);
	MFree_0D(pParam);
	this->pClassifier->param.nr_weight = 0; //because the referenced parameter-set was killed above
	this->pClassifier->param.weight = NULL;
	this->pClassifier->param.weight_label = NULL;

  if (this->pClassifier == NULL) { //don't kill the problem if a model was created, because it links to the problem's data
    killSVMproblem(pProblem);
  } else {
    this->pTraingData = pProblem;
    this->isTrained = true;
		this->classCount = getDistinctClassCount(classes, pData->Row);
  }

  return (this->pClassifier != NULL) ? SVLIB_Ok : SVLIB_Fail;
}

//====================================================================================================================
//	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
//  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
//
//  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
//====================================================================================================================
int* SC_Classifier_SVM::classify(SV_Data *pData, SV_Data* &pProbabilities) {
  int *classes = NULL;
	double *probabilities = NULL;
  SC_SVMnode *pFeatureVector = NULL;

  if (this->doScaling == true && this->pScale == NULL) {
    REPORT_ERROR(SVLIB_BadData, "Can't do scaling if no scaling parameters where found");
  }

  if (this->isTrained == true) {
    if (this->pClassifier != NULL) {
  
			MFree_0D(pProbabilities);
			if (this->pTweak->classifierSvm.probability == 1) {
				pProbabilities = new SV_Data(pData->Row, this->pClassifier->nr_class);
				MArray_1D(probabilities, this->pClassifier->nr_class, double, "SC_Classifier_SVM.classify: probabilities");
				for (long int x = 0; x < this->pClassifier->nr_class; x++) { //to intialize array in case the classifier doesn't support probability feedback
					probabilities[x] = 1.0;
				}
			}
      MArray_1D(classes, pData->Row, int, "SC_Classifier_SVM.classify: classes");

      for (long int y = 0; y < pData->Row; y++) { //classify each feature vector separately
        pFeatureVector = svDataRow2svmNode(pData, y);
				if (this->pTweak->classifierSvm.probability == 1) {
					classes[y] = (int)(this->pSVM->svm_predict_probability(this->pClassifier, pFeatureVector, probabilities));
					for (long int x = 0; x < this->pClassifier->nr_class; x++) {
						pProbabilities->Mat[y][x] = (float)(probabilities[x]);
					}
				} else {
					classes[y] = (int)(this->pSVM->svm_predict(this->pClassifier, pFeatureVector));
				}
        MFree_1D(pFeatureVector);
      }

			MFree_1D(probabilities);
    } else {
      this->isTrained = false;
    }
  }

  return classes;
}

//====================================================================================================================
//	save a trained classifier to a file
//====================================================================================================================
int SC_Classifier_SVM::saveClassifier(const char *fileName) {
  int res = SVLIB_Fail;
  char *scaleFileName = NULL;
  SV_DataIO io;

  if (strlen(fileName) > 0) {
    if (this->isTrained == true) {
      if (this->pClassifier != NULL) {
        res = this->pSVM->svm_save_model(fileName, this->pClassifier);
        if (res == SVLIB_Ok && this->doScaling == true) { //save also the scaling parameters of this training set to use it with the test-data to come
          scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
          io.OpenFile(scaleFileName, WRITE_REC);
		      res = io.PutDataRec(*this->pScale);
          res = (res == 0) ? SVLIB_Fail : SVLIB_Ok; //translate between different traditions to report errrors or success...
          io.CloseFile();
          MFree_1D(scaleFileName);
        }
      } else {
        this->isTrained = false;
      }
    }
  }

  return res;
}

//====================================================================================================================
//	load a trained classifier from a file
//====================================================================================================================
int SC_Classifier_SVM::loadClassifier(const char *fileName) {
  int res;
  char *scaleFileName = NULL;
  SV_DataIO io;

  //destroy previously trained classifier, if any
  killSVMproblem(this->pTraingData);
  if (this->pClassifier != NULL) { 
		this->pSVM->svm_destroy_param(&this->pClassifier->param);
    this->pSVM->svm_destroy_model(this->pClassifier);
  }

  //load classifier
  this->pClassifier = this->pSVM->svm_load_model(fileName);
  if (this->pClassifier != NULL) {
    res = SVLIB_Ok;
    this->isTrained = true;
  } else {
    res = SVLIB_Fail;
    this->isTrained = false;
  }

  //also try to load (if exists) the scaling parameters:
  MFree_0D(this->pScale);
  if (res == SVLIB_Ok && this->doScaling == true) {
    scaleFileName = sclib::exchangeFileExtension(fileName, ".scale");
    if (sclib::fileExists(scaleFileName) == true) {
      io.OpenFile(scaleFileName, READ_REC);
      this->pScale = io.GetAllRec();
      io.CloseFile();
    }
    MFree_1D(scaleFileName);
  }

	//set classCount
	this->classCount = this->pClassifier->nr_class;

  return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns -1 if something goes wrong
//====================================================================================================================
long int SC_Classifier_SVM::label2idx(long int label) {
	long int i, res = -1;
	
	if (this->isTrained == true) {
		for (i = 0; i < this->pClassifier->nr_class; i++) {
			if (this->pClassifier->label[i] == label) {
				res = i;
			}
		}
	}

	return res;
}

//====================================================================================================================
//	to convert between given labels and indices into the probability-parameter of the classifiy()-method
//  returns sclib::noType if something goes wrong
//====================================================================================================================
long int SC_Classifier_SVM::idx2label(long int idx) {
	long int res = sclib::noType;
	
	if (this->isTrained == true) {
		if (idx >= 0 && idx < this->pClassifier->nr_class) {
			res = this->pClassifier->label[idx];
		}
	}

	return res;
}

//====================================================================================================================
//	convert a SV_Data dataset into a format that can be subject to SVM-training; the class-labels are meant to be 
//  given in the classes-array per row of pData
//====================================================================================================================
SC_SVMproblem* SC_Classifier_SVM::svData2svmProblem(SV_Data *pData, int *classes, double *weights) {
  SC_SVMproblem *pProblem = new SC_SVMproblem;

  pProblem->l = pData->Row;
  MArray_1D(pProblem->y, pData->Row, double, "SC_Classifier_SVM.svData2svmProblem: pProblem->y");
  MArray_1D(pProblem->x, pData->Row, SC_SVMnode*, "SC_Classifier_SVM.svData2svmProblem: pProblem->x");
	if (weights != NULL) {
		MArray_1D(pProblem->W, pData->Row, double, "SC_Classifier_SVM.svData2svmProblem: pProblem->W");
	} else {
		pProblem->W = NULL;
	}

  for (long int i = 0; i < pData->Row; i++) {
    pProblem->y[i] = (double)(classes[i]);
    pProblem->x[i] = svDataRow2svmNode(pData, i);
		if (weights != NULL) {
			pProblem->W[i] = weights[i];
		}
  }

  return pProblem;
}

//====================================================================================================================
//	convert a positive and a negative example SV_Data dataset into a format that can be subject to SVM-training
//  can also be used for one-class SVM training (with pNegative==NULL)
//====================================================================================================================
SC_SVMproblem* SC_Classifier_SVM::svData2svmProblem(SV_Data *pPositive, SV_Data *pNegative, double *positiveWeights, double *negativeWeights) {
  SC_SVMproblem *pProblem = new SC_SVMproblem;

	pProblem->l = pPositive->Row + ((pNegative != NULL) ? pNegative->Row : 0);
  MArray_1D(pProblem->y, pProblem->l, double, "SC_Classifier_SVM.svData2svmProblem: pProblem->y");
  MArray_1D(pProblem->x, pProblem->l, SC_SVMnode*, "SC_Classifier_SVM.svData2svmProblem: pProblem->x");
	if (positiveWeights != NULL && negativeWeights != NULL) {
		MArray_1D(pProblem->W, pProblem->l, double, "SC_Classifier_SVM.svData2svmProblem: pProblem->W");
	} else {
		pProblem->W = NULL;
	}

  for (long int i = 0; i < pProblem->l; i++) {
    if (i < pPositive->Row) {
      pProblem->y[i] = (double)(sclib::labelPositive);
      pProblem->x[i] = svDataRow2svmNode(pPositive, i);
			if (pProblem->W != NULL) {
				pProblem->W[i] = positiveWeights[i];
			}
    } else { //only hits if pNegative!=NULL
      pProblem->y[i] = (double)(sclib::labelNegative);
      pProblem->x[i] = svDataRow2svmNode(pNegative, i-pPositive->Row);
			if (pProblem->W != NULL) {
				pProblem->W[i] = negativeWeights[i-pPositive->Row];
			}
    }
  }

  return pProblem;
}

//====================================================================================================================
//	convert a single feature-vector out of an SV_Data container into a format that can be subject to SVM-testing
//  do the scaling too, if wished (according to class-mebers)
//====================================================================================================================
SC_SVMnode* SC_Classifier_SVM::svDataRow2svmNode(SV_Data *pData, unsigned long int row) {
  SC_SVMnode *pNodes = NULL;

  MArray_1D(pNodes, pData->Col + 1, SC_SVMnode, "SC_Classifier_SVM.floatVector2svmNode: pNodes");
  for (int x = 0; x < pData->Col; x++) {
    pNodes[x].index = x + 1;
    if (this->doScaling == true && this->pScale != NULL) { //perform feature-scaling according to the previously found parameters; don't do it in an extra function for efficiency purposes: here a copy has to be created anyway
      pNodes[x].value = (pData->Mat[row][x] - this->pScale->Mat[0][x]) / (this->pScale->Mat[1][x] - this->pScale->Mat[0][x]); //(x-min)/(max-min)
    } else {
      pNodes[x].value = pData->Mat[row][x];
    }
  }
  pNodes[pData->Col].index = -1; //this tells the SVM that the vector is finished

  return pNodes;
}

//====================================================================================================================
//	creates and returns a copy of one row of SVMnodes (terminated by a node with index -1)
//====================================================================================================================
SC_SVMnode* SC_Classifier_SVM::copySVMnodes(SC_SVMnode *pNodes) {
  unsigned long int i = 0, count = 0;
  SC_SVMnode *pNewNodes = NULL, *pHook = pNodes;

  if (pNodes != NULL) {
    while (pHook[i].index != -1) {
      i++;
      count++;
    }
    if (count > 0) {
      MArray_1D(pNewNodes, count + 1, SC_SVMnode, "SC_Classifier_SVM.copyNodes: pNewNodes");
      for (i = 0; i < count; i++) {
        pNewNodes[i].index = i + 1;
        pNewNodes[i].value = pNodes[i].value;
      }
      pNewNodes[count].index = -1; //this tells the SVM that the vector is finished
    }
  }

  return pNewNodes;
}

//====================================================================================================================
//	destroy a SC_SVMproblem struct
//====================================================================================================================
void SC_Classifier_SVM::killSVMproblem(SC_SVMproblem* &pProblem) {
  SC_SVMproblem *pToKill = pProblem;
	
	if (pToKill != NULL) {
    for (long int i = 0; i < pToKill->l; i++) {
      MFree_1D(pToKill->x[i]);
    }
    MFree_1D(pToKill->x);
    MFree_1D(pToKill->y);
		MFree_1D(pToKill->W);
    MFree_0D(pToKill);
		pProblem = NULL;
  }
  
  return;
}

//====================================================================================================================
//	create a copy of a SC_SVMproblem struct
//====================================================================================================================
SC_SVMproblem* SC_Classifier_SVM::copySVMproblem(SC_SVMproblem* pProblem) {
	int i, j, count;
	SC_SVMproblem *pCopy = NULL;
	
	if (pProblem != NULL) {
		pCopy =	new SC_SVMproblem;
		pCopy->l = pProblem->l;
		MArray_1D(pCopy->x, pCopy->l, SC_SVMnode*, "SC_Classifier_SVM.copySVMproblem: pCopy->x");
		MArray_1D(pCopy->y, pCopy->l, double, "SC_Classifier_SVM.copySVMproblem: pCopy->y");
		if (pProblem->W != NULL) {
			MArray_1D(pCopy->W, pCopy->l, double, "SC_Classifier_SVM.copySVMproblem: pCopy->W");
		} else {
			pCopy->W = NULL;
		}

		for (i = 0; i < pCopy->l; i++) {
			pCopy->y[i] = pProblem->y[i];
			j = 0;
			while (pProblem->x[i][j].index != -1) {
				j++;
			}
			count = j+1; //+1 for the terminating node with index==-1
			if (count > 0) {
				MArray_1D(pCopy->x[i], count, SC_SVMnode, "SC_Classifier_SVM.copySVMproblem: pCopy->x[i]");
				for (j = 0; j < count; j++) { 
					pCopy->x[i][j] = pProblem->x[i][j];
				}
			}
			if (pCopy->W != NULL) {
				pCopy->W[i] = pProblem->W[i];
			}
		}
	}

	return pCopy;
}

//====================================================================================================================
//  create a copy of a trained SVM model
//====================================================================================================================
SC_SVMmodel* SC_Classifier_SVM::copySVMmodel(SC_SVMmodel* pModel) {
	int i, j, count;
	SC_SVMmodel *pNewModel = NULL;
	SC_SVMnode *pHook;

	if (pModel != NULL) {
		pNewModel = (SC_SVMmodel*)(malloc(sizeof(SC_SVMmodel))); //need to use malloc() here 'cause the memory is freed inside SC_SVM.svm_destroy_model() using free()

		pNewModel->param.svm_type = pModel->param.svm_type;
		pNewModel->param.kernel_type = pModel->param.kernel_type;
		pNewModel->param.degree = pModel->param.degree;
		pNewModel->param.gamma = pModel->param.gamma;
		pNewModel->param.coef0 = pModel->param.coef0;
		pNewModel->param.cache_size = pModel->param.cache_size;
		pNewModel->param.eps = pModel->param.eps;
		pNewModel->param.C = pModel->param.C;
		pNewModel->param.nr_weight = pModel->param.nr_weight;
		if (pModel->param.nr_weight > 0) {
			pNewModel->param.weight_label = (int*)(malloc(pNewModel->param.nr_weight*sizeof(int))); 
			pNewModel->param.weight = (double*)(malloc(pNewModel->param.nr_weight*sizeof(double)));
			for (i = 0; i < pNewModel->param.nr_weight; i++) {
				pNewModel->param.weight_label[i] = pModel->param.weight_label[i];
				pNewModel->param.weight[i] = pModel->param.weight[i];
			}
		} else {
			pNewModel->param.weight_label = NULL;
			pNewModel->param.weight = NULL;
		}
		pNewModel->param.nu = pModel->param.nu;
		pNewModel->param.p = pModel->param.p;
		pNewModel->param.shrinking = pModel->param.shrinking;
		pNewModel->param.probability = pModel->param.probability;

		pNewModel->nr_class = pModel->nr_class;

		pNewModel->l = pModel->l;

		pNewModel->rho = (double*)(malloc((pNewModel->nr_class*(pNewModel->nr_class-1)/2)*sizeof(double)));
		for (i = 0; i < pNewModel->nr_class*(pNewModel->nr_class-1)/2; i++) {
			pNewModel->rho[i] = pModel->rho[i];
		}

		if (pModel->label != NULL) {
			pNewModel->label = (int*)(malloc(pNewModel->nr_class*sizeof(int)));
			for (i = 0; i < pNewModel->nr_class; i++) {
				pNewModel->label[i] = pModel->label[i];
			}
		} else {
			pNewModel->label = NULL;
		}

		if (pModel->probA != NULL) {
			pNewModel->probA = (double*)(malloc((pNewModel->nr_class*(pNewModel->nr_class-1)/2)*sizeof(double)));
			for (i = 0; i < pNewModel->nr_class*(pNewModel->nr_class-1)/2; i++) {
				pNewModel->probA[i] = pModel->probA[i];
			}
		} else {
			pNewModel->probA = NULL;
		}

		if (pModel->probB != NULL) {
			pNewModel->probB = (double*)(malloc((pNewModel->nr_class*(pNewModel->nr_class-1)/2)*sizeof(double)));
			for (i = 0; i < pNewModel->nr_class*(pNewModel->nr_class-1)/2; i++) {
				pNewModel->probB[i] = pModel->probB[i];
			}
		} else {
			pNewModel->probB = NULL;
		}

		if (pModel->nSV != NULL) {
			pNewModel->nSV = (int*)(malloc(pNewModel->nr_class*sizeof(int)));
			for (i = 0; i < pNewModel->nr_class; i++) {
				pNewModel->nSV[i] = pModel->nSV[i];
			}
		} else {
			pNewModel->nSV = NULL;
		}

		pNewModel->sv_coef = (double**)(malloc(pNewModel->nr_class*sizeof(double*)));
		for (j = 0; j < pNewModel->nr_class-1; j++) {
			pNewModel->sv_coef[j] = (double*)(malloc(pNewModel->l*sizeof(double)));
			for (i = 0; i < pNewModel->l; i++) {
				pNewModel->sv_coef[j][i] = pModel->sv_coef[j][i];
			}
		}

		//in libsvm, the SVs just link to the training data; we can't do this here 'cause even if we copied it,
		//we don't no which of the training vectors are the support vectors (and don't want to check for that)
		//so, we set the freeing flag below to something that we recognize and free it ourselves before calling
		//svm_destroy_model()
		pNewModel->SV = (SC_SVMnode**)(malloc(pNewModel->l*sizeof(SC_SVMnode*)));
		for (i = 0; i < pNewModel->l; i++) {
			j = 0;
			pHook = pModel->SV[i];
			while (pHook[j].index != -1) {
				j++;
			}
			count = j;
			if (count > 0) {
				pNewModel->SV[i] = (SC_SVMnode*)(malloc((count+1)*sizeof(SC_SVMnode)));
				if (pNewModel->SV[i] == NULL) {
					REPORT_ERROR(SVLIB_NoMem, "No memory while copying SVM model");
				}
				for (j = 0; j < count; j++) {
					pNewModel->SV[i][j].index = pModel->SV[i][j].index;
					pNewModel->SV[i][j].value = pModel->SV[i][j].value;
				}
				pNewModel->SV[i][count].index = -1; //this tells the SVM that the vector is finished
			}
		}

 		pNewModel->free_sv = 42; //something that we recognize in our destructor to do freeing ourselves because we don't just link to training vectors
	}

	return pNewModel;
}

//====================================================================================================================
// add or change a class-weight mapping
//====================================================================================================================
void SC_Classifier_SVM::setClassWeight(int classLabel, double weight) {
	SC_Classifier_SVM::SC_Weight *pHook = this->pClassWeightList;
	int count = 0, idx = -1;

	if (this->pClassWeightList == NULL) {
		this->pClassWeightList = new SC_Classifier_SVM::SC_Weight(classLabel, weight);
	} else {
		while (pHook != NULL) {
			if (pHook->label == classLabel)	{
				pHook->weight = weight;
				break;
			} else if (pHook->Next == NULL) { //classLabel wasn't found until here, so add a new mapping
				pHook->Next = new SC_Classifier_SVM::SC_Weight(classLabel, weight);
				break;
			}
			pHook = pHook->Next;
		}	
	}

	return;
}

//====================================================================================================================
// remove a single class-weight mapping
//====================================================================================================================
void SC_Classifier_SVM::removeClassWeight(int classLabel) {
	SC_Classifier_SVM::SC_Weight *pHook = this->pClassWeightList;
	int count = 0, idx = -1;

	while (pHook != NULL) {
		if (pHook->label == classLabel)	{
			idx = count;
			//no break here because we need "count" to be the length of the list, thus go until its end
		}
		count++;
		pHook = pHook->Next;
	}

	if (idx >= 0 && idx < count) {
		this->pClassWeightList = sclib::removeFromList(this->pClassWeightList, idx, count);
	}

	return;
}

//====================================================================================================================
// remove all class-weight mappings
//====================================================================================================================
void SC_Classifier_SVM::removeClassWeights(void) {
	sclib::destructLinkedList(this->pClassWeightList);

	return;
}

//====================================================================================================================
//	incorporates the pClassWeight knowledge into a svm-parameterset
//====================================================================================================================
void SC_Classifier_SVM::classWeights2parameters(SC_SVMparameter *pParam) {
	SC_Classifier_SVM::SC_Weight *pHook = this->pClassWeightList;
	int cnt = 0;

	if (this->pClassWeightList == NULL) {
		pParam->nr_weight = 0;
		pParam->weight_label = NULL;
		pParam->weight = NULL;
	} else {
		pParam->nr_weight = sclib::getListCount(this->pClassWeightList);

		pParam->weight = (double*)malloc(2*pParam->nr_weight*sizeof(double)); //need to use malloc() here 'cause the memory is freed inside SC_SVM.svm_destroy_param() using free()
		pParam->weight_label = (int*)malloc(2*pParam->nr_weight*sizeof(int));
		
		while (pHook != NULL) {
			pParam->weight_label[cnt] = pHook->label;
			pParam->weight[cnt] = pHook->weight;
			cnt++;
			pHook = pHook->Next;
		}
	}

	return;
}
