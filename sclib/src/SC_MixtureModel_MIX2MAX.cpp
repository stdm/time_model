/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a abnormal GMM-IB as accidentially */
/*				found by Wei Ho 'MIX2MAX' Tsai                                   */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.02.2006																								*/
/**************************************************************************/

#include <math.h>
#include <limits>
#include <assert.h> 
#include <float.h>
#include "SC_MixtureModel_MIX2MAX.h"
#include "SC_Aux.h"
#include "SC_Gauss.h"
#include <SV_Error.h>

#define MAX_ITER_NUM 1000 //TODO: by thilo: was 100
#define MAX_B_MIXTURE_NUM 64

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel_MIX2MAX::SC_MixtureModel_MIX2MAX(SC_TweakableParameters* pTweak, SC_MixtureModel* pBackground, unsigned short int mixtureCount, unsigned short int dim) : SC_MixtureModel(pTweak) {
	this->pBackground = pBackground;
  this->pOriginalBackground	= NULL;
	this->mixtureCount = mixtureCount; 
	this->dim = dim;
  this->maxEMsteps = this->pTweak->mixtureModelMix2Max.maxEMiterations;

	if (this->pBackground != NULL) {
    this->pOriginalBackground	= new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pBackground));
	  assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
	}

	this->Hdr.ModelType = sclib::mtMIX2MAX;

	//normally, the mixtureCount must be >0, but there is one exception from this rule:
	//maybe we need only a "link" in the linked list of models, but have no data to model.
	//then, a model without any mixtures represents such a link without modelling anything.
	if (this->mixtureCount > 0) {
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX: sd");
		MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX: mean");
		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX: weight");
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX: maskLevel");

		for (int x = 0; x < this->mixtureCount; x++) {
      for (int y = 0; y < this->dim; y++) {
			  this->variance[x][y] = 0.0;
        this->sd[x][y] = 0.0;
			  this->mean[x][y] = 0.0;
      }
			this->weight[x] = 0.0;
      this->maskLevel[x] = 0.0;
		}
	} else {
		this->variance = NULL;
    this->sd = NULL;
		this->mean = NULL;
		this->weight =	NULL;
    this->maskLevel = NULL;
	}
}

//====================================================================================================================
// copy-constructor
//====================================================================================================================
SC_MixtureModel_MIX2MAX::SC_MixtureModel_MIX2MAX(const SC_MixtureModel_MIX2MAX& pParent) : SC_MixtureModel(pParent) {
	this->pBackground = pParent.pBackground;
  this->pOriginalBackground	= (pParent.pOriginalBackground != NULL) ? new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pParent.pOriginalBackground)) : NULL;

	if (this->pBackground != NULL) {
    assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
	}
}

//====================================================================================================================
// destructor
// attention: doesn't delete the background-gmm, nore the linked gmm-ib's!
// but: does delete the original background gmm, because it's a copy, not a reference/pointer!
//====================================================================================================================
SC_MixtureModel_MIX2MAX::~SC_MixtureModel_MIX2MAX() {
  this->pBackground = NULL;

  MFree_0D(this->pOriginalBackground);
  MFree_1D(this->maskLevel);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_MixtureModel_MIX2MAX& SC_MixtureModel_MIX2MAX::operator=(const SC_MixtureModel_MIX2MAX& pParent) {
	if (this != &pParent) {
		this->SC_MixtureModel::operator=(pParent);

		MFree_1D(this->maskLevel);
		MFree_0D(this->pOriginalBackground);
	  
		this->pBackground = pParent.pBackground;
		this->pOriginalBackground	= (pParent.pOriginalBackground != NULL) ? new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pParent.pOriginalBackground)) : NULL;

		if (this->pBackground != NULL) {
			assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
		}

		if (this->mixtureCount > 0) {
			MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_GMM: maskLevel");
			for (int x = 0; x < this->mixtureCount; x++) {
				this->maskLevel[x] = pParent.maskLevel[x];
			}
		}
	}
  
  return *this;
}

//====================================================================================================================
// Inherited version of training-algorithm; just calls the new one (below), which is capable of handling many differnt
// noise-models for different parts of the feature-vectors
//====================================================================================================================
int SC_MixtureModel_MIX2MAX::TrainModel(SV_Data* TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Trains the model, using segmentsToMerge segments from the linked list of features; evaluates the 
// noiseCorruptionType and chooses the suitable algorithm
//====================================================================================================================
int SC_MixtureModel_MIX2MAX::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  int res = reestimation(pData, segmentsToMerge);

  return res;
}

//train the model using the data
//has been optimized using the knowledge that using this set of formulas, only the weights get adapted
//also, own initialization and gauss/erf-computation is used now, together with possible background-model switching
//but is qualitatively the same as wesleys original version
int SC_MixtureModel_MIX2MAX::reestimation(SV_Data *pData, unsigned long int segmentsToMerge)
{
	int		i,j,k,t,m,EffectiveLength;
	float	***BSExpectation,***BSExpectation2;
	double	**BSJointProb;
	double	*SCD,*SPDF,BCD,BPDF,BSCondProb,BSJointProb_sum,MaxLogP; //,temp,B_SCondProb, E_xlz, E_xlz2, tmp; //also commented because not relevant for updating the weights...
	double	OverallScore[MAX_ITER_NUM];
	double	**Mean_numerator,**Var_numerator,*denominator;

  //by thilo
  double scaledX, **logSd, **logBgSd; //, E_xlz, E_xlz2;
  unsigned long int idx;
  long actualT	= pData->Row; //count of feature-vectors in the actual elements of the linked list
	short J = sclib::max(1, this->pBackground->getMaxMixturesInList());
	SV_Data	*pActualData = pData;	//pointer to actual element of linked list
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	SC_MixtureModel *pBackground	= this->pBackground; //don't touch the original pointer during training!
  MArray_2D(logSd, this->mixtureCount, this->dim, double, "logSd");
  MArray_2D(logBgSd, J, pBackground->getDim(), double, "logBgSd");

	if((BSJointProb=(double **)malloc(this->mixtureCount*sizeof(double *)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
	if((BSExpectation=(float ***)malloc(this->mixtureCount*sizeof(float **)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
	if((BSExpectation2=(float ***)malloc(this->mixtureCount*sizeof(float **)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
	for(i=0;i<this->mixtureCount;i++){
		if((BSJointProb[i]=(double *)malloc(pBackground->getMixtureCount()*sizeof(double)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
		if((BSExpectation[i]=(float **)malloc(pBackground->getMixtureCount()*sizeof(float *)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
		if((BSExpectation2[i]=(float **)malloc(pBackground->getMixtureCount()*sizeof(float *)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
		for(j=0;j<pBackground->getMixtureCount();j++){
			if((BSExpectation[i][j]=(float *)malloc(pCompleteData->Col*sizeof(float)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
			if((BSExpectation2[i][j]=(float *)malloc(pCompleteData->Col*sizeof(float)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
		}
	} 
	denominator=(double *)malloc(this->mixtureCount*sizeof(double));
	Mean_numerator=(double **)malloc(this->mixtureCount*sizeof(double *));
	Var_numerator=(double **)malloc(this->mixtureCount*sizeof(double *));
	for(i=0;i<this->mixtureCount;i++){
		Mean_numerator[i]=(double *)malloc(pCompleteData->Col*sizeof(double));
		Var_numerator[i]=(double *)malloc(pCompleteData->Col*sizeof(double));
	}
	SCD=(double *)malloc(pCompleteData->Col*sizeof(double));
	SPDF=(double *)malloc(pCompleteData->Col*sizeof(double));

  //by thilo
  initParameters(pCompleteData, this->pTweak->mixtureModelMix2Max.varianceLimit, this->pTweak->mixtureModelMix2Max.kMeansIterations); //by thilo
  for(m=0;m<pCompleteData->Col;m++){
    for(i=0;i<this->mixtureCount;i++){
		  logSd[i][m] = sclib::sLog(this->sd[i][m]);
    }
  }
	if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIX2MAX.txt", this, this->pTweak);}

	// maximum OverallScore search
  for(k=0;k<(long)this->maxEMsteps;k++){
		OverallScore[k]=0.0;			
		for(i=0;i<this->mixtureCount;i++){
			denominator[i]=0.0;
			for(m=0;m<pCompleteData->Col;m++)	Mean_numerator[i][m]=Var_numerator[i][m]=0.0;
      this->maskLevel[i] = 0.0; //by thilo //TODO: compute masklevel
		}
		EffectiveLength=0;

    //by thilo: for bg-model switching
    actualT	= pData->Row;
	  pActualData = pData;
	  pBackground	= this->pBackground;
    J = sclib::max(1, pBackground->getMixtureCount());

    for(t=0;t<pCompleteData->Row;t++){

      //by thilo: bg-model switching and sd-initialization
      for (m = 0; m < pCompleteData->Col; m++) {
        for(j=0;j<pBackground->getMixtureCount();j++){
		      logBgSd[j][m] = sclib::sLog(pBackground->getSd(j, m));
        }
      }

			//by thilo: handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
      //          init background-sd
		  if (t >= actualT) {
        assert(pActualData->Next != NULL);
        assert(pBackground->Next != NULL);
			  pActualData	= pActualData->Next;
			  actualT += pActualData->Row;
			  pBackground	= (SC_MixtureModel*)(pBackground->Next);
        J = sclib::max(1, pBackground->getMixtureCount());
		  }

      BSJointProb_sum=0.0; //by thilo
      MaxLogP=pCompleteData->Col*-705.0; //by thilo; was -10000000000000.0;
			for(i=0;i<this->mixtureCount;i++){
				for(m=0;m<pCompleteData->Col;m++){
					scaledX = sclib::zTransform(pCompleteData->Mat[t][m], this->mean[i][m], this->sd[i][m]);
          idx = this->gaussSolver.getIdx(scaledX);
          this->gaussSolver.tabledLogGaussianAndErf(idx, logSd[i][m], SPDF[m], SCD[m]);
				}
				for(j=0;j<J;j++){
          BSJointProb[i][j]=0.0;
					for(m=0;m<pCompleteData->Col;m++){	
						scaledX = sclib::zTransform(pCompleteData->Mat[t][m], pBackground->getMean(j, m), pBackground->getSd(j, m));
            idx = this->gaussSolver.getIdx(scaledX);
            this->gaussSolver.tabledLogGaussianAndErf(idx, logBgSd[j][m], BPDF, BCD);
						BSCondProb=BCD*exp(SPDF[m])+SCD[m]*exp(BPDF);
						//temp=BPDF-SPDF[m]; //also commented because not relevant for updating the weights...
            //the below block never gets evaluated (oh, it get's, but the results are immediately overwritten in the next block due to missing else)
            //in wesleys version, so it's commented for speed purposes here...
						/*if((fabs(temp)<744.0)&&(BCD!=0.0)){
							B_SCondProb=1.0/(1.0+(SCD[m]/BCD)*exp(temp));
              tmp = exp(SPDF[m]) / SCD[m];
              if (sclib::isFinite(tmp) == true) {
                E_xlz = this->mean[i][m] - this->variance[i][m]*tmp;
                E_xlz2 = this->mean[i][m]*this->mean[i][m] + this->variance[i][m] - (this->mean[i][m]+(double)(pCompleteData->Mat[t][m]))*this->variance[i][m]*(tmp);
              } else {
                E_xlz = pCompleteData->Mat[t][m];
                E_xlz2 = this->variance[i][m] + (pCompleteData->Mat[t][m] * pCompleteData->Mat[t][m]);
              }
              BSExpectation[i][j][m]=(float)(B_SCondProb*(double)(pCompleteData->Mat[t][m])+(1.0-B_SCondProb)*(E_xlz));
              BSExpectation2[i][j][m]=(float)((B_SCondProb*(double)(pCompleteData->Mat[t][m])*(double)(pCompleteData->Mat[t][m]))+(1.0-B_SCondProb)*(E_xlz2));
							//BSExpectation[i][j][m]=(float)(B_SCondProb*(double)(pCompleteData->Mat[t][m])+(1.0-B_SCondProb)*(this->mean[i][m]-this->variance[i][m]*exp(SPDF[m])/SCD[m]));
							//BSExpectation2[i][j][m]=(float)((B_SCondProb*(double)(pCompleteData->Mat[t][m])*(double)(pCompleteData->Mat[t][m]))+(1.0-B_SCondProb)*(this->mean[i][m]*this->mean[i][m]+this->variance[i][m]-(this->mean[i][m]+(double)(pCompleteData->Mat[t][m]))*this->variance[i][m]*(exp(SPDF[m])/SCD[m])));
               this->maskLevel[i] += 1.0 - B_SCondProb;
						}*//* if((temp<-744.0)&&(BCD!=0.0)){	//	BSCondProb=1.0;
							BSExpectation[i][j][m]=pCompleteData->Mat[t][m];
							BSExpectation2[i][j][m]=(pCompleteData->Mat[t][m])*(pCompleteData->Mat[t][m]);
              this->maskLevel[i] += 0.0;
						}else{					// BSCondProb=0.0;
							BSExpectation[i][j][m]=(float)(this->mean[i][m]);
							BSExpectation2[i][j][m]=(float)(this->mean[i][m]*this->mean[i][m]+this->variance[i][m]);
              this->maskLevel[i] += 1.0;
						}*/ //also commented because not relevant for updating the weights...
						BSJointProb[i][j] += sclib::sLog(BSCondProb);
					}
					if(BSJointProb[i][j]>MaxLogP) MaxLogP=BSJointProb[i][j]; 
				}
			}

      BSJointProb_sum=0.0;
			for(i=0;i<this->mixtureCount;i++){
				for(j=0;j<J;j++){
          if((BSJointProb[i][j]-MaxLogP)>-744.0) BSJointProb[i][j]=this->weight[i]*pBackground->getWeight(j)*exp(BSJointProb[i][j]-MaxLogP);
          else	BSJointProb[i][j]=0.0;
					BSJointProb_sum +=BSJointProb[i][j];
        }
			}
      if(BSJointProb_sum!=0.0){
				OverallScore[k] +=log(BSJointProb_sum)+MaxLogP;
				for(i=0;i<this->mixtureCount;i++){
					for(j=0;j<J;j++){
            BSJointProb[i][j] /=BSJointProb_sum; //by thilo (was extra loop above)
						denominator[i] +=BSJointProb[i][j];
						/*for(m=0;m<pCompleteData->Col;m++){
							Mean_numerator[i][m] +=BSJointProb[i][j]*(double)(BSExpectation[i][j][m]);
							Var_numerator[i][m] +=BSJointProb[i][j]*(double)(BSExpectation2[i][j][m]);
						}*/ //also commented because not relevant for updating the weights...
					}
				}
				EffectiveLength++;
			}
		}

		for(i=0;i<this->mixtureCount;i++){
			/*for(m=0;m<pCompleteData->Col;m++){
				this->mean[i][m]=Mean_numerator[i][m]/denominator[i];
				this->variance[i][m]=Var_numerator[i][m]/denominator[i]-this->mean[i][m]*this->mean[i][m];
				this->sd[i][m]=sqrt(this->variance[i][m]);
        logSd[i][m] = sclib::sLog(this->sd[i][m]);
			}*/ //also commented because not relevant for updating the weights...
			this->weight[i]=denominator[i]/(double)(EffectiveLength);
      this->maskLevel[i] /= (double)(pCompleteData->Row*pCompleteData->Col*pBackground->getMixtureCount());
		}
    //printf("\n%ld: %lf \t length-diff:%d",k,OverallScore[k],pCompleteData->Row-EffectiveLength);

		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIX2MAX.txt", this, this->pTweak);}
		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::scalarOut("MIX2MAX_likelihood.txt", OverallScore[k], this->pTweak);}

		if( (OverallScore[k]>=OverallScore[k-1])&&(OverallScore[k]-OverallScore[k-1]<this->pTweak->mixtureModelMix2Max.EMthreshold)) break;
		if( (k>10)&&(OverallScore[k]<=OverallScore[k-1]) ) break;
	}

  //block by thilo:
  this->trainingDataCount = pCompleteData->Row;
	MFree_0D(this->pOriginalBackground);
  this->pOriginalBackground =  new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)this->pBackground));
  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

  J = sclib::max(1, this->pBackground->getMaxMixturesInList());
	for(i=0;i<this->mixtureCount;i++){
		for(j=0;j<J;j++){
			free(BSExpectation[i][j]);
			free(BSExpectation2[i][j]);
		}
		free(BSJointProb[i]);
		free(BSExpectation[i]);
		free(BSExpectation2[i]);
	} 
	free(BSJointProb);
	free(BSExpectation);
	free(BSExpectation2);
	for(i=0;i<this->mixtureCount;i++){
		free(Mean_numerator[i]);
		free(Var_numerator[i]);
	}
	free(denominator);
	free(Mean_numerator);
	free(Var_numerator);
	free(SCD);
	free(SPDF);

  //by thilo
  MFree_2D(logSd);
  MFree_2D(logBgSd);

  return k;
}

//====================================================================================================================
// Test model while considering only one element of the linked list
// return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_MIX2MAX::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;					
}

//====================================================================================================================
// Test GMM-IB with a specified nr. of segments out of the linked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_MIX2MAX::TestModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  SV_Data *pScore = cal_likelihood(pData, segmentsToMerge);

  return pScore;
}

//test the model using the data
//should be equivalent to wesleys original algorithm...
//input is a model as saved by SC_Lib (i.e.: this->sd is really the standard deviation, not MIX2MAXs varConst...)
//own gauss/erf-computation is used
SV_Data* SC_MixtureModel_MIX2MAX::cal_likelihood(SV_Data* pData, unsigned long int segmentsToMerge) 
{
	int t,i,j;
	double likelihood,BPDF[MAX_B_MIXTURE_NUM],BCD[MAX_B_MIXTURE_NUM],SPDF,SCD,Pr;
	double Blikelihood,Slikelihood;

	SV_Data *pScore	= new SV_Data(1, 1);
  int	PatternNum = 20;
	int	FeatDim = this->dim;
	int FeatLen;

  //by thilo
  unsigned short int J;
	long int actualT = pData->Row; //count of feature-vectors in the actual elements of the linked list
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);	//merge selected list-elements in a linked list together
	SV_Data *pActualData = pData; //pointer to actual element of linked list
  SC_MixtureModel *pBackground;
  if (this->pTweak->mixtureModelMix2Max.bgModelCombination == true) {
    pBackground = (SC_MixtureModel*)(this->pBackground->combineModels(this->pBackground, this->pOriginalBackground, true)); //don't touch the original pointer during training!
  } else {
    pBackground = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)(this->pBackground)));
  }
  SC_MixtureModel *pBgHook = NULL;
  J = sclib::max(pBackground->getMixtureCount(), 1);

	likelihood=0.0;
  FeatLen = pCompleteData->Row;

  for(t=0;t<FeatLen;t++){

    //by thilo: handle background-model-switching and errors due to the linked-list-nature of the feature-vectors
    if (t >= actualT) { 
      assert(pActualData->Next != NULL);
      assert(pBackground->Next != NULL);
			pActualData	= pActualData->Next;
			actualT += pActualData->Row;
      pBgHook = (SC_MixtureModel*)(pBackground->Next);
      MFree_0D(pBackground);
      if (this->pTweak->mixtureModelMix2Max.bgModelCombination == true) {
        pBackground = (SC_MixtureModel*)(this->pBackground->combineModels((SC_MixtureModel*)(pBgHook), this->pOriginalBackground, true));
      } else {
        pBackground = new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)(pBgHook)));
      }
      J = sclib::max(1, pBackground->getMixtureCount());
		}

		Blikelihood=0.0;
		for(j=0;j<J;j++){
      BPDF[j] = this->gaussSolver.multivariateGaussian(pCompleteData->Mat[t], this->pBackground->getMean(j), this->pBackground->getVariance(j), this->pBackground->getSd(j), this->dim);
			BCD[j] = this->gaussSolver.multivariateAverageErf(pCompleteData->Mat[t], this->pBackground->getMean(j), this->pBackground->getSd(j), this->dim);
			Blikelihood +=BPDF[j];
		}

		Pr=0.0;
		Slikelihood=0.0;
		for(i=0;i<this->mixtureCount;i++){
      SPDF = this->gaussSolver.multivariateGaussian(pCompleteData->Mat[t], this->mean[i], this->variance[i], this->sd[i], this->dim);
			SCD = this->gaussSolver.multivariateAverageErf(pCompleteData->Mat[t], this->mean[i], this->sd[i], this->dim);
      Slikelihood +=SPDF;

      for(j=0;j<J;j++) {
				Pr +=this->weight[i]*this->pBackground->getWeight(j)*(SPDF*BCD[j]+BPDF[j]*SCD);
      }
		}

		likelihood += sclib::sLog(Pr);
	}	

	pScore->Mat[0][0] = (float)(likelihood) / (float)(FeatLen); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...;
  MFree_0D(pBackground); //by thilo
  if (pData != pCompleteData) { //by thilo
    MFree_0D(pCompleteData);
  }

	return pScore;
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_MIX2MAX::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_MixtureModel*)pSecond, false));
}
SC_MixtureModel_MIX2MAX* SC_MixtureModel_MIX2MAX::combineModels(SC_MixtureModel* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model*  SC_MixtureModel_MIX2MAX::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
   return (SC_Model*)(combineModels((SC_MixtureModel_MIX2MAX*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_MIX2MAX* SC_MixtureModel_MIX2MAX::combineModels(SC_MixtureModel_MIX2MAX* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext) {
  unsigned short int i, d, newMixtureCount = pFirst->getMixtureCount() + pSecond->getMixtureCount();
	unsigned long int newTrainingDataCount = pFirst->getTrainingDataCount() + pSecond->getTrainingDataCount();
	double *newWeight, *newMaskLevel, **newMean, **newVariance, **newSd;
	SC_MixtureModel_MIX2MAX *pNewModel;
  SC_MixtureModel *pNewOriginalBackground;
  bool tooSmall = false;
  
	if (pFirst->getDim() != pSecond->getDim()) {REPORT_ERROR(SVLIB_BadData, "Can't combine Models trained on different feature-dim!");}

	MArray_1D(newWeight, newMixtureCount, double, "SC_MixtureModel_MIX2MAX.combineModels: newWeight");
  MArray_1D(newMaskLevel, newMixtureCount, double, "SC_MixtureModel_MIX2MAX.combineModels: newMaskLevel");
	MArray_2D(newMean, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIX2MAX.combineModels: newMean");
	MArray_2D(newVariance, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIX2MAX.combineModels: newVariance");
  MArray_2D(newSd, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIX2MAX.combineModels: newSd");

  pNewModel = new SC_MixtureModel_MIX2MAX(pFirst->getTweak(), NULL, newMixtureCount, pFirst->getDim());

	//copy the mixture-components together
	for (i = 0; i < newMixtureCount; i++) {
		if (i < pFirst->getMixtureCount()) {
			newWeight[i] = (pFirst->getWeight(i) * pFirst->getTrainingDataCount()) / newTrainingDataCount; //re-normalize the weights so that their sum = 1.0
      newMaskLevel[i] = pFirst->getMaskLevel(i);
      for (d = 0; d < pFirst->getDim(); d++) {
				newMean[i][d] = pFirst->getMean(i, d);
				newVariance[i][d] = pFirst->getVariance(i, d);
        newSd[i][d] = pFirst->getSd(i, d);
			}
		} else {
			newWeight[i] = (pSecond->getWeight(i - pFirst->getMixtureCount()) * pSecond->getTrainingDataCount()) / newTrainingDataCount; //re-normalize the weights so that their sum = 1.0
      if (pSecond->Hdr.ModelType == sclib::mtMIXMAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX_ex) {
        newMaskLevel[i] = ((SC_MixtureModel_MIX2MAX*)pSecond)->getMaskLevel(i - pFirst->getMixtureCount());
      } else {
        newMaskLevel[i] = 0.0; //no gmm-ib => no noise-masking
      }
			for (d = 0; d < pFirst->getDim(); d++) {
				newMean[i][d] = pSecond->getMean(i - pFirst->getMixtureCount(), d);
				newVariance[i][d] = pSecond->getVariance(i - pFirst->getMixtureCount(), d);
        newSd[i][d] = pSecond->getSd(i - pFirst->getMixtureCount(), d);
			}
		}
    if (newWeight[i] < pFirst->getTweak()->mixtureModelMix2Max.weightLimit) {tooSmall = true;}
	}

	//make the new mixture-components the new parameters of the new model
  pNewModel->setWeight(newWeight);
  pNewModel->setMaskLevel(newMaskLevel);
  pNewModel->setMean(newMean);
  pNewModel->setVariance(newVariance);
  pNewModel->setSd(newSd);
  pNewModel->setTrainindDataCount(newTrainingDataCount);
  pNewModel->setMixtureCount(newMixtureCount);
  if (pSecond->Hdr.ModelType == sclib::mtMIXMAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX || pSecond->Hdr.ModelType == sclib::mtMIX2MAX_ex) {
    pNewOriginalBackground = (SC_MixtureModel*)(pFirst->getOriginalBackground()->combineModels(((SC_MixtureModel_MIX2MAX*)pSecond)->getOriginalBackground()));
    pNewModel->setOriginalBackground(pNewOriginalBackground);
  } else {
    SC_MixtureModel_GMM *pNewDummy = new SC_MixtureModel_GMM(pFirst->getTweak(), 0, pFirst->getDim());
    pNewOriginalBackground = (SC_MixtureModel*)(pNewDummy->combineModels(pFirst->getOriginalBackground()));
    pNewModel->setOriginalBackground(pNewOriginalBackground);
    MFree_0D(pNewDummy);
  }

  //kill too small mixtures
  if (tooSmall == true) {
    pNewModel->sortParameters();
    for (i = 0; i < pNewModel->getMixtureCount(); i++) {
      if (pNewModel->getWeight(i) < pFirst->getTweak()->mixtureModelMix2Max.weightLimit) {
        pNewModel->killMixture(i);
        i--;
      }
    }
  }

  if (keepFirstsNext == true) {
    pNewModel->Next = pFirst->Next;
  } else {
    pNewModel->Next = NULL;
  }

	return pNewModel;
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model; 
// Convert a dummy mixture (no mean, no weight, no variance) to a real mixture with equal result
//====================================================================================================================
SC_Model* SC_MixtureModel_MIX2MAX::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepDummyMixture, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel_MIX2MAX*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_MIX2MAX* SC_MixtureModel_MIX2MAX::combineModels(SC_MixtureModel_MIX2MAX* pFirst, SC_MixtureModel* pSecond, bool keepDummyMixture, bool keepFirstsNext) {
  assert("Don't use this method!!!" == 0); //TODO: why???
  
  return combineModels(pFirst, pSecond, keepFirstsNext);
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_MixtureModel_MIX2MAX::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pSecond, pSpeechFrames, segmentsToMerge, pBackgroundModels));
}
SC_MixtureModel_MIX2MAX* SC_MixtureModel_MIX2MAX::combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_MixtureModel* pBackgroundModels) {
	SC_MixtureModel_MIX2MAX* pNewModel;
	
  assert("Don't use this method!!!" == 0); //TODO: why???
  assert(this->dim == pSecond->getDim());
	
	pNewModel = new SC_MixtureModel_MIX2MAX(this->pTweak, pBackgroundModels, sclib::min(this->mixtureCount+pSecond->getMixtureCount(), this->pTweak->modelHandler.maxSpeakerModelOrder), this->dim);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

	return pNewModel;
}

//====================================================================================================================
// Sort the model-parameters according to the weights
//====================================================================================================================
void SC_MixtureModel_MIX2MAX::sortParameters() {
	unsigned short int i, d;
	SV_Data *temp;
	double *oldWeight, *oldMaskLevel, **oldMean, **oldVariance, **oldSd;
	
	temp = new SV_Data(this->mixtureCount, 2);
	MArray_1D(oldWeight, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX.sortParameters: oldWeight");
	MArray_1D(oldMaskLevel, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX.sortParameters: oldMaskLevel");
	MArray_2D(oldMean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX.sortParameters: oldMean");
	MArray_2D(oldVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX.sortParameters: oldVariance");
  MArray_2D(oldSd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX.sortParameters: oldVariance");

	for (i = 0; i < this->mixtureCount; i++) {
		temp->Mat[i][0] = (float)(this->weight[i]);
		temp->Mat[i][1] = (float)(i); //remember the old index in the 2nd column of the data-matrix
		oldWeight[i] = this->weight[i];
    oldMaskLevel[i] = this->maskLevel[i];
		for (d = 0; d < this->dim; d++) {
			oldMean[i][d] = this->mean[i][d];
			oldVariance[i][d] = this->variance[i][d];
      oldSd[i][d] = this->sd[i][d];
		}
	}

	//sort the list of weights; remember the old index in the 2nd column of the data-matrix
	//sclib::quickSort(temp, 0, this->mixtureCount - 1, 2, 0);
  sclib::quickSort(temp->Mat, 0, this->mixtureCount - 1, 2, 0);

	for (i = 0; i < this->mixtureCount; i++) {
		this->weight[i] = oldWeight[(int)(temp->Mat[i][1])];
    this->maskLevel[i] = oldMaskLevel[(int)(temp->Mat[i][1])];
		for (d = 0; d < this->dim; d++) {
			this->mean[i][d] = oldMean[(int)(temp->Mat[i][1])][d];
			this->variance[i][d] = oldVariance[(int)(temp->Mat[i][1])][d];
      this->sd[i][d] = oldSd[(int)(temp->Mat[i][1])][d];
		}
	}

	MFree_0D(temp);
	MFree_1D(oldWeight);
  MFree_1D(oldMaskLevel);
	MFree_2D(oldMean);
	MFree_2D(oldVariance);
  MFree_2D(oldSd);

	return;
}

//====================================================================================================================
// Remove the indexed mixture-component from the models mean/variance/weight-vectors/matrices
//====================================================================================================================
void SC_MixtureModel_MIX2MAX::killMixture(unsigned short int mixture) {
  unsigned short int i, i_idx = 0, d,  I = this->mixtureCount, D = this->dim;
  double *newWeight, *newMaskLevel, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I-1, double, "SC_MixtureModel_MIX2MAX.killMixture: newWeight");
  MArray_1D(newMaskLevel, I-1, double, "SC_MixtureModel_MIX2MAX.killMixture: newMaskLevel");
  MArray_2D(newVariance, I-1, D, double, "SC_MixtureModel_MIX2MAX.killMixture: newVariance");
  MArray_2D(newMean, I-1, D, double, "SC_MixtureModel_MIX2MAX.killMixture: newMean");
  MArray_2D(newSd, I-1, D, double, "SC_MixtureModel_MIX2MAX.killMixture: newSd");

  for (i = 0; i < I; i++) { //copy mixtures, leave the one to kill out
    if (i != mixture) {
      newWeight[i_idx] = this->weight[i] * (1.0  / (1.0 - this->weight[mixture])); //renormalize weight
      newMaskLevel[i_idx] = this->maskLevel[i];
      for (d = 0; d < D; d++) {
        newVariance[i_idx][d] = this->variance[i][d];
        newMean[i_idx][d] = this->mean[i][d];
        newSd[i_idx][d] = this->sd[i][d];
      }
      i_idx++;
    }
  }

  MFree_1D(this->weight);
  MFree_1D(this->maskLevel);
  MFree_2D(this->variance);
  MFree_2D(this->mean);
  MFree_2D(this->sd);

  this->weight = newWeight;
  this->maskLevel = newMaskLevel;
  this->variance = newVariance;
  this->mean = newMean;
  this->sd = newSd;
  this->mixtureCount--;

  return;
}

//====================================================================================================================
// Add one new mixture component with all means=0, variances=varianceLimit, weight=1/mixtureCount
//====================================================================================================================
void SC_MixtureModel_MIX2MAX::addMixture(void) {
  unsigned short int i, d, I = this->mixtureCount, D = this->dim;
  double *newWeight, *newMaskLevel, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I+1, double, "SC_MixtureModel_MIX2MAX.addMixture: newWeight");
  MArray_1D(newMaskLevel, I+1, double, "SC_MixtureModel_MIX2MAX.addMixture: newMaskLevel");
  MArray_2D(newVariance, I+1, D, double, "SC_MixtureModel_MIX2MAX.addMixture: newVariance");
  MArray_2D(newMean, I+1, D, double, "SC_MixtureModel_MIX2MAX.addMixture: newMean");
  MArray_2D(newSd, I+1, D, double, "SC_MixtureModel_MIX2MAX.addMixture: newSd");

  //add new mixture
  newWeight[I] = 1.0 / (this->mixtureCount+1.0);
  newMaskLevel[I] = 0.0;
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
  MFree_1D(this->maskLevel);
  MFree_2D(this->variance);
  MFree_2D(this->mean);
  MFree_2D(this->sd);

  this->weight = newWeight;
  this->maskLevel = newMaskLevel;
  this->variance = newVariance;
  this->mean = newMean;
  this->sd = newSd;
  this->mixtureCount++;

  return;
}

//====================================================================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
//====================================================================================================================
int SC_MixtureModel_MIX2MAX::SaveModel(void) {
 	int res, x, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAX model Failed!");}

  //write mixtureCount and featureDim
	bytes = io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAX Model Failed!");}

  //write weight-vector
	bytes += io.writeArray(&(this->DFile), this->weight, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAX Model Failed!");}

  //write maskLevel-vector
	bytes += io.writeArray(&(this->DFile), this->maskLevel, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAX Model Failed!");}

	//write mean and variance matrices
	for (x = 0; x < this->mixtureCount; x++) {
		bytes += io.writeArray(&(this->DFile), this->mean[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->variance[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->sd[x], this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAX Model Failed!");}
	}

  //write trainingDataCount
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAX Model Failed!");}

  //write maxEMstep
	io.writeScalar(&(this->DFile), this->maxEMsteps);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAX Failed!");}

  return bytes + MHLen; //MHLen is just an estimate of the written header-bytes...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model* SC_MixtureModel_MIX2MAX::LoadModel(void) {
	int res, x;
	unsigned short int newDim = 0, newMixtureCount = 0;
  SC_MixtureModel_GMM *pDummy = NULL;
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
	if ((DFile.good() != TRUE) || (newMixtureCount == 0)) {return(NULL);}

	//prepare this model for changing it's parameters
	if ((newDim != this->dim) || (newMixtureCount != this->mixtureCount)) {
		MFree_1D(this->weight);
    MFree_1D(this->maskLevel);
		MFree_2D(this->mean);
		MFree_2D(this->variance);
    MFree_2D(this->sd);

		this->dim = newDim;
		this->mixtureCount = newMixtureCount;

 		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX: weight");
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX: maskLevel");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX: mean");
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAX: sd");
	};

	//read weight, maskLevel, mean and variance
	io.readArray(&(this->DFile), this->weight, this->mixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAX Failed!");}
	io.readArray(&(this->DFile), this->maskLevel, this->mixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAX Failed!");}
	for (x = 0; x < this->mixtureCount; x++) {
		io.readArray(&(this->DFile), this->mean[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->variance[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->sd[x], this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAX Failed!");}
	}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAX Failed!");}

  //Set a dummy as the original background
  pDummy = new SC_MixtureModel_GMM(this->pTweak, 0, this->dim);
  this->pOriginalBackground = pDummy;

  //read maxEMstep
	io.readScalar(&(this->DFile), this->maxEMsteps, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAX Failed!");}

	return(this);
}

//====================================================================================================================
// create a (linked, not copied) SC_Signature-view on this model (for distance computation )
//====================================================================================================================
SC_Signature* SC_MixtureModel_MIX2MAX::toSignature(void) {
  double *weights;
  SC_Centroid **centroids = NULL;
  SC_Signature *pSignature;
  
  MArray_1D(centroids, this->mixtureCount, SC_Centroid*, "SC_MixtureModel_MIX2MAX.toSignature: centroids");
  MArray_1D(weights, this->mixtureCount, double, "SC_MixtureModel_MIX2MAX.toSignature: weights");

  for (unsigned short int i = 0; i < this->mixtureCount; i++) {
    weights[i] = this->weight[i] * (1.0 - this->maskLevel[i]);
		centroids[i] = new SC_Centroid_Gaussian(this->pTweak, this->dim, this->mean[i], this->variance[i]);
  }

	pSignature = new SC_Signature(centroids, weights, this->mixtureCount);

  return pSignature;
}

//====================================================================================================================
// destruct a signature created by this model's toSignature()-method
//====================================================================================================================
void SC_MixtureModel_MIX2MAX::killSignature(SC_Signature *pSignature) {
  if (pSignature != NULL) {
    MFree_0D(pSignature);
  }

  return;
}

/*

//the below functions could be replaced with the corresponding SC_Lib-versions without changing the algorithms outcome substantially...

//initialize model-parameters via k-keans algorithm
void SC_MixtureModel_MIX2MAX::K_means(SV_Data* pData)
{
	int *SampleNum,*BelongTo;
	int ClusterNo,SplitWho,TotalSampleNum;
	int	i,t,j,k;
	double *dis,TotalDistortion,e=0.02;

  if((BelongTo=(int *)malloc(pData->Row*sizeof(int)))==NULL) REPORT_ERROR(SVLIB_NoMem, "");
	if((SampleNum=(int *)malloc(this->mixtureCount*sizeof(int)))==NULL) REPORT_ERROR(SVLIB_NoMem, "");
	if((dis=(double *)malloc(this->mixtureCount*sizeof(double)))==NULL) REPORT_ERROR(SVLIB_NoMem, "");

	for(t=0;t<pData->Row;t++)
		BelongTo[t]=0;

	//************************* non-binary spliting **************************
	ClusterNo=1;
	while(1){

    // calculate centroid
		for(k=0;k<ClusterNo;k++){
			for(j=0;j<pData->Col;j++)	this->mean[k][j]=0.0;
			SampleNum[k]=0;
			for(t=0;t<pData->Row;t++){
					if(BelongTo[t]==k){
						for(j=0;j<pData->Col;j++)  this->mean[k][j] +=(double)(pData->Mat[t][j]);
						SampleNum[k]++;
					}
			}
			if(SampleNum[k]==0) printf("%d null\n",k);
			for(j=0;j<pData->Col;j++) this->mean[k][j] /=(double)SampleNum[k];
		}

		if(ClusterNo==this->mixtureCount) break;

		SplitWho=0;
		for(k=0;k<ClusterNo;k++)	if(SampleNum[k]>SampleNum[SplitWho]) SplitWho=k;
		for(j=0;j<pData->Col;j++){
			this->mean[ClusterNo][j]=this->mean[SplitWho][j]*(1.0+e);;
			this->mean[SplitWho][j] *=(1.0-e);
		}
		ClusterNo++;

    // nearest neighbor search
		TotalDistortion=0.0;
		for(t=0;t<pData->Row;t++){
			for(k=0;k<ClusterNo;k++)  dis[k]=distance(pData->Mat[t],this->mean[k],pData->Col);
			BelongTo[t]=pick_minimum(dis,ClusterNo);
			TotalDistortion +=dis[BelongTo[t]];
		}
		printf("%d %lf\n",ClusterNo,TotalDistortion);

	}
	// ********************************************************************

    // re-search and calculate centroids after cluster number has been fixed
	TotalDistortion=0.0;
	for(t=0;t<pData->Row;t++){
		for(k=0;k<ClusterNo;k++)  dis[k]=distance(pData->Mat[t],this->mean[k],pData->Col);
		BelongTo[t]=pick_minimum(dis,ClusterNo);
//      printf("%d %lf\n",BelongTo[i],dis[BelongTo[i]]);
		TotalDistortion +=dis[BelongTo[t]];
	}
	printf("%d %lf\n",ClusterNo,TotalDistortion);

	for(k=0;k<ClusterNo;k++){
		for(j=0;j<pData->Col;j++)	this->mean[k][j]=0.0;
		SampleNum[k]=0;
		for(t=0;t<pData->Row;t++){
				if(BelongTo[t]==k){
					for(j=0;j<pData->Col;j++)  this->mean[k][j] +=(double)(pData->Mat[t][j]);
					SampleNum[k]++;
				}
		}
		if(SampleNum[k]==0) printf("%d null\n",k);
		for(j=0;j<pData->Col;j++) this->mean[k][j] /=(double)SampleNum[k];
	}
	// ************************************************************************

  // calculate variance
	TotalSampleNum=0;
	for(k=0;k<this->mixtureCount;k++){
		for(j=0;j<pData->Col;j++)  this->variance[k][j]=0.0;
		SampleNum[k]=0;
		for(t=0;t<pData->Row;t++){
			if(BelongTo[t]==k){
				for(j=0;j<pData->Col;j++)  this->variance[k][j] +=((double)(pData->Mat[t][j])-this->mean[k][j])*((double)(pData->Mat[t][j])-this->mean[k][j]);
				SampleNum[k]++;
			}
		}
		if(SampleNum[k]==0) printf("%d null\n",k);
		for(j=0;j<pData->Col;j++) this->variance[k][j] /=(double)SampleNum[k];
		TotalSampleNum +=SampleNum[k];
	}

	for(i=0;i<this->mixtureCount;i++){
		for(j=0;j<pData->Col;j++){
			if(this->variance[i][j]==0.0)	this->variance[i][j]=0.000001;
			this->sd[i][j]=log(1.0/sqrt(this->variance[i][j]*6.28318));
		}
	}

	for(k=0;k<this->mixtureCount;k++) this->weight[k]=(double)SampleNum[k]/(double)TotalSampleNum;
	
	free(dis);
	free(SampleNum);
	free(BelongTo);
}


//squared distance between elements of the vectors A and B
double SC_MixtureModel_MIX2MAX::distance(float *A,double *B,int size)
{
  int i;
  double dis=0.0;
  for(i=0;i<size;i++)
	dis +=((double)(A[i])-B[i])*((double)(A[i])-B[i]);
  return dis;
}

//calculate univariate gaussian distribution
double SC_MixtureModel_MIX2MAX::LogGaussianPDF(float GX,double MeAn,double VaR,double VaRcOnSt)
{
    return (VaRcOnSt+(((double)(GX)-MeAn)*((double)(GX)-MeAn)/(-2.0*VaR)));
}

//calculate univariate cumulative gaussian distribution
double SC_MixtureModel_MIX2MAX::cumulative_distribution(float X,double MU,double SIGMA)
{
		int		ChangeSign;
		double	x,Rx;


		ChangeSign=0;
		x=((double)X-MU)/sqrt(SIGMA);
		if(x<0.0){
			x *=-1.0;
			ChangeSign=1;
		}
	
		if( (x>=0.0)&&(x<0.05) )	Rx=0.500;
		else if( (x>=0.05)&&(x<0.10) )	Rx=0.520;
		else if( (x>=0.10)&&(x<0.15) )	Rx=0.540;
		else if( (x>=0.15)&&(x<0.20) )	Rx=0.560;
		else if( (x>=0.20)&&(x<0.25) )	Rx=0.579;
		else if( (x>=0.25)&&(x<0.30) )	Rx=0.599;
		else if( (x>=0.30)&&(x<0.35) )	Rx=0.618;
		else if( (x>=0.35)&&(x<0.40) )	Rx=0.637;
		else if( (x>=0.40)&&(x<0.45) )	Rx=0.655;
		else if( (x>=0.45)&&(x<0.50) )	Rx=0.674;
		else if( (x>=0.50)&&(x<0.55) )	Rx=0.691;
		else if( (x>=0.55)&&(x<0.60) )	Rx=0.709;
		else if( (x>=0.60)&&(x<0.65) )	Rx=0.726;
		else if( (x>=0.65)&&(x<0.70) )	Rx=0.742;
		else if( (x>=0.70)&&(x<0.75) )	Rx=0.758;
		else if( (x>=0.75)&&(x<0.80) )	Rx=0.773;
		else if( (x>=0.80)&&(x<0.85) )	Rx=0.788;
		else if( (x>=0.85)&&(x<0.90) )	Rx=0.802;
		else if( (x>=0.90)&&(x<0.95) )	Rx=0.816;
		else if( (x>=0.95)&&(x<1.00) )	Rx=0.829;
		else if( (x>=1.00)&&(x<1.05) )	Rx=0.841;
		else if( (x>=1.05)&&(x<1.10) )	Rx=0.853;
		else if( (x>=1.10)&&(x<1.15) )	Rx=0.864;
		else if( (x>=1.15)&&(x<1.20) )	Rx=0.875;
		else if( (x>=1.20)&&(x<1.25) )	Rx=0.885;
		else if( (x>=1.25)&&(x<1.282) )	Rx=0.894;
		else if( (x>=1.282)&&(x<1.30) )	Rx=0.900;
		else if( (x>=1.30)&&(x<1.35) )	Rx=0.903;
		else if( (x>=1.35)&&(x<1.40) )	Rx=0.911;
		else if( (x>=1.40)&&(x<1.45) )	Rx=0.919;
		else if( (x>=1.45)&&(x<1.50) )	Rx=0.926;
		else if( (x>=1.50)&&(x<1.55) )	Rx=0.933;
		else if( (x>=1.55)&&(x<1.60) )	Rx=0.939;
		else if( (x>=1.60)&&(x<1.645) )	Rx=0.945;
		else if( (x>=1.645)&&(x<1.65) )	Rx=0.950;
		else if( (x>=1.65)&&(x<1.70) )	Rx=0.951;
		else if( (x>=1.70)&&(x<1.75) )	Rx=0.955;
		else if( (x>=1.75)&&(x<1.80) )	Rx=0.960;
		else if( (x>=1.80)&&(x<1.85) )	Rx=0.964;
		else if( (x>=1.85)&&(x<1.90) )	Rx=0.968;
		else if( (x>=1.90)&&(x<1.95) )	Rx=0.971;
		else if( (x>=1.95)&&(x<2.00) )	Rx=0.974;
		else if( (x>=2.00)&&(x<2.05) )	Rx=0.977;
		else if( (x>=2.05)&&(x<2.10) )	Rx=0.980;
		else if( (x>=2.10)&&(x<2.15) )	Rx=0.982;
		else if( (x>=2.15)&&(x<2.20) )	Rx=0.984;
		else if( (x>=2.20)&&(x<2.25) )	Rx=0.986;
		else if( (x>=2.25)&&(x<2.30) )	Rx=0.988;
		else if( (x>=2.30)&&(x<2.35) )	Rx=0.989;
		else if( (x>=2.35)&&(x<2.40) )	Rx=0.991;
		else if( (x>=2.40)&&(x<2.45) )	Rx=0.992;
		else if( (x>=2.45)&&(x<2.50) )	Rx=0.993;
		else if( (x>=2.50)&&(x<2.55) )	Rx=0.994;
		else if( (x>=2.55)&&(x<2.60) )	Rx=0.995;
		else if( (x>=2.60)&&(x<2.65) )	Rx=0.995;
		else if( (x>=2.65)&&(x<2.70) )	Rx=0.995;
		else if( (x>=2.70)&&(x<2.75) )	Rx=0.997;
		else if( (x>=2.75)&&(x<2.80) )	Rx=0.997;
		else if( (x>=2.80)&&(x<2.85) )	Rx=0.997;
		else if( (x>=2.85)&&(x<2.90) )	Rx=0.998;
		else if( (x>=2.90)&&(x<2.95) )	Rx=0.998;
		else if( (x>=2.95)&&(x<3.00) )	Rx=0.998;
		else if( (x>=3.00)&&(x<4.00) )	Rx=0.999;
		else	Rx=1.0;

		if(ChangeSign==1)	Rx=1.0-Rx;

		return Rx;
}

//return the index if the smallest element of vector q of size numb
int SC_MixtureModel_MIX2MAX::pick_minimum(double *q,int numb)
{
	int i,j;
	double temp;

	temp=q[0];
	j=0;
	for(i=1;i<numb;i++){
		if(q[i]<=temp){
			j=i;
			temp=q[i];
		}
	}

	return j;
}

double SC_MixtureModel_MIX2MAX::GaussianPDF(float *GX,double *MeAn,double *VaR,double *VaRcOnSt,int Dim)
{
  double p=0.0;
  int i;

  for(i=0;i<Dim;i++)
    p += VaRcOnSt[i]+(((double)(GX[i])-MeAn[i])*((double)(GX[i])-MeAn[i])/(-2.0*VaR[i])) ;

  return exp(p);
}

double SC_MixtureModel_MIX2MAX::cumulative_distribution(float *X,double *MU,double *SIGMA,int dim)
{
	int		i,ChangeSign;
	double	Prod,x,Rx;


	Prod=0.0;
	
	for(i=0;i<dim;i++){
		ChangeSign=0;
		x=((double)X[i]-MU[i])/sqrt(SIGMA[i]);
		if(x<0.0){
			x *=-1.0;
			ChangeSign=1;
		}
	
		if( (x>=0.0)&&(x<0.05) )	Rx=0.500;
		else if( (x>=0.05)&&(x<0.10) )	Rx=0.520;
		else if( (x>=0.10)&&(x<0.15) )	Rx=0.540;
		else if( (x>=0.15)&&(x<0.20) )	Rx=0.560;
		else if( (x>=0.20)&&(x<0.25) )	Rx=0.579;
		else if( (x>=0.25)&&(x<0.30) )	Rx=0.599;
		else if( (x>=0.30)&&(x<0.35) )	Rx=0.618;
		else if( (x>=0.35)&&(x<0.40) )	Rx=0.637;
		else if( (x>=0.40)&&(x<0.45) )	Rx=0.655;
		else if( (x>=0.45)&&(x<0.50) )	Rx=0.674;
		else if( (x>=0.50)&&(x<0.55) )	Rx=0.691;
		else if( (x>=0.55)&&(x<0.60) )	Rx=0.709;
		else if( (x>=0.60)&&(x<0.65) )	Rx=0.726;
		else if( (x>=0.65)&&(x<0.70) )	Rx=0.742;
		else if( (x>=0.70)&&(x<0.75) )	Rx=0.758;
		else if( (x>=0.75)&&(x<0.80) )	Rx=0.773;
		else if( (x>=0.80)&&(x<0.85) )	Rx=0.788;
		else if( (x>=0.85)&&(x<0.90) )	Rx=0.802;
		else if( (x>=0.90)&&(x<0.95) )	Rx=0.816;
		else if( (x>=0.95)&&(x<1.00) )	Rx=0.829;
		else if( (x>=1.00)&&(x<1.05) )	Rx=0.841;
		else if( (x>=1.05)&&(x<1.10) )	Rx=0.853;
		else if( (x>=1.10)&&(x<1.15) )	Rx=0.864;
		else if( (x>=1.15)&&(x<1.20) )	Rx=0.875;
		else if( (x>=1.20)&&(x<1.25) )	Rx=0.885;
		else if( (x>=1.25)&&(x<1.282) )	Rx=0.894;
		else if( (x>=1.282)&&(x<1.30) )	Rx=0.900;
		else if( (x>=1.30)&&(x<1.35) )	Rx=0.903;
		else if( (x>=1.35)&&(x<1.40) )	Rx=0.911;
		else if( (x>=1.40)&&(x<1.45) )	Rx=0.919;
		else if( (x>=1.45)&&(x<1.50) )	Rx=0.926;
		else if( (x>=1.50)&&(x<1.55) )	Rx=0.933;
		else if( (x>=1.55)&&(x<1.60) )	Rx=0.939;
		else if( (x>=1.60)&&(x<1.645) )	Rx=0.945;
		else if( (x>=1.645)&&(x<1.65) )	Rx=0.950;
		else if( (x>=1.65)&&(x<1.70) )	Rx=0.951;
		else if( (x>=1.70)&&(x<1.75) )	Rx=0.955;
		else if( (x>=1.75)&&(x<1.80) )	Rx=0.960;
		else if( (x>=1.80)&&(x<1.85) )	Rx=0.964;
		else if( (x>=1.85)&&(x<1.90) )	Rx=0.968;
		else if( (x>=1.90)&&(x<1.95) )	Rx=0.971;
		else if( (x>=1.95)&&(x<2.00) )	Rx=0.974;
		else if( (x>=2.00)&&(x<2.05) )	Rx=0.977;
		else if( (x>=2.05)&&(x<2.10) )	Rx=0.980;
		else if( (x>=2.10)&&(x<2.15) )	Rx=0.982;
		else if( (x>=2.15)&&(x<2.20) )	Rx=0.984;
		else if( (x>=2.20)&&(x<2.25) )	Rx=0.986;
		else if( (x>=2.25)&&(x<2.30) )	Rx=0.988;
		else if( (x>=2.30)&&(x<2.35) )	Rx=0.989;
		else if( (x>=2.35)&&(x<2.40) )	Rx=0.991;
		else if( (x>=2.40)&&(x<2.45) )	Rx=0.992;
		else if( (x>=2.45)&&(x<2.50) )	Rx=0.993;
		else if( (x>=2.50)&&(x<2.55) )	Rx=0.994;
		else if( (x>=2.55)&&(x<2.60) )	Rx=0.995;
		else if( (x>=2.60)&&(x<2.65) )	Rx=0.995;
		else if( (x>=2.65)&&(x<2.70) )	Rx=0.995;
		else if( (x>=2.70)&&(x<2.75) )	Rx=0.997;
		else if( (x>=2.75)&&(x<2.80) )	Rx=0.997;
		else if( (x>=2.80)&&(x<2.85) )	Rx=0.997;
		else if( (x>=2.85)&&(x<2.90) )	Rx=0.998;
		else if( (x>=2.90)&&(x<2.95) )	Rx=0.998;
		else if( (x>=2.95)&&(x<3.00) )	Rx=0.998;
		else if( (x>=3.00)&&(x<4.00) )	Rx=0.999;
		else	Rx=1.0;

		if(ChangeSign==1)	Rx=1.0-Rx;
		Prod +=Rx;
	}
	Prod /=(double)dim;

	return Prod;
}
*/
