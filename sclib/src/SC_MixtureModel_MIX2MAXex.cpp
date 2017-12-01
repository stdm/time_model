/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent the enhanced version of the        */
/*        abnormal GMM-IB as accidentially found by Wei Ho 'MIX2MAX' Tsai */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.06.2006																								*/
/**************************************************************************/

#include <math.h>
#include <limits>
#include <assert.h> 
#include <float.h>
#include "SC_MixtureModel_MIX2MAXex.h"
#include "SC_Aux.h"
#include "SC_Gauss.h"
#include <SV_Error.h>

#define MAX_ITER_NUM 100 
#define MAX_B_MIXTURE_NUM 64

//====================================================================================================================
// constructor
//====================================================================================================================
SC_MixtureModel_MIX2MAXex::SC_MixtureModel_MIX2MAXex(SC_TweakableParameters* pTweak, SC_MixtureModel* pBackground, unsigned short int mixtureCount, unsigned short int dim) : SC_MixtureModel(pTweak) {
	this->pBackground = pBackground;
  this->pOriginalBackground	= NULL;
	this->mixtureCount = mixtureCount; 
	this->dim = dim;
  this->maxEMsteps = this->pTweak->mixtureModelMix2Max.maxEMiterations;

	if (this->pBackground != NULL) {
    this->pOriginalBackground	= new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)pBackground));
	  assert((this->dim == this->pBackground->getDim()) || (this->pBackground->getDim() == 0));
	}

	this->Hdr.ModelType = sclib::mtMIX2MAX_ex;

	//normally, the mixtureCount must be >0, but there is one exception from this rule:
	//maybe we need only a "link" in the linked list of models, but have no data to model.
	//then, a model without any mixtures represents such a link without modelling anything.
	if (this->mixtureCount > 0) {
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex: sd");
		MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex: mean");
		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex: weight");
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex: maskLevel");

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
SC_MixtureModel_MIX2MAXex::SC_MixtureModel_MIX2MAXex(const SC_MixtureModel_MIX2MAXex& pParent) : SC_MixtureModel(pParent) {
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

//====================================================================================================================
// destructor
// attention: doesn't delete the background-gmm, nore the linked gmm-ib's!
// but: does delete the original background gmm, because it's a copy, not a reference/pointer!
//====================================================================================================================
SC_MixtureModel_MIX2MAXex::~SC_MixtureModel_MIX2MAXex() {
  this->pBackground = NULL;

  MFree_0D(this->pOriginalBackground);
  MFree_1D(this->maskLevel);
}

//====================================================================================================================
// overloaded assignment-operator
//====================================================================================================================
SC_MixtureModel_MIX2MAXex& SC_MixtureModel_MIX2MAXex::operator=(const SC_MixtureModel_MIX2MAXex& pParent) {
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
int SC_MixtureModel_MIX2MAXex::TrainModel(SV_Data* TrainData) {
	return TrainModel(TrainData, 1);
}

//====================================================================================================================
// Trains the model, using segmentsToMerge segments from the linked list of features; evaluates the 
// noiseCorruptionType and chooses the suitable algorithm
//====================================================================================================================
int SC_MixtureModel_MIX2MAXex::TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  int res = reestimation(pData, segmentsToMerge);

  return res;
}

//train the model using the data
//here, the algorithm of the SC_MixtureModel_GMM get's used now...
int SC_MixtureModel_MIX2MAXex::reestimation(SV_Data *pData, unsigned long int segmentsToMerge)
{
	int		i,j,k,t,m,EffectiveLength;
	double	**BSJointProb;
	double	*SCD,*SPDF,BCD,BPDF,BSCondProb,BSJointProb_sum,MaxLogP;
	double	OverallScore[MAX_ITER_NUM];
	double	*denominator;
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	SC_MixtureModel *pBackground	= this->pBackground; //don't touch the original pointer during training!
	//double *p_i, *Sum_p_i, **Sum_p_i_x, **Sum_p_i_xx, *pb, Sum_pb; //GMM

	if((BSJointProb=(double **)malloc(this->mixtureCount*sizeof(double *)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
	for(i=0;i<this->mixtureCount;i++){
		if((BSJointProb[i]=(double *)malloc(pBackground->getMixtureCount()*sizeof(double)))==NULL)  REPORT_ERROR(SVLIB_NoMem, "");
	} 
	denominator=(double *)malloc(this->mixtureCount*sizeof(double));
	SCD=(double *)malloc(pCompleteData->Col*sizeof(double));
	SPDF=(double *)malloc(pCompleteData->Col*sizeof(double));
	/*MArray_1D(p_i, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.TrainModel: p_i"); //GMM
	MArray_1D(pb, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.TrainModel: pb"); //GMM
	MArray_1D(Sum_p_i, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.TrainModel: Sum_p_i"); //GMM
	MArray_2D(Sum_p_i_x, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.TrainModel: Sum_p_i_x"); //GMM
	MArray_2D(Sum_p_i_xx, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.TrainModel: Sum_p_i_xx"); //GMM*/

  //by thilo
  double scaledX, **logSd, **logBgSd; //, E_xlz, E_xlz2;
  unsigned long int idx;
  MArray_2D(logSd, this->mixtureCount, this->dim, double, "logSd");
  MArray_2D(logBgSd, pBackground->getMixtureCount(), pBackground->getDim(), double, "logBgSd");
	//double mBPDF[MAX_B_MIXTURE_NUM],mBCD[MAX_B_MIXTURE_NUM],mSPDF,mSCD; //mdim

  //by thilo
  initParameters(pCompleteData, this->pTweak->mixtureModelMix2Max.varianceLimit, this->pTweak->mixtureModelMix2Max.kMeansIterations); //by thilo
  for(m=0;m<pCompleteData->Col;m++){
    for(j=0;j<pBackground->getMixtureCount();j++){
		  logBgSd[j][m] = sclib::sLog(pBackground->getSd(j, m));
    }
    for(i=0;i<this->mixtureCount;i++){
		  logSd[i][m] = sclib::sLog(this->sd[i][m]);
    }
  }
	if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIX2MAXex.txt", this, this->pTweak);}

	// maximum OverallScore search
	for(k=0;k<(long)this->maxEMsteps;k++){
		OverallScore[k]=0.0;			

		for(i=0;i<this->mixtureCount;i++){
			denominator[i]=0.0;
			/*Sum_p_i[i] = 0.0; //GMM
      for(m=0;m<pCompleteData->Col;m++) { //GMM
				Sum_p_i_x[i][m] = 0.0; //GMM
				Sum_p_i_xx[i][m] = 0.0; //GMM
      } //GMM*/
      this->maskLevel[i] = 0.0; //by thilo //TODO: compute masklevel
		}
		EffectiveLength=0;

		for(t=0;t<pCompleteData->Row;t++){

      BSJointProb_sum=0.0; //by thilo
			MaxLogP=-10000000000000.0;
			for(i=0;i<this->mixtureCount;i++){
        //pb[i] = log(this->weight[i]); //GMM
				for(m=0;m<pCompleteData->Col;m++){
					scaledX = sclib::zTransform(pCompleteData->Mat[t][m], this->mean[i][m], this->sd[i][m]);
          idx = this->gaussSolver.getIdx(scaledX);
          this->gaussSolver.tabledLogGaussianAndErf(idx, logSd[i][m], SPDF[m], SCD[m]);
          //pb[i] += SPDF[m]; //GMM
				}
        //Sum_pb = (i == 0) ? pb[i] : sclib::sLogAdd(Sum_pb, pb[i]); //GMM

				//mSPDF = this->gaussSolver.multivariateGaussian(pCompleteData->Mat[t], this->mean[i], this->variance[i], this->sd[i], this->dim); //mdim
				//mSCD = this->gaussSolver.multivariateAverageErf(pCompleteData->Mat[t], this->mean[i], this->sd[i], this->dim); //mdim
				for(j=0;j<pBackground->getMixtureCount();j++){
					/*if (i == 0) { //mDim
						mBPDF[j] = this->gaussSolver.multivariateGaussian(pCompleteData->Mat[t], this->pBackground->getMean(j), this->pBackground->getVariance(j), this->pBackground->getSd(j), this->dim); //mdim
						//mBCD[j] = this->gaussSolver.multivariateAverageErf(pCompleteData->Mat[t], this->pBackground->getMean(j), this->pBackground->getSd(j), this->dim); //mdim
					} //mDim*/
  				BSJointProb[i][j]=0.0;
					for(m=0;m<pCompleteData->Col;m++){	
						scaledX = sclib::zTransform(pCompleteData->Mat[t][m], pBackground->getMean(j, m), pBackground->getSd(j, m));
            idx = this->gaussSolver.getIdx(scaledX);
            this->gaussSolver.tabledLogGaussianAndErf(idx, logBgSd[j][m], BPDF, BCD);
						BSCondProb=BCD*exp(SPDF[m])+SCD[m]*exp(BPDF); //original version
						//BSCondProb=exp(SPDF[m]); //removed CDs & bg-PDF for test-purposes => same performance as GMM
						//BSCondProb=exp(SPDF[m]) + exp(BPDF); //removed CDs for test-purposes
						//BSCondProb=exp(SPDF[m]) + exp(BPDF) - exp(SPDF[m])*exp(BPDF); //removed CDs for test-purposes, "or"-idea from eyke huellermeier
						
						BSJointProb[i][j] +=log(BSCondProb);
					}
					//BSJointProb[i][j] = sclib::sLog(mSPDF + mBPDF[j] - mSPDF*mBPDF[j]); //mdim //sclib::sLog(mSPDF*mBCD[j]+mBPDF[j]*mSCD);
					if(BSJointProb[i][j]>MaxLogP) MaxLogP=BSJointProb[i][j];
				}
			}
			
      BSJointProb_sum=0.0;
			for(i=0;i<this->mixtureCount;i++){
        /*p_i[i] = sclib::sExp(pb[i] - Sum_pb); //GMM
        Sum_p_i[i] += p_i[i]; //GMM
        for (m = 0; m < pCompleteData->Col; m++) { //GMM
					Sum_p_i_x[i][m] += p_i[i] * pCompleteData->Mat[t][m]; //GMM
          Sum_p_i_xx[i][m] += p_i[i] * pCompleteData->Mat[t][m] * pCompleteData->Mat[t][m]; //GMM
				} //GMM*/
				for(j=0;j<pBackground->getMixtureCount();j++){
          if((BSJointProb[i][j]-MaxLogP)>-744.0)	BSJointProb[i][j]=this->weight[i]*pBackground->getWeight(j)*exp(BSJointProb[i][j]-MaxLogP);
					else	BSJointProb[i][j]=0.0;
					BSJointProb_sum +=BSJointProb[i][j];
				}
			}
      if(BSJointProb_sum!=0.0){
				OverallScore[k] +=log(BSJointProb_sum)+MaxLogP;
				for(i=0;i<this->mixtureCount;i++){
					for(j=0;j<pBackground->getMixtureCount();j++){
            BSJointProb[i][j] /=BSJointProb_sum; //by thilo (was extra loop above)
						denominator[i] +=BSJointProb[i][j];
					}
				}
				EffectiveLength++;
			}
		}

		for(i=0;i<this->mixtureCount;i++){
			for(m=0;m<pCompleteData->Col;m++){
				/*this->mean[i][m] = Sum_p_i_x[i][m] / Sum_p_i[i]; //GMM
				this->variance[i][m] = (Sum_p_i_xx[i][m] / Sum_p_i[i]) - (this->mean[i][m] * this->mean[i][m]); //GMM
				if (this->variance[i][m] < this->pTweak->mixtureModelGmm.varianceLimit) {this->variance[i][m] = this->pTweak->mixtureModelMix2MaxEx.varianceLimit;} //GMM
        this->sd[i][m] = sqrt(this->variance[i][m]); //GMM
				logSd[i][m] = sclib::sLog(this->sd[i][m]); //GMM*/
			}
			this->weight[i]=denominator[i]/(double)(EffectiveLength);
      this->maskLevel[i] /= (double)(pCompleteData->Row*pCompleteData->Col*pBackground->getMixtureCount());
		}

		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::classOut("MIX2MAXex.txt", this, this->pTweak);}
		if (this->pTweak->debug.debugMode & sclib::dbModelCreation) {sclib::scalarOut("MIX2MAXex_likelihood.txt", OverallScore[k], this->pTweak);}

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

	for(i=0;i<this->mixtureCount;i++){
		free(BSJointProb[i]);
	} 
	free(BSJointProb);
	free(denominator);
	free(SCD);
	free(SPDF);

  //by thilo
  MFree_2D(logSd);
  MFree_2D(logBgSd);
  /*MFree_1D(p_i); //GMM
	MFree_1D(pb); //GMM
	MFree_1D(Sum_p_i); //GMM
	MFree_2D(Sum_p_i_x); //GMM
	MFree_2D(Sum_p_i_xx); //GMM*/

  return k;
}

//====================================================================================================================
// Test model while considering only one element of the linked list
// return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_MIX2MAXex::TestModel(SV_Data *TestData) {
	SV_Data	*pScore = TestModel(TestData, 1);

	return pScore;					
}

//====================================================================================================================
// Test GMM-IB with a specified nr. of segments out of the linked list, return likelihood of the data given the model
//====================================================================================================================
SV_Data* SC_MixtureModel_MIX2MAXex::TestModel(SV_Data *pData, unsigned long int segmentsToMerge) {
  SV_Data *pScore = cal_likelihood(pData, segmentsToMerge);

  return pScore;
}

//test the model using the data
//should be equivalent to wesleys original algorithm...
//input is a model as saved by SC_Lib (i.e.: this->sd is really the standard deviation, not MIX2MAXs varConst...)
SV_Data* SC_MixtureModel_MIX2MAXex::cal_likelihood(SV_Data* pData, unsigned long int segmentsToMerge) 
{
	int t,i,j;
	double likelihood,BPDF[MAX_B_MIXTURE_NUM],BCD[MAX_B_MIXTURE_NUM],SPDF,SCD,Pr;
	double Blikelihood,Slikelihood;

	/*
	int d;
	double b[MAX_B_MIXTURE_NUM][MAX_B_MIXTURE_NUM], s[MAX_B_MIXTURE_NUM], p;//*/

  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	SV_Data *pScore	= new SV_Data(1, 1);
  int	PatternNum = 20;
	int	FeatDim = this->dim;
	int FeatLen = pCompleteData->Row;

	likelihood=0.0;

  for(t=0;t<FeatLen;t++){
		
		Blikelihood=0.0;
		for(j=0;j<this->pBackground->getMixtureCount();j++){
      BPDF[j] = this->gaussSolver.multivariateGaussian(pCompleteData->Mat[t], this->pBackground->getMean(j), this->pBackground->getVariance(j), this->pBackground->getSd(j), this->dim);
			BCD[j] = this->gaussSolver.multivariateAverageErf(pCompleteData->Mat[t], this->pBackground->getMean(j), this->pBackground->getSd(j), this->dim);
			Blikelihood +=BPDF[j];
		
			/*
			for (d = 0; d < pCompleteData->Col; d++) {
				b[j][d] = this->gaussSolver.gaussian(pCompleteData->Mat[t][d], this->pBackground->getMean(j, d), this->pBackground->getVariance(j, d), this->pBackground->getSd(j, d));
			}//*/
		}

		Pr=0.0;
		Slikelihood=0.0;
		for(i=0;i<this->mixtureCount;i++){
      SPDF = this->gaussSolver.multivariateGaussian(pCompleteData->Mat[t], this->mean[i], this->variance[i], this->sd[i], this->dim);
			SCD = this->gaussSolver.multivariateAverageErf(pCompleteData->Mat[t], this->mean[i], this->sd[i], this->dim);
      Slikelihood +=SPDF;
		
			///*
      for(j=0;j<this->pBackground->getMixtureCount();j++) {
				Pr +=this->weight[i]*this->pBackground->getWeight(j)*(SPDF*BCD[j]+BPDF[j]*SCD); //original
				//Pr +=this->weight[i]*this->pBackground->getWeight(j)*(SPDF + BPDF[j]); //suggestion by reviewer 2 => same result as above => the difference in order of magnitude between the PDFs and CDs is so great that the formula is completely dominated by the PDF terms
				//Pr +=this->weight[i]*this->pBackground->getWeight(j)*(SPDF + BPDF[j] - SPDF*BPDF[j]); //suggestion by eyke huellermeier
				
				/*
				p = 1.0;
				for (d = 0; d < pCompleteData->Col; d++) {
					s[d] = this->gaussSolver.gaussian(pCompleteData->Mat[t][d], this->mean[i][d], this->variance[i][d], this->sd[i][d]);
					//p *= s[d] + b[j][d] - s[d]*b[j][d];
					p *= sclib::max(s[d], b[j][d]); //2 gmms per-frame per-dim compensation
				}
				Pr += this->weight[i] * this->pBackground->getWeight(j) * p;//*/
      }
		}


		likelihood +=log(Pr); //original
		//likelihood += log(sclib::max(Slikelihood, Blikelihood)); //simple scheme of reviewer 2: frame is either full noise or full voice
	}	

	pScore->Mat[0][0] = (float)(likelihood) / (float)(FeatLen); //divide by T as suggested in "Speaker Verification Using Adapted Gaussian Mixture Models"...;
  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }

	return pScore;
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of this and the pSecond and return a new model
//====================================================================================================================
SC_Model*	SC_MixtureModel_MIX2MAXex::combineModels(SC_Model* pSecond) {
  return (SC_Model*)(combineModels(this, (SC_MixtureModel*)pSecond, false));
}
SC_MixtureModel_MIX2MAXex* SC_MixtureModel_MIX2MAXex::combineModels(SC_MixtureModel* pSecond) {
  return combineModels(this, pSecond, false);
}

//====================================================================================================================
// Combine 2 Models by adding the mixtures of pFirst and the pSecond and return a new model
//====================================================================================================================
SC_Model*  SC_MixtureModel_MIX2MAXex::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext) {
   return (SC_Model*)(combineModels((SC_MixtureModel_MIX2MAXex*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_MIX2MAXex* SC_MixtureModel_MIX2MAXex::combineModels(SC_MixtureModel_MIX2MAXex* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext) {
  unsigned short int i, d, newMixtureCount = pFirst->getMixtureCount() + pSecond->getMixtureCount();
	unsigned long int newTrainingDataCount = pFirst->getTrainingDataCount() + pSecond->getTrainingDataCount();
	double *newWeight, *newMaskLevel, **newMean, **newVariance, **newSd;
	SC_MixtureModel_MIX2MAXex *pNewModel;
  SC_MixtureModel *pNewOriginalBackground;
  bool tooSmall = false;
  
	if (pFirst->getDim() != pSecond->getDim()) {REPORT_ERROR(SVLIB_BadData, "Can't combine Models trained on different feature-dim!");}

	MArray_1D(newWeight, newMixtureCount, double, "SC_MixtureModel_MIX2MAXex.combineModels: newWeight");
  MArray_1D(newMaskLevel, newMixtureCount, double, "SC_MixtureModel_MIX2MAXex.combineModels: newMaskLevel");
	MArray_2D(newMean, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIX2MAXex.combineModels: newMean");
	MArray_2D(newVariance, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIX2MAXex.combineModels: newVariance");
  MArray_2D(newSd, newMixtureCount, pFirst->getDim(), double, "SC_MixtureModel_MIX2MAXex.combineModels: newSd");

  pNewModel = new SC_MixtureModel_MIX2MAXex(pFirst->getTweak(), NULL, newMixtureCount, pFirst->getDim());

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
        newMaskLevel[i] = ((SC_MixtureModel_MIX2MAXex*)pSecond)->getMaskLevel(i - pFirst->getMixtureCount());
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
    pNewOriginalBackground = (SC_MixtureModel*)(pFirst->getOriginalBackground()->combineModels(((SC_MixtureModel_MIX2MAXex*)pSecond)->getOriginalBackground()));
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
SC_Model* SC_MixtureModel_MIX2MAXex::combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepDummyMixture, bool keepFirstsNext) {
  return (SC_Model*)(combineModels((SC_MixtureModel_MIX2MAXex*)pFirst, (SC_MixtureModel*)pSecond, keepFirstsNext));
}
SC_MixtureModel_MIX2MAXex* SC_MixtureModel_MIX2MAXex::combineModels(SC_MixtureModel_MIX2MAXex* pFirst, SC_MixtureModel* pSecond, bool keepDummyMixture, bool keepFirstsNext) {
  assert("Don't use this method!!!" == 0); //TODO: why???
  
  return combineModels(pFirst, pSecond, keepFirstsNext);
}

//====================================================================================================================
// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
//====================================================================================================================
SC_Model* SC_MixtureModel_MIX2MAXex::combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels) {
  return (SC_Model*)(combineModels((SC_MixtureModel*)pSecond, pSpeechFrames, segmentsToMerge, pBackgroundModels));
}
SC_MixtureModel_MIX2MAXex* SC_MixtureModel_MIX2MAXex::combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_MixtureModel* pBackgroundModels) {
	SC_MixtureModel_MIX2MAXex* pNewModel;
	
  assert("Don't use this method!!!" == 0); //TODO: why???
  assert(this->dim == pSecond->getDim());
	
	pNewModel = new SC_MixtureModel_MIX2MAXex(this->pTweak, pBackgroundModels, sclib::min(this->mixtureCount+pSecond->getMixtureCount(), this->pTweak->modelHandler.maxSpeakerModelOrder), this->dim);
	pNewModel->TrainModel(pSpeechFrames, segmentsToMerge);

	return pNewModel;
}

//====================================================================================================================
// Sort the model-parameters according to the weights
//====================================================================================================================
void SC_MixtureModel_MIX2MAXex::sortParameters() {
	unsigned short int i, d;
	SV_Data *temp;
	double *oldWeight, *oldMaskLevel, **oldMean, **oldVariance, **oldSd;
	
	temp = new SV_Data(this->mixtureCount, 2);
	MArray_1D(oldWeight, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.sortParameters: oldWeight");
	MArray_1D(oldMaskLevel, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.sortParameters: oldMaskLevel");
	MArray_2D(oldMean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.sortParameters: oldMean");
	MArray_2D(oldVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.sortParameters: oldVariance");
  MArray_2D(oldSd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.sortParameters: oldVariance");

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
void SC_MixtureModel_MIX2MAXex::killMixture(unsigned short int mixture) {
  unsigned short int i, i_idx = 0, d,  I = this->mixtureCount, D = this->dim;
  double *newWeight, *newMaskLevel, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I-1, double, "SC_MixtureModel_MIX2MAXex.killMixture: newWeight");
  MArray_1D(newMaskLevel, I-1, double, "SC_MixtureModel_MIX2MAXex.killMixture: newMaskLevel");
  MArray_2D(newVariance, I-1, D, double, "SC_MixtureModel_MIX2MAXex.killMixture: newVariance");
  MArray_2D(newMean, I-1, D, double, "SC_MixtureModel_MIX2MAXex.killMixture: newMean");
  MArray_2D(newSd, I-1, D, double, "SC_MixtureModel_MIX2MAXex.killMixture: newSd");

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
void SC_MixtureModel_MIX2MAXex::addMixture(void) {
  unsigned short int i, d, I = this->mixtureCount, D = this->dim;
  double *newWeight, *newMaskLevel, **newVariance, **newMean, **newSd;

  MArray_1D(newWeight, I+1, double, "SC_MixtureModel_MIX2MAXex.addMixture: newWeight");
  MArray_1D(newMaskLevel, I+1, double, "SC_MixtureModel_MIX2MAXex.addMixture: newMaskLevel");
  MArray_2D(newVariance, I+1, D, double, "SC_MixtureModel_MIX2MAXex.addMixture: newVariance");
  MArray_2D(newMean, I+1, D, double, "SC_MixtureModel_MIX2MAXex.addMixture: newMean");
  MArray_2D(newSd, I+1, D, double, "SC_MixtureModel_MIX2MAXex.addMixture: newSd");

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
int SC_MixtureModel_MIX2MAXex::SaveModel(void) {
 	int res, x, bytes;
	SV_DataIO io;

  res = SaveHdr();
	if (res == SVLIB_Fail) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex model Failed!");}

  //write mixtureCount and featureDim
	bytes = io.writeScalar(&(this->DFile), this->dim);
	bytes += io.writeScalar(&(this->DFile), this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex Model Failed!");}

  //write weight-vector
	bytes += io.writeArray(&(this->DFile), this->weight, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex Model Failed!");}

  //write maskLevel-vector
	bytes += io.writeArray(&(this->DFile), this->maskLevel, this->mixtureCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex Model Failed!");}

	//write mean and variance matrices
	for (x = 0; x < this->mixtureCount; x++) {
		bytes += io.writeArray(&(this->DFile), this->mean[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->variance[x], this->dim);
		bytes += io.writeArray(&(this->DFile), this->sd[x], this->dim);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex Model Failed!");}
	}

  //write trainingDataCount
	bytes += io.writeScalar(&(this->DFile), this->trainingDataCount);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex Model Failed!");}

  //write maxEMstep
	io.writeScalar(&(this->DFile), this->maxEMsteps);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Save SC_MixtureModel_MIX2MAXex Failed!");}

  return bytes + MHLen; //MHLen is just an estimate of the written header-bytes...
}

//====================================================================================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//====================================================================================================================
SV_Model* SC_MixtureModel_MIX2MAXex::LoadModel(void) {
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

 		MArray_1D(this->weight, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex: weight");
    MArray_1D(this->maskLevel, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex: maskLevel");
    MArray_2D(this->mean, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex: mean");
		MArray_2D(this->variance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex: variance");
    MArray_2D(this->sd, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex: sd");
	};

	//read weight, maskLevel, mean and variance
	io.readArray(&(this->DFile), this->weight, this->mixtureCount, codeSizes, fileSizes);

	//TODO: remove output
	//sclib::vectorOut("weights.txt", this->weight, this->mixtureCount, false, this->pTweak);

	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAXex Failed!");}
	io.readArray(&(this->DFile), this->maskLevel, this->mixtureCount, codeSizes, fileSizes);
	if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAXex Failed!");}
	for (x = 0; x < this->mixtureCount; x++) {
		io.readArray(&(this->DFile), this->mean[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->variance[x], this->dim, codeSizes, fileSizes);
		io.readArray(&(this->DFile), this->sd[x], this->dim, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAXex Failed!");}
	}

  //read trainingDataCount
	io.readScalar(&(this->DFile), this->trainingDataCount, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAXex Failed!");}

  //Set a dummy as the original background
  pDummy = new SC_MixtureModel_GMM(this->pTweak, 0, this->dim);
  this->pOriginalBackground = pDummy;

  //read maxEMstep
	io.readScalar(&(this->DFile), this->maxEMsteps, codeSizes, fileSizes);
  if (DFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Load SC_MixtureModel_MIX2MAXex Failed!");}

	return(this);
}

//====================================================================================================================
// create a (linked, not copied) SC_Signature-view on this model (for distance computation )
//====================================================================================================================
SC_Signature* SC_MixtureModel_MIX2MAXex::toSignature(void) {
  double *weights;
  SC_Centroid **centroids = NULL;
  SC_Signature *pSignature;
  
  MArray_1D(centroids, this->mixtureCount, SC_Centroid*, "SC_MixtureModel_MIX2MAXex.toSignature: centroids");
  MArray_1D(weights, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.toSignature: weights");

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
void SC_MixtureModel_MIX2MAXex::killSignature(SC_Signature *pSignature) {
  if (pSignature != NULL) {
    MFree_0D(pSignature);
  }

  return;
}

/*
int SC_MixtureModel_MIX2MAXex::reestimation(SV_Data *pData, unsigned long int segmentsToMerge) {
	double Log_p_X = 0.0, Log_p_X_old; //the log-likelihood of the data given the actual/last model: p(X|lambda)
	double *p_i; //class-probability for acoustic classes (mixture-component) i: p(i|x_t, lambda) for specific t and all i
	double *Sum_p_i; //SUM_t=1..T p(i|x_t, lambda) for all i 
	double **Sum_p_i_x; //SUM_t=1..T p(i|x_t, lambda)*x_t for all i 
	double **Sum_p_i_xx; //SUM_t=1..T p(i|x_t, lambda)*(x_t)^2 for all i 
	double *pb; //values of gauss-densities multiplied with the class-weights for a specific feature-vector and all mixtures: p_i * b_i(x_t)
	double Sum_pb; //SUM_i=1..M p_i * b_i(x_t) for a specific t
	double Log_Prod_pb; //log(PROD_t=1..T SUM_i=1..M p_i * b_i(x_t))
	unsigned long int t; //0..T over all feature-vectors
	unsigned short int i; //0..M over all mixture-components
	unsigned short int d; //0..D over all components of the feature-vectors
	unsigned long int iterationCount = 0;
  SV_Data *pCompleteData = (segmentsToMerge == 1 || sclib::getListCount(pData) == 1) ? pData : pData->MergeData(segmentsToMerge);
	SC_MixtureModel *pBackground	= this->pBackground; //don't touch the original pointer during training!
  double **scaledVariance;

	MArray_1D(p_i, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.TrainModel: p_i");
	MArray_1D(pb, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.TrainModel: pb");
	MArray_1D(Sum_p_i, this->mixtureCount, double, "SC_MixtureModel_MIX2MAXex.TrainModel: Sum_p_i");
	MArray_2D(Sum_p_i_x, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.TrainModel: Sum_p_i_x");
	MArray_2D(Sum_p_i_xx, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.TrainModel: Sum_p_i_xx");
  MArray_2D(scaledVariance, this->mixtureCount, this->dim, double, "SC_MixtureModel_MIX2MAXex.TrainModel: scaledVariance");
	
	//The model-parameters must be initialized some way prior to the em-algorithm
	initParameters(pCompleteData, this->pTweak->mixtureModelMix2MaxEx.varianceLimit, this->pTweak->mixtureModelMix2MaxEx.kMeansIterations);
	//if (this->pTweak->debug.debugMode & SCLIB_DB_GMM) {sclib::classOut("gmm.txt", this, this->pTweak);}

  //initialize the scaled variance to speed up gaussian computation
	for (i = 0; i < this->mixtureCount; i++) {
  	for (d = 0; d < this->dim; d++) {
      scaledVariance[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
		}			
  }

	//now we can start... ladys 'n gentleman: the EM-Algorithm ;-)
	do {
		Log_p_X_old = Log_p_X;
		iterationCount++;

		//calculate all SUM/PROD terms of the formulas a-priori in some loops:
		Log_Prod_pb = 0.0;
		for (i = 0; i < this->mixtureCount; i++) { //1. initialize with zero... needed by the third loop
			Sum_p_i[i] = 0.0;
			for (d = 0; d < this->dim; d++) {
				Sum_p_i_x[i][d] = 0.0;
				Sum_p_i_xx[i][d] = 0.0;
			}
      this->maskLevel[i] = 0.0; //new compared to GMM algorithm
		}

		for (t = 0; t < (unsigned long int)pCompleteData->Row; t++) {
			for (i = 0; i < this->mixtureCount; i++) { //2. compute the sum_pb over all mixtures i needed for p_X and the third loop
        pb[i] = log(this->weight[i]);
				for (d = 0; d < this->dim; d++) { //the multivariate gaussian pdf may be evalutated this way (component-wise), because we have diagonal covariances (=> features are assumed uncorrelated)!
          pb[i] += scaledVariance[i][d] + (-0.5 * (pCompleteData->Mat[t][d] - this->mean[i][d]) * (pCompleteData->Mat[t][d] - this->mean[i][d]) / this->variance[i][d]);
        }
        Sum_pb = (i == 0) ? pb[i] : sclib::sLogAdd(Sum_pb, pb[i]);
			}//i=0..M

      Log_Prod_pb += Sum_pb;

			for (i = 0; i < this->mixtureCount; i++) { //3. compute the p_i's needed by the reestimation formulas
        p_i[i] = sclib::sExp(pb[i] - Sum_pb);
        Sum_p_i[i] += p_i[i];
        for (d = 0; d < this->dim; d++) {
					Sum_p_i_x[i][d] += p_i[i] * pCompleteData->Mat[t][d];
          Sum_p_i_xx[i][d] += p_i[i] * pCompleteData->Mat[t][d] * pCompleteData->Mat[t][d];
				}
			} //i=0..M
		} //t=0..T

		Log_p_X = Log_Prod_pb;

		//let's start with the reestimation process:
		for (i = 0; i < this->mixtureCount; i++) {
			this->weight[i] = Sum_p_i[i] / pCompleteData->Row; //p_i* = 1/T SUM_t=1..T p(i|x_t, lambda)
      //assert(this->weight[i] > 0.0);
      //assert(Sum_p_i[i] != 0.0);
      if (this->weight[i] > 0.0) {
  		  for (d = 0; d < this->dim; d++) {
				  this->mean[i][d] = Sum_p_i_x[i][d] / Sum_p_i[i]; //mu_i* = (SUM_t=1..T p(i|x_t, lambda)*x_t) / (SUM_t=1..T p(i|x_t, lambda))
				  this->variance[i][d] = (Sum_p_i_xx[i][d] / Sum_p_i[i]) - (this->mean[i][d] * this->mean[i][d]); //sigma_i* = ((SUM_t=1..T p(i|x_t, lambda)*(x_t)^2) / (SUM_t=1..T p(i|x_t, lambda))) - (mu_i*)^2
				  if (this->variance[i][d] < this->pTweak->mixtureModelGmm.varianceLimit) {this->variance[i][d] = this->pTweak->mixtureModelMix2MaxEx.varianceLimit;} //do variance-limiting
          this->sd[i][d] = sqrt(this->variance[i][d]);
          scaledVariance[i][d] = log(1.0 / (this->sd[i][d] * sclib::sqrt_2pi)); //precompute the first factor in the normal density so the fast gaussian calculation can be used
			  }			
      } else {
        killMixture(i);
        printf("\n     mixture %i killed in run %i, %i remaining", i, iterationCount, this->mixtureCount);
      }
  	}
  	
  	//if (this->pTweak->debug.debugMode & SCLIB_DB_GMM) {sclib::classOut("gmm.txt", this, this->pTweak);}
		//if (this->pTweak->debug.debugMode & SCLIB_DB_GMM) {sclib::scalarOut("gmm_likelihood.txt", Log_p_X, this->pTweak);}

    //assert(Log_p_X >= Log_p_X_old || iterationCount == 1);
    if (!(Log_p_X >= Log_p_X_old || iterationCount == 1)) {
      printf("\n     abnormal log-likelihood: %f (old: %f)", Log_p_X, Log_p_X_old);
    }
    
    if (iterationCount == this->maxEMsteps) {
      break;
    }

	} while (fabs(Log_p_X - Log_p_X_old) > fabs(this->pTweak->mixtureModelMix2MaxEx.EMthreshold));

	this->trainingDataCount = pCompleteData->Row;
  this->pOriginalBackground =  new SC_MixtureModel_GMM(*((SC_MixtureModel_GMM*)this->pBackground)); //new compared to GMM algorithm
  
  //kill too small mixtures
  this->sortParameters();
  for (i = 0; i < this->mixtureCount; i++) {
    if (this->weight[i] < this->pTweak->mixtureModelMix2MaxEx.weightLimit) { //kill too small mixtures
      this->killMixture(i);
      i--;
    }
  }

  if (pData != pCompleteData) {
    MFree_0D(pCompleteData);
  }
	MFree_1D(p_i);
	MFree_1D(pb);
	MFree_1D(Sum_p_i);
	MFree_2D(Sum_p_i_x);
	MFree_2D(Sum_p_i_xx);
  MFree_2D(scaledVariance);
	
	return iterationCount;
}
*/
