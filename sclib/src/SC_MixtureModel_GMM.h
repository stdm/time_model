/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a GMM as described in						  */
/*				'Robust Text-Independent Speaker Identification Using Gaussian  */
/*				'Mixture Speaker Models', D.Reynolds, R.C.Rose, 1995 (IEEE)			*/
/*																																				*/
/*		Some Issues:																												*/
/*			- The Model Uses only diagonal Covariance-Matrices, because this	*/
/*				has proven to yield best results in the above paper. To draw		*/
/*				maximum reward from this simplification, we use only a vector		*/
/*				to store the variances. This simplifies some calculations as		*/
/*				well as the memory-usage																				*/
/*			- Variance-limiting is applied.																		*/
/*			- Initialization of the model prior to estimating parameters			*/
/*			  isn't well researched; here, the means are build randomly, and  */
/*				some iterations of k-means algorithm privides somehow				    */
/*				meaningfull distribution of the data along the mixtures 				*/
/*				accordning to the paper.																				*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 06.04.2004																								*/
/**************************************************************************/

#ifndef __SC_MixtureModel_GMM_H__
#define __SC_MixtureModel_GMM_H__

#include "SC_MixtureModel.h"

class SC_MixtureModel_GMM : public SC_MixtureModel {

  private :

  protected :

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
		SC_MixtureModel_GMM(SC_TweakableParameters* pTweak, unsigned short int mixtureCount, unsigned short int dim);
    SC_MixtureModel_GMM(SC_TweakableParameters *pTweak, unsigned short int mixtureCount, unsigned short int dim, double *weight, double **mean, double **variance);
		SC_MixtureModel_GMM(const SC_MixtureModel_GMM& pParent);
	  virtual ~SC_MixtureModel_GMM();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_MixtureModel_GMM& operator=(const SC_MixtureModel_GMM& pParent);

	  //====================================================================================================================
	  // methods to estimate and test a model
	  //====================================================================================================================
	  virtual int         TrainModel(SV_Data *TrainData);
		virtual int				  TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
	  virtual SV_Data*		TestModel(SV_Data *TestData);
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge);

    //====================================================================================================================
	  // get & set (new) background model; (wo do not have a background-model in a gmm, so return NULL/do nothing)
    // this is just for compatibility with the functions in SC_Cluster, which have to work also with an GMM-IB!
	  //====================================================================================================================
    virtual SC_Model*		getBackground() {return NULL;}
    virtual void		    setBackground(SC_Model* pBackground) {return;}

    virtual double getWeightLimit(void) {return this->pTweak->mixtureModelGmm.weightLimit;}

    //====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond);
    virtual SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);
		SC_MixtureModel_GMM* combineModels(SC_MixtureModel* pSecond);
    SC_MixtureModel_GMM* combineModels(SC_MixtureModel* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext = false);

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);
    SC_MixtureModel_GMM* combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) {return this->mixtureCount*this->dim*2 + this->mixtureCount;} //mean & variance & weight per dimension and mixture

	  //====================================================================================================================
	  // load/save a model from/to file
	  //====================================================================================================================
	  virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);
};

#endif
