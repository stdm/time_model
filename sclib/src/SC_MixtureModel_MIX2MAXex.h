/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent the enhanced version of the        */
/*        abnormal GMM-IB as accidentially found by Wei Ho 'MIX2MAX' Tsai  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 21.06.2006																								*/
/**************************************************************************/

#ifndef __SC_MixtureModel_MIX2MAXex_H__
#define __SC_MixtureModel_MIX2MAXex_H__

#include "SC_MixtureModel.h"
#include "SC_MixtureModel_GMM.h"

class SC_MixtureModel_MIX2MAXex : public SC_MixtureModel {

  private :

  protected :
	  
	  //====================================================================================================================
	  // For each mixture, the masLevel tells to which percentage (0..1) this mixture is masked by noise
	  //====================================================================================================================
    double*             maskLevel;

	  //====================================================================================================================
	  // Model's parameters
	  //====================================================================================================================
		SC_MixtureModel*   				pBackground;						//the actual pre-estimated background-model
		SC_MixtureModel*   				pOriginalBackground;		//a copy of the background-model with which this model was originally trained

    //wesleys test-algorithms and auxiliary functions:
    SV_Data* cal_likelihood(SV_Data* pData, unsigned long int segmentsToMerge);
    //double cumulative_distribution(float *X,double *MU,double *SIGMA,int dim);
    //double GaussianPDF(float *GX,double *MeAn,double *VaR,double *VaRcOnSt,int Dim);

    //wesleys training-algorithms and auxiliary functions:
    int reestimation(SV_Data *pData, unsigned long int segmentsToMerge);
    //void K_means(SV_Data* pData);
    //int pick_minimum(double *q,int numb);
    //double distance(float *A,double *B,int size);
    //double LogGaussianPDF(float GX,double MeAn,double VaR,double VaRcOnSt);
    //double cumulative_distribution(float X,double MU,double SIGMA);

public :

	  //====================================================================================================================
	  // constructor
	  //====================================================================================================================
		SC_MixtureModel_MIX2MAXex(SC_TweakableParameters* pTweak, SC_MixtureModel* pBackground, unsigned short int mixtureCount, unsigned short int dim);
    SC_MixtureModel_MIX2MAXex(const SC_MixtureModel_MIX2MAXex& pParent);
	  
		//====================================================================================================================
		// destructor
		// attention: doesn't delete the integrated background-gmm!
		//====================================================================================================================
		virtual ~SC_MixtureModel_MIX2MAXex();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_MixtureModel_MIX2MAXex& operator=(const SC_MixtureModel_MIX2MAXex& pParent);

	  //====================================================================================================================
	  // getters and setters
	  //====================================================================================================================
		SC_Model*				    getBackground(void) {return this->pBackground;}
		void								setBackground(SC_Model* pBackground) {this->pBackground = (SC_MixtureModel*)(pBackground); return;}
		SC_Model*				    getOriginalBackground(void) {return this->pOriginalBackground;}
		void								setOriginalBackground(SC_MixtureModel* pBackground) {MFree_0D(this->pOriginalBackground); this->pOriginalBackground = pBackground; return;}

	  //====================================================================================================================
	  // get & set other parameters not common to other model-types
	  //====================================================================================================================
		inline double*			getMaskLevel(void) {return this->maskLevel;}
	  inline double				getMaskLevel(unsigned short int mixture) {return this->maskLevel[mixture];}
    inline void         setMaskLevel(double* newMaskLevel) {MFree_1D(this->maskLevel); this->maskLevel = newMaskLevel;};
    inline void         setMaskLevel(short mixture, double value) {this->maskLevel[mixture] = value;}
    virtual double getWeightLimit(void) {return this->pTweak->mixtureModelMix2Max.weightLimit;}

	  //====================================================================================================================
	  // methods to estimate and test a model
	  //====================================================================================================================
	  virtual int  				TrainModel(SV_Data *TrainData);
	  virtual int  				TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
	  virtual SV_Data*		TestModel(SV_Data *TestData);
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge);

		//====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond);
    virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);
    virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepDummyMixture, bool keepFirstsNext);
		SC_MixtureModel_MIX2MAXex* combineModels(SC_MixtureModel* pSecond);
    SC_MixtureModel_MIX2MAXex* combineModels(SC_MixtureModel_MIX2MAXex* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext = false);
    SC_MixtureModel_MIX2MAXex* combineModels(SC_MixtureModel_MIX2MAXex* pFirst, SC_MixtureModel* pSecond, bool keepDummyMixture, bool keepFirstsNext);

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);
		SC_MixtureModel_MIX2MAXex* combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_MixtureModel* pBackgroundModels);

		//====================================================================================================================
	  // Sort the model-parameters according to the weights
	  //====================================================================================================================
    virtual void				sortParameters(void);

    //====================================================================================================================
	  // Remove the indexed mixture-component from the models mean/variance/weight-vectors/matrices
	  //====================================================================================================================
    virtual void        killMixture(unsigned short int mixture);
    
    //====================================================================================================================
    // Add one new mixture component with all means=0, variances=1.0, weight=1/mixtureCount
    //====================================================================================================================
    virtual void        addMixture(void);

	  //====================================================================================================================
	  // load/save a model from/to file
	  //====================================================================================================================
	  virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
    virtual SC_Signature* toSignature(void);
    virtual void          killSignature(SC_Signature *pSignature);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) {return this->mixtureCount*this->dim*2 + this->mixtureCount*2;} //mean & variance & weight & maskLevel per dimension and mixture
};

#endif
