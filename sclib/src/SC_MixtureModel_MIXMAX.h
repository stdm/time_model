/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a GMM-IB as described in		  			*/
/*				'Integrated Models of Signal and Background with Application to */
/*				 Speaker Identification in Noise', R.C.Rose, E.M.Hofstetter,		*/
/*			   D.A.Reynolds, 1994 (IEEE)																			*/
/*				 &																															*/	
/*				'Speech Recognition Using Noise-Adaptive Prototypes', A.Nadas,	*/
/*				 D.Nahamoo, M.A.Picheny, 1989 (IEEE)														*/
/*																																				*/
/*		Some Issues:																												*/
/*			- The Model Uses only diaginal Covariance-Matrices, because this	*/
/*				has proven to yield best results in the above paper. To draw		*/
/*				maximum reward from this simplification, we use only a vector		*/
/*				to store the variances. This simplifies some calculations as		*/
/*				well as the memory-usage																				*/
/*			- Variance-limiting is applied.																		*/
/*			- Initialization of the model prior to estimating parameters			*/
/*			  isn't well researched; here, the means are build randomly, and  */
/*				some iterations of k-means algorithm privides somehow			    	*/
/*				meaningfull distribution of the data along the mixtures 				*/
/*				according to the paper.																				  */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 06.04.2004																								*/
/**************************************************************************/

#ifndef __SC_MixtureModel_MIXMAX_H__
#define __SC_MixtureModel_MIXMAX_H__

#include "SC_MixtureModel.h"
#include "SC_MixtureModel_GMM.h"

class SC_MixtureModel_MIXMAX : public SC_MixtureModel {

  private :

  protected :
	  
	  //====================================================================================================================
	  // For each mixture, the masLevel tells to which percentage (0..1) this mixture is masked by noise
	  //====================================================================================================================
    double*             maskLevel;

	  //====================================================================================================================
	  // The type of signal-noise-interaction: additive or max (as an approximation to log-additive) as defined by the above 
    // #defines
	  //====================================================================================================================
    unsigned char       noiseCorruptionType;

	  //====================================================================================================================
	  // Model's parameters
	  //====================================================================================================================
		SC_MixtureModel*   				pBackground;						//the actual pre-estimated background-model
		SC_MixtureModel*   				pOriginalBackground;		//a copy of the background-model with which this GMM-IB was originally trained

	  //====================================================================================================================
	  // The traning/testing-routines for a max signal-noise-interaction function 
	  //====================================================================================================================
    int                 TrainWithMaxNoise(SV_Data *pData, unsigned long int segmentsToMerge);
    int                 TrainWithMaxNoise_LogArithmetic(SV_Data *pData, unsigned long int segmentsToMerge);
    SV_Data*            TestWithMaxNoise(SV_Data *pData, unsigned long int segmentsToMerge);
    SV_Data*            TestWithMaxNoise_LogArithmetic(SV_Data *pData, unsigned long int segmentsToMerge);

	  //====================================================================================================================
	  // The traning/testing-routines for an additive signal-noise-interaction function 
	  //====================================================================================================================
    //int                 TrainWithAdditiveNoise(SV_Data *pData, unsigned long int segmentsToMerge);
    //SV_Data*            TestWithAdditiveNoise(SV_Data *pData, unsigned long int segmentsToMerge);

public :

	  //====================================================================================================================
	  // constructor
	  //====================================================================================================================
		SC_MixtureModel_MIXMAX(SC_TweakableParameters* pTweak, SC_MixtureModel* pBackground, unsigned short int mixtureCount, unsigned short int dim);
    SC_MixtureModel_MIXMAX(const SC_MixtureModel_MIXMAX& pParent);
	  
		//====================================================================================================================
		// destructor
		// attention: doesn't delete the integrated background-gmm!
		//====================================================================================================================
		virtual ~SC_MixtureModel_MIXMAX();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_MixtureModel_MIXMAX& operator=(const SC_MixtureModel_MIXMAX& pParent);

	  //====================================================================================================================
	  // getters and setters
	  //====================================================================================================================
		SC_Model*				    getBackground(void) {return this->pBackground;}
		void								setBackground(SC_Model* pBackground) {this->pBackground = (SC_MixtureModel*)(pBackground); return;}
		SC_Model*				    getOriginalBackground(void) {return this->pOriginalBackground;}
		void								setOriginalBackground(SC_MixtureModel* pBackground) {MFree_0D(this->pOriginalBackground); this->pOriginalBackground = pBackground; return;}
    unsigned char       getNoiseCurruptionType(void) {return this->noiseCorruptionType;}
    void                setNoiseCorruptionType(unsigned char newType) {this->noiseCorruptionType = newType;}

	  //====================================================================================================================
	  // get & set other parameters not common to other model-types
	  //====================================================================================================================
		inline double*			getMaskLevel(void) {return this->maskLevel;}
	  inline double				getMaskLevel(unsigned short int mixture) {return this->maskLevel[mixture];}
    inline void         setMaskLevel(double* newMaskLevel) {MFree_1D(this->maskLevel); this->maskLevel = newMaskLevel;};
    inline void         setMaskLevel(short mixture, double value) {this->maskLevel[mixture] = value;}
    virtual double getWeightLimit(void) {return this->pTweak->mixtureModelMixMax.weightLimit;}

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
		virtual SC_Model* combineModels(SC_Model* pSecond);
    virtual SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);
    virtual SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepDummyMixture, bool keepFirstsNext);
		SC_MixtureModel_MIXMAX* combineModels(SC_MixtureModel* pSecond);
    SC_MixtureModel_MIXMAX* combineModels(SC_MixtureModel_MIXMAX* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext = false);
    SC_MixtureModel_MIXMAX* combineModels(SC_MixtureModel_MIXMAX* pFirst, SC_MixtureModel* pSecond, bool keepDummyMixture, bool keepFirstsNext);

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);
		SC_MixtureModel_MIXMAX* combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_MixtureModel* pBackgroundModels);

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
