/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a Baggenstoss'-GMM as described in	*/
/*				'Robust Text-Independent Speaker Identification Using Gaussian  */
/*				'Mixture Speaker Models', D.Reynolds, R.C.Rose, 1995 (IEEE)			*/
/*				'Statistical Modeling Using Gaussian Mixtures and HMM with      */
/*         Matlab', P.M.Baggenstoss, 2002 (web)														*/
/*																																				*/
/*		Some Issues:																												*/
/*		  -	This model represents a GMM using P.M.Baggenstoss' algorithms		*/
/*				for training etc. in the full-covariance- and altered version   */
/*        for diagonal covariances. It is capable of deleting, splitting  */
/*        and merging mixtures during training/combining :-)							*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 29.03.2005																								*/
/**************************************************************************/

#ifndef __SC_MixtureModel_bGMM_H__
#define __SC_MixtureModel_bGMM_H__

#include "SC_MixtureModel.h"
#include "SC_BaggenstossEM.h"
#include "SC_BaggenstossEMex.h"

class SC_MixtureModel_bGMM : public SC_MixtureModel {

  private :

  protected :

    unsigned short int  maxMixtureCount; //upper bound for mixture-splitting
    double* minVariance; //minimum variance in contrast to min. standard-deviation in basic class.
		bool fullCovariance;
		bool verbose;

		SC_BaggenstossEM *pEM; //the real thing

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
		SC_MixtureModel_bGMM(SC_TweakableParameters* pTweak, unsigned short int mixtureCount, unsigned short int dim, bool fullCovariance, bool verbose);
    SC_MixtureModel_bGMM(const SC_MixtureModel_bGMM& pParent);
	  virtual ~SC_MixtureModel_bGMM();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_MixtureModel_bGMM& operator=(const SC_MixtureModel_bGMM& pParent);

	  //====================================================================================================================
	  // methods to estimate and test a model
	  //====================================================================================================================
	  virtual int         TrainModel(SV_Data *pData);
    virtual int         TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
		virtual int         TrainModel(SV_Data *pData, unsigned long int segmentsToMerge, long int samplesPerMode, bool bias, double maxCloseness, bool addModes, double splitThresh);
	  virtual SV_Data*		TestModel(SV_Data *pData);
		virtual SV_Data*		TestModel(SV_Data *pData, unsigned long int segmentsToMerge);

    //====================================================================================================================
	  // get & set (new) background model; (wo do not have a background-model in a gmm, so return NULL/do nothing)
    // this is just for compatibility with the functions in SC_Cluster, which have to work also with an GMM-IB!
	  //====================================================================================================================
    virtual SC_Model*		getBackground() {return NULL;}
    virtual void        setBackground(SC_Model* pBackground) {return;}

    virtual double getWeightLimit(void) {return this->pTweak->mixtureModelBgmm.weightLimit;}

    //====================================================================================================================
    // Set the minimum-variance from outside
    //====================================================================================================================
    void                setMinVariance(double* newMinVariance);
    void                setMinVariance(double newMinVariance);
    
    //====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond);
    virtual SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);
    SC_MixtureModel_bGMM* combineModels(SC_MixtureModel_bGMM* pSecond);
    SC_MixtureModel_bGMM* combineModels(SC_MixtureModel_bGMM* pFirst, SC_MixtureModel_bGMM* pSecond, bool keepFirstsNext = false);

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);
    SC_MixtureModel_bGMM* combineModels(SC_MixtureModel_bGMM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);

		//====================================================================================================================
		// Two new getters/setters...
		//====================================================================================================================
		double** getCovar(unsigned short int mixture) {return (this->pEM != NULL) ? this->pEM->getCovar(mixture) : NULL;} //computes and returns the covariance-matrix of the given mixture
    void setMaxMixtureCount(unsigned short int maxMixtureCount) {this->maxMixtureCount = maxMixtureCount;}
		bool isFullCovar(void) {return this->fullCovariance;}
		bool isVerbose(void) {return this->verbose;}

	  //====================================================================================================================
	  // load/save a model from/to file
	  //====================================================================================================================
	  virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void);

    //====================================================================================================================
	  // to adapt to the interface
	  //====================================================================================================================
		virtual void				sortParameters(void) {return;} //not available
		virtual void        killMixture(unsigned short int mixture) {return;} //not available
		virtual void        addMixture(void) {return;} //not available

		virtual inline double*			getWeight(void) {return this->pEM->getWeight();}
	  virtual inline double				getWeight(unsigned short int mixture) {return this->pEM->getWeight()[mixture];}
		virtual inline double**			getVariance(void) {return (this->fullCovariance == false) ? ((SC_BaggenstossEMex*)this->pEM)->getVariance() : NULL;}
		virtual inline double*			getVariance(unsigned short int mixture) {return (this->fullCovariance == false) ? ((SC_BaggenstossEMex*)this->pEM)->getVariance()[mixture] : NULL;}
		virtual inline double				getVariance(unsigned short int mixture, unsigned short int dim) {return (this->fullCovariance == false) ? ((SC_BaggenstossEMex*)this->pEM)->getVariance()[mixture][dim] : -1.0;}
		virtual inline double**			getSd(void) {return NULL;}
    virtual inline double*			getSd(unsigned short int mixture) {return NULL;}
    virtual inline double				getSd(unsigned short int mixture, unsigned short int dim) {return -1.0;}
		virtual inline double**			getMean(void) {return this->pEM->getMean();}
    virtual inline double*			getMean(unsigned short int mixture) {return this->pEM->getMean()[mixture];}
	  virtual inline double				getMean(unsigned short int mixture, unsigned short int dim) {return this->pEM->getMean()[mixture][dim];}
    virtual inline void         setVariance(double** newVariance) {if (this->fullCovariance == false) {((SC_BaggenstossEMex*)this->pEM)->setVariance(newVariance);}}
    virtual inline void         setVariance(short mixture, short dim, double value) {if (this->fullCovariance == false) {double **emVar = ((SC_BaggenstossEMex*)this->pEM)->getVariance(); emVar[mixture][dim] = value;}}
		virtual inline void         setVariance(short mixture, double value) {if (this->fullCovariance == false) {double **emVar = ((SC_BaggenstossEMex*)this->pEM)->getVariance(); for (short d = 0; d < this->dim; d++) {emVar[mixture][d] = value;}}}
    virtual inline void         setSd(double** newSd) {return;}
    virtual inline void         setSd(short mixture, short dim, double value) {return;}
    virtual inline void         setSd(short mixture, double value) {return;}
    virtual inline void         setMean(double** newMean) {this->pEM->setMean(newMean);}
    virtual inline void         setMean(short mixture, short dim, double value) {this->pEM->getMean()[mixture][dim] = value;}
    virtual inline void         setMean(short mixture, double value) {for (short d = 0; d < this->dim; d++) {this->pEM->getMean()[mixture][d] = value;}}
    virtual inline void         setWeight(double* newWeight) {this->pEM->setWeight(newWeight);}
    virtual inline void         setWeight(short mixture, double value) {this->pEM->getWeight()[mixture] = value;}

    virtual SC_Signature* toSignature(void);
    virtual void          killSignature(SC_Signature *pSignature);

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
    // this method follows the instructions in "A Monte-Carlo Method for Score Normalization in Automatic Speaker 
    // Verification using Kullback-Leibler Distances" by Ben, Blouet, Bimbot on ICASSP 2002
	  //====================================================================================================================
    virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count);
};
#endif
