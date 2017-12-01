/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Model to represent a quasi-GMM as described in		  					*/
/*				'Real-Time Unsupervised Speaker Change Detection',							*/
/*				 L.Lu, H.J.Zhang, 2002 (IEEE)																		*/
/*																																				*/
/*    How it works:																												*/
/*			A GMM-1 is estimated from a initial chunk of trainig data; as			*/
/*			more training-data gets available, the estimation of the cov of		*/
/*			GMM-1 can be imporved; if after updating no improvement is found,	*/
/*			the model evolves to GMM-2, that means, a new cov for a new				*/
/*			mixture is created. this emerges till MAX_MIXTURES mixtures are		*/
/*			reached or no more trainingdata are available.										*/
/*																																				*/
/*			because the model is meant to be used for speaker-change-					*/
/*			detection, where few data is available, it yields methods to			*/
/*			compute the match-score between a model and a cov instead of			*/
/*			a model and some feature-vectors.																	*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 13.03.2004																								*/
/**************************************************************************/

#ifndef __SC_Model_qGMM_H__
#define __SC_Model_qGMM_H__

//#include "SC_TweakableParameters.h"
#include "SC_DistanceMeasures.h"
#include "SC_Model.h"
#include "SC_MatrixFunctions.h"
//#include <SV_Model.h>

class SC_Model_qGMM : public SC_Model {

  private :

  protected :
	  
		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

	  //====================================================================================================================
	  // Model's parameters
	  //====================================================================================================================
		unsigned short int maxMixtures;       //max. number of mixture-components
    unsigned short int dim;							  //dimension of feature-vectors
	  unsigned short int mixtureCount;		  //count of mixture-components
	  double* weight;		                    //mixture weight vector
	  double** mean;			                  //array of mean-vectors
    unsigned long int* trainingDataCount; //count of feature-vectors used to estimate cov
	  double*** cov;	                  		//array of covariance-matrices of model
	  bool complete;								        //true, if MAX_MIXTURES mixtures are estimated good enough so that the model needs no further refinement

    //====================================================================================================================
	  // For dealing with the cov-matrices
	  //====================================================================================================================
    SC_MatrixFunctions* pMatrixFunc;

    //====================================================================================================================
	  // There are some tweakable parameters in the SC_Lib library; they can be centraly managed in this class.
	  //====================================================================================================================
    SC_TweakableParameters* pTweak;

    //====================================================================================================================
	  // For computing divergence shape distance between cov's
	  //====================================================================================================================
    SC_DistanceMeasures* pDist;

    //====================================================================================================================
	  // return the sum of the trainingDataCount[i]
	  //====================================================================================================================
    unsigned long int getCompleteTrainingDataCount(void);

	  //====================================================================================================================
	  // combine a given cov-matrix with new data to the cov resulting from estimating a cov out of the combination of 
	  // the original and the new training data (using a method by p.m.baggenstoss)
	  //====================================================================================================================
    double** combineCov(double** cov1, double** cov2, double* mean1, double* mean2, double weight1, double weight2, unsigned short int dim);
    
    //====================================================================================================================
	  // combine two covariance-matrices using the method of lu and zhang
	  //====================================================================================================================
    double** combineCov(double** cov_1, double** cov_2, unsigned long int n_1, unsigned long int n_2, unsigned short int dim);

	  //====================================================================================================================
	  // the model-training-method needs a threshold to decide wether the actual mixture is estimated well enough or not;
	  // this method computes the divergenceShape-distance between the input cov and a another cov, which's values differ
	  // from the first one's about percentSimilarity percent. this can serve as the threshold.
	  //====================================================================================================================
	  double getThreshold(double** cov, double percentSimilarity);

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_Model_qGMM(SC_TweakableParameters* pTweak, unsigned short int dim, unsigned short int maxMixtures);
	  virtual ~SC_Model_qGMM();

	  //====================================================================================================================
	  // give access to protected members
    // attention: no tests here wether the field indices really exist!
	  //====================================================================================================================
		inline unsigned short int	getMixtureCount(void) {return this->mixtureCount;}
		inline unsigned short int	getDim(void) {return this->dim;}
		inline double*						getWeight(void) {return this->weight;}
	  inline double							getWeight(unsigned short int mixture) {return this->weight[mixture];}
		inline double**						getMean(void) {return this->mean;}
    inline double*						getMean(unsigned short int mixture) {return this->mean[mixture];}
	  inline double							getMean(unsigned short int mixture, unsigned short int dim) {return this->mean[mixture][dim];}

    inline void         setMean(double** newMean) {MFree_2D(this->mean); this->mean = newMean;}
    inline void         setMean(short mixture, short dim, double value) {this->mean[mixture][dim] = value;}
    inline void         setMean(short mixture, double value) {for (short d = 0; d < this->dim; d++) {this->mean[mixture][d] = value;}}
    inline void         setWeight(double* newWeight) {MFree_1D(this->weight); this->weight = newWeight;}
    inline void         setWeight(short mixture, double value) {this->weight[mixture] = value;}
    inline void         setMixtureCount(unsigned short int newMixtureCount) {this->mixtureCount = newMixtureCount;}

	  //====================================================================================================================
	  // methods to estimate a model and to compute scores between models
	  //====================================================================================================================
	  virtual int TrainModel(SV_Data *TrainData);
    virtual int	TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
	  double TestModel(double** testCov, double* testMean, unsigned long int testDataCount);
    virtual SV_Data* TestModel(SV_Data *TestData); //use the (2) above method!
    virtual SV_Data* TestModel(SV_Data *TestData, unsigned long int segmentsToMerge); //use the (3) above method!

	  //====================================================================================================================
	  // method to provide acces to private member
	  //====================================================================================================================
    double** getLastCov(void) {return this->cov[this->mixtureCount-1];}
    double* getLastMean(void) {return this->mean[this->mixtureCount-1];}
    unsigned long int getLastTrainingDataCount(void) {return this->trainingDataCount[this->mixtureCount-1];}

	  //====================================================================================================================
	  // output stuff
	  //====================================================================================================================
	  virtual int SaveModel(void);
	  virtual SV_Model* LoadModel(void);

    //====================================================================================================================
	  // get & set (new) background model; 
	  //====================================================================================================================
		virtual SC_Model*		getBackground(void) {return NULL;}
		virtual void				setBackground(SC_Model* pBackground) {return;};

		//====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond) {return NULL;};
		virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false) {return NULL;};

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL) {return NULL;}

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
		virtual SC_Signature* toSignature(void) {return NULL;}
		virtual void          killSignature(SC_Signature *pSignature) {return;}

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
		virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count) {return NULL;}

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) {return this->mixtureCount*this->dim + this->mixtureCount*this->dim*this->dim + this->mixtureCount;} //mean & covariance & weight per dimension and mixture

	  //====================================================================================================================
	  // generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
		// method needs to return false;
	  //====================================================================================================================
		virtual bool scoreIsAverage(void) {return false;}
};

#endif
