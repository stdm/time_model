/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_Model as a base-class for new mixture models									*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.04.2004																								*/
/**************************************************************************/

#ifndef __SC_MixtureModel_H__
#define __SC_MixtureModel_H__

#include <iostream>
#include "SC_Model.h"
#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
#include "SC_Centroid_Gaussian.h"
#include "SC_Gauss.h"
#include <SV_Error.h>
#include <SV_Model.h>
#include <SV_Data.h>

class SC_MixtureModel : public SC_Model {

  private :

  protected :

		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

	  //====================================================================================================================
	  // Model's parameters
	  //====================================================================================================================
	  unsigned short int	mixtureCount;						//count of mixture-components
		unsigned short int	dim;										//dimension of feature-vectors
	  double*							weight;		              //mixture weight vector
	  double**						variance;	          		//array of variance-vectors (see above for explanation)
    double**            sd;                     //the standard-deviation = sqrt(variance)
	  double**						mean;			              //array of mean-vectors
    unsigned long int   maxEMsteps;             //maximum number of EM iterations

		//====================================================================================================================
	  // The model-parameters must be initialized some way prior to the em-algorithm
	  //====================================================================================================================
    void								initParameters(SV_Data *pData, double varianceLimit, unsigned int kMeansCount = 10);

		//====================================================================================================================
	  // To compute the sd (srqt(variance)) out of the given variance
	  //====================================================================================================================
    double   variance2sd(double variance);
    double*  variance2sd(double *variance, unsigned int dim);
    double** variance2sd(double **variance, unsigned int len, unsigned int dim);

		//====================================================================================================================
	  // Constants for dummy-mixture
	  //====================================================================================================================
		const double dummyWeight; //one mixture has the weight 1.0 ;-)
		const double dummyMean; //mean for a dummy mixture (should be very small)
		const double dummyVariance; //variance for a dummy mixture (should be near zero)

  public :

		//====================================================================================================================
	  // For the calculation of gaussian functions
	  //====================================================================================================================
    static SC_Gauss gaussSolver;

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_MixtureModel(SC_TweakableParameters *pTweak = NULL);
		SC_MixtureModel(const SC_MixtureModel &pParent);
	  virtual ~SC_MixtureModel();

	  //====================================================================================================================
	  // assignment operator
	  //====================================================================================================================
		SC_MixtureModel& operator=(const SC_MixtureModel& pParent);

	  //====================================================================================================================
	  // give access to protected members
    // attention: no tests here wether the field indices really exist!
	  //====================================================================================================================
		virtual inline unsigned short int	getMixtureCount(void) {return this->mixtureCount;}
		virtual inline unsigned short int	getDim(void) {return this->dim;}
		virtual inline double*						getWeight(void) {return this->weight;}
	  virtual inline double							getWeight(unsigned short int mixture) {return this->weight[mixture];}
		virtual inline double**						getVariance(void) {return this->variance;}
    virtual inline double*						getVariance(unsigned short int mixture) {return this->variance[mixture];}
    virtual inline double							getVariance(unsigned short int mixture, unsigned short int dim) {return this->variance[mixture][dim];}
		virtual inline double**						getSd(void) {return this->sd;}
    virtual inline double*						getSd(unsigned short int mixture) {return this->sd[mixture];}
    virtual inline double							getSd(unsigned short int mixture, unsigned short int dim) {return this->sd[mixture][dim];}
		virtual inline double**						getMean(void) {return this->mean;}
    virtual inline double*						getMean(unsigned short int mixture) {return this->mean[mixture];}
	  virtual inline double							getMean(unsigned short int mixture, unsigned short int dim) {return this->mean[mixture][dim];}
    virtual inline unsigned long int  getMaxEMsteps(void) {return this->maxEMsteps;}
    virtual double										getWeightLimit(void) = 0;

    virtual inline void         setVariance(double** newVariance) {MFree_2D(this->variance); this->variance = newVariance;}
    virtual inline void         setVariance(short mixture, short dim, double value) {this->variance[mixture][dim] = value;}
    virtual inline void         setVariance(short mixture, double value) {for (short d = 0; d < this->dim; d++) {this->variance[mixture][d] = value;}}
    virtual inline void         setSd(double** newSd) {MFree_2D(this->sd); this->sd = newSd;}
    virtual inline void         setSd(short mixture, short dim, double value) {this->sd[mixture][dim] = value;}
    virtual inline void         setSd(short mixture, double value) {for (short d = 0; d < this->dim; d++) {this->sd[mixture][d] = value;}}
    virtual inline void         setMean(double** newMean) {MFree_2D(this->mean); this->mean = newMean;}
    virtual inline void         setMean(short mixture, short dim, double value) {this->mean[mixture][dim] = value;}
    virtual inline void         setMean(short mixture, double value) {for (short d = 0; d < this->dim; d++) {this->mean[mixture][d] = value;}}
    virtual inline void         setWeight(double* newWeight) {MFree_1D(this->weight); this->weight = newWeight;}
    virtual inline void         setWeight(short mixture, double value) {this->weight[mixture] = value;}
    virtual inline void         setMixtureCount(unsigned short int newMixtureCount) {this->mixtureCount = newMixtureCount;}
    virtual inline void         setMaxEMsteps(unsigned long int maxSeps) {this->maxEMsteps = maxSeps; return;}
    
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
	  // get & set (new) background model; 
	  //====================================================================================================================
		virtual SC_Model*		getBackground(void) = 0;
		virtual void				setBackground(SC_Model* pBackground) = 0;

		//====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond) = 0;
    virtual SC_Model*   combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false) = 0;

		//====================================================================================================================
		// Combine 2 Models by estimating it on the combined training-feature-vectors of both parents
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL) = 0;
    
	  //====================================================================================================================
	  // methods to train/test a model (maybe with specified nr. of segments in the linked list of feature-vectors)
		// all mixture models return _average_ log likelihood values as suggested in "Speaker Verification Using Adapted 
		// Gaussian Mixture Models" to compensate for length effects
	  //====================================================================================================================
	  virtual int         TrainModel(SV_Data *TrainData) = 0;
		virtual int				  TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) = 0;
    virtual SV_Data*		TestModel(SV_Data *TestData) = 0;
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) = 0;
		virtual double			testModel(SV_Data *TestData, unsigned long int segmentsToMerge = 1) {SV_Data *pScore = TestModel(TestData, segmentsToMerge); double res = pScore->Mat[0][0]; MFree_0D(pScore); return res;}

		//====================================================================================================================
	  // Step trough a linked list of models (starting from this one) and return the greatest nr. of mixture-components,
		// that one single model in the list holds 
	  //====================================================================================================================
		unsigned int				getMaxMixturesInList(void);

	  //====================================================================================================================
	  // load/save a model from/to file
	  //====================================================================================================================
	  virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);
	  virtual SV_Model*		oldLoadModel(void); //uses the old load-method (machine-dependant, but there are still some files stored this way...)

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
    virtual SC_Signature* toSignature(void);
    virtual void          killSignature(SC_Signature *pSignature);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) = 0;

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
    // this method follows the instructions in "A Monte-Carlo Method for Score Normalization in Automatic Speaker 
    // Verification using Kullback-Leibler Distances" by Ben, Blouet, Bimbot on ICASSP 2002
	  //====================================================================================================================
    virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count);
};

#endif
