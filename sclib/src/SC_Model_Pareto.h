/**************************************************************************/
/*    Intended to model a speaker's voice's probability density via the   */
/*    non-parametric Pareto Density Estimation technique by Alfred Ultsch */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.02.2006																								*/
/**************************************************************************/

#ifndef __SC_Model_Pareto_H__
#define __SC_Model_Pareto_H__

#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
#include "SC_Model.h"
#include "SC_PDE.h"

class SC_Model_Pareto : public SC_Model {

  private :

  protected :

		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

    SC_Model_Pareto *pBackground;
    SC_PDE *pPDEsolver; //to compute the PDE
    
    SV_Data **pDensities; //1 SV_Data-container per dim: the last col contains the density-estimation, all previous cols contain the original data
    double *paretoRadius; //1 entry per dim
    unsigned int dim; //=1 if useMarginalDistributions=false, otherwise tells how many independant cols there are
    bool useMarginalDistributions; //if true, each col of the training-/test-data is handled independantly; then dim=col, otherwise dim=1
    
  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_Model_Pareto(SC_TweakableParameters *pTweak = NULL, SC_Model_Pareto *pBackgound = NULL, bool useMarginalDistributions = true);
    SC_Model_Pareto(const SC_Model_Pareto& pParent);
	  virtual ~SC_Model_Pareto();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_Model_Pareto& operator=(const SC_Model_Pareto& pParent);

    //====================================================================================================================
    // get & set (new) background model; 
    //====================================================================================================================
    virtual SC_Model*   getBackground(void) {return this->pBackground;};
    virtual void		    setBackground(SC_Model* pBackground) {this->pBackground = (SC_Model_Pareto*)(pBackground); return;};
  
    inline unsigned int getDim(void) {return this->dim;}
    inline bool         getMarginalDistributionUsage(void) {return this->useMarginalDistributions;}
    inline double*      getParetoRadius(void) {return this->paretoRadius;}
    inline double       getParetoRadius(unsigned int dim) {return (dim < this->dim) ? this->paretoRadius[dim] : SVLIB_Fail;}
    inline SV_Data**    getDensities(void) {return this->pDensities;}
    inline SV_Data*     getDensities(unsigned int dim) {return (dim < this->dim) ? this->pDensities[dim] : NULL;}
    SV_Data*            getTrainingData(void);

    //====================================================================================================================
    // Combine 2 Models by adding the mixtures of this and the second one and return a new model
    //====================================================================================================================
    virtual SC_Model*   combineModels(SC_Model* pSecond);
    virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);
    SC_Model_Pareto*    combineModels(SC_Model_Pareto* pSecond);
    SC_Model_Pareto*	  combineModels(SC_Model_Pareto* pFirst, SC_Model_Pareto* pSecond, bool keepFirstsNext = false);

    //====================================================================================================================
    // Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
    //====================================================================================================================
    virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL);
    SC_Model_Pareto*    combineModels(SC_Model_Pareto* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model_Pareto* pBackgroundModels = NULL);

    //====================================================================================================================
    // methods to train & test a model with specified nr. of segments in the linked list of feature-vectors
    //====================================================================================================================
	  virtual int         TrainModel(SV_Data *TrainData);
		virtual int		      TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
	  virtual SV_Data*		TestModel(SV_Data *TestData);
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge);

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
    virtual unsigned int getFreeParameterCount(void) {return (this->dim > 0) ? this->dim * this->pDensities[0]->Row : 0;}

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
    virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count);
};

#endif
