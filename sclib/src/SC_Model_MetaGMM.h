/**************************************************************************/
/*    This is a "meta-GMM": it holds one special GMM per dimension, so    */
/*    that each dimension can have a speratae mixture count. Maybe this   */
/*    helps to save some parameters.                                      */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 07.09.2009																								*/
/**************************************************************************/

#ifndef __SC_Model_MetaGMM_H__
#define __SC_Model_MetaGMM_H__

#include "SC_Model.h"
#include "SC_MixtureModel_GMM.h"

class SC_Model_MetaGMM : public SC_Model {

  private :

  protected :

	  //====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

		unsigned short int dim;
		unsigned short int maxMixtureCount;
		unsigned short int *mixtureCount;
		SC_MixtureModel_GMM **pGMM;
		double **orthTransMat;

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
		SC_Model_MetaGMM(SC_TweakableParameters* pTweak, unsigned short int maxMixtureCount, unsigned short int dim);
		SC_Model_MetaGMM(const SC_Model_MetaGMM& pParent);
	  virtual ~SC_Model_MetaGMM();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_Model_MetaGMM& operator=(const SC_Model_MetaGMM& pParent);

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

    //====================================================================================================================
    // Combine 2 Models
    //====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond) {return NULL;} //this model-type needs reestimation!
		virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!
		SC_Model_MetaGMM*  	combineModels(SC_Model_MetaGMM* pSecond) {return NULL;} //this model-type needs reestimation!
		SC_Model_MetaGMM*	  combineModels(SC_Model_MetaGMM* pFirst, SC_Model_MetaGMM* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL);
    SC_Model_MetaGMM* combineModels(SC_Model_MetaGMM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model_MetaGMM* pBackgroundModels = NULL);

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
		virtual SC_Signature* toSignature(void);
		virtual void killSignature(SC_Signature *pSignature);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void);

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
		virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count); 

	  //====================================================================================================================
	  // generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
		// method needs to return false;
	  //====================================================================================================================
		virtual bool scoreIsAverage(void);

	  //====================================================================================================================
	  // load/save a model from/to file
	  //====================================================================================================================
	  virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);
};

#endif
