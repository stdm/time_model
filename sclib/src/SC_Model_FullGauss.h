/**************************************************************************/
/*    A wrapper around SV_Lib's SV_Model_Gaus to implement a full-        */
/*    covariance multivariate gaussian model (no mixture-model!).         */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.05.2008																								*/
/**************************************************************************/

#ifndef __SC_Model_FullGauss_H__
#define __SC_Model_FullGauss_H__

#include "SC_Model.h"
#include <SV_Model_Gaus.h>

class SC_Model_FullGauss : public SC_Model {

  private :

  protected :

		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

		SV_Model_Gaus *pGauss; //the real thing

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_Model_FullGauss(SC_TweakableParameters *pTweak);
    SC_Model_FullGauss(SC_Model_FullGauss& pParent);
	  virtual ~SC_Model_FullGauss();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_Model_FullGauss& operator=(SC_Model_FullGauss& pParent);

    //====================================================================================================================
    // get & set (new) background model (not existent here, so give dummies to fullfill the interface)
    //====================================================================================================================
    virtual SC_Model*   getBackground(void) {return NULL;};
    virtual void		    setBackground(SC_Model* pBackground) {return;};

    //====================================================================================================================
    // Combine 2 Models
    //====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond) {return NULL;} //this model-type needs reestimation!
		virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!
		SC_Model_FullGauss*	combineModels(SC_Model_FullGauss* pSecond) {return NULL;} //this model-type needs reestimation!
		SC_Model_FullGauss*	combineModels(SC_Model_FullGauss* pFirst, SC_Model_FullGauss* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!

    //====================================================================================================================
    // Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
    //====================================================================================================================
    virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL);
    SC_Model_FullGauss*	combineModels(SC_Model_FullGauss* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model_FullGauss* pBackgroundModels = NULL);

    //====================================================================================================================
    // methods to train & test a model with specified nr. of segments in the linked list of feature-vectors
		// here, the segments to use (merge) are NOT merged but used as distinct temporal patterns to train the hmm; 
		// segmentsToMerge patterns are actually used.
    //====================================================================================================================
	  virtual int         TrainModel(SV_Data *TrainData);
		virtual int		      TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
	  virtual SV_Data*		TestModel(SV_Data *TestData);

		//====================================================================================================================
		// Test the model with specified nr. of segments in the linked list of feature-vectors.
		// The return value is a matrix with one column and T+1 rows where T is the number of test feature vectors: the 0th
		// row holds the score as defined in Yegnanarayana & Kishore's paper on AANN as an alternative to GMMs, the following 
		// rows include each individual vector's confidence score as proposed in Dhanajaya & Yegnanarayana's paper on speaker 
		// change detection.
		//====================================================================================================================
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge);

    //====================================================================================================================
    // load/save a model from/to file
    //====================================================================================================================
    virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
		virtual SC_Signature* toSignature(void) {return NULL;} //a full-gauss model has no sound signature representation yet
		virtual void          killSignature(SC_Signature *pSignature) {return;}

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void);

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
		virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count); 
};

#endif
