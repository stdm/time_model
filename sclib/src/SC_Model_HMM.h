/**************************************************************************/
/*    A wrapper around SV_Lib's SV_Model_CHMM to implement a continous    */
/*    density Hidden Markov Model for kinds of layouts (i.e. ergodic/     */
/*    fully-connected, left-to-right/Bakis, you name it). The code in     */
/*    SV_Model_CHMM was slightly changed to incorporate the initial state */
/*    probability to allow for this flexibility. Now, the implementation  */
/*    closely follows Rabiner, "A Tutorial on Hidden Markov Models and    */
/*    Selected Applications in Speech Recognition", 1989. All additional  */
/*    sugar is implemented here: compliance with the SC_Model interface,  */
/*    including sampling from the model, and a useable interface to best- */
/*    path estimation via the Viterbi algorithm.                          */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 30.04.2008																								*/
/**************************************************************************/

#ifndef __SC_Model_HMM_H__
#define __SC_Model_HMM_H__

#include "SC_Model.h"
#include <SV_Model_CHMM.h>

class SC_Model_HMM : public SC_Model {

  private :

  protected :

		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

		//====================================================================================================================
		// Parses a string that represents a boolean matrix (columns delimited by a space ' ', rows delimited by a semicolon
		// ';', true by '1', false by '0'). The entries in row i, column j in this matrix mean if the hmm-state i should have
		// a link (state transition probability possibly > 0) to state j or not. Returned is the reconstructed matrix with
		// integer datatype as needed by the SV_Model_CHMM class for the ConfMat member. if the string is invalid or "" or 
		// NULL, a fully-connected model is assumed or too short or missing rows are filled up with '1'.
		//====================================================================================================================
		int** parseTransitionStructureString(const char *transitionStructure, unsigned int stateCount);

		//====================================================================================================================
		// The other way round
		//====================================================================================================================
		char* constructTransitionStructureString(int **transitionMatrix, unsigned int stateCount);
		
		SV_Model_CHMM *pHMM; //the real thing

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_Model_HMM(SC_TweakableParameters *pTweak, unsigned int stateCount, const char *transitionStructure, unsigned int mixturesPerState, bool useOrthogonalTransform = false, bool leftToRight = false, unsigned int maxIterations = 100, bool verbose = false);
    SC_Model_HMM(SC_Model_HMM& pParent);
	  virtual ~SC_Model_HMM();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_Model_HMM& operator=(SC_Model_HMM& pParent);

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
		SC_Model_HMM*				combineModels(SC_Model_HMM* pSecond) {return NULL;} //this model-type needs reestimation!
		SC_Model_HMM*				combineModels(SC_Model_HMM* pFirst, SC_Model_HMM* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!

    //====================================================================================================================
    // Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
    //====================================================================================================================
    virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL);
    SC_Model_HMM*				combineModels(SC_Model_HMM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model_HMM* pBackgroundModels = NULL);

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
		virtual SC_Signature* toSignature(void) {return NULL;} //a HMM has no sound signature representation
		virtual void          killSignature(SC_Signature *pSignature) {return;}

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void);

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
		virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count); 

	  //====================================================================================================================
	  // uses the viterbi-algorithm to find the best state sequence through the previously trained model. an int-array with 
		// rows corresponding to the pFeatures-rows is returned, where each entry gives a state-number between 0 and 
		// stateCount-1
	  //====================================================================================================================
		int* getBestStateSequence(SV_Data *pFeatures, double &logLikelihood);
};

#endif
