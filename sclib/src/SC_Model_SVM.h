/**************************************************************************/
/*    A model that doesn't solve the (hard) problem of estimating the pdf */
/*    of the data but instead the simpler problem of estimating the level */
/*    set (support or center of mass of the pdf) via a one-class nu-SVM   */
/*    as described in                                                     */
/*      "An Online Kernel Change Detection Algorithm", Desobry, Davy,     */
/*      Doncarli, IEEE Trans. Sig. Proc., Vol 53, No. 8, 2005,            */
/*      pp. 2961-2974                                                     */
/*    and                                                                 */
/*      "Unsupervised Speaker Indexing Using One-Class Support Vector     */
/*      Machines", Fergani, Davy, Houacine, 2006                          */
/*                                                                        */
/*    The model is typically evaluated via computing a distance (proposed */
/*    in the same papers) to a second model of same type; such, when      */
/*    evaluated giving only a dataset, this is first modeled.             */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 17.04.2007																								*/
/**************************************************************************/

#ifndef __SC_Model_SVM_H__
#define __SC_Model_SVM_H__

#include "SC_Model.h"
#include "SC_Classifier_SVM.h"

class SC_Model_SVM : public SC_Model {

  private :

  protected :

		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

    SC_Classifier_SVM *pSVM; //to compute and store the SVM model
    unsigned int dim; //dimensionality of training feature vectors
		bool verbose; //shall the model be talkative, i.e. print things to stdout?
		bool distanceBasedTesting; //if true, the arc-distance (SC_DistanceMeasures::svmArcDistance()) is used in TestModel() instead of SVM classification result
		bool doParameterSearch; //if true, the SVM searches for a good gamma parameter using the training data; not a good choice if distanceBasedScoring is true!
    
  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_Model_SVM(SC_TweakableParameters *pTweak = NULL, bool distanceBasedTesting = true, bool doParameterSearch = false, bool verbose = false);
    SC_Model_SVM(const SC_Model_SVM& pParent, bool justLink = false);
	  virtual ~SC_Model_SVM();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_Model_SVM& operator=(const SC_Model_SVM& pParent);

    //====================================================================================================================
    // get & set (new) background model (not existent here, so give dummies to fullfill the interface)
    //====================================================================================================================
    virtual SC_Model*   getBackground(void) {return NULL;};
    virtual void		    setBackground(SC_Model* pBackground) {return;};
  
    inline unsigned int getDim(void) {return this->dim;}
		inline SC_SVMmodel* getSVMmodel(void) {return this->pSVM->getClassifier();}
		inline SC_Classifier_SVM* getSVM(void) {return this->pSVM;}
		inline void setDistanceBasedTesting(bool distanceBasedTesting) {this->distanceBasedTesting = distanceBasedTesting; return;}

    //====================================================================================================================
    // Combine 2 Models by adding the mixtures of this and the second one and return a new model
    //====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond) {return NULL;} //this model-type needs reestimation!
		virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!
		SC_Model_SVM*				combineModels(SC_Model_SVM* pSecond) {return NULL;} //this model-type needs reestimation!
		SC_Model_SVM*				combineModels(SC_Model_SVM* pFirst, SC_Model_SVM* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!

    //====================================================================================================================
    // Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
    //====================================================================================================================
    virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL);
    SC_Model_SVM*				combineModels(SC_Model_SVM* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model_SVM* pBackgroundModels = NULL);

    //====================================================================================================================
    // methods to train & test a model with specified nr. of segments in the linked list of feature-vectors
    //====================================================================================================================
	  virtual int         TrainModel(SV_Data *TrainData);
		virtual int		      TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
		virtual int					TrainModel(SC_SVMproblem *pProblem);
	  virtual SV_Data*		TestModel(SV_Data *TestData);
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge);
		virtual SV_Data*    TestModel(SC_Model_SVM *pModel); //scoring is always accomplished here by giving a distance between this model and a second model (possibly temporary created from a TestData dataset)

    //====================================================================================================================
    // load/save a model from/to file
    //====================================================================================================================
    virtual SV_Model*		LoadModel(void);
    virtual int					SaveModel(void);

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
		virtual SC_Signature* toSignature(void) {return NULL;} //a SVM has no sound signature representation
		virtual void          killSignature(SC_Signature *pSignature) {return;}

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) {return 2;} //free parameters: left and right boundary of the region of mapped feature-vectors to SVM kernel space

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
		virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count);

	  //====================================================================================================================
	  // generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
		// method needs to return false;
	  //====================================================================================================================
		virtual bool scoreIsAverage(void) {return false;}
};

#endif
