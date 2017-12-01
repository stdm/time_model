/**************************************************************************/
/*    This is the testbed for a new kind of model that also regards the   */
/*    time order of feature vectors in the given feature sets             */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.02.2009																								*/
/**************************************************************************/
 
#ifndef __SC_Model_Time_H__
#define __SC_Model_Time_H__

#include "SC_Model.h"
#include <map>
#include <string>

class SC_TestResultMap { //needed so there are no warnings for template members without dll interface below... (C4251)
	public:
		std::map<std::string, double> map;
};

class SC_Model_Time : public SC_Model {

  private :

  protected :

		//====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os);

		SC_Model *pSubModel; //the real model
		char subModelType; //model typ of the real model
		bool verbose; //shall the model be talkative, i.e. print things to stdout?
		unsigned int dim; //expected dimensionality of feature vectors
		unsigned int syllableLength; //length in [ms] of a "syllable", i.e. the local time horizon (trajectory) to reflect in the model
		unsigned int trajectoryStep; //frame shift between trajectories in [#frames]
		SV_Data *pNorm; //normalization matrix for pre-trajectory (i.e. original, i.i.d.) feature set
		bool removeTiming; //if true, timing information is removed in the succession of frames by clustering them into #templateCount templates, and replcing each frame with the nearest template, thereby discarding frames that fall unto the same template than their predecessor; this way, timing is removed by retaining spectral individuality of this feature-set
		unsigned int templateCount; //nr. of templates the data is clustered into when removing timing
		unsigned int clusteringIterations; //nr. of iterations of kMeans in clustering to find timing-removal templates
		bool replaceTrainingData; //if true, the given (list of) training data is replaced with its trajectory-version
		bool checkForTrajectorization; //if true, the given training-/test-data is first checked if already in trajectory-form; if it is (Hdr->signature[2]>0), it is directly used without trajectorizing it again
		bool cacheResults; //if true, results on datasets are cached and the scoring is speeded up; then, no detailed results are returned but just the overall likelihood of the dataset
		SC_TestResultMap checksumCache; //tested complete data sets get a checksum that is associated with the their likelihood that can be used to do faster scoring in conjunction the e.g. the CLR measure
		static SC_Model *pWorldModel; //world-model of type subModelType to remove common parts of speech
		static int instanceCounter; //nr. of instances currently created; to kill world model on last instance destruction
		
  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
		SC_Model_Time(SC_TweakableParameters *pTweak, unsigned int dim, unsigned int syllableLength = 100, unsigned int trajectoryStep = 1, char subModelType = sclib::mtSVM, bool removeTiming = true, unsigned int templateCount = 0, unsigned int clusteringIterations = 100, bool replaceTrainingData = false, bool checkForTrajectorization = true, const char *worldModelFile = "", const char *normalizationFile = "", bool cacheResults = true, bool verbose = true);
    SC_Model_Time(SC_Model_Time& pParent, bool justLink = false);
	  virtual ~SC_Model_Time();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_Model_Time& operator=(SC_Model_Time& pParent);

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
		SC_Model_Time*	combineModels(SC_Model_Time* pSecond) {return NULL;} //this model-type needs reestimation!
		SC_Model_Time*	combineModels(SC_Model_Time* pFirst, SC_Model_Time* pSecond, bool keepFirstsNext = false) {return NULL;} //this model-type needs reestimation!

    //====================================================================================================================
    // Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
    //====================================================================================================================
    virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL);
    SC_Model_Time*	combineModels(SC_Model_Time* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model_Time* pBackgroundModels = NULL);

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
		virtual SC_Signature* toSignature(void) {return NULL;} //no sound signature representation yet
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
	  // generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
		// method needs to return false;
	  //====================================================================================================================
		virtual bool scoreIsAverage(void);

	  //====================================================================================================================
	  // getter for protected members
	  //====================================================================================================================
		SC_Model* getSubModel(void) {return this->pSubModel;}
};

#endif
