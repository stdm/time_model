/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates an SVM-classifier based on the libsvm-SVM-         */
/*        implementation v2.81 (with additions by Hsuan-Tien Lin to       */
/*        support weighted examples for c-svc and epsilon-svr)            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Classifier_SVM_H__
#define __SC_Classifier_SVM_H__

#include "SC_ClassifierWithWeights.h"
#include "SC_SVM.h"
#include "SC_Aux.h"

class SCLIB_API SC_Classifier_SVM : public SC_ClassifierWithWeights {
  
  private :

		class SC_Weight { //to form linked lists of weight-class-mappings
			public:
				SC_Weight(int label = sclib::labelPositive, double weight = 0.5) {
					this->label = label;
					this->weight = weight;
					this->Next = NULL;
				}
				SC_Weight(const SC_Weight& pParent) {
					this->label = pParent.label;
					this->weight = pParent.weight;
					this->Next = NULL;
				}
				int label;
				double weight;
				SC_Weight *Next;
				int Valid(void) {return 1;}
		};

  protected :
 
    SC_SVM *pSVM; //the SVM algorithms
    SC_SVMmodel *pClassifier; //a trained SVM-model
    SC_SVMproblem *pTraingData; //the original training data in SC_SVM-readable format (need to be stored as long as the trained model persists because it links to it)
    bool doParameterSearch; //controls whether to use default parameters or do a parameter search
		SC_Classifier_SVM::SC_Weight *pClassWeightList; //so that the user can specifiy weights per class
		bool justLinked; //true if the trained SVM and the training data is just linked during copying from a parent model

    //====================================================================================================================
    //	reduce a training-set by returning a new one with maxSize vectors drawn at random
    //====================================================================================================================
		SC_SVMproblem* reduceTrainingSet(SC_SVMproblem *pProblem, unsigned long int maxSize);

		//====================================================================================================================
    //	find suitable parameter settings according to the tweakable parameters
    //====================================================================================================================
    SC_SVMparameter* getParameters(SC_SVMproblem *pProblem = NULL);
  
    //====================================================================================================================
    //	check the given parameter-set for correctness/applicability
    //====================================================================================================================
    bool checkParameters(SC_SVMproblem *pPoblem, SC_SVMparameter *pParam);

    //====================================================================================================================
    //	incorporates the pClassWeight knowledge into a svm-parameterset
    //====================================================================================================================
		void classWeights2parameters(SC_SVMparameter *pParam);

  public :

    SC_Classifier_SVM(SC_TweakableParameters* pTweak, bool doParameterSearch = true, bool doScaling = true, bool verbose = true);
		SC_Classifier_SVM(const SC_Classifier_SVM& pParent, bool justLink = false);
    virtual ~SC_Classifier_SVM();

    //====================================================================================================================
    //	setter
    //====================================================================================================================
		void setDoParameterSearch(bool doParameterSearch) {this->doParameterSearch = doParameterSearch; return;}

    //====================================================================================================================
    //	train a one-class nu-SVM (more of a model than a classifier)
    //====================================================================================================================
		virtual int trainOneClass(SV_Data *pData);
		virtual int trainOneClass(SC_SVMproblem *pProblem);

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative);

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //  in the two-class-case, the constants SCIB_LABEL_POSITIVE/sclib::labelNegative should be used to indicate the 
    //  classes, and they normaly should evaluate to +1/-1, respectively
		//  the weights in the two arrays together (corresponding to rows in the SV_Data objects) are meant to sum up to 1
    //====================================================================================================================
		virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative, double *positiveWeights, double *negativeWeights);

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes);

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
		//  the weights (corresponding to rows in the SV_Data object) are meant to sum up to 1
    //====================================================================================================================
		virtual int trainMultiClass(SV_Data *pData, int *classes, double *weights);

    //====================================================================================================================
    //	classifiy previously unseen test-data; returned is an array of classlabels, each entry corresponding to the 
    //  respective row in pData; if available, the probabilities for each class-decision are given in the pProbabilities 
		//  parameter: the rows therein correspond to the pData-rowes, and the columns correspond to the classes
    //
    //  ATTENTION: only a previously trained (or loaded) classifier can be used for classification!
    //====================================================================================================================
    virtual int* classify(SV_Data *pData, SV_Data* &pProbabilities);
    
    //====================================================================================================================
    //	save and load a trained classifier to/from a file
    //====================================================================================================================
    virtual int saveClassifier(const char *fileName);
    virtual int loadClassifier(const char *fileName);

    //====================================================================================================================
    //	to convert between given labels and indices into the probability-parameter of the classifiy()-method
    //====================================================================================================================
		virtual long int label2idx(long int label);
		virtual long int idx2label(long int idx);

		//====================================================================================================================
		// below are methods to make the interface and data-exchange to and with libsvm easier...
		//====================================================================================================================

    //====================================================================================================================
    //	access to protected member (for the SC_Model_SVM class and if distances between two one-class SVMs should be 
		//  computetd)
    //====================================================================================================================
		inline SC_SVM* getSVM(void) {return this->pSVM;}
		inline SC_SVMmodel* getClassifier(void) {return this->pClassifier;}

    //====================================================================================================================
    //	convert a SV_Data dataset into a format that can be subject to SVM-training; the class-labels are meant to be 
    //  given in the classes-array per row in pData
    //====================================================================================================================
    SC_SVMproblem* svData2svmProblem(SV_Data *pData, int *classes, double *weights = NULL);

    //====================================================================================================================
    //	convert a positive and a negative example SV_Data dataset into a format that can be subject to SVM-training
    //  can also be used for one-class SVM training (with pNegative==NULL)
    //====================================================================================================================
    SC_SVMproblem* svData2svmProblem(SV_Data *pPositive, SV_Data *pNegative = NULL, double *positiveWeights = NULL, double *negativeWeights = NULL);

    //====================================================================================================================
    //	convert a single feature-vector out of an SV_Data container into a format that can be subject to SVM-testing
    //  do the scaling too, if wished
    //====================================================================================================================
    SC_SVMnode* svDataRow2svmNode(SV_Data *pData, unsigned long int row);

    //====================================================================================================================
    //	creates and returns a copy of one row of SVMnodes (termionated by a node with index -1)
    //====================================================================================================================
    SC_SVMnode* copySVMnodes(SC_SVMnode *pNodes);

    //====================================================================================================================
    //	destroy a SC_SVMproblem struct
    //====================================================================================================================
    static void killSVMproblem(SC_SVMproblem* &pProblem);

    //====================================================================================================================
    //	create a copy of a SC_SVMproblem struct
    //====================================================================================================================
		SC_SVMproblem* copySVMproblem(SC_SVMproblem* pProblem);

		//====================================================================================================================
		//  create a copy of a trained SVM model
		//====================================================================================================================
		SC_SVMmodel* copySVMmodel(SC_SVMmodel* pModel);

		//====================================================================================================================
		//  return a gamma-parameter used in one-class SVM computation that is optimal with respect to the given data
		//  this also creates a scale-matrix if scaling is switched on for this classifier
		//====================================================================================================================
		double getOneClassGammaParameter(SV_Data *pData);

		//====================================================================================================================
		//  return a gamma-parameter used in one-class SVM computation that is optimal with respect to the given data
		//====================================================================================================================
		double getOneClassGammaParameter(SC_SVMproblem *pProblem, SC_SVMparameter * pParam);

		//====================================================================================================================
		//  some methods to set weights per class; it is stored in an internal data structure and used whenever the 
		//  class-weightsd could be brought in
		//====================================================================================================================
		void setClassWeight(int classLabel, double weight);
		void removeClassWeight(int classLabel);
		void removeClassWeights(void);
};

#endif
