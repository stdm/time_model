/**************************************************************************/
/*    Responsibility:																											*/
/*      - conatins a maximum-likelihood classifier based on one of the    */
/*        implemented SC_Model classes                                    */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Classifier_ML_H__
#define __SC_Classifier_ML_H__

#include "SC_Classifier.h"
#include "SC_Conversion.h"
#include "SC_Model.h"

class SCLIB_API SC_Classifier_ML : public SC_Classifier {
  
  private :

  protected :

    class SC_ClassMapping {
      public: 
        SC_ClassMapping(void) {this->Next = NULL;}
        int classLabel;
        int modelIndex;
        SC_ClassMapping *Next;
        int Valid(void) {return 1;}
    };

    SC_Classifier_ML::SC_ClassMapping *pMapping; //for mapping between model-index and original class-label
    SC_Model **pModels; //array of statstical models for the trained classes; two entrys in the two-class case
    unsigned int modelCount;

    //proved mappings between internal calssifiers (models) and class-labels
    bool addMapping(int classLabel, int modelIndex);

		//====================================================================================================================
    //	mappings between class labels and indices (into the probabiloty-container, e.g.)
    //====================================================================================================================
    int getMappedLabel(int modelIndex);
    int getMappedIndex(int classLabel);

		//====================================================================================================================
    //	free the memory of the models
    //====================================================================================================================
    void freeModels(void);

  public :

    SC_Classifier_ML(SC_TweakableParameters* pTweak, bool doScaling = false);
    virtual ~SC_Classifier_ML();

    //====================================================================================================================
    //	train a classifier for distinguishing between two classes, for which examples are given in the two SV_Data objects
    //====================================================================================================================
    virtual int trainTwoClass(SV_Data *pPositive, SV_Data *pNegative);

    //====================================================================================================================
    //	train a classifier for distinguishing between several classes
    //  the complete training-data (for all classes) is given in the SV_Data container, while the class-labes are given in 
    //  the classes-array, which has as many entrys as there are rows in pData, each entry corresponding with the 
    //  respective row of pData.
    //====================================================================================================================
    virtual int trainMultiClass(SV_Data *pData, int *classes);

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
};

#endif
