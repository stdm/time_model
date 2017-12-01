/**************************************************************************/
/*    Derived from:																												*/
/*      - SV_Model as a base-class for new speaker models									*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.02.2006																								*/
/**************************************************************************/

#ifndef __SC_Model_H__
#define __SC_Model_H__

#include <iostream>
#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
#include "SC_Signature.h"
#include "SC_Api.h"
#include <SV_Error.h>
#include <SV_Model.h>
#include <SV_Data.h>

class SC_Model : public SV_Model {

  private :

  protected :
	  
	  //====================================================================================================================
	  // Last fileName passed successfully to OpenFile
	  //====================================================================================================================
		char lastUsedFileName[sclib::bufferSize];

	  //====================================================================================================================
	  // Model's parameters
	  //====================================================================================================================
		unsigned long int		trainingDataCount;			//count of feature-vectors used for training

    //====================================================================================================================
	  // There are some tweakable parameters in the SC_Lib library; they can be centraly managed in this class.
	  //====================================================================================================================
    SC_TweakableParameters* pTweak;

	  //====================================================================================================================
	  // Dump model's parameter in ASCII in a virtual fashion
	  //====================================================================================================================
	  virtual ostream& modelOut(ostream& os) = 0;

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
	  SC_Model(SC_TweakableParameters *pTweak = NULL);
		SC_Model(const SC_Model& pParent);
	  virtual ~SC_Model();

	  //====================================================================================================================
	  // assignment operator
	  //====================================================================================================================
		SC_Model& operator=(const SC_Model& pParent);

		//====================================================================================================================
		// opens the file just as the svlib version, but manages to remember the given filename in a new class member
		//====================================================================================================================
		void OpenFile(const char *FName, int Mode);

	  //====================================================================================================================
	  // give access to protected members
    // attention: no tests here wether the field indices really exist!
	  //====================================================================================================================
    inline unsigned long int  getTrainingDataCount(void) {return this->trainingDataCount;}
    inline void               setTrainindDataCount(unsigned long int newTrainingDataCount) {this->trainingDataCount = newTrainingDataCount;}
    inline SC_TweakableParameters* getTweak(void) {return this->pTweak;}
    
    //====================================================================================================================
	  // get & set (new) background model; 
	  //====================================================================================================================
		virtual SC_Model*		getBackground(void) = 0;
		virtual void				setBackground(SC_Model* pBackground) = 0;

		//====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond) = 0;
    virtual SC_Model*	  combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false) = 0;

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model*   combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge = 0, SC_Model* pBackgroundModels = NULL) = 0;
    
	  //====================================================================================================================
	  // methods to train/test a model (maybe with specified nr. of segments in the linked list of feature-vectors)
	  //====================================================================================================================
	  virtual int         TrainModel(SV_Data *TrainData) = 0;
		virtual int				  TrainModel(SV_Data *pData, unsigned long int segmentsToMerge) = 0;
	  virtual SV_Data*		TestModel(SV_Data *TestData) = 0;
		virtual SV_Data*		TestModel(SV_Data *TestData, unsigned long int segmentsToMerge) = 0;

	  //====================================================================================================================
	  // load/save a model from/to file
	  //====================================================================================================================
	  virtual SV_Model*		LoadModel(void) = 0;
    virtual int					SaveModel(void) = 0;

	  //====================================================================================================================
	  // create a (linked, not copied) SC_Signature-view on this model (for distance computation ), and destruct it
	  //====================================================================================================================
    virtual SC_Signature* toSignature(void) = 0;
    virtual void          killSignature(SC_Signature *pSignature) = 0;

	  //====================================================================================================================
	  // Dump model's parameter in ASCII 
	  //====================================================================================================================
	  friend  ostream& operator<<(ostream& os, SC_Model& Data);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) = 0;

	  //====================================================================================================================
	  // draw count samples distributed according to the pdf modeled by this model using monte-carlo methods
	  //====================================================================================================================
    virtual SV_Data*    drawSamplesFromDistribution(unsigned long int count) = 0;

	  //====================================================================================================================
	  // Give access to the internal file-variable used for loading/saving (necessary for seekg() calls if more than one 
		// model resides in a file)
	  //====================================================================================================================
		fstream* getDFile(void) {return &(this->DFile);}

	  //====================================================================================================================
	  // generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
		// method needs to return false;
	  //====================================================================================================================
		virtual bool scoreIsAverage(void) {return true;}
};

#endif
