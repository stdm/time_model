/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the change        */
/*        detction (speaker and acoustic) process													*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.11.2006																								*/
/**************************************************************************/

#ifndef __SC_Score_ChangeDetection_H__
#define __SC_Score_ChangeDetection_H__

#include "SC_Score.h"

class SCLIB_API SC_Score_ChangeDetection : public SC_Score {
	private :
    
		class SC_Boundary {
			public:
				SC_Boundary(unsigned long int boundary, SC_Boundary* Next = NULL) {
					this->boundary = boundary; 
					this->Next = Next;
				};
				unsigned long int boundary;
				SC_Boundary *Next;
				int Valid(void) {return 1;}
		};

  protected :
    
		//members to hold scores of the classification process
		unsigned long int ca; //number or count of classes
		long int ***scatterMatrix[4];
		double *recall[4];
		double *precision[4];
		double *missRate[4];
		double *falseAlarmRate[4];
		double *errorRate[4];
		double *specificity[4];
		double *accuracy[4];

		//storage for the raw data boundary data according to it's classification as gt-based or TP, FP or FN (TN is uninteresting...)
		SC_Score_ChangeDetection::SC_Boundary *gtBoundaryList[2];
		SC_Score_ChangeDetection::SC_Boundary *hypoTPlist[2];
		SC_Score_ChangeDetection::SC_Boundary *hypoFPlist[2];
		SC_Score_ChangeDetection::SC_Boundary *hypoFNlist[2];

		//====================================================================================================================
		//  Add a boundary to the given list (see above), so that it can be outputted when printing the report
		//====================================================================================================================
		void addBoundary(SC_Score_ChangeDetection::SC_Boundary* &boundaryList, unsigned long int boundary);

		//====================================================================================================================
		//  To make operator<<() kind of virtual...
		//====================================================================================================================
		virtual ostream& output(ostream& OutS); 

		//====================================================================================================================
		//  This method returns for a given class index the types (or typesNot, according to "typesNot"-parameter) that 
		//  correspond with the analyzed parts of the signal for this class
		//  TODO: this should probably be placed into the algorithm's class itself and given to this class as a parameter for 
		//  better decoupling...
		//====================================================================================================================
		long int getRelatedTypes(long int idx, bool typesNot = false);

	public :

    //====================================================================================================================
    //	The constructor
    //====================================================================================================================
    SC_Score_ChangeDetection(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT);

    //====================================================================================================================
    //	The destructor
    //====================================================================================================================
    virtual ~SC_Score_ChangeDetection();

    //====================================================================================================================
    //	Fills the internal score-variables by computing their values according to the frameList in pGT, so that the
    //  get*()-Functions return reasonable values (before calling calcScores(), they return all 0)
		//  "start" and "end" refer to sample-numbers so that the area of the frameList for which scores shall be computed can
		//  be specified; this way, scores can be calculated only for parts of the video/corpus, e.g. for a scene.
		//  In algorithmicUncertaintyDiameter a value [in samples!] can be given which describes the precision with which the
		//  specific algorithm responsible for the results (and only the algorithm, not the gt...) can predict the place of 
		//  event-on- and -offsets
    //====================================================================================================================
    virtual void calcScores(unsigned long int start = 0, unsigned long int end = 0, unsigned long int algorithmicUncertaintyDiameter = 0); 

    //====================================================================================================================
    //	These get*()-functions give access to the results computed by calcScores()
		//  There are different names for the same measure, so there are different methods providing the same result
    //====================================================================================================================
		virtual double getRecall(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->recall[ce] != NULL) ? this->recall[ce][classIdx] : 0.0;}
		virtual double getPrecision(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->precision[ce] != NULL) ? this->precision[ce][classIdx] : 0.0;}
		virtual double getF1(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx<this->ca && this->recall[ce]!=NULL && this->precision[ce]!=NULL && this->recall[ce][classIdx]*this->precision[ce][classIdx]>0.0) ? (2.0*this->recall[ce][classIdx]*this->precision[ce][classIdx])/(this->recall[ce][classIdx]+this->precision[ce][classIdx]) : 0.0;}
		virtual double getMissRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->missRate[ce] != NULL) ? this->missRate[ce][classIdx] : 0.0;}
		virtual double getMDR(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getMissRate(classIdx,countedEntity);}
		virtual double getFalseAlarmRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->falseAlarmRate[ce] != NULL) ? this->falseAlarmRate[ce][classIdx] : 0.0;}
		virtual double getFAR(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getFalseAlarmRate(classIdx, countedEntity);}
		virtual double getErrorRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->errorRate[ce] != NULL) ? this->errorRate[ce][classIdx] : 0.0;}
		virtual double getSpecificity(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->specificity[ce] != NULL) ? this->specificity[ce][classIdx] : 0.0;}
		virtual double getAccuracy(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->accuracy[ce] != NULL) ? this->accuracy[ce][classIdx] : 0.0;}
		virtual double getSensitivity(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getRecall(classIdx, countedEntity);}
		virtual double getTruePositiveRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getRecall(classIdx, countedEntity);}
		virtual double getOmission(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return 1.0 - getRecall(classIdx, countedEntity);}
		virtual double getProducersAccuracy(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return 1.0 - getRecall(classIdx, countedEntity);}
		virtual double getCommission(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return 1.0 - getPrecision(classIdx, countedEntity);}
		virtual double getUsersAccuracy(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return 1.0 - getPrecision(classIdx, countedEntity);}
		virtual double getFallout(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getSpecificity(classIdx, countedEntity);}
		virtual double getFalseNegativeRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getMissRate(classIdx, countedEntity);}
		virtual double getFalsePositiveRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return getFalseAlarmRate(classIdx, countedEntity);}
		virtual double getFidelity(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return 1.0 - getErrorRate(classIdx, countedEntity);}
		virtual long int** getScatterMatrix(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {return this->scatterMatrix[countedEntity][classIdx];}

    //====================================================================================================================
    //	Convert class-tags (sclib::atSpeakerBoundary etc.) to indices into the scatter-matrix/result-vectors or vice 
		//  versa; SVLIB_Fail is returned if the mapping can't be established
    //====================================================================================================================
		long int class2idx(unsigned long int classTag); 
		long int idx2class(unsigned long int classIdx, char* className = NULL, bool shortName = false);

    //====================================================================================================================
    //	Returns an or-concatenated list of the audio-types this class is responsible for scoring 
		//  (e.g. sclib::atSpeech|sclib::atNoise) or sclib::noType if there is no such type (e.g. in case of speaker id)
		//  This is useful e.g. to now which types can be fed into class2idx()
    //====================================================================================================================
		long int responsibility(void) {return sclib::atSpeakerBoundary|sclib::atNoiseBoundary;}
};

#endif
