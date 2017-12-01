/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the audio-type    */
/*        classification process																					*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 27.09.2006																								*/
/**************************************************************************/

#ifndef __SC_Score_AudioTypeClassification_H__
#define __SC_Score_AudioTypeClassification_H__

#include "SC_Score.h"

class SCLIB_API SC_Score_AudioTypeClassification : public SC_Score {
	private :

  protected :
    
		//members to hold scores of the classification process
		unsigned long int ca; //numer or count of classes
		long int **scatterMatrix[4];
		double averageOmission[4];
		double averageCommission[4];
		double missclassificationRate[4];
		double kappaStatistic[4];
		double *recall[4];
		double *precision[4];
		double *missRate[4];
		double *falseAlarmRate[4];
		double *errorRate[4];
		double *specificity[4];
		double *accuracy[4];

		//special members to hold extra speech/non-speech results
		long int **sn_scatterMatrix[4];
		double sn_averageOmission[4];
		double sn_averageCommission[4];
		double sn_missclassificationRate[4];
		double sn_kappaStatistic[4];

		//====================================================================================================================
		//  To make operator<<() kind of virtual...
		//====================================================================================================================
		virtual ostream& output(ostream& OutS); 

	public :

    //====================================================================================================================
    //	The constructor
    //====================================================================================================================
    SC_Score_AudioTypeClassification(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT);

    //====================================================================================================================
    //	The destructor
    //====================================================================================================================
    virtual ~SC_Score_AudioTypeClassification();

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
		virtual double getAverageOmission(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->averageOmission[ce];}
		virtual double getAverageCommission(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->averageCommission[ce];}
		virtual double getMissclassificationRate(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->missclassificationRate[ce];}
		virtual double getKappaStatistic(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->kappaStatistic[ce];}
		virtual double getKHatIndex(unsigned int countedEntity = sclib::ceSample) {return getKappaStatistic(countedEntity);};
		virtual double getSNaverageOmission(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->sn_averageOmission[ce];}
		virtual double getSNaverageCommission(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->sn_averageCommission[ce];}
		virtual double getSNmissclassificationRate(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->sn_missclassificationRate[ce];}
		virtual double getSNkappaStatistic(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->sn_kappaStatistic[ce];}
		virtual double getSNkHatIndex(unsigned int countedEntity = sclib::ceSample) {return getKappaStatistic(countedEntity);};
		virtual double getRecall(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->recall[ce] != NULL) ? this->recall[ce][classIdx] : 0.0;}
		virtual double getPrecision(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->precision[ce] != NULL) ? this->precision[ce][classIdx] : 0.0;}
		virtual double getMissRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->missRate[ce] != NULL) ? this->missRate[ce][classIdx] : 0.0;}
		virtual double getFalseAlarmRate(unsigned long int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->falseAlarmRate[ce] != NULL) ? this->falseAlarmRate[ce][classIdx] : 0.0;}
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
		virtual long int** getScatterMatrix(unsigned int countedEntity = sclib::ceSample) {return this->scatterMatrix[countedEntity];}
		virtual long int** getSNscatterMatrix(unsigned int countedEntity = sclib::ceSample) {return this->sn_scatterMatrix[countedEntity];}

    //====================================================================================================================
    //	Convert class-tags (sclib::atSpeech etc.) to indices into the scatter-matrix/result-vectors or vice versa
		//  SVLIB_Fail is returned if the mapping can't be established
    //====================================================================================================================
		long int class2idx(unsigned long int classTag);
		long int idx2class(unsigned long int classIdx, char* className = NULL, bool shortName = false);

    //====================================================================================================================
    //	Returns an or-concatenated list of the audio-types this class is responsible for scoring 
		//  (e.g. sclib::atSpeech|sclib::atNoise) or sclib::noType if there is no such type (e.g. in case of speaker id)
		//  This is useful e.g. to now which types can be fed into class2idx()
    //====================================================================================================================
		long int responsibility(void) {return sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atMusic|sclib::atBackground|sclib::atSilence|sclib::atBreath|sclib::atAction|sclib::atUndefined|sclib::atSpeech|sclib::atNoise;}
};

#endif
