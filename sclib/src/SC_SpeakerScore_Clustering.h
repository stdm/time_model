/**************************************************************************/
/*    Responsibility:																											*/
/*		  - Computes scores to measure the performance of the created       */
/*        clustering according to the ground truth                  			*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 13.09.2005																								*/
/**************************************************************************/

#ifndef __SC_SpeakerScore_Clustering_H__
#define __SC_SpeakerScore_Clustering_H__

#include "SC_Api.h"
#include "SC_SpeakerScore.h"
#include "SC_Partition.h"
#include "SC_GroundTruth.h"
#include "SC_TweakableParameters.h"
#include <SV_Error.h>

class SCLIB_API SC_SpeakerScore_Clustering : public SC_SpeakerScore {
	private :

  protected :

		unsigned long int cu; //nr. of rows (cu=clusters, i.e. hypothesized speakers) in the scatter matrix (here this can be >=ca)
    long int** scatterMatrix[4]; 
		long int** mergeHistory;
		double overallRecall[4];
		double overallPrecision[4];
		double missclassificationRate[4];
		double averageClusterPurity[4];
		double averageClassPurity[4];
		double average2precision; //exists independent from the current CE
		double purity[4];
		double randIndex[4];
		double BBNmetric[4]; 
		double DER; //exists independant form the current CE (i.e. it is always computed sample-based)
		double* averagePrecision; //exists independent from the current CE
		double* classPurity[4];
		double* clusterPurity[4];
		double* recall[4];
		double* precision[4];
		double* missRate[4];
		double* falseAlarmRate[4];
		double* errorRate[4];
		double* specificity[4];
		double* accuracy[4];
		bool mergeHistoryIsValid; //==false, if the final partition is not a direct successor of previous partitions all present in the partition-list (maybe because it was altered by SC_SpeakerIdentificator), which leads to a meaningless merge-history and average-precision measures.
		bool computeMergeHistory; //merge history gets computed only if this is true;

    SC_Partition *pPartitionList;
    SC_Partition *pFinalPartition;

    //====================================================================================================================
		//  To make operator<<() kind of virtual...
		//====================================================================================================================
		virtual ostream& output(ostream& OutS); 

    //====================================================================================================================
		//	Exchanges for two arbitrary cluster-tags the indexes to which they are mapped
		//====================================================================================================================
		void exchangeClusterMappings(unsigned long int clusterIdx1, unsigned long int clusterIDx2, int countedEntity = sclib::ceSample);

		//====================================================================================================================
		//	Creates (recursively) the merge-history for the given clusterID and stores the findings in the given resultVector
		//====================================================================================================================
		void createMergeHistory(long int clusterID, SC_Partition *pPartitionList, SC_Partition *pCurrentPartition, long int finalOwnerID, long int* resultVector);

		//====================================================================================================================
		//	if a global criterion has been computed during clustering, but hasn't been used as the termination criterion (i.e. 
		//  sclib::tcTrue has been used), we have a final partition with only one cluster but due to the global criterion this
		//  class thinks that "cu" distinct speakers have been found; this method picks and returns the partition from the
		//  list that has "cu" clusters and thus is the wanted "virtual" final partition.
		//====================================================================================================================
		SC_Partition* findVirtualFinalPartition(SC_Partition *pPartitionList, SC_Partition *pFinalPartition, unsigned int virtualSpeakerCount);

		//====================================================================================================================
		//	debug output of the internal mapping of hypo-ids to scatter-matrix rows
		//====================================================================================================================
		void clusterMappingOut(const char *fileName, unsigned int ce);

  public :

    //====================================================================================================================
    //	The constructor
    //  the final partition must be a pointer to an onject within the partition-list in order to get destructed!
    //====================================================================================================================
    SC_SpeakerScore_Clustering(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT, SC_Partition *pPartitionList = NULL, SC_Partition *pFinalPartition = NULL, bool computeMergeHistory = true);

    //====================================================================================================================
    //	The destructor
    //  destructs all the linked partitions, too!!!
    //====================================================================================================================
    virtual ~SC_SpeakerScore_Clustering();

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
		virtual double getOverallRecall(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->overallRecall[ce];}
		virtual double getOverallPrecision(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->overallPrecision[ce];}
		virtual double getAverageClusterPurity(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->averageClusterPurity[ce];}
		virtual double getAverageClassPurity(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->averageClassPurity[ce];}
		virtual double getAverage2Precision(void) {return this->average2precision;}
		virtual double getAverageAveragePrecision(void) {return getAverage2Precision();}
		virtual double getPurity(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->purity[ce];}
		virtual double getRandIndex(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->randIndex[ce];}
		virtual double getBBNmetric(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->BBNmetric[ce];}
		virtual double getAveragePrecision(unsigned int clusterIdx) {return (clusterIdx < this->cu) ? this->averagePrecision[clusterIdx] : 0.0;}
		virtual double getClassPurity(unsigned int classIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (classIdx < this->ca && this->classPurity[ce] != NULL) ? this->classPurity[ce][classIdx] : 0.0;}
		virtual double getClusterPurity(unsigned int clusterIdx, unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return (clusterIdx < this->cu && this->clusterPurity[ce] != NULL) ? this->clusterPurity[ce][clusterIdx] : 0.0;}
		virtual double getDiarizationErrorRate(void) {return this->DER;}
		virtual double getDER(void) {return getDiarizationErrorRate();}
		virtual double getMissclassificationRate(unsigned int countedEntity = sclib::ceSample) {unsigned int ce = (countedEntity < 4) ? countedEntity : 0; return this->missclassificationRate[ce];}
    SC_Partition* getFinalPartition(void) {return this->pFinalPartition;}
    unsigned int getClusterCount(void) {return sclib::getListCount(this->pFinalPartition->getClusters());}

    //====================================================================================================================
    //	Functions to set certain values during the clustering-process to prepare score-calculation
    //====================================================================================================================
    void setFinalPartition(SC_Partition *pPartition) {this->pFinalPartition = pPartition; return;}
		void setPartitionList(SC_Partition *pPartitionList) {this->pPartitionList = pPartitionList;  return;}
		void setComputeMergeHistory(bool computeMergeHistory) {this->computeMergeHistory = computeMergeHistory; return;}

    //====================================================================================================================
    //	Convert cluster-tags (hypo-speaker-id's) to indices into the scatter-matrix/result-vectors or vice versa; 
		//  SVLIB_Fail is returned if the mapping can't be established
    //====================================================================================================================
		long int cluster2idx(unsigned long int clusterTag, int countedEntity = sclib::ceSample);
		long int idx2cluster(unsigned long int clusterIdx, int countedEntity = sclib::ceSample, char* clusterName = NULL, bool shortName = false);

    //====================================================================================================================
    //	Returns an or-concatenated list of the audio-types this class is responsible for scoring 
		//  (e.g. sclib::atSpeech|sclib::atNoise) or sclib::noType if there is no such type (e.g. in case of speaker id)
		//  This is useful e.g. to now which types can be fed into class2idx()
    //====================================================================================================================
		long int responsibility(void) {return sclib::noType;}
};

#endif
