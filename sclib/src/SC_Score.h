/**************************************************************************/
/*    Responsibility:																											*/
/*		  - base class for computing scores (missmatch between groundtruth  */
/*        and algorithmic results) according to different measures and    */
/*        different tasks                                                 */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Score_H__
#define __SC_Score_H__

#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_GroundTruth.h"
#include "SC_TweakableParameters.h"

class SCLIB_API SC_Score {
	private :

  protected :

		unsigned long int start; //start-sample of score-claculation
		unsigned long int end; //end-sample of score-calculation
		char purpose[sclib::bufferSize]; //a string containing a short statement what a derived class does score ("Audio-type classification assessment" etc.) for report printing

		SC_TweakableParameters *pTweak;
		SC_GroundTruth *pGT;

    //====================================================================================================================
		//  To make operator<<() kind of virtual...
		//====================================================================================================================
		virtual ostream& output(ostream& OutS); 

    //====================================================================================================================
		//  Regarding the 2 extra rows in the scatter-matrices, population & possible set:
		//    Population is used for: 
		//      - accuracy [detection measure]
		//    Possible is used for:
		//      - missclassification rate [s/us classification score]
		//      - kappa statistic [s/us classification score]
		//      - overall recall [unsupervised classification score]
		//      - average class purity [unsupervised classification score]
		//      - average cluster purity [unsupervised classification score]
		//      - purity [unsupervised classification score]
		//====================================================================================================================

    //====================================================================================================================
    //	This is rather simple: Based on the number of overall indetifyable items and the number of correct 
		//  identifications, the id-accuracy (or recognition rate) is returned as the percentage of correctItems; the return
		//  value may indicvate error by returning SVLIB_Fail, otherwise SVLIB_Ok
    //====================================================================================================================
		int calcIdentificationScores(long int overallItems, long int correctItems, double &accuracy);

    //====================================================================================================================
    //	Receives a scatter matrix (4x2) of the following format:
    //
    //               GroundTruth
    //  Hypothesized [ TP | FP ] (number of true positives, false positives of the current CE)
    //               [ FN | TN ] (number of false negatives, true negatives of the current CE)
    //               [PPop|NPop] (number of positives/negatives according to ground-truth in the population)
    //               [PPos|NPos] (number of positives/negatives according to ground-truth in the possible set)
    //
    //  the "population" is the count of overall available CEs in the area under consideration
    //  the "possible set" is the count of CEs actually subject to the detection task (e.g. for change detection, only 
    //  segment boundaries are actually searched for changes, where the overall population would also include the parts in 
    //  the segments)
    //
    //  Returns recall (=sensitivity =true positive rate), precision, specificity (=fallout), miss rate (=false negative 
		//  rate=MDR), false alarm rate (=FAR=false positive rate), error rate (=1-fidelity) and accuracy as parameters and an 
		//  integer indicating an error (if <> SVLIB_OK)
    //====================================================================================================================
    int calcDetectionScores(long int** scatterMatrix, double& recall, double& precision, double& missRate, double& falseAlarmRate, double& errorRate, double& specificity, double& accuracy);

		//====================================================================================================================
		//	Receives a scatter matrix ((cu+2) x ca) of the following format:
		//
		//               GroundTruth (classes: 1..ca)
		//  Hypothesized [ x_1_1  | x_1_2 | ... | x_1_ca  ] (the correct classified CEs appear on the main diagonal of the 
		//  (classes:    [ x_2_1  | x_2_2 |     | ...     ]  matrix, i.e. x_i_i for i=1..ca; the rest in this upper part shows
		//   1..ca)      [ ...    | ...   |     |         ]  the spread of errors)
		//               [ x_ca_1 |       |     | x_ca_ca ] 
		//               [Pop_1   | ...   |     | Pop_ca  ] (overall number of CEs belonging to this class, maybe not all where clustered)
		//               [Pos_1   | ...   |     | Pos_ca  ] (overall number of CEs belonging to this class in the possible set, maybe more than where clustered)
		//
		//  the "population" is the count of overall available CEs in the area under consideration
		//  the "possible" is the count of CEs actually subject to the classification task (e.g. for change detection, only 
		//  segment boundaries are actually searched for changes, where the overall population would also include the parts in
		//  the segments)
		//
		//  The correspondence betweeen correct clusters and classes is provided by the arrangement of the rows and cols in
		//  the scatter matrix: It is such that cluster/row i belongs to class/col i; 
		//
		//  Returned are average omission & commission, the missclassification rate and the kappa statistic (=k-hat index) to 
		//  characterize the overall result; for ech class separately, recall (=1-omission =1-producer's accuracy) & precision 
		//  (=1-commission =1-user's accuracy), miss & FA & error rate, specificity and accuracy (for further details, see 
		//  calcDetectionScores() comment) are provided
		//
		//  If additionalVectorCols>0, the per-class result vectors are made such an amount bigger (may be needed by the 
		//  calling function)
		//
		//	The 2nd available version doesn't do the per-class detection scoring
		//
		//  An integer indicating an error (if <> SVLIB_OK) is returned
		//====================================================================================================================
    int calcSupervisedClassificationScores(long int** scatterMatrix, unsigned long int ca, double& averageOmission, double& averageCommission, double& missclassificationRate, double& kappaStatistic, double*& recall, double*& precision, double*& missRate, double*& falseAlarmRate, double*& errorRate, double*& specificity, double*& accuracy, unsigned int additionalVectorFields = 0);
    int calcSupervisedClassificationScores(long int** scatterMatrix, unsigned long int ca, double& averageOmission, double& averageCommission, double& missclassificationRate, double& kappaStatistic);

    //====================================================================================================================
    //	Receives a scatter matrix ((cu+2) x ca) of the following format:
    //
    //               GroundTruth (classes: 1..ca)
    //  Hypothesized [ x_1_1  | x_1_2 | ... | x_1_ca  ] (the correct clustered CEs appear on the main diagonal of the 
    //  (clusters:   [ x_2_1  | x_2_2 |     | ...     ]  matrix, i.e. x_i_i for i=1..ca; the rest in this upper part shows
    //   1..cu)      [ ...    | ...   |     |         ]  the spread of errors)
    //               [ x_cu_1 |       |     | x_cu_ca ] (cu >= ca)
    //               [Pop_1   | ...   |     | Pop_ca  ] (overall number of CEs belonging to this class, maybe not all where clustered)
    //               [Pos_1   | ...   |     | Pos_ca  ] (overall number of CEs belonging to this class in the possible set, maybe more than where clustered)
    //
    //  the "population" is the count of overall available CEs in the area under consideration
    //  the "possible" is the count of CEs actually subject to the classification task (e.g. for change detection, only 
    //  segment boundaries are actually searched for changes, where the overall population would also include the parts in
    //  the segments)
		//
		//  The correspondence betweeen correct clusters and classes is provided by the arrangement of the rows and cols in
    //  the scatter matrix: It is such that cluster/row i belongs to class/col i; the problem for the just fitting 
    //  clusters (redundant clusters, where CEs of one class are in the majority and therefore this class assigns the 
    //  cluster's label, but there is a bigger one which is labeled as correct) is solved by ordering them this way:
    //  after row ca comes the the biggest cluster fitting (but not being correct) to class 1, followed by the 2nd biggest
    //  etc. after all fitting clusters of class 1 comes the biggest cluster fitting to class 2 and so on.
    //
		//  In addition, some information concerning the order of the classification process, if it was hierarchical, is 
		//  expected to compute average precision as a measure of quality of the cluster merging process in an information 
		//  retrieval sense. This information is included in the array mergeHistory, where each cell i holds an array j=1..J 
		//  of 1/0 content telling if j'th merge was good with respect to cluster i; because each merge-histroy can have a 
		//  different length J, each array is terminated after an entry with the value of "-1". a merge is considered as 
		//  "good/1", if the  merged clusters both where fitting to the final id. After the "-1" entry there follows one 
		//  last element in the array containing the (gt-based) total number of correct segments for this cluster, to include 
		//  recall into the average-precision measure
		//  
		//                       1.. j ..J
		//  1..    cluster_1  -> 1|0|0|1|1|-1|3
		//    i    cluster_2  -> 1|1|-1|2
		//    ..cu cluster_cu -> 1|0|1|1|1|1|1|0|-1|7
		//
		//  Returned are overall recall & precision, the missclassification rate (=error rate), average class & cluster 
		//  purity, average average precision, purity, the rand index and the BBN metric as measures of the overall goodness.
    //  To assess goodness in every single class, the following measures are provided per class (sometimes also per 
		//  cluster), in the order the classes/clusters appear in the scatter matrix, respectively: average precision, class &
    //  cluster purity, all values returned by calcDetectionScores() (as each class on its own can be viewed as a
    //  detection problem).
    //
    //  An integer indicating an error (if <> SVLIB_OK) is returned
    //====================================================================================================================
    int calcUnsupervisedClassificationScores(long int** scatterMatrix, unsigned long int cu, unsigned long int ca, long int** mergeHistory, long int maxMergeHistoryLength, double& overallRecall, double& overallPrecision, double& missclassificationRate, double& averageClusterPurity, double& averageClassPurity, double& average2precision, double& purity, double& randIndex, double& BBNmetric, double*& averagePrecision, double*& classPurity, double*& clusterPurity, double*& recall, double*& precision, double*& missRate, double*& falseAlarmRate, double*& errorRate, double*& specificity, double*& accuracy);

    //====================================================================================================================
    //	Flips the elements of a 4x2 detection scatter matrix such that the other class' values are now the positives
    //====================================================================================================================
		void flipDetectionScatterMatrix(long int** scatterMatrix);

		//====================================================================================================================
		//	ranking1 and ranking2 represent one classification experiment per line, with the score for different clases in 
		//  the columns; both experiments are meant to are identical, except that the models/classifiers used to create the
		//  scores shall differ. The concordance score gives now a value in the interval [0..1] of how much both 
		//  classifiers/models agree in their ranking (classification) on all possible pairs of classes (columns) for each
		//  row. 1 means perfect concordance (match), 0 means total discordance
		//====================================================================================================================
		double getConcordance(double** ranking1, double** ranking2, int len, int dim);

  public :

    //====================================================================================================================
    //	The constructor
    //====================================================================================================================
    SC_Score(SC_TweakableParameters *pTweak, SC_GroundTruth *pGT);

    //====================================================================================================================
    //	The destructor; doesn't destruct the linked list in which this object might be a member
    //====================================================================================================================
    virtual ~SC_Score();

    //====================================================================================================================
    //	Fills the internal score-variables by computing their values according to the frameList in pGT, so that the
    //  get*()-Functions return reasonable values (before calling calcScores(), they return all 0)
		//  "start" and "end" refer to sample-numbers so that the area of the frameList for which scores shall be computed can
		//  be specified; this way, scores can be calculated only for parts of the video/corpus, e.g. for a scene.
		//  In algorithmicUncertaintyDiameter a value [in samples!] can be given which describes the precision with which the
		//  specific algorithm responsible for the results (and only the algorithm, not the gt...) can predict the place of 
		//  event-on- and -offsets
    //====================================================================================================================
		virtual void calcScores(unsigned long int start = 0, unsigned long int end = 0, unsigned long int algorithmicUncertaintyDiameter = 0) = 0;

		//====================================================================================================================
		//	Copies a string containing a description of the given CE or "" if the CE is invalid into the elsewhere alocated 
		//  buffer
		//====================================================================================================================
		virtual void ce2text(unsigned int countedEntity, char *buffer);

		//====================================================================================================================
		//	Presents the results, also care for the output of the complete framelist when wished in the tweakable parameters
		//  this method works bestthe following way: a linked list of score-objects in the order of processing should be 
		//  created (e.g. ATC->V/Uv->CD->Clustering->ID for speaker diarization) and then the printReport() metod of the first
		//  list element should be invoked with parameters such sthat the introduxctory header, the complete framelist and the 
		//  other scores are printed. This will pwoduce a complete estimony of the overall analyzation run and all scoreable
		//  algorithms in a suitable form. This procedure can on wish be repeated for different parts of the analyzed data: 
		//  For the whole thing at once (start=0, end=pGT->getAudioSampleCount()-1) and for each scene separately...
		//====================================================================================================================
    virtual void printReport(const char* fileName, bool introductoryHeader = true, unsigned long int elapsedSeconds = 0, bool allowFrameListOut = true, bool ignoreNextScore = false);
		friend ostream&	operator<<(ostream& os, SC_Score& pScore);

    //====================================================================================================================
    //	To build a linked list for complete report printing
    //====================================================================================================================
		SC_Score *Next;

    //====================================================================================================================
    //	Returns an or-concatenated list of the audio-types this class is responsible for scoring 
		//  (e.g. sclib::atSpeech|sclib::atNoise) or sclib::noType if there is no such type (e.g. in case of speaker id)
		//  This is useful e.g. to now which types can be fed into class2idx()
    //====================================================================================================================
		virtual long int responsibility(void) = 0;
};

#endif
