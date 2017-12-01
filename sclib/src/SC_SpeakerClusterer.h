/**************************************************************************/
/*    Responsibility:																											*/
/*      - Performs the speaker clustering process                         */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 03.03.2006																								*/
/**************************************************************************/

#ifndef __SC_SpeakerClusterer_H__
#define __SC_SpeakerClusterer_H__

#include <map>
#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Cluster.h"
#include "SC_SpeakerScore_Clustering.h"
#include "SC_DistanceMeasures.h"

class SCLIB_API SC_SpeakerClusterer {
	
  private:

	protected:
	
    //====================================================================================================================
		// There are some tweakable parameters in the SC_Lib library; they can be centraly managed in this class.
		//====================================================================================================================
    SC_TweakableParameters* pTweak;

		//====================================================================================================================
		//	Claculates the global criterion evaluating the partition according to the termination criterion given in the
    //  tweakable parameters
    //  Parameters: - pPartition is the pointer to the first element of the linked list of clusters forming the partition
    //              - pParent1 and pParent2 are the last merged clusters
    //              - pChild is the result of this mergence
    //              - the numerical result of the algorithm is placed in exactValue, if != NULL
    //              - clusterCount is the current number of clusters in the partition
    //              - speakerCount is the real number of speakers according to the ground-truth
		//====================================================================================================================
    double getClusterGlobalCriterion(SC_Partition *pPartition, SC_Cluster *pParent1, SC_Cluster *pParent2, SC_Cluster *pChild);

		//====================================================================================================================
		//	Decide if it is good to proceed cluster-mergence according to the criterion given in the tweakable parameters
		//====================================================================================================================
    bool getClusterContinuation(double globalCriterion, unsigned int clusterCount = 0, unsigned int speakerCount = 0);

		//====================================================================================================================
		//	Pick and return a pointer to the partition form the linked list of all created partitions, which can be regarded 
    //  as the final result according to the choosen global criterion
		//====================================================================================================================
    SC_Partition* getClusterFinalPartition(SC_GroundTruth *pGT, SC_Partition *pFirstPartition);

		//====================================================================================================================
		//	Returns the distance between the two clusters, based on the precomputed distances between the single segments in 
		//  the clusters, stored in the distance matrix according to the choosen linkage criterion (doesn't apply to 
		//  linkageMerge though, this is computed directly in getClusterDistanceMatrix(); the segmentMap maps segments 
		//  (segmentStart) to rows/cols of the distance matrix
		//====================================================================================================================
    double getClusterDistance(SC_Cluster *pFirst, SC_Cluster *pSecond, double **segmentDistanceMatrix, std::map<unsigned long int, unsigned long int> segmentMap);

		//====================================================================================================================
		//	Compute the distance between to single segment models of 2 clusters
		//====================================================================================================================
    double getModelDistance(SC_Model *pModel1, SC_Model *pBackground1, SV_Data *pFeatures1, unsigned long int segmentCount1, SC_Model *pModel2, SC_Model *pBackground2, SV_Data *pFeatures2, unsigned long int segmentCount2);

		//====================================================================================================================
		//	Compute all pairwise distances between two single segment models, return a distance matrix (upper triangle) and
		//  a map that maps segmentStart samples to rows/cols of the distance-matrix
		// the new distance matrix is returned together with the x/y row/col-indexes of the closest pair of clusters
		//====================================================================================================================
		double** getSegmentDistanceMatrix(SC_Partition *pPartition, std::map<unsigned long int, unsigned long int> &segmentMap);

		//====================================================================================================================
		// (re-)compute the distance-matrix needed by the clustering algorithm based on the merged models in the clusters;
		// the new distance matrix is returned together with the x/y row/col-indexes of the closest pair of clusters and a 
		// mapping
		//====================================================================================================================
		double** getClusterDistanceMatrix(SC_Partition *pPartition, unsigned long int minDistanceIndex[2], double** &segmentBasedDistances, std::map<unsigned long int, unsigned long int> &segmentMap, double **oldDistanceMatrix = NULL);

		//====================================================================================================================
		// Splits the given list of cluster-objetcs: all clusters containing altogether a too short speach portion (as 
		// compared to the tweakable parameter speakerClusterer.speechSegLengthThreshold) are removed from the list and 
		// returned as a newly linked list of "invalid" clusters
		//====================================================================================================================
		SC_Cluster* splitClusterList(SC_GroundTruth *pGT, SC_Cluster *&pAllClusters);

		//====================================================================================================================
		//  This method returns, for a given partition-list (that depicts the succession of a clustering process), the alpha
		//  (BIC penalty factor) that, for the given values of an objective function per list element that must be minimized, 
		//  yields the best result when BIC is used as a termination criterion for clustering. In finalPartitionIdx, the index 
		//  of the partition that would be choosen as final by BIC with the computed alpha is returned.
		//====================================================================================================================
		double tuneBIC(SC_Partition *pPartitionList, double *objectiveFunction, unsigned long int &finalPartitionIdx);

		//====================================================================================================================
		//  This method returns, for a given partition-list (that depicts the succession of a clustering process), the ICR 
		//  threshold that, for the given values of an objective function per list element that must be minimized, 
		//  yields the best result when ICR is used as a termination criterion for clustering. In finalPartitionIdx, the index 
		//  of the partition that would be choosen as final by ICR with the computed threshold is returned.
		//====================================================================================================================
		double tuneICR(SC_Partition *pPartitionList, double *objectiveFunction, unsigned long int &finalPartitionIdx);

		//====================================================================================================================
		//  This method returns, for a given partition-list (that depicts the succession of a clustering process), the WCD 
		//  penalty factor that, for the given values of an objective function per list element that must be minimized, 
		//  yields the best result when WCD is used as a termination criterion for clustering. In finalPartitionIdx, the index 
		//  of the partition that would be choosen as final by WCD with the computed penalty factor is returned.
		//====================================================================================================================
		double tuneWCD(SC_Partition *pPartitionList, double *objectiveFunction, unsigned long int &finalPartitionIdx);

		SC_DistanceMeasures *pDist;
    bool verbose;

	public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_SpeakerClusterer(SC_TweakableParameters* pTweak, bool verbose = true);
		virtual ~SC_SpeakerClusterer();

 		//====================================================================================================================
		// Do the clustering with the initial segmentation found in pPartition; 
		// return a score-object including the final partition and the list of all created partitions
    // This function also alters the framelist in the groundtruth-object
		// Invalid clusters (those who's segments are too short for clustering) are removed from pAllClusters in the 
		// beginning and returned as a separate list in pInvalidClusters; they can be combined with the invalidClusters from
		// the modeling process and used for identification...
		//====================================================================================================================
		SC_SpeakerScore_Clustering* clusterSpeakers(SC_GroundTruth *pGT, SC_Cluster *&pAllClusters, SC_Cluster *&pInvalidClusters, long int progressCount = 0, long int progressMax = 0, bool (*doProceed)(long int &curr, long int max, unsigned int cnt) = NULL);

    //====================================================================================================================
		// Implements the idea from Kwon, Narayanan, "Robust Speaker Identification based on Selective use of Feature 
		// Vectors", 2006: Having all speaker-models and recspective training-sets at hand, rebuild the models based solely on
		// those vectors that have the maximum likelihood (among all given models) with the model they where used to 
		// construct. Thus, models are hoped to be non-overlapping and short utterances might be scored more reliably.
		//====================================================================================================================
		void constructNonOverlappingClusters(SC_Cluster *&pAllClusters);

    //====================================================================================================================
		// Perform speaker clustering manually in an interactive manner, return the final partition
		//====================================================================================================================
		SC_Partition* interactiveClustering(SC_Partition* &pFirstPartition);

    //====================================================================================================================
		// Perform #repetitions runs of random clustering and print the results to the file; "random" here means (thats not
		// so obvious because its a 2-stage process, first choosing cluster-count/-sizes and then assignung clusters to 
		// individual segments, and, when compared to a human, a human will not act totally radnom on the first task): first
		// a number of singel segments for a cluster is choosen according to the given distribution, then segments are 
		// assigned to this cluster randomly. this process is iterated until no unassigned segments remain. the distribution 
		// is expected as an array where the content at index i gives the probability that a cluster of size i is present in
		// the clustering; maxClusterSize, then, is the maximal index into the array (yes, 0-size is included but should be
		// 0); the repetitions exist to find a true random result by performing a kind of monte carlo simulation experiment 
		// and averaging the random scores reached in each radnom run
		//====================================================================================================================
		void randomClustering(const char *fileName, SC_GroundTruth *pGT, SC_Cluster *pClusterList, double *clusterSizeDistribution, unsigned int maxClusterSize, unsigned int repetitions = 1000);
};

#endif

