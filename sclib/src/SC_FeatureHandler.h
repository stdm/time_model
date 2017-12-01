/**************************************************************************/
/*    Responsibility:																											*/
/*		  - provides the possbility to extract all implemented features     */
/*        in one framework                                                */
/*      - helps to save/load features to a file                           */
/*      - helps aggregating features per sub-clip or to a combined vector */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.02.2006																								*/
/**************************************************************************/

#ifndef __SC_FeatureHandler_H__
#define __SC_FeatureHandler_H__

#include <map>
#include <string>
#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Corpus.h"
#include "SC_Signal.h"
#include "SC_MatrixFunctions.h"
#include "SC_SVM.h"
#include <SV_Data.h>
#include <SV_Signal.h>

class SCLIB_API SC_FeatureHandler {
	private :

  protected :

		typedef struct {
       char* fileName;
			 int featureLen;
       int sampleStep;
       short sampleSize;
       short sampleType;
       int featureDim;
       int sampleCount;
		} SC_HTK_HeaderType;

    SC_TweakableParameters *pTweak;
    unsigned long int featureCount;
		bool verbose;

		//====================================================================================================================
		//	Auxiliary method used by extractFeatures() that saves original pTweak entries, overwrites them with special 
		//  parameters, extracts the corresponding feature and restores the original parameters
		//====================================================================================================================
		template<class T> void extractWithSpecialParameters(SV_Data** &pFeatures, unsigned long int featureType, T *pSpecialParameters, T *pOriginalParameters, T *pEffectiveParameters, SC_Corpus *pCorpus, SC_Signal *pSignal, unsigned long int segmentStart, unsigned long int segmentEnd, bool extractOnlySpecialFeatures = true) {
			bool hasAllocated = false;

			if (*(T*)(pSpecialParameters)!=*pEffectiveParameters || extractOnlySpecialFeatures==true) { //only extract a second version of the feature-set (there is already one for standard parameters!) if the parameters really differ
				if (pOriginalParameters == NULL) { //save the original parameters, if not already done
					pOriginalParameters = new T();
					*(T*)(pOriginalParameters) = *pEffectiveParameters;
					hasAllocated = true;
				}

				*pEffectiveParameters = *(T*)(pSpecialParameters); //copy special parameters to the location where extractFeatures() expects the parameters used for feature extraction
				extractFeatures(pFeatures, pCorpus, pSignal, segmentStart, segmentEnd, featureType); //only extract the current feature
				*pEffectiveParameters = *(T*)(pOriginalParameters); //restore original parameters

				if (hasAllocated == true) {
					MFree_0D(pOriginalParameters);
				}
			}

			return;			
		}

		//====================================================================================================================
		//	This method is called by extractFeatured() to extract a single feature with a given extractor and care for proper
		//  debug output and progress report
		//====================================================================================================================
		SV_Data* extractFeature(SC_Corpus *pCorpus, SC_Signal *pSignal, SV_Feature* pExtractor, SC_TweakableParameters::SC_FeaturePar *pParameters, unsigned long int segmentStart, unsigned long int segmentEnd);

    //====================================================================================================================
		//  this method receives an already allocated pFeatures array of proper (this->featureCount) size and fills it with 
		//  the respective selected feature-sets according to the parameters in pTweak->feature*. if an array-entry already 
		//  includes a feature-set, the new one is added to the end of the linked list it represents.
    //====================================================================================================================
    void extractFeatures(SV_Data** &pFeatures, SC_Corpus *pCorpus, SC_Signal *pSignal, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int featureTypes);

		//====================================================================================================================
		//	this is called by branchBoundFeatureSelection() to do the recursive branching with backtracking; the scatter 
		//  matrixes are used to compute the performance measer inside, see also branchBoundFeatureSelection(). the colKey
		//  is an array with a cell per column of the original feature set, having a 'x' character at those locations where
		//  the corresponding column belongs to the current sub-featureset. the result is the best performing colKey.
		//====================================================================================================================
		void branch(double **withinScatter, double **betweenScatter, char *colKey, unsigned int originalDim, unsigned int level, unsigned int maxLevel, std::map<std::string, int> &visited, double &bestPerformance, char* &bestColKey, SC_MatrixFunctions *pMatFunc);

  public :
		
    SC_FeatureHandler(SC_TweakableParameters *pTweak, bool verbose = true);
    virtual ~SC_FeatureHandler();

    //====================================================================================================================
    //	Change the parameter-container to extract features with different parameter-settings
    //====================================================================================================================
    void setTweak(SC_TweakableParameters *pTweak) {this->pTweak = pTweak; return;}

    //====================================================================================================================
    //	Get the number of different features extracted by the extractFeatures method (= dimensionality of the returned 
    //  array of SV_Data objects)
    //====================================================================================================================
    unsigned long int getFeatureCount(void) {return this->featureCount;}

    //====================================================================================================================
    //	return a (new) string containing the feature's name for the given feature-nr (a SCLIB_FEATURE_* constant or a 
		//  concatenation of those)
    //====================================================================================================================
    char* getFeatureName(unsigned long int featureNr);
    
    //====================================================================================================================
		//	Method to control feature-extraction procedure: An this->featureCount-dimensional array of feature-vector classes 
		//  is returned, whith the following content:
		//    - pFeatures[0]:  MFCCs or NULL
		//    - pFeatures[1]:  NULL (was: Energy&ZCR)
		//    - pFeatures[2]:  FbE's or NULL
		//    - pFeatures[3]:  Spectrum or NULL
		//    - pFeatures[4]:  Band Periodicity or NULL
		//    - pFeatures[5]:  Brightness&Bandwidth or NULL
		//    - pFeatures[6]:  NFR or NULL
		//    - pFeatures[7]:  Spectrum Flux or NULL
		//    - pFeatures[8]:  Subband Power or NULL
		//    - pFeatures[9]:  ZCR or NULL
		//    - pFeatures[10]: STE or NULL                  //Bing
		//		//- pFeatures[11]: CepstralPeak or NULL         //Bing
		//    //- pFeatures[12]: Wavelet Energy Distribution  //Nan
		//    - pFeatures[13]: LPC or NULL;                 //Jun
		//    - pFeatures[14]: LPCresidual or NULL          //Jun
		//    //- pFeatures[15]: AAP or NULL                  //Basti
		//    //- pFeatures[16]: BFDAC or NULL                //Basti
		//    - pFeatures[17]: SDP or NULL                  //Bing
		//    - pFeatures[18]: Pitch or NULL
		//    - pFeatures[19]: LSPs or NULL
		//    - pFeatures[20]: Formants or NULL
		//    - pFeatures[21]: Samples or NULL
		//  NULL means that this class of feature-vectors wasn't selected by the featureTypes parameter; each entry might 
		//  include a linked list of features of the same class, but extracted using different sets of parameters. in the 
		//  standard case, only the features for the standard parameter sets (pTweak->feature*) are extracted, but additional 
		//  parametersets can be given in the linked list pSpecialParameters.
    //====================================================================================================================
		SV_Data** extractFeatures(SC_Corpus *pCorpus, SC_Signal *pSignal, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int featureTypes, SC_TweakableParameters::SC_FeaturePar *pSpecialParameters = NULL, bool extractOnlySpecialFeatures = false);

    //====================================================================================================================
    //	Save an array of (linked lists of) features as returned by extractFeatures()
    //====================================================================================================================
    bool saveFeatures(const char *fileName, SV_Data **pFeatures, unsigned long int featureTypes = sclib::featureAllFeatures);

    //====================================================================================================================
    //	Save a single feature-set (this function is able to append data to an existing feature-file); linked lists of 
    //  feature-sets are not supported
    //====================================================================================================================
    bool saveFeature(const char *fileName, SV_Data *pFeature);

    //====================================================================================================================
    //	Save a linked list of features to a file 
    //====================================================================================================================
    bool saveFeatureList(const char *fileName, SV_Data *pList);

    //====================================================================================================================
    //	Load all features from the file and return them in one big (merged) SV_Data object instead as a linked list
    //====================================================================================================================
		SV_Data* loadAndMergeFeatures(const char *fileName, long int firstCol = -1, long int lastCol = -1, bool pivot = false);    
    
		//====================================================================================================================
		//	Load all features from the file and return them in one big (merged) SVMproblem object instead as a linked list
		//  that can be used for one-class SVM training (i.e. labels always positive, no weights)
		//====================================================================================================================
		SC_SVMproblem* loadAndMergeSvmFeatures(const char *fileName, long int firstCol = -1, long int lastCol = -1);

		//====================================================================================================================
		//	Merges a linked list of feature-vectors to one big matrix not in main memory but on hard-disk using a temporary
		//  file; the original list is replaced by the merged dataset to save memory (at all time there is only one version
		//  of the data in memory)
		//  if a filename is given explicitly, the data is stored therein and NOT deleted after merging (debug-dir isn't used)
		//====================================================================================================================
		bool mergeOnDisk(SV_Data* &pDataList, const char *fileName = "");

    //====================================================================================================================
    //	Load an array of features as returned by extractFeatures(); if there are more than one dataset of each featureType
    //  in the file, the corresponding entry in the feature-array is a linked list
    //====================================================================================================================
    SV_Data** loadFeatures(const char *fileName, unsigned long int featureTypes = 0);
    
    //====================================================================================================================
		//	Load a single feature, if featureType > 0 (if there are more than 1 entries with given featureType, a linked list 
		//  is returned; if featureType == 0, return a linked list of all featuresets in the file
    //====================================================================================================================
    SV_Data* loadFeature(const char *fileName, unsigned long int featureType = 0);

		//====================================================================================================================
		//	Normalize a set of feature-vectors by rescaling each vector x to the interval [0..1] by using the min/max values 
		//  for each dimension (best obtained from a big set of training-data): x' = (x-min)/(max-min)
		//====================================================================================================================
    SV_Data* normalize(SV_Data *pFeatures, double *min, double *max, bool replace = true);

		//====================================================================================================================
		//	Does normalization as above; the min thereby comes from the first row of pNorm while the max comes from the second 
		//  row
		//====================================================================================================================
    SV_Data* normalize(SV_Data *pFeatures, SV_Data *pNorm, bool replace = true);

		//====================================================================================================================
		//	normalizes the length of each column vector to unity (1)
		//====================================================================================================================
		void normalizeColumns(SV_Data *pFeatures);

		//====================================================================================================================
		//	normalizes the length of each row vector to unity (1)
		//====================================================================================================================
		void normalizeRows(SV_Data *pFeatures);

    //====================================================================================================================
    //	Removes the effect of normalization
    //====================================================================================================================
		SV_Data* unNormalize(SV_Data *pFeatures, SV_Data *pNorm, bool replace = true);

    //====================================================================================================================
    //	Create a normalization matrix as expected by normalize()
    //====================================================================================================================
		SV_Data* createNormalizationMatrix(SV_Data *pFeatures);

		//====================================================================================================================
		//	Standardize a set of feature-vectors by rescaling each vector x by the mean m and standard deviation sd of this 
		//  feature's distribution (best obtained from a big set of training-data): x' = (x-m)/sd
		//====================================================================================================================
		SV_Data* standardize(SV_Data *pFeatures, double *mean, double *sd, bool replace = true);

		//====================================================================================================================
		//	Does standardization as above; the mean thereby comes from the first row of pNorm while the sd comes from the 
		//  second row
		//====================================================================================================================
		SV_Data* standardize(SV_Data *pFeatures, SV_Data *pStd, bool replace = true);

    //====================================================================================================================
    //	Create a standardization matrix as expected by standardize()
    //====================================================================================================================
		SV_Data* createStandardizationMatrix(SV_Data *pFeatures);

    //====================================================================================================================
    //	Concatenate the columns of all non-NULL array-entries of features (normalization should be done before or 
		//  afterwards to ensure that the scales of the features are similar; otherwise, concatenation in not meaningful) in 
		//  the array going from 0 to featureCount-1 (featureCount is necessary is a parameter to be able to handle feature-
		//  sets of other cardinality as returned by extractFeatures(), as is the case in the LZL ATC algorithm)
		//  If selectionMap != 0, only those features are concatenated that have their index set as a sclib::bit-flag in this map
		//  If individualColSelectionMap != NULL, for each previously selected (non-NULL, selected in selectionMap or 
		//  selectionMap == 0) featureset it is lookad at it's index-position: Only these columns are concatenated which's 
		//  bitflag (first features has index 1!) is set to 1
		//  The featureMap provides, indexed by the featureType, start- and end-index of the columns of this featureType
    //====================================================================================================================
    SV_Data* combineFeatureVectors(SV_Data** pFeatures, unsigned long int featureCount, unsigned long int selectionMap = 0, unsigned long int *individualColSelectionMap = NULL);
    SV_Data* combineFeatureVectors(SV_Data** pFeatures, unsigned long int featureCount, std::map<unsigned long int, std::pair<unsigned int, unsigned int> > &featureMap, unsigned long int selectionMap = 0, unsigned long int *individualColSelectionMap = NULL);

    //====================================================================================================================
    //	Can be called prior to combineFeatureVectors to bring all feature sets to a common frame-rate/-size (smallest 
		//  frame-rate is used as the common ground)
    //====================================================================================================================
		void equalizeFrameParameters(SV_Data** pFeatures, unsigned long int featureCount, unsigned long int selectionMap = 0);

    //====================================================================================================================
    //	Load and save features in HTK (Hidden Markov Model Tolkit) format
    //====================================================================================================================
    SV_Data* loadHTKfeatures(const char *fileName);
    int saveHTKfeatures(const char *fileName, SV_Data *pFeatures);

    //====================================================================================================================
    //	Returns the overall number of rows in alle the SV_Data objects in the (possibly) linked list of objects
    //  count only the rows in the first segmentsToCount objects or all if this is =0    
    //====================================================================================================================
    unsigned long int countFeaturesInList(SV_Data *pFeatureList, unsigned long int segmentsToCount = 0);

    //====================================================================================================================
    //	Destruct array of (linked lists of) features as returned by extractFeatures() method
    //====================================================================================================================
    void freeFeatures(SV_Data** &pFeatures);

    //====================================================================================================================
    //	Check if the values of the features are alright (finite)
		//  The variable parameters provide information on where th (first) error occured in case of defective features
    //====================================================================================================================
    bool checkFeatures(SV_Data **pFeatures, unsigned long int &feature, unsigned long int &row, unsigned long int &column);
    bool checkFeatures(SV_Data *pFeatures, unsigned long int &row, unsigned long int &column);

    //====================================================================================================================
    //	Convert to another frame-size/-step by weighted averaging all given frames in the range of a new sized frame
		//  originalSampleRange is the number of samples used to construct the originalFrameSet (used in order to deduce 
		//  correct amount of new frames (there might be one missing in the original set that fits in the new one)).
		//  "ignoreValues" may point to an array of values resembling a valid row of pOriginalFrameSet, which will be excluded
		//  from bulding the linear weighted mean (e.g. if there are values of special meaning, like 0.0 as a value describing 
		//  that there is no pitch in an pitch-feature-set)
    //====================================================================================================================
    SV_Data* convertFrameRate(SV_Data *pOriginalFrameSet, unsigned long int originalSampleRange, unsigned int newFrameSize, unsigned int newFrameStep, double *ignoreValues = NULL);

    //====================================================================================================================
		//	Each Algorithm class gets as input the array pFeatures of all extracted features, one array-entry for a feature.
		//  It expects that if it needs sclib::featureXYZ (a feature ID) for processing, it will find it at position 
		//  pFeatures[sclib::bit(sclib::featureXYZ]. But there is only room for one such feature set, what happens if 
		//  different instances of the same feature where extracted using different parameters (each algorithm can specifiy 
		//  special parameter sets for the features it uses in case they deviate from the standard püarameters)? they are 
		//  deposited as linked lists at the specified location in the array; but, further on, how does the algorithm class 
		//  now find "its" feature-set (with it's special parameters) in the linked list? it doesn't have to find it, it just
		//  always takes the first element of the linked list. THIS method, called just before the call to the algorithm 
		//  itself, cares for putting the right linked-list-elements at the frontal position, given a linked list of the 
		//  special parameters an algorithm wants. so, after calling this method, all non-specific feature-sets (those 
		//  extracted using the standrad parameters) reside in the front position of the linked lists in the pFeatures-array-
		//  entries, except for those feature-sets for which this method got a special parameter-set; in this case, those 
		//  feature-sets corresponding to the given parameter-set are pushed to the front (coherence between feature-sets and 
		//  parameter-sets is established via the "client"-member of a parameter-set and the first character in the 
		//  "Hdr.Signature"-enry of a feature-set, which should both hold a non-zero algorithm ID).
		//====================================================================================================================
		void prepareFeatureSet(SV_Data **pFeatures, SC_TweakableParameters::SC_FeaturePar *pSpecialParameterList = NULL);

		//====================================================================================================================
		//  returns true if all wanted feature sets (and those with special parameters) are existent, false if one is NULL
		//====================================================================================================================
		bool checkFeatureSet(SV_Data **pFeatures, unsigned long int featureTypes, SC_TweakableParameters::SC_FeaturePar *pSpecialParameterList);

    //====================================================================================================================
		//	randomizes the order of rows in the given frameset; in fact, not frames but blocks of size #preserveBlockSize 
		//  frames are re-ordered, i.e. if e.g. preserveBlockSize==5, the first 5 frames will be consecutive also in the 
		//  original, than comes the next block of 5 originally consecutive samples, and so forth
		//====================================================================================================================
		void splice(SV_Data *pFeatures, unsigned int preserveBlockSize = 1);

    //====================================================================================================================
		//	creates and returns a new set of features with #intermediateFrameCount frames between two successive frames of the
		//  original frame-set, which interpolate between them via a sigmoid function of given steepness; blendStep 
		//  corresponds with the preserveBlockSize-parameter of splice(), so that glueing frames could only be inserted all 
		//  #blendStep frames.
		//====================================================================================================================
		SV_Data* blend(SV_Data *pFeatures, unsigned int intermediateFrameCount, double steepness, unsigned int blendStep = 1);

		//====================================================================================================================
		//	same algorithm as above, but different storage container: featurs is an array with pointer to individual array-
		//  rows, one for each row (i.e. its *not* a matrix that can be deallocated by MFree_2D()!). the result is alike.
		//====================================================================================================================
		double** blend(double **features, unsigned int length, unsigned int dim, unsigned int intermediateFrameCount, double steepness, unsigned int blendStep, unsigned int &resultingRows);

    //====================================================================================================================
		//	takes a single feature set (ignoring further linked sets) and returns a new linked list of shorter feature-sets
		//  (#splitCount elements) that have the same time-order as the frames in pFeatureSet and an overlap of 
    //  #overlapInFrames frames between two consecutive sets
		//====================================================================================================================
		SV_Data* splitFeatureSet(SV_Data *pFeatureSet, unsigned int splitCount, unsigned int overlapInFrames);

		//====================================================================================================================
		//	the feature selection algorithm presented in Kotti, Benetos, Kotropoulos, "Computationally Efficient and Robust 
		//  BIC-Based Speaker Segmentation", 2008. Input is a linked list of feature-sets with correspondig boolean labels to
		//  assign each list element to one of two classes (the performance measure uses the relationship if inter- and intra-
		//  class scatter), and the number of feature columns to select. The algorithm uses a depth-first-search branch-and-
		//  bound startegy to find the most effective feature set. Returned is an array of charcters, wehere each bit 
		//  corresponds with the respective column of the original feature set (see mfcc-feature-class for more on this).
		//  see also: "Classification, Parameter Estimation and State Estimation: An Engineering Approach Using MATLAB", 
		//  chapter 6
		//====================================================================================================================
		unsigned char* branchBoundFeatureSelection(SV_Data *pFeatureList, bool *labels, unsigned int resultingColumnCount);

		//====================================================================================================================
		//	return a subset of columns of the given feature set
		//====================================================================================================================
		SV_Data* getSubSet(SV_Data *pFeatures, unsigned int startCol, unsigned int endCol);

		//====================================================================================================================
		//	replaces outliers in the data by the median of the resprective dimension; regard 1.5*IQR (inter quartile range) as
		//  mild outliers, 3.0*IQR as extreme outliers. http://in.answers.yahoo.com/question/index?qid=1006052718450
		//  the percentage of outliers is returned
		//====================================================================================================================
		double removeOutliers(SV_Data *pFeatures, const char mode = sclib::outlierRemoveExtreme); 

		//====================================================================================================================
		//	set all values below minPitch to zero (important to clean pitch sampled from a GMM where zero isn't zero)
		//====================================================================================================================
		void cleanPitch(SV_Data *pPitch, int column, double minPitch);

		//====================================================================================================================
		//	concatenates #framesPerTrajectory frames to a feature trajectory in order to incorporate time-dependencies into a
		//  feature vector; move #trajectoryStep frames forward to begin the next trajectory.
		//  if removeTiming==true, the original data is clustered via #clusteringIterations of kMeans (k=templateCount; if 0, 
		//  66% of #frames) and each vector is the replaced by the nearest template; then, successive identical templates are 
		//  reduced to just one instance before the trajectories are formed. this way it is hoped to remove timing differencs 
		//  in trajectories (e.g. faster or slower articulation of syllables) by retaining the original spectral content
		//====================================================================================================================
		SV_Data* createTrajectories(SV_Data *pFeature, unsigned int framesPerTrajectory, unsigned int trajectoryStep = 1, bool removeTiming = true, unsigned int templateCount = 32, unsigned int clusteringIterations = 10);

		//====================================================================================================================
		//	convert trajectories back to single frames; just what was removed when removing timing misses
		//====================================================================================================================
		SV_Data* unwindTrajectories(SV_Data *pTrajectories, unsigned int framesPerTrajectory, unsigned int trajectoryStep);

		//====================================================================================================================
		//	load data stored in libSVM format
		//====================================================================================================================
		SV_Data* loadLibSVMfeatures(const char *fileName, int* &labels);

		//====================================================================================================================
		//	returns the MD5 checksum of the data matrix (not regarding the header or linked list nature or stuff)
		//====================================================================================================================
		char* getChecksum(SV_Data *pFeature, unsigned long int segmentsToMerge = 0);

		//====================================================================================================================
		//	shrink the given dataset by percent% via setting the row-count to a smaller value (the original data remains the 
		//  same; it can be restored by setting pFeature->Row to its original value)		
		//====================================================================================================================
		int reduceByPercent(SV_Data *pFeature, double percent = 5.0);
};

#endif
