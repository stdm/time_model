/**************************************************************************/
/*	This class is used to implement some common tasks that would normally */
/*  by solved by writing a small program that uses the sclib to do some   */
/*  work, such as training a specific algorithm or building specified     */
/*  models. In order to minimize the amount of needed projects of this    */
/*  type, the same functionality can be coded in one single "main"        */
/*  function here that can be feeded with parameters from a uniform       */
/*  external script that can such accomplish several tasks. In the long   */
/*  run this class could exchange programs like trainer, modeller, scivo, */
/*  ...                                                                   */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 04.03.2008																								*/
/**************************************************************************/

#ifndef __SC_MainTasks_H__
#define __SC_MainTasks_H__

#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
class SCLIB_API SC_MainTasks {

	private :
    
		//some parameters that are assumed common to the majority of tasks
    SC_TweakableParameters *pTweak;
		char audioFile[sclib::bufferSize]; //could also be a corpus file name
		char sceneFile[sclib::bufferSize];
		char segmentationFile[sclib::bufferSize];
		double videoFrameRate;

	protected:

	public :

		//====================================================================================================================
		//	(de)-construct this thingy
		//====================================================================================================================
	  SC_MainTasks(SC_TweakableParameters* pTweak = NULL, const char *audioFileName = "", const char *segmentationFileName = "", const char *sceneFileName = "", double videoFrameRate = 25.0);
		virtual ~SC_MainTasks();
	
		//====================================================================================================================
		//	(re)-initialize the class members
		//====================================================================================================================
		void initParameters(SC_TweakableParameters* pTweak, const char *audioFileName, const char *segmentationFileName = "", const char *sceneFileName = "", double videoFrameRate = 25.0);

		//====================================================================================================================
		//	training script for the LZL audio type classifier
		//====================================================================================================================
		bool trainAudioTypeLZL();
		
		//====================================================================================================================
		//	the enhancer script to enhance (cleanse) the speech in the audio file (SCiVo corpus)
		//====================================================================================================================
		bool enhancer();

		//====================================================================================================================
		//	the modeller script to build individual models of given type, feature and audioportion
		//====================================================================================================================
		bool modeller(const char *modelFile, const char *featureFile, const char *normalizationFile, unsigned long int audioTypes, unsigned long int audioTypesNot, unsigned long int featureType, unsigned long int modelType, unsigned short int modelOrder);

		//====================================================================================================================
		//	the combineModels script from modeller
		//====================================================================================================================
		bool combineModels(const char *modelFile1, const char *modelFile2, const char *combinedModelFile, unsigned long int modelType);

		//====================================================================================================================
		//	the statistician script; prints densities of speaker duration distributions and detectability of change points
		//====================================================================================================================
		bool statistician();

		//====================================================================================================================
		//	Bing's testscript for the v/uv stuff
		//====================================================================================================================
		bool vuvDetection();

		//====================================================================================================================
		//	the Wesley script
		//====================================================================================================================
		bool wesley(long int whichExperiment, unsigned long int noiseRobustModelType, char *rawDataDir);

		//====================================================================================================================
		//	the SCiVo script
		//====================================================================================================================
		bool scivo(char *explicitModelFile = NULL, char ***explicitModels = NULL);

		//====================================================================================================================
		//	simple testscript for the ESPS pitch tracker
		//====================================================================================================================
		bool pitchTest();

		//====================================================================================================================
		//	"training script" (prints the 2 thresholds to the console) for the LZL silence detector
		//====================================================================================================================
		bool trainSilenceLZL();

		//====================================================================================================================
		//	test script for a audio type classifier
		//====================================================================================================================
		bool testAudioType();

		//====================================================================================================================
		//	tests various feature-re-synthesis approaches like mfcc2wav etc.
		//====================================================================================================================
		bool synthesisTest();

		//====================================================================================================================
		//	tests the method called by java/videana
		//====================================================================================================================
		bool javaTest(int *shotList = NULL, int shotListLength = 0);

		//====================================================================================================================
		//	testscript for the algorithm in Kotti, Benetos, Kotropoulos, "Computationally Efficient and Robust BIC-Based 
		//  Speaker Segmentation", 2008
		//====================================================================================================================
		bool kottiBICtest(char *bicTrainFile, char*durationModelFile, char *fsClass1File, char *fsClass2File);

		//====================================================================================================================
		//	model building part of the experiments from Reynolds, "Speaker Identification and Verification using Gaussian 
		//  Mixture Speaker Models", 1995; Reynolds, Rose, "Robust Text-Independent Speaker Identification using Gaussian
		//  Mixture Speaker Models", 1995
		//====================================================================================================================
		bool buildReynoldsSpeakerModels(double reduceByPercent = 0.0);

		//====================================================================================================================
		//	model test part of the experiments from Reynolds, "Speaker Identification and Verification using Gaussian 
		//  Mixture Speaker Models", 1995; Reynolds, Rose, "Robust Text-Independent Speaker Identification using Gaussian
		//  Mixture Speaker Models", 1995
		//====================================================================================================================
		bool testReynoldsSpeakerModels(char *modelListFile, double reduceByPercent = 0.0);

		//====================================================================================================================
		//	testcontainer for Jun's HHT-code
		//====================================================================================================================
		//use this method for HHT analysis with 16000 SampleRate
		bool HHTtest(double SampleRate);	
		//bool HHTtest(unsigned long int segStart, unsigned long int segEnd);

		//====================================================================================================================
		//	resynthesizes the speaker-clustering-relevant audio (feature2wav & model2wav) in the given corpus
		//====================================================================================================================
		bool synthesizeCorpus();

		//====================================================================================================================
		//	test task for a webservice: gets a wav-file (name) as input, builds a GMM from all samples inside and returns a
		//  new wav file name of the resynthesized, inverted GMM
		//====================================================================================================================
		bool wav2gmm(char* &resultFileName);

		//====================================================================================================================
		//	all test stuff for the time model
		//====================================================================================================================
		bool timeModelTest();

		//====================================================================================================================
		//	test container to determine length of voice-dependent speech "sounds"
		//====================================================================================================================
		bool splicer(int blockSizeMs = 100, int olaMaxIterations = 150, double olaErrorTarget = 4.0);

		//====================================================================================================================
		//	extracts feature(s) for given file; skips sampleFeed samples at the beginning of the signal in order to evaluate 
		//  what effect small displacements in frame position may have on the feature vectors.
		//  switch off pre-emphasis, resampling and low-/highpass-filtering for the choosen (signle) feature-set (via the 
		//  speaker-model-feature in the modelHandlers tweakable parameters) in order to get comparable results 
		//====================================================================================================================
		bool featureDisplacemetTest(int sampleFeed = 0, bool filePerFrame = true);

		//====================================================================================================================
		//	receives the filename of a saved feature set as "audioFile", loads it and then performs k-means clustering and 
		//  saves the found centroids
		//====================================================================================================================
		bool createTemplates(int k);

		//====================================================================================================================
		//	creates the context vectors (trajectories) for each speaker in the given training set and creates a matrix of 
		//  those vectors (along with a file specifiying the speaker id for each column)
		//====================================================================================================================
		bool createContextVectorMatrix(char *matrixFileName, char *idFileName);

		//====================================================================================================================
		//	finds a sparse representation of the audioFile's context vectors against the trained matrix via Matlab's l1-Magic
		//  and evaluates it
		//====================================================================================================================
		bool testContextVectorMatrix(char *matrixFileName, char *idFileName);

		//====================================================================================================================
		//	produces synthetic data to build a GMM from and visualize it
		//====================================================================================================================
		bool gmmVisualizationTest();

		//====================================================================================================================
		//	as the webvoice methods for time models, but to call via C++ instead of Java
		//====================================================================================================================
		bool wav2timeModel(char *resultPath, int syllableLength = 130, int trajectoryStep = 1, char subModelType = 15, bool removeTiming = true, int templateCount = 0, int clusteringIterations = 100, int frameSize = 20, int frameStep = 10, bool useMFCC = true, bool addPitch = true, int mfccOrder = 20, int lpcOrder = 10, double preemphasis = 0.97, int window = 1, int fftSize = 512, int frequencyScale = 1, int filterbankSize = 24, double minFilterbankFrequency = 0.0, double maxFilterbankFrequency = 7400.0, double olaErrorTarget = 4.0, int olaMaxIterationCount = 150, int pitchFrameSize = 10, int pitchFrameStep = 10, double pitchCandidateThresh = 0.3, double pitchLagWeight = 0.3, double pitchFrequencyWeight = 0.02, double pitchTransitionCost = 0.005, double pitchTransitionAmpModCost = 0.5, double pitchTransitionSpecModCost = 0.5, double pitchVoiceBias = 0.0, double pitchDoublingCost = 0.35, double pitchMinF0 = 10.0, double pitchMaxF0 = 400.0, int pitchCandidateCount = 20, double pitchCandidateWindowSize = 0.0075);

		//====================================================================================================================
		//	audio feature extraction for markus (audio-"sift"-idea)
		//====================================================================================================================
		bool videoFeatureExtractor(char *fileName, unsigned long int featureTypes);

		//====================================================================================================================
		//	creates the context vectors (trajectories) for each speaker in the given TIMIT corpus file, and saves them in a 
		//  new text file under the name fileName
		//  uses the timeModel normalization matrix and speakerModelFeature
		//====================================================================================================================
		bool createTIMITcontextVectors(char *fileName);

		//====================================================================================================================
		//	training the time world model & normalization matrix for ACM MM 2009 experiments, stripped by anything not needed 
		//  to work on the used TIMIT data; 
		//  uses timeModel normalization matrix and worldModel filenames
		//  based on modeller() script
		//====================================================================================================================
		bool trainTimeWorldModelAndNormMatrix(const char *tmpFeaturesFile);
};

#endif
