/**************************************************************************/
/*    All tweakable parameters in all the	SC_* lib                        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 22.03.2004																								*/
/**************************************************************************/

/*
How and when to add a new parameter - a short introduction and manual
----------------------------------------------------------------------

Q: When to add a parameter?
A: Always! Each tiny tweakable parameter in every class shold appear here so that this class
   is the central repository for managing all user-definable values. This includes thresholds,
   small epsilons, alphas, deltas (or whatever), choices of algorithms to use etc.

Q: How to add a parameter?
A: All parameters of a class are grouped together in a new struct named SC_classnamePar, where 
   classname is the name of the class without the prefix "SC_". If you want to add a parameter 
   for a class that has not already some parameters in here, you must also add such a struct and
   then create a new member of that type. The naming-conventions for this members are just 
   "classname" (as above), with the little change that successive upper case letters are not allowed 
   (so convert  them to lower case).
   When there exists the struct and the member, you can add the parametername to the struct; try to
   give it a meaningful name (no matter if it gets long) and also a comment describing exactly what 
   it does.
   Then, you have to alter 3 methods of this class by adding this new parameter there (in all 
   methods here, the paremers are ordered alphabetically according to there struct-name and 
   parameter-name; try to do likewise!):
    1. the constructor: give the  new parameter a meaningful initialization there that manages to 
       serve as a good default parameter (don't expect someone to change it/make up his mind for 
       another setting when he just uses your algorithm!)
    2. operator<<: add the new parameter to both branches (verbose and not verbose) of the output 
       poerator; as above, but the parameter in alphabetical order to the other ones.
    3. setByName(): include a block in the setByName()-method for the new parameter; use the correct
       string-to-parametertype-function to convert the given value to your type; consider checking
       the value if it's content is valid prior to setting the parameter, as is done with the 
       parameter "debug.debugDir"

 Q: How to encode choices?
 A: All choices (which algorithm to use, what modeltype to use, ...) should be defined via constants
    in SC_Aux.h in the namespace sclib. They should begin with a short abbreviation indicating to 
		which group of values the constant belongs (e.g. 'at_' for "audio type"), followed by a meaningful 
		name to describe the choice to make. include a comment after the new constant describing its meaning

 that's it. 
 thilo
*/

#ifndef __SC_TWEAKABLEPARAMETERS_H__
#define __SC_TWEAKABLEPARAMETERS_H__

#include <iostream>
#include <map>
#include <string.h>
#include "SC_Api.h"

struct ltstr { //	a "const char*" comparison operator for maps as used below
	bool operator() (const char *a, const char *b) const {
		return strcmp(a, b) < 0;
	}
};

class SC_ParameterMap { //needed so there are no warnings for template members without dll interface below... (C4251)
	public:
		std::map<const char*, int, ltstr> map;
};

class SCLIB_API SC_TweakableParameters {
  private:

		SC_ParameterMap parameters; //maps parameter-names to integers for setByValue()

	protected:
		
		//for toggeling output-mode
		bool verbose;
    char* originalDebugDir; //here is the original debugDir without prefix stored, to make unsetting the prefix possible
	
	public:

    //Parameters for SC_Classifier_AdaBoost
    typedef struct {
      int weakClassifierType;
			unsigned int maxWeakClassifiers; //upper bound on the number of weak classifiers
    } SC_ClassifierAdaBoostPar;
    SC_TweakableParameters::SC_ClassifierAdaBoostPar classifierAdaBoost;

    //Parameters for SC_Classifier_DecisionStump
    typedef struct {
      bool splitMultiway; //if true, each feature can be split in more than 2 intervals
    } SC_ClassifierDecisionStumpPar;
    SC_TweakableParameters::SC_ClassifierDecisionStumpPar classifierDecisionStump;

    //Parameters for SC_Classifier_ML
    typedef struct {
      unsigned long int modelType;
      unsigned short int minOrder;
      unsigned short int maxOrder;
    } SC_ClassifierMLPar;
    SC_TweakableParameters::SC_ClassifierMLPar classifierMl;

    //Parameters for SC_Classifier_SVM
    typedef struct { //look inside the SC_SVM class for parameter meaning, or refer directly to the svmlib-documentation
	    int svm_type;
	    int kernel_type;
	    double degree;
	    double gamma; //in one-class SVM: th greater the better the fit to the data (danger of overfitting); <<0.1 advised, but for a simple 4-modal 2d-dataset (a cluster in each corner of a square) i needed >20.0 for 4 differnt modes and not just one big circle in the middle
	    double coef0;
	    double nu; //lower bound on the on the fraction of SVs, and an upper bound on the fraction of outliers 
	    double cache_size;
	    double C;
	    double eps;
	    double p;
	    int shrinking;
	    int probability;
      bool doCV; //here, parametersearch utilizing cross validation can be switchen on/off
      int cvFolds; //cross validation will be carried out cvFold-fold
      long cvMaxDatasetSize; //limit cv training set to this size; if the current problem is larger, randomly a subset of size cvMaxDatasetSize is selected
      double cvCoarseCMin; //the following parameters control the aread of the 2-step (coarse/fine) parameter search for C and gamma
      double cvCoarseCMax;
      double cvCoarseCStep;
      double cvCoarseGammaMin;
      double cvCoarseGammaMax;
      double cvCoarseGammaStep;
      double cvFineCStep;
      double cvFineCRadius; //radius around the previously found best working parameter during coarse search
      double cvFineGammaStep;
      double cvFineGammaRadius;
			int oneClassGammaSearchMaxIterations; //max nr of pairwise distances evaluated to find a suitable gamma for one-class svms
			int oneClassGammaSearchRepeats; //if vectors for pairwise distances are choosen randomly, this process will be repeated for a number of times to assure it converges
    } SC_ClassifierSVMPar;
    SC_TweakableParameters::SC_ClassifierSVMPar classifierSvm;

    //Parameters for SC_Cluster
    typedef struct {
      unsigned int mergeMode; //how should two speaker-models (describing the cluster) be merged: retrained or added up? ->  sclib::mergeRetrain or sclib::mergeAddUp
    } SC_ClusterPar;
    SC_TweakableParameters::SC_ClusterPar cluster;

    //Parameters for debugging output
    typedef struct {
		  unsigned long int debugMode; //toggles debug-output for different Algorithms: sclib::dbSpeakerModels, sclib::dbNoiseModels, sclib::dbClustering, SCLIB_DB_MIXMAX, SCLIB_DB_GMM, sclib::dbCompleteResults, sclib::dbWav
		  char* debugDir; //directory to store debug-output and final results to
      bool useDebugPrefix; //if true, the filename of the analyzed file is put as a prefix before all debugging output filenames
    } SC_DebugPar;
    SC_TweakableParameters::SC_DebugPar debug;

    //Parameters for distance measures
    typedef struct {
      double BICpenaltyFactor; //the (de-)emphasize parameter for the model-complexity penalty factor during delta-BIC computation
			unsigned int mergeMode; //how should two speaker-models (describing the cluster) be merged: retrained or added up? -> sclib::mergeRetrain or sclib::mergeAddUp
			double WCDpenaltyFactor; //penalty factor within the WCD global criterion to (de-)emphasize the want of big clusters: the greater this value (>1), the	more big clusters are favoured
			unsigned short int groundDistance; //the ground- or basis distance measure for distances between distributions like the EMD or Beigi distance
			double ICRthreshold; //comparison threshold for ICR; greater threshold produces less clusters; needs to be tuned to actual data: compute ICR measure for several homogeneous clusters, take mean+std of the measurements as threshold
    } SC_DistanceMeasurePar;
    SC_TweakableParameters::SC_DistanceMeasurePar distanceMeasure;

		//for earth mover's distance
		typedef struct {
			unsigned long int maxSigSize;
			unsigned long int maxIterations;
			unsigned short int debugLevel;
		} SC_EmdPar;
		SC_TweakableParameters::SC_EmdPar emd;

    //Parameters for SC_Enhancement
    typedef struct {
			bool doEnhancement; //switches speech enhancement on/off
      unsigned long int minNoiseDuration; //min. number of ms of noise to build a reliable noise model [ms]
      unsigned long int noiseModelUpdateRate; //period in [ms] after which the noise-model should be updated (only in environment-mode and if possible)
			char* speechModelFile; //fileName for the universal clean speech model to load
    } SC_EnhancementPar;
    SC_TweakableParameters::SC_EnhancementPar enhancement;

		class SCLIB_API SC_FeaturePar {
			public:
				int frameSize; //in [ms]
				int frameStep; //in [ms]
				double lowCut; //lowcut-frequency	[Hz]
				double highCut; //highcut-frequency [Hz]
				double preEmphasizeFactor; //pre-emphasis, usually 0.97
				double sampleRate; //!=0.0 only if a special forced sampleRate is needed for a feature
				SC_FeaturePar *Next; //to form linked lists of parameters
				char client; //an id pointing to the algorithm that "ordered" this specific feature set (0 in case of a general feature set)
				unsigned long int featureType; //the feature type id
				SC_FeaturePar();
				int Valid(void) {
					return 1;
				}
				bool operator==(SC_FeaturePar const &second) { //don't use Next and client as comparison elements!
          return this->sampleRate==second.sampleRate && this->frameSize==second.frameSize && this->frameStep==second.frameStep && this->highCut==second.highCut && this->lowCut==second.lowCut && this->preEmphasizeFactor==second.preEmphasizeFactor;
				}
				bool operator!=(SC_FeaturePar const &second) { //don't use Next and client as comparison elements!
          return this->sampleRate!=second.sampleRate || this->frameSize!=second.frameSize || this->frameStep!=second.frameStep || this->highCut!=second.highCut || this->lowCut!=second.lowCut || this->preEmphasizeFactor!=second.preEmphasizeFactor;
				}
		};

    //Parameters for extracting Band Peridiocity feature
		class SCLIB_API SC_FeatureBandPeriodicityPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				SC_FeatureBandPeriodicityPar();
		}; 
    SC_TweakableParameters::SC_FeatureBandPeriodicityPar featureBandPeriodicity;

    //Parameters for extracting Brightness & Bandwidth feature
		class SCLIB_API SC_FeatureBrightnessBandwidthPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				SC_FeatureBrightnessBandwidthPar();
    };
    SC_TweakableParameters::SC_FeatureBrightnessBandwidthPar featureBrightnessBandwidth;

    //Parameters for extracting filter bank energys / MFCC's
		class SCLIB_API SC_FeatureFbEPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int filterBankSize; //count of filters in filter-bank for log filterbank energy features
				unsigned char frequencyScale; //frequency-scale to use for building the filterbank in feature-computation
  			unsigned short int MFCCorder; //feature-order of MFCC's
				//TODO: MFCCcoeffSelection not used in sclib!
  			unsigned char MFCCcoeffSelection[8]; //select specific coefficients for MFCC's; first bit of first element corresponds with first MFCC and so forth
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				bool CMN; //should cepstral mean subtraction be performed?
				bool addDeltas; //should delta features be added?
				bool addDeltaDeltas; //should delta-delta features be added (only works if deltas are added also)
				unsigned char smoothing; //shall spectrum smoothing be performed? sclib::smooth*
				bool dEnergy; //true => first MFCC-Coeff. represents delta feature
				double minFilterBankFrequency; //lower end of the filterbank (in [Hz])
				double maxFilterBankFrequency; //upper end of the filterbank (in [Hz])
				unsigned char resultType; //result type (linear fbe's, log fbe's, mfccs)
				SC_FeatureFbEPar();
    };
    SC_TweakableParameters::SC_FeatureFbEPar featureFbe;

		//Parameters for extracting STE feature  /Bing Shi
		class SCLIB_API SC_FeatureSTEPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				bool useButterworth; //true if the special butterworth filter should be used (only possible for 8kHz signals)
				bool scaleResult; //true if the result should be energy per ms, not energy per frame
				SC_FeatureSTEPar();
    };
    SC_TweakableParameters::SC_FeatureSTEPar featureSte;

    //Parameters in the SC_FeatureHandler class
    typedef struct {
    } SC_FeatureHandlerPar;
    SC_TweakableParameters::SC_FeatureHandlerPar featureHandler;

    //Parameters for extracting MFCC's as in the sv_lib
		class SCLIB_API SC_FeatureMFCCPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int fftSize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				unsigned short int filterBankSize; //count of filters in mel-filter-bank for MFCC's         
  			unsigned short int MFCCorder; //feature-order of MFCC's
				unsigned char coeffSelection[8]; //select specific coefficients for MFCC's; first bit of first element corresponds with first MFCC and so forth
				bool CMN; //cepstral mean normalization
				bool dEnergy; //true => first MFCC-Coeff. represents delta feature
				bool addDeltas; //should delta features be added?
				bool addDeltaDeltas; //should delta-delta features be added (only works if deltas are added also)
				int method; //which algorithm to use? svlib or sclib implementation
				unsigned char sclib_frequencyScale; //frequency-scale to use for building the filterbank in feature-computation
				unsigned char sclib_smoothing; //shall spectrum smoothing be performed? sclib::smooth*
				double sclib_minFilterBankFrequency; //lower end of the filterbank (in [Hz])
				double sclib_maxFilterBankFrequency; //upper end of the filterbank (in [Hz])
				SC_FeatureMFCCPar();
    };
    SC_TweakableParameters::SC_FeatureMFCCPar featureMfcc;

    //Parameters for extracting Noise Frame Ratio feature
		class SCLIB_API SC_FeatureNFRPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				double NFRthreshold; //threshold for the peak in the norm. correlation function between two adjacent magnitude spectra for a frame to be considered as noise (0..1)
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				SC_FeatureNFRPar();
    };
    SC_TweakableParameters::SC_FeatureNFRPar featureNfr;

    //Parameters for extracting spectrum feature
		class SCLIB_API SC_FeatureSpectrumPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				bool logarithmize; //if true, the log-power-spectrum is returned (generally adviceable)
				bool createPhase; //if true, a second SV_Data object is returned (linked to by the Next pointer, identifyable by Hdr.Signature[1]==1) that holds the corresponding phase spectra to the (log) power spectra
				SC_FeatureSpectrumPar();
    };
    SC_TweakableParameters::SC_FeatureSpectrumPar featureSpectrum;

    //Parameters for extracting Spectrum Flux feature
		class SCLIB_API SC_FeatureSpectrumFluxPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				SC_FeatureSpectrumFluxPar();
    };
    SC_TweakableParameters::SC_FeatureSpectrumFluxPar featureSpectrumFlux;

    //Parameters for extracting subband power feature
		class SCLIB_API SC_FeatureSubBandPowerPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short int FFTsize; //needs to be a power of 2; size of FFT
  			unsigned short int window; //which window function to use during signal processing? -> sclib::wnd*
				SC_FeatureSubBandPowerPar();
    };
    SC_TweakableParameters::SC_FeatureSubBandPowerPar featureSubBandPower;

    //by Nan
    //Parameters for extracting ZCR feature
		class SCLIB_API SC_FeatureZCRPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				bool useChebyshev; //true if the special chebyshev filter should be used (only possible for 8kHz signals)
				bool scaleResult; //true if the result should be ZCR/ms and not just #zero-crossings
				SC_FeatureZCRPar();
    };
    SC_TweakableParameters::SC_FeatureZCRPar featureZcr;

		//by Jun
		//Parameters for LPC
		class SCLIB_API SC_FeatureLPCPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short LPCorder; //coef for LPC, usually 16
				unsigned int window; //which window function to use during signal processing? -> sclib::wnd*
				bool computeGain; //by thilo: if true, the each LPC filter's gain is stored in the last (additional) column of the feature vector
				SC_FeatureLPCPar();
		};
		SC_TweakableParameters::SC_FeatureLPCPar featureLpc;

    //by Jun
		//Parameters for LPCresidual
		class SCLIB_API SC_FeatureLPCresidualPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short order; //LPC coef for residual, usually 4
				SC_FeatureLPCresidualPar();
		};
		SC_TweakableParameters::SC_FeatureLPCresidualPar featureLpcResidual;

		//Parameters for LSPs
		class SCLIB_API SC_FeatureLSPPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				unsigned short LPCorder; //order of underlying LPC analysis
				unsigned int window; //which window function to use during signal processing? -> sclib::wnd*
				int method; //choice of algorithm to derive LSPs from LPCs: SPEEX or MELP (SPEEX recommended)
				double delta; //delta for root-searching
				int bisections; //number of bisections during root-search
				double minSeparation; //minimum separation delta for LSPs and MELP method 
				int maxLoops; //maximum number of loops of separation check and MELP method
				SC_FeatureLSPPar();
		};
		SC_TweakableParameters::SC_FeatureLSPPar featureLsp;

    //Parameters for extracting pitch
		class SCLIB_API SC_FeaturePitchPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				int method; //choice of algorithm to derive Pitch: ESPS or SVlib (ESPS recommended)
				bool sqrt; //true if the squareroot of pitch shall be delivered as suggested by Rose [2002]
				float esps_cand_thresh; //crosscorrelation peak threshold: Determines crosscorrelation peak height required for a peak to be considered a pitch-peak candidate. [.01 .99]
				float esps_lag_weight; //weight to pitch interval shortness: Amount of weight given to the shortness of the proposed pitch interval. Higher numbers make high F0 estimates more likely [0 1]
				float esps_freq_weight; //weight to F0 continuity: Strength of F0 continuity.  Higher numbers impose smoother contours. [0 1]
				float esps_trans_cost; //cost of voicing-state transition: Fixed cost of making a voicing-state transition.  Higher numbers discourage state changes. [0 1]
				float esps_trans_amp; //trans_cost modulated by local rate of amp change: Voicing-state transition cost modulated by the local rate of amplitude change.  Higher numbers discourage transitions EXCEPT when the rate of amplitude change is great. [0 100]
				float esps_trans_spec; //trans_cost modulated by local rate of spectral change: Voicing-state transition cost modulated by the local rate of spectral change.  Higher numbers discourage transitions EXCEPT when the rate of spectral change is great. [0 100]
				float esps_voice_bias; //weight for voice state preference: Determines fixed preference for voiced or unvoiced state.  Positive numbers encourage the voiced hypothesis, negative numbers the unvoiced. [-1 1]
				float esps_double_cost; //cost of a rapid one-octave F0 change: The cost of a rapid one-octave (up or down) F0 change.  High numbers discourage any jumps, low numbers permit octave jumps. [0 10]
				float esps_min_f0; //minimum F0 to search for (Hz): Minimum F0 to search for. Note that computational cost grows as 1/min_f0. [10 < min_F0]
				float esps_max_f0; //maximum F0 to search for (Hz): Maximum F0 to search for. [>min_F0 <Fs/2]
				int esps_n_cands; //max number of candidate F0-peaks in any frame: The maximum number of correlation peaks considered as possible F0-peak candidates in any frame.  At most, the top n-cands candidates are considered in each frame. The computational cost grows approximately as n_cands SQUARED! [3 100]
				float esps_wind_dur; //correlation window duration (sec): Size of correlation window.  Computation increases as wind_dur. [10/Fs .1]
				SC_FeaturePitchPar();
    };
    SC_TweakableParameters::SC_FeaturePitchPar featurePitch;

    //Parameters for extracting pitch
		class SCLIB_API SC_FeatureFormantPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				int esps_lpc_ord; //lpc prediction order
				int esps_lpc_type; //use bsa's stabilized covariance if != 0
				int esps_w_type; //window type: 0=rectangular; 1=Hamming; 2=cos**4; 3=Hanning
				double esps_ds_freq; //internally downsample to this samplerate
				double esps_wdur; //frame size for LPC analysis in [s]
				double esps_nom_f1; //nominal F1 frequqncy (to override base settings of hypothesized formant positions)
				double esps_cor_wdur; //window-length for crosscorrelation F0 estimator
				double esps_frame_int; //frame step in [s]
				int esps_nform; //nr. of formants to track (max. MAXFORMANTS)
				SC_FeatureFormantPar();
    };
    SC_TweakableParameters::SC_FeatureFormantPar featureFormant;

    //Parameters for the SC_Feature_SDP class //by bing
		class SCLIB_API SC_FeatureSDPPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				int m;
				int lag;
				int color;
				int n;
				int pictureSize;
				int tau;
				double preEmphasizeFactor;
				SC_FeatureSDPPar();
    };
    SC_TweakableParameters::SC_FeatureSDPPar featureSdp;

    //Parameters for the SC_Feature_Samples
		class SCLIB_API SC_FeatureSamplesPar : public SC_TweakableParameters::SC_FeaturePar {
			public:
				SC_FeatureSamplesPar();
    };
    SC_TweakableParameters::SC_FeatureSamplesPar featureSamples;

    //General parameters
    typedef struct {
		  unsigned int firstScene; //first scene to process
		  unsigned int lastScene; //last scene to process
		  unsigned long int sceneSelection; //bitflags to select specific (maybe not all) scenes between first- and last scene; not considered if =0
	    unsigned long int shortSpeechThreshold; //defines a lower limit to the speech-segment-length: shorter segments won't be considered for model-building (before clustering, during SCD, ...)! zero means: all segments will be considered. in [ms]
      unsigned int pauseSilenceThreshold; //in ms; shorter segments will be regarded as pauses within speech, longer segments will be regarded as speech-separating silence
			char* preClusteringResultsPrefix; //fileName for the pre-clustering results (groundTruth and (short) single clusters); can be saved/loaded if given and existent
			char* featurePrefix; //fileName to save extracted (atc) features to (per scene) so that multiple runs of the same algorithms on the same video dont have to load trhe signal and do feature extraction all the time; can be saved/loaded if given and existent
    } SC_GeneralPar;
    SC_TweakableParameters::SC_GeneralPar general;

    //Parameters in SC_GroundTruth
    typedef struct {
		  unsigned int internalFrameSize; //length of a frame in the internal frameList [ms]
      unsigned long int pseudoSceneLength; //if no scene-boundary-file is given, pseudo-scene-boundarys of approximately this length (in ms) are generated
      short int videoFrameMachineOffset; //on different machines, different video-codecs might be used, so that the video-frames given in the ground-truth may vary by a constant amount from the ones computed on the machine the program runs on; this factor compensates for this
			bool storeProbabilityInformation; //if memory is short, this can be turned off if not needed (e.g. for classifier-training)
    } SC_GroundTruthPar;
    SC_TweakableParameters::SC_GroundTruthPar groundTruth;

    //Parameters for SC_MixtureModel_GMM
    typedef struct {
      double weightLimit; //mixtures with less weights will be killed
		  double varianceLimit; //too little variances may distort the model, so this will serve as a lower threshold according to reynolds/rose
		  double EMthreshold; //threshold to tell the EM-algorithm to stop (level at which we assume no further improvement from further iterations)
      unsigned long int maxEMiterations; //max. number of EM iterations
      unsigned short int kMeansIterations; //how many rounds should the kMeans algorithms iterate for model-initialization?
    } SC_MixtureModelGMMPar;
    SC_TweakableParameters::SC_MixtureModelGMMPar mixtureModelGmm;

    //Parameters for SC_MixtureModel_bGMM
    typedef struct {
      double weightLimit; //mixtures with less weights will be killed
		  double varianceLimit; //too little variances may distort the model, so this will serve as a lower threshold according to reynolds/rose
		  double EMthreshold; //threshold to tell the EM-algorithm to stop (level at which we assume no further improvement from further iterations)
      unsigned long int maxEMiterations; //max. number of EM iterations
			bool fullCovariance; //if true, the full covariance model is used, otherwise diagonal covariance
			bool randomInitialization; //if false, the model is initialized in a deterministic way (for debugging)
    } SC_MixtureModelBGMMPar;
    SC_TweakableParameters::SC_MixtureModelBGMMPar mixtureModelBgmm;

    //Parameters for MIXMAX models
    typedef struct {
      double weightLimit; //mixtures with less weights will be killed
		  double varianceLimit; //too little variances may distort the model, so this will serve as a lower threshold according to reynolds/rose
		  double EMthreshold; //threshold to tell the EM-algorithm to stop (level at which we assume no further improvement from further iterations)
      unsigned long int maxEMiterations; //max. number of EM iterations
      unsigned long int maxERFloops; //max. number of loops too compute approx. cumulative gaussian (error function, erf)
      unsigned char noiseCorruptionType; //additive (unsupported, sclib::nctAdditive) or log-additive (sclib::nctMax) noise in the MIXMAX; also switches between logarithmic and normal arithmetic (postfix _LOG or without postfix)
      bool bgModelCombination; //shall the background model during MIXMAX training/testing be combined with the original background?
      unsigned short int  kMeansIterations; //how many rounds should the kMeans algorithms iterate for model-initialization?
    } SC_MixtureModelMixMaxPar;
    SC_TweakableParameters::SC_MixtureModelMixMaxPar mixtureModelMixMax;

    //Paramerts for UBM models
    typedef struct {
      bool adaptMeans; //should the means of the UBM be adapted or not?
      bool adaptVariances; //should the variances of the UBM be adapted or not?
      bool adaptWeights; //should the weights of the UBM be adapted or not?
      double relevanceFactor; //relevance factor during adaption: it controls how much new data should be observed before a new parameter begins replacing an old one
      unsigned int scoringMethod; //controls the scoring (testing) method: GMM-UBM scoring (fast, returns ormalized likelihood ratio), GMM-UM scoring with top-mixture-cache usage, or standard GMM scoring (omitting the likelihood ration/score normalization): SCLIB_SCORING_*
      unsigned int topCmixtures; //only the C best-scoring mixtures from the UBM are used for scoring the GMM 
      char* ubmFileName; //fileName of the UBM to load/use
      double varianceLimit; //limit variances (sse below/above for explanation)
    } SC_MixtureModelGMMUBMPar;
    SC_TweakableParameters::SC_MixtureModelGMMUBMPar mixtureModelGmmubm;

    //Parameters for Mix2Max-Models
    typedef struct {
      double weightLimit; //mixtures with less weights will be killed
		  double varianceLimit; //too little variances may distort the model, so this will serve as a lower threshold according to reynolds/rose
		  double EMthreshold; //threshold to tell the EM-algorithm to stop (level at which we assume no further improvement from further iterations)
      unsigned long int maxEMiterations; //max. number of EM iterations
      unsigned long int maxERFloops; //max. number of loops too compute approx. cumulative gaussian (error function, erf)
      unsigned short int kMeansIterations; //how many rounds should the kMeans algorithms iterate for model-initialization?
      bool bgModelCombination; //shall the background model during Mix2Max training/testing be combined with the original background?
    } SC_MixtureModelMix2MaxPar;
    SC_TweakableParameters::SC_MixtureModelMix2MaxPar mixtureModelMix2Max;

    //Parameters for Mix2Max-Ex-Models
    typedef struct {
      double weightLimit; //mixtures with less weights will be killed
		  double varianceLimit; //too little variances may distort the model, so this will serve as a lower threshold according to reynolds/rose
		  double EMthreshold; //threshold to tell the EM-algorithm to stop (level at which we assume no further improvement from further iterations)
      unsigned long int maxEMiterations; //max. number of EM iterations
      unsigned long int maxERFloops; //max. number of loops too compute approx. cumulative gaussian (error function, erf)
      bool bgModelCombination; //shall the background model during training/testing be combined with the original background?
      unsigned short int kMeansIterations; //how many rounds should the kMeans algorithms iterate for model-initialization?
    } SC_MixtureModelMix2MaxExPar;
    SC_TweakableParameters::SC_MixtureModelMix2MaxExPar mixtureModelMix2MaxEx;

		//Parameters for the Hidden Markov Model
		typedef struct {
			unsigned int stateCount; //number of states in the HMM
			char *transitionStructure; //the string represents a boolean matrix (columns delimited by ' ', rows by ';', true is '1', false is '0') where a '1' in cell i,j means: hmm-state i can have a transition probability to state j of >0
			unsigned int mixturesPerState; //number of mixture components in the GMM-like densities in each state; equal for all states
			bool useOrthogonalTransform; //if true, features are transformed to be orthogonal prior to training/testing
			bool leftToRight; //if true, a left-to-right model is build (i.e. initial state probability is 1 for state 0 and zero for all others; take care that the transitionStructure also represents a left-to-right linkage!
			unsigned int maxIterations; //upper bound on the number of training iterations for Baum-Welch algorithm
			bool verbose; //if true, the model outputs probabilities during training
		} SC_ModelHmmPar;
		SC_TweakableParameters::SC_ModelHmmPar modelHmm;

		//Parameters for the Vector Quantisation Model
		typedef struct {
			unsigned int codebookSize; //size of the vq codebook
			unsigned int splitMethod; //method to train codebook
			unsigned int maxIterations; //upper bound for iterations during training
		} SC_ModelVqPar;
		SC_TweakableParameters::SC_ModelVqPar modelVq;

		//Parameters for the One-Class SVM Model
		typedef struct {
			bool distanceBasedTesting; //if true, the arc-distance (SC_DistanceMeasures::svmArcDistance()) is used in TestModel() instead of SVM classification result
			bool doParameterSearch; //if true, the SVM searches for a good gamma parameter using the training data; not a good choice if distanceBasedScoring is true!
		} SC_ModelSvmPar;
		SC_TweakableParameters::SC_ModelSvmPar modelSvm;

		//Parameters for the Time Context Model
		typedef struct {
			unsigned int syllableLength; //length in [ms] of a "syllable", i.e. the local time horizon (trajectory) to reflect in the model
			unsigned int trajectoryStep; //frame shift between trajectories in [#frames]
			char subModelType; //model typ of the real model
			bool removeTiming; //if true, timing information is removed in the succession of frames by clustering them into #templateCount templates, and replcing each frame with the nearest template, thereby discarding frames that fall unto the same template than their predecessor; this way, timing is removed by retaining spectral individuality of this feature-set
			unsigned int templateCount; //nr. of templates the data is clustered into when removing timing
			unsigned int clusteringIterations; //nr. of iterations of kMeans in clustering to find timing-removal templates
			bool replaceTrainingData; //if true, the given (list of) training data is replaced with its trajectory-version
			bool checkForTrajectorization; //if true, the given training-/test-data is first checked if already in trajectory-form; if it is (Hdr->signature[2]>0), it is directly used without trajectorizing it again
			char *worldModelFile; //fileName of a saved worldModel (to remove common syllables) of type subModelType
			char *normalizationFile; //fileName for a normalization matrix for the un-trajectorized features
		} SC_ModelTimePar;
		SC_TweakableParameters::SC_ModelTimePar modelTime;

    //Parameters for the SC_ModelHandler class
    typedef struct {
      double SNRthreshold; //in segments with a SNR above this threshold noise is not considered a problem (e.g. they aren't modeled by an GMM-IB)
      unsigned long int foregroundModelType; //which model to use for speaker-modelling? (sclib::mt*)
      unsigned long int backgroundModelType; //which model to use for background-modelling, if used? (sclib::mt*)
      unsigned short int maxNoiseModelOrder; //max. order for noise-GMM's
      unsigned short int maxSpeakerModelOrder; //max. order for speaker-GMM-IB's
		  unsigned long	int msPerGaussian; //for each msPerGaussian milliseconds a model gets one gaussian to represent the feature-vectors
      unsigned long int orderGuessEMsteps; //number of EM stpes for mixture-models during order-guessing (very little... 3-5?)
      unsigned short int orderGuessMode; //how to guess the order of mixture-models?
			char* onlyThisSpeaker; //if given, only segments of a speaker who's name in the groundtruth exactly matching this one are avaluated for pre-clustering model-building
			unsigned long int speakerModelFeature; //A sclibb::feature* constant that defines which (single) feature to use for building the pre-clustering speaker models
			char outlierRemovalMode; //can be set to a sclib::outlierRemove* constant to remove outliers from speech prior to model training in buildSpeakerModels()
    } SC_ModelHandlerPar;
    SC_TweakableParameters::SC_ModelHandlerPar modelHandler;

    //Parameters for Pareto-Models
    typedef struct {
			bool useMarginalDistributions; //model all dimensions separately (=true) or jointly (=false)
		} SC_ModelPareto;
    SC_TweakableParameters::SC_ModelPareto modelPareto;

    //Parameters for quasi-GMM model by lu/zhang
    typedef struct {
      double deltaBIClambda; //amplifying factor in the deltaBIC-formula
      unsigned short int maxMixtures; //maximum mixture-components in qGMM-model
      double percentDifference; //a threshold is needed to decide whether an updated matrix is significantly different from it's old version.
                                //so, a new matrix is computed with values which ar percentDiffernce different from the old one, and their distance is computed
                                //if the distance between the updated and old matrix is greater than the distance between the syntetically altered and the old matrix, the update is regarded as beeing significant.
    } SC_ModelQGMMPar;
    SC_TweakableParameters::SC_ModelQGMMPar modelQgmm;

    //Parameters for SC_Score_SpeakerClustering
    typedef struct {
      double BBNmetricLambda; //how much are a few big clusters favoured against high cluster purity
    } SC_ScorePar;
    SC_TweakableParameters::SC_ScorePar score;

    //Parameters for the LZL audio-segmentation/classification algorithm
    typedef struct {
			char* actionModelFileName; //path to and fileName of the saved, trained ACTION GMM
      char* classifierFileName; //path to and fileName of the saved, trained classifier-tree
      char* featureFileName; //path to and fileName of the saved, aggregated features
      char* normalizationFileName; //path to and fileName of the saved normalization parameters
			double speechSpecificity; //user-specifiyable parameter (videana) to detect more (->1.0) or less (->0.0) speech (instead of music/background/action)
			double musicSpecificity; //user-specifiyable parameter (videana) to detect more (->1.0) or less (->0.0) music (instead of background/action)
      unsigned long int subClipLength; //length (in [ms]) of a subclip to decide its audio-type
			SC_TweakableParameters::SC_FeatureMFCCPar mfccParameters; 
			SC_TweakableParameters::SC_FeatureZCRPar zcrParameters; 
			SC_TweakableParameters::SC_FeatureSpectrumFluxPar sfParameters; 
			SC_TweakableParameters::SC_FeatureSubBandPowerPar sbpParameters; 
			SC_TweakableParameters::SC_FeatureBrightnessBandwidthPar bbParameters; 
			SC_TweakableParameters::SC_FeatureBandPeriodicityPar bpParameters; 
			SC_TweakableParameters::SC_FeatureNFRPar nfrParameters; 
			SC_TweakableParameters::SC_ClassifierSVMPar svmParameters;
    } SC_Segmentation_AudioTypeLZLPar;
    SC_TweakableParameters::SC_Segmentation_AudioTypeLZLPar segmentationAudioTypeLzl;

    //Parameters for the LZL silence detector
    typedef struct {
      double zcrSilenceThreshold; //threshold on the ZCR feature: less indicates silence
      double energySilenceThreshold; //threshold on the energy feature: less indicates silence
			double specificity; //user-specifiyable parameter (videana) to detect more (->1.0) or less (->0.0) silence
			SC_TweakableParameters::SC_FeatureZCRPar zcrParameters;
			SC_TweakableParameters::SC_FeatureSTEPar steParameters;
    } SC_Segmentation_SilenceLZLPar;
    SC_TweakableParameters::SC_Segmentation_SilenceLZLPar segmentationSilenceLzl;

    //Parameters for the change detection algorithm by lu/zhang
    typedef struct {
      unsigned short int lastNdistances; //count of distances used to compute adaptive threshold
      unsigned long int detectorWindowLength; //size for the sliding window in the lu-zhang-window-algorithm, where the changepoint is assumed between this and the next window (in [ms])
      unsigned long int detectorWindowStep; //step for the sliding window in the lu-zhang-window-algorithm (in [ms])
      double adaptiveThresholdAlpha; //amplifying factor for the last N distances
			double bayesianThreshold; //threshold to compare bayesian likelihood ration with
			char* priorsFileName; //fileName (+ absolute or relative path) for saving/loading prior probabilities
			char* time2ChangeModelFileName; //fileName (+ absolute or relative path) for saving/loading the density model for the time2change
    } SC_Segmentation_ChangesLZPar;
    SC_TweakableParameters::SC_Segmentation_ChangesLZPar segmentationChangesLz;


    //Parameters for the change detection algorithm by kotti/benetos/kotropoulos
    typedef struct {
      double r; //resolution at which BIC tests are performed (half mean utterance duration) in [s]
			double lambda; //data-dependant BIC penalty weighting factor
			double tolerance; //range (in [s]) around a groundtruth change point in which a hypothesized change point must lie (+- tolerance/2) in order to be considered a correct detection
			SC_TweakableParameters::SC_FeatureMFCCPar mfccParameters;
    } SC_Segmentation_ChangesKBKar;
    SC_TweakableParameters::SC_Segmentation_ChangesKBKar segmentationChangesKbk;

    //Parameters for SC_SegmentationHandler, the class for calling/handling all segmentation algorithms
    typedef struct {
      unsigned short int audioTypeMode;
      unsigned short int changeDetectorMode; //controls which algorithm to use as a change-detector    
      unsigned short int silenceDetectorMode;
      unsigned short int vUvDetectorMode;
    } SC_SegmentationHandlerPar;
    SC_TweakableParameters::SC_SegmentationHandlerPar segmentationHandler;

    //Parameters in the Adaptive Silence Detector by li/Narayanan/Kuo
    typedef struct {
      unsigned short int energyQuantizationLevel; //count of classes to determine otsu-threshold
    } SC_Segmentation_SilenceLNKPar;
    SC_TweakableParameters::SC_Segmentation_SilenceLNKPar segmentationSilenceLnk;

    //Parameters in the V/Uv detector based on the Adaptive Silence Detector by li/Narayanan/Kuo
    typedef struct {
      unsigned short int energyQuantizationLevel; //count of classes to determine otsu-threshold
    } SC_Segmentation_VUvLNKPar;
    SC_TweakableParameters::SC_Segmentation_VUvLNKPar segmentationVUvLnk;

		typedef struct {
			SC_TweakableParameters::SC_FeaturePitchPar pitchParameters;
		} SC_Segmentation_VUvESPSPar;
		SC_TweakableParameters::SC_Segmentation_VUvESPSPar segmentationVUvEsps;

		//parameters for the signal-handler
		typedef struct {
			int forceSampleRate; //if >0, all loaded signals will be forced (resampled) to have this samplerate (in [Hz])
		} SC_SignalHandlerPar;
		SC_TweakableParameters::SC_SignalHandlerPar signalHandler;

		//parameters for the MPEG decoder class
		typedef struct {
			bool fastSeeking; //true if fast seeking is wished (experimentally at the moment, NOT advised!)
      bool hqResampling; //true if higher quality resampling (than provided by ffmpeg) is wished
			int outputChannelCount; //number of channels to take (possibly after resampling) from an MPEG file
			int outputSampleRate; //sample rate to take (possibly after resampling) from an MPEG file
		} SC_Signal_MPEGPar;
		SC_TweakableParameters::SC_Signal_MPEGPar signalMpeg;

    //Parameters for the SC_SpeakerClusterer class
    typedef struct {
			bool doClustering; //switches speaker-clustering on/off
			bool constructNonOverlappingClusters; //if true, all cluster models are retrained to be maximally apart from each other; after an idea by Kwon and Narayanan
	    unsigned long int speechSegLengthThreshold; //defines a lower limit to the speech-segment-length: shorter segments won't be clustered! zero means: all segments will be considered. in [ms]
  		unsigned short int distanceMeasure; //distance-measure for clustering: sclib::dmCLR or sclib::dmGLR or sclib::dmBeigi at the moment
      unsigned short int globalCriterion; //global criterion for clustering: sclib::gcBIC, sclib::gcWCD or sclib::gcNone at the moment
		  unsigned short int terminationCriterion; //termination-criterion for clustering: sclib::tcGc, sclib::tcTrue, sclib::tcFalse or sclib::tcKnowledge at the moment
			unsigned char linkageMode; //decide if single/complete/average/"merge" linkage is used: sclib::linkage*
			double specificity; //user-specifiyable parameter (videana) to detect more (->1.0) or less (->0.0) speakers (if >1, it is taken as the wanted number of speakers directly)
			char *firstDistanceMatrixPrefix; //fileName for the first distance matrix during speaker clustrering; can be saved/loaded if given and existent
    } SC_SpeakerClustererPar;
    SC_TweakableParameters::SC_SpeakerClustererPar speakerClusterer;

    //Parameters for the SC_SpeakerIdentificator class
    typedef struct {
      unsigned long int dNormSampleCount; //number of samples drawn from the models for monte-carlo estimation of the KL-distances during dNorm computation
			bool doIdentification; //switches subsequent speaker-id on/off
      unsigned short int normalizationMode; //which normalization-mode to use
			bool useUBMs; //if true, all speaker-models are converted to GMM-UBMs before identification to allow better scoring of short segments
    } SC_SpeakerIdentificatorPar;
    SC_TweakableParameters::SC_SpeakerIdentificatorPar speakerIdentification;

    //Parameters for fft, dct, ... in SC_Transform
    typedef struct {
      double taperingLength; //how many percent of the frame shall be affected by the window function (between 0 and 1)?
    } SC_TransformPar;
    SC_TweakableParameters::SC_TransformPar transform;

	  //Parameters for the Clusterer  //by bing
    typedef struct {
			int maxIterations;
			int numCluster;
			int iterationsToIniKMeanList;
    } SC_ClustererPar;
    SC_TweakableParameters::SC_ClustererPar clusterer;

    //Parameters for the resampling class
    typedef struct {
      bool fastConversion;
    } SC_ResamplingPar;
    SC_TweakableParameters::SC_ResamplingPar resampling;

    void setDebugDir(const char* dir);		//brings 'dir' to the desired form
    void setDebugPrefix(const char* prefix); //adds the prefix to the debugDir
    void unsetDebugPrefix(void); //removes the previously setted prefix from debugDir by restoring the value in originalDebugDir

    //====================================================================================================================
    //	The constructor
    //====================================================================================================================
    SC_TweakableParameters(const char *iniFileName = "", bool verbose = true);

    //====================================================================================================================
    //	The destructor
    //====================================================================================================================
  	virtual ~SC_TweakableParameters();

    //====================================================================================================================
    //	Set the parameter speciffied by it's name to the given value, which is converted implicitly to the needed type
    //  The return-value indicates if the parameter existed or the value was valid (false, if one of this failed, else 
    //  true)
    //====================================================================================================================
    bool setByName(const char* parameterName, const char* value);

    //====================================================================================================================
    //	Some functions for printing out the parameters
    //====================================================================================================================
		void toggleVerboseMode(bool verbose) {this->verbose = verbose; return;}
		friend SCLIB_API std::ostream& operator<<(std::ostream& os, SC_TweakableParameters& pTweak);
};

#endif
