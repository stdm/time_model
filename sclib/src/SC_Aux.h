/**************************************************************************/
/*    Some auxiliary functions needed by SC_*															*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.02.2004																								*/
/**************************************************************************/

#ifndef __SC_AUX_H__
#define __SC_AUX_H__

#include <map>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <limits>
#include <time.h>
#include <math.h>
#include <vector>
#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_MatrixFunctions.h"
#include <SV_Data.h>
#include <SV_DataIO.h>
#include <SV_Error.h>
#include <GN_Rand.h>

namespace sclib { //put all these standard functions and declarations into a common namespace to un-pollute the environment

	//====================================================================================================================
	//  typedefs
	//====================================================================================================================
	//there are some differences between gcc and vc here which can be handled using this simple typedef
	#ifndef _MSC_VER
		typedef std::_Ios_Openmode OpenMode;
	#else
		typedef ios_base::open_mode OpenMode;
		//typedef std::ios_base::openmode sclib::OpenMode;
	#endif

	//====================================================================================================================
	//  constant values
	//====================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
	//  common random seed (and stateful random number generator) to have repeatable pseudo-random results
	//-------------------------------------------------------------------------------------------------------------------
	const int randomSeed = 1;
	GN_Rand* getRandomizer(void);

	//-------------------------------------------------------------------------------------------------------------------
	//  some mathematical constans or pre-calculated common values
	//-------------------------------------------------------------------------------------------------------------------
	const double ln_2 = 0.69314718055994530941723212145818; //natural logarithm of 2
	const unsigned long int pow_2[32] = {1, 2, 4, 8, 16, 
																			 32, 64, 128, 256, 512, 
																			 1024, 2048, 4096, 8192, 16384,
																			 32768, 65536, 131072, 262144, 524288,
																			 1048576, 2097152, 4194304, 8388608, 16777216,
																			 33554432, 67108864, 134217728, 268435456, 536870912,
																			 1073741824, 2147483648UL}; //values of 2^i at index-position i from[0..31]
	const double pi = 3.1415926535897932384626433832795; //=pi
	const double sqrt_pi = 1.7724538509055160272981674833411; //=sqrt(pi)
	const double sqrt_2 = 1.4142135623730950488016887242097; //=sqrt(2)
	const double sqrt_2pi = 2.506628274631000502415765284811; //=sqrt(2*pi)
	const double one_div_sqrt_2pi = 0.39894228040143267793994605993438; //=1/sqrt(2*pi)
	const double two_pi = 6.283185307179586476925286766559; //=2*pi
  const double onehundredeighty_div_pi = 57.295779513082320876798154814105; //=180/pi
  const double sqrt_2_div_2 = 0.70710678118654752440084436210485; //=sqrt(2)/2
	const double pi_div_2 = 1.5707963267948966192313216916398; //=pi/2
  const double log_2pi = 1.8378770664093454835606594728112; //=log(2*pi)

	//-------------------------------------------------------------------------------------------------------------------
	//  Bitflags to represent the content of an audioFrame in the frameList of a groundTruth object; because they are 
	//  used as bitflags, they must be a power of 2!
	//  Several of these audioTypes have to be mutually exclusive, both in the groundtruth and in the hypothesized 
	//  results; the method checkConsistency() checks if all the rules are fullfilled (and it states the rules, too :-)
	//-------------------------------------------------------------------------------------------------------------------
	const long int atNoise = 1; //0: all non-speech-frames of specific minimum energy, e.g. background-noise, music, babble
  
	const long int atPureSpeech = 2; //1: for the audio-segmenter: clean speech (sub-group for speech)
	const long int atNoisySpeech = 4; //2: for the audio-segmenter: noisy speech (sub-group for speech)
	const long int atBackground = 8; //3: for the audio-segmenter: background noise including machine and nature sounds, laughter, babble, etc. (sub-group for noise)
	const long int atMusic = 16; //4: for the audio-segmenter: music (sub-group for noise)
	const long int atAction = 32; //5: for the audio-segmenter: actionsounds such as screams, shots, explosions and other loud transients (sub-group for noise)
	const long int atBreath = 64; //6: for the audio-segmenter: breathing sounds (sub-group for noise)
	const long int atUndefined = 512; //9: 

	const long int atSpeech = 1024; //10: all speech frames of specific minimum energy
	const long int atPause = 2048; //11: short silence regions within speech (speech with too low energy)
	const long int atVoiced = 4096; //12: further details for speech frames: voiced speech
	const long int atUnvoiced = 8192; //13: further details for speech frames: unvoiced speech
	const long int atMaleVoice = 16384; //14: further details for speech frames: male voice
	const long int atFemaleVoice = 32768; //15: further details for speech frames: female voice

	const long int atSilence = 131072; //17: speech and noise regions of low energy are regarded as silence
	const long int atShort = 262144; //18: speech or noise region is too short to analyze (set by some algorithms if they fail)

	const long int atSceneBoundary = 524288; //19: beginning of a new video-scene
	const long int atSpeakerBoundary = 1048576; //20: beginning of a speech-segment of a new speaker
	const long int atNoiseBoundary = 2097152; //21: beginning of a new background noise type
	const long int atShotBoundary = 4194304; //22: beginning of a new video shot
	const long int atArtificialBoundary = 8388608; //23: marks the co-occuring standard-boundary as not being real, but artificially introduced by the readGroundTruth() method for technical reasons (only appears in the groundTruth-column of the frameList)

	const long int atSpeechSegmentStart = 16777216; //24: for the final results: starting point of a evaluated, speaker-specific speech-sgement
	const long int atSpeechSegmentEnd = 33554432; //25: as above, but end-point

	const long int atAllTypes = 2147483647; //31: all types

	//-------------------------------------------------------------------------------------------------------------------
	//  phone codes found in the TIMIT transcription files (phones are acoustically distinguishble units; phonemes in 
	//  contrast are the smallest units that make a difference in meaning of a word; allophones are the different phones 
	//  that can stand for a phoneme)
	//-------------------------------------------------------------------------------------------------------------------
	//stops
	const char phone_b = 1;
	const char phone_d = 2;
	const char phone_g = 3;
	const char phone_p = 4;
	const char phone_t = 5;
	const char phone_k = 6;
	const char phone_dx = 7;
	const char phone_q = 8;

	//closure symbols
	const char phone_bcl = 9;
	const char phone_dcl = 10;
	const char phone_gcl = 11;
	const char phone_pcl = 12;
	const char phone_tck = 13;
	const char phone_kcl = 14;
	const char phone_tcl = 15;

	//affricates
	const char phone_jh = 16;
	const char phone_ch = 17;

	//fricatives
	const char phone_s = 18;
	const char phone_sh = 19;
	const char phone_z = 20;
	const char phone_zh = 21;
	const char phone_f = 22;
	const char phone_th = 23;
	const char phone_v = 24;
	const char phone_dh = 25;

	//nasals
	const char phone_m = 26;
	const char phone_n = 27;
	const char phone_ng = 28;
	const char phone_em = 29;
	const char phone_en = 30;
	const char phone_eng = 31;
	const char phone_nx = 32;

	//semivowels/glides
	const char phone_l = 33;
	const char phone_r = 24;
	const char phone_w = 35;
	const char phone_y = 36;
	const char phone_hh = 37;
	const char phone_hv = 38;
	const char phone_el = 39;

	//vowels
	const char phone_iy = 40;
	const char phone_ih = 41;
	const char phone_eh = 42;
	const char phone_ey = 43;
	const char phone_ae = 44;
	const char phone_aa = 45;
	const char phone_aw = 46;
	const char phone_ay = 47;
	const char phone_ah = 48;
	const char phone_ao = 49;
	const char phone_oy = 50;
	const char phone_ow = 51;
	const char phone_uh = 52;
	const char phone_uw = 53;
	const char phone_ux = 53;
	const char phone_er = 55;
	const char phone_ax = 56;
	const char phone_ix = 57;
	const char phone_axr = 58;
	const char phone_ax_h = 59;

	//other symbols
	const char phone_pau = 60;
	const char phone_epi = 61;
	const char phone_h_sharp = 62;
	const char phone_stress1 = 63;
	const char phone_stress2 = 64;

	//-------------------------------------------------------------------------------------------------------------------
	// constants to define modes (in the gt-functions an elsewhere)
	//-------------------------------------------------------------------------------------------------------------------
	//modes in the gt
	const char modeLabelAdd = 0;
	const char modeLabelRemove = 1;
	const char modeHypothesized = 3;
	const char modeGroundtruth = 4;

	//modes in speaker-id scoring
	const char modeIncludShort = 0; 
	const char modeIgnoreShort = 1;

	//modes in fft-calculation
	const char modeFourierCoefficients = 0; //fft yields fourier-coefficients
	const char modeFourierTransform = 1; //fft yields fourier-transform
	const char modeInverseFourierTransform = 2; //fft yields inverse fourier-transform

	//for LSP extraction
	const char modeSpeex = 1;
	const char modeMELP = 2;

	//for Pitch extraction
	const char modeESPS = 1;
	const char modeSVlib = 2;
	const char modeAS = 3; //by the AS v/uv algorithm

	//for MFCC-extraction
	const char modeSClib = 1;
	const char modeSlaney = 3;
	//const char modeSVlib = 2; //already defined for pitch extraction

	//for change-detection
	const char modeSpeakerChange = 1;
	const char modeAcousticChange = 2;

	//for model-building
	const char modeForeground = 0;
	const char modeBackground = 1;

	//for mfcc-extraction and model-loading
	//const char modeSClib = 1; //already defined for MFCC extraction
	const char modeForeign = 2; //Use foreign code s for several things (mfcc-extraction, model-loading)

	//for vq-model training
	const char modeLBG = 0; //train vq-codebook by LBG algorithm
	const char modeSplit = 1; //train vq-codebook by splitting

	//for outlier removel
	const char outlierRemoveNone = 0;    //no outliers are removed
	const char outlierRemoveExtreme = 1; //only extreme outliers are removed
	const char outlierRemoveMild = 2;    //even mild outliers are removed

	//for cluster-merging
	const unsigned int mergeNone = 0; //don't merge models in child clusters (used when single- or complete linkage is wanted to save some time)
	const unsigned int mergeAddUp = 1; //merge clusters by adding up the mixtures of the two models
	const unsigned int mergeRetrain = 2; //merge clusters by retraining the models on the complete data

	//for clustering linkage selection
	const char linkageMerge = 0; //a newly trained model based on all segments in the cluster is used for distance comparison
	const char linkageSingle = 1; //the two nearest single segment models in the cluster are used for distance computation
	const char linkageComplete = 2; //the two farthest apart single segment models...
	const char linkageAverage = 3; //the average distance from all pairwise single segment distances between the two clusters is used
	const char linkageInteractive = 4; //interactive clustering linkage, to evaluate manual clustering done by humans

	//-------------------------------------------------------------------------------------------------------------------
	//  ground-truth types
	//-------------------------------------------------------------------------------------------------------------------
	const char gtStandard = 0;
	const char gtTIMIT = 1;
	const char gtMPEG7 = 2;
	const char gtSCiVo = 3;
	const char gtWesley = 4;
	const char gtMAC = 5;

	//-------------------------------------------------------------------------------------------------------------------
	//  for signal type determination
	//-------------------------------------------------------------------------------------------------------------------
	const char stGuess = 0; //guess the correct audio-type from the filename's extension (or, if this is ambigous, from the files content)
	const char stWave = 1; //MS RIFF WAVE
	const char stMPEG = 2; //MPEG video
	const char stMP2 = 3; //MP2 audio
	const char stMP3 = 4; //MP3 audio
	const char stNIST = 5; //NIST waveform with SPHERE heder
	const char stOther = 6; //unknown format, will try to open it via ffmpeg (SC_Signal_MPEG)...
	const char stJWave = 7; //WAVE read over a stream from java

	//-------------------------------------------------------------------------------------------------------------------
	//  for classifier-type selection (SVM/ML/...)
	//-------------------------------------------------------------------------------------------------------------------
	const char ctUndefined = -1;
	const char ctML = 0;
	const char ctSVM = 1;
	const char ctDecisionStump = 2;
	const char ctAdaBoost = 3;
	//const char ctNN = 4;
	const char ctTree = 5;

	//-------------------------------------------------------------------------------------------------------------------
	// for (mixture) models
	//-------------------------------------------------------------------------------------------------------------------
	const char mtQGMM = 6; //quasi GMM for LZH SCD
	const char mtGMM_new = 7; //new (own) GMM implementation
	const char mtMIXMAX = 8; //MIXMAX model with integrated background
	const char mtBGMM = 9; //model type for baggenstoss-gmm
	const char mtMIX2MAX_ex = 10; //extended (modified) MIX2MAX model
	const char mtMIX2MAX = 11; //MIX2MAX model
	const char mtPareto = 12; //Pareto Density Estimation in a model
	const char mtGMM_UBM = 13; //GMM with Universal Background Model
	//const char mtPictureCentroids = 14; //SC_Model_PictureCentroids //by Bing
	const char mtSVM = 15; //one-class SVM model
	//const char mtANN = 16; //neural network model
	const char mtHMM = 17; //(continuous density) hidden markov model
	const char mtFullGauss = 18; //full covariance gaussian model
	const char mtVQ = 19; //vector quatization model
	//const char mtGroup = 20; //a meta-model that build separete models per group/class and scores only features with fitting group-models
	const char mtTime = 21; //a (meta?)model that also models the time order of feature vectors
	const char mtMetaGMM = 22; //a GMM that builds a seperate GMM per dimension (to possibly save mixtures for less complicated dimensions)

	//-------------------------------------------------------------------------------------------------------------------
	// error types (in v/uv scoring)
	//-------------------------------------------------------------------------------------------------------------------
	const char etV2Uv = 1;
	const char etUv2V = 2;

	//-------------------------------------------------------------------------------------------------------------------
	// constants for the Counted Entity (ce): count by sample, segment, shot or scene
	//-------------------------------------------------------------------------------------------------------------------
	const char ceSample = 0;
	const char ceSegment = 1;
	const char ceShot = 2;
	const char ceScene = 3;
	
	//-------------------------------------------------------------------------------------------------------------------
	//  distance measures
	//-------------------------------------------------------------------------------------------------------------------
	const unsigned short int dmCLR = 1; //distance-mesaure: Cross Likelihood Ratio
	const unsigned short int dmGLR = 2; //distance-measure: Generalized Likelihood Ratio
	const unsigned short int dmBeigi = 4; //distance-measure: Beigi-Distance between two distributions
	const unsigned short int dmEMD = 8; //distance measure: The Earth Mover's Distancve between two distibutions
	const unsigned short int dmEuclid = 16; //distance-measure: Euclidean ground-distance for sclib::dmEMD or sclib::dmBeigi (combine by e.g. sclib::dmEMD|sclib::dmEuclid)
	const unsigned short int dmMahalanobis = 32; //distance-measure: Mahalanobis ground-distance for sclib::dmEMD or sclib::dmBeigi (combine by e.g. sclib::dmEMD|sclib::dmMahalanobis)
	const unsigned short int dmKullbackLeibler = 64; //distance-measure: Kullback-Leibler ground-distance for sclib::dmEMD or sclib::dmBeigi (combine by e.g. sclib::dmEMD|sclib::dmKullbackLeibler)
	const unsigned short int dmSvmArcDistance = 128; //distance measure: SVM arc distance between one-class SVMs
	const unsigned short int dmBhattacharyya = 256; //distance-measure: Bhattacharyya ground-distance for sclib::dmEMD or sclib::dmBeigi (combine by e.g. sclib::dmEMD|sclib::dmBhattacharyya)

	//-------------------------------------------------------------------------------------------------------------------
	//  termination criterions
	//-------------------------------------------------------------------------------------------------------------------
	const char tcGc = 0; //termination criterion: According to the choosen global criterion
	const char tcTrue = 1; //termination criterion: Always true (one big cluster to see the complete dendrogram), but GC is calculated
	const char tcFalse = 2; //termination criterion: Always false (no clustering, just the results of SCD)
	const char tcKnowledge = 3; //termination criterion: stops after reaching the real speaker count known a priori from the nr. of ground-truth speaker names
	const char tcOptimal = 4; //termination criterion: chooses the partition that minimizes the diarization error rate when compared to ground truth
	const char tcGcOptimal = 5; //termination criterion: chooses a value for the selected global criterion that results in the minimum reachable diarization error rate for this GC

	//-------------------------------------------------------------------------------------------------------------------
	//  global clustering evaluation criterions
	//-------------------------------------------------------------------------------------------------------------------
	const char gcBIC = 0; //gloabl criterion to evaluate a partition/clustering: Bayesian Information Criterion
	const char gcWCD = 1; //global criterion to evaluate a partition/clustering: Within Cluster Dispersion
	const char gcNone = 2; //global criterion to evaluate a partition/clustering: no global criterion computation
	const char gcICR = 3; //global criterion to evaluate a partition/clustering: Information Change Rate by Han/Narayanan

	//-------------------------------------------------------------------------------------------------------------------
	//  global debug mode (must be powers of 2)
	//-------------------------------------------------------------------------------------------------------------------
	const unsigned long int dbNone = 0;	//no debug output
	const unsigned long int dbSpeakerModels = 1; //all final speaker-models: speakerModels.txt
	const unsigned long int dbNoiseModels = 2; //all final noise-models: noiseModels.txt
	const unsigned long int dbClustering = 4; //distance-matrix and value of termination-criterion: dist.txt, globalCriterion.txt
	const unsigned long int dbModelCreation = 8; //the creation-process of a gmm-ib: MIXMAX.txt, MIXMAX_likelihood.txt
	const unsigned long int dbCompleteResults = 32; //show not only speaker-turns-frames, but information for all frames of analyzed audio-file: frameList.txt
	const unsigned long int dbWav = 64; //the frames leading to the speaker- and noise-models to listen to: speech_XXXX.wav, noise_XXXX.wav, mplete_XXXX.wav (where XXXX is the current scene- and speaker-sgement-label
	const unsigned long int dbSegmentStatistics = 128; //statistics for brutto- and netto-speech-segment-length after silence- and unvoiced-removal: segStat.txt
	const unsigned long int dbSNR = 256; //the SNR between the features used to nuild speaker- and corresponding noise modles are printed to SNR.txt
	const unsigned long int dbFeatures = 512; //the exracted features are printed in plain text
	const unsigned long int dbConsistencyCheck = 1024; //are the groundtruth and algorithmic results consistent?
	const unsigned long int dbModelGenerationReport = 2048; //if models are saved via saveModels-Method, an report is generated, too   //nan
	const unsigned long int dbSCD = 4096; //debug-ouput for speaker change detection goes to scd_debug.txt
	const unsigned long int dbClassifierTraining = 8192; //some output during classifier-training
	const unsigned long int dbSpeechTrainData = 16384; //speaker model training data is put out in matlab format in "speechTrainData.txt"

	//-------------------------------------------------------------------------------------------------------------------
	//  for signal processing: which window-function to use?
	//-------------------------------------------------------------------------------------------------------------------
	const char wndNone = 0; //no window
	const char wndRectangle = 0; //just the same as applying no window
	const char wndHamming = 1; //hamming-window
	const char wndHanning = 2; //hanning-window
	const char wndHann = 2; //sometimes the Hanning-window is also called Hann-window
	const char wndVonHann = 2; //..or von-Hann-window, which suites best the name of it's inventor
	const char wndBartlett = 3; //bartlett-window

	//-------------------------------------------------------------------------------------------------------------------
	//  for log-fbe/mfcc-extraction
	//-------------------------------------------------------------------------------------------------------------------
	const char scaleLinear = 0;
	const char scaleMel = 1;
	const char scaleExpoLog = 2;
	const char scaleModMel = 3; 
	const char scaleBark = 4;

	//for result-type in fbe-extraction
	const char resultLinear = 0;
	const char resultLog = 1;
	const char resultCepstrum = 2;

	//for spectrum-smoothing in SC_Feature_FbE
	const char smoothNone = 0; //deactivates smoothing
	const char smoothLight = 1; //cuts highest 10 percent of the frequencys of the powerSpectrum
	const char smoothMiddle = 2; //holds lowest one third of the frequencys of the powerspectrum
	const char smoothHeavy = 3; //holds only the lowest 6 frequencys of the powerSpectrum

	//-------------------------------------------------------------------------------------------------------------------
	//  for model order estimation
	//-------------------------------------------------------------------------------------------------------------------
	const char guessHeuristic = 1; //guess mixture-model order by a heuristic
	const char guessParameterSearch = 2; //guess mixture-model parameters by conducting a parameter search

	//-------------------------------------------------------------------------------------------------------------------
	//  for noise-corruption type definition in the gmm-ib
	//-------------------------------------------------------------------------------------------------------------------
	const char nctAdditive = 0; //the algorithms for additive noise are unsupported/discontinued!
	const char nctAdditiveLog = 1; //there is no log-arithmetic for additive noise
	const char nctMax = 2; //max()-interaction (as an approximation to log-additive noise)
	const char nctMaxLog = 3; //as max(), but the algorithms use logarithmic arithmetic

	//-------------------------------------------------------------------------------------------------------------------
	//  no-X constants
	//-------------------------------------------------------------------------------------------------------------------
	const char noNode = -1; //for the classifierTree class
	const char noPhone = 0; //no phone code
	const char noSpeaker = -1; //default label in the frameList for "no valid spaker-id"
	const char noType = 0; //default label for methods dealing with the framelist indicating that no type is specified
	const char noSegment = -1; //return value of getNextBoundary/getNextSegment indicating that ther is no next segment 
	const char noPeak = -1; //Bing
	const char noPitch = 0; //Bing 

	//-------------------------------------------------------------------------------------------------------------------
	//for SC_Centroid_*
	//-------------------------------------------------------------------------------------------------------------------
	const char centroidGaussian = 1;
	const char centroidPoint = 2;
	const char centroidSignature = 3;

	//-------------------------------------------------------------------------------------------------------------------
	// Some constants used for 'getNextBoundary'
	//-------------------------------------------------------------------------------------------------------------------
	const char searchForward = 0;
	const char searchBackward = 1;
	const char searchMiddle = 2;
	const char searchWithin = 3;

	//-------------------------------------------------------------------------------------------------------------------
	// definitions for labels (for classifiers)
	//-------------------------------------------------------------------------------------------------------------------
	const char labelNegative = -1;
	const char labelPositive = +1;

	//-------------------------------------------------------------------------------------------------------------------
	//  used for timer-measurement conversion
	//-------------------------------------------------------------------------------------------------------------------
	const char alignmentStart = 0;
	const char alignmentEnd = 1;

	//-------------------------------------------------------------------------------------------------------------------
	//  for noise-model update during speech enhancement
	//-------------------------------------------------------------------------------------------------------------------
	const char updateGuess = 0; //the initial value: take the noise-model update method best fitting to the given data and models
	const char updateEnvironment = 1; //take surrounding frames till the minimum noise duration is reached
	const char updateBoundary = 2; //take surrounding frames till the next audiotype-boundarys (if found, otherwise switch to the environment-method)
	const char updateNoUpdate = 3; //don't update the noise models during speech enhancement (useful if a predefined model was supplied)

	//-------------------------------------------------------------------------------------------------------------------
	//  for feature-extraction (SC_FeatureHandler) and -recognition
	//-------------------------------------------------------------------------------------------------------------------
	const long int featureNoFeature = 0;
	const long int featureMFCC = 1; 
	//const long int featureEnergyZCR = 2;
	const long int featureFbE = 4;
	const long int featureSpectrum = 8;
	const long int featureBandPeriodicity = 16;
	const long int featureBrightnessBandwidth = 32;
	const long int featureNFR = 64;
	const long int featureSpectrumFlux = 128;
	const long int featureSubbandPower = 256;
	const long int featureZCR = 512;
	const long int featureSTE = 1024;
	//const long int featureCepstralPeak = 2048;
	//const long int featureWaveletEnergyDistribution = 4096;
	const long int featureLPC = 8192;
	const long int featureLPCresidual = 16384;
	//const long int featureAAP = 32768;
	//const long int featureBFDAC = 65536;
	const long int featureSDP = 131072;
	const long int featurePitch = 262144; //this isn't necessarily extracted by the feature-handler's extractFeatures() method, but can also be a meta-feature that is a by-product of v/uv classification that is later on pasted into the results of the feature-handler
	const long int featureLSP = 524288;
	const long int featureFormant = 1048576;
	const long int featureSamples = 2097152;
	const long int featureAllFeatures = 4194303; //all above constans concatenated by "|"

	//-------------------------------------------------------------------------------------------------------------------
	//  for algorithm selection when there are more than one possible implementations
	//-------------------------------------------------------------------------------------------------------------------
	const char algorithm_nothing = 0; //no algorithm selected -> omit this stage/phase/whatever...
	const char algorithm_cd_LZW = 1; //change detector: Lu/Zhang (sliding-)Window-based
	//const char algorithm_cd_LZS = 2; //change detector: Lu/Zhang Segment-based
	const char algorithm_sd_LNK = 3; //silende detector: Li/Narayanan/Kuo
	const char algorithm_vud_LNK = 4; //voiced/unvoiced detector: Li/Narayanan/Kuo
	const char algorithm_ac_LZL = 5; //acoustic classification: Lu/Zhang/Li
	const char algorithm_sd_LZL = 6; //silence detector: Lu/Zhang/Li
	//const char algorithm_vud_AS = 7; //Bing
	//const char algorithm_vud_JM = 8; //nan
	//const char algorithm_vud_AAP = 9; // Basti
	//const char algorithm_vud_SAM = 10; // Basti
	//const char algorithm_cd_TST = 11; //change detector: Thilo STadelmann
	const char algorithm_vud_ESPS = 12; //voiced/unvoiced detector: using the ESPS pitch estimates
	const char algorithm_cd_KBK = 13; //change detector by Kotti/Benetos/Kotropoulos
	const char algorithm_cd_Std = 14; //standard change detector: puts a change point at the beginning of each new segment

	//-------------------------------------------------------------------------------------------------------------------
	//  for SC_MixtureModel_GMM_UBM.TestModel()
	//-------------------------------------------------------------------------------------------------------------------
	const char scoringGMM_UBM = 1; //use the fast mehtod (using only the C top mixtures) by Reynolds et al.
	const char scoringGMM_UBM_CACHE = 2; //as above, but also try to use the top-mixture-cache (can be turned off only for speed reasons if it is known in advance that there will be no cache hits=?
	const char scoringGMM = 3; //score the GM-UBM as a normal GMM, omitting the UBM and don't doing score-normalization (a log-likelihood is returned, not a ratio of two such things)

	//-------------------------------------------------------------------------------------------------------------------
	//  for SC_SpeakerIdentifiactor
	//-------------------------------------------------------------------------------------------------------------------
	const char normalizationNone = 0; //use no score normalization during speaker identification
	const char normalizationDNorm = 1; //use DNorm score normalization

	//-------------------------------------------------------------------------------------------------------------------
	// for resampling
	//-------------------------------------------------------------------------------------------------------------------
	const char resampleAverage = 0;
	const char resampleFirstChannel = 1; //must have the value "1"!
	const char resampleSecondChannel = 2; //must have the value "2"!
	const char resampleThirdChannel = 3; //must have the value "3"!
	const char resampleFourthChannel = 4; //must have the value "4"!
	const char resampleFifthChannel = 5; //must have the value "5"!

	//-------------------------------------------------------------------------------------------------------------------
	// for change detection
	//-------------------------------------------------------------------------------------------------------------------
	const char refpointStart = 0;
	const char refpointMiddle = 1;

	//-------------------------------------------------------------------------------------------------------------------
	// some upper bounds
	//-------------------------------------------------------------------------------------------------------------------
	const unsigned int maxSpeakers = 768; //the maximum number of speakers in the the ground-truth
	const char maxRectangles = 4; //max. nr of rectangles in a filter in sc_metafeature_rectangle

	//-------------------------------------------------------------------------------------------------------------------
	//  some other constants 
	//-------------------------------------------------------------------------------------------------------------------
	const int bufferSize = 1023; //size of string-buffers to fetch lines from ground-truth files or contain speaker-names etc.
	const char matlabSyntax[] = "MATLAB"; //when used as the separator for matrix-/vector-out methods, matlab-syntax is outputted
	const int encodingPCM = 5; //PCM encoding for RIFF-WAVE signals

	//====================================================================================================================
	//  auxiliary mathematical functions
	//====================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
	//  returns the sum of the integer numbers between (and including) high- and low-end
	//-------------------------------------------------------------------------------------------------------------------
	inline unsigned long int gaussSum(unsigned long int highEnd, unsigned long int lowEnd = 1) {unsigned long int n = highEnd-lowEnd+1; return (n%2 == 0) ? (n/2 * (lowEnd+highEnd)) : ((n-1)/2 * (lowEnd+highEnd-1) + highEnd);}

	//-------------------------------------------------------------------------------------------------------------------
	//  calculates the logarithm to base 2, sclib::log2(x), as log(x)/log(2)
	//-------------------------------------------------------------------------------------------------------------------
	inline double log2(double x) {return log(x) / sclib::ln_2;}
	inline double ld(double x) {return log(x) / sclib::ln_2;} //the same as log2(), just if one looks for the different name

	//-------------------------------------------------------------------------------------------------------------------
	//  calculates the logarithm to base e, sclib::ln(x), which already is log(x) in the standard library
	//-------------------------------------------------------------------------------------------------------------------
	inline double ln(double x) {return log(x);}

	//-------------------------------------------------------------------------------------------------------------------
	//  add two numbers in the log-domain where a addition in the linear domain is needed (involves one log() and one 
	//  exp() function call)
	//-------------------------------------------------------------------------------------------------------------------
	inline double logAdd(double logX, double logY) {return logX + log(1.0 + exp(logY-logX));}

	//-------------------------------------------------------------------------------------------------------------------
	//  subtracts two numbers in the log-domain where a subtraction in the linear domain is needed (involves one log() and 
	//  one exp() function call)
	//-------------------------------------------------------------------------------------------------------------------
	inline double logSub(double logX, double logY) {return logX + log(1.0 - exp(logY - logX));}

	//-------------------------------------------------------------------------------------------------------------------
	//  secure natural logarithm: if the argument is zero or negativ, a very small (logarithmic) number is returned 
	//  instead of a non-finite number
	//-------------------------------------------------------------------------------------------------------------------
	inline double sLog(double x) {return (x > 0.0) ? log(x) : -705.0;}

	//-------------------------------------------------------------------------------------------------------------------
	//  secure exponential function: if the argument is too small/big, zero/4e306 is returned instead of a denaturated 
	//  (DEN) number
	//-------------------------------------------------------------------------------------------------------------------
	inline double sExp(double x) {return (fabs(x) < 706.0) ? exp(x) : ((x >= 706.0) ? 4.0917e306 : 0.0);}

	//-------------------------------------------------------------------------------------------------------------------
	//  add two numbers in the log-domain where a addition in the linear domain is needed (involves one log() and one 
	//  exp() function call)
	//  if the difference between the two arguments is so great that problems may arise during exp()-computation, the
	//  max()-approximation is rather used
	//-------------------------------------------------------------------------------------------------------------------
	inline double sLogAdd(double logX, double logY) {return (fabs(logY-logX) < 706.0) ? sclib::logAdd(logX, logY) : ((logX>logY) ? logX : logY);}

	//-------------------------------------------------------------------------------------------------------------------
	//  subtracts two numbers in the log-domain where a subtraction in the linear domain is needed (involves one log() and 
	//  one exp() function call)
	//  if the difference between the two arguments is so great that problems may arise during exp()-computation, the
	//  standard computation is rather used
	//-------------------------------------------------------------------------------------------------------------------
	inline double sLogSub(double logX, double logY) {return (logY-logX<0.0 && logY-logY>-706.0) ? sclib::logSub(logX, logY) : exp(logX) - exp(logY);}

	//-------------------------------------------------------------------------------------------------------------------
	// The signum-functions returns the sign of it's argument
	//-------------------------------------------------------------------------------------------------------------------
	inline int sg(double x) {return (x < 0) ? -1 : 1;}

	//-------------------------------------------------------------------------------------------------------------------
	//  the sigmoid function in its basic form has an s-shape, centered around 0 (where it has the value of 0.5), below
	//  this it drops fastly towards zero, above it aproaching 1 (p is a scaling parameter that makes the curve steeper 
	//  as it is increased):
	//
	//                  y ^
	//                    |          
	//                  1 |  __ _ _----------
	//                    | /                                        1
	//                    ||                           sig(x) = ------------    (inversion: 1-sig(x))
	//                 .5 /                                      1 + e^(-p*x)
	//                   ||
	//  __________- - --/ |
	//  -------------------------------------> 
	//  -inf              0              inf x
	//
	//  the 2nd version lets one define the x-position of the infliction-point (pointFiveAt), the point where it shall 
	//  reach 1 (in fact, where it is >=1-eps) and where it shall reach zero (i.e. <=0+eps); if zer0At>oneAt, the 
	//  function is inversed such that it comes from 1 and approaches zero with rising x
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API inline double sigmoid(double x, double p = 1.0) {return 1.0 / (1.0 + exp(-p * x));}
	SCLIB_API double sigmoid(double x, double zeroAt, double pointFiveAt, double oneAt, double eps = 0.000001, double p = 1.0);

	//====================================================================================================================
	//  this is the inverted version of the sigmoidal function: give a sigmoid value [0..1] and get the x value according
	//  to the parameters: minAt describes the x-value for y=0, midAt the x-value for y=0.5, maxAt the value for
	//  y=1
	//====================================================================================================================
	SCLIB_API double invSigmoid(double y, double eps = 0.000001, double p = 1.0);
	SCLIB_API double invSigmoid(double y, double minAt, double midAt, double maxAt, double eps = 0.000001, double p = 1.0);

	//-------------------------------------------------------------------------------------------------------------------
	//  calculates the binomial coefficient, i.e.  /n\, with an iterative programm based on one by B.R.Preiss, see 
	//                                             \k/
	//  http://www.brpreiss.com/books/opus4/html/page467.html; the return value is a double to avoid overflow because the
	//  numbers grow big quick here, the downside is: rounding off errors!
	//-------------------------------------------------------------------------------------------------------------------
	double binomi(unsigned long int n, unsigned long int k);

	//-------------------------------------------------------------------------------------------------------------------
	//  generates a random number between min and max
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API double rand(double min, double max);

	//-------------------------------------------------------------------------------------------------------------------
	//  generates a random index between 0 and maxIdx
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API unsigned int rand(unsigned int maxIdx);

	//-------------------------------------------------------------------------------------------------------------------
	//  generates n random numbers according to the given multivariate gaussian distribution (mean and covariance matrix)
	//-------------------------------------------------------------------------------------------------------------------
	double** randN(unsigned int n, unsigned int dim, double *mean, double **covar);

	//-------------------------------------------------------------------------------------------------------------------
	//  given an array with probabilities (i.e. single values between 0..1 that altogether sum up to 1), an index into 
	//  this array is drawn according to the given distribution. this is useful e.g. for randomly choosing a mixture 
	//  component of a gmm (or a state in a hmm) according to the mixture weights (state transition probabilities).
	//-------------------------------------------------------------------------------------------------------------------
	unsigned int drawIndexFromDistribution(double *weight, unsigned int weightCount);

	//-------------------------------------------------------------------------------------------------------------------
	//	rounds the given double value to the next integer: return the greater int if the part after the comma is greater
	//  or equal 0.5, otherwise the smaller int
	//-------------------------------------------------------------------------------------------------------------------
	inline long int round(double value) {return (long int)((value > 0) ? value+0.5 : value-0.5);}

	//-------------------------------------------------------------------------------------------------------------------
	//  z-transformation of a normal distributed variable to fit the standard normal distribution
	//-------------------------------------------------------------------------------------------------------------------
	inline double zTransform(double x, double mean, double sd) {return (x - mean) / sd;}

	//-------------------------------------------------------------------------------------------------------------------
	//	checks (and returns an altered, definitively valid version) of the given fftSize; it must be a power of 2 and
	//  greater or equal to the given signal length
	//-------------------------------------------------------------------------------------------------------------------
	int checkFftSize(int givenFftSize, int givenSignalLength);

	//-------------------------------------------------------------------------------------------------------------------
	//	returns the "next" (next representable) double that is greater than x
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API double incrementDouble(double x);

	//-------------------------------------------------------------------------------------------------------------------
	//	returns the "next" (next representable) double that is smaller than x
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API double decrementDouble(double x);

	//-------------------------------------------------------------------------------------------------------------------
	//	Scale the given value, which belongs to the interval [oldIntervalMin..oldIntervalMax], to the new interval
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> double scaleToInterval(T value, T oldIntervalMin, T oldIntervalMax, T newIntervalMin, T newIntervalMax) {
		double tmp = (value-oldIntervalMin) / (double)(oldIntervalMax-oldIntervalMin); //scale to [0..1]

		return tmp * (newIntervalMax-newIntervalMin) + newIntervalMin; //scale from [0..1] to [newIntervalMin..newIntervalMax]
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	return the ratio of value1 and value2, but such that (a) both values are taken without sign, i.e. their absolute 
	//  valuale, and (b), the greater absolute value is always the nominator, i.e. the result is >= 1.0
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> double absRatio(T value1, T value2) {
		double abs1 = fabs((double)(value1)), abs2 = fabs((double)(value2));

		return (abs1 > abs2) ? abs1/abs2 : abs2/abs1;
	}

	//-------------------------------------------------------------------------------------------------------------------
	// Tests if the given value is a power of 2
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> bool isPowerOfTwo(T value) {
		double tmp = sclib::log2(value);
		return (tmp > floor(tmp)) ? false : true;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	checks a given value for if it is finite or not
	// 
	//	the idea is from the following newsgroup-posting:
	//		From:Igor Tandetnik (itandetnik@mvps.org)
	//		Subject:Re: Detecting NaN and Infinity
	//		Newsgroups:microsoft.public.vc.language
	//		Date:2003-04-29 08:49:56 PST
	//		
	//	http://groups.google.de/groups?hl=de&lr=&ie=UTF-8&threadm=uz9LFbmDDHA.2704%40TK2MSFTNGP11.phx.gbl&rnum=4&prev=/groups%3Fhl%3Dde%26lr%3D%26ie%3DUTF-8%26q%3D_finite(%2Bportable
	//	2004-06-14
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> bool isFinite(T value) {
		bool res = false;
	  
		res = (value == value); //check for +-1.#IND
		res = res && !(!(0. <= value) && !(value <= 0.)); //check for not NaN
		res = res && !(value == std::numeric_limits<T>::infinity()); //check for not 1.#INF
		res = res && !(value == -1 * std::numeric_limits<T>::infinity()); //check for not -1.#INF

		return res;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	return the minimum of two values ('<'-opetator must exist!); the 2nd value is converted to the type of the 1st one
	//-------------------------------------------------------------------------------------------------------------------
	template<class T, class C> inline T min(T x, C y) {
		return ((x < (T)(y)) ? x : (T)(y));
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	return the maximum of two values ('>'-opetator must exist!); the 2nd value is converted to the type of the 1st one
	//-------------------------------------------------------------------------------------------------------------------
	template<class T, class C> inline T max(T x, C y) {
		return ((x > (T)(y)) ? x : (T)(y));
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	return x or the lower or upper bound, if x exceeds it (<, > operators must exist)
	//-------------------------------------------------------------------------------------------------------------------
	template<class T, class C> inline T getBetween(T lowerBound, C x, T upperBound) {
		return sclib::max((T)(sclib::min(x, upperBound)), lowerBound);
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	add a new item to an already calculated mean value
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> double addToMean(T mean, unsigned long int count, T toAdd) {
		double newMean;
		
		newMean = (double)(count) * (double)(mean); //remove old normalization
		newMean += (double)(toAdd); //add new item
		newMean /= (double)(count+1); //renormalize

		return newMean;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	subtract an item from an already calculated mean value
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> double subtractFromMean(T mean, unsigned long int count, T toSubtract) {
		double newMean;
		
		newMean = (double)(count) * (double)(mean); //remove old normalization
		newMean -= (double)(toSubtract); //subtract item
		newMean /= (double)(count-1); //renormalize

		return newMean;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	returns true if the value lies in the interval between min & max
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> inline bool isBetween(T minValue, T value, T maxValue) {
		return (value >= minValue && value <= maxValue) ? true : false;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//  given are the bounaries of two segments; returns amount of intersection if the two segments intersect
	//-------------------------------------------------------------------------------------------------------------------
	unsigned long int intersect(long int start1, long int end1, long int start2, long int end2);

	//-------------------------------------------------------------------------------------------------------------------
	//	template-version of quickSort for any 1D-Array with '<', '>' operators defined; sorts the array in ascending 
	//  order; l is lower bound, r is upper bound of the array (so array is from i=l;i<=r;i++)
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> void quickSort(T* a, int l, int r) {
		int i, j;
		T	x;		//storage for comparison-element
		T	temp; //storage for temp-element for triangle-exchange
		
		if ( l>=r ) { //nothing to do here
			return;
		} 

		i=l;
		j=r;
		x = a[(l+r)/2];

		//start of division-step
		while (i <= j) {
			while (a[i] < x) {i++;}
			while (a[j] > x) {j--;}
			if (i<=j)	{ //do the triangle-exchange
				temp = a[i];
				a[i] = a[j];
				a[j] = temp;

				i++; 
				j--;
			}
		}

		sclib::quickSort (a, l, j);
		sclib::quickSort (a, i, r);

		return;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	template-version of quickSort for any 2D-Array with '<', '>' operators defined; sorts the matrix in ascending 
	//  order; l is lower bound, r is upper bound of the array's 1st dimension (rows, from i=l;i<=r;i++)
	//  dim is nr. of cols (from d=0;d<dim;d++ note the less-than-operator, not less-or-equal here !!!)
	//  sortBy is the col which gives the reference for sorting
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> void quickSort(T** a, int l, int r, int dim, int sortBy) {
		int i, j, k;
		T x;		//storage for comparison-element
		T temp; //storage for temp-element for triangle-exchange
		
		if ( l>=r ) { //nothing to do here
			return;
		} 

		i=l;
		j=r;
		x = a[(l+r)/2][sortBy];

		//start of division-step
		while (i <= j) {
			while (a[i][sortBy] < x) {i++;}
			while (a[j][sortBy] > x) {j--;}
			if (i<=j)	{ //do the triangle-exchange
				for (k = 0; k < dim; k++) {
					temp = a[i][k];
					a[i][k] = a[j][k];
					a[j][k] = temp;
				}

				i++; 
				j--;
			}
		}

		sclib::quickSort (a, l, j, dim, sortBy);
		sclib::quickSort (a, i, r, dim, sortBy);

		return;
	}

	/*
	void mergeCombine(int *array, const int left, const int right, const int pivot) {
		int temp, i = left, j = pivot+1;

		while (j != right+1 && i != j) { //continue until either list runs out
			if (array[j] <= array[i]) { //Move the jth element in front of the ith element
				temp = array[j];
				for (int k = j; k > i; k--) {
					array[k] = array[k-1]; //Shifts elements from i to j one step forward
				}
				array[i] = temp; //Puts final element in the right place
				i++;
				j++;
			} else {
				i++; // Skip to the next element
			}
		}

		return;
	}
	void mergeSort(int *array, int left, int right) {
		int pivot = left + ((right - left) / 2); //middle of left and right, complicated to avoid overflow (see http://googleresearch.blogspot.com/2006/06/extra-extra-read-all-about-it-nearly.html for details)

		if (left != right) {
			mergeSort(array, left, pivot); //First half
			mergeSort(array, pivot + 1, right); //Second half
			mergeCombine(array, left, right, pivot); //Combines two sorted arrays (left -> pivot) and (pivot+1 -> right)
		}

		return;
	}
	*/

	//-------------------------------------------------------------------------------------------------------------------
	//	selects the k-th smallest element in the array; therefor, a modified sclib::quickSort is used
	//  after its application, the k-th smallest ellement is data[k]
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> void selectK(T* data, int l, int r, unsigned long int k) { // name "select" conflicts with unix/linux syscall
		int i, j;
		T	x;		//storage for comparison-element
		T	temp; //storage for temp-element for triangle-exchange
		
		if ( l>=r ) { //nothing to do here
			return;
		} 

		i=l;
		j=r;
		x = data[(l+r)/2];

		//start of division-step
		while (i <= j) {
			while (data[i] < x) {i++;}
			while (data[j] > x) {j--;}
			if (i<=j)	{ //do the triangle-exchange
				temp = data[i];
				data[i] = data[j];
				data[j] = temp;

				i++; 
				j--;
			}
		}

		if (l<=(int)(k) && (int)(k)<=j) {
			selectK (data, l, j, k);
		} else if (i<=(int)(k) && (int)(k)<=r) {
			selectK (data, i, r, k);
		}

		return;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	returns the p-th percentile of the elements of the array; if inPlace == true, the algorithms will be faster, but 
	//  afterwards the elements of data will be rearranged (all indexes below the percentile will be <=, all above will 
	//  be >=)
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T percentile(T* data, unsigned long int dim, unsigned long int whichPercentile, bool inPlace = false) {
		T p, *temp = NULL;
		unsigned long int idx = sclib::max(0, sclib::round(dim * (double)(whichPercentile / 100.0)) - 1);

		if (inPlace == true) {
			temp = data;
		} else {
			MArray_1D(temp, dim, T, "sclib::percentile: temp");
			for (unsigned long int d = 0; d < dim; d++) {
				temp[d] = data[d];
			}
		}

		sclib::selectK(temp, 0, dim-1, idx);
		p = temp[idx];

		if (inPlace == false) {
			MFree_1D(temp);
		}

		return p;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	analyze the data-array and return a new array having 101 elements (0 to 100), and element i being the i-th 
	//  percentile of the data
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T* percentiles(T* data, unsigned long int dim, bool inPlace = false) {
		T *percentiles, *temp;
		unsigned int i;

		MArray_1D(percentiles, 101, T, "sclib::percentiles: percentiles");

		if (inPlace == true) {
			temp = data;
		} else {
			MArray_1D(temp, dim, T, "percentile: temp");
			for (unsigned long int d = 0; d < dim; d++) {
				temp[d] = data[d];
			}
		}

		sclib::quickSort(temp, 0, dim-1);

		percentiles[0] = temp[0] - std::numeric_limits<T>::epsilon();
		for (i = 1; i < 101; i++) {
			percentiles[i] = temp[((i * dim) / 100) - 1];
		}

		if (inPlace == false) {
			MFree_1D(temp);
		}

		return percentiles;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	returns the median of the elements of the array; if inPlace == true, the algorithms will be faster, but 
	//  afterwards the elements of data will be rearranged (all indexes below the median will be <=, all above will be >=)
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T median(T* data, unsigned long int dim, bool inPlace = false) {
		return sclib::percentile(data, dim, 50, inPlace);
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	inverts the order of elements in the array and returns the pointer to the array (no new space needed, pointer 
	//  stays the same!)
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T* invertArray(T* array, unsigned int dim) {
		T tmp;
		
		for (unsigned int i = 0; i < dim/2; i++) {
			tmp = array[i];
			array[i] = array[dim-1-i];
			array[dim-1-i] = tmp;
		}
	
		return array;
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	rotates the contents of an array by shiftIdx positions (positive: to the right, negative: to the left), treating
	//  the array as a circular buffer. no new space is allocated within and the pointer to the original (and final) 
	//  array is returned; idea from http://dsalgo.blogspot.com/2006/07/rotate-array.html
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T* rotateArray(T* array, unsigned int dim, int shiftIdx) {
		//consider the following example: array = [0, 1, 2, 3, 4, 5], dim = 6, shiftIdx = 2 (or, identically, -4), => expected result = [4, 5, 0, 1, 2, 3]
		sclib::invertArray(array, dim); //invert the hole array: '01|2345' => '54|3210'
		
		//convert multiple rounds of shifting into only one round, and convert left shift (negative shiftIdx) to appropriate right-shift with same result
		if (shiftIdx < 0) {
			if (shiftIdx < -(int)(dim)) {
				shiftIdx = (shiftIdx % (int)(dim)) * -1;
			}
			shiftIdx = dim + shiftIdx; //shiftIdx is negative here, so addition will do			
		} else {
			if ((unsigned int)(shiftIdx) > dim) {
				shiftIdx %= (int)(dim);
			}
		}
		
		sclib::invertArray(array, shiftIdx); //invert first sublist: '54|3210' => '45|3210'
		sclib::invertArray(&(array[shiftIdx]), dim-shiftIdx); //invert second sublist, e.g. '45|3210' => '45|0123'
		
		return array;
	}

	//===================================================================================================================
	//  auxiliary bit functions
	//  a good ressource for further stuff is http://graphics.stanford.edu/~seander/bithacks.html 
	//===================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
	//	shortcut to (int)pow(2,position) to set the bit# position in an integer
	//  eg: bit(3)|bit(5)|bit(7) is an int with bits 3,5 and 7 set, all others are 0.
	//-------------------------------------------------------------------------------------------------------------------
	inline unsigned long int bitValue(unsigned int position) {return (position < 32) ? sclib::pow_2[position] : sclib::round(pow(2.0, (double)(position)));}
	inline unsigned long int bit(unsigned int position) {return sclib::bitValue(position);} //just another name for bitValue()

	//-------------------------------------------------------------------------------------------------------------------
	//	shortcut to (int)log2((double)value) to get the setted sclib::bit# of the integer "value" that has to be a power 
	//  of 2 eg: bitPosition(64) = 6
	//-------------------------------------------------------------------------------------------------------------------
	inline unsigned int bitPosition(unsigned long int value) {unsigned int c = 0;	while (value >>= 1) c++; return c;}	
	inline unsigned int invBit(unsigned long int value) {return sclib::bitPosition(value);} //just the old name for bitPosition()

	//-------------------------------------------------------------------------------------------------------------------
	//	shortcut to ((value & test) == test) to test if similar bits are set in the integers "test" and "value"
	//  eg: bitTest(64, 6) = true, bitTest(64, 5) = false
	//-------------------------------------------------------------------------------------------------------------------
	inline bool bitTest(unsigned long int value, unsigned long int test) {return ((value&test) == test) ? true : false;}

	//-------------------------------------------------------------------------------------------------------------------
	// Returns the number of bits equal to "1" in the given integer; optimal for sparse ones; 
	// See http://www-db.stanford.edu/~manku/bitcount/bitcount.html (2006-10-03) for more algorithms and explanations
	//-------------------------------------------------------------------------------------------------------------------
	unsigned int bitCount(unsigned long int n);

	//====================================================================================================================
	//  auxiliary string functions
	//====================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
	//	a "const char*" comparison operator for maps as used in SC_TweakableParameters
	//-------------------------------------------------------------------------------------------------------------------
	struct ltstr {
		bool operator() (const char *a, const char *b) const {
			return strcmp(a, b) < 0;
		}
	};

	//-------------------------------------------------------------------------------------------------------------------
	//	read a complete, single line from an ascii-file
	//-------------------------------------------------------------------------------------------------------------------
	unsigned int readline(FILE* file, char* buffer, unsigned int maxlen);

	//-------------------------------------------------------------------------------------------------------------------
	//	return the number of lines in the given file; 0 in case of error or empty file
	//-------------------------------------------------------------------------------------------------------------------
	unsigned long int countLines(const char* fileName);

	//-------------------------------------------------------------------------------------------------------------------
	//	check if a character contains a number (didn't find this in the stdlib :-)
	//-------------------------------------------------------------------------------------------------------------------
	bool isNum(char test);

	//-------------------------------------------------------------------------------------------------------------------
	//	converts argument to bool (didn't find this in the stdlib :-)
	//-------------------------------------------------------------------------------------------------------------------
	inline bool atob(const char* value) {return (atoi(value)!=0) ? true : false;}

	//-------------------------------------------------------------------------------------------------------------------
	//	cut the frontmost integer out of a string containing numbers separated by whitespaces
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API int getNextIntFromString(char* buffer, int size, const char *separators = " ,;\t\n\0\r");

	//-------------------------------------------------------------------------------------------------------------------
	//	extract the frontmost integer (containing numbers separated by whitespaces) out of a string starting at startPos, 
	//  return last read position to provide a new starting-pos
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API int getNextIntFromString(const char* buffer, int size, int &res, int startPos = 0, const char *separators = " ,;\t\n\0\r");

	//-------------------------------------------------------------------------------------------------------------------
	//	cut the frontmost double out of a string containing numbers separated by whitespaces
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API double getNextDoubleFromString(char* buffer, int size, const char *separators = " ,;\t\n\0\r");

	//-------------------------------------------------------------------------------------------------------------------
	//	extract the frontmost double (containing numbers separated by whitespaces) out of a string starting at startPos, 
	//  return last read position to provide a new starting-pos
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API int getNextDoubleFromString(const char* buffer, int size, double &res, int startPos = 0, const char *separators = " ,;\t\n\0\r");

	//-------------------------------------------------------------------------------------------------------------------
	//	cut the frontmost string out of a string containing strings separated by whitespaces
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API char* getNextStringFromString(char *buffer, int size, const char *separators = " ,;\t\n\0\r");

	//-------------------------------------------------------------------------------------------------------------------
	//	returns the remainder of inLine after a commentPrefix-sign, if any; otherwise returns NULL
	//-------------------------------------------------------------------------------------------------------------------
	char* extractComment(const char* inLine, int size, const char commentPrefix = '#');

	//-------------------------------------------------------------------------------------------------------------------
	//  returns strIn without heading/trailing whitespaces (ASCII <= 32); a pointer to a new string is returned
	//-------------------------------------------------------------------------------------------------------------------
	char* lTrim(const char* strIn);
	char* rTrim(const char* strIn);
	char* trim(const char* strIn);
	char* trimInPlace(char* strIn, bool left = true, bool right = true);

	//-------------------------------------------------------------------------------------------------------------------
	// Converts all lower-case lettes in 'text' to upper-case
	//-------------------------------------------------------------------------------------------------------------------
	char* uCase(char *text);

	//-------------------------------------------------------------------------------------------------------------------
	// Converts all upper-case lettes in 'text' to lower-case
	//-------------------------------------------------------------------------------------------------------------------
	char* lCase(char *text);

	//-------------------------------------------------------------------------------------------------------------------
	//	A simple 'like'-function (and an auxiliary function) written originally by Stefan Hogedal. It is supposed to mimic 
	//  the behaviour of MS SQL Server 1.1
	// 
	//  The code comes from the newsgroup-thread "C++ has not a Like Operator?", 9 Apr. 1997 09:00, from 
	//  microsoft.public.win32.programmer.gdi 
	//  http://groups.google.de/group/microsoft.public.win32.programmer.gdi/browse_frm/thread/c1285fe694e5c900?page=end&q=like+operator+%22c%2B%2B%22&hl=de&
	//  24.09.2005
	//
	//  Additions by thilo: Character-escaping within the likeExpr with the '\'-character
	//-------------------------------------------------------------------------------------------------------------------
	bool like(const char *text, const char *likeExpr); 
	bool inSet(const char *likeExpr, char c); 

	//-------------------------------------------------------------------------------------------------------------------
	// Compares the two strings on max. len characters ignoring case; return values are the same as for strncmp()
	//-------------------------------------------------------------------------------------------------------------------
	int strincmp(const char *str1, const char *str2, int len);

	//-------------------------------------------------------------------------------------------------------------------
	//  searches 'strIn' and replaces all occurences of 'in' with 'out'; changes the value of the argument, doesn't 
	//  allocate new space or touch the pointer
	//-------------------------------------------------------------------------------------------------------------------
	void strReplace(const char *strIn, char in, char out);

	//-------------------------------------------------------------------------------------------------------------------
	//	this function takes the fileName and changes it's extension (if anyone, else: appends the new one) to the new one.
	//  the space for the new fileName is allocated by this function and the pointer is returned, the 2 parameters remain
	//  untouched. The new extension is meant to be given with the heading '.'
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API char* exchangeFileExtension(const char *fileName, const char *newExtension);

	//-------------------------------------------------------------------------------------------------------------------
	//	returns a new string (space is allocated within, original pointer remains untouched) containing only a filename
	//  without path information
	//-------------------------------------------------------------------------------------------------------------------
	char* extractFileName(const char *fullPath);
	
	//-------------------------------------------------------------------------------------------------------------------
	//	returns a new string (space is allocated within, original pointer remains untouched) containg only a path
	//  without trailing fileName; if giveLastSlash is true (the default), the last charcter will be the '/', so that a 
	//  filename can immediately be added to this path and no information is lost by splitting a full path into its 
	//  ingredients, filename and path.
	//-------------------------------------------------------------------------------------------------------------------
	char* extractPath(const char *fullPath, bool giveLastSlash = true);

	//-------------------------------------------------------------------------------------------------------------------
	//	extract and return the extension of the given filename in a new string
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API char* extractExtension(const char *fileName);

	//-------------------------------------------------------------------------------------------------------------------
	//	add a postfix to a given filename by inserting it directly before the dot introducing the extension (or directly 
	//  at the end, if there is no extension); a new string is returned
	//-------------------------------------------------------------------------------------------------------------------
	char* addPostfixToFilename(const char *fileName, const char *postfix);

	//-------------------------------------------------------------------------------------------------------------------
	//	takes something that should be a path and returns it in a from that can be used as a path throughout this 
	//  library; this means: convert '\' and '\\' to '/' and adda trailing '/'
	//-------------------------------------------------------------------------------------------------------------------
	char* makePath(const char *pathCandidate);

	//-------------------------------------------------------------------------------------------------------------------
	//	takes a string of the form "key=value" and return the key and value in newly allocated buffers
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API bool extractKeyValue(const char *pair, char* &key, char* &value);

	//====================================================================================================================
	//  auxiliary file handling functions
	//====================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
	//	checks if a given filename corresponds with a real existing file
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API bool fileExists(const char* fileName);

	//-------------------------------------------------------------------------------------------------------------------
	//	checks if a given pathname corresponds with a real existing path which is accessible for writing
	//-------------------------------------------------------------------------------------------------------------------
	SCLIB_API bool pathExists(const char* pathName);

	//====================================================================================================================
	//  auxiliary memory allocation/de-allocation functions and macros
	//====================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
	//	3-d array initialization
	//-------------------------------------------------------------------------------------------------------------------
	#define MArray_3D(Target, Dim1, Dim2, Dim3, Type, Msg) {  \
		Type *_Pnt, **_Pnt2; int _Cnt; \
		if ( (Dim1)<=0 || (Dim2) <= 0 || (Dim3) <= 0 || (_Pnt = new Type[(Dim1)*(Dim2)*(Dim3)])==NULL) {(*FunPoint)(SVLIB_NoMem, Msg, __FILE__, __LINE__);} \
		if ( (_Pnt2  = new Type* [Dim1*Dim2]) == NULL) {(*FunPoint)(SVLIB_NoMem, Msg, __FILE__, __LINE__);} \
		for(_Cnt = 0; _Cnt < Dim1*Dim2; _Cnt++) {_Pnt2[_Cnt] = _Pnt + _Cnt * (Dim3);} \
		if ( (Target = new Type**[Dim1]) == NULL) {(*FunPoint)(SVLIB_NoMem, Msg, __FILE__, __LINE__);} \
		for(_Cnt = 0; _Cnt < Dim1; _Cnt++) {Target[_Cnt] = _Pnt2 + _Cnt * (Dim2);} \
	}

	//-------------------------------------------------------------------------------------------------------------------
	//	3-d array destruction
	//-------------------------------------------------------------------------------------------------------------------
	#define MFree_3D(Target) { if (Target != NULL) {delete [] (**Target); delete [] (*Target);  delete [] (Target); Target = NULL;}}

	//-------------------------------------------------------------------------------------------------------------------
	//	0-d array destruction
	//-------------------------------------------------------------------------------------------------------------------
	#define MFree_0D(Target) { if (Target != NULL) {delete (Target); Target = NULL;}}

	//-------------------------------------------------------------------------------------------------------------------
	//	array-in-array destruction (like a matrix, but not allocated in one pice)
	//-------------------------------------------------------------------------------------------------------------------
	#define MFree_2Dex(Target, length) { if (Target != NULL) {for (unsigned long int i = 0; i < (unsigned long int)(length); i++) {delete[] Target[i];}  delete[] (Target); Target = NULL;}}

	//====================================================================================================================
	//  functions to make work with various linked list data-types easier (types are assumed to have a ->Valid() function 
	//  and a ->Next Pointer)
	//====================================================================================================================

	//--------------------------------------------------------------------------------------------------------------------
	//	returns a pointer to the last element of a linked list of objects;
	//	the object must be a pointer to a class containing a 'Next'-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T getLastInList(T pFirst) {
		T pHook = pFirst;

		if (pHook != NULL) {
			while (pHook->Next != NULL) {
				pHook	= (T)(pHook->Next);
			}
		}

		return pHook;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	counts the elements in a linked list of objects
	//	the object must be a pointer to a class containing a 'Next'-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> unsigned long int getListCount(T pFirst) {
		unsigned long int count = 1;
		T pHook = pFirst;
		
		if (pHook != NULL) {
			while (pHook->Next != NULL) {
				count++;
				pHook	= (T)(pHook->Next);
			}
		} else {
			count = 0;
		}

		return count;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	perform a kind of random-access on a linked list of objects (slow if list is long, though)
	//	the object must be a pointer to a class containing a 'Next'-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T getListWithIndex(T pFirst, unsigned long int index) {
		unsigned long int	count = 1;
		T pHook = pFirst;
		
		if (pHook != NULL) {
			while ((pHook->Next != NULL) && (count-1 < index)) {
				count++;
				pHook	= (T)(pHook->Next);
			}
		}

		return pHook;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	returns the index (first element has index 0, 2nd has index 1, ...) of the item in the given list or SVLIB_Fail, 
	//  if it is not in the list
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> long int getListItemIndex(T pFirst, T pItem) {
		long int res = SVLIB_Fail, count = 0;
		T pHook = pFirst;
		
		while (pHook != NULL) {
			if (pHook == pItem) {
				res = count;
				break;
			}
			count++;
			pHook	= (T)(pHook->Next);
		}

		return res;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	returns a pointer to the list-item just before the given one
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T getPreviousListItem(T pFirst, T pItem) {
		long int res = SVLIB_Fail, count = 0;
		T pHook = pFirst, pPrevious = NULL;
		
		if (pHook != NULL) {
			while (pHook->Next != NULL) {
				if (pHook->Next == pItem) {
					pPrevious = pHook;
					break;
				}
				pHook	= (T)(pHook->Next);
			}
		}

		return pPrevious;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	remove an object from a linked list of objects and return the (maybe new) beginning of the list
	//	the object must be a pointer to a class containing a 'Next'-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T removeFromList(T pFirst, unsigned long int indexToRemove, unsigned long int listCount) {
		unsigned long int	count = 1;
		T pPrevious;
		T pNext;
		T pToRemove = sclib::getListWithIndex(pFirst, indexToRemove);
		T pNewStart = pFirst;
		
		//get pointers to the previous and next object
		if (indexToRemove > 0) {
			if (indexToRemove < listCount) {
				pPrevious = sclib::getListWithIndex(pFirst, indexToRemove-1);
				pNext = pToRemove->Next;
			} else {
				pPrevious = sclib::getListWithIndex(pFirst, indexToRemove-1);
				pNext = NULL;
			}
		} else {
			pPrevious = NULL;
			pNext = pToRemove->Next;
		}

		//re-link the list
		pToRemove->Next = NULL;
		if (pPrevious == NULL) {
			if (pNext == NULL) {
				pNewStart = NULL;
			} else {
				pNewStart = pNext;
			}
		} else {
			pPrevious->Next = pNext;
		}

		return pNewStart;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	remove an object from a linked list of objects and return the (maybe new) beginning of the list; in pRemoved the
	//  pointer to the removed elemt is returned
	//	the object must be a pointer to a class containing a 'Next'-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T removeFromList(T pFirst, unsigned long int indexToRemove, unsigned long int listCount, T &pRemoved) {
		unsigned long int	count = 1;
		T pPrevious;
		T pNext;
		T pToRemove = sclib::getListWithIndex(pFirst, indexToRemove);
		T pNewStart = pFirst;
		
		//get pointers to the previous and next object
		if (indexToRemove > 0) {
			if (indexToRemove < listCount) {
				pPrevious = sclib::getListWithIndex(pFirst, indexToRemove-1);
				pNext = pToRemove->Next;
			} else {
				pPrevious = sclib::getListWithIndex(pFirst, indexToRemove-1);
				pNext = NULL;
			}
		} else {
			pPrevious = NULL;
			pNext = pToRemove->Next;
		}

		//re-link the list
		pToRemove->Next = NULL;
		if (pPrevious == NULL) {
			if (pNext == NULL) {
				pNewStart = NULL;
			} else {
				pNewStart = pNext;
			}
		} else {
			pPrevious->Next = pNext;
		}

		pRemoved = pToRemove;

		return pNewStart;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	destructs a whole linked list of objects, if maxLength = 0; otherwise only maxLength items are destructed and a 
	//  pointer to remaining items (or NULL if none) is returned.
	//  the class must have a Valid()-function and a Next-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T destructLinkedList(T& pFirst, unsigned long int maxLength = 0) {
		T pHook;
		T pLast;
		unsigned long int count = 0;

		pHook = pFirst;
		if ((pHook != NULL) && (pHook->Valid() != 0)) {
			do {
 				pLast = pHook;
				pHook = (T)(pHook->Next);
				MFree_0D(pLast);
				count++;
			} while (pHook != NULL && (count < maxLength || maxLength == 0));
		}

		pFirst = NULL;

		return pHook;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	creates a copy of a linked list of objects; the class must have a copy-constructor and a Next-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T* copyLinkedList(T* pFirst) {
		T* pHook;
		T* pNewFirst = NULL;
		T* pNewHook;

		if (pFirst != NULL) {
			pNewFirst = new (T)(*pFirst);
			pHook = pFirst;
			pNewHook = pNewFirst;

			while (pHook->Next != NULL) {
				pNewHook->Next = new (T)(*(pHook->Next));
				pHook = pHook->Next;
				pNewHook = pNewHook->Next;
			}

			pNewHook->Next = NULL;
		}

		return pNewFirst;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	creates a copy of a linked list of objects; the class must have a copy-constructor and a Next-pointer, and the 
	//  copy constructor receives 2 booleans as 2nd and 3rd arguments
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T* copyLinkedList(T* pFirst, bool arg2, bool arg3) {
		T* pHook;
		T* pNewFirst = NULL;
		T* pNewHook;

		if (pFirst != NULL) {
			pNewFirst = new (T)(*pFirst, arg2, arg3);
			pHook = pFirst;
			pNewHook = pNewFirst;

			while (pHook->Next != NULL) {
				pNewHook->Next = new (T)(*(pHook->Next), arg2, arg3);
				pHook = pHook->Next;
				pNewHook = pNewHook->Next;
			}

			pNewHook->Next = NULL;
		}

		return pNewFirst;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	destroys the connections between members of a linked list; 
	//  the class must have a Valid()-function and a Next-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> void unchainLinkedList(T* pFirst) {
		T pHook;
		T pLast;

		pHook = pFirst;
		if ((pHook != NULL) && (pHook->Valid() != 0)) {
			do {
 				pLast = pHook;
				pHook = (T)(pHook->Next);
				pLast->Next = NULL;
			} while (pHook != NULL);
		}

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	adds the new element to the end of the list given by it's first element; if the list is empty, pNew is the new 
	//  first element. the class must have a Valid()-function and a Next-pointer
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> void addToList(T* &pFirst, T* pNew) {
		if (pFirst == NULL) {
			pFirst = pNew;
		} else {
			getLastInList(pFirst)->Next = pNew;
		}

		return;
	}

	//====================================================================================================================
	//	functions to facilitate (debug) output: accept a data-type and write it into an ASCII-file
	//====================================================================================================================

	//--------------------------------------------------------------------------------------------------------------------
	//	print a string to a file
	//--------------------------------------------------------------------------------------------------------------------
	SCLIB_API const char* nullFilter(const char *text);

	//--------------------------------------------------------------------------------------------------------------------
	//	print a string to a file
	//--------------------------------------------------------------------------------------------------------------------
	SCLIB_API void stringOut(const char *fileName, const char *string, SC_TweakableParameters *pTweak = NULL, const char* separator = "\n");

	//--------------------------------------------------------------------------------------------------------------------
	//	print a scalar to a file
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void scalarOut(const char* fileName, T scalar, SC_TweakableParameters *pTweak = NULL, bool fullDetail = false, const char* separator = "\n") {
		fstream	fileOut;
		char* fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);
		if (fullDetail == false) {
			fileOut << setw(8) << scalar << separator;
		} else {
			fileOut << setprecision(32) << scalar << separator;
		}
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a scalar, annotated by a \0-terminated string, to a file
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void scalarOutEx(const char* fileName, T scalar, const char* annotation, SC_TweakableParameters *pTweak = NULL) {
		fstream	fileOut;
		char* fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);
		fileOut << annotation << ": " << setw(8) << scalar << "\n";
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print the content of a map (std::vector) to a file
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void mapOut(const char* fileName, std::vector<T>* map, SC_TweakableParameters *pTweak = NULL, const char* separator = "\n") {
		fstream	fileOut;
		char* fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);
		for (long int t = 0; t < (long int)(map->size()); t++) {
			fileOut << (*map)[t] << separator;
		}
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a tupel to a file, using two kinds of separators: first one between the tupel-elements, second one between
	//  different tupels
	//--------------------------------------------------------------------------------------------------------------------
	template<class T, class C> void tupelOut(const char* fileName, T value1, C value2, SC_TweakableParameters *pTweak = NULL, const char* separator1 = "; ", const char* separator2 = "\n") {
		fstream	fileOut;
		char* fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);
		fileOut << value1 << separator1 << value2 << separator2;
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a tripel to a file, using two kinds of separators: first one between the tripel-elements, second one between
	//  different tripels
	//--------------------------------------------------------------------------------------------------------------------
	template<class T, class C, class D>	void tripelOut(const char* fileName, T value1, C value2, D value3, SC_TweakableParameters *pTweak = NULL, const char* separator1 = "; ", const char* separator2 = "\n") {
		fstream	fileOut;
		char* fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);
		fileOut << value1 << separator1 << value2 << separator1 << value3 << separator2;
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a quadrupel to a file, using two kinds of separators: first one between the tripel-elements, second one 
	//  between different quadrupels
	//--------------------------------------------------------------------------------------------------------------------
	template<class T, class C, class D, class E>	void quadrupelOut(const char* fileName, T value1, C value2, D value3, E value4, SC_TweakableParameters *pTweak = NULL, const char* separator1 = "; ", const char* separator2 = "\n") {
		fstream	fileOut;
		char* fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);
		fileOut << value1 << separator1 << value2 << separator1 << value3 << separator1 << value4 << separator2;
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a vector to a file
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void vectorOut(const char* fileName, T* vector, unsigned long int dim, bool printVertical = false, SC_TweakableParameters *pTweak = NULL, const char* separator = " ") {
		unsigned long int x;
		fstream fileOut;
		char* fName;
		bool matlab = (strncmp(separator, sclib::matlabSyntax, strlen(sclib::matlabSyntax)) == 0) ? true : false;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);

		if (matlab == true) {
			fileOut << "[";
		}

		for (x = 0; x < dim; x++) {
			if (matlab == true) {
				fileOut << vector[x] << " ";
			} else {
				fileOut << setprecision(32) << vector[x] << (printVertical == true ? "\n" : separator);//setw(32) << vector[x] << (printVertical == true ? "\n" : " ");
			}
		}

		if (matlab == true) {
			fileOut << "]";
		}

		fileOut << "\n";

		fileOut.close();
		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a vector to a file
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void vectorOutEx(const char* fileName, T* vector, unsigned long int dim, bool printVertical = false, SC_TweakableParameters *pTweak = NULL, const char* separator = " ") {
		unsigned long int x;
		fstream fileOut;
		char* fName;
		bool matlab = (strncmp(separator, sclib::matlabSyntax, strlen(sclib::matlabSyntax)) == 0) ? true : false;
		double figure;
		unsigned long int idx;
		SC_MatrixFunctions mFunc;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);

		if (matlab == true) {
			fileOut << "[";
		}

		for (x = 0; x < dim; x++) {
			if (matlab == true) {
				fileOut << vector[x] << " ";
			} else {
				fileOut << setprecision(32) << vector[x] << (printVertical == true ? "\n" : separator);//setw(32) << vector[x] << (printVertical == true ? "\n" : " ");
			}
		}

		if (matlab == true) {
			fileOut << "]";
		}
		fileOut << "\n";

		//additional info
		figure = mFunc.min(vector, dim, &idx);
		fileOut << setw(8) << figure << separator << " (min @ " << idx << ")\n";
		figure = mFunc.max(vector, dim, &idx);
		fileOut << setw(8) << figure << separator << " (max @ " << idx << ")\n";
		figure = mFunc.mean(vector, dim);
		fileOut << setw(8) << figure << separator << " (mean)\n";
		figure = mFunc.variance(vector, dim, &figure);
		fileOut << setw(8) << figure << separator << " (variance)\n";
		figure = median(vector, dim, false);
		fileOut << setw(8) << figure << separator << " (median)\n";

		fileOut.close();
		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a matrix to a file, accompanied by additional information: min/max/mean/variance/median
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> void matrixOut(const char* fileName, T** matrix, unsigned long int len, unsigned long int dim, SC_TweakableParameters *pTweak = NULL, unsigned long int minLen = 0, unsigned long int maxLen = 0, unsigned long int minDim = 0, unsigned long int maxDim = 0, const char* separator = ";", bool printLastNewline = true) {
		unsigned long int x, y;
		fstream fileOut;
		char* fName;
		unsigned long int mxDim = maxDim, mxLen = maxLen, length, dimension;
		bool matlab = (strncmp(separator, sclib::matlabSyntax, strlen(sclib::matlabSyntax)) == 0) ? true : false;

		if (mxLen == 0) {mxLen = len;}
		if (mxDim == 0) {mxDim = dim;}
		length = mxLen - minLen; 
		dimension = mxDim - minDim;

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);

		if (matlab == true) {
			fileOut << "[";
		}

		for ( y = minLen; y < mxLen; y++) { //len == dim_y
			for (x = minDim; x < mxDim; x++) { //dim == dim_x
				if (matlab == true) {
					fileOut << matrix[y][x] << " ";
				} else {
					fileOut << setw(8) << matrix[y][x] << separator;
				}
			}
			if (matlab == true) {
				fileOut << "; ";
			} else {
				fileOut << "\n";
			}
		}

		if (matlab == true) {
			fileOut << "]";
		}

		if (printLastNewline == true) {
			fileOut << "\n";
		}

		fileOut.close();
		MFree_1D(fName);

		return;	
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a matrix to a file, accompanied by additional information for each column: min/max/mean/variance
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> void matrixOutEx(const char* fileName, T** matrix, unsigned long int len, unsigned long int dim, SC_TweakableParameters *pTweak = NULL, const char* separator = ";") {
		unsigned long int x, y;
		fstream fileOut;
		char* fName;
		double *mean = NULL, *variance = NULL, *min = NULL, *max = NULL;
		SC_MatrixFunctions mFunc;
		bool matlab = (strncmp(separator, sclib::matlabSyntax, strlen(sclib::matlabSyntax)) == 0) ? true : false;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		
		fileOut.open(fName, ios_base::out|ios_base::app);

		if (matlab == true) {
			fileOut << "[";
		}

		for ( y = 0; y < len; y++) { //len == dim_y
			for (x = 0; x < dim; x++) { //dim == dim_x
				if (matlab == true) {
					fileOut << matrix[y][x] << " ";
				} else {
					fileOut << setw(8) << matrix[y][x] << separator;
				}
			}
			if (matlab == true) {
				fileOut << "; ";
			} else {
				fileOut << "\n";
			}
		}

		if (matlab == true) {
			fileOut << "]";
		}

		fileOut << "\n";

		//additional info
		min = mFunc.min(matrix,len, dim);
		for (x = 0; x < dim; x++) { //dim == dim_x
			fileOut << setw(8) << min[x] << (matlab==false?separator:"\n");
		}
		fileOut << " (min)\n";
		MFree_1D(min);
		max = mFunc.max(matrix, len, dim);
		for (x = 0; x < dim; x++) { //dim == dim_x
			fileOut << setw(8) << max[x] << (matlab==false?separator:"\n");
		}
		fileOut << " (max)\n";
		MFree_1D(max);
		mean = mFunc.mean(matrix, len, dim);
		for (x = 0; x < dim; x++) { //dim == dim_x
			fileOut << setw(8) << mean[x] << (matlab==false?separator:"\n");
		}
		fileOut << " (mean)\n";
		variance = mFunc.variance(matrix, len, dim, mean);
		for (x = 0; x < dim; x++) { //dim == dim_x
			fileOut << setw(8) << variance[x] << (matlab==false?separator:"\n");
		}
		fileOut << " (variance)\n";
		MFree_1D(mean);
		MFree_1D(variance);
		
		fileOut.close();
		MFree_1D(fName);

		return;	
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a class to a file (class must have overloaded '<<'-operator)
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void classOut(const char* fileName, T* pClass, SC_TweakableParameters *pTweak = NULL, sclib::OpenMode mode = ios_base::out|ios_base::app, const char *additionalInfo = "") {
		fstream fileOut;
		char*		fName;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		fileOut.open(fName, mode);
		fileOut << additionalInfo;
		fileOut << *pClass;
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//--------------------------------------------------------------------------------------------------------------------
	//	print a linked list of objects (with a "Next" pointer member) to a file (class must have overloaded '<<'-operator)
	//--------------------------------------------------------------------------------------------------------------------
	template<class T>	void listOut(const char* fileName, T* pList, SC_TweakableParameters *pTweak = NULL, sclib::OpenMode mode = ios_base::out|ios_base::app, const char *additionalInfo = "") {
		fstream fileOut;
		char*		fName;
		T* pHook = pList;
		unsigned long int count = 0;
		
		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}
		fileOut.open(fName, mode);
		fileOut << additionalInfo;
		while (pList != NULL) {
			fileOut << "----------- list element nr. " << count++ << "-----------" << endl;
			fileOut << *pList;
			pList = pList->Next;
		}
		fileOut.close();

		MFree_1D(fName);

		return;
	}

	//====================================================================================================================
	//	functions to facilitate binary input/output
	//====================================================================================================================

	//-------------------------------------------------------------------------------------------------------------------
  //  save a matrix like a model; if the file exists, it will be overwritten
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> int saveMatrix(char *fileName, T **matrix, unsigned long int len, unsigned long int dim, SC_TweakableParameters *pTweak = NULL) {
		SV_DataIO io;
		int bytes;
		char *fName;
		fstream	matFile;
		SV_DataIO::SV_DatatypeSizes codeSizes;
		io.getCurrentDatatypeSizes(codeSizes);

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}

		matFile.open(fName, ios_base::out|ios_base::binary);
		bytes = io.writeMachineHeader(&matFile, codeSizes);
		bytes += io.writeScalar(&matFile, len);
		bytes += io.writeScalar(&matFile, dim);
		bytes += io.writeMatrix(&matFile, matrix, len, dim);
		if (matFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Saving matrix failed!");}
		matFile.clear();
		matFile.close();

		MFree_1D(fName);

		return bytes;
	}

	//-------------------------------------------------------------------------------------------------------------------
  //  load a matrix like a model
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T** loadMatrix(char *fileName, T matType, unsigned long int &len, unsigned long int &dim, SC_TweakableParameters *pTweak = NULL) {
		T **mat = NULL;
		SV_DataIO io;
		int bytes;
		char *fName;
		fstream	matFile;
		SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
		io.getCurrentDatatypeSizes(codeSizes);

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}

		matFile.open(fName, ios_base::in|ios_base::binary);
		bytes = io.readMachineHeader(&matFile, fileSizes);
		bytes += io.readScalar(&matFile, len, codeSizes, fileSizes);
		bytes += io.readScalar(&matFile, dim, codeSizes, fileSizes);
		MArray_2D(mat, (int)(len), (int)(dim), T, "SC_MatrixFunctions.loadMatrix: mat");
		bytes += io.readMatrix(&matFile, mat, len, dim, codeSizes, fileSizes);
		if (matFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Loading matrix failed!");
			MFree_2D(mat);
			len = 0;
			dim = 0;
		}
		matFile.clear();
		matFile.close();

		MFree_1D(fName);

		return mat;
	}

	//-------------------------------------------------------------------------------------------------------------------
  //  save a vector like a model; if the file exists, it will be overwritten
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> int saveVector(char *fileName, T *vector, unsigned long int dim, SC_TweakableParameters *pTweak = NULL) {
		SV_DataIO io;
		int bytes;
		char *fName;
		fstream	vecFile;
		SV_DataIO::SV_DatatypeSizes codeSizes;
		io.getCurrentDatatypeSizes(codeSizes);

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}

		vecFile.open(fName, ios_base::out|ios_base::binary);
		bytes = io.writeMachineHeader(&vecFile, codeSizes);
		bytes += io.writeScalar(&vecFile, dim);
		bytes += io.writeArray(&vecFile, vector, dim);
		if (vecFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Saving vector failed!");}
		vecFile.clear();
		vecFile.close();

		MFree_1D(fName);

		return bytes;
	}

	//-------------------------------------------------------------------------------------------------------------------
  //  load a vector like a model
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> T* loadVector(char *fileName, T vecType, unsigned long int &dim, SC_TweakableParameters *pTweak = NULL) {
		T *vec = NULL;
		SV_DataIO io;
		int bytes;
		char *fName;
		fstream	vecFile;
		SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
		io.getCurrentDatatypeSizes(codeSizes);

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}

		vecFile.open(fName, ios_base::in|ios_base::binary);
		bytes = io.readMachineHeader(&vecFile, fileSizes);
		bytes += io.readScalar(&vecFile, dim, codeSizes, fileSizes);
		MArray_1D(vec, dim, T, "SC_MatrixFunctions.loadVector: vec");
		bytes += io.readArray(&vecFile, vec, dim, codeSizes, fileSizes);
		if (vecFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Loading matrix failed!");
			MFree_1D(vec);
			dim = 0;
		}
		vecFile.clear();
		vecFile.close();

		MFree_1D(fName);

		return vec;
	}

	//-------------------------------------------------------------------------------------------------------------------
  //  save a hash map of scalar keys and values like a model; if the file exists, it will be overwritten
	//-------------------------------------------------------------------------------------------------------------------
	template<class T, class C> int saveScalarMap(char *fileName, std::map<T, C> scalarMap, SC_TweakableParameters *pTweak = NULL) {
		SV_DataIO io;
		int bytes;
		char *fName;
		fstream	mapFile;
		SV_DataIO::SV_DatatypeSizes codeSizes;
		io.getCurrentDatatypeSizes(codeSizes);
		typename std::map<T, C>::iterator i;

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}

		mapFile.open(fName, ios_base::out|ios_base::binary);
		bytes = io.writeMachineHeader(&mapFile, codeSizes);
		bytes += io.writeScalar(&mapFile, (int)(scalarMap.size()));
		for (i = scalarMap.begin(); i != scalarMap.end(); ++i) {
			bytes += io.writeScalar(&mapFile, i->first);
			bytes += io.writeScalar(&mapFile, i->second);
		}
		if (mapFile.good() != TRUE) {REPORT_ERROR(SVLIB_Fail, "Saving scalar map failed!");}
		mapFile.clear();
		mapFile.close();

		MFree_1D(fName);

		return bytes;		
	}
	
	//-------------------------------------------------------------------------------------------------------------------
  //  load a hash map of scalar keys and values like a model
	//-------------------------------------------------------------------------------------------------------------------
	template<class T, class C> bool loadScalarMap(char *fileName, std::map<T, C> scalarMap, SC_TweakableParameters *pTweak = NULL) {
		SV_DataIO io;
		int bytes, size;
		char *fName;
		bool success = true;
		fstream	mapFile;
		SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes;
		io.getCurrentDatatypeSizes(codeSizes);

		if (pTweak != NULL) {
			fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
			sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
		} else {
			fName = new char[strlen(fileName) + 1];
			sprintf(fName, "%s\0", fileName);
		}

		scalarMap.clear();
		mapFile.open(fName, ios_base::in|ios_base::binary);
		bytes = io.readMachineHeader(&mapFile, fileSizes);
		bytes += io.readScalar(&mapFile, size, codeSizes, fileSizes);
		for (int i = 0; i < size; i++) {
			T first;
			C second;
			bytes += io.readScalar(&mapFile, first, codeSizes, fileSizes);
			bytes += io.readScalar(&mapFile, second, codeSizes, fileSizes);
			scalarMap[first] = second;
		}
		if (mapFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Loading matrix failed!");
			scalarMap.clear();
			success = false;
		}
		mapFile.clear();
		mapFile.close();

		MFree_1D(fName);

		return success;
	}

	//====================================================================================================================
	//  other auxiliary functions
	//====================================================================================================================

	//--------------------------------------------------------------------------------------------------------------------
	//	this function sets or retrieves the status of what the error handler does in release mode: throwing an exception
	//  (as needed when the library is called via JNI) or exiting with error code
	//--------------------------------------------------------------------------------------------------------------------
	bool errorHandlerThrows(bool justGet = true, bool newValue = true);

	//--------------------------------------------------------------------------------------------------------------------
	//	new svlib-like error-handler, which allows resuming 
	//--------------------------------------------------------------------------------------------------------------------
	void errorHandler(int ErrorCode, const char* ErrorMsg, const char* FName, int LNum);

	//--------------------------------------------------------------------------------------------------------------------
	// Generate synthetic data to test the models
	//--------------------------------------------------------------------------------------------------------------------
	void generateTestData(SV_Data* &pNoise, SV_Data* &pSpeech, unsigned long int D, unsigned long int J, unsigned long int I);
	
	//-------------------------------------------------------------------------------------------------------------------
	// computes the numer of rows that a feature-set has when a sliding window of given size is moved at given step-size 
	// over the data of given size; windowSize might be a float (as in SC_Feature_PitchESPS), but should normally be 
	// unsigned long int, too
	//-------------------------------------------------------------------------------------------------------------------
	template<class T> unsigned long int getRowCount(unsigned long int dataSize, T windowSize, unsigned long int windowStep) {
		unsigned long int rows = 0;

		//the old formula "(dataSize / windowStep) - (windowSize / windowStep) + 1" is mathematical equivalent to this one,
		//but there is the chance of rounding errors because of implicitly used rounding via integer-arithmetic because 
		//rounding occurs independantly in two terms and not just (as is correct) in one. an example is dataSize=300,
		//windowSize=186, windowStep=5

		if ((T)(dataSize) >= windowSize) {
			rows = (unsigned long int)(dataSize - windowSize + 0.5) / windowStep; //nr. of times the window can be moved while completely fitting into the dataset (integer-arithmetic is correct here!)
			rows++; //even if the window can't be moved, the window fits 1 time into the data-set (the enclosing if(dataSize >= windowSize) ensures that)
		}

		return rows;
	}

	//--------------------------------------------------------------------------------------------------------------------
	// assuming that the currentValue is between 1 and maxValue, the percentage of how much of the maxValue is already 
	// reached is printed on stdOut if there is at least minIncrease percentage points increase as compared to the 
	// lastPercentage; set firstPrint=true if this is the first percentage to be printed in this position so that 
	// indentation is handled correctly; the last printed percentage is returned (for the next lastPercentage parameter)
	//--------------------------------------------------------------------------------------------------------------------
	double printPercentage(unsigned long int maxValue, unsigned long int currentValue, double lastPercentage, double minIncrease, bool firstPrint = false);

	//--------------------------------------------------------------------------------------------------------------------
	// ofton, a method returns a new object that should directly be stored in a pointer to one of its inputs, but because
	// the object pointed to by the input pointer must first be destroyed, a temporary pointer is needed to catch the 
	// output, free the old input and reassign the output (now residing in the tmp pointer) to the input pointer; this 
	// function does it in one line, freeing the caller from the need to declare a temporary pointer. the pointers may 
	// mean pointers to a single object or starting points of arrays (isArray==true)
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T* directOutput(T* output, T* &input, bool isArray = false) {
		if (isArray == true) {
			MFree_1D(input);
		} else {
			MFree_0D(input);
		}

		input = output;

		return input;
	}

	//--------------------------------------------------------------------------------------------------------------------
	// variant of the above function for matrixes and arrays of object-pointer.
	//--------------------------------------------------------------------------------------------------------------------
	template<class T> T** directOutput(T** output, T** &input, bool isMatrix = false) {
		if (isMatrix == true) {
			MFree_2D(input);
		} else {
			MFree_1D(input);
		}

		input = output;

		return input;
	}

	//--------------------------------------------------------------------------------------------------------------------
	// toggle the value of flag between newValue (reset=false) and it's original value (reset=true); only works if flag's
	// address doesn't change (because it is used as a key to find the original value)!
	//--------------------------------------------------------------------------------------------------------------------
	/*template<class T> void toggle(T &flag, bool reset = false, T newValue = (T)(0)) {
		static std::map<T*, T> originalValue;

		if (reset == false) {
			originalValue[&flag] = flag;
			flag = newValue;
		} else {
			flag = originalValue[&flag];
			originalValue.erase(&flag);
		}

		return;
	}*/

} //end of namespace sclib

#endif
