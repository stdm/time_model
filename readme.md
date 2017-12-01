## Release 1.1b of sclib

In mid 2017, I decided to clean up sclib in order to publish code able to reproduce the ACM MM 2009 results (Unfolding Speaker Clustering Potential - a Biomimetic Approach). I removed certain experimental classes, fixed some typos in comments, and tried to make working with the time model more accessible. This process is not finished yet (no single function currently will flawlessly repeat the experiment from the ACM MM 2009 paper), but at least normalization matrix creation and trajectory extraction works directly (and doing the clustering experiment just needs some glue code in SC_MainTasks). While doing so, I converted the windows VC project files to x64 (it compiled & ran without changes to the code); this changes havn't been reflected in the linux GCC project files yet. Other important changes are reflected inline below (e.g., commented-out classes no longer available).

Winterthur, July 2017
Thilo Stadelmann
stdm@zhaw.ch


## Release 1.0b of sclib

This is the first release of the speaker classification library "sclib". It has been developed in the course of my Ph.D. research, and this work is now completed. The software is not mature enough to be called a product, usable by end users/developers. Yet it possibly is of enormous benefit for other researchers working in the field and having (or are striving for) a deep understanding of the corresponding audio analysis techniques. The software is therefore classifed as havin "beta" status. More on this below in Section "status". 

As is the case with other "research prototypes" (I deem the sclib to be more, but lets stick with this label for a moment as a term for software emerging from research projects rather than from product- and user-driven commercial efforts), the sclib's "soft skills" are regrettably underdeveloped. Besides good code, there is not much documentation or installation instructions. In fact, the best overview of its structure, the overall big picture of ist development, its content, aimed audience/task and descriptions of the included algorithms together with their background is given in my thesis (T. Stadelmann, "Voice Modeling Methods for Automatic Speaker Recognition under Adverse Conditions", Marburg, Germany, 2010). Then there is this document that adds some statements about usage and utilizability from a developer point of view. The details are packed into lots of explaining comments in the code itself.

The software is provided as is, without any warranty, expressed or implied. Do with it as you like (but don't violate the copyright of other projects linked/used by sclib), but please cite my thesis (or specific papers) and send me a note as you do so.

The remainder of this readme is structured as follows: "Content of Package" gives a very coarse overview of the directories and their content; "Status" explains the maturity of each class in the sclib; "Windows" explains how to build the system on a windows machine; "Linux" explains how to build it under linux.

Marburg, April 2010
Thilo Stadelmann
stadelmann@informatik.uni-marburg.de


## Content of Package

- /__build_first: code of libraries used by sclib; given for reference (these versions compile and run under linux and windows and work together with the sclib), but maybe its a good idea to use the current versions from the original developers
- /data: i/o directory of the sclib: /data/trained contains trained classifiers/models, /data/test is for debug output; /data contains ini-files with parameter settings for all algorithms in the sclib as well as example files for loading groundtruth and corpora of speech
- /dlls: under windows, all used libraries and headers are located here, and the reesults of compilation end up here, too (v1.1b: libfann is no longer needed)
- /Gardener: a GUI written in Visual Basic; not up to date, but helpful to visualize the speaker clustering results
- /General: the Visual Studio solution files for all related projects
- /sclib: code an project files for the sclib
- /svlib: code an project files for the svlib that is used as a basis for the sclib
- /testbed: a small C++ test program used as a test environment for SC_MainTasks methods


### Status

Parts of the sclib are well tested and used and represent the state of the art in audio processing research. Others mark more or less dead ends of the development, ideas that haven't been carried out until the end or just did not work. Others just did not work with the used parameters and the given datasets, but the code might be ok. Nobody wants to use a software with such mixed up maturities. So, this section gives an overview about the status of each contained class so that you know what you can rely on out of the box and which parts you have to look into very closely before using.

| Class | Status |
| ----- | ------ |
| SC_Api | well tested and applied |
|	SC_Aux | well tested and applied |
|	SC_BaggenstossEM | well tested and applied |
|	SC_BaggenstossEMex | well tested and applied |
|	SC_Centroid | well tested and applied |
|	SC_Centroid_Gaussian | well tested and applied |
|	SC_Centroid_Point | well tested and applied |
|	SC_Centroid_Signature | not tested, should work |
|	SC_Classifier | well tested and applied |
|	SC_Classifier_AdaBoost | tested and applied |
|	SC_Classifier_DecisionStump | tested and applied |
|	SC_Classifier_ML | tested and applied |
|	SC_Classifier_SVM | well tested and applied |
|	SC_ClassifierHandler | well tested and applied |
|	SC_ClassifierTree | well tested and applied |
|	SC_ClassifierWithWeights | well tested and applied |
|	SC_Cluster | well tested and applied |
|	SC_Clusterer | parted: the SV_Data-returning kMeans and its initialization methods are well tested and applied, rest is implemented by student (seems to be dead end) |
|	SC_Conversion | well tested and applied |
|	SC_Corpus | well tested and applied |
|	SC_Corpus_MAC | well tested and applied |
|	SC_Corpus_MPEG7 | well tested and applied |
|	SC_Corpus_SCiVo | well tested and applied |
|	SC_Corpus_TIMIT | well tested and applied |
|	SC_Corpus_Wesley | well tested and applied |
|	SC_DistanceMeasures | well tested and applied |
|	SC_EMD | well tested and applied (ported code from Y. Rubner) |
|	SC_Enhancement | experimental - stopped developing before extensive testing |
|	SC_Feature_BandPeriodicity | well tested and applied |
|	SC_Feature_BrightnessBandwidth | well tested and applied |
|	SC_Feature_FbE | well tested and applied |
|	SC_Feature_Formant | not tested, seems to work (ported ESPS code) |
|	SC_Feature_LPC | tested and applied |
|	SC_Feature_LPCresidual | tested and applied |
|	SC_Feature_LSP | well tested and applied (ported MELP/Speex code) |
|	SC_Feature_MFCC | well tested and applied |
|	SC_Feature_NFR | well tested and applied |
|	SC_Feature_Pitch | well tested and applied (ported ESPS code) |
|	SC_Feature_SDP | tested and applied |
|	SC_Feature_STE | well tested and applied |
|	SC_Feature_Samples | well tested and applied |
|	SC_Feature_Spectrum | well tested and applied |
|	SC_Feature_SpectrumFlux | well tested and applied |
|	SC_Feature_SubBandPower | well tested and applied |
|	SC_Feature_ZCR | well tested and applied |
|	SC_FeatureHandler | well tested and applied |
|	SC_Gauss | well tested and applied |
|	SC_GroundTruth | well tested and applied |
|	SC_GroundTruth_MAC | well tested and applied |
|	SC_GroundTruth_MPEG7 | well tested and applied |
|	SC_GroundTruth_SCiVo | well tested and applied |
|	SC_GroundTruth_TIMIT | well tested and applied |
|	SC_GroundTruth_Wesley | well tested and applied |
|	SC_HHT | implemented by student (haven't verified reliability, but seems to work) |
|	SC_Ini | well tested and applied |
|	SC_JStreamReader | well tested and applied |
|	SC_Lib | well tested and applied |
|	SC_LDA | implemented by student (seems to work) |
|	SC_MainTasks | melting pot of well tested&applied as well as experimental test scripts |
|	SC_Matlab | interface works, application of l1magic did not yield good results |
|	SC_MatrixFunctions | well tested and applied |
|	SC_MD5 | well tested and applied |
|	SC_MixtureModel | well tested and applied |
|	SC_MixtureModel_bGMM | tested and applied |
|	SC_MixtureModel_GMM | well tested and applied |
|	SC_MixtureModel_GMM_UBM | not tested much, but seems to work |
|	SC_MixtureModel_MIX2MAX | tested and applied (in the thesis this called the "AMU") |
|	SC_MixtureModel_MIX2MAXex | tested and applied (in the thesis this are the extensions to the "AMU") |
|	SC_MixtureModel_MIXMAX | well tested and applied |
|	SC_Model | well tested and applied |
|	SC_Model_FullGauss | tested and applied (wrapper for svlib code) |
|	SC_Model_HMM | not tested much (wrapper for svlib code) |
|	SC_Model_MetaGMM | well tested and applied |
|	SC_Model_Pareto | tested and applied |
|	SC_Model_qGMM | applied but not really tested |
|	SC_Model_SVM | tested and applied (arc-distance is experimental) |
|	SC_Model_Time | tested and applied, but only prototype of greater idea |
|	SC_Model_VQ | not tested (wrapper for svlib code) |
|	SC_ModelHandler | well tested and applied |
|	SC_Partition | well tested and applied |
|	SC_PDE | tested and applied |
|	SC_Resampling | well tested and applied |
|	SC_Score | tested an applied |
|	SC_Score_AudioTypeClassification | tested and applied |
|	SC_Score_ChangeDetection | tested and applied |
|	SC_Score_VUv | tested and applied |
|	SC_SDP | well tested and applied |
|	SC_Segmentation_AudioType | well tested and applied |
|	SC_Segmentation_AudioType_LZL | well tested and applied |
|	SC_Segmentation_Changes | all change detectors work unstatisfactory, partly due to no good parameter values found, but bugs in code may remain... |
|	SC_Segmentation_Changes_KBK | well tested, but does not yield what can be expected - bugs or parameter-settings? |
|	SC_Segmentation_Changes_LZW | works so-so |
|	SC_Segmentation_Changes_Std | well tested and applied |
|	SC_Segmentation_Silence | well tested and applied |
|	SC_Segmentation_Silence_LZL | well tested and applied |
|	SC_Segmentation_Silence_LNK | works so-so |
|	SC_Segmentation_VUv | all voiced/unvoiced detectors are somewhat untested (student projects) because in the end it worked better to use all speech for speaker modeling instead of only voiced speech |
|	SC_Segmentation_VUv_ESPS | well tested and applied |
|	SC_Segmentation_VUv_LNK | works so-so |
|	SC_SegmentationHandler | well tested and applied |
|	SC_Signal | well tested and applied |
|	SC_SignalHandler | well tested and applied |
|	SC_Signal_jWAVE | well tested and applied |
|	SC_Signal_MPEG | well tested and applied |
|	SC_Signal_NIST | well tested and applied |
|	SC_Signal_WAVE | well tested and applied |
|	SC_Signature | well tested and applied |
|	SC_SpeakerClusterer | well tested and applied |
|	SC_SpeakerIdentificator | implemented but not much tested/applied |
|	SC_SpeakerScore | tested and applied |
|	SC_SpeakerScore_Classification | tested and applied |
|	SC_SpeakerScore_Clustering | tested and applied |
|	SC_SpeakerScore_Identification | tested and applied |
|	SC_SVM | well tested and applied (ported libsvm code) |
|	SC_Synthesis | well tested and applied |
|	SC_Timer | well tested and applied |
|	SC_Transform | well tested and applied |
|	SC_TweakableParameters | well tested and applied |
	

### Windows

The sclib compiles using Microsoft Visual Studio 2008 (v9.0.x), already using the Express Edition (v1.1b: Visual Studio 2010 Professional v10.0.40219.1 SP1Rel). Tested Windows versions are XP, Vista, 7 (v1.1b: Win10). 

To build the software, first the /dlls directory has to be added to the PATH environment variable in order for windows to find the dll-files of the loaded libraries. The libraries are prebuild for 32bit Win7 using VS8. If they cannont be linked, build the libraries from scratch using for example the code in the /__build_first directory. All libraries but ffmpeg employ visual studio project files, so that building should be no problem after maybe altering the output directories. Ffmeg is a special case - it can only be compiled under windows using MSys and a special procedure described in the file init-windows.sh; additionally, because sclib uses functions from the ffmpeg libraries that are deprecated there, it might not be a good idea to use the latest ffmpeg code. After ffmpeg, libfann (v1.1b: not needed anymore), easyBMP and libsamplerate are compiled, the sclib general solution file can be opened. It contains the projects sclib, svlib and testbed. If one of the linked libraries is not available, its usage can be switched off by (a) undefining the corresponding constant in SC_Api.h and (b) removing their filename from "Additional Dependencies" key in "Configuration Properties/Linker/Input" of the sclib project properties.

Using testbed, a great number of examples can be explored that are stored as methods in SC_MainTasks.


### Linux

The sclib code has been tested to compile using GCC v4.3.2 and v4.4.3, automake v1.10.1 and v1.11.1, and autoconf v2.61 and v2.63. All "project settings" are to be found in the three files init.sh, configure.ac and Makefile.am. If anything does not work (e.g., a library is not found, or a header), look into these files and maybe change them. A good idea after a system update may be to delete the file aclocal.m4 (it is created again by calling init.sh). 

As a prerequisite, the used libraries have to be compiled and installed. To this end, they all contain (at least in the version provided here) a file named "init.sh" in their base directories. It can be edited to set the correct install path (PREFIX), and then the combination init, make, make install should do everything. The location (path) of the libraries must be stored in the environment variable LIBRARY_PATH and LD_LIBRARY_PATH, and the necessary header files in C_INCLUDE_PATH and CPLUS_INCLUDE_PATH. The ffmpeg libraries maybe copy their headers into a subdirectory /ffmpeg inside the include-path, but they have to be moved to the include-path directly in order to be loeaded correctly.

Then, sclib can be build. This is done by first building the svlib: init, make, make install (additional information can be found in /svlib/src/init.sh. Afterwards, the same is done for the sclib. Which files are to be compiled and which libraries are to be linked can be configured in /sclib/src/Makefile.am. Finally, testbed can be compiled using the same process (without make install) in /testbed/src.