/**************************************************************************/
/*	This is the general include file for the SC_Lib library, which should */
/*  be included by all it's users instead of the individual SC_*.h files  */
/*                                                                        */
/*  It also manages (in the correspondig .cpp file) dll-management for    */
/*  windows and provides the "main" function for calls via JNI (probably  */
/*  only important for videana)																						*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 05.03.2004																								*/
/**************************************************************************/

#ifndef __SC_Lib_H__
#define __SC_Lib_H__

//====================================================================================================================
//  All headers necessary in external programs (not all headers in the lib...)
//====================================================================================================================
#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_Cluster.h"
#include "SC_Clusterer.h"
#include "SC_Classifier.h"
#include "SC_Classifier_ML.h"
#include "SC_Classifier_SVM.h"
#include "SC_ClassifierTree.h"
#include "SC_Centroid.h"
#include "SC_Centroid_Point.h"
#include "SC_Corpus.h"
#include "SC_Corpus_MPEG7.h"
#include "SC_Corpus_MAC.h"
#include "SC_Corpus_TIMIT.h"
#include "SC_Corpus_SCiVo.h"
#include "SC_Corpus_Wesley.h"
#include "SC_DistanceMeasures.h"
#include "SC_Enhancement.h"
#include "SC_Feature_FbE.h"
#include "SC_Feature_Formant.h"
#include "SC_Feature_LPC.h"
#include "SC_Feature_LPCresidual.h"
#include "SC_Feature_MFCC.h"
#include "SC_Feature_Pitch.h"
#include "SC_FeatureHandler.h"
#include "SC_GroundTruth.h"
#include "SC_GroundTruth_MAC.h"
#include "SC_GroundTruth_MPEG7.h"
#include "SC_GroundTruth_SCiVo.h"
#include "SC_GroundTruth_TIMIT.h"
#include "SC_GroundTruth_Wesley.h"
#include "SC_HHT.h"
#include "SC_Ini.h"
#include "SC_MainTasks.h"
#include "SC_Matlab.h"
#include "SC_MixtureModel_bGMM.h"
#include "SC_MixtureModel_GMM.h"
#include "SC_MixtureModel_GMM_UBM.h"
#include "SC_Model_HMM.h"
#include "SC_Model_FullGauss.h"
#include "SC_Model_Pareto.h"
#include "SC_Model_SVM.h"
#include "SC_Model_Time.h"
#include "SC_ModelHandler.h"
#include "SC_Partition.h"
#include "SC_Signal.h"
#include "SC_Score.h"
#include "SC_Score_AudioTypeClassification.h"
#include "SC_Score_ChangeDetection.h"
#include "SC_SpeakerScore_Clustering.h"
#include "SC_SpeakerScore_Identification.h"
#include "SC_SpeakerScore_Classification.h"
#include "SC_Score_VUv.h"
#include "SC_Segmentation_AudioType_LZL.h"
#include "SC_Segmentation_Changes_KBK.h"
#include "SC_Segmentation_Changes_LZW.h"
#include "SC_Segmentation_Silence_LZL.h"
#include "SC_Segmentation_Silence_LNK.h"
#include "SC_SegmentationHandler.h"
#include "SC_SignalHandler.h"
#include "SC_Signature.h"
#include "SC_SpeakerClusterer.h"
#include "SC_SpeakerIdentificator.h"
#include "SC_Synthesis.h"
#include "SC_TweakableParameters.h"
#include "SC_Timer.h"

//====================================================================================================================
//  "Main" function to be called by Videana via JNI (does in principle what SCiVo does...)
//  Returned are lists with the following information per video-frame:
//    - segmentation-results with bitflags giving the labels belonging to each frame (as parameter)
//    - speaker-ids (as parameter)
//    - probabilities of segmentation results in columns 0-31 and of speaker-id in column 32 (as parameter), if not 
//      switched off in the tweakable parameters
//    - probabilityCols is the number of columns in the probability-matrix
//    - length of each of the three lists (=nr. of video-frames -1)
//  the shotList contains videoFrames-numbers at which shot transition occur
//====================================================================================================================
SCLIB_API_C unsigned long int audioSegmentation(const char *fileName, SC_TweakableParameters *pTweak, double videoFrameRate, long int* segmentationResults, long int* speakerIDs, long int* probabilities, int* probabilityCols, int *shotList = NULL, int shotListLength = 0);

//====================================================================================================================
//  Returns the value of the "progress" global variable, the progress (in percent) of the current stage (quite coarse: 
//  segmentation or speaker-id) of the audioSegmentation() function
//====================================================================================================================
SCLIB_API_C double getProgress(void);

//====================================================================================================================
//  Returns the value of the "progressStage" global variable, the name of the stage of the audioSegmentation() 
//  function; memory for this variable is allocated and freed inside the library, so don't touch it!
//====================================================================================================================
SCLIB_API_C char* getProcessingStage(void);

//====================================================================================================================
//  Called from outside to tell the audioSegmentation() function to immediately return
//====================================================================================================================
SCLIB_API_C void doAbort(void);

//====================================================================================================================
//  A bunch of memory-deallocation methods to be called from a (java) application and an associated enumeration to 
//  tell the methods which type the given object has; at the moment, only those types currently passed to videana are
//  recognized in the enumeration; TODO: enhance it with all exported types
//====================================================================================================================
enum SC_Types {scLong=0, scDouble=1, scTweakableParameters=2, svData=3, scModelHandler=4, scModel=5, scCentroid=6, scSignature=7, scDistanceMeasure=8}; //we want those elements to have the same numbers even when we later on enhance the list (possibly changing the order), that's why we give the corresponding numbers explicitly
SCLIB_API_C void deleteScalar(void *scalar, SC_Types baseType);
SCLIB_API_C void deleteArray(void *array, SC_Types baseType);
SCLIB_API_C void deleteMatrix(void **matrix, SC_Types baseType);

//====================================================================================================================
//  Methods to work with array of pointers (to objects) from within java applications
//====================================================================================================================
SCLIB_API_C void* constructPointerArray(SC_Types baseType, int dim);
SCLIB_API_C void* getPointerArrayElement(void **array, int idx);
SCLIB_API_C int setPointerArrayElement(void **array, int idx, void *newPointer);
SCLIB_API_C void deletePointerArray(void **array, SC_Types baseType);

//====================================================================================================================
//  Get instant printf()'s from e.g. Java
//====================================================================================================================
SCLIB_API_C void getInstantPrintfs(void);
typedef int SC_Bool; //java dispatch.dll cannot handle bools, so we give ints instead (==0 <=> false, ==1 <=> true) and rename it in order to distinguish it from ordinary ints

//====================================================================================================================
//  Interface to the SC_TweakableParameters class
//====================================================================================================================
SCLIB_API_C SC_TweakableParameters* scTweakableParameters_construct(const char *fileName = "", SC_Bool verbose = 1);
SCLIB_API_C SC_Bool scTweakableParameters_setByName(SC_TweakableParameters *pTweak, const char *parameterName, const char *value);
SCLIB_API_C SC_TweakableParameters* createTweak(const char *fileName = ""); //deprectated, just exists for compatibility reasons
SCLIB_API_C SC_Bool setParameterByName(SC_TweakableParameters *pTweak, char *parameterName, char *value); //deprectated, just exists for compatibility reasons

//====================================================================================================================
//  Interface to the SV_Data class
//====================================================================================================================
SCLIB_API_C SV_Data* svData_construct(int rows = 0, int cols = 0);
SCLIB_API_C SC_Bool svData_setCol(SV_Data *pData, int col);
SCLIB_API_C int svData_getCol(SV_Data *pData);
SCLIB_API_C SC_Bool svData_setRow(SV_Data *pData, int row);
SCLIB_API_C int svData_getRow(SV_Data *pData);
SCLIB_API_C SC_Bool svData_alloc(SV_Data *pData);
SCLIB_API_C SC_Bool svData_setMat(SV_Data *pData, int row, int col, float value);
SCLIB_API_C float svData_getMat(SV_Data *pData, int row, int col);
SCLIB_API_C SC_Bool svData_setNext(SV_Data *pData, SV_Data *pNext);
SCLIB_API_C SV_Data* svData_getNext(SV_Data *pData);
SCLIB_API_C SV_Data* svData_mergeData(SV_Data *pData, int maxSegments = 0);

//====================================================================================================================
//  Interface to the SC_ModelHandler class
//====================================================================================================================
SCLIB_API_C SC_ModelHandler* scModelHandler_construct(SC_TweakableParameters *pTweak, SC_Bool verbose = 1);
SCLIB_API_C int scModelHandler_guessModelOrder(SC_ModelHandler *pModelHandler, SV_Data *pFeatures, SC_Model *pBackground = NULL, int modelType = sclib::mtGMM_new, int minOrder = 1, int maxOrder = 2048, int segmentsToMerge = 1);
SCLIB_API_C SC_Model* scModelHandler_buildModel(SC_ModelHandler *pModelHandler, SV_Data *pFeatures, SC_Model *pBackground = NULL, int modelOrder = 16, int modelType = sclib::mtGMM_new, int segmentsToMerge = 1);
SCLIB_API_C double scModelHandler_testModel(SC_ModelHandler *pModelHandler, SC_Model *pModel, SV_Data *pFeatures, int segmentsToMerge = 1, SC_Bool forceNoNormalization = 0, SC_Model *pBackground = NULL);
SCLIB_API_C int scModelHandler_saveModel(SC_ModelHandler *pModelHandler, const char *fileName, SC_Model *pModel, SC_Bool useDebugDir = 1);
SCLIB_API_C SC_Model* scModelHandler_loadModel(SC_ModelHandler *pModelHandler, const char *fileName, int modelType = sclib::mtGMM_new);
SCLIB_API_C SC_Bool scModelHandler_setTweak(SC_ModelHandler *pModelHandler, SC_TweakableParameters *pTweak);

//====================================================================================================================
//  Interface to SC_Centroid_Point class
//====================================================================================================================
SCLIB_API_C SC_Centroid_Point* scCentroidPoint_construct(SC_TweakableParameters *pTweak = NULL, int dim = 0, double *coordinate = NULL, SC_Bool justLink = 1);
SCLIB_API_C double* scCentroidPoint_getCoordinate(SC_Centroid_Point* pCentroid);
SCLIB_API_C int scCentroidPoint_getDim(SC_Centroid_Point* pCentroid);
SCLIB_API_C int scCentroidPoint_setCoordinate(SC_Centroid_Point* pCentroid, double *newCoordinate);
SCLIB_API_C int scCentroidPoint_setDim(SC_Centroid_Point* pCentroid, int newDim);
SCLIB_API_C double scCentroidPoint_getDistance(SC_Centroid_Point* pCentroid, SC_Centroid *secondCentroid);

//====================================================================================================================
//  Interface to SC_Signature class
//====================================================================================================================
SCLIB_API_C SC_Signature* scSignature_construct(SC_Centroid **centroids = NULL, double *weights = NULL, int n = 0, SC_Bool justLinkWeights = 1, double smallestUnnormalizedWeight = 1.0);

//====================================================================================================================
//  Interface to SC_DistanceMeasures class
//====================================================================================================================
SCLIB_API_C SC_DistanceMeasures* sc_DistanceMeasures_construct(SC_TweakableParameters* pTweak, SC_MatrixFunctions* pMatrixFunc = NULL, SC_Bool verbose = 1);
SCLIB_API_C double scDistanceMeasures_EMD(SC_Signature* pSignature1, SC_Signature* pSignature2, SC_TweakableParameters *pTweak);

#endif
