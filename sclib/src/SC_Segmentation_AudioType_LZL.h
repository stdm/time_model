/**************************************************************************/
/*    Responsibility:																											*/
/*      - Implements the audio-segmenter published in "Content-based Audio*/
/*        Classification and Segmentation by Using Support Vector         */
/*        Machines", Lu/Zhang/Li 2003                                     */
/*      - Provides possibilities to train & use the classifier            */
/*      - Capable to separate speech/non-pure speech/music/noise from     */
/*        audiofeatures of a non-silent signal                            */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 10.03.2006																								*/
/**************************************************************************/

#ifndef __SC_Segmentation_AudioType_LZL_H__
#define __SC_Segmentation_AudioType_LZL_H__

#include "SC_Segmentation_AudioType.h"
#include "SC_TweakableParameters.h"
#include "SC_GroundTruth.h"
#include "SC_Api.h"

class SCLIB_API SC_Segmentation_AudioType_LZL : public SC_Segmentation_AudioType {
  
  private :

    //====================================================================================================================
    //  To construct a linked list of labels easily
    //====================================================================================================================
    class SC_Labels {
      protected:
        int *labels;
        unsigned long int count;

      public:
        SC_Labels *Next;

        SC_Labels(int *labels = NULL, unsigned long int count = 0) { //the given labels get copied!
          this->Next = NULL;
          this->labels = NULL;
          this->count = 0;
          setLabels(labels, count);
        }

        SC_Labels(int label, unsigned long int count) { //a list is created with "count" entrys of the value "label"
          this->Next = NULL;
          this->labels = NULL;
          this->count = 0;
          setLabels(label, count);
        }

        virtual ~SC_Labels() {
          MFree_1D(this->labels);
        }

        void setLabels(int *labels, unsigned long int count) {
          MFree_1D(this->labels);
          this->count = count;
          if (this->count > 0) {
            MArray_1D(this->labels, count, int, "SC_Labels.SC_Labels: labels");
            for (unsigned long int x = 0; x < count; x++) {
              this->labels[x] = labels[x];
            }
          }
          return;
        }

        void setLabels(int label, unsigned long int count) {
          MFree_1D(this->labels);
          this->count = count;
          if (this->count > 0) {
            MArray_1D(this->labels, count, int, "SC_Labels.SC_Labels: labels");
            for (unsigned long int x = 0; x < count; x++) {
              this->labels[x] = label;
            }
          }
          return;
        }

        SC_Labels* addList(int label, unsigned long int count) { //add a new labels-array to the list (Next-pointer) and return Next-pointer
          this->Next = new SC_Labels(label, count);
          return this->Next;
        }

        int* getLabels(void) {
          return this->labels;
        }

        unsigned long int getCount(void) {
          return this->count;
        }

        int Valid(void) {
          return 1;
        }

        int* merge(unsigned long int &count) { //return a single array with all labels in the linked list of label-arrays
          int *completeLabels = NULL;
          long int x, overallCount = 0;
          SC_Labels *pHook = this;
          
          while (pHook != NULL) { //count overall number of labels 
            overallCount += pHook->getCount();
            pHook = pHook->Next;
          }
          count = overallCount;

          if (overallCount > 0) { //construct concatenated list
            MArray_1D(completeLabels, overallCount, int, "SC_Labels.merge: completeLabels");
            pHook = this;
            overallCount = 0;
            while (pHook != NULL) { //count overall number of labels 
              for (x = 0; x < (long int)(pHook->getCount()); x++) {
                completeLabels[overallCount + x] = pHook->labels[x];
              }
              overallCount += pHook->getCount();
              pHook = pHook->Next;
            }
          }

          return completeLabels;
        }
    };

		//TODO: testtesttest
		//int actionThreshClassifyer(SC_GroundTruth *pGT, SV_Data *pMFCC, unsigned long int segmentStart, unsigned long int segmentEnd);

  protected :

    //====================================================================================================================
    //	Convert the per-frame-features to per-subclip-features and return a new SV_Data object containing:
    //   - The following values for the columns sztarting with startCol and ending with endCol (or the last available col,
    //     if endCol=0)
    //   - The mean (if divByLengthMinusOne=false; otherwise mean*(T/(T-1)) with T number of frames in a subclip)
    //   - The standard deviation of the features in the subclip, if addSd=true
    //  Subclip-length has to be given in samples
    //  At the end there may be some frames missing that didn't fully fit into a subclip; we don't car about them if 
		//  createUnderpopulatedFrames==false
    //====================================================================================================================
    SV_Data* aggregatePerSubClip(SV_Data *pFeatures, unsigned long int subClipLength, unsigned int startCol = 0, unsigned int endCol = 0, bool addSd = false, bool divByLengthMinusOne = false, bool createUnderpopulatedFrames = false);

    //====================================================================================================================
    //	This is the realy classification/segmentation algorithm as described in the paper "Content-based Audio 
    //  Classification and Segmentation by Using Support Vector Machines", Lu/Zhang/Li 2003. Though most if it's work is
    //  done in the feature-extraction/feature-handler/classifier-classes, it is rather short & simple
    //  Because the classifier is trained on signals of a specific sampleRate and the testdata needs to look alike, this
    //  Algorithm possibly receives features which have a different SR as stated in the groundTruth (but match the 
    //  classifier). This is compensated for here (and only here throughout the lib except feature extraction, at the 
    //  moment).
    //====================================================================================================================
    int lzlAlgorithm(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
    //	This is the real training algorithm as described in the paper "Content-based Audio Classification and Segmentation 
    //  by Using Support Vector Machines", Lu/Zhang/Li 2003. 
    //  The feature-vectors should be classified into silence/non-silence beforehand because this information is used 
    //  inside
    //====================================================================================================================
    int trainLzlAlgorithm(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
    //	Check if the provided features are suitable for the purpose of this class
    //====================================================================================================================
    bool checkFeatures(SV_Data **pFeatures);

    SV_Data *pTrainingData[32]; //to temporarily store the training-data (one audio-type per cell), so it can be loaded in small parts

    bool cacheFeatures; //if true, the features are saved to disk during partially loading and the loaded at once -> saves memory because no merging needed
    char *cacheFilePrefix; //prefix for the caching files; the suffix determines the contained autio-type, i.e. the numbers from 1-32

    //====================================================================================================================
    //	Loads all features from the disk-cache in ready-to-use form (no more merging needed, all done on disk) togehter
    //  with the respective labels; returns true on success, false otherwise
    //====================================================================================================================
    bool loadFromCache(long int *neededTypes, int typeCount, SV_Data* &pFeatures, int* &labels);

  public :

    SC_Segmentation_AudioType_LZL(SC_TweakableParameters* pTweak, bool cacheFeatures2Disk = true);
    virtual ~SC_Segmentation_AudioType_LZL();

    //====================================================================================================================
		// If a classification-algorithm needs training, this can be handled using this function; otherwise it doesn't need
    // implementation
    // The features needed for training are meant to reside in an array of feature-containers as returned by 
    // SC_FeatureHandler.extractFeatures(), where theire respective labels are properly stored in the groundtruth
		//====================================================================================================================
    virtual int trainClassifier(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
		// classifiy the given audio-segment according to the underlying audio-type (speech/noise/music/...)
    // pFeatures must be an array of feature-sets as returned by the SC_FeatureHandler->extractFeatures() method (with 
    // the log of the feature-set constants SCLIB_FEATURE_* as indices into the array)
		//====================================================================================================================
    virtual int classifyAudioType(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures);

    //====================================================================================================================
    //  If the training-database is too big to hold in memory completely (not forgetting the number of copys for 
    //  copying-together, aggregation, normalization, and inside the classifier), this function can be called several 
    //  times prior to call the trainClassifier(): Each time, a different subpart, representing one needed audio-type,
    //  can be provided (even the same audioType can be provided in parts); it's feature-vectors get aggregated per 
    //  subclip and afterwards stored for training, so the original features can be destroyed. After providing 
    //  feature-sets for all relevant audio-types, the training-algorithm can be called.
    //  pFeatures is meant to be an array as returned by SC_FeatureHandler.extractFeatures(), and must contain frames 
    //  of the specified (single, not or-concatenated) audio-type (but can also contain others, though this gives no 
    //  sense...)
    //====================================================================================================================
    int partiallyLoadTrainingData(SC_GroundTruth *pGT, unsigned long int segmentStart, unsigned long int segmentEnd, SV_Data **pFeatures, unsigned long int audioType);

		//====================================================================================================================
		// Returns the width [in ms!!!] of the region in which the found segment-boundaries may lie for a given exact 
		// position; this uncertainty region is due to frame-based anlysis and internal windowsizes, e.g.
		// Only the factors due to this specific algorithm are taken into account (factors regarding thr ground-truth class
		// are handled therein)
		//====================================================================================================================
		virtual unsigned long int getUncertaintyRegionWidth(void);

		//====================================================================================================================
		// Returns a or-concatenated list of SCLIB_FEATURE_* constants for all feature-types used by this algorithm
		//====================================================================================================================
		virtual unsigned long int getUsedFeatures(void) {return sclib::featureMFCC|sclib::featureZCR|sclib::featureSubbandPower|sclib::featureBrightnessBandwidth|sclib::featureSpectrumFlux|sclib::featureBandPeriodicity|sclib::featureNFR;}

		//====================================================================================================================
		// Returns linked list of feature-parameter objects in case the actual algorithm needs other than the standard 
		// parameters
		//====================================================================================================================
		virtual SC_TweakableParameters::SC_FeaturePar* getSpecialFeatureParameters(void);

		//====================================================================================================================
		// Returns the algorithm's name
		//====================================================================================================================
		virtual const char* getName(void) {return "Segmentation_AudioType_LZL";};
};

#endif
