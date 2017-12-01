/**************************************************************************/
/*    Derived from:																												*/
/*      - SC_MixtureModel to represent a GMM-UBM as described in				  */
/*				'Speaker Verification Using Adapted Gaussian Mixture Models',   */
/*				D.A.Reynolds, T.F.Quatieri, R.B.Dunn, 2000 (Academic Press)     */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 14.04.2006																								*/
/**************************************************************************/

#ifndef __SC_MixtureModel_GMM_UBM_H__
#define __SC_MixtureModel_GMM_UBM_H__

#include "SC_MixtureModel.h"
#include "SC_MixtureModel_GMM.h"

class SC_MixtureModel_GMM_UBM : public SC_MixtureModel {

  private :

  protected :

 	  //====================================================================================================================
	  // this nested class provides the possibility of having an global top-scoring-mixture-list to speed up successive 
    // model-test-runs on the same feature-set
    //====================================================================================================================
    class SC_MixtureCache {
      private:
        unsigned int topListSize;
        double **topList;
        SV_Data *pLastFeatureset;
    
      public:
        SC_MixtureCache(unsigned int topListSize = 0) {this->topList = NULL; this->pLastFeatureset = NULL; this->topListSize = topListSize;}; //no memory is allocated at construction-time!
        virtual ~SC_MixtureCache() {MFree_2D(this->topList); MFree_0D(this->pLastFeatureset);};

	      //====================================================================================================================
	      // initialize the internal topList holding a list of the tp-scoring mixtures for each frame in the feature-set
	      //====================================================================================================================
        void initTopList(SV_Data *pFeatures, unsigned int topListSize) {
          MFree_2D(this->topList);
          this->topListSize = topListSize;
          MArray_2D(this->topList, (long int)(pFeatures->Row*this->topListSize), 2, double, "SC_MixtureCache.initTopList: topList");
          
          MFree_0D(this->pLastFeatureset);
          this->pLastFeatureset = new SV_Data((SV_Data&)(*pFeatures), true);

          return;
        }

	      //====================================================================================================================
	      // copies the list of listSize topScoring Mixtures for the given frame-index into the internal topList for all frames
        // featureVectorIndex (the frame-index) starts with 0, listSize really means the size, not the index of the last 
        // element in newTopList; it is assumed that the internal topList was initialized with initTopList prior to calling 
        // this method so that it has the right size
        // (the topList per frame is assumed to be arranged in ascending order with the best-scoring mixture being the last 
        //  entry for that frame)
	      //====================================================================================================================
        void copyTopList(unsigned long int featureVectorIndex, double **newTopList) {
          unsigned long int count = 0;

          for (unsigned long int y = featureVectorIndex*this->topListSize; y < (featureVectorIndex+1)*this->topListSize; y++) {
            for (int x = 0; x < 2; x++) {
              this->topList[y][x] = newTopList[count][x];
            }
            count++;
          }

          return;
        }

	      //====================================================================================================================
	      // returns the mixture-index and likelihood of the top-scoring mixture with given rank for the givn frame-index
        // (the topList per frame is arranged in ascending order with the best-scoring mixture being the last entry for that 
        // frame)
	      //====================================================================================================================
        short int getTopMixture(unsigned long int featureVectorIndex, unsigned short int rank, double &likelihood) {
          unsigned long int idx = ((featureVectorIndex+1) * this->topListSize) - rank;
          
          likelihood = this->topList[idx][1];
          return (short int)sclib::round(this->topList[idx][0]);
        }

	      //====================================================================================================================
	      // returns true if the given feature-set matches the one for which the cache stores the topList; returns false 
        // otherwise; the two feature-sets are sayed to match if the header-information (nr. of rows/cols, frame-size and 
        // -step, ID) and a random entry in the feature-matrix match.
	      //====================================================================================================================
        bool cacheHit(SV_Data *pFeatures) {
          bool res = false;
          unsigned long int idx_x, idx_y;

          if (this->pLastFeatureset != NULL &&
              pFeatures->Row == this->pLastFeatureset->Row &&
              pFeatures->Col == this->pLastFeatureset->Col &&
              pFeatures->Hdr.frameSize != this->pLastFeatureset->Hdr.frameSize &&
              pFeatures->Hdr.frameStep != this->pLastFeatureset->Hdr.frameStep &&
              pFeatures->Hdr.ID == this->pLastFeatureset->Hdr.ID) {
						idx_x = (unsigned long int)(floor(sclib::rand(0.0, (double)(pFeatures->Col-1))));
            idx_y = (unsigned long int)(floor(sclib::rand(0.0, (double)(pFeatures->Row-1))));
            if (pFeatures->Mat[idx_y][idx_x]-this->pLastFeatureset->Mat[idx_y][idx_x] < std::numeric_limits<float>::epsilon()*2) {
              res = true;
            }
          }

          return res;
        }
    };

    static SC_MixtureModel_GMM_UBM::SC_MixtureCache topMixtureCache; //hold a global cache for top-scoring mixtures to speed u successive model-tests on the same feature-set with different models having the same UBM
		friend void libDestruct(void); //allows the library destructor (in SC_Lib.cpp) to destruct the topMixtureCache
    double **pTopMixtureList; //a list to remember the index and likelihood-value of the C (see pTweak) top-scoring mixtures of the UBM ((double)idx|likelihood)

    bool insertTopMixture(unsigned short int mixtureNum, double likelihood); //add a mixture/likelihood-pair to the list of top-scoring ones; returns true if it was added (if it was top-scoring), otherwise dosn't add and returns false
    short int getTopMixture(unsigned short int rank, double &likelihood); //returns the mixture-number of the top-scoring mixture with given rank (1..C) or -1 in the case of error

    SC_MixtureModel_GMM *pUBM; //pointer to the UBM from which the parameters are adopted

  public :

	  //====================================================================================================================
	  // constructor/destructor
	  //====================================================================================================================
		SC_MixtureModel_GMM_UBM(SC_TweakableParameters* pTweak, SC_MixtureModel_GMM *pUBM);
    SC_MixtureModel_GMM_UBM(const SC_MixtureModel_GMM_UBM& pParent);
	  ~SC_MixtureModel_GMM_UBM();

    //====================================================================================================================
    // overloaded assignment-operator
    //====================================================================================================================
    SC_MixtureModel_GMM_UBM& operator=(const SC_MixtureModel_GMM_UBM& pParent);

	  //====================================================================================================================
	  // methods to estimate and test a model
	  //====================================================================================================================
	  virtual int TrainModel(SV_Data *TrainData);
		virtual int	TrainModel(SV_Data *pData, unsigned long int segmentsToMerge);
	  virtual SV_Data* TestModel(SV_Data *TestData);
		virtual SV_Data* TestModel(SV_Data *TestData, unsigned long int segmentsToMerge);

    //====================================================================================================================
	  // get & set (new) background model; (wo do not have a background-model in a gmm, so return NULL/do nothing)
    // this is just for compatibility with the functions in SC_Cluster, which have to work also with an GMM-IB!
	  //====================================================================================================================
    virtual SC_Model* getBackground() {return NULL;}
    virtual void setBackground(SC_Model* pBackground) {return;}
    void setUBM(SC_MixtureModel_GMM *pNewUBM) {this->pUBM = pNewUBM; return;}
    SC_MixtureModel* getUBM(void) {return this->pUBM;}

    virtual double getWeightLimit(void) {return 0.0;} //no weight-limiting in the GMM-UBM

    //====================================================================================================================
		// Combine 2 Models by adding the mixtures of this and the second one and return a new model
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond);
    virtual SC_Model* combineModels(SC_Model* pFirst, SC_Model* pSecond, bool keepFirstsNext = false);
		SC_MixtureModel_GMM_UBM* combineModels(SC_MixtureModel* pSecond);
    SC_MixtureModel_GMM_UBM* combineModels(SC_MixtureModel* pFirst, SC_MixtureModel* pSecond, bool keepFirstsNext = false);

		//====================================================================================================================
		// Combine 2 Models by training a new one with the combination of the training-data of this and pSecond.
		//====================================================================================================================
		virtual SC_Model* combineModels(SC_Model* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);
    SC_MixtureModel_GMM_UBM* combineModels(SC_MixtureModel* pSecond, SV_Data* pSpeechFrames, unsigned long int segmentsToMerge, SC_Model* pBackgroundModels);

	  //====================================================================================================================
	  // for computing BIC etc.
	  //====================================================================================================================
    virtual unsigned int getFreeParameterCount(void) {return this->mixtureCount*this->dim*2 + this->mixtureCount;} //mean & variance & weight per dimension and mixture

	  //====================================================================================================================
	  // generally, all child classes return averaged scores (i.e. divided by the number of test patterns); if not, this 
		// method needs to return false;
	  //====================================================================================================================
		virtual bool scoreIsAverage(void) {return (this->pTweak->mixtureModelGmmubm.scoringMethod == sclib::scoringGMM);}
};

#endif
