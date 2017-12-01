/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle the TIMIT corpus              */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 09.04.2006																								*/
/**************************************************************************/

#include "SC_Corpus_TIMIT.h"
#include "SC_SignalHandler.h"
#include "SC_Signal_NIST.h"
#include "SC_FeatureHandler.h"
#include "SC_Aux.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Corpus_TIMIT::SC_Corpus_TIMIT(SC_TweakableParameters* pTweak, const char *corpusFileName, bool verbose, bool setSceneBoundaries, bool setAllSpeakerBoundaries, unsigned int sbSkipCount, bool sceneBoundariesOnlyAtSpeakerBoundaries) : SC_Corpus(pTweak) {
	this->verbose = verbose;
  this->pGT = new SC_GroundTruth_TIMIT(this->pTweak, corpusFileName, this->verbose, (setSceneBoundaries && !sceneBoundariesOnlyAtSpeakerBoundaries), setAllSpeakerBoundaries, sbSkipCount, (setSceneBoundaries && sceneBoundariesOnlyAtSpeakerBoundaries));
	this->pSignalHandler->setSignalType(((sclib::like(corpusFileName, "%conTIMIT%")==true) ? sclib::stWave : sclib::stNIST)); //conTIMIT data has RIFF/Wave files, TMIT is NIST format
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Corpus_TIMIT::~SC_Corpus_TIMIT() {

}

//====================================================================================================================
// provides the possibility to load the desired samples even if they are located in different files
// due to the huge size of the TIMIT corpus it is not recommended to load it into memory in one single part, but
// to operate on it scene-wise ("scenes" are here "files").
// if the given segment-borders nearly fit the real borders of a file (nearly=within fewer samples than 1 GT-internal
// frameSize), the segment-borders are aligned to the file-borders (that's the reason for the reference parameters)
//====================================================================================================================
SC_Signal* SC_Corpus_TIMIT::loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries) {
  unsigned long int sample, offset = 0, sampleCount = 0, count = 0, start, end, x, internalFrameSize = this->pGT->getInternalFrameSize();
  char *fileName = NULL;
  short *pSamples = NULL, *Buf_L;
	double lastPercentage = 0.0;
  SC_Signal *pSignal = new SC_Signal_NIST((SC_Signal_NIST*)this->pGT->getSignalPrototype()), *pTempSignal = NULL;

	if (unchangeableBoundaries == false) {
		//here comes a heuristic that analyses whether the given segment-borders (in samples) fall within less than one GT-internal frameLength
		//of the real file-borders (stored in GT and returned by whichFile() in the offset/sampleCount parameter), but are not exactly as these
		//if this happens because scene-borders (which are ofton used as segment-start- and -end) are only stored frame- and not sample-precise.
		//so, if we find the above condition true, we guess that the user wanted to load a file from it's real beginning or till it's real end 
		//rather than omit the first or last few samples (remember: less than one GT-internal frameSize) and alter the given borders accordingly
		fileName = ((SC_GroundTruth_TIMIT*)this->pGT)->whichFile(segmentStart + internalFrameSize, offset, sampleCount);
		if (abs((long)segmentStart - (long)offset) < (long)internalFrameSize) {
			segmentStart = offset; //move segmentStart to real first sample of file
		}
		MFree_1D(fileName);
		fileName = ((SC_GroundTruth_TIMIT*)this->pGT)->whichFile(segmentEnd - internalFrameSize, offset, sampleCount); //so we are really near the end of the wanted file, not near the beginning of the next one...
		if (abs((long)segmentEnd-(long)offset-(long)sampleCount) < (long)internalFrameSize) {
			segmentEnd = offset + sampleCount - 1; //move the segmentEnd to the real last sample of the file
		}
		MFree_1D(fileName);
	}

  //now we have the real borders, we can allocate memory for the whole buffer
  MArray_1D(pSamples, segmentEnd-segmentStart+1, short, "SC_Corpus_TIMIT.loadSignal: pSamples");

  if (this->verbose == true) {
    printf("Loading TIMIT signals between %d and %d (%0.1f min.): ", segmentStart, segmentEnd, this->pGT->getConverter()->sample2ms(segmentEnd-segmentStart+1)/60000.0);
  }

  sample = segmentStart;
  while (sample < segmentEnd) {
    fileName = ((SC_GroundTruth_TIMIT*)this->pGT)->whichFile(sample, offset, sampleCount);
    start = ((long)segmentStart - (long)offset < 0) ? 0 : segmentStart - offset;
    end = (segmentEnd - offset > sampleCount) ? sampleCount - 1 : segmentEnd - offset;

    pTempSignal = this->pSignalHandler->loadSignal(fileName, start, end, this->pTweak->signalHandler.forceSampleRate);

    Buf_L = pTempSignal->GetBuf_L();
    for (x = 0; x < (unsigned long)pTempSignal->GetLen(); x++) {
      pSamples[count+x] = Buf_L[x];  
    }
    count += pTempSignal->GetLen();
    sample += pTempSignal->GetLen();
    
    if (this->verbose == true) {
			lastPercentage = sclib::printPercentage(segmentEnd, sample, lastPercentage, 2.5, !(count-pTempSignal->GetLen()));
    }

    MFree_1D(fileName);
    MFree_0D(pTempSignal);
  }
 
  pSignal->setBuf_L(pSamples, count);
	pSignal->setSampleCount(this->getGT()->getAudioSampleCount());

  if (this->verbose == true) {
		sclib::printPercentage(1, 1, 0.0, 0.0, false);
    printf("done!\n");
  }

  return pSignal;
}

/*
//====================================================================================================================
// don't use this method!!! reason:
// at the end of each file there are some frames missing because they don't fit completely into the loaded signal
// part (only the samples of one file are loaded at one moment). as more and more files are loaded and features are 
// extracted, the frame-numbers get more and more out of sync with the way frames are mapped to samples (and vice 
// versa) in the ground-truth object.
// so better use scene-based analysis as in SCiVo or the TIMIT data (scene-boundarys are set at each file's first 
// sample) for recognition/dectection involving the groundtruth, and then maybe copy samples together (using the
// copyFramesTogether()-method of the ground-truth object) to build bigger models etc.
// nevertheless, this method's implementation can serve as an example for how to concatenate the feature-sets after
// all ground-truth work has been done.
//====================================================================================================================
// there might be a memory problem when first loading a huge part of the corpus by calling loadSignal and the doing
// feature extraction (using the SC_FeatureHandler class) on this large sample-array (feature extraction even 
// involves min. 1 copy of the comlete samples...).
// so this method provides a way of specifiying which features should be extracted on which part of the signal 
// (can be located in different files as above), and the features get extracted file-by-file and are then concatenated
//====================================================================================================================
SV_Data** SC_Corpus_TIMIT::extractFeatures(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int featureTypes) {
  unsigned long int sample = segmentStart, offset = 0, sampleCount = 0, start, end, x, featureCount;
  char *fileName = NULL;
  double currentPercent, percent = 0.0;
  SV_Data **pTempFeatures, **pFeaturesStart = NULL, **pFeaturesEnd = NULL, **pFeatures = NULL;
  SC_Signal *pTempSignal = NULL;
  SC_SignalHandler signalHandler(this->pTweak, sclib::stNIST);
  SC_FeatureHandler featureHandler(this->pTweak);

  featureCount = featureHandler.getFeatureCount();
  MArray_1D(pFeaturesEnd, featureCount, SV_Data*, "SC_Corpus_TIMIT.extractFeatures: pFeaturesEnd");

  if (this->verbose == true) {
    printf("Extracting TIMIT features between %d and %d (%0.1f min.): ", segmentStart, segmentEnd, this->pGT->sample2ms(segmentEnd-segmentStart+1)/60000.0);
  }

  while (sample < segmentEnd) {
    fileName = ((SC_GroundTruth_TIMIT*)this->pGT)->whichFile(sample, offset, sampleCount);
    start = ((long)segmentStart - (long)offset < 0) ? 0 : segmentStart - offset;
    end = (segmentEnd - offset >= sampleCount) ? sampleCount - 1 : segmentEnd - offset;

    pTempSignal = signalHandler.loadSignal(fileName, start, end, this->pTweak->signalHandler.forceSampleRate);
    MFree_1D(fileName);

    pTempFeatures = featureHandler.extractFeatures(this->pGT, pTempSignal, offset+start, offset+end, featureTypes);
    MFree_0D(pTempSignal);

    if (sample == segmentStart) { //first iteration, mark start of linked list
      pFeaturesStart = pTempFeatures;
      for (x = 0; x < featureCount; x++) {
        pFeaturesEnd[x] = pFeaturesStart[x];
      }
    } else {
      for (x = 0; x < featureCount; x++) { //construct linked list, move pointer to end of list
        if (pTempFeatures[x] != NULL) { //there might be empty slots in the array if not all features where selected for extraction (the normal case)
          pFeaturesEnd[x]->Next = pTempFeatures[x];
          pFeaturesEnd[x] = pFeaturesEnd[x]->Next;
        }
      }
      MFree_1D(pTempFeatures);
    }

    if (this->verbose == true) {
      printf(".");
      currentPercent = ((double)sample/(double)segmentEnd) * 100.0;
      if (currentPercent > percent) {
        printf("%02.1f%%", currentPercent);
        percent = max(currentPercent + 2.5, percent + 2.5);
      }
    }

    sample += end - start + 1;
  }

  MFree_1D(pFeaturesEnd);
  MArray_1D(pFeatures, featureCount, SV_Data*, "SC_Corpus_TIMIT.extractFeatures: pFeatures");
  for (x = 0; x < featureCount; x++) {
    if (pFeaturesStart[x] != NULL) {
      pFeatures[x] = pFeaturesStart[x]->MergeData();
    } else {
      pFeatures[x] = NULL;
    }
  }
  MFree_1D(pFeaturesStart);

  if (this->verbose == true) {
    printf("done!\n");
  }

  return pFeatures;
}
*/
