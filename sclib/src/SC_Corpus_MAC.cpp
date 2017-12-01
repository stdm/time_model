/**************************************************************************/
/*    Responsibility:																											*/
/*      - encapsulates algorithms zu handle Thilo's MAC (movie audio      */
/*        classification) corpus                                          */
/*      - some assumptions about the corpus:                              */
/*        * each file contains only one audio-type (music/background/     */
/*          noisy speech/pure speech/breath/action) and maybe silence-    */
/*          regions                                                       */
/*        * the groundtruth about the audio-type is extracted from the    */
/*          name of the first subdirectory under the corpus' base path;   */
/*          there is no groundtruth about silence regions!                */
/*        * all files are in the format PCM/16bit/8kHz/mono               */
/*        * typically, all sound-files mentioned in a corpus-file should  */
/*          be ordered by type, so that after a change from directory     */
/*          /MUSIC/ to /BACKGROUND/, no more files of type MUSIC are      */
/*          expected; under this circumstances, the second loadSignal()   */
/*          function works best...                                        */
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 12.06.2006																								*/
/**************************************************************************/

#include "SC_Corpus_MAC.h"
#include "SC_SignalHandler.h"
#include "SC_Signal_WAVE.h"
#include "SC_FeatureHandler.h"
#include "SC_Aux.h"

//====================================================================================================================
//	constructor
//====================================================================================================================
SC_Corpus_MAC::SC_Corpus_MAC(SC_TweakableParameters* pTweak, const char *corpusFileName, bool verbose) : SC_Corpus(pTweak) {
  this->verbose = verbose;
  this->pGT = new SC_GroundTruth_MAC(this->pTweak, corpusFileName, this->verbose);
	this->pSignalHandler->setSignalType(sclib::stWave);
}

//====================================================================================================================
//	destructor
//====================================================================================================================
SC_Corpus_MAC::~SC_Corpus_MAC() {

}

//====================================================================================================================
// provides the possibility to load the desired samples even if they are located in different files
// to due the huge size of the TIMIT corpus it is nor recommended to load it into memory in one single part, but
// to operate on in scene-wise
// if the given segment-borders nearly fit the real borders of a file (nearly=within fewer samples than 1 GT-internal
// frameSize), the segment-borders are aligned to the file-borders (that's the reason for the reference paramerts)
//====================================================================================================================
SC_Signal* SC_Corpus_MAC::loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries) {
  unsigned long int sample, offset = 0, sampleCount = 0, count = 0, start, end, x, internalFrameSize = this->pGT->getInternalFrameSize();
  long int types;
  char *fileName = NULL;
  short *pSamples = NULL, *Buf_L;
	double lastPercentage = 0.0;
  SC_Signal *pSignal = new SC_Signal_WAVE((SC_Signal_WAVE*)this->pGT->getSignalPrototype()), *pTempSignal = NULL;

	if (unchangeableBoundaries == false) {
		//here comes a heuristic that analyses whether the given segment-borders (in samples) fall within less than one GT-internal frameLength
		//of the real file-borders (stored in GT and returned by whichFile() in the offset/sampleCount parameter), but are not exactly as these
		//if this happens because scene-borders (which are ofton used as segment-start- and -end) are only stored frame- and not sample-precise.
		//so, if we find the above condition true, we guess that the user wanted to load a file from it's real beginning or till it's real end 
		//rather than omit the first or last few samples (remember: less than one GT-internal frameSize) and alter the given borders accordingly
		fileName = ((SC_GroundTruth_MAC*)this->pGT)->whichFile(segmentStart + internalFrameSize, offset, sampleCount, types);
		if (abs((long)segmentStart - (long)offset) < (long)internalFrameSize) {
			segmentStart = offset; //move segmentStart to real first sample of file
		}
		MFree_1D(fileName);
		fileName = ((SC_GroundTruth_MAC*)this->pGT)->whichFile(segmentEnd - internalFrameSize, offset, sampleCount, types); //so we are really near the and of the wanted file, not near the beginning of the next one...
		if (abs((long)segmentEnd-(long)offset-(long)sampleCount) < (long)internalFrameSize) {
			segmentEnd = offset + sampleCount - 1; //move the segmentEnd to the real last sample of the file
		}
		MFree_1D(fileName);
	}

  //now we have the real borders, we can allocate memory for the whole buffer
  MArray_1D(pSamples, segmentEnd-segmentStart+1, short, "SC_Corpus_MAC.loadSignal: pSamples");

  if (this->verbose == true) {
    printf("Loading MAC signals between %d and %d (%0.1f min.): ", segmentStart, segmentEnd, this->pGT->getConverter()->sample2ms(segmentEnd-segmentStart+1)/60000.0);
  }

  sample = segmentStart;
  while (sample < segmentEnd) {
		fileName = ((SC_GroundTruth_MAC*)this->pGT)->whichFile(sample, offset, sampleCount, types);
    start = ((long)segmentStart - (long)offset < 0) ? 0 : segmentStart - offset;
    end = (segmentEnd - offset > sampleCount) ? sampleCount - 1 : segmentEnd - offset;

    pTempSignal = this->pSignalHandler->loadSignal(fileName, start, end, this->pTweak->signalHandler.forceSampleRate);

    Buf_L = pTempSignal->GetBuf_L();
    for (x = 0; x < (unsigned long)pTempSignal->GetLen(); x++) {
      pSamples[count+x] = Buf_L[x];  
    }
    count += pTempSignal->GetLen();
    sample += pTempSignal->GetLen();
    
    MFree_1D(fileName);
    MFree_0D(pTempSignal);

    if (this->verbose == true) {
			lastPercentage = sclib::printPercentage(segmentEnd, sample, lastPercentage, 2.5, !start);
    }
  }
 
  pSignal->setBuf_L(pSamples, count);

  if (this->verbose == true) {
    printf("done!\n");
  }

  return pSignal;
}

//====================================================================================================================
// Another way to load data, more suitable for this corpus and it's usage: Provide just the start-sample and how many
// samples of which type should approximately be loaded; the function than returns the endsample of the actually
// loaded segment and the loaded signal itself, which's length depends on the following rules:
//  - the loaded segment is homogenious according to the given type
//  - if there is a little sclib::bit more or less content of the given type than wanted, the load-length is adapted
//  - always loads files completely, doesn't stop in the midst of one; so it loads one complete file as a minimum
//  - if possible, the algorithm prefers to load more than less as wished
//  - it is guaranteed that everything between segmentStart and segmentEnd gets loaded, so that the resultant signal
//    and it's boundarys can be used in conjunction with the frameList
//====================================================================================================================
SC_Signal* SC_Corpus_MAC::loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, unsigned long int approximateLoadLength, long int type) {
  unsigned long int sample, offset = 0, sampleCount = 0, count = 0, x, internalFrameSize = this->pGT->getInternalFrameSize();
  long int types, start, end;
  char *fileName = NULL;
  short *pSamples = NULL, *Buf_L;
  bool finished = false;
	double lastPercentage = 0.0;
  SC_Signal *pSignal = NULL, *pTempSignal = NULL;
  SC_SignalHandler signalHandler(this->pTweak, sclib::stWave);

  //align the startsample with the next fitting one of the wanted type
  end = segmentStart - 1;
  do {
    start = end + 1;
    fileName = ((SC_GroundTruth_MAC*)this->pGT)->whichFile(start + internalFrameSize, offset, sampleCount, types);
    start = offset; //move segmentStart to real first sample of file
    end = start + sampleCount - 1;
    MFree_1D(fileName);
  } while (sclib::bitTest(types, type) != true && sampleCount > 0);
  segmentStart = start;
  
  if (sampleCount > 0) { //otherwise there where no more samples of the wanted type, so nothing to load...
    //decide how much/which portions to load
    end = segmentStart - 1;
    sampleCount = 0;
    while (finished != true) {
      start = end + 1;
      fileName = ((SC_GroundTruth_MAC*)this->pGT)->whichFile(start, offset, sampleCount, types);
      MFree_1D(fileName);
      //we don't need to set start here...
      end = offset + sampleCount - 1;
      if (sclib::bitTest(types, type) == true) {
        count += sampleCount;
        segmentEnd = sclib::min(end, (long)this->pGT->getAudioSampleCount() - 1);
        if (count >= approximateLoadLength || end >= (long)this->pGT->getAudioSampleCount() - 1) {
          finished = true;
        }
      } else {
        finished = true;
      }
    }

    //now we have the real borders, we can allocate memory for the whole buffer
    MArray_1D(pSamples, count, short, "SC_Corpus_MAC.loadSignal: pSamples");

    if (this->verbose == true) {
      printf("Loading MAC signals between %d and %d (%0.1f min.): ", segmentStart, segmentEnd, this->pGT->getConverter()->sample2ms(segmentEnd-segmentStart+1)/60000.0);
    }

    //now really load the samples between the previously determined boundarys
    count = 0;
    sample = segmentStart;
    while (sample < segmentEnd) {
      fileName = ((SC_GroundTruth_MAC*)this->pGT)->whichFile(sample, offset, sampleCount, types);
      start = 0; //start- and end-sample relating to the concrete file: simple here (in contrast to the other loadSignal() function) because we only want complete files, not parts of their content
      end = sampleCount - 1;

      pTempSignal = signalHandler.loadSignal(fileName, start, end, this->pTweak->signalHandler.forceSampleRate);
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
   
    pSignal = new SC_Signal_WAVE((SC_Signal_WAVE*)this->pGT->getSignalPrototype());
    pSignal->setBuf_L(pSamples, count);

    if (this->verbose == true) {
      sclib::printPercentage(1, 1, 0.0, 0.0, false);
      printf("done!\n");
    }
  }

  return pSignal;
}
