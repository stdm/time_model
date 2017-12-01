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

#ifndef __SC_Corpus_MAC_H__
#define __SC_Corpus_MAC_H__

#include "SC_Corpus.h"
#include "SC_Signal.h"
#include "SC_GroundTruth_MAC.h"
#include <SV_Data.h>

class SCLIB_API SC_Corpus_MAC : public SC_Corpus {
	
  private:

  protected:
	
    bool verbose; //controls some console-output

  public :
 		
    //====================================================================================================================
		// Constructor, destructor
		//====================================================================================================================
    SC_Corpus_MAC(SC_TweakableParameters* pTweak, const char *corpusFileName, bool verbose = true);
		virtual ~SC_Corpus_MAC();

    //====================================================================================================================
		// provides the possibility to load the desired samples even if they are located in different files
    // to due the huge size of the TIMIT corpus it is nor recommended to load it into memory in one single part, but
    // to operate on in scene-wise
    // if the given segment-borders nearly fit the real borders of a file (nearly=within fewer samples than 1 GT-internal
    // frameSize), the segment-borders are aligned to the file-borders (that's the reason for the reference paramerts)
		//====================================================================================================================
    SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, bool unchangeableBoundaries = false);

    //====================================================================================================================
    // Another way to load data, more suitable for this corpus and it's usage: Provide just the start-sample and how many
    // samples of which type should approximately be loaded; the function than returns the endsample of the actually
    // loaded segment and the loaded signal itself, which's length depends on the following rules:
    //  - the loaded segment is homogenious according to the given type
    //  - if there is a little sclib::bit more or less content of the given type than wanted, the load-length is adapted
    //  - always loads files completely, doesn't stop in the midst of one; so it loads one complete file as a minimum
    //  - if prossible, the algorithm prefers to load more than less as wished
    //  - it is guaranteed that everything between segmentStart and segmentEnd gets loaded, so that the resultant signal
    //    and it's boundarys can be used in conjunction with the frameList
    //====================================================================================================================
    SC_Signal* loadSignal(unsigned long int &segmentStart, unsigned long int &segmentEnd, unsigned long int approximateLoadLength, long int type);

    void setVerbose(bool verbose) {this->verbose = verbose; return;}

		//====================================================================================================================
    // Return the fileName of the audio file containg the given startSample
		//====================================================================================================================
		virtual char* getCurrentAudioFileName(unsigned long int startSample) {unsigned long int offset, sampleCount; long int types; return ((SC_GroundTruth_MAC*)this->getGT())->whichFile(startSample, offset, sampleCount, types);}
};

#endif
