/**************************************************************************/
/*    Responsibility:																											*/
/*		  - organizes file-io (audio-stream, scene-list, segement-list)			*/
/*			- organizes timeline in terms of video-frames											*/
/*			- maps milliseconds and audio-samples to video-frames							*/
/*			- privides access to the audio-segments														*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 16.02.2004																								*/
/**************************************************************************/

#ifndef __SC_GroundTruth_H__
#define __SC_GroundTruth_H__

#include <iostream>
#include "SC_Aux.h"
#include "SC_Api.h"
#include "SC_TweakableParameters.h"
#include "SC_Signal.h"
#include "SC_Conversion.h"
#include <SV_Lib.h>

//====================================================================================================================
// ATTENTION: Basically, this class provides read/write access to a data-structure called the "frameList", which holds
// some information about small portions of the analyzed audio-file called herein a "frame". The algorithms providing
// the or relying on this information also process not the samples of the audiofile directly, but operate on frames, 
// but the length and progression/step of these frames has not necessarily to coincide with the frames in the frame-
// list. 
//
// For this reason, *all values* describing positions in the audiofile *given to or returned by methods* of this 
// class are *sample-based*. 
//
// They are internally converted to the FLI (frame list index, index into the frameList-array) measure and can 
// externally be converted to the current audioframe-measure using methods provided by the SC_Conversion class
//====================================================================================================================
class SCLIB_API SC_GroundTruth {

  private:

	protected :

    class SC_SpeakerMapping {
      public:
				SC_SpeakerMapping(bool correct = false, long int groundTruthID = sclib::noSpeaker, long int hypothesizedID = sclib::noSpeaker, SC_SpeakerMapping *Next = NULL) {
					this->correct = correct;
					this->groundTruthID = groundTruthID;
					this->hypothesizedID = hypothesizedID;
					this->Next = Next;
				}
        bool correct; //true, if the mapping is correct (not only fitting), i.e. if the hypothesihzed speaker corresponds to the cluster containing the most speech of this groundtruth speaker
        long int groundTruthID;
        long int hypothesizedID;
        SC_SpeakerMapping *Next;
        int Valid(void) {return 1;}
    };

    SC_Signal *pSignalPrototype; //holds parameter-settings of originally opened file
    
    //====================================================================================================================
		// The framelist: 
    //  - The 0. component holds an entry for every audioFrame describing it's content with the bitflag-constants 
    //    defined above, as stated in the ground-truth. 
    //  - The 1. component holds the same, but as hypothesized by the algorithms of this lib. 
    //  - The 2. component holds the ground-truth Speaker-ID for this frame (if available)
    //  - The 3. component holds the hypothesized Speaker-ID as computed by the classification algorithms
		//====================================================================================================================
		long int** frameList; //must be signed!!!

    //====================================================================================================================
		// Stores information about how likely the results in the hypo-col. of the frameList are. for each bitFlag in the
		// frameList there is a col in the probabilityList including a (log-)likelihood score derived by the algorithm which 
		// set that flag or 1.0 as a standard value; there is an additional col holding the probability-information for the
		// speaker-decision
		//====================================================================================================================
		double** probabilityList;

    //====================================================================================================================
		// Holds pointers to the speaker-names given in the ground-truth as comments behind speech-segments; their index is 
    // the corresponding ground-truth Speaker-ID as in the second component of the frameList;
    // This array is filled by insertSpeakerNames() and queried by getSpeakerIDfromName()
		//====================================================================================================================
		char** speakerNames;

    SC_TweakableParameters *pTweak;
    SC_GroundTruth::SC_SpeakerMapping *pSpeakerMapping;

		//====================================================================================================================
		// Some classdata to remember important parameters of the signal, the video and the feature-vectors
		//====================================================================================================================
		char*								audioFileName;
		unsigned long int   internalFrameSize;  //in samples
    unsigned long int   internalFrameCount; 
    unsigned long int   audioSampleRate;    //in Hz
    unsigned long int   audioSampleCount;

		//====================================================================================================================
		// To decide which incarnation of an ground-truth object is available
		//====================================================================================================================
    int gtType;

		//====================================================================================================================
		// This number tells the width (in samples) of the region in which an object might really be located when the 
		// groundtruth tells it is located at sample x; the uncertainty measured here is due to the groundtruthing- (i.e. 
		// hand-labeling) process; the total uncertainty is further increased by the frame-steps/window-steps in the feature-
		// extraction- and algorithm-stage and the information-storing (i.e. frame-quantisation in the frameList)
		// the complete uncertainty-area is returned by 
		//====================================================================================================================
    unsigned long int uncertaintyRegion;
		
		//====================================================================================================================
		// Some constants to model the relationship between types
		//====================================================================================================================
		static const long int speechRelatedTypes = sclib::atSpeech|sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atMaleVoice|sclib::atFemaleVoice|sclib::atVoiced|sclib::atUnvoiced|sclib::atPause|sclib::atSpeakerBoundary|sclib::atSpeechSegmentStart|sclib::atSpeechSegmentEnd;
		static const long int noiseRelatedTypes = sclib::atNoise|sclib::atAction|sclib::atBackground|sclib::atBreath|sclib::atMusic|sclib::atUndefined|sclib::atNoiseBoundary;
		static const long int detailedAudioTypes = sclib::atMusic|sclib::atAction|sclib::atBackground|sclib::atBreath|sclib::atPureSpeech|sclib::atNoisySpeech|sclib::atUndefined;

    //====================================================================================================================
    //  Inserts a speaker-name into the speakerNames[]-array and returns it's speaker-id (index into the array). If the
    //  name already exists, the existing index is returned; if anything goes wrong (e.g. speakerName is NULL or '\0' or 
    //  MAX_SPEAKRS speakers are reached), -1 (sclib::noSpeaker) is returned
    //====================================================================================================================
    int insertSpeakerName(const char* speakerName, int size);

		//====================================================================================================================
		// Initialize the frameList
		//====================================================================================================================
		virtual void initFrameList(void);
    
    //====================================================================================================================
    // The maximum speech segment length in the specified scene will be set to maxLength audioFrames; longer segments
    // will be split, with speaker-boundarys before each newly created segment!
    //====================================================================================================================
		void setMaxSpeechSegLength(unsigned long int segmentStart, unsigned long int segmenteEnd, unsigned long int maxLength, int origin = sclib::modeHypothesized);
		
		//====================================================================================================================
		// The minimum speech segment length in the specified segment will be set to minLength samples; shorter segments
		// will be erased
		//====================================================================================================================
		void setMinSpeechSegLength(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int minLength);

    //====================================================================================================================
    // method to generate pseudo scene boundarys at suitable positions
    // this is helpful if no scene-file is available to read real detected scene boundarys, but if one wishes to part the
    // audio into little pieces anyway (to reduce computational load and memory requirements)
    // this function analyzes the framelist and generates scene-boundarys approximately all 'approximateSceneLength' 
    // samples, but cares for not cutting homogenious speech segments into pieces. care also that each scene has some 
    // noise (max. 250ms after an extension because of a connected speech segment) for background estimation at it's start 
    // and end
    // ATTENTION: The generated boundaries are set for both the groundtruth- and hypothesized data
    //====================================================================================================================
    unsigned long int generatePseudoSceneBoundarys(unsigned long int approximateSceneLength);

		//====================================================================================================================
		// Methods to convert samples and internal audioFrames
		//====================================================================================================================
    unsigned long int   sample2FLI(unsigned long int sample);
    unsigned long int   FLI2sample(unsigned long int frameListIndex, unsigned int alignment = sclib::alignmentStart);

		//====================================================================================================================
		// Returns true if the questionable frame (internal frame, this method can't be called from outside!) has labels 
		//  according to this rules:
		//  - types:
		//    - uniteTypes==true : every label in types (concatenated with OR: e.g. sclib::atSpeech|sclib::atPause) must apply
		//                         to the frame
		//    - uniteTypes==false: min. one of the labels in types must apply to the frame
		//  - typesNot:
		//    - uniteTypesNot==t.: all of the labels in typesNot together must not apply to the frame (a single one is ok, though)
		//    - uniteTypesNot==f.: no single one of the labels in typesNot must apply to the frame
		//====================================================================================================================
    bool testFrame(unsigned long int FLI, unsigned int types = sclib::noType, bool uniteTypes = false, unsigned int typesNot = sclib::noType, bool uniteTypesNot = false, int origin = sclib::modeHypothesized);
		
		//====================================================================================================================
		//	Method to label the single frame FLI with specific marker(s) (type(s))
		//   - types: signle audioType-constant (SCLIB_AT_*) or OR-concatenated list of them
		//   - uniteTypes: if true and the labelsshould be removed, they are removed only if they all apply to the frame, 
		//                 otherwise each one is removed separately
		//	 - speakerID: the frame is given this speaker-id if it contains SPEECH; if speakerID==sclib::noSpeaker, [...]
		//   - forceNoSpeaker: [...] has to be true to set (or better: remove) it anyhow
		//   - action: sclib::modeLabelRemove or sclib::modeLabelAdd
		//   - orighin: sclib::modeGroundtruth or sclib::modeHypothesized
		//   - playByTheRules: if true, a set of heuristics is applied to make the action not violate any rule/assumption 
		//                     associated with the organisation of labels in the framelist; e.g. for label-adding it is 
		//                     checked then if the given types are consistent with mutual exclusiveness rules one by one, and 
		//                     already attached labels get removed if problems exist.
		//====================================================================================================================
    void setFrame(unsigned long int FLI, long int types, bool uniteTypes = false, int speakerID = sclib::noSpeaker, int action = sclib::modeLabelAdd, int origin = sclib::modeHypothesized, bool forceNoSpeaker = false, bool playByTheRules = false);

		//====================================================================================================================
		//  returns the next occasion of the given types/no-types combination (regarding the unite* flags as in testFrame())
		//  in the given direction; helpful for e.g. sclib::atSpeechSegment*
		//====================================================================================================================
		long int getNextFrame(unsigned long int initFLI, long int types = sclib::noType, bool uniteTypes = false, unsigned int typesNot = sclib::noType, bool uniteTypesNot = false, unsigned int direction = sclib::searchForward, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// Return the speaker-id for the given internal frame
		//====================================================================================================================
		long int getSpeakerID(unsigned long int FLI, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// Read Groundtruth from file(s) and initialize internal data structures (frameList etc.)
		//====================================================================================================================
    virtual bool readGroundTruth(void) = 0;

		//====================================================================================================================
		// To convert between different time measurement units
		//====================================================================================================================
		SC_Conversion *pConverter;

	public :
		
		//====================================================================================================================
		// Constructor, Destructor
		//====================================================================================================================
    SC_GroundTruth(SC_TweakableParameters *pTweak, const char* audioFileName);
		virtual ~SC_GroundTruth();
	
		//====================================================================================================================
		// Returns a SC_Signal object that includes the same parameter-settings as the originally opened file, but not the
    // samples
		//====================================================================================================================
    SC_Signal* getSignalPrototype(void) {return this->pSignalPrototype;};

		//====================================================================================================================
		// Methods to convert between the different time measurement units
		//====================================================================================================================
		unsigned long int sample2scene(unsigned long int sample, unsigned long int lastSample = 0, unsigned long int lastScene = 0);
		unsigned long int sample2shot(unsigned long int sample, unsigned long int lastSample = 0, unsigned long int lastShot = 0);

		//====================================================================================================================
		// Returns the number of columns of the probabilityList: Number of bitflags per FLI + 1 additional col for the 
		// probability of the speaker-id decision
		//====================================================================================================================
		unsigned long int getProbabilityListDim(void) {return (this->frameList != NULL) ? (sizeof(this->frameList[0][0])*8 + 1) : 0;}

		//====================================================================================================================
		//	Method to get the start and end of a segment (speech, noise, silence, ...) in terms of samples out of the 
		//  framelist. The behaviour of this method is dependant on the flag 'direction':
		//    - sclib::searchForward:  starting from 'initPosition', the next start- and endpoint of a segment of the 
		//                             desired type is returned (the standart case)
		//    - sclib::searchBackward: searching from 'initPosition' backward, the last segement completely before the 
		//                             'initPosition' is returned
		//    - sclib::searchMiddle:   it is searched for a startpoint before- and a endpoint after 'initPosition'; if a 
		//                             boundary lies on the initposition, the following segment (starting at 
		//                             'initPosition') is returned
		//    - sclib::searchWithin:   as sclib::searchForward, but if the startFrame already has the desired 
		//                             segmentType, it is returned as the segmentStart negelecting the fact that the
		//                             segment started earlier
		//  If the speaker-id changes within a segment and speakerBordersSegment==true, the changepoint is also considered as 
		//  the segment-end; the same is for other boundaries (speaker, scene or shot, may be OR-concatenated in the OR sense)
		//  when boundaryBordersSegment==true;
		//	ATTENTION:	segments are marked by continous labeling; e.g. if one frame is labeled as silence and the next is 
		//              not, this means the end of this segment!
		//====================================================================================================================
		void getNextSegment(unsigned long int initPosition, long int &segmentStart, long int &segmentEnd, long int type, unsigned int direction = sclib::searchForward, int origin = sclib::modeHypothesized, bool speakerBordersSegment = false, bool boundaryBordersSegment = false, unsigned long int boundaryType = sclib::atSceneBoundary);

		//====================================================================================================================
		//  If andTypes==true, this method just calls getNextSegment() and returns it's result; new functionality is added 
		//  for the other case: Then, the closest segment borders (minimum distance between appropriate segment-border and 
		//  initPosition) having one of the types in 'types' is returned together with it's type; it is assumed that the given
		//  types only appear mutually exclusive, a fact that can also be checked for if checkMutualTypeExclusiveness==true
		//====================================================================================================================
		unsigned long int getClosestSegment(unsigned long int initPosition, long int &segmentStart, long int &segmentEnd, long int types, bool uniteTypes = false, unsigned int direction = sclib::searchForward, int origin = sclib::modeHypothesized, bool speakerBordersSegment = false, bool boundaryBordersSegment = false, unsigned long int boundaryType = sclib::atSceneBoundary, bool checkMutualTypeExclusiveness = true);

		//====================================================================================================================
		//  returns the next occasion of the given type in the given direction; helpful for e.g. sclib::atSpeechSegment*
		//====================================================================================================================
		long int getNextOccasion(unsigned long int initPosition, long int type, unsigned int direction = sclib::searchForward, int origin = sclib::modeHypothesized, unsigned int alignment = sclib::alignmentStart);

		//====================================================================================================================
		//	Method to get the start and end of a boundary (scene, speaker, ...) in terms of samples out of the frameList
		//  the behaviour of this method is dependant on the flag 'direction':
		//    - sclib::searchForward:  starting from 'initPosition', the next start- and endpoint of a segment of the 
		//                             desired type is returned (the standart case)
		//    - sclib::searchBackward: searching from 'initPosition' backward, the last segement completely before the 
		//                             'initPosition' is returned
		//    - sclib::searchMiddle:   it is searched for a startpoint before- and a endpoint after 'initPosition'; if a 
		//                             boundary lies on the initposition, the following segment (starting at 
		//                             'initPosition') is returned
		//	ATTENTION:	boundarys in the frameList are only marked by their beginning, the start of the next one means the 
		//              end+1 of the last one!
		//====================================================================================================================
    void getNextBoundary(unsigned long int initPosition, long int &segmentStart, long int &segmentEnd, long int type, unsigned int direction = sclib::searchForward, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		//	method to test whether a given type exists at least once in the frameList between start and end (or, if end is bad 
		//  [<start, >frameCount], internalFrameCount) given in terms of samples.
		//====================================================================================================================
		bool existsSegmentType(unsigned long int segmentStart, unsigned long int segmentEnd, long int type, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		//	Method to label a segment (given by start- and end-point in samples) with specific marker(s) (types) by processing
		//  each frame in it with setFrame().
		//====================================================================================================================
		void setSegment(unsigned long int segmentStart, unsigned long int segmentEnd, long int types, bool uniteTypes = false, int speakerID = sclib::noSpeaker, int action = sclib::modeLabelAdd, int origin = sclib::modeHypothesized, bool forceNoSpeaker = false, bool playByTheRules = false);

		//====================================================================================================================
		//	Method to label a segment (given by start- and end-point in samples) with specific marker(s) (type(s));
		//  only those frames in the segment get processed by setFrame() which have the desired types; the number of actually 
		//  processed (by setFrame()) frames is returned
		//====================================================================================================================
		unsigned long int setSegmentIf(unsigned long int segmentStart, unsigned long int segmentEnd, long int ifTypes = sclib::noType, bool uniteIfTypes = false, long int ifTypesNot = sclib::noType, bool uniteIfTypesNot = false, long int types = sclib::noType, bool uniteTypes = false, int speakerID = sclib::noSpeaker, int action = sclib::modeLabelAdd, int origin = sclib::modeHypothesized, bool forceNoSpeaker = false, bool playByTheRules = false);

		//====================================================================================================================
		//	set the probability for the occurence of the given type(s) in the framelist from segmentStart to segmentEnd; 
		//  if type==sclib::noType, the probability for the speakerID-decision is set
		//====================================================================================================================
		void setProbability(unsigned long int segmentStart, unsigned long int segmentEnd, long int types, double probability);

		//====================================================================================================================
		// tests all frames in the given sample-borders and returns the number of samples meeting the requirements.
		// andFrames means: all frames have to meet the requirements, otherwise 0 is returned; for other parameters, see 
		// method testFrame() for their respective meaning.
		//====================================================================================================================
		unsigned long int testSegment(unsigned long int segmentStart, unsigned long int segmentEnd, bool andFrames = true, unsigned int types = sclib::noType, bool uniteTypes = false, unsigned int typesNot = sclib::noType, bool uniteTypesNot = false, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// can be used instead of testSegment() when - after each internal frame, which's size and concept is hidden outside 
		// this class - the result if the frame is positive or negative according to the given criteria is needed; this result
		// is returned in terms of nr. of positive samples in the positiveSamples parameter (>0 => FLI was positive), and 
		// positiveSamplePosition holds the startSample of this FLI; as long as the function-return-value is true, the tested 
		// FLI was between segmentStart and segmentEnd
		//====================================================================================================================
		bool testSegmentCallback(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned short int &positiveSamples, unsigned long int &positiveSamplePosition, unsigned int types = sclib::noType, bool uniteTypes = false, unsigned int typesNot = sclib::noType, bool uniteTypesNot = false, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// Test if the given segment (given by boundaries in samples) is speaker-homogenious (all frames are labeled with the 
		// same speaker-id)
		//====================================================================================================================
		bool isSpeakerHomogenious(long int segmentStart, long int segmentEnd, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// Test if the given segment (given by boundaries in samples) is speaker-homogenious (all frames are labeled with the 
		// same speaker-id) and if it's speaker is the one specified by speaker-name; does only make sense on ground-truth 
		// data, because only there speaker names may be available
		//====================================================================================================================
		bool isSpeakerHomogenious(long int segmentStart, long int segmentEnd, const char* speakerName);

		//====================================================================================================================
		// Return the speaker-id for the given sample-nr
		//====================================================================================================================
		long int getSamplesSpeakerID(unsigned long int sample, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		//  Returns the speaker-name for a given ground-truth speaker-id; returns "" if the given speaker-id is unknown
		//====================================================================================================================
		const char* getSpeakerName(int speakerGID);
	
    //====================================================================================================================
		// Returns the index (= ground-truth speaker-id) of this speaker-name into the speakerNames[]-array
    // returns -1 (sclib::noSpeaker) if the given speaker-name is not found
		//====================================================================================================================
    int getSpeakerIDfromName(const char* speakerName, int size);

		//====================================================================================================================
		// Returns the number of all samples having the sclib::atSpeech label within the speaker-boundary the given segment is 
		// belonging to (needed if the segment-length-threshold for clustering should be evaluated: A new cluster contains all 
		// speech belonging to one speaker-segment, so not the length of a single speech-segment has to be compared to the
		// segLengthThreshold, but all speech in it's speakers segment).
		// Here, origin==sclib::modeGroundtruth makes perfect sense as a default value 'cause the method is meant to be used
		// to forecast the number of speakers/ce's to be recognized by the algorithms depending on the groundtruth
		// ATTENTION: Relys on the fact that a new cluster always contains the speech of one speaker's segment, i.e. it is 
		//            created by SC_ModelHandler.buildSpeakerModels()
		//====================================================================================================================
		unsigned long int getSurroundingSpeakersSegmentSize(unsigned long int segmentStart, unsigned long int segmentEnd, int origin = sclib::modeGroundtruth);

		//====================================================================================================================
		//  To distinguish at runtime between different implementations
		//====================================================================================================================
		int getGTtype(void) {return this->gtType;}

		//====================================================================================================================
		// Returns the total width of the area in which an event may really have been occured of it is told in the 
		// groundtruth to have been occured at exact position X (the middle of that region);
    // The value of additionalLag (in samples) is added to the final result
		//====================================================================================================================
		unsigned long int getUncertaintyRegionWidth(bool includeGroundtruthUncertainty = true);

		//====================================================================================================================
		// The most common audio-type in the given segment is returned, as well as (in the percent parameter) the percentage 
		// with which it occurs; in the statistics-parameter, an array is returned in which each cell tells how many FLIs of 
		// the type of pow(2, index) are present in the segment (the last (32th) entry holds the sum of the previous ones); 
		// wantedTypes specifies which audio-types should be considered, default is all. 
		// This method exploits that the audio-types are just bitflags
		//====================================================================================================================
    unsigned long int getPrevailingAudioType(unsigned long int segmentStart, unsigned long int segmentEnd, double &percent, unsigned long int *&statistics, unsigned long int wantedTypes = sclib::atAllTypes, int origin = sclib::modeHypothesized);

    //====================================================================================================================
    //	Returns the ID of the speaker that appears most ofton in the given segment or sclib::noSpeaker if there is none
    //====================================================================================================================
		long int getMajorSpeakerID(unsigned long int segmentStart, unsigned long int segmentEnd, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// Go through the given segment and return:
		//  - overallGT: #samples having the type(s) according to groundTruth
		//  - overallHypo: #samples having the type(s) according to algorithmic results
		//  - overallContradiction: #samples that have the type(s) according to origin but not according to the other column
		//  - overallAgreement: #samples havin the type(s) in both GT and hypo columns
		//====================================================================================================================
		void getTypeStatistics(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int &overallGT, unsigned long int &overallHypo, unsigned long int &overallErrors, unsigned long int &overallAgreement, unsigned long int types, bool uniteTypes = false, int origin = sclib::modeGroundtruth);

    //====================================================================================================================
    //	Method that labels very short silence-segments (within speech) as 'pause' (instead of silence) and longer 
    //	sclib::atSilence segments as sclib::atSilence alone (without the tag sclib::atSpeech)
    //	thereby, 'very short' is defined by the parameter 'threshold', which is expressed in samples
    //====================================================================================================================
    void silence2pause(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned int threshold, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		//	This function marks speech segments (given by their boundarys in samples), which are shorter than 
		//  'segmentLengthThreshold' samples, with the sclib::atShort tag
		//====================================================================================================================
    void markShortSpeechSegments(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int segmentLengthThreshold);

    //====================================================================================================================
    // The speech segment length in this segment (given by borders in samples) will be set to the given parameters: 
    // Longer segments will be split, shorter segments will be erased. All values are sample-based 
    //====================================================================================================================
		void setSpeechSegLength(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int minLength, unsigned long int maxLength);

		//====================================================================================================================
		// get the count of audio-frames/the nr of the last frame fitting completely into a specified segment
		//====================================================================================================================
    unsigned long int getAudioFrameCountInSegment(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned int audioFrameSize, unsigned int audioFrameStep);
		unsigned long int getLastAudioFrameNrInSegment(unsigned long int segmentStart, unsigned long int segmentEnd, unsigned int audioFrameSize, unsigned int audioFrameStep);

		//====================================================================================================================
		// Speaker-specific segments in the frameList are modeled only by a flag at the beginning of a new speaker's segment,
		// so the proposed end of the last speech-segment is the beginning of the new one - 1. To get the real end of the last 
		// segment, this function counts backwards from the proposed end on till it finds the first sclib::atSpeech-frame. 
		// This is then the real segment-end; input. and output-values are sample-based
		//====================================================================================================================
    unsigned long int getSpeechEnd(unsigned long int proposedSpeechEnd, int origin = sclib::modeHypothesized);

 		//====================================================================================================================
		// Methods to do some debugging output
		//====================================================================================================================
		void frameListOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types = sclib::noType, unsigned long int typesNot = sclib::noType, sclib::OpenMode mode = ios_base::out|ios_base::app, int origin = sclib::modeHypothesized, bool suppressHeader = false);
    virtual void analyzedOut(const char* fileName, sclib::OpenMode mode = ios_base::out|ios_base::app);
		virtual void segmentStatisticsOut(const char* fileName, unsigned long int segmentStart, unsigned long int segmentEnd);
		virtual ostream& output(ostream& OutS, unsigned long int segmentStart, unsigned long int segmentEnd, unsigned long int types = sclib::noType, unsigned long int typesNot = sclib::noType, int origin = sclib::modeHypothesized, bool suppressHeader = false); 
		friend ostream&	operator<<(ostream& os, SC_GroundTruth& pGT);

    //====================================================================================================================
    // Copys all frames between start and end together in one new SV_Data-object, which has one of the frame-types
    // defined by 'types' (types may be multiple types added together with OR: e.g. sclib::atSpeech || sclib::atNoise || 
    // sclib::atPause) and definitely not one of the frame-types defined by 'typesNot'.
    // If needed, one can specifiy the parameters additionalCols and startCol: The new SV_Data-object will then have 
    // additionalCols more columns than it's parent, with the original data copied from startCol on.
		// If setStartSample is true and additionalCols is at least 1 and startCol is at least 1, the number of the first 
		// sample belonging to this frame (overall) is inserted in the zeroth col of each vector. If startCol or 
		// additionalCols is 0 but setting startSample is wished anyway, a second SV_Data-object is created with only 1 column
		// and rows corresponding to the rows of the returned object; this 2nd oject will have the sample-numbers in its 
		// single column and will be accessible via the ->Next pointer of the first SV_Data-object.
    // 'start' and 'end' must be real sample-numbers, not relative to 'offset'!!!
    // The offset means which sample is the first one in the first frame (first row) of the pComplete object
    //====================================================================================================================
    SV_Data* copyFramesTogether(SV_Data* pComplete, unsigned long int offset, unsigned long int start, unsigned long int end, unsigned long int types = sclib::noType, unsigned long int typesNot = sclib::noType, unsigned int additionalCols = 0, unsigned int startCol = 0, bool setFrameNr = false);

 		//====================================================================================================================
		// access to protected members
		//====================================================================================================================
    unsigned long int getAudioSampleRate(void) {return this->audioSampleRate;}
    unsigned long int getAudioSampleCount(void) {return this->audioSampleCount;}
    unsigned long int getInternalFrameSize(void) {return this->internalFrameSize;} //in samples
    char* getAudioFileName(void) {return this->audioFileName;}
		SC_Conversion* getConverter(void) {return this->pConverter;}

    //====================================================================================================================
    //  Read an ASCII-file containing a list when to load which explicit background model. It returns the nr of rows of 
    //  the table explicitModels, which has the form desired  by SC_ModelHandler::buildSpeakerModels():
    //    - explicitModels[x][0][0]: sceneNumber or 0 (then the modell applies to all scenes)
    //    - explicitModels[x][1][0]: validSegmentNumber or 0 (then the modell applies for all valid segments in the scene)
    //    - explicitModels[x][2]   : the full path/filename to the model-file to load
    //  Here's a constraint: the scene-/segment-numbers must be <=127 and the filename has to be <=255 characters in 
    //  length.
    //====================================================================================================================
    virtual unsigned short int readExplicitModelList(const char* modelListFileName, char*** &explicitModels);

		//====================================================================================================================
		// returns the numer of different names in the speakerNames[]-array; this corresponds to the nr. of different 
		// speakers if ground-truthing is carryed out the right way
		// counts only those speakers which appear in the scenes to evaluate: >= firstScene, <= lastScene, in sceneSelection
		// count only those whose segments where long enough to model/cluster
		// count only those segments between start & end (if end > 0, otherwise all)
		//====================================================================================================================
    unsigned int getRealSpeakerCount(bool considerSegmentLengthThresholds = true, unsigned long int start = 0, unsigned long int end = 0);

    //====================================================================================================================
    // counts the number of scenes to evaluate aacording to firstScene, lastScene and sceneSelection
    //====================================================================================================================
		unsigned long int getRealSceneCount(void);

		//====================================================================================================================
    // establishes (and after that returns) a mapping between speaker-id's from the groundtruth and from the clustering-
    // process; if the mapping already exists, just the 'correct'-flag is altered
    //====================================================================================================================
    void addSpeakerMapping(long int groundTruthID, long int hypothesizedID, bool correct = false);
    long int getSpeakerGIDFromHID(long int hypothesizedID, bool &correct);
    long int getSpeakerHIDFromGID(long int groundTruthID, bool &correct);
    void removeAllSpeakerMappings(void);

 		//====================================================================================================================
		// check the ground-truth (or hypothesized results) for consistency; this includes the check for mutually exclusive
		// tupels of audio-types (e.g. no internal frame can be speech and noise together); returns true if everything is ok
		// violations of rules are reported in a file if debug output SCLIB_BD_GT_INCOSISTENCIES is turned on
		//====================================================================================================================
		virtual bool checkConsistency(unsigned long int segmentStart, unsigned long int segmentEnd, int origin = sclib::modeHypothesized);

		//====================================================================================================================
		// Returns true if the first parameter represents a subtype of the second one (e.g. MUSIC is a subtype of NOISE), 
		// false if it isn't (e.g. because the first parameter is never a subtype)
		//====================================================================================================================
		static bool isSubType(long int questionableType, long int type);

    //====================================================================================================================
    // Returns true if the first parameter is semantically the opposite of the second one (e.g. SPEECH and NOISE, meaning 
    // the absence of one means the presence of the other), false if it isn't (e.g. BREATH and MUSIC, they don't have such 
    // a struct relationship although the< are mutually exclusive)
    //====================================================================================================================
    static bool isOppositeType(long int questionableType, long int type);

		//====================================================================================================================
		// Returns true if the first parameter and the second parameters (can be OR-concatenated) are mutually exclusive (e.g. 
		// VOICED and UNVOICED speech); in the parameter problem a OR-concatenated list of those types actally violating the
		// mutual exclusiveness is given
		//====================================================================================================================
		static bool areMutualExclusiveTypes(long int questionableType, long int types, long int &problem);

		//====================================================================================================================
		// Returns a string constant representing the Name (or abbreviation) of/for the given type
		//====================================================================================================================
		static const char* getAudioTypeName(long int audioType, bool shortName = false);

		//====================================================================================================================
		// Fills the buffer (and returns a pointer to it) with the concatenated (separated by the given separator) names of 
		// the types present in audioTypes (may be OR-concatenated)
		//====================================================================================================================
		static char* getAudioTypesNames(long int audioTypes, char *buffer, const char *separator = ", ", bool shortNames = false);

		//====================================================================================================================
		// Give access to semantically grouped audio-types
		//====================================================================================================================
		static const long int getSpeechRelatedTypes(void) {return SC_GroundTruth::speechRelatedTypes;};
		static const long int getNoiseRelatedTypes(void) {return SC_GroundTruth::noiseRelatedTypes;};
		static const long int getDetailedAudioTypes(void) {return SC_GroundTruth::detailedAudioTypes;};

		//====================================================================================================================
		// File-I/O for this class, so that the sate of a groundtruth-object can be saved and loaded from/to a file
		//====================================================================================================================
		virtual bool save(const char *fileName);
		virtual SC_GroundTruth* load(const char *fileName);

		//====================================================================================================================
		// Two methods that provide a way to load the complete signal of one speaker, for all speakes: The returns a (gt) 
		// speaker-id each time it is called until it returns false when all speaker-ids are returned; the second gives, for a
		// given (gt) speaker-id, all segments belonging to this speaker; with a nested loop, over all speakers, therein over 
		// all segments, therein doing loading, feature extraction and modeling, speaker models for all speakers in the corpus
		// can be built.
		//====================================================================================================================
		bool getSpeakersCallback(long int &speakerId);
		bool getSpeakersSegmentsCallback(long int speakerId, long int &segmentStart, long int &segmentEnd);

		//====================================================================================================================
		//  Calculate the Diarization error Rate (DER) according to Huang, Marchertet, Visweswariah, Potamianos, "The IBM RT07 
		//  Evaluation System for Speaker Diarization on Lecutre Meetings", 2007:
		//    "In accordance to NIST scoring, results are reported in terms of diarization error rate (DER). DER is calculated 
		//     by first finding the optimal one-to-one mapping between reference speakers and the hypothesized ones, and then 
		//     computing the percentage of time that is wrongly assigned according to the optimal mapping. DER includes 
		//     speaker error time, missed speaker time, and false alarm speaker time, thus also taking SAD errors into 
		//     account" (SAD=speech activity detection)
		//  It is placed here to save lots of computational time when the getNext* methods are used as compared to just
		//  traversing the internal data strucuture once herein
		//====================================================================================================================
		double calcDER(unsigned long int softBoundaryDiameter);
};

#endif
