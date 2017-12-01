/**************************************************************************/
/*    Responsibility:																											*/
/*      - container to do speech enhancement	  													*/
/*                                                                        */
/*    Does speech enhacement: It operates on a loaded buffer of time-     */
/*    discrete audio-samples and with the help of an externally created   */
/*    universal clean speech model (UCSM). By calling enhance(), the user */
/*    initiates the division of the samples into frames, followed by      */
/*    an enhancement of each single frame. After calling enhance(), the   */
/*		reconstructed clean speech resides in the buffer 'pEnhancedSignal'  */
/*		and can be copied or linked to a SV_Signal class to write it back  	*/
/*		into an audio-file or just base feature-extraction upon it				  */
/*																																				*/
/*		This class needs the frames defined in pTweak to be of length of a 	*/
/*		power of 2 with 50% overlap!  																			*/
/*																																				*/
/*    Author  : Thilo Stadelmann															            */
/*    Date    : 07.03.2005																								*/
/**************************************************************************/

#ifndef __SC_Enhancement_H__
#define __SC_Enhancement_H__

#include <stdlib.h>
#include "SC_Aux.h"
#include "SC_Api.h"
#include "SC_GroundTruth.h"
#include "SC_Transform.h"
#include "SC_MixtureModel.h"
#include "SC_ModelHandler.h"
#include "SC_Feature_Spectrum.h"
#include "SC_TweakableParameters.h"
#include <SV_Data.h>

class SCLIB_API SC_Enhancement {
  
  private :

  protected :

    short* pOriginalSignal; //just a pointer to the old signal (will not be modified, must be managed externally on the caller side)
    short* pEnhancedSignal; //the new enhanced signal
    unsigned long int signalLength; //length of - and pEnhancedSignal-buffer
    unsigned short int frameLength; //in samples, must be a power of 2; the frames will overlap 50%
    unsigned short int fftLength; //the next power of two starting from frameLength
    unsigned long int segmentStart; //first sample of pOriginalSignal's frames
    unsigned short int bewareLastFrame; //set by signal2frames() to indicate that the last frame was incomplete (the amount tells how many zeros where padded); used by frames2samples to reconstruct the signal
    short noiseModelUpdateMethod; //could be ENVIRONMENT or BOUNDARY, dependant on whethter the pGT-class holds information about noise-changes (boundaries) -> see aove #defines
    unsigned long int noiseModelLastFrame; //stores the right boundary of the frames specifying the current noise model (only relevant for the environment-update-rule in updateNoiseModel()); it is stored as a frame-index into the feature-set, not in absolute samples!
    unsigned long int noiseModelFirstFrame; //stores the left boundary...
    unsigned long int noiseModelLastUpdate; //the frame-index of evalutaed frame (currentFrameIdx) when the noise model was last updated
    bool verbose; //to control some console output...

    SC_ModelHandler* pModelHandler;
    SC_MixtureModel* pCleanSpeechModel; //the universal clean speech model (externally created, just loaded from a file)
    SC_MixtureModel* pNoiseModel; //built & updated in this class
    SC_GroundTruth* pGT; //for speech<->nonspeech-information
    SC_TweakableParameters* pTweak;
    SC_Transform *pTrans; //for the spectral transformation stuff
    SC_Feature_Spectrum *pExtractor; //used during noise-model-update

		//====================================================================================================================
		// converts the discrete audio samples into a matrix of successive, 50% overlapping frames.
		//====================================================================================================================
    SV_Data* signal2frames(void);

		//====================================================================================================================
		// converts the the matrix of succesive, 50% overlapping frames back into a discrete audio signal with exactly the 
		// same length as the original one.
		//====================================================================================================================
		bool frames2signal(SV_Data* pFrames);

		//====================================================================================================================
		// build the optimal noise-model for the background surrounding the frame with given index (relative to 
		// this->segmentStart)
		//====================================================================================================================
		bool updateNoiseModel(SV_Data* pFrames, unsigned long int currentFrameIdx);

		//====================================================================================================================
		// this is the real enhancer: gets a single frame (in samples) and returns an enhanced one.
		// it is assumed that the frameLength is a power of 2 as is controlled in the constructor
		// information about speech-nonspeech of frames needs to be present in the pGT-class
		// the UCSM and the noiseModel must be loaded.
		// frame and resFrame (and theire power-spectra) differ because of the influence of windowing prior to fft; this 
		// effect is compensated in the enhancement-framework by only taking the middle samples of each frame to construct the 
		// new signal.
		//====================================================================================================================
		double* mixmax(double* frame, bool usePostprocessing = true, bool useEnergyAdjustment = true);

  public :

		//====================================================================================================================
		// the constructor: the frameLength as specified in the pTweak class must be a power of 2 with 50% overlap
		//====================================================================================================================
    SC_Enhancement(SC_GroundTruth* pGT, SC_TweakableParameters* pTweak, const char* speechModelFile = "", const char* noiseModelFile = "", unsigned long int segmentStart = 0, short* pSignal = NULL, unsigned long int signalLength = 0, bool verbose = true);
    virtual ~SC_Enhancement();

		//====================================================================================================================
		// loads the universal clean speech model (UCSM) out of the specified file; returns true if successful
		// the old UCSM (if any) is freed anyway
		//====================================================================================================================
    bool setCleanSpeechModel(const char* fileName);

		//====================================================================================================================
		// loads an explicigt noise model (UCSM) out of the specified file; returns true if successful
		// the old noise model (if any) is freed anyway
		//====================================================================================================================
    bool setNoiseModel(const char* fileName);

		//====================================================================================================================
		// bend the internal pOriginalSignal-pointer to pSignal, set internal signalLength accordingly
		// the internal buffer is NOT freed before bending the pointer, because it is assumed that is stored & and managed
		// outside this class and this->pOriginalSignal is just a pointer to it
		//====================================================================================================================
		void setOriginalSignalPointer(short* pSignal, unsigned long int length);

		//====================================================================================================================
		// returns a pointer to the (internally stored, not just linked) enhanced signal and it's length (which is the same as
		// the length of the original signal)
		//====================================================================================================================
		void getEnhancedSignalPointer(short* &pSignal, unsigned long int &length);
		
		short* getEnhancedSignalPointer(void) {return this->pEnhancedSignal;}
		unsigned long int getEnhancedSignalLength(void) {return this->signalLength;}
		void setSegmentStart(unsigned long int newStart) {this->segmentStart = newStart; return;}
		
		//====================================================================================================================
		// the pointer to the enhaced signal is set to NULL without freeing the memory, so the enhanced signal remains after 
		// destruction of this class; this is usefull if it is still used by another class (SC_Signal*) without making a copy
		// before destructing the enhancer
		//====================================================================================================================
		void forgetEnhancedSignalPointer(void);

		//====================================================================================================================
		// this is the framework for speech-enhancement: it divides the signal into frames, calls enhancement for each frame,
		// and reassembles a new signal from the enhanced frames, which is stored in the class' pEnhancedSignal-buffer.
		// false is returned if a single frame fails to be enhanced or some preliminarys are missing (UCSM, signal)
		// frames need to overlap each other by 50% as controlled by the constructor
		//====================================================================================================================
    bool enhance(bool usePostprocessing = true, bool useEnergyAdjustment = true);
};

#endif
