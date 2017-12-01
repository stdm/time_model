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

#include <math.h>
#include <assert.h>
#include <limits.h>
#include "SC_Enhancement.h"
#include "SC_Aux.h"
#include <SV_Error.h>

//====================================================================================================================
// the constructor: the frameLength as specified in the pTweak class must be a power of 2 with 50% overlap
//====================================================================================================================
SC_Enhancement::SC_Enhancement(SC_GroundTruth* pGT, SC_TweakableParameters* pTweak, const char* speechModelFile, const char* noiseModelFile, unsigned long int segmentStart, short* pSignal, unsigned long int signalLength, bool verbose) {  
  this->verbose = verbose;
  this->pGT = pGT;
  if (this->pGT == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "SC_Enhancement needs the pGT-Argument. It mustn't be NULL!");
	}
  this->pTweak = pTweak;
  if (this->pTweak == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "SC_Enhancement needs the ptweak-Argument. It mustn't be NULL!");
	}
  this->pModelHandler = new SC_ModelHandler(this->pTweak, false);
	if (this->pModelHandler->isMixtureModel(this->pTweak->modelHandler.foregroundModelType) != true || this->pModelHandler->isMixtureModel(this->pTweak->modelHandler.backgroundModelType) != true) {
		REPORT_ERROR(SVLIB_BadArg, "SC_Enhancement needs simple mixture models for speech and noise");
  }

  //set and test frame-length and -overlap parameters
  this->segmentStart = segmentStart;
  this->frameLength = (unsigned short int)(this->pGT->getConverter()->ms2sample(this->pTweak->featureSpectrum.frameSize));
  if (sclib::isPowerOfTwo(this->frameLength) == false) {
		REPORT_ERROR(SVLIB_BadArg, "SC_Enhancement needs a frameLength of 2^n samples!");
	}
  if (this->pGT->getConverter()->ms2sample(this->pTweak->featureSpectrum.frameStep) != this->frameLength/2) {REPORT_ERROR(SVLIB_BadArg, "SC_Enhancement needs a frame-overlap of 50%!");}
  
  //initialize the fft-related stuff
  this->fftLength = this->frameLength;
  if (this->pTweak->featureSpectrum.FFTsize != this->frameLength) {
		REPORT_ERROR(SVLIB_BadArg, "SC_Enhancement needs a fft-size equal to the frame-length!");
	}
  this->pTrans = new SC_Transform(this->fftLength, sclib::wndHanning, this->pTweak->transform.taperingLength);

  //initialize the feature-extraction stuff
  this->pExtractor = new SC_Feature_Spectrum(this->pGT->getAudioSampleRate(), this->frameLength, this->frameLength/2, this->pTweak->featureSpectrum.preEmphasizeFactor, this->fftLength, this->pTweak->featureSpectrum.window, true, false); //we don't need phase here, but the log-spectrum!

  //load universal clean speech model
  this->pCleanSpeechModel = NULL;
  if (speechModelFile != NULL && strcmp(speechModelFile, "") != 0) {
    if (!setCleanSpeechModel(speechModelFile)) {
      REPORT_ERROR(SVLIB_Fail, "Universal clean speech model for the MIXMAX speech enhancer couldn't be loaded!");
    }
  }

  //load explicit noise model (if available)
  this->pNoiseModel = NULL;
  this->noiseModelUpdateMethod = sclib::updateGuess;
  if (noiseModelFile != NULL && strcmp(noiseModelFile, "") != 0) {
    if (!setNoiseModel(noiseModelFile)) {
      REPORT_ERROR(SVLIB_Fail, "Explicit noise model for the MIXMAX speech enhancer couldn't be loaded!");
    } else {
      this->noiseModelUpdateMethod = sclib::updateNoUpdate;
    }
  }
	
  this->pOriginalSignal = pSignal;
  this->signalLength = signalLength;
  this->pEnhancedSignal = NULL;
  this->pNoiseModel = NULL;

  return;
}

//====================================================================================================================
// the destructor: the new encanced signal is freed; if you need it any further, make a copy prior to destructing!
//====================================================================================================================
SC_Enhancement::~SC_Enhancement(void) {
  MFree_0D(this->pModelHandler);
  MFree_0D(this->pTrans);
  MFree_0D(this->pCleanSpeechModel);
  MFree_0D(this->pNoiseModel);
  MFree_0D(this->pExtractor);
  MFree_1D(this->pEnhancedSignal);
  
  return;
}

//====================================================================================================================
// loads the universal clean speech model (UCSM) out of the specified file; returns true if successful
// the old UCSM (if any) is freed anyway
//====================================================================================================================
bool SC_Enhancement::setCleanSpeechModel(const char* fileName) {
  MFree_0D(this->pCleanSpeechModel);
  this->pCleanSpeechModel = (SC_MixtureModel*)this->pModelHandler->loadModel(fileName, this->pTweak->modelHandler.foregroundModelType);

  return (this->pCleanSpeechModel != NULL) ? true : false;
}

//====================================================================================================================
// loads an explicigt noise model (UCSM) out of the specified file; returns true if successful
// the old noise model (if any) is freed anyway
//====================================================================================================================
bool SC_Enhancement::setNoiseModel(const char* fileName) {
  MFree_0D(this->pNoiseModel);
  this->pNoiseModel = (SC_MixtureModel*)this->pModelHandler->loadModel(fileName, this->pTweak->modelHandler.backgroundModelType);

  return (this->pNoiseModel != NULL) ? true : false;
}

//====================================================================================================================
// this is the framework for speech-enhancement: it divides the signal into frames, calls enhancement for each frame,
// and reassembles a new signal from the enhanced frames, which is stored in the class' pEnhancedSignal-buffer.
// false is returned if a single frame fails to be enhanced or some preliminarys are missing (UCSM, signal)
// frames need to overlap each other by 50% as controlled by the constructor
//====================================================================================================================
bool SC_Enhancement::enhance(bool usePostprocessing, bool useEnergyAdjustment) {
  double *frame, *resFrame;
  unsigned long int frameCount, sampleCount;
  SV_Data* pFrames;
  bool res;
	double lastPercentage = 0.0;
 
  if (this->pOriginalSignal == NULL) {return false;}
  if (this->pCleanSpeechModel == NULL) {return false;}
  if ((pFrames = signal2frames()) == NULL) {return false;}
    
  //some init stuff
  MArray_1D(frame, this->frameLength, double, "SC_Enhancement.enhance: frame");
  this->noiseModelLastFrame = 0; //neeed (perhaps) in updateNoiseModel()
  this->noiseModelFirstFrame = 0;
  this->noiseModelLastUpdate = 0;
  
  //loop over all frames of the original old noisy signal
  for (frameCount = 0; frameCount < (unsigned)(pFrames->Row); frameCount++) {

    if (this->pGT->testSegment(this->segmentStart+pGT->getConverter()->audioFrame2sample(frameCount, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentStart), this->segmentStart+pGT->getConverter()->audioFrame2sample(frameCount, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentEnd), true, sclib::atSpeech) > 0) { //enhance only speech frames
      if (updateNoiseModel(pFrames, frameCount) == true) { //update the noise-model; if it fails, there is no valid ground to base enhancing on, so just leave this frame as it is
        //convert to double :-(
        for (sampleCount = 0; sampleCount < this->frameLength; sampleCount++) {
          frame[sampleCount] = (double)(pFrames->Mat[frameCount][sampleCount]);
        }

        //enhance!
        resFrame = mixmax(frame, usePostprocessing, useEnergyAdjustment);
        if (resFrame == NULL) {  //enhancing should work
          MFree_0D(pFrames);
          MFree_1D(frame);
          return false;
        }

        //convert back to float :-(
        for (sampleCount = 0; sampleCount < this->frameLength; sampleCount++) {
          pFrames->Mat[frameCount][sampleCount] = (float)(resFrame[sampleCount]);
        }

        if (this->verbose == true) {
					lastPercentage = sclib::printPercentage(pFrames->Row, frameCount, lastPercentage, 1.0, !frameCount);
        }

        MFree_1D(resFrame);
      } //if noise-model-update succeeded
    } //if frame conatins speech, not noise

  } //for frameCount

  //convert the frames back into a time-discrete audio-signal
  res = frames2signal(pFrames);
  
  //sclib::vectorOut("orig.txt", this->pOriginalSignal, this->signalLength, true, this->pTweak);
  //sclib::vectorOut("new.txt", this->pEnhancedSignal, this->signalLength, true, this->pTweak);

  if (this->verbose == true) {
		sclib::printPercentage(1, 1, 0.0, 0.0, false);
  }

  MFree_0D(pFrames);
  MFree_1D(frame);
  return res;
}

//====================================================================================================================
// converts the discrete audio samples into a matrix of successive, 50% overlapping frames.
//====================================================================================================================
SV_Data* SC_Enhancement::signal2frames(void) {
  unsigned long int stepSize, frameCount, framesOverall, sampleCount;
  SV_Data* pFrames;
  
  if (this->pOriginalSignal == NULL) {return NULL;}
  
  //some init stuff
  stepSize = this->frameLength / 2; //50% overlap hard coded
  framesOverall = this->signalLength / (this->frameLength/2);
  this->bewareLastFrame = (unsigned short)(((framesOverall-1)*stepSize + (this->frameLength-1) + 1) - this->signalLength); //if this is non-zero, the last frame doesn't fit entirely into the originalSignal; it overlaps about the amount of this variable
  pFrames = new SV_Data(framesOverall, this->frameLength);
  pFrames->Hdr.frameSize = this->frameLength;
  pFrames->Hdr.frameStep = stepSize;
  pFrames->Hdr.sampleRate = this->pGT->getAudioSampleRate();

  //loop over all frames of the original old noisy signal
  for (frameCount = 0; frameCount < framesOverall; frameCount++) {
    //fill the actual frame; if it is longer than the originalSignal, fill it with zeros.
    for (sampleCount = 0; sampleCount < this->frameLength; sampleCount++) {
      if ((frameCount * stepSize + sampleCount) < this->signalLength) {
        pFrames->Mat[frameCount][sampleCount] = (float)(this->pOriginalSignal[frameCount * stepSize + sampleCount]);
      } else {
        pFrames->Mat[frameCount][sampleCount] = 0.0;
      }
    } //for sampleCount
  } //for frameCount

  pFrames->Hdr.frameSize = this->frameLength;
  pFrames->Hdr.frameStep = stepSize;
  pFrames->Hdr.sampleRate = this->pGT->getAudioSampleRate();

  return pFrames;
}

//====================================================================================================================
// converts the the matrix of succesive, 50% overlapping frames back into a discrete audio signal with exactly the 
// same length as the original one.
//====================================================================================================================
bool SC_Enhancement::frames2signal(SV_Data* pFrames) {
  unsigned long int frameCount, sampleCount, startSample, endSample, offset;

  if (pFrames == NULL) {return false;}
  
  //allocate space for the new enhanced signal within this class
  MFree_1D(this->pEnhancedSignal);
  MArray_1D(this->pEnhancedSignal, this->signalLength, short, "SC_Enhancement.enhance: this->pEnhancedSignal");

  //loop over all frames of the original old noisy signal
  for (frameCount = 0; frameCount < (unsigned)(pFrames->Row); frameCount++) {
    //reassemble the complete new enhanced signal
    //because there is 50% overlap in the frames, we will only take the frameLength/2 samples in the middle of each enhanced frame to
    //build the new signal. exception: the first frame conducts 3/4 from it's beginning, the last frame conducts it's last 3/4 samples
    //because of the lack of an successor/predecessor
    if (frameCount == 0) { //the first frame is taken form the beginning till 3/4 of it's length
      offset = 0;
      startSample = 0;
      endSample = this->frameLength*3/4;
    } else if (frameCount == pFrames->Row-1) { //the last frame is taken from 1/4 of it's lenth to it's end, if it wasn't padded with zeros
      offset = frameCount * this->frameLength/2 + this->frameLength/4;
      startSample = this->frameLength/4;
      endSample = this->frameLength - this->bewareLastFrame;
    } else { //normal frames are taken from 1/4 to 3/4 of their length
      offset = frameCount * this->frameLength/2 + this->frameLength/4;
      startSample = this->frameLength/4;
      endSample = this->frameLength*3/4;
    }

    for (sampleCount = startSample; sampleCount < endSample; sampleCount++) {
      this->pEnhancedSignal[offset+sampleCount-startSample] = (short)(floor(pFrames->Mat[frameCount][sampleCount]));
    } //for sampleCount
  } //for frameCount

  return true;
}

//====================================================================================================================
// build the optimal noise-model for the background surrounding the frame with given index (relative to 
// this->segmentStart)
//====================================================================================================================
bool SC_Enhancement::updateNoiseModel(SV_Data* pFrames, unsigned long int currentFrameIdx) {
  unsigned long int framesNeeded = this->pGT->getConverter()->ms2audioFrame(this->pTweak->enhancement.minNoiseDuration, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep); //we need at least xxx ms audio to build a reliable noise model (xxx specified in the tweakable parameters)
  long int noiseStart, noiseEnd, step;
  unsigned long int sample, frameNr, frameCount, maxFrameNr, minFrameNr;
  bool sign;
  SV_Data *pNoiseFrames = NULL, *pNoiseSpectrum = NULL;
  
  //no update wanted because of an explicitely loaded background model?
  if (this->noiseModelUpdateMethod == sclib::updateNoUpdate && this->pNoiseModel != NULL) {
    return true;
  }

  //if this hasn't been done yet, determine the best method for noise-model updating
  if (this->noiseModelUpdateMethod == sclib::updateGuess) {
    this->noiseModelUpdateMethod = (this->pGT->existsSegmentType(0, this->pGT->getAudioSampleCount(), sclib::atNoiseBoundary)) ? sclib::updateBoundary : sclib::updateEnvironment;
  }
  
  if (this->noiseModelUpdateMethod == sclib::updateBoundary) {

    if ((this->pGT->testSegment(this->segmentStart+this->pGT->getConverter()->audioFrame2sample(currentFrameIdx, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentStart), this->segmentStart+this->pGT->getConverter()->audioFrame2sample(currentFrameIdx, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentEnd), true, sclib::atNoiseBoundary) > 0) || (this->pNoiseModel == NULL)) { //if the background-model needs an update
      //select the correct set of frames for the new model
      this->pGT->getNextBoundary(this->segmentStart+this->pGT->getConverter()->audioFrame2sample(currentFrameIdx, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentStart), noiseStart, noiseEnd, sclib::atNoiseBoundary, sclib::searchMiddle);
      if (noiseStart == sclib::noSegment || noiseEnd == sclib::noSegment) { //this should not happen!
        REPORT_ERROR(SVLIB_BadData, "Expexted noise-boundarys during speech enhancement/noise-model-update, but found none");
        return false;
      }
      noiseStart = sclib::max(this->segmentStart, noiseStart); //the start must lie in this segment
      noiseEnd = sclib::min(noiseEnd, this->segmentStart+this->pGT->getConverter()->audioFrame2sample(pFrames->Row, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep)); //the end also
      pNoiseFrames = this->pGT->copyFramesTogether(pFrames, this->segmentStart, noiseStart, noiseEnd, sclib::atNoise, 0); //sclib::atUnvoiced|sclib::atPause|sclib::atSilence
      if (pNoiseFrames == NULL || pNoiseFrames->Row < (long)framesNeeded) {
        MFree_0D(pNoiseFrames);
        return false;
      }
      //this->pGT->storeSignal(this->pTweak, this->segmentStart, this->segmentStart + pFrames->Row, "noise.wav", sclib::atNoise);
    } else { //no need to update the model because there is no change in the background noise
      return true;
    }

  }  else { // if (this->noiseModelUpdateMethod == sclib::updateEnvironment) {
    
    if (this->pNoiseModel == NULL || (currentFrameIdx - this->noiseModelLastUpdate) >= this->pGT->getConverter()->ms2audioFrame(this->pTweak->enhancement.noiseModelUpdateRate, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep)) {
		//if (currentFrameIdx >= noiseModelLastFrame) { //update the noise-model if the possibility of new relevant noise-samples is exists
      //this is the normal case when no ECD (environment change detector) is available:
      //noise from before and after the frame# 'currentFrame' is used to build the model; 
      //the total amount of frames should reach xxx ms of audio (xxx specified in the tweakable parameters).
      step = 0;
      sign = true;
      frameNr = currentFrameIdx;
      frameCount = 0;
      maxFrameNr = 0;
      minFrameNr = std::numeric_limits<unsigned long int>::max();

      //search with alternating direction and increasing distance around the 'currentFrame' till the search passes
      //beyond the borders of the current segment or the xxx ms are full (xxx as in pTweak->enhancement.minNoiseDuration)
      pNoiseFrames = new SV_Data(framesNeeded, this->frameLength);
      while (frameCount < framesNeeded-1 && frameNr >= 0 && frameNr <= (unsigned long)(pFrames->Row-1)) {
        //copy the current frame if it consists of noise
        if (this->pGT->testSegment(this->segmentStart+this->pGT->getConverter()->audioFrame2sample(frameNr, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentStart), this->segmentStart+this->pGT->getConverter()->audioFrame2sample(frameNr, pFrames->Hdr.frameSize, pFrames->Hdr.frameStep, sclib::alignmentEnd), true, sclib::atNoise) > 0) { //sclib::atUnvoiced|sclib::atPause|sclib::atSilence
          for (sample = 0; sample < this->frameLength; sample++) {
            pNoiseFrames->Mat[frameCount][sample] = pFrames->Mat[frameNr][sample];
          }
          maxFrameNr = sclib::max(maxFrameNr, frameNr);
          minFrameNr = sclib::min(minFrameNr, frameNr);
          frameCount++; //this one must reach framesNeeded-1
        }
        
        //compute the next frameNr to look at
        step = (sign == true) ? abs(step)+1 : step*(-1);
        sign = !sign;
        frameNr = currentFrameIdx + step;

        //if 'frameNr' is outside the current segment, look in the other direction
        if (frameNr < 0) { //before this segment
          frameNr += 2*abs(step) + 1; //the next frame in the forward direction
          step = abs(step) + 1; //so that the next round doesn't yield the same frame
        }
        if (frameNr > (unsigned long)(pFrames->Row-1)) { //after this segment
          frameNr -= 2*abs(step); //the next frame in the backward direction
          step++;
          sign = false; //so that the next round yields the next frame in backward direction
        }
      } //while...

      if (frameCount < framesNeeded-1) { //not enough frames to build a model //TODO: why (-1)?
        MFree_0D(pNoiseFrames);      
        return false;
      }

      this->noiseModelLastFrame = maxFrameNr; //store the right boundary of the features building the noise model (as a frame-idx into the feature-set, not in samples!!!)
      this->noiseModelLastFrame = minFrameNr; //store the left boundary...
      this->noiseModelLastUpdate = currentFrameIdx; //store the frame when it was last updated
    } else {//noiseModel != NULL
      return true;
    }
  } //ENVIRONMENT-method

  //convert the sample-frames to log-spectrum features (as in SC_Feature_Spectrum)
  pNoiseSpectrum = this->pExtractor->ExtractFeature(pNoiseFrames);
  //this->pGT->storeSignal(this->pTweak, this->segmentStart, this->segmentStart + pFrames->Row, "noise.wav", sclib::atNoise, 0);
  MFree_0D(pNoiseFrames);
  
  //everything was ok till here, so (re-)build the model
  MFree_0D(this->pNoiseModel);
  this->pNoiseModel = (SC_MixtureModel*)this->pModelHandler->buildModel(pNoiseSpectrum, NULL, sclib::atNoise, 1);
  //sclib::classOut("features.txt", pNoiseSpectrum, this->pTweak);
  MFree_0D(pNoiseSpectrum);

  //sclib::classOut("noisemodel.txt", this->pNoiseModel, this->pTweak);
  //sclib::classOut("cleanspeechmodel.txt", this->pCleanSpeechModel, this->pTweak);

  return (this->pNoiseModel != NULL) ? true : false;
}
 
//====================================================================================================================
// this is the real enhancer: gets a single frame (in samples) and returns an enhanced one.
// it is assumed that the frameLength is a power of 2 as is controlled in the constructor
// information about speech-nonspeech of frames needs to be present in the pGT-class
// the UCSM and the noiseModel must be loaded.
// frame and resFrame (and theire power-spectra) differ because of the influence of windowing prior to fft; this 
// effect is compensated in the enhancement-framework by only taking the middle samples of each frame to construct the 
// new signal.
//====================================================================================================================
double* SC_Enhancement::mixmax(double* frame, bool usePostprocessing, bool useEnergyAdjustment) {
  unsigned short int k, K, i, I, j, J;
  double *powerSpectrum = NULL, *phaseAngle = NULL, *resFrame = NULL;
  double **hij, h, **q, f, F, *g, *G, **Ri, **Rj, p, qComplete;
  double Ez, z, m, v, wi, deltaK, newZ, sd;
  double qMax = -1.0 * numeric_limits<double>::max();
  double exponent, log_eps = log(numeric_limits<double>::epsilon());
  double oldEnergy, newEnergy, deltaEnergy;

  if (this->pNoiseModel == NULL || this->pCleanSpeechModel == NULL) {return NULL;}

  //get the log-power-spectrum
  powerSpectrum = pTrans->powerSpectrum(frame, this->frameLength, phaseAngle, false, true);

  //init-stuff
  K = this->frameLength/2 +1;
  I = this->pCleanSpeechModel->getMixtureCount();
  J = this->pNoiseModel->getMixtureCount();
  MArray_2D(Ri, I, K, double, "Ri");
  MArray_2D(Rj, J, K, double, "Rj");
  MArray_2D(hij, I, J, double, "h");
  MArray_2D(q, I, J, double, "q");
  MArray_1D(g, J, double, "g");
  MArray_1D(G, J, double, "G");

  //calculate Ri,Rj (the ratios f/F,g/G) and hij (the pdf for mixtures i,j) (f,F,g,G and hij are in the log-domain)
  for (k = 0; k < K; k++) {
    z = powerSpectrum[k]; 
    for (i = 0; i < I; i++) {
      m = this->pCleanSpeechModel->getMean(i, k); 
      v = this->pCleanSpeechModel->getVariance(i, k);
      sd = this->pCleanSpeechModel->getSd(i, k);
      f = this->pCleanSpeechModel->gaussSolver.logGaussian(z, m, v, sd);
      F = this->pCleanSpeechModel->gaussSolver.logErf(z, m, sd);
      if (F < f) {
        F = sclib::sLog(this->pCleanSpeechModel->gaussSolver.approxErf(z, exp(f), m, v, v));
      }
      Ri[i][k] = exp(f-F); //f/F;
      assert(sclib::isFinite(Ri[i][k]));
      for (j = 0; j < J; j++) {
        if (i == 0) {
          m = this->pNoiseModel->getMean(j, k); 
          v = this->pNoiseModel->getVariance(j, k);
          sd = this->pNoiseModel->getSd(j, k);
          g[j] = this->pNoiseModel->gaussSolver.logGaussian(z, m, v, sd);
          //assert(g[j] > -600);
          G[j] = this->pNoiseModel->gaussSolver.logErf(z, m, sd);
          if (G[j] < g[j]) {
            G[j] = sclib::sLog(this->pNoiseModel->gaussSolver.approxErf(z, exp(g[j]), m, v, v));
          }
          Rj[j][k] = exp(g[j]-G[j]); //g[j]/G[j];
          assert(sclib::isFinite(Rj[j][k]));
        }
        if (k == 0) {
          hij[i][j] = 0.0;
        }
        
        exponent = f + G[j] - F - g[j];
        assert(sclib::isFinite(exponent));
        if (exponent > -log_eps) { //1 + x == x => log(1+exp(x)) == x
          hij[i][j] += f + G[j];
        } else if (exponent < log_eps) { //1 + x == 1 => log(1+exp(x)) == log(1) == 0
          hij[i][j] += F + g[j];
        } else { //the standard case
          hij[i][j] += log(1.0 + exp(exponent)) + F + g[j];
        }
        assert(sclib::isFinite(hij[i][j]));
      } //j
    } //i
  } //k

  //calculate h (the complete pdf h(z)) (also in the log-domain)
  for (i = 0; i < I; i++) {
    wi = log(this->pCleanSpeechModel->getWeight(i));
    for (j = 0; j < J; j++) {
      q[i][j] = wi + log(this->pNoiseModel->getWeight(j)) + hij[i][j]; //in the paper: c_i * h_i(z)
      assert(sclib::isFinite(q[i][j]));
      if (q[i][j] > qMax) {
        qMax = q[i][j];
      }
    } //j
  } //i
  h = qMax;
  if (I>1 || J>1) {
    for (i = 0; i < I; i++) {
      for (j = 0; j < J; j++) {
        h += exp(q[i][j] - qMax); //in the paper: SUM_j{c_j * h_j(z)}
      } //j
    } //i
  } //I>1, J>1
  assert(sclib::isFinite(h));

  //calculate q (the class conditioned probability, NOT logarithmic!!!)
  qComplete = 0;
  for (i = 0; i < I; i++) {
    for (j = 0; j < J; j++) {
      q[i][j] = exp(q[i][j] - h);
      qComplete += q[i][j];
      assert(sclib::isFinite(q[i][j]));
    } //j
  } //i

  //now really do the enhancement!!!
  oldEnergy = 0.0;
  newEnergy = 0.0;
  for (k = 0; k < K; k++) {
    deltaK = (k < 37) ? log(0.35) : log(0.18);
    newZ = 0.0;
    for (i = 0; i < I; i++) {
      Ez = this->pCleanSpeechModel->getMean(i, k) - this->pCleanSpeechModel->getVariance(i, k) * Ri[i][k];
      for (j = 0; j < J; j++) {
        if (Ri[i][k] <= 0.0) {
          p = 0.0; 
        } else {
          p = 1.0 / (1.0 + (Rj[j][k] / Ri[i][k]));
        }
        newZ += (q[i][j] * (1.0/qComplete)) * (p*powerSpectrum[k] + (1.0-p)*Ez);
      } //j
    } //i

    oldEnergy += powerSpectrum[k];

    //nonlinear postprocessing
    if (usePostprocessing == true) {
      powerSpectrum[k] = sclib::max(powerSpectrum[k] + deltaK, newZ);
    } else {
      powerSpectrum[k] = newZ;
    }

    newEnergy += powerSpectrum[k];
  } //k

  //adjust total frame energy so that the enhanced signal doesn't appear so much quieter than the original one
  if (useEnergyAdjustment == true) {
    deltaEnergy = fabs(newEnergy - oldEnergy) / (double)(K);
    for (k = 0; k < K; k++) {
      powerSpectrum[k] += deltaEnergy;
    }
  }

  MFree_2D(Ri);
  MFree_2D(Rj);
  MFree_2D(hij);
  MFree_2D(q);
  MFree_1D(g);
  MFree_1D(G);

  //transform back to time-based signal
  pTrans->powerPhase2ft(powerSpectrum, phaseAngle, false, this->frameLength, false, true);
  resFrame = pTrans->ifft(powerSpectrum, phaseAngle, this->frameLength);

  MFree_1D(phaseAngle);
  MFree_1D(powerSpectrum);

  return resFrame;
}

//====================================================================================================================
// bend the internal pOriginalSignal-pointer to pSignal, set internal signalLength accordingly
// the internal buffer is NOT freed before bending the pointer, because it is assumed that is stored & and managed
// outside this class and this->pOriginalSignal is just a pointer to it
//====================================================================================================================
void SC_Enhancement::setOriginalSignalPointer(short* pSignal, unsigned long int length) {
  this->pOriginalSignal = pSignal; 
  this->signalLength = length;
  return;
}

//====================================================================================================================
// returns a pointer to the (internally stored, not just linked) enhanced signal and it's length (which is the same as
// the length of the original signal)
//====================================================================================================================
void SC_Enhancement::getEnhancedSignalPointer(short* &pSignal, unsigned long int &length) {
  pSignal = this->pEnhancedSignal; 
  length = this->signalLength;
  return;
}

//====================================================================================================================
// the pointer to the enhaced signal is set to NULL without freeing the memory, so the enhanced signal remains after 
// destruction of this class; this is usefull if it is still used by another class (SC_Signal*) without making a copy
// before destructing the enhancer
//====================================================================================================================
void SC_Enhancement::forgetEnhancedSignalPointer(void) {
  this->pEnhancedSignal = NULL; 
  this->signalLength = 0;
  return;
}
