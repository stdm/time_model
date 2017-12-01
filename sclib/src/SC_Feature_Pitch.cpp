/**************************************************************************/
/*    Wrapper around different pitch extraction (in [Hz]) algorithms:     */
/*     - the one from the SV_Lib (SV_Feature_Pitch)                       */
/*     - ESPS/Talkin pitch detector, Copyright (c) 1990-1996 Entropic     */
/*       Research Laboratory, Inc., which is described as follows:        */
/*                                                                        */
/*    Brief description: The pitch estimation algorithm is the ESPS       */
/*    get_f0 tool. The ESPS/waves+ get_f0 algorithm is the same as RAPT   */
/*    (Robust Algorithm for Pitch Tracking). RAPT is described by David   */
/*    Talkin in chapter 14 of the book "Speech Coding and Synthesis",     */
/*    edited by W.B. Kleijn and K.K. Paliwal, ISBN: 0-444-82169-4.        */
/*    Basically, RAPT uses normalised cross-correlation of the waveform   */
/*    to obtain f0 candidates. Dynamic programming is then used to find   */
/*    the optimal f0 track.                                               */
/*                                                                        */
/*    Malcolm Slaney (Author of the Auditory Toolbox for Matlab) writes   */
/*    on http://www.auditory.org/mhonarc/2003/msg00407.html about it in   */
/*    July 2003:                                                          */
/*                                                                        */
/*    I think the ESPS pitch tracker has the best combination of wide     */
/*    distribution and very careful testing and benchmarking.  It is the  */
/*    benchmark algorithm I use to test the performance of all new pitch  */
/*    trackers. I'm really happy that the code is available. (I expect    */
/*    better pitch trackers to be developed now.  I for one, will not     */
/*    accept a pitch paper for publication unless they demonstrate better */
/*    performance than Talkin's approach--or a similarly good tracker. No */
/*    excuses now!! :-).                                                  */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 18.12.2007																								*/
/**************************************************************************/

#include <float.h>
#include <limits.h>
#include "SC_Feature_Pitch.h"
#include "SC_Aux.h"
#include <SV_Feature_Pitch.h>
#include "SC_FeatureHandler.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_Pitch::SC_Feature_Pitch(int sampleRate, int frameLength, int frameStep, int method, bool takeSqrt, float esps_cand_thresh, float esps_lag_weight, float esps_freq_weight, float esps_trans_cost, float esps_trans_amp, float esps_trans_spec, float esps_voice_bias, float esps_double_cost, float esps_min_f0, float esps_max_f0, int esps_n_cands, float esps_wind_dur, bool verbose) : SV_Feature() {
  //general initializations
	this->Para.WinSz = frameLength;
	this->Para.StpSz = frameStep;
	this->Para.SRate = sampleRate; 
	this->verbose = verbose;
	this->method = method;
	this->takeSqrt = takeSqrt;

	//ESPS parameters
#ifdef SC_USE_ESPSPITCH
	this->get_f0_par.cand_thresh = esps_cand_thresh;
	this->get_f0_par.lag_weight = esps_lag_weight;
	this->get_f0_par.freq_weight = esps_freq_weight;
	this->get_f0_par.trans_cost = esps_trans_cost;
	this->get_f0_par.trans_amp = esps_trans_amp;
	this->get_f0_par.trans_spec = esps_trans_spec;
	this->get_f0_par.voice_bias = esps_voice_bias;
	this->get_f0_par.double_cost = esps_double_cost;
	this->get_f0_par.min_f0 = esps_min_f0;
	this->get_f0_par.max_f0 = esps_max_f0;
	this->get_f0_par.frame_step = this->Para.StpSz / (float)(this->Para.SRate); //in [s]
	this->get_f0_par.wind_dur = esps_wind_dur;
	this->get_f0_par.n_cands = esps_n_cands;
	this->get_f0_par.mean_f0 = 200; //Does nothing yet.  Will allow biasing toward a speaker's characteristic F0.
	this->get_f0_par.mean_f0_weight = (float)(0.0); //Does nothing yet.  Used with "mean_F0".
	this->get_f0_par.conditioning = 0; //Does nothing yet.

	//ESPS-specific initializations
	this->headF = NULL;
	this->tailF = NULL;
	this->cmpthF = NULL;
	this->pcands = NULL;
	this->cir_buff_growth_count = 0;
	this->locs = NULL;
	this->wReuse = 0;
	this->f0p = NULL;
	this->vuvp = NULL;
	this->rms_speech = NULL;
	this->acpkp = NULL;
	this->peaks = NULL;
	this->windstat = NULL;
	this->first_time = 1;
	this->debug_level = (this->verbose == true) ? 1 : 0;
	this->stat = NULL;
	this->mem = NULL;
	MArray_1D(this->ProgName, sclib::bufferSize, char, "SC_Feature_Pitch.SC_Feature_Pitch: ProgName");
	sprintf(this->ProgName, "ESPS Pitch extractor in sclib");

	//former static variables in methods
	this->static_memsize = 0;
	this->static_nframes_old = 0;
	this->static_foutput = NULL;
	this->static_ncoeff = 127;
	this->static_ncoefft = 0;
  this->static_co = NULL; 
	this->static_mem = NULL;
  this->static_fsize = 0;
	this->static_resid = 0;
	this->static_din = NULL;
	this->static_n0 = 0;
  this->static_wsize = 0;
	this->static_wind = NULL;
  this->static_wsize_hwindow = 0;
	this->static_wind_hwindow = NULL;
  this->static_wsize_hnwindow = 0;
	this->static_wind_hnwindow = NULL;
  this->static_nwind = 0;
	this->static_dwind = NULL;
  this->static_nwind_lpc = 0;
	this->static_dwind_lpc = NULL;
  this->static_dbdata = NULL;
	this->static_dbsize = 0;
  this->static_dbdata_crossfi = NULL;
	this->static_dbsize_crossfi = 0;
#endif
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_Pitch::~SC_Feature_Pitch() {
#ifdef SC_USE_ESPSPITCH
	MFree_1D(this->ProgName);
	destructESPSmembers();
#endif
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data* SC_Feature_Pitch::ExtractFeature(void) {
	SV_Data *pTmp = NULL, *pPitch = NULL;
	int frameSize = 0, frameStep = 0;
	SC_FeatureHandler handler(NULL, false); //the method we call doesn't need tweakable parameters
	double ignoreValues[1] = {0.0};

	switch (this->method) {
		case sclib::modeESPS:
#ifdef SC_USE_ESPSPITCH
			frameSize = this->Para.WinSz; 
			frameStep = this->Para.StpSz;
			this->Para.WinSz = 10 * this->Para.SRate/1000; //the standard for the get_f0() method is 10 ms
			this->Para.StpSz = this->Para.WinSz; //the ESPS pitch tracker needs equal frame-size and -step, so we set frameSize=frameStep
			this->get_f0_par.frame_step = this->Para.StpSz / (float)(this->Para.SRate); //in [s]

			pTmp = get_f0();

			if (pTmp != NULL && (frameSize != this->Para.WinSz || frameStep != this->Para.StpSz)) {
				pPitch = handler.convertFrameRate(pTmp, this->Len, frameSize, frameStep, ignoreValues);
				MFree_0D(pTmp);
			} else {
				pPitch = pTmp;
			}

			this->Para.WinSz = frameSize;
			this->Para.StpSz = frameStep;
			this->get_f0_par.frame_step = this->Para.StpSz / (float)(this->Para.SRate); //in [s]
#else
			REPORT_ERROR(SVLIB_Fail, "ESPS code not available");
#endif
			break;
		case sclib::modeSVlib:
			pPitch = svlibPitch();
			break;
		default:
			REPORT_ERROR(SVLIB_BadArg, "Choosen pitch detector unknown!\n");
			break;
	}

	//take the sqrt of the pitch values as suggested by Rose [2002] (forensics book)
	if (this->takeSqrt == true) {
		for (int y = 0; y < pPitch->Row; y++) {
			pPitch->Mat[y][0] = sqrt(pPitch->Mat[y][0]);
		}
	}

	return pPitch;
}

//====================================================================================================================
// the function to call from ExtractFeatures() to get the result of the SV_Lib pitch tracker
//====================================================================================================================
SV_Data* SC_Feature_Pitch::svlibPitch(void) {
	SV_Data *pPitch = NULL;
  SV_Feature_Pitch *pPitchExtractor = new SV_Feature_Pitch();
	
	//set parameters
	pPitchExtractor->Para.SRate = this->Para.SRate;
	pPitchExtractor->Para.WinSz = this->Para.WinSz;
	pPitchExtractor->Para.StpSz = this->Para.StpSz;
	pPitchExtractor->Para.RmvSilence = 0;
	pPitchExtractor->Para.Smooth = 1; //smooth pitch contour with 3-point median filter

	//get Pitch
	pPitchExtractor->setSignal(this->Sig, this->Len, false);
	pPitch = pPitchExtractor->ExtractFeature();
	pPitchExtractor->setSignal(NULL, 0, false);
	MFree_0D(pPitchExtractor);

	//set header
	pPitch->Hdr.ID = sclib::featurePitch;
	pPitch->Hdr.frameSize = this->Para.WinSz;
	pPitch->Hdr.frameStep = this->Para.StpSz;
	pPitch->Hdr.sampleRate = this->Para.SRate;
	pPitch->Hdr.Signature[1] = sclib::modeSVlib; //encode which extractor was used

	//convert pitch period in [ms] to pitch frequency in [Hz]
	for (int y = 0; y < pPitch->Row; y++) {
		pPitch->Mat[y][0] = (float)(1000.0) / pPitch->Mat[y][0];
	}

  return pPitch;
}

#ifdef SC_USE_ESPSPITCH
//====================================================================================================================
// the function to call from ExtractFeatures() to get the result of the ESPS pitch tracker
//====================================================================================================================
SV_Data* SC_Feature_Pitch::get_f0(void) {
	SV_Data *pPitch = NULL;
	F0_params *par = &(this->get_f0_par);
	long int buff_size, actsize, total_samps, sdstep = 0, start;
	int i, done, vecsize, rowCount;
	float *fSamples = NULL;
	float *f0p, *vuvp, *rms_speech, *acpkp;

  if (this->Len >= ((this->Para.StpSz/(float)(this->Para.SRate) * 2.0) + this->get_f0_par.wind_dur) * this->Para.SRate) {
		//Initialize variables in get_f0.c; allocate data structures; determine length and overlap of input frames to read.
		destructESPSmembers(); //kill what has been potentially initialized by a previous call to init_dp_f0()
		if (init_dp_f0(this->Para.SRate, par, &buff_size, &sdstep)==0 && buff_size<=INT_MAX && sdstep<=INT_MAX) {
			if (this->verbose == true) {
				printf("  SC_Feature_Pitch.init_dp_f0 returned buff_size %ld, sdstep %ld.\n", buff_size, sdstep);
			}

			total_samps = this->Len;
			if (buff_size > total_samps) {
				buff_size = total_samps;
			}

			//actsize = get_feasd_recs(sd_rec, 0L, buff_size, ihd, ifile);
			//actsize = getBuffer(fSamples, 0, buff_size, this->Sig, this->Len);
			start = 0;
			actsize = sclib::min(buff_size, this->Len);
			MArray_1D(fSamples, sclib::max(buff_size, sdstep), float, "SC_Feature_Pitch.get_f0: fSamples");

			//create final data container
			rowCount = sclib::getRowCount(this->Len, this->Para.WinSz, this->Para.StpSz);
			pPitch = new SV_Data(rowCount, 1); //pitch [Hz]
			pPitch->Hdr.frameSize = this->Para.WinSz;
			pPitch->Hdr.frameStep = this->Para.StpSz;
			pPitch->Hdr.ID = sclib::featurePitch;
			pPitch->Hdr.sampleRate = this->Para.SRate;
			pPitch->Hdr.Signature[1] = sclib::modeESPS; //encode which extractor was used

			//this->debug_level = 1000;
			rowCount = 0;
			while (true) {
				done = (actsize < buff_size) || (total_samps == buff_size);
				getBuffer(fSamples, start, actsize, this->Sig, this->Len);

				if (dp_f0(fSamples, (int)(actsize), (int)(sdstep), this->Para.SRate, par, &f0p, &vuvp, &rms_speech, &acpkp, &vecsize, done)) {
					REPORT_ERROR(SVLIB_Fail, "problem in SC_Feature_PitcESPS.dp_f0()");
					MFree_0D(pPitch);
					break;
				}

				//store findings (pitch and probability of voiced speech) in final feature matrix
				for (i = vecsize - 1; i >= 0; i--) {
					pPitch->Mat[rowCount][0] = f0p[i];
					//pPitch->Mat[rowCount][1] = vuvp[i]; ignore voicing 'cause it is implicitly encoded in (voiced <=> pitch>0)
					//ignore root mean squared measurement in rms_speech[i]
					//ignore peak normalized cross-correlation value used to find the pitch in acpkp[i]
					rowCount++;
				}

				if (done) {
					break;
				}

				start += sdstep; 
				actsize = sclib::min(buff_size, this->Len-start);
				total_samps -= sdstep;
				if (actsize > total_samps) {
					actsize = total_samps;
				}
			}

			MFree_1D(fSamples);

			//somehow, get_f0() yields some frames less than expected, so initialize those entries with 0
			for (i = rowCount; i < pPitch->Row; i++) {
				pPitch->Mat[i][0] = (float)(0.0);
				//pPitch->Mat[i][1] = (float)(0.0);
			}
		} else {
			REPORT_ERROR(SVLIB_Fail, "problem in SC_FeaturePitchESPS.init_dp_f0()\n");
		}
	}// else {
	//	REPORT_ERROR(SVLIB_BadData, "input range too small for analysis by get_f0\n");
	//}

	destructESPSmembers(); //free memory of members
	
	return pPitch;
}

//====================================================================================================================
// mimic a macro from the ESPS package that defines an assertion with printed error message
//====================================================================================================================
void SC_Feature_Pitch::spsassert(void *pointer, const char *msg) {
	if (!pointer) {
		REPORT_ERROR(SVLIB_Fail, msg);
	}

	return;
}

//====================================================================================================================
// get the signal window-wise, return nr. of samples in buffer (without padded zeros if "samples" is over)
//====================================================================================================================
int SC_Feature_Pitch::getBuffer(float *buffer, int start, int len, float *samples, int overallLen) {
	int x;
	int actualLen = 0;
	int loopLen = (len+start > overallLen) ? overallLen-start : len;

	for (x = 0; x < loopLen; x++) {
		buffer[x] = (float)(samples[start+x]);
		actualLen++;
	}
	for (x = actualLen; x < len; x++) { //pad buffer with 0th if samples is shorter than wanted
		buffer[x] = (float)(0.0);
	}

	return actualLen;
}

//====================================================================================================================
// free all ESPS data members
//====================================================================================================================
void SC_Feature_Pitch::destructESPSmembers(void) {
	Frame *pHook, *pNext;

	free(this->pcands);
	this->pcands = NULL;
	free(this->f0p);
	this->f0p = NULL;
	free(this->vuvp);
	this->vuvp = NULL;
	free(this->acpkp);
	this->acpkp = NULL;
	free(this->rms_speech);
	this->rms_speech = NULL;
	free(this->peaks);
	this->peaks = NULL;
	free(this->locs);
	this->locs = NULL;
	free(this->windstat);
	this->windstat = NULL;

	pHook = this->headF;
	this->headF = NULL;
	this->tailF = NULL;
	if (pHook != NULL) {
		pHook->prev->next = NULL;
		while (pHook != NULL) {
			pNext = pHook->next;

			free(pHook->cp->correl);
			free(pHook->cp);
			free(pHook->dp->locs);
			free(pHook->dp->pvals);
			free(pHook->dp->mpvals);
			free(pHook->dp->prept);
			free(pHook->dp->dpvals);
			free(pHook->dp);
			free(pHook);
			pHook = NULL;

			pHook = pNext;
		}
	}

	if (this->stat != NULL) {
		free(this->stat->stat);
		free(this->stat->rms);
		free(this->stat->rms_ratio);
		free(this->stat);
		this->stat = NULL;
	}

  free(this->mem);
  this->mem = NULL;

	this->static_memsize = 0;
	this->static_nframes_old = 0;
	free(this->static_foutput);
	this->static_foutput = NULL;
	this->static_ncoeff = 127;
	this->static_ncoefft = 0;
  free(this->static_co);
	this->static_co= NULL; 
	free(this->static_mem);
	this->static_mem= NULL;
  this->static_fsize = 0;
	this->static_resid = 0;
	free(this->static_din);
	this->static_din= NULL;
	this->static_n0 = 0;
  this->static_wsize = 0;
	free(this->static_wind);
	this->static_wind= NULL;
  this->static_wsize_hwindow = 0;
	free(this->static_wind_hwindow);
	this->static_wind_hwindow= NULL;
  this->static_wsize_hnwindow = 0;
	free(this->static_wind_hnwindow);
	this->static_wind_hnwindow= NULL;
  this->static_nwind = 0;
	free(this->static_dwind);
	this->static_dwind= NULL;
  this->static_nwind_lpc = 0;
	free(this->static_dwind_lpc);
	this->static_dwind_lpc= NULL;
  free(this->static_dbdata);
	this->static_dbdata= NULL;
	this->static_dbsize = 0;
  free(this->static_dbdata_crossfi);
	this->static_dbdata_crossfi= NULL;
	this->static_dbsize_crossfi = 0;

	return;
}

//====================================================================================================================
// Below is the code taken and adapted from the ESPS project
//====================================================================================================================

/*
 * This material contains unpublished, proprietary software of 
 * Entropic Research Laboratory, Inc. Any reproduction, distribution, 
 * or publication of this work must be authorized in writing by Entropic 
 * Research Laboratory, Inc., and must bear the notice: 
 *
 *    "Copyright (c) 1990-1996 Entropic Research Laboratory, Inc. 
 *                   All rights reserved"
 *
 * The copyright notice above does not evidence any actual or intended 
 * publication of this source code.     
 *
 * Written by:  David Talkin
 * Checked by:
 * Revised by:  Derek Lin, David Talkin
 *
 * Brief description:  Estimate speech fundamental frequency.
 *
 */

//static char *sccs_id = "@(#)dp_f0.c	1.14	10/21/96	ERL"; //dp_f0.c

/* A fundamental frequency estimation algorithm using the normalized
   cross correlation function and dynamic programming.  The algorithm
   implemented here is similar to that presented by B. Secrest and
   G. Doddington, "An integrated pitch tracking algorithm for speech
   systems", Proc. ICASSP-83, pp.1352-1355.  It is fully described
   by D. Talkin, "A robust algorithm for ptich tracking (RAPT)", in
   W. B. Kleijn & K. K. Paliwal (eds.) Speech Coding and Synthesis,
   (New York: Elsevier, 1995). */

/* For each frame, up to par->n_cands cross correlation peaks are
   considered as F0 intervals.  Each is scored according to its within-
   frame properties (relative amplitude, relative location), and
   according to its connectivity with each of the candidates in the
   previous frame.  An unvoiced hypothesis is also generated at each
   frame and is considered in the light of voicing state change cost,
   the quality of the cross correlation peak, and frequency continuity. */

/* At each frame, each candidate has associated with it the following
   items:
	its peak value
	its peak value modified by its within-frame properties
	its location
	the candidate # in the previous frame yielding the min. err.
		(this is the optimum path pointer!)
	its cumulative cost: (local cost + connectivity cost +
		cumulative cost of its best-previous-frame-match). */

/* Dynamic programming is then used to pick the best F0 trajectory and voicing
   state given the local and transition costs for the entire utterance. */

/* To avoid the necessity of computing the full crosscorrelation at
   the input sample rate, the signal is downsampled; a full ccf is
   computed at the lower frequency; interpolation is used to estimate the
   location of the peaks at the higher sample rate; and the fine-grained
   ccf is computed only in the vicinity of these estimated peak
   locations. */

//extern int  debug_level;
//extern char *ProgName;
  
/*
 * READ_SIZE: length of input data frame in sec to read
 * DP_CIRCULAR: determines the initial size of DP circular buffer in sec
 * DP_HIST: stored frame history in second before checking for common path 
 *      DP_CIRCULAR > READ_SIZE, DP_CIRCULAR at least 2 times of DP_HIST 
 * DP_LIMIT: in case no convergence is found, DP frames of DP_LIMIT secs
 *      are kept before output is forced by simply picking the lowest cost
 *      path
 */

#define READ_SIZE 0.2
#define DP_CIRCULAR 1.5
#define DP_HIST 0.5
#define DP_LIMIT 1.0

/* 
 * stationarity parameters -
 * STAT_WSIZE: window size in sec used in measuring frame energy/stationarity
 * STAT_AINT: analysis interval in sec in measuring frame energy/stationarity
 */
#define STAT_WSIZE 0.030
#define STAT_AINT 0.020

/*--------------------------------------------------------------------*/
int SC_Feature_Pitch::get_Nframes(long buffsize, int pad, int step)
{
  if (buffsize < pad)
    return (0);
  else
    return ((buffsize - pad)/step);
}


/*--------------------------------------------------------------------*/
int SC_Feature_Pitch::init_dp_f0(double freq, F0_params *par, long *buffsize, long *sdstep)
{
  int nframes;
  int i;
  int stat_wsize, agap, ind, downpatch;

/*
 * reassigning some constants 
 */

  this->tcost = par->trans_cost;
  this->tfact_a = par->trans_amp;
  this->tfact_s = par->trans_spec;
  this->vbias = par->voice_bias;
  this->fdouble = par->double_cost;
  this->frame_int = par->frame_step;
  
  this->step = round(this->frame_int * freq);
  this->size = round(par->wind_dur * freq);
  this->frame_int = (float)(((float)step)/freq);
  this->wdur = (float)(((float)size)/freq);
  this->start = round(freq / par->max_f0);
  this->stop = round(freq / par->min_f0);
  this->nlags = this->stop - this->start + 1;
  this->ncomp = this->size + this->stop + 1; /* # of samples required by xcorr
			      comp. per fr. */
  this->maxpeaks = 2 + (this->nlags/2);	/* maximum number of "peaks" findable in ccf */
  this->ln2 = log((float)(2.0));
  this->size_frame_hist = (int) (DP_HIST / frame_int);
  this->size_frame_out = (int) (DP_LIMIT / frame_int);

/*
 * SET UP THE D.P. WEIGHTING FACTORS:
 *      The intent is to make the effectiveness of the various fudge factors
 *      independent of frame rate or sampling frequency.                
 */
  
  /* Lag-dependent weighting factor to emphasize early peaks (higher freqs)*/
  this->lagwt = par->lag_weight/this->stop;
  
  /* Penalty for a frequency skip in F0 per frame */
  this->freqwt = par->freq_weight/this->frame_int;
  
  i = (int) (READ_SIZE *freq);
  if(this->ncomp >= this->step) nframes = ((i-this->ncomp)/this->step ) + 1;
  else nframes = i / this->step;

  /* *buffsize is the number of samples needed to make F0 computation
     of nframes DP frames possible.  The last DP frame is patched with
     enough points so that F0 computation on it can be carried.  F0
     computaion on each frame needs enough points to do

     1) xcross or cross correlation measure:
           enough points to do xcross - ncomp

     2) stationarity measure:
           enough to make 30 msec windowing possible - ind

     3) downsampling:
           enough to make filtering possible -- downpatch
 
     So there are nframes whole DP frames, padded with pad points
     to make the last frame F0 computation ok.

  */

  /* last point in data frame needs points of 1/2 downsampler filter length 
     long, 0.005 is the filter length used in downsampler */
  downpatch = (((int) (freq * 0.005))+1) / 2;

  stat_wsize = (int) (STAT_WSIZE * freq);
  agap = (int) (STAT_AINT * freq);
  ind = ( agap - stat_wsize ) / 2;
  i = stat_wsize + ind;
  pad = downpatch + ((i>ncomp) ? i:ncomp);
  *buffsize = nframes * this->step + this->pad;
  *sdstep = nframes * this->step;
  
  /* Allocate space for the DP storage circularly linked data structure */

  size_cir_buffer = (int) (DP_CIRCULAR / frame_int);

  /* creating circularly linked data structures */
  this->tailF = alloc_frame(this->nlags, par->n_cands);
  this->headF = this->tailF;

  /* link them up */
  for(i=1; i<this->size_cir_buffer; i++){
    this->headF->next = alloc_frame(this->nlags, par->n_cands);
    this->headF->next->prev = this->headF;
    this->headF = this->headF->next;
  }
  this->headF->next = this->tailF;
  this->tailF->prev = this->headF;

  this->headF = this->tailF;

  /* Allocate sscratch array to use during backtrack convergence test. */
  if( ! this->pcands ) {
    this->pcands = (int *) malloc( par->n_cands * sizeof(int));
    spsassert(this->pcands,"can't allocate pathcands");
  }

  /* Allocate arrays to return F0 and related signals. */

  /* Note: remember to compare *vecsize with size_frame_out, because
     size_cir_buffer is not constant */
  this->output_buf_size = this->size_cir_buffer;
  this->rms_speech = (float*)malloc(sizeof(float) * this->output_buf_size);
  spsassert(this->rms_speech,"rms_speech malloc failed");
  this->f0p = (float*)malloc(sizeof(float) * this->output_buf_size);
  spsassert(this->f0p,"f0p malloc failed");
  this->vuvp = (float*)malloc(sizeof(float)* this->output_buf_size);
  spsassert(this->vuvp,"vuvp malloc failed");
  this->acpkp = (float*)malloc(sizeof(float) * this->output_buf_size);
  spsassert(this->acpkp,"acpkp malloc failed");

  /* Allocate space for peak location and amplitude scratch arrays. */
  this->peaks = (float*)malloc(sizeof(float) * this->maxpeaks);
  spsassert(this->peaks,"peaks malloc failed");
  this->locs = (int*)malloc(sizeof(int) * this->maxpeaks);
  spsassert(this->locs, "locs malloc failed");
  
  /* Initialise the retrieval/saving scheme of window statistic measures */
  this->wReuse = agap / step;
  if (this->wReuse){
      this->windstat = (Windstat *) malloc( this->wReuse * sizeof(Windstat));
      spsassert(this->windstat, "windstat malloc failed");
      for(i=0; i<this->wReuse; i++){
	  this->windstat[i].err = 0;
	  this->windstat[i].rms = 0;
      }
  }

  if(this->debug_level){
    Fprintf(stderr, "%s: done with initialization:\n", this->ProgName);
    Fprintf(stderr,
	    " size_cir_buffer:%d  xcorr frame size:%d start lag:%d nlags:%d\n",
	    this->size_cir_buffer, this->size, this->start, this->nlags);
  }

  this->num_active_frames = 0;
  this->first_time = 1;

  return(0);
}
  

/*--------------------------------------------------------------------*/
int SC_Feature_Pitch::dp_f0(float *fdata, int buff_size, int sdstep, double freq, F0_params *par, float **f0p_pt, float **vuvp_pt, float **rms_speech_pt, float **acpkp_pt, int *vecsize, int last_time)
{
  float  maxval, engref, *sta, *rms_ratio, *dsdata;
  register float ttemp, ftemp, ft1, ferr, err, errmin;
  register int  i, j, k, loc1, loc2;
  int   nframes, maxloc, ncand, ncandp, minloc,
        decimate, samsds;

	Stat *stat = NULL;

  nframes = get_Nframes((long) buff_size, pad, this->step); /* # of whole frames */

  if(debug_level)
    Fprintf(stderr,
	    "%s: ******* Computing %d dp frames ******** from %d points\n",
	    ProgName, nframes, buff_size);

  /* Now downsample the signal for coarse peak estimates. */

  decimate = (int)(freq/2000.0);	/* downsample to about 2kHz */
  if (decimate <= 1)
    dsdata = fdata;
  else {
    samsds = ((nframes-1) * step + ncomp) / decimate;
    dsdata = downsample(fdata, buff_size, sdstep, freq, &samsds, decimate, 
			this->first_time, last_time);
    if (!dsdata) {
      Fprintf(stderr, "%s: can't get downsampled data.\n", ProgName);
      return 1;
    }
  }

  /* Get a function of the "stationarity" of the speech signal. */

  stat = get_stationarity(fdata, freq, buff_size, nframes, step, first_time);
  if (!stat) { 
    Fprintf(stderr, "%s: can't get stationarity\n", ProgName);
    return(1);
  }
  sta = stat->stat;
  rms_ratio = stat->rms_ratio;

  /***********************************************************************/
  /* MAIN FUNDAMENTAL FREQUENCY ESTIMATION LOOP */
  /***********************************************************************/
  if(!first_time && nframes > 0) headF = headF->next;

  for(i = 0; i < nframes; i++) {
 
    /* NOTE: This buffer growth provision is probably not necessary.
       It was put in (with errors) by Derek Lin and apparently never
       tested.  My tests and analysis suggest it is completely
       superfluous. DT 9/5/96 */
    /* Dynamically allocating more space for the circular buffer */
    if(headF == tailF->prev){
			Frame *frm;

			if(cir_buff_growth_count > 5){
				Fprintf(stderr,
					"%s: too many requests (%d) for dynamically allocating space.\n   There may be a problem in finding converged path.\n",
					ProgName, cir_buff_growth_count);
				return(1);
      }
      if(debug_level) 
				Fprintf(stderr, "%s: allocating %d more frames for DP circ. buffer.\n",
					ProgName, size_cir_buffer);
      frm = alloc_frame(nlags, par->n_cands);
      headF->next = frm;
      frm->prev = headF;
      for(k=1; k<size_cir_buffer; k++){
				frm->next = alloc_frame(nlags, par->n_cands);
				frm->next->prev = frm;
				frm = frm->next;
      }
      frm->next = tailF;
      tailF->prev = frm;
      cir_buff_growth_count++;
		}

    headF->rms = stat->rms[i];
    get_fast_cands(fdata, dsdata, i, step, size, decimate, start,
		   nlags, &engref, &maxloc,
		   &maxval, headF->cp, peaks, locs, &ncand, par);
    

    /*    Move the peak value and location arrays into the dp structure */
		{
      register float *ftp1, *ftp2;
      register short *sp1;
      register int *sp2;
      
      for(ftp1 = headF->dp->pvals, ftp2 = peaks,
					sp1 = headF->dp->locs, sp2 = locs, j=ncand; j--; ) {
				*ftp1++ = *ftp2++;
				*sp1++ = *sp2++;
      }
      *sp1 = -1;		/* distinguish the UNVOICED candidate */
      *ftp1 = maxval;
      headF->dp->mpvals[ncand] = vbias+maxval; /* (high cost if cor. is high)*/
    }

    /* Apply a lag-dependent weight to the peaks to encourage the selection
       of the first major peak.  Translate the modified peak values into
       costs (high peak ==> low cost). */
    for(j=0; j < ncand; j++){
      ftemp = (float)(1.0 - ((float)locs[j] * lagwt));
      headF->dp->mpvals[j] = (float)(1.0 - (peaks[j] * ftemp));
    }
    ncand++;			/* include the unvoiced candidate */
    headF->dp->ncands = ncand;

    /*********************************************************************/
    /*    COMPUTE THE DISTANCE MEASURES AND ACCUMULATE THE COSTS.       */
    /*********************************************************************/

    ncandp = headF->prev->dp->ncands;
    for(k=0; k<ncand; k++){	/* for each of the current candidates... */
      minloc = 0;
      errmin = FLT_MAX;
			if((loc2 = headF->dp->locs[k]) > 0) { /* current cand. is voiced */
				for(j=0; j<ncandp; j++){ /* for each PREVIOUS candidate... */
					/*    Get cost due to inter-frame period change. */
					loc1 = headF->prev->dp->locs[j];
					if (loc1 > 0) { /* prev. was voiced */
						ftemp = log((float)(((double) loc2) / loc1));
						ttemp = fabs(ftemp);
						ft1 = fdouble + fabs(ftemp + ln2);
						if (ttemp > ft1)
							ttemp = ft1;
						ft1 = fdouble + fabs(ftemp - ln2);
						if (ttemp > ft1)
							ttemp = ft1;
						ferr = ttemp * freqwt;
					} else {		/* prev. was unvoiced */
						ferr = tcost + (tfact_s * sta[i]) + (tfact_a / rms_ratio[i]);
					}
					/*    Add in cumulative cost associated with previous peak. */
					err = ferr + headF->prev->dp->dpvals[j];
					if(err < errmin){	/* find min. cost */
						errmin = err;
						minloc = j;
					}
				}
      } else {			/* this is the unvoiced candidate */
				for(j=0; j<ncandp; j++){ /* for each PREVIOUS candidate... */
	  
					/*    Get voicing transition cost. */
					if (headF->prev->dp->locs[j] > 0) { /* previous was voiced */
						ferr = tcost + (tfact_s * sta[i]) + (tfact_a * rms_ratio[i]);
					}
					else
						ferr = 0.0;
					/*    Add in cumulative cost associated with previous peak. */
					err = ferr + headF->prev->dp->dpvals[j];
					if(err < errmin){	/* find min. cost */
						errmin = err;
						minloc = j;
					}
				}
			}
      /* Now have found the best path from this cand. to prev. frame */
      if (first_time && i==0) {		/* this is the first frame */
				headF->dp->dpvals[k] = headF->dp->mpvals[k];
				headF->dp->prept[k] = 0;
      } else {
				headF->dp->dpvals[k] = errmin + headF->dp->mpvals[k];
				headF->dp->prept[k] = minloc;
			}
		} /*    END OF THIS DP frame */

    if (i < nframes - 1)
      headF = headF->next;
    
    if (debug_level >= 2) {
      Fprintf(stderr,"%d engref:%10.0f max:%7.5f loc:%4d\n",
	      i,engref,maxval,maxloc);
    }
    
  } /* end for (i ...) */

  /***************************************************************/
  /* DONE WITH FILLING DP STRUCTURES FOR THE SET OF SAMPLED DATA */
  /*    NOW FIND A CONVERGED DP PATH                             */
  /***************************************************************/

  *vecsize = 0;			/* # of output frames returned */

  num_active_frames += nframes;

  if( num_active_frames >= size_frame_hist  || last_time ){
		Frame *frm;
    int  num_paths, best_cand, frmcnt, checkpath_done = 1;
    float patherrmin;
      
    if(debug_level)
      Fprintf(stderr, "%s: available frames for backtracking: %d\n", 
	      ProgName, num_active_frames);
      
    patherrmin = FLT_MAX;
    best_cand = 0;
    num_paths = headF->dp->ncands;

    /* Get the best candidate for the final frame and initialize the
       paths' backpointers. */
    frm = headF;
    for(k=0; k < num_paths; k++) {
      if (patherrmin > headF->dp->dpvals[k]){
				patherrmin = headF->dp->dpvals[k];
				best_cand = k;	/* index indicating the best candidate at a path */
      }
      pcands[k] = frm->dp->prept[k];
    }

    if(last_time){     /* Input data was exhausted. force final outputs. */
      this->cmpthF = this->headF;		/* Use the current frame as starting point. */
    } else {
      /* Starting from the most recent frame, trace back each candidate's
				 best path until reaching a common candidate at some past frame. */
      frmcnt = 0;
      while (1) {
				frm = frm->prev;
				frmcnt++;
				checkpath_done = 1;
				for(k=1; k < num_paths; k++){ /* Check for convergence. */
					if(pcands[0] != pcands[k])
						checkpath_done = 0;
				}
				if( ! checkpath_done) { /* Prepare for checking at prev. frame. */
					for(k=0; k < num_paths; k++){
						pcands[k] = frm->dp->prept[pcands[k]];
					}
				} else {	/* All paths have converged. */
					cmpthF = frm;
					best_cand = pcands[0];
					if(debug_level)
						Fprintf(stderr,
							"%s: paths went back %d frames before converging\n",
							ProgName, frmcnt);
					break;
				}
				if(frm == tailF){	/* Used all available data? */
					if( num_active_frames < size_frame_out) { /* Delay some more? */
						checkpath_done = 0; /* Yes, don't backtrack at this time. */
						cmpthF = NULL;
					} else {		/* No more delay! Force best-guess output. */
						checkpath_done = 1;
						cmpthF = headF;
						Fprintf(stderr,
							"%s: WARNING: no converging path found after going back %d frames, will use the lowest cost path\n",
							ProgName, num_active_frames);
					}
					break;
				} /* end if (frm ...) */
      }	/* end while (1) */
		} /* end if (last_time) ... else */

    /*************************************************************/
    /* BACKTRACKING FROM cmpthF (best_cand) ALL THE WAY TO tailF    */
    /*************************************************************/
    i = 0;
    frm = cmpthF;	/* Start where convergence was found (or faked). */
		while( frm != tailF->prev && checkpath_done){
      if( i == output_buf_size ){ /* Need more room for outputs? */
				output_buf_size *= 2;
				if(debug_level)
					Fprintf(stderr,
						"%s: reallocating space for output frames: %d\n",
						ProgName, output_buf_size);
				rms_speech = (float *) realloc((char *) rms_speech,
										sizeof(float) * output_buf_size);
				spsassert(rms_speech, "rms_speech realloc failed in dp_f0()");
				f0p = (float *) realloc((char *) f0p,
							sizeof(float) * output_buf_size);
				spsassert(f0p, "f0p realloc failed in dp_f0()");
				vuvp = (float *) realloc(vuvp, sizeof(float) * output_buf_size);
				spsassert(vuvp, "vuvp realloc failed in dp_f0()");
				acpkp = (float *) realloc(acpkp, sizeof(float) * output_buf_size);
				spsassert(acpkp, "acpkp realloc failed in dp_f0()");
			}
      rms_speech[i] = frm->rms;
      acpkp[i] =  frm->dp->pvals[best_cand];
      loc1 = frm->dp->locs[best_cand];
      vuvp[i] = 1.0;
      best_cand = frm->dp->prept[best_cand];
      ftemp = (float)(loc1);
      if(loc1 > 0) {		/* Was f0 actually estimated for this frame? */
				if (loc1 > start && loc1 < stop) { /* loc1 must be a local maximum. */
					float cormax, cprev, cnext, den;
					  
					j = loc1 - start;
					cormax = frm->cp->correl[j];
					cprev = frm->cp->correl[j+1];
					cnext = frm->cp->correl[j-1];
					den = (float)(2.0 * ( cprev + cnext - (2.0 * cormax) ));
					/*
					* Only parabolic interpolate if cormax is indeed a local 
					* turning point. Find peak of curve that goes though the 3 points
					*/
					  
					if (fabs(den) > 0.000001)
						ftemp += (float)(2.0 - ((((5.0*cprev)+(3.0*cnext)-(8.0*cormax))/den)));
				}
				f0p[i] = (float)(freq/ftemp);
			} else {		/* No valid estimate; just fake some arbitrary F0. */
				f0p[i] = 0;
				vuvp[i] = 0.0;
      }
      frm = frm->prev;
	  
      if (debug_level >= 2)
				Fprintf(stderr," i:%4d%8.1f%8.1f\n",i,f0p[i],vuvp[i]);
      /* f0p[i] starts from the most recent one */ 
      /* Need to reverse the order in the calling function */
      i++;
    } /* end while() */
    if (checkpath_done){
      *vecsize = i;
      tailF = cmpthF->next;
      num_active_frames -= *vecsize;
    }
  } /* end if() */

  if (debug_level)
    Fprintf(stderr, "%s: writing out %d frames.\n", ProgName, *vecsize);
  
  *f0p_pt = f0p;
  *vuvp_pt = vuvp;
  *acpkp_pt = acpkp;
  *rms_speech_pt = rms_speech;
  *acpkp_pt = acpkp;
  
  if(first_time) first_time = 0;
  return(0);
}


/*--------------------------------------------------------------------*/
SC_Feature_Pitch::Frame* SC_Feature_Pitch::alloc_frame(int nlags, int ncands)
{
	Frame *frm;
  int j;

  frm = (Frame*)malloc(sizeof(Frame));
  frm->dp = (Dprec *) malloc(sizeof(Dprec));
  spsassert(frm->dp,"frm->dp malloc failed in alloc_frame");
  frm->dp->ncands = 0;
  frm->cp = (Cross *) malloc(sizeof(Cross));
  spsassert(frm->cp,"frm->cp malloc failed in alloc_frame");
  frm->cp->correl = (float *) malloc(sizeof(float) * nlags);
  spsassert(frm->cp->correl, "frm->cp->correl malloc failed");
  /* Allocate space for candidates and working arrays. */
  frm->dp->locs = (short*)malloc(sizeof(short) * ncands);
  spsassert(frm->dp->locs,"frm->dp->locs malloc failed in alloc_frame()");
  frm->dp->pvals = (float*)malloc(sizeof(float) * ncands);
  spsassert(frm->dp->pvals,"frm->dp->pvals malloc failed in alloc_frame()");
  frm->dp->mpvals = (float*)malloc(sizeof(float) * ncands);
  spsassert(frm->dp->mpvals,"frm->dp->mpvals malloc failed in alloc_frame()");
  frm->dp->prept = (short*)malloc(sizeof(short) * ncands);
  spsassert(frm->dp->prept,"frm->dp->prept malloc failed in alloc_frame()");
  frm->dp->dpvals = (float*)malloc(sizeof(float) * ncands);
  spsassert(frm->dp->dpvals,"frm->dp->dpvals malloc failed in alloc_frame()");
    
  /*  Initialize the cumulative DP costs to zero */
  for(j = ncands-1; j >= 0; j--)
    frm->dp->dpvals[j] = 0.0;

  return(frm);
}


/*--------------------------------------------------------------------*/
/* push window stat to stack, and pop the oldest one */

int SC_Feature_Pitch::save_windstat(float *rho, int order, float err, float rms)
{
    int i,j;

    if(wReuse > 1){               /* push down the stack */
	for(j=1; j<wReuse; j++){
	    for(i=0;i<=order; i++) windstat[j-1].rho[i] = windstat[j].rho[i];
	    windstat[j-1].err = windstat[j].err;
	    windstat[j-1].rms = windstat[j].rms;
	}
	for(i=0;i<=order; i++) windstat[wReuse-1].rho[i] = rho[i]; /*save*/
	windstat[wReuse-1].err = err;
	windstat[wReuse-1].rms = rms;
	return 1;
    } else if (wReuse == 1) {
	for(i=0;i<=order; i++) windstat[0].rho[i] = rho[i];  /* save */
	windstat[0].err = err;
	windstat[0].rms = rms;
	return 1;
    } else 
	return 0;
}


/*--------------------------------------------------------------------*/
int SC_Feature_Pitch::retrieve_windstat(float *rho, int order, float *err, float *rms)
{
    Windstat wstat;
    int i;
	
    if(wReuse){
	wstat = windstat[0];
	for(i=0; i<=order; i++) rho[i] = wstat.rho[i];
	*err = wstat.err;
	*rms = wstat.rms;
	return 1;
    }
    else return 0;
}


/*--------------------------------------------------------------------*/
float SC_Feature_Pitch::get_similarity(int order, int size, float *pdata, float *cdata, float *rmsa, float *rms_ratio, float pre, float stab, int w_type, int init)
{
  float rho3[BIGSORD+1], err3, rms3, rmsd3, b0, t, a2[BIGSORD+1], 
      rho1[BIGSORD+1], a1[BIGSORD+1], b[BIGSORD+1], err1, rms1, rmsd1;
  //float itakura(), wind_energy();

/* (In the lpc() calls below, size-1 is used, since the windowing and
   preemphasis function assumes an extra point is available in the
   input data array.  This condition is apparently no longer met after
   Derek's modifications.) */

  /* get current window stat */
  lpc(order, stab, size-1, cdata,
      a2, rho3, (float *) NULL, &err3, &rmsd3, pre, w_type);
  rms3 = wind_energy(cdata, size, w_type);
  
  if(!init) {
      /* get previous window stat */
      if( !retrieve_windstat(rho1, order, &err1, &rms1)){
	  lpc(order, stab, size-1, pdata,
	      a1, rho1, (float *) NULL, &err1, &rmsd1, pre, w_type);
	  rms1 = wind_energy(pdata, size, w_type);
      }
      a_to_aca(a2+1,b,&b0,order);
      t = itakura(order,b,&b0,rho1+1,&err1) - (float)(.8);
      if(rms1 > 0.0)
	  *rms_ratio = (float)((0.001 + rms3)/rms1);
      else
	  if(rms3 > 0.0)
	      *rms_ratio = 2.0;	/* indicate some energy increase */
	  else
	      *rms_ratio = 1.0;	/* no change */
  } else {
      *rms_ratio = 1.0;
      t = 10.0;
  }
  *rmsa = rms3;
  save_windstat( rho3, order, err3, rms3);
  return((float)(0.2/t));
}


/* -------------------------------------------------------------------- */
/* This is an ad hoc signal stationarity function based on Itakura
 * distance and relative amplitudes.
 */
/* 
  This illustrates the window locations when the very first frame is read.
  It shows an example where each frame step |  .  | is 10 msec.  The
  frame step size is variable.  The window size is always 30 msec.
  The window centers '*' is always 20 msec apart.
  The windows cross each other right at the center of the DP frame, or
  where the '.' is.

                          ---------*---------   current window

              ---------*---------  previous window

  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |
              ^           ^  ^
              ^           ^  ^
              ^           ^  fdata
              ^           ^
              ^           q
	      p

                          ---
                          ind

  fdata, q, p, ind, are variables used below.
   
*/

SC_Feature_Pitch::Stat* SC_Feature_Pitch::get_stationarity(float *fdata, double freq, int buff_size, int nframes, int frame_step, int first_time)
{
  //static Stat *stat;
  //static int nframes_old = 0, memsize;
  //static float *mem;
  float preemp = (float)(0.4), stab = (float)(30.0);
  float *p, *q, *r, *datend;
  int ind, i, j, m, size, order, agap, w_type = 3;

  agap = (int) (STAT_AINT *freq);
  size = (int) (STAT_WSIZE * freq);
  ind = (agap - size) / 2;

  if( this->static_nframes_old < nframes || !stat || first_time){
    /* move this to init_dp_f0() later */
    this->static_nframes_old = nframes;
    if(stat){ 
      free((char *) stat->stat);
      free((char *) stat->rms);
      free((char *) stat->rms_ratio);
      free((char *) stat);
    }
    stat = (Stat *) malloc(nframes *sizeof(Stat));
    spsassert(stat,"stat malloc failed in get_stationarity");
    stat->stat = (float*)malloc(sizeof(float)*nframes);
    spsassert(stat->stat,"stat->stat malloc failed in get_stationarity");
    stat->rms = (float*)malloc(sizeof(float)*nframes);
    spsassert(stat->rms,"stat->rms malloc failed in get_stationarity");
    stat->rms_ratio = (float*)malloc(sizeof(float)*nframes);
    spsassert(stat->rms_ratio,"stat->ratio malloc failed in get_stationarity");
    this->static_memsize = (int) (STAT_WSIZE * freq) + (int) (STAT_AINT * freq);
    mem = (float *) malloc( sizeof(float) * this->static_memsize);
    spsassert(mem, "mem malloc failed in get_stationarity()");
    for(j=0; j<this->static_memsize; j++) mem[j] = 0;
  }
  
  if(nframes == 0) return(stat);

  q = fdata + ind;
  datend = fdata + buff_size;

  if((order = (int)(2.0 + (freq/1000.0))) > BIGSORD) {
    Fprintf(stderr,
	    "%s: Optimim order (%d) exceeds that allowable (%d); reduce Fs\n",
	    ProgName, order, BIGSORD);
    order = BIGSORD;
  }

  /* prepare for the first frame */
  for(j=this->static_memsize/2, i=0; j<this->static_memsize; j++, i++) mem[j] = fdata[i];

  /* never run over end of frame, should already taken care of when read */

  for(j=0, p = q - agap; j < nframes; j++, p += frame_step, q += frame_step){
      if( (p >= fdata) && (q >= fdata) && ( q + size <= datend) )
	  stat->stat[j] = get_similarity(order,size, p, q, 
					     &(stat->rms[j]),
					     &(stat->rms_ratio[j]),preemp,
					     stab,w_type, 0);
      else {
	  if(first_time) {
	      if( (p < fdata) && (q >= fdata) && (q+size <=datend) )
		  stat->stat[j] = get_similarity(order,size, NULL, q,
						     &(stat->rms[j]),
						     &(stat->rms_ratio[j]),
						     preemp,stab,w_type, 1);
	      else{
		  stat->rms[j] = 0.0;
		  stat->stat[j] = (float)(0.01 * 0.2);   /* a big transition */
		  stat->rms_ratio[j] = 1.0;   /* no amplitude change */
	      }
	  } else {
	      if( (p<fdata) && (q+size <=datend) ){
		  stat->stat[j] = get_similarity(order,size, mem, 
						     mem + (this->static_memsize/2) + ind,
						     &(stat->rms[j]),
						     &(stat->rms_ratio[j]),
						     preemp, stab,w_type, 0);
		  /* prepare for the next frame_step if needed */
		  if(p + frame_step < fdata ){
		      for( m=0; m<(this->static_memsize-frame_step); m++) 
			  mem[m] = mem[m+frame_step];
		      r = q + size;
		      for( m=0; m<frame_step; m++) 
			  mem[this->static_memsize-frame_step+m] = *r++;
		  }
	      }
	  }
      }
  }

  /* last frame, prepare for next call */
  for(j=(this->static_memsize/2)-1, p=fdata + (nframes * frame_step)-1; j>=0 && p>=fdata; j-- ) //by thilo: "&& p>=fdata" from snack code
    mem[j] = *p--;
  return(stat);
}


/* -------------------------------------------------------------------- */
/*	Round the argument to the nearest integer.			*/

int SC_Feature_Pitch::round(double flnum)
{
  return((flnum >= 0.0) ? (int)(flnum + 0.5) : (int)(flnum - 0.5));
}

//static char *sccs_id = "@(#)get_cands.c	1.5	9/9/96	ERL"; //get_cands.c

#define TRUE 1
#define FALSE 0

//static void get_cand(), peak(), do_ffir();
//static int lc_lin_fir(), downsamp();

/* ----------------------------------------------------------------------- */
void SC_Feature_Pitch::get_fast_cands(float *fdata, float *fdsdata, int ind, int step, int size, int dec, int start, int nlags, float *engref, int *maxloc, float *maxval, Cross *cp, float *peaks, int *locs, int *ncand, F0_params *par)
{
  int decind, decstart, decnlags, decsize, i, j, *lp;
  float *corp, xp, yp, lag_wt;
  register float *pe;

  lag_wt = par->lag_weight/nlags;
  decnlags = 1 + (nlags/dec);
  if((decstart = start/dec) < 1) decstart = 1;
  decind = (ind * step)/dec;
  decsize = 1 + (size/dec);
  corp = cp->correl;
    
  crossf(fdsdata + decind, decsize, decstart, decnlags, engref, maxloc,
	maxval, corp);
  cp->maxloc = *maxloc;	/* location of maximum in correlation */
  cp->maxval = *maxval;	/* max. correlation value (found at maxloc) */
  cp->rms = sqrt(*engref/size); /* rms in reference window */
  cp->firstlag = decstart;

  get_cand(cp,peaks,locs,decnlags,ncand,par->cand_thresh); /* return high peaks in xcorr */

  /* Interpolate to estimate peak locations and values at high sample rate. */
  for(i = *ncand, lp = locs, pe = peaks; i--; pe++, lp++) {
    j = *lp - decstart - 1;
    peak(&corp[j],&xp,&yp);
    *lp = (*lp * dec) + (int)(0.5+(xp*dec)); /* refined lag */
    *pe = (float)(yp*(1.0 - (lag_wt* *lp))); /* refined amplitude */
  }
  
  if(*ncand >= par->n_cands) {	/* need to prune candidates? */
    register int *loc, *locm, lt;
    register float smaxval, *pem;
    register int outer, inner, lim;
    for(outer=0, lim = par->n_cands-1; outer < lim; outer++)
      for(inner = *ncand - 1 - outer,
	  pe = peaks + (*ncand) -1, pem = pe-1,
	  loc = locs + (*ncand) - 1, locm = loc-1;
	  inner--;
	  pe--,pem--,loc--,locm--)
	if((smaxval = *pe) > *pem) {
	  *pe = *pem;
	  *pem = smaxval;
	  lt = *loc;
	  *loc = *locm;
	  *locm = lt;
	}
    *ncand = par->n_cands-1;  /* leave room for the unvoiced hypothesis */
  }
  crossfi(fdata + (ind * step), size, start, nlags, 7, engref, maxloc,
	  maxval, corp, locs, *ncand);

  cp->maxloc = *maxloc;	/* location of maximum in correlation */
  cp->maxval = *maxval;	/* max. correlation value (found at maxloc) */
  cp->rms = sqrt(*engref/size); /* rms in reference window */
  cp->firstlag = start;
  get_cand(cp,peaks,locs,nlags,ncand,par->cand_thresh); /* return high peaks in xcorr */
    if(*ncand >= par->n_cands) {	/* need to prune candidates again? */
    register int *loc, *locm, lt;
    register float smaxval, *pe, *pem;
    register int outer, inner, lim;
    for(outer=0, lim = par->n_cands-1; outer < lim; outer++)
      for(inner = *ncand - 1 - outer,
	  pe = peaks + (*ncand) -1, pem = pe-1,
	  loc = locs + (*ncand) - 1, locm = loc-1;
	  inner--;
	  pe--,pem--,loc--,locm--)
	if((smaxval = *pe) > *pem) {
	  *pe = *pem;
	  *pem = smaxval;
	  lt = *loc;
	  *loc = *locm;
	  *locm = lt;
	}
    *ncand = par->n_cands - 1;  /* leave room for the unvoiced hypothesis */
  }
}

/* ----------------------------------------------------------------------- */
float* SC_Feature_Pitch::downsample(float *input, int samsin, int state_idx, double freq, int *samsout, int decimate, int first_time, int last_time)
{
  //static float	b[2048];
  //static float *foutput;
  float	beta = 0.0;
  //static int	ncoeff = 127, ncoefft = 0;
  int init;

  if(input && (samsin > 0) && (decimate > 0) && *samsout) {
    if(decimate == 1) {
      return(input);
    }

    if(first_time){
      int nbuff = (samsin/decimate) + (2*this->static_ncoeff);

      this->static_ncoeff = ((int)(freq * .005)) | 1;
      beta = (float)(.5/decimate);
      this->static_foutput = (float*)malloc(sizeof(float) * nbuff);
      spsassert(this->static_foutput, "Can't allocate foutput in downsample");
      for( ; nbuff > 0 ;)
	this->static_foutput[--nbuff] = 0.0;

      if( !lc_lin_fir(beta,&this->static_ncoeff,this->static_b)) {
	fprintf(stderr,"\nProblems computing interpolation filter\n");
	free(this->static_foutput);
	return(NULL);
      }
      this->static_ncoefft = (this->static_ncoeff/2) + 1;
    }		    /*  endif new coefficients need to be computed */

    if(first_time) init = 1;
    else if (last_time) init = 2;
    else init = 0;
    
    if(downsamp(input,this->static_foutput,samsin,samsout,state_idx,decimate,this->static_ncoefft,this->static_b,init)) {
      return(this->static_foutput);
    } else
      Fprintf(stderr,"Problems in downsamp() in downsample()\n");
  } else
    Fprintf(stderr,"Bad parameters passed to downsample()\n");
  
  return(NULL);
}

/* ----------------------------------------------------------------------- */
/* Get likely candidates for F0 peaks. */
void SC_Feature_Pitch::get_cand(Cross *cross, float *peak, int *loc, int nlags, int *ncand, float cand_thresh)
{
  register int i, lastl, *t;
  register float o, p, q, *r, *s, clip;
  int start, ncan, maxl;

  clip = cand_thresh * cross->maxval;
  maxl = cross->maxloc;
  lastl = nlags - 2;
  start = cross->firstlag;

  r = cross->correl;
  o= *r++;			/* first point */
  q = *r++;	                /* middle point */
  p = *r++;
  s = peak;
  t = loc;
  ncan=0;
  for(i=1; i < lastl; i++, o=q, q=p, p= *r++){
    if((q > clip) &&		/* is this a high enough value? */
      (q >= p) && (q >= o)){ /* NOTE: this finds SHOLDERS and PLATEAUS
				      as well as peaks (is this a good idea?) */
	*s++ = q;		/* record the peak value */
	*t++ = i + start;	/* and its location */
	ncan++;			/* count number of peaks found */
      }
  }
/*
  o = q;
  q = p;
  if( (q > clip) && (q >=0)){
    *s++ = q;
    *t++ = i+start;
    ncan++;
  }
*/
  *ncand = ncan;
}

/* ----------------------------------------------------------------------- */
/* buffer-to-buffer downsample operation */
/* This is STRICTLY a decimator! (no upsample) */
int SC_Feature_Pitch::downsamp(float *in, float *out, int samples, int *outsamps, int state_idx, int decimate, int ncoef, float fc[], int init)
{
  if(in && out) {
    do_ffir(in, samples, out, outsamps, state_idx, ncoef, fc, 0, decimate, init);
    return(TRUE);
  } else
    printf("Bad signal(s) passed to downsamp()\n");
  return(FALSE);
}

/*      ----------------------------------------------------------      */
void SC_Feature_Pitch::do_ffir(register float *buf, register int in_samps, register float *bufo, register int *out_samps, int idx, register int ncoef, float *fc, register int invert, register int skip, register int init)
/* fc contains 1/2 the coefficients of a symmetric FIR filter with unity
    passband gain.  This filter is convolved with the signal in buf.
    The output is placed in buf2.  If(invert), the filter magnitude
    response will be inverted.  If(init&1), beginning of signal is in buf;
    if(init&2), end of signal is in buf.  out_samps is set to the number of
    output points placed in bufo. */
{
  register float *dp1, *dp2, *dp3, sum, integral;
  //static float *co=NULL, *mem=NULL;
  //static float state[1000];
  //static int fsize=0, resid=0;
  register int i, j, k, l;
  register float *sp;
  register float *buf1;

  buf1 = buf;
  if(ncoef > this->static_fsize) {/*allocate memory for full coeff. array and filter memory */
    if(this->static_co)
      free(this->static_co);
    if(this->static_mem)
      free(this->static_mem);
    this->static_fsize = 0;
    i = (ncoef+1)*2;
    if(!((this->static_co = (float *)malloc(sizeof(float)*i)) &&
	 (this->static_mem = (float *)malloc(sizeof(float)*i)))) {
      fprintf(stderr,"allocation problems in do_fir()\n");
      exit(-1);
    }
    this->static_fsize = ncoef;
  }

  /* fill 2nd half with data */
  for(i=ncoef, dp1=this->static_mem+ncoef-1; i-- > 0; )  *dp1++ = *buf++;  

  if(init & 1) {	/* Is the beginning of the signal in buf? */
    /* Copy the half-filter and its mirror image into the coefficient array. */
    for(i=ncoef-1, dp3=fc+ncoef-1, dp2=this->static_co, dp1 = this->static_co+((ncoef-1)*2),
	integral = 0.0; i-- > 0; )
      if(!invert) *dp1-- = *dp2++ = *dp3--;
      else {
	integral += (sum = *dp3--);
	*dp1-- = *dp2++ = -sum;
      }
    if(!invert)  *dp1 = *dp3;	/* point of symmetry */
    else {
      integral *= 2;
      integral += *dp3;
      *dp1 = integral - *dp3;
    }

    for(i=ncoef-1, dp1=this->static_mem; i-- > 0; ) *dp1++ = 0;
  }
  else
    for(i=ncoef-1, dp1=this->static_mem, sp=this->static_state; i-- > 0; ) *dp1++ = *sp++;

  i = in_samps;
  this->static_resid = 0;

  k = (ncoef << 1) -1;	/* inner-product loop limit */

  if(skip <= 1) {       /* never used */
/*    *out_samps = i;	
    for( ; i-- > 0; ) {	
      for(j=k, dp1=this->static_mem, dp2=this->static_co, dp3=this->static_mem+1, sum = 0.0; j-- > 0;
	  *dp1++ = *dp3++ )
	sum += *dp2++ * *dp1;

      *--dp1 = *buf++;	
      *bufo++ = (sum < 0.0)? sum -0.5 : sum +0.5; 
    }
    if(init & 2) {	
      for(i=ncoef; i-- > 0; ) {
	for(j=k, dp1=this->static_mem, dp2=this->static_co, dp3=this->static_mem+1, sum = 0.0; j-- > 0;
	    *dp1++ = *dp3++ )
	  sum += *dp2++ * *dp1;
	*--dp1 = 0.0;
	*bufo++ = (sum < 0)? sum -0.5 : sum +0.5; 
      }
      *out_samps += ncoef;
    }
    return;
*/
  } 
  else {			/* skip points (e.g. for downsampling) */
    /* the buffer end is padded with (ncoef-1) data points */
    for( l=0 ; l < *out_samps; l++ ) {
      for(j=k-skip, dp1=this->static_mem, dp2=this->static_co, dp3=this->static_mem+skip, sum=0.0; j-- >0;
	  *dp1++ = *dp3++)
	sum += *dp2++ * *dp1;
      for(j=skip; j-- >0; *dp1++ = *buf++) /* new data to memory */
	sum += *dp2++ * *dp1;
      *bufo++ = (float)((sum<0.0) ? sum -0.5 : sum +0.5);
    }
    if(init & 2){
      this->static_resid = in_samps - *out_samps * skip;
      for(l=this->static_resid/skip; l-- >0; ){
	for(j=k-skip, dp1=this->static_mem, dp2=this->static_co, dp3=this->static_mem+skip, sum=0.0; j-- >0;
	    *dp1++ = *dp3++)
	    sum += *dp2++ * *dp1;
	for(j=skip; j-- >0; *dp1++ = 0.0)
	  sum += *dp2++ * *dp1;
	*bufo++ = (float)((sum<0.0) ? sum -0.5 : sum +0.5);
	(*out_samps)++;
      }
    }
    else
      for(dp3=buf1+idx-ncoef+1, l=ncoef-1, sp=this->static_state; l-- >0; ) *sp++ = *dp3++;
  }
}

/*      ----------------------------------------------------------      */
int SC_Feature_Pitch::lc_lin_fir(register float fc, int *nf, float *coef)
/* create the coefficients for a symmetric FIR lowpass filter using the
   window technique with a Hanning window. */
{
    register int	i, n;
    register double	twopi, fn, c;

    if(((*nf % 2) != 1))
	*nf = *nf + 1;
    n = (*nf + 1)/2;

    /*  Compute part of the ideal impulse response (the sin(x)/x kernel). */
		twopi = sclib::pi * 2.0;
    coef[0] = (float)(2.0 * fc);
    c = sclib::pi;
    fn = twopi * fc;
    for(i=1;i < n; i++) coef[i] = (float)(sin(i * fn)/(c * i));

    /* Now apply a Hanning window to the (infinite) impulse response. */
    /* (Probably should use a better window, like Kaiser...) */
    fn = twopi/(double)(*nf);
    for(i=0;i<n;i++) 
	coef[n-i-1] *= (float)(.5 - (.5 * cos(fn * ((double)i + 0.5))));
    
    return(TRUE);
}


/* ----------------------------------------------------------------------- */
/* Use parabolic interpolation over the three points defining the peak
 * vicinity to estimate the "true" peak. */
void SC_Feature_Pitch::peak(float *y, float *xp, float *yp) /* y: vector of length 3 defining peak;  x,y: values of parabolic peak fitting the input points. */
{
  register float a, c;
  
  a = (y[2]-y[1])+(float)(.5*(y[0]-y[2]));
  if(fabs(a) > .000001) {
    *xp = c = (y[0]-y[2])/(float)(4.0*a);
    *yp = y[1] - (a*c*c);
  } else {
    *xp = 0.0;
    *yp = y[1];
  }
}

//static char *sccs_id = "@(#)sigproc.c	1.4	9/9/96	ERL";

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Return a time-weighting window of type type and length n in dout.
 * Dout is assumed to be at least n elements long.  Type is decoded in
 * the switch statement below.
 */
int SC_Feature_Pitch::get_window(register float *dout, register int n, register int type)
{
  //static float *din = NULL;
  //static int n0 = 0;
  float preemp = 0.0;

  if(n > this->static_n0) {
    register float *p;
    register int i;
    
    if(this->static_din) free(this->static_din);
    this->static_din = NULL;
    if(!(this->static_din = (float*)malloc(sizeof(float)*n))) {
      Fprintf(stderr,"Allocation problems in get_window()\n");
      return(FALSE);
    }
    for(i=0, p=this->static_din; i++ < n; ) *p++ = 1;
    this->static_n0 = n;
  }
  return(window(this->static_din, dout, n, preemp, type));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Apply a rectangular window (i.e. none).  Optionally, preemphasize. */
void SC_Feature_Pitch::rwindow(register float *din, register float *dout, register int n, register float preemp)
{
  register float *p;
 
  // If preemphasis is to be performed,  this assumes that there are n+1 valid  samples in the input buffer (din).
  if(preemp != 0.0) {
    for( p=din+1; n-- > 0; )
      *dout++ = (float)(*p++) - (preemp * *din++);
  } else {
    for( ; n-- > 0; )
      *dout++ =  *din++;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generate a cos^4 window, if one does not already exist. */
void SC_Feature_Pitch::cwindow(register float *din, register float *dout, register int n, register float preemp)
{
  register int i;
  register float *p;
  //static int wsize = 0;
  //static float *wind=NULL;
  register float *q, co;
 
  if(this->static_wsize != n) {		//Need to create a new cos**4 window?
    register double arg, half=0.5;
    
    if(this->static_wind) this->static_wind = (float*)realloc(this->static_wind,n*sizeof(float));
    else this->static_wind = (float*)malloc(n*sizeof(float));
    this->static_wsize = n;
    for(i=0, arg=3.1415927*2.0/(this->static_wsize), q=this->static_wind; i < n; ) {
      co = (float)(half*(1.0 - cos((half + (double)i++) * arg)));
      *q++ = co * co * co * co;
    }
  }
	//If preemphasis is to be performed,  this assumes that there are n+1 valid samples in the input buffer (din).
  if(preemp != 0.0) {
    for(i=n, p=din+1, q=this->static_wind; i--; )
      *dout++ = *q++ * ((float)(*p++) - (preemp * *din++));
  } else {
    for(i=n, q=this->static_wind; i--; )
      *dout++ = *q++ * *din++;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generate a Hamming window, if one does not already exist. */
void SC_Feature_Pitch::hwindow(register float *din, register float *dout, register int n, register float preemp)
{
  register int i;
  register float *p;
  //static int wsize = 0;
  //static float *wind=NULL;
  register float *q;

  if(this->static_wsize_hwindow != n) {		//Need to create a new Hamming window?
    register double arg, half=0.5;
    
    if(this->static_wind_hwindow) this->static_wind_hwindow = (float*)realloc(this->static_wind_hwindow,n*sizeof(float));
    else this->static_wind_hwindow = (float*)malloc(n*sizeof(float));
    this->static_wsize_hwindow = n;
    for(i=0, arg=3.1415927*2.0/(this->static_wsize_hwindow), q=this->static_wind_hwindow; i < n; )
      *q++ = (float)(.54 - .46 * cos((half + (double)i++) * arg));
  }
// If preemphasis is to be performed,  this assumes that there are n+1 valid samples in the input buffer (din).
  if(preemp != 0.0) {
    for(i=n, p=din+1, q=this->static_wind_hwindow; i--; )
      *dout++ = *q++ * ((float)(*p++) - (preemp * *din++));
  } else {
    for(i=n, q=this->static_wind_hwindow; i--; )
      *dout++ = *q++ * *din++;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generate a Hanning window, if one does not already exist. */
void SC_Feature_Pitch::hnwindow(register float *din, register float *dout, register int n, register float preemp)
{
  register int i;
  register float *p;
  //static int wsize = 0;
  //static float *wind=NULL;
  register float *q;

  if(this->static_wsize_hnwindow != n) {		//Need to create a new Hanning window?
    register double arg, half=0.5;
    
    if(this->static_wind_hnwindow) this->static_wind_hnwindow = (float*)realloc(this->static_wind_hnwindow,n*sizeof(float));
    else this->static_wind_hnwindow = (float*)malloc(n*sizeof(float));
    this->static_wsize_hnwindow = n;
    for(i=0, arg=3.1415927*2.0/(this->static_wsize_hnwindow), q=this->static_wind_hnwindow; i < n; )
      *q++ = (float)(half - half * cos((half + (double)i++) * arg));
  }
  //If preemphasis is to be performed,  this assumes that there are n+1 valid samples in the input buffer (din).
  if(preemp != 0.0) {
    for(i=n, p=din+1, q=this->static_wind_hnwindow; i--; )
      *dout++ = *q++ * ((float)(*p++) - (preemp * *din++));
  } else {
    for(i=n, q=this->static_wind_hnwindow; i--; )
      *dout++ = *q++ * *din++;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Apply a window of type type to the short PCM sequence of length n
 * in din.  Return the floating-point result sequence in dout.  If preemp
 * is non-zero, apply preemphasis to tha data as it is windowed.
 */
int SC_Feature_Pitch::window(register float *din, register float *dout, register int n, register float preemp, int type)
{
  switch(type) {
  case 0:			// rectangular
    rwindow(din, dout, n, preemp);
    break;
  case 1:			// Hamming
    hwindow(din, dout, n, preemp);
    break;
  case 2:			// cos^4
    cwindow(din, dout, n, preemp);
    break;
  case 3:			// Hanning
    hnwindow(din, dout, n, preemp);
    break;
  default:
    Fprintf(stderr,"Unknown window type (%d) requested in window()\n",type);
    return(FALSE);
  }
  return(TRUE);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Compute the pp+1 autocorrelation lags of the windowsize samples in s.
 * Return the normalized autocorrelation coefficients in r.
 * The rms is returned in e.
 */
void SC_Feature_Pitch::autoc(register int windowsize, register float *s, register int p, register float *r, register float *e)
{
  register int i, j;
  register float *q, *t, sum, sum0;

  for( i=windowsize, q=s, sum0=0.0; i--;) {
    sum = *q++;
    sum0 += sum*sum;
  }
  *r = 1.;			// r[0] will always =1.
  if(sum0 == 0.0) {		// No energy: fake low-energy white noise.
    *e = 1.;			// Arbitrarily assign 1 to rms.
    // Now fake autocorrelation of white noise.
    for ( i=1; i<=p; i++){
      r[i] = 0.;
    }
    return;
  }
  *e = (float)(sqrt((double)(sum0/windowsize)));
  sum0 = (float)(1.0/sum0);
  for( i=1; i <= p; i++){
    for( sum=0.0, j=windowsize-i, q=s, t=s+i; j--; )
      sum += (*q++) * (*t++);
    *(++r) = sum*sum0;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Using Durbin's recursion, convert the autocorrelation sequence in r
 * to reflection coefficients in k and predictor coefficients in a.
 * The prediction error energy (gain) is left in *ex.
 * Note: durbin returns the coefficients in normal sign format.
 *	(i.e. a[0] is assumed to be = +1.)
 */
void SC_Feature_Pitch::durbin(register float *r, register float *k, register float *a, register int p, register float *ex) //p: analysis order
{
  float  bb[BIGSORD];
  register int i, j;
  register float e, s, *b = bb;

  e = *r;
  *k = -r[1]/e;
  *a = *k;
  e *= (float)((1. - (*k) * (*k)));
  for ( i=1; i < p; i++){
    s = 0;
    for ( j=0; j<i; j++){
      s -= a[j] * r[i-j];
    }
    k[i] = ( s - r[i+1] )/e;
    a[i] = k[i];
    for ( j=0; j<=i; j++){
      b[j] = a[j];
    }
    for ( j=0; j<i; j++){
      a[j] += k[i] * b[i-j-1];
    }
    e *= (float)(( 1. - (k[i] * k[i]) ));
  }
  *ex = e;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Compute the autocorrelations of the p LP coefficients in a. 
 *  (a[0] is assumed to be = 1 and not explicitely accessed.)
 *  The magnitude of a is returned in c.
 *  2* the other autocorrelation coefficients are returned in b.
 */
void SC_Feature_Pitch::a_to_aca (float *a, float *b, float *c, register int p)
{
  register float  s, *ap, *a0;
  register int  i, j;

  for ( s=1., ap=a, i = p; i--; ap++ )
    s += *ap * *ap;

  *c = s;
  for ( i = 1; i <= p; i++){
    s = a[i-1];
    for (a0 = a, ap = a+i, j = p-i; j--; )
      s += (*a0++ * *ap++);
    *b++ = (float)(2. * s);
  }

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Compute the Itakura LPC distance between the model represented
 * by the signal autocorrelation (r) and its residual (gain) and
 * the model represented by an LPC autocorrelation (c, b).
 * Both models are of order p.
 * r is assumed normalized and r[0]=1 is not explicitely accessed.
 * Values returned by the function are >= 1.
 */
float SC_Feature_Pitch::itakura (register int p, register float *b, register float *c, register float *r, register float *gain)
{
  register float s;

  for( s= *c; p--; )
    s += *r++ * *b++;

  return (s/ *gain);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Compute the time-weighted RMS of a size segment of data.  The data
 * is weighted by a window of type w_type before RMS computation.  w_type
 * is decoded above in window().
 */
float SC_Feature_Pitch::wind_energy(register float *data, register int size, register int w_type) /* data: input PCM data; size: size of window; w_type: window type */
{
  //static int nwind = 0;
  //static float *dwind = NULL;
  register float *dp, sum, f;
  register int i;

  if(this->static_nwind < size) {
    if(this->static_dwind) this->static_dwind = (float*)realloc(this->static_dwind,size*sizeof(float));
    else this->static_dwind = (float*)malloc(size*sizeof(float));
    if(!this->static_dwind) {
      Fprintf(stderr,"Can't allocate scratch memory in wind_energy()\n");
      return(0.0);
    }
  }
  if(this->static_nwind != size) {
    get_window(this->static_dwind, size, w_type);
    this->static_nwind = size;
  }
  for(i=size, dp = this->static_dwind, sum = 0.0; i-- > 0; ) {
    f = *dp++ * (float)(*data++);
    sum += f*f;
  }
  return((float)sqrt((double)(sum/size)));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generic autocorrelation LPC analysis of the short-integer data
 * sequence in data.
 */
int SC_Feature_Pitch::lpc(int lpc_ord, float lpc_stabl, int wsize, float *data, float *lpca, float *ar, float *lpck, float *normerr, float *rms, float preemp,int type)
			 //int lpc_ord,		/* Analysis order */
       //wsize,			/* window size in points */
       //type;		/* window type (decoded in window() above) */
       //float lpc_stabl,	/* Stability factor to prevent numerical problems. */
       //*lpca,		/* if non-NULL, return vvector for predictors */
       //*ar,		/* if non-NULL, return vector for normalized autoc. */
       //*lpck,		/* if non-NULL, return vector for PARCOR's */
       //*normerr,		/* return scaler for normalized error */
       //*rms,		/* return scaler for energy in preemphasized window */
       //preemp;
       //float *data;	/* input data sequence; assumed to be wsize+1 long */
{
  //static float *dwind=NULL;
  //static int nwind=0;
  float rho[BIGSORD+1], k[BIGSORD], a[BIGSORD+1],*r,*kp,*ap,en,er,wfact = (float)(1.0); //by thilo: initilization of wfact as in snack

  if((wsize <= 0) || (!data) || (lpc_ord > BIGSORD)) return(FALSE);
  
  if(this->static_nwind_lpc != wsize) {
    if(this->static_dwind_lpc) this->static_dwind_lpc = (float*)realloc(this->static_dwind_lpc,wsize*sizeof(float));
    else this->static_dwind_lpc = (float*)malloc(wsize*sizeof(float));
    if(!this->static_dwind_lpc) {
      Fprintf(stderr,"Can't allocate scratch memory in lpc()\n");
      return(FALSE);
    }
    this->static_nwind_lpc = wsize;
  }
  
  window(data, this->static_dwind_lpc, wsize, preemp, type);
  if(!(r = ar)) r = rho;	/* Permit optional return of the various */
  if(!(kp = lpck)) kp = k;	/* coefficients and intermediate results. */
  if(!(ap = lpca)) ap = a;
  autoc( wsize, this->static_dwind_lpc, lpc_ord, r, &en );
  if(lpc_stabl > 1.0) {	/* add a little to the diagonal for stability */
    int i;
    float ffact;
    ffact =(float)(1.0/(1.0 + exp((-lpc_stabl/20.0) * log(10.0))));
    for(i=1; i <= lpc_ord; i++) rho[i] = ffact * r[i];
    *rho = *r;
    r = rho;
    if(ar)
      for(i=0;i<=lpc_ord; i++) ar[i] = r[i];
  }
  durbin ( r, kp, &ap[1], lpc_ord, &er);
  switch(type) {		/* rms correction for window */
  case 0:
    wfact = 1.0;		/* rectangular */
    break;
  case 1:
    wfact = (float)(.630397);		/* Hamming */
    break;
  case 2:
    wfact = (float)(.443149);		/* (.5 - .5*cos)^4 */
    break;
  case 3:
    wfact = (float)(.612372);		/* Hanning */
    break;
  }
  *ap = 1.0;
  if(rms) *rms = en/wfact;
  if(normerr) *normerr = er;
  return(TRUE);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Return a sequence based on the normalized crosscorrelation of the signal
   in data.
 *
  data is the input speech array
  size is the number of samples in each correlation
  start is the first lag to compute (governed by the highest expected F0)
  nlags is the number of cross correlations to compute (set by lowest F0)
  engref is the energy computed at lag=0 (i.e. energy in ref. window)
  maxloc is the lag at which the maximum in the correlation was found
  maxval is the value of the maximum in the CCF over the requested lag interval
  correl is the array of nlags cross-correlation coefficients (-1.0 to 1.0)
 *
 */
void SC_Feature_Pitch::crossf(float *data, int size, int start, int nlags, float *engref, int *maxloc, float *maxval, float *correl)
{
  //static float *dbdata=NULL;
  //static int dbsize = 0;
  register float *dp, *ds, sum, st;
  register int j;
  register  float *dq, t, *p, engr, *dds, amax;
  register  double engc;
  int i, iloc, total;
  int sizei, sizeo, maxsize;

  /* Compute mean in reference window and subtract this from the
     entire sequence.  This doesn't do too much damage to the data
     sequenced for the purposes of F0 estimation and removes the need for
     more principled (and costly) low-cut filtering. */
  if((total = size+start+nlags) > this->static_dbsize) {
    if(this->static_dbdata)
      free(this->static_dbdata);
    this->static_dbdata = NULL;
    this->static_dbsize = 0;
    if(!(this->static_dbdata = (float*)malloc(sizeof(float)*total))) {
      Fprintf(stderr,"Allocation failure in crossf()\n");
      exit(-1);
    }
    this->static_dbsize = total;
  }
  for(engr=0.0, j=size, p=data; j--; ) engr += *p++;
  engr /= size;
  for(j=size+nlags+start, dq = this->static_dbdata, p=data; j--; )  *dq++ = *p++ - engr;

  maxsize = start + nlags;
  sizei = size + start + nlags + 1;
  sizeo = nlags + 1;
 
  /* Compute energy in reference window. */
  for(j=size, dp=this->static_dbdata, sum=0.0; j--; ) {
    st = *dp++;
    sum += st * st;
  }

  *engref = engr = sum;
  if(engr > 0.0) {    /* If there is any signal energy to work with... */
    /* Compute energy at the first requested lag. */  
    for(j=size, dp=this->static_dbdata+start, sum=0.0; j--; ) {
      st = *dp++;
      sum += st * st;
    }
    engc = sum;

    /* COMPUTE CORRELATIONS AT ALL OTHER REQUESTED LAGS. */
    for(i=0, dq=correl, amax=0.0, iloc = -1; i < nlags; i++) {
      for(j=size, sum=0.0, dp=this->static_dbdata, dds = ds = this->static_dbdata+i+start; j--; )
	sum += *dp++ * *ds++;
      *dq++ = t = (float)(sum/sqrt((double)(engc*engr))); /* output norm. CC */
      engc -= (double)(*dds * *dds); /* adjust norm. energy for next lag */
      if((engc += (double)(*ds * *ds)) < 1.0)
	engc = 1.0;		/* (hack: in case of roundoff error) */
      if(t > amax) {		/* Find abs. max. as we go. */
	amax = t;
	iloc = i+start;
      }
    }
    *maxloc = iloc;
    *maxval = amax;
  } else {	/* No energy in signal; fake reasonable return vals. */
    *maxloc = 0;
    *maxval = 0.0;
    for(p=correl,i=nlags; i-- > 0; )
      *p++ = 0.0;
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Return a sequence based on the normalized crosscorrelation of the
   signal in data.  This is similar to crossf(), but is designed to
   compute only small patches of the correlation sequence.  The length of
   each patch is determined by nlags; the number of patches by nlocs, and
   the locations of the patches is specified by the array locs.  Regions
   of the CCF that are not computed are set to 0. 
 *
  data is the input speech array
  size is the number of samples in each correlation
  start0 is the first (virtual) lag to compute (governed by highest F0)
  nlags0 is the number of lags (virtual+actual) in the correlation sequence
  nlags is the number of cross correlations to compute at each location
  engref is the energy computed at lag=0 (i.e. energy in ref. window)
  maxloc is the lag at which the maximum in the correlation was found
  maxval is the value of the maximum in the CCF over the requested lag interval
  correl is the array of nlags cross-correlation coefficients (-1.0 to 1.0)
  locs is an array of indices pointing to the center of a patches where the
       cross correlation is to be computed.
  nlocs is the number of correlation patches to compute.
 *
 */
void SC_Feature_Pitch::crossfi(float *data, int size, int start0, int nlags0, int nlags, float *engref, int *maxloc, float *maxval, float *correl, int *locs, int nlocs)
{
  //static float *dbdata=NULL;
  //static int dbsize = 0;
  register float *dp, *ds, sum, st;
  register int j;
  register  float *dq, t, *p, engr, *dds, amax;
  register  double engc;
  int i, iloc, start, total;

  /* Compute mean in reference window and subtract this from the
     entire sequence. */
  if((total = size+start0+nlags0) > this->static_dbsize_crossfi) {
    if(this->static_dbdata_crossfi)
      free(this->static_dbdata_crossfi);
    this->static_dbdata_crossfi = NULL;
    this->static_dbsize_crossfi = 0;
    if(!(this->static_dbdata_crossfi = (float*)malloc(sizeof(float)*total))) {
      Fprintf(stderr,"Allocation failure in crossfi()\n");
      exit(-1);
    }
    this->static_dbsize_crossfi = total;
  }
  for(engr=0.0, j=size, p=data; j--; ) engr += *p++;
  engr /= size;
/*  for(j=size+nlags0+start0, t = -2.1, amax = 2.1, dq = this->static_dbdata_crossfi, p=data; j--; ) {
    if(((smax = *p++ - engr) > t) && (smax < amax))
      smax = 0.0;
    *dq++ = smax;
  } */
  for(j=size+nlags0+start0, dq = this->static_dbdata_crossfi, p=data; j--; ) {
    *dq++ = *p++ - engr;
  }

  /* Zero the correlation output array to avoid confusing the peak
     picker (since all lags will not be computed). */
  for(p=correl,i=nlags0; i-- > 0; )
    *p++ = 0.0;

  /* compute energy in reference window */
  for(j=size, dp=this->static_dbdata_crossfi, sum=0.0; j--; ) {
    st = *dp++;
    sum += st * st;
  }

  *engref = engr = sum;
   amax=0.0;
  iloc = -1;
  if(engr > 0.0) {
    for( ; nlocs > 0; nlocs--, locs++ ) {
      start = *locs - (nlags>>1);
      if(start < start0)
	start = start0;
      dq = correl + start - start0;
      /* compute energy at first requested lag */  
      for(j=size, dp=this->static_dbdata_crossfi+start, sum=0.0; j--; ) {
	st = *dp++;
	sum += st * st;
      }
      engc = sum;

      /* COMPUTE CORRELATIONS AT ALL REQUESTED LAGS */
      for(i=0; i < nlags; i++) {
	for(j=size, sum=0.0, dp=this->static_dbdata_crossfi, dds = ds = this->static_dbdata_crossfi+i+start; j--; )
	  sum += *dp++ * *ds++;
	if(engc < 1.0)
	  engc = 1.0;		/* in case of roundoff error */
	*dq++ = t = (float)(sum/sqrt((double)(10000.0 + (engc*engr))));
	engc -= (double)(*dds * *dds);
	engc += (double)(*ds * *ds);
	if(t > amax) {
	  amax = t;
	  iloc = i+start;
	}
      }
    }
    *maxloc = iloc;
    *maxval = amax;
  } else {
    *maxloc = 0;
    *maxval = 0.0;
  }
}
#endif
