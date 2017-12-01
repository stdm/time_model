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

#ifndef __SC_Feature_Pitch_H__
#define __SC_Feature_Pitch_H__

#include "SC_TweakableParameters.h"
#include "SC_Api.h"
#include <SV_Feature.h>

class SC_Feature_Pitch : public SV_Feature {

	private :

	protected :

    //====================================================================================================================
		// declarations for the ESPS/"Talkin" pitch detector
    //====================================================================================================================
#ifdef SC_USE_ESPSPITCH
		#define Fprintf (void)fprintf
		#define BIGSORD 100

		typedef struct f0_params {
		float cand_thresh,	/* only correlation peaks above this are considered */
					lag_weight,	/* degree to which shorter lags are weighted */
					freq_weight,	/* weighting given to F0 trajectory smoothness */
					trans_cost,	/* fixed cost for a voicing-state transition */
					trans_amp,	/* amplitude-change-modulated VUV trans. cost */
					trans_spec,	/* spectral-change-modulated VUV trans. cost */
					voice_bias,	/* fixed bias towards the voiced hypothesis */
					double_cost,	/* cost for octave F0 jumps */
					mean_f0,		/* talker-specific mean F0 (Hz) */
					mean_f0_weight,	/* weight to be given to deviations from mean F0 */
					min_f0,		/* min. F0 to search for (Hz) */
					max_f0,		/* max. F0 to search for (Hz) */
					frame_step,	/* inter-frame-interval (sec) */
					wind_dur;		/* duration of correlation window (sec) */
		int   n_cands,		/* max. # of F0 cands. to consider at each frame */
					conditioning;     /* Specify optional signal pre-conditioning. */
		} F0_params;

		typedef struct cross_rec { /* for storing the crosscorrelation information */
			float	rms;	/* rms energy in the reference window */
			float	maxval;	/* max in the crosscorr. fun. q15 */
			short	maxloc; /* lag # at which max occured	*/
			short	firstlag; /* the first non-zero lag computed */
			float	*correl; /* the normalized corsscor. fun. q15 */
		} Cross;

		typedef struct dp_rec { /* for storing the DP information */
			short	ncands;	/* # of candidate pitch intervals in the frame */
			short	*locs; /* locations of the candidates */
			float	*pvals; /* peak values of the candidates */
			float	*mpvals; /* modified peak values of the candidates */
			short	*prept; /* pointers to best previous cands. */
			float	*dpvals; /* cumulative error for each candidate */
		} Dprec;

		typedef struct windstat_rec {  /* for lpc stat measure in a window */
				float rho[BIGSORD+1];
				float err;
				float rms;
		} Windstat;

		typedef struct sta_rec {  /* for stationarity measure */
			float *stat;
			float *rms;
			float *rms_ratio;
		} Stat;

		typedef struct frame_rec{
			Cross *cp;
			Dprec *dp;
			float rms;
			struct frame_rec *next;
			struct frame_rec *prev;
		} Frame;

		/*
		* headF points to current Frame in the circular buffer, 
		* tailF points to the Frame where tracks start
		* cmpthF points to starting Frame of converged path to backtrack
		*/
		Frame *headF, *tailF, *cmpthF;

		int *pcands;	/* array for backtracking in convergence check */
		int cir_buff_growth_count;

		int size_cir_buffer,	/* # of frames in circular DP buffer */
				size_frame_hist,	/* # of frames required before convergence test */
				size_frame_out,	/* # of frames before forcing output */
				num_active_frames,	/* # of frames from tailF to headF */
				output_buf_size;	/* # of frames allocated to output buffers */

		/* 
		* DP parameters
		*/
		float tcost, tfact_a, tfact_s, frame_int, vbias, fdouble, wdur, ln2, freqwt, lagwt;
		int step, size, nlags, start, stop, ncomp, *locs;
		short maxpeaks;

		int wReuse;  /* number of windows seen before resued */
		Windstat *windstat;

		float *f0p, *vuvp, *rms_speech, *acpkp, *peaks;
		int first_time, pad;

		int debug_level;
		char *ProgName;

		Stat *stat;
		float *mem;

		int get_Nframes(long buffsize, int pad, int step);
		int init_dp_f0(double freq, F0_params *par, long *buffsize, long *sdstep);
		int dp_f0(float *fdata, int buff_size, int sdstep, double freq, F0_params *par, float **f0p_pt, float **vuvp_pt, float **rms_speech_pt, float **acpkp_pt, int *vecsize, int last_time);
		Frame *alloc_frame(int nlags, int ncands);
		int save_windstat(float *rho, int order, float err, float rms);
		int retrieve_windstat(float *rho, int order, float *err, float *rms);
		float get_similarity(int order, int size, float *pdata, float *cdata, float *rmsa, float *rms_ratio, float pre, float stab, int w_type, int init);
		Stat* get_stationarity(float *fdata, double freq, int buff_size, int nframes, int frame_step, int first_time);
		int round(double flnum);
		void get_fast_cands(float *fdata, float *fdsdata, int ind, int step, int size, int dec, int start, int nlags, float *engref, int *maxloc, float *maxval, Cross *cp, float *peaks, int *locs, int *ncand, F0_params *par);
		float *downsample(float *input, int samsin, int state_idx, double freq, int *samsout, int decimate, int first_time, int last_time);
		void get_cand(Cross *cross, float *peak, int *loc, int nlags, int *ncand, float cand_thresh);
		int downsamp(float *in, float *out, int samples, int *outsamps, int state_idx, int decimate, int ncoef, float fc[], int init);
		void do_ffir(register float *buf, register int in_samps, register float *bufo, register int *out_samps, int idx, register int ncoef, float *fc, register int invert, register int skip, register int init);
		int lc_lin_fir(register float fc, int *nf, float *coef);
		void peak(float *y, float *xp, float *yp);
		int get_window(register float *dout, register int n, register int type);
		void rwindow(register float *din, register float *dout, register int n, register float preemp);
		void cwindow(register float *din, register float *dout, register int n, register float preemp);
		void hwindow(register float *din, register float *dout, register int n, register float preemp);
		void hnwindow(register float *din, register float *dout, register int n, register float preemp);
		int window(register float *din, register float *dout, register int n, register float preemp, int type);
		void autoc(register int windowsize, register float *s, register int p, register float *r, register float *e);
		void durbin(register float *r, register float *k, register float *a, register int p, register float *ex);
		void a_to_aca(float *a, float *b, float *c, register int p);
		float itakura(register int p, register float *b, register float *c, register float *r, register float *gain);
		float wind_energy(register float *data, register int size, register int w_type);
		int lpc(int lpc_ord, float lpc_stabl, int wsize, float *data, float *lpca, float *ar, float *lpck, float *normerr, float *rms, float preemp,int type);
		void crossf(float *data, int size, int start, int nlags, float *engref, int *maxloc, float *maxval, float *correl);
		void crossfi(float *data, int size, int start0, int nlags0, int nlags, float *engref, int *maxloc, float *maxval, float *correl, int *locs, int nlocs);

		//by thilo to mimic functions from the ESPS package
		void spsassert(void *pointer, const char *msg); //mimic a macro from the ESPS package that defines an assertion with printed error message
		int getBuffer(float *buffer, int start, int len, float *samples, int overallLen); //get the signal window-wise, return nr. of samples in buffer (without padded zeros if "samples" is over)
		void destructESPSmembers(void); //free all ESPS data members

		//former static variables
		int static_memsize, static_nframes_old;
		float static_b[2048];
		float *static_foutput;
		int static_ncoeff, static_ncoefft;
    float *static_co, *static_mem;
    float static_state[1000];
    int static_fsize, static_resid;
		float *static_din;
		int static_n0;
    int static_wsize;
		float *static_wind;
    int static_wsize_hwindow;
		float *static_wind_hwindow;
    int static_wsize_hnwindow;
		float *static_wind_hnwindow;
		int static_nwind;
		float *static_dwind;
    int static_nwind_lpc;
		float *static_dwind_lpc;
    float *static_dbdata;
		int static_dbsize;
    float *static_dbdata_crossfi;
		int static_dbsize_crossfi;
#endif

		//====================================================================================================================
		// the function to call from ExtractFeatures() to get the result of the ESPS pitch tracker
		//====================================================================================================================
		SV_Data* get_f0(void);

		//====================================================================================================================
		// the function to call from ExtractFeatures() to get the result of the SV_Lib pitch tracker
		//====================================================================================================================
		SV_Data* svlibPitch(void);

    //====================================================================================================================
		// general declarations
    //====================================================================================================================
		bool verbose;
		int method;
		bool takeSqrt;
#ifdef SC_USE_ESPSPITCH
		F0_params get_f0_par;
#endif

	public :

    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_Pitch(int sampleRate, int frameLength, int frameStep, int method, bool takeSqrt, float esps_cand_thresh, float esps_lag_weight, float esps_freq_weight, float esps_trans_cost, float esps_trans_amp, float esps_trans_spec, float esps_voice_bias, float esps_double_cost, float esps_min_f0, float esps_max_f0, int esps_n_cands, float esps_wind_dur, bool verbose = true);
		virtual ~SC_Feature_Pitch();

    //====================================================================================================================
		// override base class method, return pitch sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
};

#endif
