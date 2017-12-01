/**************************************************************************/
/*    Wrapper (i.e. as little as possible is changed from the original    */
/*    code, just enough to put in a class and get rid of tcl/tk stuff)    */
/*    around the ESPS/Talkin formant tracker from libsnack, described     */
/*    there as follows:                                                   */
/*                                                                        */
/*    "A formant tracker based on LPC polynomial roots and dynamic        */
/*    programming. At each frame, the LPC poles are ordered by increasing */
/*    frequency.  All "reasonable" mappings of the poles to F1, F2, ...   */
/*    are performed. The cost of "connecting" each of these mappings with */
/*    each of the mappings in the previous frame is computed.  The lowest */
/*    cost connection is then chosen as the optimum one.  At each frame,  */
/*    each mapping has associated with it a cost based on the formant     */
/*    bandwidths and frequencies.  This "local" cost is finally added to  */
/*    the cost of the best "connection."  At end of utterance (or after a */
/*    reasonable delay like .5sec) the best mappings for the entire       */
/*    utterance may be found by retracing back through best candidate     */
/*    mappings, starting at end of utterance (or current frame)."         */
/*                                                                        */
/*    Original copyright notice:                                          */
/*                                                                        */
/*    "This software has been licensed to the Centre of Speech Technology,*/
/*    KTH by AT&T Corp. and Microsoft Corp. with the terms in the         */
/*    accompanying file BSD.txt, which is a BSD style license.            */
/*    "Copyright (c) 1987-1990  AT&T, Inc.                                */
/*    "Copyright (c) 1986-1990  Entropic Speech, Inc.                     */
/*    "Copyright (c) 1990-1994  Entropic Research Laboratory, Inc.        */
/*                   All rights reserved"                                 */
/*    Written by:  David Talkin                                           */
/*    Revised by: John Shore"                                             */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 06.06.2008																								*/
/**************************************************************************/

#ifndef __SC_Feature_Formant_H__
#define __SC_Feature_Formant_H__

#include "SC_TweakableParameters.h"
#include "SC_Api.h"
#include <SV_Feature.h>

class SC_Feature_Formant : public SV_Feature {

	private :

	protected :

    //====================================================================================================================
		// declarations for the ESPS/Talkin formant tracker
    //====================================================================================================================
#ifdef SC_USE_ESPSFORMANT
		//Structure definitions for the formant tracker.. 
		typedef struct form_latt { /* structure of a DP lattice node for formant tracking */
			short ncand; /* # of candidate mappings for this frame */
			short **cand;      /* pole-to-formant map-candidate array */
			short *prept;	 /* backpointer array for each frame */
			double *cumerr; 	 /* cum. errors associated with each cand. */
		} FORM;
		typedef struct pole_array {   /* structure to hold raw LPC analysis data */
			double rms;    /* rms for current LPC analysis frame */
			double rms2;    /* rms for current F0 analysis frame */
			double f0;     /* fundamental frequency estimate for this frame */
			double pv;		/* probability that frame is voiced */
			double change; /* spec. distance between current and prev. frames */
			short npoles; /* # of complex poles from roots of LPC polynomial */
			double *freq;  /* array of complex pole frequencies (Hz) */
			double *band;  /* array of complex pole bandwidths (Hz) */
		} POLE;

		//soome hardcoded maxima, for formants and lpc-order, respectively
		#define MAXFORMANTS 7
		#define MAXORDER 30

		//debug levels
		#define DEB_PAUSE	8	
		#define DEB_LPC_PARS	4
		#define DEB_PARAMS	2
		#define DEB_ENTRY	1
				
		//Here are the major fudge factors for tweaking the formant tracker.
		//by thilo: all those where formerly "static"
		#define MAXCAN	300 //maximum number of candidate mappings allowed
		double MISSING; //equivalent delta-Hz cost for missing formant
		double NOBAND; //equivalent bandwidth cost of a missing formant;
		double DF_FACT; //cost for proportional frequency changes; with good "stationarity" function: DF_FACT =  80.0
		double DFN_FACT; //cost for proportional dev. from nominal freqs.
		double BAND_FACT; //cost per Hz of bandwidth in the poles
		double F_BIAS; //bias toward selecting low-freq. poles
		double F_MERGE; //cost of mapping f1 and f2 to same frequency
		double	*fre;
		double fnom[7]; //"nominal" freqs.
		double fmins[7]; //frequency bounds
		double fmaxs[7]; //for 1st 5 formants
		int	maxp; //number of poles to consider
		int maxf; //number of formants to find
	  int ncan, domerge;
		short **pc;

		//outer parameters
		int lpc_ord; //lpc prediction order
		int lpc_type; //use bsa's stabilized covariance if != 0
		int w_type; //window type: 0=rectangular; 1=Hamming; 2=cos**4; 3=Hanning
		double ds_freq; //downsample to this samplerate
		double wdur; //frame size for LPC analysis in [s]
		double nom_f1; //nominal F1 frequqncy (to override base settings of hypothesized formant positions)
		double cor_wdur; //window-length for crosscorrelation F0 estimator
		double frame_int; //frame step in [s]
		double preemp; //preemphasis-factor
		int nform; //nr. of formants to track (max. MAXFORMANTS)
		int debug; //debug-level

		//methods for formant tracking
		SV_Data* formantCmd(void); //the "entry point" to the formant tracking stuff
		int canbe(int pnumb, int fnumb); //can this pole be this freq.?
		void candy(int cand, int pnumb, int fnumb); //This does the real work of mapping frequencies to formants.
		void get_fcand(int npole, double *freq, double *band, int nform, short **pcan); //Given a set of pole frequencies and allowable formant frequencies for nform formants, calculate all possible mappings of pole frequencies to formants, including, possibly, mappings with missing formants.
		void set_nominal_freqs(double f1);
		double get_stat_max(register POLE **pole, register int nframes); //find the maximum in the "stationarity" function (stored in rms)
		float** dpform(float **ps, int psLen, int psDim, int psSampleRate, POLE **pole, int nform, double nom_f1);

		//auxiliary functions
		double integerize(register double time, register double freq);
		int eround(register double flnum); //Round the argument to the nearest integer.
		float** lpc_poles(float *sp, int length, int sampleRate, double wdur, double frame_int, int lpc_ord, double preemp, int lpc_type, int w_type, int &resultLength, int &resultSampleRate, POLE** &resultPoles);
		double frand();
		int lpcbsa(int np, double lpc_stabl, int wind, short *data, double *lpc, double *rho, double *nul1, double *nul2, double *energy, double preemp);
		int lc_lin_fir(register double fc, int *nf, double coef[]); //create the coefficients for a symmetric FIR lowpass filter using the window technique with a Hanning window.
		void do_fir(short *buf, int in_samps, short *bufo, int ncoef, short ic[], int invert); //ic contains 1/2 the coefficients of a symmetric FIR filter with unity passband gain.  This filter is convolved with the signal in buf. The output is placed in buf2.  If invert != 0, the filter magnitude response will be inverted.
		int get_abs_maximum(register short *d, register int n);
		int dwnsamp(short *buf, int in_samps, short **buf2, int *out_samps, int insert, int decimate, int ncoef, short ic[], int *smin, int *smax);
		int ratprx(double a, int *k, int *l, int qlim);
		float* Fdownsample(float *s, int sampleRate, double freq2, int start, int end, int &newLength, int &newSampleRate);
		float* highpass(float *s, int length);

		//auxiliary signal processing stuff
		/*int*/long long dlpcwtd(double *s, int *ls, double *p, int *np, double *c, double *phi, double *shi, double *xl, double *w); //pred anal subroutine with ridge reg
		int w_covar(short *xx, int *m, int n, int istrt, double *y, double *alpha, double *r0, double preemp, int w_type); //covariance LPC analysis; originally from Markel and Gray (a translation from the fortran)
		void w_window(register short *din, register double *dout, register int n, register double preemp, int type);
		void rwindow(register short *din, register double *dout, register int n, register double preemp);
		void cwindow(register short *din, register double *dout, register int n, register double preemp);
		void hwindow(register short *din, register double *dout, register int n, register double preemp);
		void hnwindow(register short *din, register double *dout, register int n, register double preemp);
		int formant(int lpc_order, double s_freq, double *lpca, int *n_form, double *freq, double *band, int init); //Find the roots of the LPC denominator polynomial and convert the z-plane zeros to equivalent resonant frequencies and bandwidths. The complex poles are then ordered by frequency.
		int qquad(double a, double b, double c, double *r1r, double *r1i, double *r2r, double *r2i); //find x, where a*x**2 + b*x + c = 0
		int lbpoly(double *a, int order, double *rootr, double *rooti); //Rootr and rooti are assumed to contain starting points for the root search on entry to lbpoly(); return FALSE on error
		void dcwmtrx(double *s, int* ni, int *nl, int *np, double *phi, double *shi, double *ps, double *w); //cov mat for wtd lpc
		int dchlsky(double *a, int *n, double *t, double *det); //performs cholesky decomposition
		void dlwrtrn(double *a, int *n, double *x, double *y); //routine to solve ax=y with cholesky
		/*int*/long long dcovlpc(double *p, double *s, double *a, int *n, double *c); //solve p*a=s using stabilized covariance method
		void dreflpc(double *c, double *a, /*int*/long long *n); //convert ref to lpc
		int lpc(int lpc_ord, double lpc_stabl, int wsize, short *data, double *lpca, double *ar, double *lpck, double *normerr, double *rms, double preemp, int type);
		void autoc(register int windowsize, register double *s, register int p, register double *r, register double *e); //Compute the pp+1 autocorrelation lags of the windowsize samples in s. Return the normalized autocorrelation coefficients in r. The rms is returned in e.
		void durbin (register double *r, register double *k, register double *a, register int p, register double *ex); //Compute the AR and PARCOR coefficients using Durbin's recursion. Note: Durbin returns the coefficients in normal sign format. (i.e. a[0] is assumed to be = +1.)

		//former static local variables in functions, turned into members of the class to ensure proper encapsulation
    double *static_x;
		int static_nold;
		int static_mold;
		double *static_b, *static_beta, *static_grc, *static_cc, static_gam, static_s;
		double *static_rr, *static_ri;
    int static_wsize;
    double *static_wind;
    double *static_dwind;
    int static_nwind;
		int static_i, static_owind, static_wind1;
		/*int*/long long static_mm;
		double static_w[1000];
		double static_Fdownsample_beta, static_Fdownsample_b[256];
		int static_ncoeff, static_ncoefft, static_nbits;
		short static_ic[256];
		short *static_lcf;
		int static_len;
    		
		//these (maybe former static variables) are for the auxiliary signal processing functions; they don't need initialization and are not local only for speed purposes, I guess...
		double *psl,*pp2,*ppl2,*pc2,*pcl,*pph1,*pph2,*pph3,*pphl;
		double *static_pdl1,*static_pdl2,*static_pdl3,*static_pdl4,*static_pdl5,*static_pdl6,*static_pdll;
		double *pa_1,*pa_2,*pa_3,*pa_4,*pa_5,*pal,*pt;
		double *pxl,*pa,*py,*pyl,*pa1,*px;
		double *static_pp,*static_ppl,*static_pa;
		double *static_pa1,*static_pa2,*static_pa3,*static_pa4,*static_pa5,*static_pc;

		//by thilo to mimic functions from the ESPS package
		void destructESPSmembers(void); //free all ESPS data members
#endif

    //====================================================================================================================
		// general declarations
    //====================================================================================================================
		bool verbose;

	public :

    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_Formant(int sampleRate, int frameLength, int frameStep, int esps_lpc_ord, int esps_lpc_type, int esps_w_type, double esps_ds_freq, double esps_wdur, double esps_nom_f1, double esps_cor_wdur, double esps_frame_int, double esps_preemp, int esps_nform, bool verbose = true);
		virtual ~SC_Feature_Formant();

    //====================================================================================================================
		// override base class method, return formant sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
};

#endif
