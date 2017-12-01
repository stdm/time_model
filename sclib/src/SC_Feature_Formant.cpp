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

#include "SC_Feature_Formant.h"
#include "SC_Conversion.h"
#include "SC_FeatureHandler.h"

//====================================================================================================================
// constructor
//====================================================================================================================
SC_Feature_Formant::SC_Feature_Formant(int sampleRate, int frameLength, int frameStep, int esps_lpc_ord, int esps_lpc_type, int esps_w_type, double esps_ds_freq, double esps_wdur, double esps_nom_f1, double esps_cor_wdur, double esps_frame_int, double esps_preemp, int esps_nform, bool verbose) : SV_Feature() {
  //general initializations
	this->Para.WinSz = frameLength;
	this->Para.StpSz = frameStep;
	this->Para.SRate = sampleRate; 
	this->Para.Alpha = esps_preemp;
	this->verbose = verbose;

	//ESPS parameters
#ifdef SC_USE_ESPSFORMANT
	//here come external parameters
	this->lpc_ord = esps_lpc_ord; //12;
	this->lpc_type = esps_lpc_type; //0;
	this->w_type = esps_w_type; //in snack code: 2; in wavesurfer: 1
	this->ds_freq = esps_ds_freq; //10000.0;
	this->wdur = esps_wdur; //.049;
	this->cor_wdur = esps_cor_wdur; //.01;
	this->frame_int = esps_frame_int; //.01;
	this->preemp = esps_preemp; //.7;
	this->nform = esps_nform; //4;
	this->nom_f1 = esps_nom_f1; //-10.0;
		
	//here come internal parameters
	this->debug = (this->verbose == true) ? DEB_ENTRY|DEB_PARAMS|DEB_LPC_PARS|DEB_PAUSE  : 0;
	this->MISSING = 1;
	this->NOBAND = 1000;
	this->DF_FACT =  20.0; //with good "stationarity" function: DF_FACT =  80.0
	this->DFN_FACT = 0.3;
	this->BAND_FACT = .002;
	this->F_BIAS	= 0.000; // 0.0004;
	this->F_MERGE = 2000.0;
	
	//here come internal constants
	this->fnom[0] = 500; this->fnom[1] = 1500; this->fnom[2] = 2500; this->fnom[3] = 3500; this->fnom[4] = 4500; this->fnom[5] = 5500; this->fnom[6] = 6500;
	this->fmins[0] = 50; this->fmins[1] = 400; this->fmins[2] = 1000; this->fmins[3] = 2000; this->fmins[4] = 2000; this->fmins[5] = 3000; this->fmins[6] = 3000;
	this->fmaxs[0] = 1500; this->fmaxs[1] = 3500; this->fmaxs[2] = 4500; this->fmaxs[3] = 5000; this->fmaxs[4] = 6000; this->fmaxs[5] = 6000; this->fmaxs[6] = 8000;

	//here come starting values 
	this->domerge = TRUE;	

	//here come former function-static variables with their respective initializations	
  this->static_x = NULL;
	this->static_nold = 0;
	this->static_mold = 0;
	this->static_b = NULL; 
	this->static_beta = NULL;
	this->static_grc = NULL;
	this->static_cc = NULL;
	MArray_1D(this->static_rr, MAXORDER, double, "SC_Feature_Formant_ this->static_rr");
	MArray_1D(this->static_ri, MAXORDER, double, "SC_Feature_Formant_ this->static_rr");
  this->static_wsize = 0;
  this->static_wind = NULL;
  this->static_dwind = NULL;
  this->static_nwind = 0;
	this->static_owind = 0;
	this->static_Fdownsample_beta = 0.0;
	this->static_ncoeff = 127;
	this->static_ncoefft = 0;
	this->static_nbits = 15;
  this->static_lcf = NULL;
  this->static_len = 0;

	this->fre = NULL;
	this->pc = NULL;
	this->psl = NULL;
	this->pp2 = NULL;
	this->ppl2 = NULL;
	this->pc2 = NULL;
	this->pcl = NULL; 
	this->pph1 = NULL;
	this->pph2 = NULL;
	this->pph3 = NULL;
	this->pphl = NULL;
	this->static_pdl1 = NULL;
	this->static_pdl2 = NULL;
	this->static_pdl3 = NULL;
	this->static_pdl4 = NULL;
	this->static_pdl5 = NULL;
	this->static_pdl6 = NULL;
	this->static_pdll = NULL;
	this->pa_1 = NULL;
	this->pa_2 = NULL;
	this->pa_3 = NULL;
	this->pa_4 = NULL;
	this->pa_5 = NULL;
	this->pal = NULL;
	this->pt = NULL;
	this->pxl = NULL;
	this->pa = NULL;
	this->py = NULL;
	this->pyl = NULL;
	this->pa1 = NULL;
	this->px = NULL;
	this->static_pp = NULL;
	this->static_ppl = NULL;
	this->static_pa = NULL;
	this->static_pa1 = NULL;
	this->static_pa2 = NULL;
	this->static_pa3 = NULL;
	this->static_pa4 = NULL;
	this->static_pa5 = NULL;
	this->static_pc = NULL;
#endif
}

//====================================================================================================================
// default destructor
//====================================================================================================================
SC_Feature_Formant::~SC_Feature_Formant() {
#ifdef SC_USE_ESPSFORMANT
	destructESPSmembers();
#endif
}

//====================================================================================================================
// free all ESPS data members
//====================================================================================================================
void SC_Feature_Formant::destructESPSmembers(void) {
	MFree_1D(this->static_x);
	MFree_1D(this->static_b)
	MFree_1D(this->static_beta)
	MFree_1D(this->static_grc)
	MFree_1D(this->static_cc)
	MFree_1D(this->static_rr);
	MFree_1D(this->static_ri);
	MFree_1D(this->static_wind);
	MFree_1D(this->static_dwind);
	free((void *)(this->static_lcf));

	/*
  free(this->fre);
	free(this->pc);
	free(this->psl);
	free(this->pp2);
	free(this->ppl2);
	free(this->pc2);
	free(this->pcl); 
	free(this->pph1);
	free(this->pph2);
	free(this->pph3);
	free(this->pphl);
	free(this->static_pdl1);
	free(this->static_pdl2);
	free(this->static_pdl3);
	free(this->static_pdl4);
	free(this->static_pdl5);
	free(this->static_pdl6);
	free(this->static_pdll);
	free(this->pa_1);
	free(this->pa_2);
	free(this->pa_3);
	free(this->pa_4);
	free(this->pa_5);
	free(this->pal);
	free(this->pt);
	free(this->pxl);
	free(this->pa);
	free(this->py);
	free(this->pyl);
	free(this->pa1);
	free(this->px);
	free(this->static_pp);
	free(this->static_ppl);
	free(this->static_pa);
	free(this->static_pa1);
	free(this->static_pa2);
	free(this->static_pa3);
	free(this->static_pa4);
	free(this->static_pa5);
	free(this->static_pc);
	*/

	return;
}

//====================================================================================================================
// This is the engine of deriving features
//====================================================================================================================
SV_Data* SC_Feature_Formant::ExtractFeature(void) {
	SV_Data *pTmp = NULL, *pFormants = NULL;
	int originalFrameSize = 0, originalFrameStep = 0;
	SC_FeatureHandler handler(NULL, false); //the method we call doesn't need tweakable parameters
	SC_Conversion converter(this->Para.SRate);
	double ignoreValues[1] = {0.0};

#ifdef SC_USE_ESPSFORMANT
	pTmp = formantCmd(); //uses class members for parameters

	if (pTmp != NULL) {
		originalFrameSize = converter.ms2sample(converter.sample2ms(this->Para.WinSz, this->Para.SRate), sclib::alignmentStart, pTmp->Hdr.sampleRate); //formant tracking maybe involves downsampling, be prepared for that!
		originalFrameStep = converter.ms2sample(converter.sample2ms(this->Para.StpSz, this->Para.SRate), sclib::alignmentStart, pTmp->Hdr.sampleRate); 

		//check if the desired ("original") frame parameters deviate from the one internaly used in the formant tracker, and interpolated to the original ones, if necessary
		if (pTmp->Hdr.frameSize != originalFrameSize || pTmp->Hdr.frameStep != originalFrameStep) {
			pFormants = handler.convertFrameRate(pTmp, this->Len, originalFrameSize, originalFrameStep, ignoreValues);
			MFree_0D(pTmp);
		} else {
			pFormants = pTmp;
		}
	}
#else
	REPORT_ERROR(SVLIB_Fail, "ESPS code not available");
#endif

	return pFormants;
}

#ifdef SC_USE_ESPSFORMANT
//====================================================================================================================
// the function to call from ExtractFeatures() to get the result of the ESPS formant tracker
//====================================================================================================================
SV_Data* SC_Feature_Formant::formantCmd(void)
{
  float *dssnd = NULL, *hpsnd = NULL, **polesnd = NULL;
  float **formantsnd = NULL, *hpsrcsnd, *polesrcsnd;
  int poleLength, poleSRate, newLen = this->Len, newSRate = this->Para.SRate; //by thilo: initialize with original values
  POLE **poles;
  SV_Data *pFormants = NULL;
	SC_Conversion converter;
  
  //Check for errors in specifying parameters
  if(this->nform > (this->lpc_ord-4)/2){
		REPORT_ERROR(SVLIB_BadArg, "Number of formants must be <= (lpc order - 4)/2");
    return NULL;
  }
  if(this->nform > MAXFORMANTS){
    REPORT_ERROR(SVLIB_BadArg, "A maximum of 7 formants are supported at this time");
    return NULL;
  }

	if(this->ds_freq < this->Para.SRate) {
    dssnd = Fdownsample(this->Sig, this->Para.SRate, this->ds_freq, 0, this->Len-1, newLen, newSRate);
  }

  hpsrcsnd = (dssnd ? dssnd : this->Sig);
    
  if (preemp < 1.0) { //be sure DC and rumble are gone!
    hpsnd = highpass(hpsrcsnd, newLen);
  }
  if (hpsrcsnd != this->Sig) { //block by thilo
		MFree_1D(hpsrcsnd);
  }	
  
	polesrcsnd = (hpsnd ? hpsnd : this->Sig);

	if(!(polesnd = lpc_poles(polesrcsnd, newLen, newSRate, this->wdur, this->frame_int, this->lpc_ord, this->preemp, this->lpc_type, this->w_type, poleLength, poleSRate, poles))) {
    REPORT_ERROR(SVLIB_Fail, "Problems in lpc_poles()");
    //TODO: release poles?!
    if (polesrcsnd != this->Sig) { //block by thilo
			MFree_1D(polesrcsnd);
    }
    return NULL;
  }
  if (polesrcsnd != this->Sig) { //block by thilo
		MFree_1D(polesrcsnd);
  }	  

  //LPC poles are now available for the formant estimator.
  if (!(formantsnd = dpform(polesnd, poleLength, this->lpc_ord, poleSRate, poles, this->nform, this->nom_f1))) {
    REPORT_ERROR(SVLIB_Fail, "Problems in dpform()");
    //TODO: release poles!?
    MFree_2D(polesnd);
    return NULL;
  }
  MFree_2D(polesnd); //Snack_DeleteSound(polesnd);

	if (poleLength > 0) {
		pFormants = new SV_Data();
		pFormants->Row = poleLength;
		pFormants->Col = this->nform * 2; //first come the this->nform formants (F1, F2, ..., Fn), then their respective bandwidth' (bw-F1, bw-F2, ..., bw-Fn)
		pFormants->Mat = formantsnd;
		pFormants->Hdr.ID = sclib::featureFormant;
		pFormants->Hdr.sampleRate = sclib::min(this->Para.SRate, this->ds_freq);
		converter.setAudioSampleRate(pFormants->Hdr.sampleRate);
		pFormants->Hdr.frameSize = converter.ms2sample(sclib::round(this->wdur*1000.0));
		pFormants->Hdr.frameStep = converter.ms2sample(sclib::round(this->frame_int*1000.0));
	}
	
  return pFormants;
}

//====================================================================================================================
// Below is the code taken and adapted from the ESPS project
//====================================================================================================================

//can this pole be this freq.?
int SC_Feature_Formant::canbe(int pnumb, int fnumb) 
{
return((this->fre[pnumb] >= this->fmins[fnumb])&&(this->fre[pnumb] <= this->fmaxs[fnumb]));
}

//This does the real work of mapping frequencies to formants.
//cand: candidate number being considered
//pnumb: pole number under consideration
//fnumb: formant number under consideration
void SC_Feature_Formant::candy(int cand, int pnumb, int fnumb)
{
  int i,j;

  if(fnumb < this->maxf) this->pc[cand][fnumb] = -1;
  if((pnumb < this->maxp)&&(fnumb < this->maxf)){
    /*   printf("\ncan:%3d  pnumb:%3d  fnumb:%3d",cand,pnumb,fnumb); */
    if(canbe(pnumb,fnumb)){
      this->pc[cand][fnumb] = pnumb;
      if(this->domerge&&(fnumb==0)&&(canbe(pnumb,fnumb+1))){ /* allow for f1,f2 merger */
	this->ncan++;
	this->pc[this->ncan][0] = this->pc[cand][0];
	candy(this->ncan,pnumb,fnumb+1); /* same pole, next formant */
      }
      candy(cand,pnumb+1,fnumb+1); /* next formant; next pole */
      if(((pnumb+1) < this->maxp) && canbe(pnumb+1,fnumb)){
	/* try other frequencies for this formant */
	this->ncan++;			/* add one to the candidate index/tally */
	/*		printf("\n%4d  %4d  %4d",ncan,pnumb+1,fnumb); */
	for(i=0; i<fnumb; i++)	/* clone the lower formants */
	  this->pc[this->ncan][i] = this->pc[cand][i];
	candy(this->ncan,pnumb+1,fnumb);
      }
    } else {
      candy(cand,pnumb+1,fnumb);
    }
  }
  /* If all pole frequencies have been examined without finding one which
     will map onto the current formant, go on to the next formant leaving the
     current formant null. */
  if((pnumb >= this->maxp) && (fnumb < this->maxf-1) && (this->pc[cand][fnumb] < 0)){
    if(fnumb){
      j=fnumb-1;
      while((j>0) && this->pc[cand][j] < 0) j--;
      i = ((j=this->pc[cand][j]) >= 0)? j : 0;
    } else i = 0;
    candy(cand,i,fnumb+1);
  }
}

//Given a set of pole frequencies and allowable formant frequencies
//for nform formants, calculate all possible mappings of pole frequencies
//to formants, including, possibly, mappings with missing formants.
//freq, band: poles ordered by increasing FREQUENCY
void SC_Feature_Formant::get_fcand(int npole, double *freq, double *band, int nform, short **pcan)
{	
  this->ncan = 0;
  this->pc = pcan;
  this->fre = freq;
  this->maxp = npole;
  this->maxf = nform;
  candy(this->ncan, 0, 0);
  this->ncan++;	/* (converts ncan as an index to ncan as a candidate count) */
}

void SC_Feature_Formant::set_nominal_freqs(double f1)
{
  int i;
  for(i=0; i < MAXFORMANTS; i++) {
    this->fnom[i] = ((i * 2) + 1) * f1;
    this->fmins[i] = this->fnom[i] - ((i+1) * f1) + 50.0;
    this->fmaxs[i] = this->fnom[i] + (i * f1) + 1000.0;
  }
}

/*      ----------------------------------------------------------      */

//find the maximum in the "stationarity" function (stored in rms)
double SC_Feature_Formant::get_stat_max(register POLE **pole, register int nframes)
{
  register int i;
  register double amax, t;

  for(i=1, amax = (*pole++)->rms; i++ < nframes; )
    if((t = (*pole++)->rms) > amax) amax = t;

  return(amax);
}

//by thilo: changed "Sound" datatype to "float**"
float** SC_Feature_Formant::dpform(float **ps, int psLen, int psDim, int psSampleRate, POLE **pole, int nform, double nom_f1)
{
  double pferr, conerr, minerr, dffact, ftemp, berr, ferr, bfact, ffact,
         rmsmax, fbias, **fr, **ba, rmsdffact, merger=0.0, merge_cost,
         FBIAS;
  register int	i, j, k, l, ic, ip, mincan=0;
  short	**pcan;
  FORM	**fl;
  //POLE	**pole; /* raw LPC pole data structure array */
  float **fbs; //Sound *fbs;
  int dmaxc,dminc,dcountc,dcountf;
  
  if(ps) {
    if(nom_f1 > 0.0)
      set_nominal_freqs(nom_f1);
    //pole = (POLE**)ps->extHead;
    rmsmax = get_stat_max(pole, psLen);
    FBIAS = F_BIAS /(.01 * psSampleRate);
    /* Setup working values of the cost weights. */
    dffact = (this->DF_FACT * .01) * psSampleRate; /* keep dffact scaled to frame rate */
    bfact = this->BAND_FACT /(.01 * psSampleRate);
    ffact = this->DFN_FACT /(.01 * psSampleRate);
    merge_cost = this->F_MERGE;
    if(merge_cost > 1000.0) this->domerge = FALSE;

    /* Allocate space for the formant and bandwidth arrays to be passed back. */
    if(this->debug & DEB_ENTRY){
      printf("Allocating formant and bandwidth arrays in dpform()\n");
    }
    fr = (double**)malloc(sizeof(double*) * nform * 2);
    ba = fr + nform;
    for(i=0;i < nform*2; i++){
      fr[i] = (double*)malloc(sizeof(double) * psLen);
    }
    /*    cp = new_ext(ps->name,"fb");*/
    /*    if((fbs=new_signal(cp,SIG_UNKNOWN,dup_header(ps->header),fr,ps->length,		       ps->samprate, nform * 2))) {*/
    if (1) {
      /* Allocate space for the raw candidate array. */
      if(this->debug & DEB_ENTRY){
	printf("Allocating raw candidate array in dpform()\n");
      }
      pcan = (short**)malloc(sizeof(short*) * MAXCAN);
      for(i=0;i<MAXCAN;i++) pcan[i] = (short*)malloc(sizeof(short) * nform);

      /* Allocate space for the dp lattice */
      if(this->debug & DEB_ENTRY){
	printf("Allocating DP lattice structure in dpform()\n");
      }
      fl = (FORM**)malloc(sizeof(FORM*) * psLen);
      for(i=0;i < psLen; i++)
	fl[i] = (FORM*)malloc(sizeof(FORM));

      /*******************************************************************/
      /* main formant tracking loop */
      /*******************************************************************/
      if(this->debug & DEB_ENTRY){
	printf("Entering main computation loop in dpform()\n");
      }
      for(i=0; i < psLen; i++){	/* for all analysis frames... */

	ncan = 0;		/* initialize candidate mapping count to 0 */

	/* moderate the cost of frequency jumps by the relative amplitude */
	rmsdffact = pole[i]->rms;
	rmsdffact = rmsdffact/rmsmax;
	rmsdffact = rmsdffact * dffact;

	/* Get all likely mappings of the poles onto formants for this frame. */
	if(pole[i]->npoles){	/* if there ARE pole frequencies available... */
	  get_fcand(pole[i]->npoles,pole[i]->freq,pole[i]->band,nform,pcan);

	  /* Allocate space for this frame's candidates in the dp lattice. */
	  fl[i]->prept =  (short*)malloc(sizeof(short) * this->ncan);
	  fl[i]->cumerr = (double*)malloc(sizeof(double) * this->ncan);
	  fl[i]->cand =   (short**)malloc(sizeof(short*) * this->ncan);
	  for(j=0;j<this->ncan;j++){	/* allocate cand. slots and install candidates */
	    fl[i]->cand[j] = (short*)malloc(sizeof(short) * nform);
	    for(k=0; k<nform; k++)
	      fl[i]->cand[j][k] = pcan[j][k];
	  }
	}
	fl[i]->ncand = this->ncan;
	/* compute the distance between the current and previous mappings */
	for(j=0;j<this->ncan;j++){	/* for each CURRENT mapping... */
	  if( i ){		/* past the first frame? */
	    minerr = 0;
	    if(fl[i-1]->ncand) minerr = 2.0e30;
	    mincan = -1;
	    for(k=0; k < fl[i-1]->ncand; k++){ /* for each PREVIOUS map... */
	      for(pferr=0.0, l=0; l<nform; l++){
		ic = fl[i]->cand[j][l];
		ip = fl[i-1]->cand[k][l];
		if((ic >= 0)	&& (ip >= 0)){
		  ftemp = 2.0 * fabs(pole[i]->freq[ic] - pole[i-1]->freq[ip])/
		           (pole[i]->freq[ic] + pole[i-1]->freq[ip]);
    /*		  ftemp = pole[i]->freq[ic] - pole[i-1]->freq[ip];
		  if(ftemp >= 0.0)
		    ftemp = ftemp/pole[i-1]->freq[ip];
		  else
		    ftemp = ftemp/pole[i]->freq[ic]; */
		  /* cost prop. to SQUARE of deviation to discourage large jumps */
		  pferr += ftemp * ftemp;
		}
		else pferr += this->MISSING;
	      }
	      /* scale delta-frequency cost and add in prev. cum. cost */
	      conerr = (rmsdffact * pferr) + fl[i-1]->cumerr[k]; 
	      if(conerr < minerr){
		minerr = conerr;
		mincan = k;
	      }
	    }			/* end for each PREVIOUS mapping... */
	  }	else {		/* (i.e. if this is the first frame... ) */
	    minerr = 0;
	  }

	  fl[i]->prept[j] = mincan; /* point to best previous mapping */
	  /* (Note that mincan=-1 if there were no candidates in prev. fr.) */
	  /* Compute the local costs for this current mapping. */
	  for(k=0, berr=0, ferr=0, fbias=0; k<nform; k++){
	    ic = fl[i]->cand[j][k];
	    if(ic >= 0){
	      if( !k ){		/* F1 candidate? */
		ftemp = pole[i]->freq[ic];
		merger = (this->domerge &&
			  (ftemp == pole[i]->freq[fl[i]->cand[j][1]]))?
			  merge_cost: 0.0;
	      }
	      berr += pole[i]->band[ic];
	      ferr += (fabs(pole[i]->freq[ic]-fnom[k])/fnom[k]);
	      fbias += pole[i]->freq[ic];
	    } else {		/* if there was no freq. for this formant */
	      fbias += this->fnom[k];
	      berr += this->NOBAND;
	      ferr += this->MISSING;
	    }
	  }

	  /* Compute the total cost of this mapping and best previous. */
	  fl[i]->cumerr[j] = (FBIAS * fbias) + (bfact * berr) + merger +
	                     (ffact * ferr) + minerr;
	}			/* end for each CURRENT mapping... */

	if(this->debug & DEB_LPC_PARS){
	  printf("\nFrame %4d  # candidates:%3d stat:%f prms:%f",i,ncan,rmsdffact,pole[i]->rms);
	  for (j=0; j<this->ncan; j++){
	    printf("\n	");
	    for(k=0; k<nform; k++)
	      if(pcan[j][k] >= 0)
		printf("%6.0f ",pole[i]->freq[fl[i]->cand[j][k]]);
	      else
		printf("  NA   ");
	    printf("  cum:%7.2f pp:%d",fl[i]->cumerr[j], fl[i]->prept[j]);
	  }
	}
      }				/* end for all analysis frames... */	
      /**************************************************************************/

      /* Pick the candidate in the final frame with the lowest cost. */
      /* Starting with that min.-cost cand., work back thru the lattice. */
      if(this->debug & DEB_ENTRY){
	printf("Entering backtrack loop in dpform()\n");
      }
      dmaxc = 0;
      dminc = 100;
      dcountc = dcountf = 0;
      for(mincan = -1, i=psLen - 1; i>=0; i--){
	if(this->debug & DEB_LPC_PARS){
	  printf("\nFrame:%4d mincan:%2d ncand:%2d ",i,mincan,fl[i]->ncand);
	}
	if(mincan < 0)		/* need to find best starting candidate? */
	  if(fl[i]->ncand){	/* have candidates at this frame? */
	    minerr = fl[i]->cumerr[0];
	    mincan = 0;
	    for(j=1; j<fl[i]->ncand; j++)
	      if( fl[i]->cumerr[j] < minerr ){
		minerr = fl[i]->cumerr[j];
		mincan = j;
	      }
	  }
	if(mincan >= 0){	/* if there is a "best" candidate at this frame */
	  if((j = fl[i]->ncand) > dmaxc) dmaxc = j;
	  else
	    if( j < dminc) dminc = j;
	  dcountc += j;
	  dcountf++;
	  for(j=0; j<nform; j++){
	    k = fl[i]->cand[mincan][j];
	    if(k >= 0){
	      fr[j][i] = pole[i]->freq[k];
	      if(this->debug & DEB_LPC_PARS){
		printf("%6.0f",fr[j][i]);
	      }
	      ba[j][i] = pole[i]->band[k];
	    } else {		/* IF FORMANT IS MISSING... */
	      if(i < psLen - 1){
		fr[j][i] = fr[j][i+1]; /* replicate backwards */
		ba[j][i] = ba[j][i+1];
	      } else {
		fr[j][i] = this->fnom[j]; /* or insert neutral values */
		ba[j][i] = NOBAND;
	      }
	      if(this->debug & DEB_LPC_PARS){
		printf("%6.0f",fr[j][i]);
	      }
	    }
	  }
	  mincan = fl[i]->prept[mincan];
	} else {		/* if no candidates, fake with "nominal" frequencies. */
	  for(j=0; j < nform; j++){
	    fr[j][i] = this->fnom[j];
	    ba[j][i] = this->NOBAND;
	    if(debug & DEB_LPC_PARS){
	      printf("%6.0f",fr[j][i]);
	    }
	  } 
	}			/* note that mincan will remain =-1 if no candidates */
      }				/* end unpacking formant tracks from the dp lattice */
      /* Deallocate all the DP lattice work space. */
      /*if(debug & DEB_ENTRY){
	printf("%s complete; max. cands:%d  min. cands.:%d average cands.:%f\n",
	     fbs->name,dmaxc,dminc,((double)dcountc)/dcountf);
	printf("Entering memory deallocation in dpform()\n");
      }*/
      for(i=psLen - 1; i>=0; i--){
	if(fl[i]->ncand){
	  if(fl[i]->cand) {
	    for(j=0; j<fl[i]->ncand; j++) free((void *)fl[i]->cand[j]);
	    free((void *)fl[i]->cand);
	    free((void *)fl[i]->cumerr);
	    free((void *)fl[i]->prept);
	  }
	}
      }
      for(i=0; i<psLen; i++)	free((void *)fl[i]);
      free((void *)fl);
      fl = 0;
      
      for(i=0; i<psLen; i++) {
	free((void *)pole[i]->freq);
	free((void *)pole[i]->band);
	free((void *)pole[i]);
      }
      free((void *)pole);

      /* Deallocate space for the raw candidate aray. */
      for(i=0;i<MAXCAN;i++) free((void *)pcan[i]);
      free((void *)pcan);

      //fbs = Snack_NewSound(ps->samprate, SNACK_FLOAT, nform * 2);
      //Snack_ResizeSoundStorage(fbs, ps->length);
      MArray_2D(fbs, psLen, nform*2, float, "SC_Feature_Formant.dpform: fbs");
      for (i = 0; i < psLen; i++) {
	for (j = 0; j < nform * 2; j++) {
	  fbs[i][j] = (float)fr[j][i]; //Snack_SetSample(fbs, j, i, (float)fr[j][i]);
	}
      }
      //fbs->length = psLen;

      for(i = 0; i < nform*2; i++) free((void *)fr[i]);
      free((void *)fr);

      return(fbs);
    } else
      printf("Can't create a new Signal in dpform()\n");
  } else
    printf("Bad data pointers passed into dpform()\n");
  return(NULL);
}

/* lpc_poles.c */

/* computation and I/O routines for dealing with LPC poles */

/*************************************************************************/
double SC_Feature_Formant::integerize(register double time, register double freq)
{
  register int i;

  i = (int) (.5 + (freq * time));
  return(((double)i)/freq);
}

/*	Round the argument to the nearest integer.			*/
int eround(register double flnum)
{
	return((flnum >= 0.0) ? (int)(flnum + 0.5) : (int)(flnum - 0.5));
}

/*************************************************************************/
//by thilo: changed "Sound*" data type to "float*(*)"
float** SC_Feature_Formant::lpc_poles(float *sp, int length, int sampleRate, double wdur, double frame_int, int lpc_ord, double preemp, int lpc_type, int w_type, int &resultLength, int &resultSampleRate, POLE** &resultPoles)
{
  int i, j, size, step, nform, init, nfrm;
  POLE **pole;
  double lpc_stabl = 70.0, energy, lpca[MAXORDER], normerr,
         *bap=NULL, *frp=NULL, *rhp=NULL;
  short *datap, *dporg;
  float **lp; //Sound *lp;

  if(lpc_type == 1) { /* force "standard" stabilized covariance (ala bsa) */
    wdur = 0.025;
    preemp = exp(-62.831853 * 90. / sampleRate); /* exp(-1800*pi*T) */
  }
  if((lpc_ord > MAXORDER) || (lpc_ord < 2)/* || (! ((short**)sp->data)[0])*/)
    return(NULL);
  /*  np = (char*)new_ext(sp->name,"pole");*/
  wdur = integerize(wdur,(double)sampleRate);
  frame_int = integerize(frame_int,(double)sampleRate);
  nfrm= 1 + (int) (((((double)length)/sampleRate) - wdur)/(frame_int));
  if(nfrm >= 1/*lp->buff_size >= 1*/) {
    size = (int) (.5 + (wdur * sampleRate));
    step = (int) (.5 + (frame_int * sampleRate));
    pole = (POLE**)malloc(nfrm/*lp->buff_size*/ * sizeof(POLE*));
    datap = dporg = (short *) malloc(sizeof(short) * length);
    for (i = 0; i < length; i++) {
      datap[i] = (short) sp[i]; //Snack_GetSample(sp, 0, i);
    }
    for(j=0, init=TRUE/*, datap=((short**)sp->data)[0]*/; j < nfrm/*lp->buff_size*/;j++, datap += step){
      pole[j] = (POLE*)malloc(sizeof(POLE));
      pole[j]->freq = frp = (double*)malloc(sizeof(double)*lpc_ord);
      pole[j]->band = bap = (double*)malloc(sizeof(double)*lpc_ord);

      switch(lpc_type) {
      case 0:
	if(! lpc(lpc_ord,lpc_stabl,size,datap,lpca,rhp,NULL,&normerr,
		 &energy, preemp, w_type)){
	  printf("Problems with lpc in lpc_poles()");
	  break;
	}
	break;
      case 1:
	if(! lpcbsa(lpc_ord,lpc_stabl,size,datap,lpca,rhp,NULL,&normerr,
		    &energy, preemp)){
          printf("Problems with lpc in lpc_poles()");
	  break;
	}
	break;
      case 2:
	{
	  int Ord = lpc_ord;
	  double alpha, r0;

	  w_covar(datap, &Ord, size, 0, lpca, &alpha, &r0, preemp, 0);
	  if((Ord != lpc_ord) || (alpha <= 0.0))
	    printf("Problems with covar(); alpha:%f  Ord:%d\n",alpha,Ord);
	  energy = sqrt(r0/(size-Ord));
	}
	break;
      }
      pole[j]->change = 0.0;
       /* don't waste time on low energy frames */
       if((pole[j]->rms = energy) > 1.0){
	 formant(lpc_ord,(double)sampleRate, lpca, &nform, frp, bap, init);
	 pole[j]->npoles = nform;
	 init=FALSE;		/* use old poles to start next search */
       } else {			/* write out no pole frequencies */
	 pole[j]->npoles = 0;
	 init = TRUE;		/* restart root search in a neutral zone */
       }
/*     if(debug & 4) {
	 printf("\nfr:%4d np:%4d rms:%7.0f  ",j,pole[j]->npoles,pole[j]->rms);
	 for(k=0; k<pole[j]->npoles; k++)
	   printf(" %7.1f",pole[j]->freq[k]);
	 printf("\n                   ");
	 for(k=0; k<pole[j]->npoles; k++)
	   printf(" %7.1f",pole[j]->band[k]);
	 printf("\n");
	 }*/
     } /* end LPC pole computation for all lp->buff_size frames */
    /*     lp->data = (caddr_t)pole;*/
    free((void *)dporg);
    //lp = Snack_NewSound((int)(1.0/frame_int), LIN16, lpc_ord);
    //Snack_ResizeSoundStorage(lp, nfrm);
    resultSampleRate = (int)(1.0/frame_int);
    MArray_2D(lp, nfrm, lpc_ord, float, "SC_Feature_Formant.lpc_poles: lp");
    for (i = 0; i < nfrm; i++) {
      for (j = 0; j < lpc_ord; j++) {
        lp[i][j] = (float)pole[i]->freq[j]; //Snack_SetSample(lp, j, i, (float)pole[i]->freq[j]);
      }
    }
    //lp->length = nfrm;
    resultLength = nfrm;
    //lp->extHead = (char *)pole;
    resultPoles = pole;
    return(lp);
  } else {
    printf("Bad buffer size in lpc_poles()\n");
  }
  return(NULL);
}

/**********************************************************************/
double SC_Feature_Formant::frand()
{
  return (((double)rand())/(double)RAND_MAX);
}
    
/**********************************************************************/
/* a quick and dirty interface to bsa's stabilized covariance LPC */
#define NPM	30	/* max lpc order		*/

int SC_Feature_Formant::lpcbsa(int np, double lpc_stabl, int wind, short *data, double *lpc, double *rho, double *nul1, double *nul2, double *energy, double preemp)
{
  //static int i, mm, owind=0, wind1;
  //static double w[1000];
  double rc[NPM],phi[NPM*NPM],shi[NPM],sig[1000];
  double xl = .09, fham, amax;
  register double *psp1, *psp3, *pspl;

  if(this->static_owind != wind) {		/* need to compute a new window? */
    fham = 6.28318506 / wind;
    for(psp1=this->static_w,this->static_i=0;this->static_i<wind;this->static_i++,psp1++)
      *psp1 = .54 - .46 * cos(this->static_i * fham);
    this->static_owind = wind;
  }
  wind += np + 1;
  this->static_wind1 = wind-1;

  for(psp3=sig,pspl=sig+wind; psp3 < pspl; )
    *psp3++ = (double)(*data++) + .016 * frand() - .008;
  for(psp3=sig+1,pspl=sig+wind;psp3<pspl;psp3++)
    *(psp3-1) = *psp3 - preemp * *(psp3-1);
  for(amax = 0.,psp3=sig+np,pspl=sig+this->static_wind1;psp3<pspl;psp3++)
    amax += *psp3 * *psp3;
  *energy = sqrt(amax / (double)this->static_owind);
  amax = 1.0/(*energy);
	
  for(psp3=sig,pspl=sig+this->static_wind1;psp3<pspl;psp3++)
    *psp3 *= amax;
  if((this->static_mm=dlpcwtd(sig,&this->static_wind1,lpc,&np,rc,phi,shi,&xl,this->static_w))!=np) {
    printf("LPCWTD error mm<np %d %d\n",this->static_mm,np);
    return(FALSE);
  }
  return(TRUE);
}

/*	Copyright (c) 1987, 1988, 1989 AT&T	*/
/*	  All Rights Reserved	*/

/*	THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF AT&T	*/
/*	The copyright notice above does not evidence any	*/
/*	actual or intended publication of such source code.	*/

/* downsample.c */
/* a quick and dirty downsampler */

#ifndef TRUE
# define TRUE 1
# define FALSE 0
#endif

//ifndef by thilo
#ifndef PI
# define PI 3.1415927
#endif

/*      ----------------------------------------------------------      */

//create the coefficients for a symmetric FIR lowpass filter using the
//window technique with a Hanning window.
int SC_Feature_Formant::lc_lin_fir(register double fc, int *nf, double coef[])
{
    register int	i, n;
    register double	twopi, fn, c;

    if(((*nf % 2) != 1) || (*nf > 127)) {
	if(*nf <= 126) *nf = *nf + 1;
	else *nf = 127;
    }
    n = (*nf + 1)/2;

    /*  compute part of the ideal impulse response */
    twopi = PI * 2.0;
    coef[0] = 2.0 * fc;
    c = PI;
    fn = twopi * fc;
    for(i=1;i < n; i++) coef[i] = sin(i * fn)/(c * i);

    /* Now apply a Hanning window to the (infinite) impulse response. */
    fn = twopi/((double)(*nf - 1));
    for(i=0;i<n;i++) 
	coef[i] *= (.5 + (.5 * cos(fn * ((double)i))));
    
    return(TRUE);
}

/*      ----------------------------------------------------------      */

//ic contains 1/2 the coefficients of a symmetric FIR filter with unity
//passband gain.  This filter is convolved with the signal in buf.
//The output is placed in buf2.  If invert != 0, the filter magnitude
//response will be inverted.
void SC_Feature_Formant::do_fir(short *buf, int in_samps, short *bufo, int ncoef, short ic[], int invert)
{
    register short  *buft, *bufp, *bufp2, stem;
    short co[256], mem[256];
    register int i, j, k, l, m, sum, integral;
    
    for(i=ncoef-1, bufp=ic+ncoef-1, bufp2=co, buft = co+((ncoef-1)*2),
	integral = 0; i-- > 0; )
      if(!invert) *buft-- = *bufp2++ = *bufp--;
      else {
	integral += (stem = *bufp--);
	*buft-- = *bufp2++ = -stem;
      }
    if(!invert)  *buft-- = *bufp2++ = *bufp--; /* point of symmetry */
    else {
      integral *= 2;
      integral += *bufp;
      *buft-- = integral - *bufp;
    }
/*         for(i=(ncoef*2)-2; i >= 0; i--) printf("\n%4d%7d",i,co[i]);  */
    for(i=ncoef-1, buft=mem; i-- > 0; ) *buft++ = 0;
    for(i=ncoef; i-- > 0; ) *buft++ = *buf++;
    l = 16384;
    m = 15;
    k = (ncoef << 1) -1;
    for(i=in_samps-ncoef; i-- > 0; ) {
      for(j=k, buft=mem, bufp=co, bufp2=mem+1, sum = 0; j-- > 0;
	  *buft++ = *bufp2++ )
	sum += (((*bufp++ * *buft) + l) >> m);

      *--buft = *buf++;		/* new data to memory */
      *bufo++ = sum; 
    }
    for(i=ncoef; i-- > 0; ) {	/* pad data end with zeros */
      for(j=k, buft=mem, bufp=co, bufp2=mem+1, sum = 0; j-- > 0;
	  *buft++ = *bufp2++ )
	sum += (((*bufp++ * *buft) + l) >> m);
      *--buft = 0;
      *bufo++ = sum; 
    }
}

/* ******************************************************************** */

int SC_Feature_Formant::get_abs_maximum(register short *d, register int n)
{
  register int i;
  register short amax, t;

  if((t = *d++) >= 0) amax = t;
  else amax = -t;
  
  for(i = n-1; i-- > 0; ) {
    if((t = *d++) > amax) amax = t;
    else {
      if(-t > amax) amax = -t;
    }
  }
  return((int)amax);
}

/* ******************************************************************** */

int SC_Feature_Formant::dwnsamp(short *buf, int in_samps, short **buf2, int *out_samps, int insert, int decimate, int ncoef, short ic[], int *smin, int *smax)
{
  register short  *bufp, *bufp2;
  short	*buft;
  register int i, j, k, l, m;
  int imax, imin;

  if(!(*buf2 = buft = (short*)malloc(sizeof(short)*insert*in_samps))) {
    perror("malloc() in dwnsamp()");
    return(FALSE);
  } 

  k = imax = get_abs_maximum(buf,in_samps);
  if (k == 0) k = 1;
  if(insert > 1) k = (32767 * 32767)/k;	/*  prepare to scale data */
  else k = (16384 * 32767)/k;
  l = 16384;
  m = 15;
    

  /* Insert zero samples to boost the sampling frequency and scale the
     signal to maintain maximum precision. */
  for(i=0, bufp=buft, bufp2=buf; i < in_samps; i++) {
    *bufp++ = ((k * (*bufp2++)) + l) >> m ; 
    for(j=1; j < insert; j++) *bufp++ = 0;
  }
    
  do_fir(buft,in_samps*insert,buft,ncoef,ic,0);
    
  /*	Finally, decimate and return the downsampled signal. */
  *out_samps = j = (in_samps * insert)/decimate;
  k = decimate;
  for(i=0, bufp=buft, imax = imin = *bufp; i < j; bufp += k,i++) {
    *buft++ = *bufp;
    if(imax < *bufp) imax = *bufp;
    else
      if(imin > *bufp) imin = *bufp;
  }
  *smin = imin;
  *smax = imax;
  *buf2 = (short*)realloc((void *) *buf2, sizeof(short) * (*out_samps));
  return(TRUE);
}

/*      ----------------------------------------------------------      */

int SC_Feature_Formant::ratprx(double a, int *k, int *l, int qlim)
{
    double aa, af, q, em, qq = 0, pp = 0, ps, e;
    int	ai, ip, i;
    
    aa = fabs(a);
    ai = (int) aa;
/*    af = fmod(aa,1.0); */
    i = (int) aa;
    af = aa - i;
    q = 0;
    em = 1.0;
    while(++q <= qlim) {
	ps = q * af;
	ip = (int) (ps + 0.5);
	e = fabs((ps - (double)ip)/q);
	if(e < em) {
	    em = e;
	    pp = ip;
	    qq = q;
	}
    };
    *k = (int) ((ai * qq) + pp);
    *k = (a > 0)? *k : -(*k);
    *l = (int) qq;
    return(TRUE);    
}

/* ----------------------------------------------------------------------- */

//by thilo: changed data type "Sound*" to float*", see comments for what was previously written
float* SC_Feature_Formant::Fdownsample(float *s, int sampleRate, double freq2, int start, int end, int &newLength, int &newSampleRate) //Sound* Fdownsample(Sound *s, double freq2, ...
{
  short	*bufin, **bufout;
  //static double	beta = 0.0, b[256];
  double	ratio_t, maxi, ratio, beta_new, freq1;
  //static int	ncoeff = 127, ncoefft = 0, nbits = 15;
  //static short	ic[256];
  int	insert, decimate, out_samps, smin, smax;
  float *so; //Sound *so;s

  register int i, j;

  freq1 = sampleRate; //s->samprate;
  
  if((bufout = (short**)malloc(sizeof(short*)))) {
    bufin = (short *) malloc(sizeof(short) * (end - start + 1));
    for (i = start; i <= end; i++) {
      bufin[i-start] = (short) s[i]; //(short) Snack_GetSample(s, 0, i);
    }

    ratio = freq2/freq1;
    ratprx(ratio,&insert,&decimate,10);
    ratio_t = ((double)insert)/((double)decimate);

    if(ratio_t > .99) return(s);
  
    freq2 = ratio_t * freq1;
    beta_new = (.5 * freq2)/(insert * freq1);

    if(this->static_Fdownsample_beta != beta_new){
      this->static_Fdownsample_beta = beta_new;
      if( !lc_lin_fir(this->static_Fdownsample_beta,&this->static_ncoeff,this->static_Fdownsample_b)) {
	printf("\nProblems computing interpolation filter\n");
	return(FALSE);
      }
      maxi = (1 << this->static_nbits) - 1;
      j = (this->static_ncoeff/2) + 1;
      for(this->static_ncoefft = 0, i=0; i < j; i++){
	this->static_ic[i] = (int) (0.5 + (maxi * this->static_Fdownsample_b[i]));
	if(this->static_ic[i]) this->static_ncoefft = i+1;
      }
    }				/*  endif new coefficients need to be computed */

    if(dwnsamp(bufin,end-start+1,bufout,&out_samps,insert,decimate,this->static_ncoefft,this->static_ic,
	       &smin,&smax)){
      /*      so->buff_size = so->file_size = out_samps;*/
      //so = Snack_NewSound(0, LIN16, s->nchannels);
      //Snack_ResizeSoundStorage(so, out_samps);
      MArray_1D(so, out_samps, float, "SC_Feature_Formant.Fdownsample: so");
      for (i = 0; i < out_samps; i++) {
				so[i] = (float)(*bufout)[i]; //Snack_SetSample(so, 0, i, (float)(*bufout)[i]);
      }
      newLength = out_samps; //so->length = out_samps;
      newSampleRate = (int)freq2; //so->samprate = (int)freq2;
      free((void *)*bufout);
      free((void *)bufout);
      free((void *)bufin);
      return(so);
    } else
      printf("Problems in dwnsamp() in downsample()\n");
  } else
       printf("Can't create a new Signal in downsample()\n");
  
  return(NULL);
}

/*      ----------------------------------------------------------      */

//by thilo changed data type from "Sound*" to "float*"
float* SC_Feature_Formant::highpass(float *s, int length) //Sound* highpass(Sound *s)
{

  short *datain, *dataout;
  //static short *lcf;
  //static int len = 0;
  double scale, fn;
  register int i;
  float *so; //Sound *so;

  /*  Header *h, *dup_header();*/
  
#define LCSIZ 101
  /* This assumes the sampling frequency is 10kHz and that the FIR
     is a Hanning function of (LCSIZ/10)ms duration. */

  datain = (short *) malloc(sizeof(short) * length); //s->length
  dataout = (short *) malloc(sizeof(short) * length); //s->length
  for (i = 0; i < length; i++) { //Snack_GetLength(s)
    datain[i] = (short)s[i]; //(short) Snack_GetSample(s, 0, i);
  }

  if(!this->static_len) {		/* need to create a Hanning FIR? */
    this->static_lcf = (short*)malloc(sizeof(short) * LCSIZ);
    this->static_len = 1 + (LCSIZ/2);
    fn = PI * 2.0 / (LCSIZ - 1);
    scale = 32767.0/(.5 * LCSIZ);
    for(i=0; i < this->static_len; i++) 
      this->static_lcf[i] = (short) (scale * (.5 + (.4 * cos(fn * ((double)i)))));
  }
  do_fir(datain,length,dataout,this->static_len,this->static_lcf,1); /* in downsample.c */ //s->length
  //so = Snack_NewSound(s->samprate, LIN16, s->nchannels);
  MArray_1D(so, length, float, "SC_Feature_Formant.highpass: so");
  if (so == NULL) return(NULL);
  //Snack_ResizeSoundStorage(so, s->length);
  for (i = 0; i < length; i++) { //s->length
    so[i]= (float)dataout[i]; //Snack_SetSample(so, 0, i, (float)dataout[i]);
  }
  //so->length = s->length;
  free((void *)dataout);
  free((void *)datain);
  return(so);
}

//----
//below are help routines from ESPS' sigproc2.c part of the library
//----

//pred anal subroutine with ridge reg
//s - speech
//ls - length of s
//p - pred coefs
//np - polyn order
//c - ref coers
//phi - cov matrix
//shi - cov vect
/*int*/long long SC_Feature_Formant::dlpcwtd(double *s, int *ls, double *p, int *np, double *c, double *phi, double *shi, double *xl, double *w)
{
/*int*/long long m,np1,mm;
double d,pre,pre3,pre2,pre0,pss,pss7,thres;
double ee;
np1  =  *np  +  1;
dcwmtrx(s,np,ls,np,phi,shi,&pss,w);
if(*xl>=1.0e-4)
	{
	pph1 = phi; ppl2 = p + *np;
	for(pp2=p;pp2<ppl2;pp2++){
		*pp2 = *pph1;
		pph1 += np1;
	}
	*ppl2 = pss;
	pss7 = .0000001 * pss;
	mm = dchlsky(phi,np,c,&d);
	if(mm< *np)fprintf(stderr,"LPCHFA error covariance matric rank %d \n",mm);
	dlwrtrn(phi,np,c,shi);
	ee = pss;
	thres = 0.;
	pph1 = phi; pcl = c + mm;
	for(pc2=c;pc2<pcl;pc2++)
		{
		if(*pph1<thres)break;
		ee = ee - *pc2 * *pc2;
		if(ee<thres)break;
		if(ee<pss7)
			fprintf(stderr,"LPCHFA is losing accuracy\n");
		}
	m = pc2 - c;
	if(m != mm)
		fprintf(stderr,"*W* LPCHFA error - inconsistent value of m %d \n",m);
	pre = ee * *xl;
	pphl = phi + *np * *np;
	for(pph1=phi+1;pph1<pphl;pph1+=np1)
		{
		pph2 = pph1;
		for(pph3=pph1+ *np-1;pph3<pphl;pph3+= *np)
			{
			*pph3 = *(pph2++);
			}
		}
	pp2 = p; pre3 = .375 * pre; pre2 = .25 * pre; pre0 = .0625 * pre;
	for(pph1=phi;pph1<pphl;pph1+=np1)
		{
		*pph1 = *(pp2++) + pre3;
		if((pph2=pph1- *np)>phi)
			*pph2 = *(pph1-1) = *pph2 - pre2;
		if((pph3=pph2- *np)>phi)
			*pph3 = *(pph1-2) = *pph3 + pre0;
		}
	*shi -= pre2;
	*(shi+1) += pre0;
	*(p+ *np) = pss + pre3;
	}
m = dcovlpc(phi,shi,p,np,c);
return(m);
}

//covariance LPC analysis; originally from Markel and Gray
//(a translation from the fortran)
int SC_Feature_Formant::w_covar(short *xx, int *m, int n, int istrt, double *y, double *alpha, double *r0, double preemp, int w_type)
{
  //static double *x=NULL;
  //static int nold = 0;
  //static int mold = 0;
  //static double *b = NULL, *beta = NULL, *grc = NULL, *cc = NULL, gam,s;
 int ibeg, ibeg1, ibeg2, ibegmp, np0, ibegm1, msq, np, np1, mf, jp, ip,
     mp, i, j, minc, n1, n2, n3, npb, msub, mm1, isub, m2;
  int mnew = 0;

  if((n+1) > this->static_nold) {
    if(this->static_x) free((void *)this->static_x);
    this->static_x = NULL;
    if(!(this->static_x = (double*)malloc((n+1)*sizeof(double)))) {
      printf("Allocation failure in w_covar()\n");
      return(FALSE);
    }
    memset(this->static_x, 0, (n+1) * sizeof(double));
    this->static_nold = n+1;
  }

  if(*m > this->static_mold) {
    if(this->static_b) free((void *)this->static_b); if(this->static_beta) free((void *)this->static_beta); if (this->static_grc) free((void *)this->static_grc); if (this->static_cc) free((void *)this->static_cc);
    this->static_b = this->static_beta = this->static_grc = this->static_cc = NULL;
    mnew = *m;
    
    if(!((this->static_b = (double*)malloc(sizeof(double)*((mnew+1)*(mnew+1)/2))) && 
       (this->static_beta = (double*)malloc(sizeof(double)*(mnew+3)))  &&
       (this->static_grc = (double*)malloc(sizeof(double)*(mnew+3)))  &&
       (this->static_cc = (double*)malloc(sizeof(double)*(mnew+3)))))   {
      printf("Allocation failure in w_covar()\n");
      return(FALSE);
    }
    this->static_mold = mnew;
  }

  w_window(xx, this->static_x, n, preemp, w_type);  

 ibeg = istrt - 1;
 ibeg1 = ibeg + 1;
 mp = *m + 1;
 ibegm1 = ibeg - 1;
 ibeg2 = ibeg + 2;
 ibegmp = ibeg + mp;
 i = *m;
 msq = ( i + i*i)/2;
 for (i=1; i <= msq; i++) this->static_b[i] = 0.0;
 *alpha = 0.0;
 this->static_cc[1] = 0.0;
 this->static_cc[2] = 0.0;
 for(np=mp; np <= n; np++) {
   np1 = np + ibegm1;
   np0 = np + ibeg;
   *alpha += this->static_x[np0] * this->static_x[np0];
   this->static_cc[1] += this->static_x[np0] * this->static_x[np1];
   this->static_cc[2] += this->static_x[np1] * this->static_x[np1];
 }
 *r0 = *alpha;
 this->static_b[1] = 1.0;
 this->static_beta[1] = this->static_cc[2];
 this->static_grc[1] = -this->static_cc[1]/this->static_cc[2];
 y[0] = 1.0;
 y[1] = this->static_grc[1];
 *alpha += this->static_grc[1]*this->static_cc[1];
 if( *m <= 1) return(FALSE);		/* need to correct indices?? */
 mf = *m;
 for( minc = 2; minc <= mf; minc++) {
   for(j=1; j <= minc; j++) {
     jp = minc + 2 - j;
     n1 = ibeg1 + mp - jp;
     n2 = ibeg1 + n - minc;
     n3 = ibeg2 + n - jp;
     this->static_cc[jp] = this->static_cc[jp - 1] + this->static_x[ibegmp-minc]*this->static_x[n1] - this->static_x[n2]*this->static_x[n3];
   }
   this->static_cc[1] = 0.0;
   for(np = mp; np <= n; np++) {
     npb = np + ibeg;
     this->static_cc[1] += this->static_x[npb-minc]*this->static_x[npb];
   }
   msub = (minc*minc - minc)/2;
   mm1 = minc - 1;
   this->static_b[msub+minc] = 1.0;
   for(ip=1; ip <= mm1; ip++) {
     isub = (ip*ip - ip)/2;
     if(this->static_beta[ip] <= 0.0) {
       *m = minc-1;
       return(TRUE);
     }
     this->static_gam = 0.0;
     for(j=1; j <= ip; j++)
       this->static_gam += this->static_cc[j+1]*this->static_b[isub+j];
     this->static_gam /= this->static_beta[ip];
     for(jp=1; jp <= ip; jp++)
       this->static_b[msub+jp] -= this->static_gam*this->static_b[isub+jp];
   }
   this->static_beta[minc] = 0.0;
   for(j=1; j <= minc; j++)
     this->static_beta[minc] += this->static_cc[j+1]*this->static_b[msub+j];
   if(this->static_beta[minc] <= 0.0) {
     *m = minc-1;
     return(TRUE);
   }
   this->static_s = 0.0;
   for(ip=1; ip <= minc; ip++)
     this->static_s += this->static_cc[ip]*y[ip-1];
   this->static_grc[minc] = -this->static_s/this->static_beta[minc];
   for(ip=1; ip < minc; ip++) {
     m2 = msub+ip;
     y[ip] += this->static_grc[minc]*this->static_b[m2];
   }
   y[minc] = this->static_grc[minc];
   this->static_s = this->static_grc[minc]*this->static_grc[minc]*this->static_beta[minc];
   *alpha -= this->static_s;
   if(*alpha <= 0.0) {
     if(minc < *m) *m = minc;
     return(TRUE);
   }
 }
 return(TRUE);
}

void SC_Feature_Formant::w_window(register short *din, register double *dout, register int n, register double preemp, int type)
{
  switch(type) {
  case 0:
    rwindow(din, dout, n, preemp);
    return;
  case 1:
    hwindow(din, dout, n, preemp);
    return;
  case 2:
    cwindow(din, dout, n, preemp);
    return;
  case 3:
    hnwindow(din, dout, n, preemp);
    return;
  default:
    printf("Unknown window type (%d) requested in w_window()\n",type);
  }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void SC_Feature_Formant::rwindow(register short *din, register double *dout, register int n, register double preemp)
{
  register short *p;
 
/* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
  if(preemp != 0.0) {
    for( p=din+1; n-- > 0; )
      *dout++ = (double)(*p++) - (preemp * *din++);
  } else {
    for( ; n-- > 0; )
      *dout++ =  *din++;
  }
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void SC_Feature_Formant::cwindow(register short *din, register double *dout, register int n, register double preemp)
{
  register int i;
  register short *p;
  //static int wsize = 0;
  //static double *wind=NULL;
  register double *q, co;
 
  if(this->static_wsize != n) {		/* Need to create a new cos**4 window? */
    register double arg, half=0.5;
    
    if(this->static_wind) this->static_wind = (double*)realloc((void *)this->static_wind,n*sizeof(double));
    else this->static_wind = (double*)malloc(n*sizeof(double));
    this->static_wsize = n;
    for(i=0, arg=3.1415927*2.0/(this->static_wsize), q=this->static_wind; i < n; ) {
      co = half*(1.0 - cos((half + (double)i++) * arg));
      *q++ = co * co * co * co;
    }
  }
/* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
  if(preemp != 0.0) {
    for(i=n, p=din+1, q=this->static_wind; i-- > 0; )
      *dout++ = *q++ * ((double)(*p++) - (preemp * *din++));
  } else {
    for(i=n, q=this->static_wind; i-- > 0; )
      *dout++ = *q++ * *din++;
  }
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void SC_Feature_Formant::hwindow(register short *din, register double *dout, register int n, register double preemp)
{
  register int i;
  register short *p;
  //static int wsize = 0;
  //static double *wind=NULL;
  register double *q;

  if(this->static_wsize != n) {		/* Need to create a new Hamming window? */
    register double arg, half=0.5;
    
    if(this->static_wind) this->static_wind = (double*)realloc((void *)this->static_wind,n*sizeof(double));
    else this->static_wind = (double*)malloc(n*sizeof(double));
    this->static_wsize = n;
    for(i=0, arg=3.1415927*2.0/(this->static_wsize), q=this->static_wind; i < n; )
      *q++ = (.54 - .46 * cos((half + (double)i++) * arg));
  }
/* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
  if(preemp != 0.0) {
    for(i=n, p=din+1, q=this->static_wind; i-- > 0; )
      *dout++ = *q++ * ((double)(*p++) - (preemp * *din++));
  } else {
    for(i=n, q=this->static_wind; i-- > 0; )
      *dout++ = *q++ * *din++;
  }
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void SC_Feature_Formant::hnwindow(register short *din, register double *dout, register int n, register double preemp)
{
  register int i;
  register short *p;
  //static int wsize = 0;
  //static double *wind=NULL;
  register double *q;

  if(this->static_wsize != n) {		/* Need to create a new Hamming window? */
    register double arg, half=0.5;
    
    if(this->static_wind) this->static_wind = (double*)realloc((void *)this->static_wind,n*sizeof(double));
    else this->static_wind = (double*)malloc(n*sizeof(double));
    this->static_wsize = n;
    for(i=0, arg=3.1415927*2.0/(this->static_wsize), q=this->static_wind; i < n; )
      *q++ = (half - half * cos((half + (double)i++) * arg));
  }
/* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
  if(preemp != 0.0) {
    for(i=n, p=din+1, q=this->static_wind; i-- > 0; )
      *dout++ = *q++ * ((double)(*p++) - (preemp * *din++));
  } else {
    for(i=n, q=this->static_wind; i-- > 0; )
      *dout++ = *q++ * *din++;
  }
}

//Find the roots of the LPC denominator polynomial and convert the z-plane
//zeros to equivalent resonant frequencies and bandwidths.
//The complex poles are then ordered by frequency.
//lpc_order: order of the LP model
//n_form: number of COMPLEX roots of the LPC polynomial
//init: preset to true if no root candidates are available
//s_freq: the sampling frequency of the speech waveform data
//lpca: linear predictor coefficients
//freq: returned array of candidate formant frequencies
//band: returned array of candidate formant bandwidths
int SC_Feature_Formant::formant(int lpc_order, double s_freq, double *lpca, int *n_form, double *freq, double *band, int init)
{
  double  x, flo, pi2t, theta;
  //static double  rr[MAXORDER], ri[MAXORDER];
  int	i,ii,iscomp1,iscomp2,fc,swit;

  if(init){ /* set up starting points for the root search near unit circle */
    x = sclib::pi/(lpc_order + 1); //M_PI
    for(i=0;i<=lpc_order;i++){
      flo = lpc_order - i;
      this->static_rr[i] = 2.0 * cos((flo + 0.5) * x);
      this->static_ri[i] = 2.0 * sin((flo + 0.5) * x);
    }
  }
  if(! lbpoly(lpca,lpc_order,this->static_rr,this->static_ri)){ /* find the roots of the LPC polynomial */
    *n_form = 0;		/* was there a problem in the root finder? */
    return(FALSE);
  }

  pi2t = sclib::pi * 2.0 /s_freq; //M_PI

  /* convert the z-plane locations to frequencies and bandwidths */
  for(fc=0, ii=0; ii < lpc_order; ii++){
    if((this->static_rr[ii] != 0.0)||(this->static_ri[ii] != 0.0)){
      theta = atan2(this->static_ri[ii],this->static_rr[ii]);
      freq[fc] = fabs(theta / pi2t);
      if((band[fc] = 0.5 * s_freq *
	  log(((this->static_rr[ii] * this->static_rr[ii]) + (this->static_ri[ii] * this->static_ri[ii])))/sclib::pi) < 0.0) //M_PI
	band[fc] = -band[fc];
      fc++;			/* Count the number of real and complex poles. */

      if((this->static_rr[ii] == this->static_rr[ii+1])&&(this->static_ri[ii] == -this->static_ri[ii+1]) /* complex pole? */
	 && (this->static_ri[ii] != 0.0)) ii++; /* if so, don't duplicate */
    }
  }


  /* Now order the complex poles by frequency.  Always place the (uninteresting)
     real poles at the end of the arrays. 	*/
  theta = s_freq/2.0;		/* temporarily hold the folding frequency. */
  for(i=0; i < fc -1; i++){	/* order the poles by frequency (bubble) */
    for(ii=0; ii < fc -1 -i; ii++){
      /* Force the real poles to the end of the list. */
      iscomp1 = (freq[ii] > 1.0) && (freq[ii] < theta);
      iscomp2 = (freq[ii+1] > 1.0) && (freq[ii+1] < theta);
      swit = (freq[ii] > freq[ii+1]) && iscomp2 ;
      if(swit || (iscomp2 && ! iscomp1)){
	flo = band[ii+1];
	band[ii+1] = band[ii];
	band[ii] = flo;
	flo = freq[ii+1];
	freq[ii+1] = freq[ii];
	freq[ii] = flo;
      }
    }
  }
  /* Now count the complex poles as formant candidates. */
  for(i=0, theta = theta - 1.0, ii=0 ; i < fc; i++)
    if( (freq[i] > 1.0) && (freq[i] < theta) ) ii++;
  *n_form = ii;

  return(TRUE);
}

/*		lbpoly.c		*/
/*					*/
/* A polynomial root finder using the Lin-Bairstow method (outlined
	in R.W. Hamming, "Numerical Methods for Scientists and
	Engineers," McGraw-Hill, 1962, pp 356-359.)		*/

#define MAX_ITS	100	/* Max iterations before trying new starts */
#define MAX_TRYS	100	/* Max number of times to try new starts */
#define MAX_ERR		1.e-6	/* Max acceptable error in quad factor */

//find x, where a*x**2 + b*x + c = 0
//r1r, r2r, r1i, r2i: return real and imag. parts of roots
int SC_Feature_Formant::qquad(double a, double b, double c, double *r1r, double *r1i, double *r2r, double *r2i) 
{
double  numi;
double  den, y;

	if(a == 0.0){
		if(b == 0){
		   printf("Bad coefficients to _quad().\n");
		   return(FALSE);
		}
		*r1r = -c/b;
		*r1i = *r2r = *r2i = 0;
		return(TRUE);
	}
	numi = b*b - (4.0 * a * c);
	if(numi >= 0.0) {
		/*
		 * Two forms of the quadratic formula:
		 *  -b + sqrt(b^2 - 4ac)           2c
		 *  ------------------- = --------------------
		 *           2a           -b - sqrt(b^2 - 4ac)
		 * The r.h.s. is numerically more accurate when
		 * b and the square root have the same sign and
		 * similar magnitudes.
		 */
		*r1i = *r2i = 0.0;
		if(b < 0.0) {
			y = -b + sqrt(numi);
			*r1r = y / (2.0 * a);
			*r2r = (2.0 * c) / y;
		}
		else {
			y = -b - sqrt(numi);
			*r1r = (2.0 * c) / y;
			*r2r = y / (2.0 * a);
		}
		return(TRUE);
	}
	else {
		den = 2.0 * a;
		*r1i = sqrt( -numi )/den;
		*r2i = -*r1i;
		*r2r = *r1r = -b/den;
		return(TRUE);
	}
}

//block by thilo:
#ifndef _MSC_VER
	#ifndef DBL_MAX
		#define DBL_MAX std::numeric_limits<double>::max()
	#endif
#endif
//end by thilo

//Rootr and rooti are assumed to contain starting points for the root search on entry to lbpoly().
//return FALSE on error
//a: coeffs. of the polynomial (increasing order)
//order: the order of the polynomial
//rootr, rooti: the real and imag. roots of the polynomial
int SC_Feature_Formant::lbpoly(double *a, int order, double *rootr, double *rooti) 
{
    int	    ord, ordm1, ordm2, itcnt, i, k, mmk, mmkp2, mmkp1, ntrys;
    double  err, p, q, delp, delq, b[MAXORDER], c[MAXORDER], den;
    double  lim0 = 0.5*sqrt(DBL_MAX);

    for(ord = order; ord > 2; ord -= 2){
	ordm1 = ord-1;
	ordm2 = ord-2;
	/* Here is a kluge to prevent UNDERFLOW! (Sometimes the near-zero
	   roots left in rootr and/or rooti cause underflow here...	*/
	if(fabs(rootr[ordm1]) < 1.0e-10) rootr[ordm1] = 0.0;
	if(fabs(rooti[ordm1]) < 1.0e-10) rooti[ordm1] = 0.0;
	p = -2.0 * rootr[ordm1]; /* set initial guesses for quad factor */
	q = (rootr[ordm1] * rootr[ordm1]) + (rooti[ordm1] * rooti[ordm1]);
	for(ntrys = 0; ntrys < MAX_TRYS; ntrys++)
	{
	    int	found = FALSE;

	    for(itcnt = 0; itcnt < MAX_ITS; itcnt++)
	    {
		double	lim = lim0 / (1 + fabs(p) + fabs(q));

		b[ord] = a[ord];
		b[ordm1] = a[ordm1] - (p * b[ord]);
		c[ord] = b[ord];
		c[ordm1] = b[ordm1] - (p * c[ord]);
		for(k = 2; k <= ordm1; k++){
		    mmk = ord - k;
		    mmkp2 = mmk+2;
		    mmkp1 = mmk+1;
		    b[mmk] = a[mmk] - (p* b[mmkp1]) - (q* b[mmkp2]);
		    c[mmk] = b[mmk] - (p* c[mmkp1]) - (q* c[mmkp2]);
		    if (b[mmk] > lim || c[mmk] > lim)
			break;
		}
		if (k > ordm1) { /* normal exit from for(k ... */
		    /* ????	b[0] = a[0] - q * b[2];	*/
		    b[0] = a[0] - p * b[1] - q * b[2];
		    if (b[0] <= lim) k++;
		}
		if (k <= ord)	/* Some coefficient exceeded lim; */
		    break;	/* potential overflow below. */

		err = fabs(b[0]) + fabs(b[1]);

		if(err <= MAX_ERR) {
		    found = TRUE;
		    break;
		}

		den = (c[2] * c[2]) - (c[3] * (c[1] - b[1]));  
		if(den == 0.0)
		    break;

		delp = ((c[2] * b[1]) - (c[3] * b[0]))/den;
		delq = ((c[2] * b[0]) - (b[1] * (c[1] - b[1])))/den;  

		/* printf("\nerr=%f  delp=%f  delq=%f  p=%f  q=%f",
		   err,delp,delq,p,q); */

		p += delp;
		q += delq;

	    } /* for(itcnt... */

	    if (found)		/* we finally found the root! */
		break;
	    else { /* try some new starting values */
		p = ((double)rand() - 0.5*RAND_MAX)/(double)RAND_MAX;
		q = ((double)rand() - 0.5*RAND_MAX)/(double)RAND_MAX;
		/* fprintf(stderr, "\nTried new values: p=%f  q=%f\n",p,q); */
	    }

	} /* for(ntrys... */
	if((itcnt >= MAX_ITS) && (ntrys >= MAX_TRYS)){
	  /*	    printf("Exceeded maximum trial count in _lbpoly.\n");*/
	    return(FALSE);
	}

	if(!qquad(1.0, p, q,
		  &rootr[ordm1], &rooti[ordm1], &rootr[ordm2], &rooti[ordm2]))
	    return(FALSE);

	/* Update the coefficient array with the coeffs. of the
	   reduced polynomial. */
	for( i = 0; i <= ordm2; i++) a[i] = b[i+2];
    }

    if(ord == 2){		/* Is the last factor a quadratic? */
	if(!qquad(a[2], a[1], a[0],
		  &rootr[1], &rooti[1], &rootr[0], &rooti[0]))
	    return(FALSE);
	return(TRUE);
    }
    if(ord < 1) {
	printf("Bad ORDER parameter in _lbpoly()\n");
	return(FALSE);
    }

    if( a[1] != 0.0) rootr[0] = -a[0]/a[1];
    else {
	rootr[0] = 100.0;	/* arbitrary recovery value */
	printf("Numerical problems in lbpoly()\n");
    }
    rooti[0] = 0.0;

    return(TRUE);
}

//cov mat for wtd lpc
void SC_Feature_Formant::dcwmtrx(double *s, int* ni, int *nl, int *np, double *phi, double *shi, double *ps, double *w)
{
	double sm;
	int i,j;
	*ps = 0;
	for(this->static_pdl1=s+*ni,this->static_pdl2=w,this->static_pdll=s+*nl;this->static_pdl1<this->static_pdll;this->static_pdl1++,this->static_pdl2++)
		*ps += *this->static_pdl1 * *this->static_pdl1 * *this->static_pdl2;

	for(this->static_pdl3=shi,this->static_pdl4=shi+*np,this->static_pdl5=s+*ni;this->static_pdl3<this->static_pdl4;this->static_pdl3++,this->static_pdl5--){
		*this->static_pdl3 = 0.;
		for(this->static_pdl1=s+*ni,this->static_pdl2=w,this->static_pdll=s+*nl,this->static_pdl6=this->static_pdl5-1;
			this->static_pdl1<this->static_pdll;this->static_pdl1++,this->static_pdl2++,this->static_pdl6++)
			*this->static_pdl3 += *this->static_pdl1 * *this->static_pdl6 * *this->static_pdl2;

	}

	for(i=0;i<*np;i++)
		for(j=0;j<=i;j++){
			sm = 0.;
			for(this->static_pdl1=s+*ni-i-1,this->static_pdl2=s+*ni-j-1,this->static_pdl3=w,this->static_pdll=s+*nl-i-1;
				this->static_pdl1<this->static_pdll;)
				sm += *this->static_pdl1++ * *this->static_pdl2++ * *this->static_pdl3++;

			*(phi + *np * i + j) = sm;
			*(phi + *np * j + i) = sm;
		}
}

//performs cholesky decomposition
//a - l * l(transpose)
//l - lower triangle
//det det(a)
//a - nxn matrix
//return - no of reduced elements
//results in lower half + diagonal
//upper half undisturbed.
int SC_Feature_Formant::dchlsky(double *a, int *n, double *t, double *det)
{
	double sm;
	int m;
	*det = 1.;
	m = 0;
	pal = a + *n * *n;
	for(this->pa_1=a;this->pa_1<pal;this->pa_1+= *n){
		this->pa_3=this->pa_1;
		this->pt = t;
		for(this->pa_2=a;this->pa_2<=this->pa_1;this->pa_2+= *n){
			sm = *this->pa_3;	/*a(i,j)*/
			this->pa_5 = this->pa_2;
			for(this->pa_4=this->pa_1;this->pa_4<this->pa_3;this->pa_4++)
				sm =  sm - *this->pa_4 * *(this->pa_5++);
			if(pa_1==pa_2){
				if(sm<=0.)return(m);
				*this->pt = sqrt(sm);
				*det = *det * *this->pt;
				*(this->pa_3++) = *this->pt;
				m++;
				*this->pt = 1. / *this->pt;
				this->pt++;
			}
			else
				*(this->pa_3++) = sm * *(this->pt++);
		}
	}
	return(m);
}

//routine to solve ax=y with cholesky
//a - nxn matrix
//x,y -vectors
void SC_Feature_Formant::dlwrtrn(double *a, int *n, double *x, double *y)
{
	double sm;
	*x = *y / *a;
	this->pxl = x + 1;
	this->pyl = y + *n;
	this->pa = a + *n;
	for(this->py=y+1;this->py<this->pyl;this->py++,this->pxl++){
		sm = *this->py;
		this->pa1 = this->pa;
		for(this->px=x;this->px<this->pxl;this->px++)
			sm = sm - *(this->pa1++) * *this->px;
		pa += *n;
		*this->px = sm / *this->pa1;
	}
}

//solve p*a=s using stabilized covariance method
//p - cov nxn matrix
//s - corrvec
//a lpc coef *a = 1.
//c - ref coefs
/*int*/long long SC_Feature_Formant::dcovlpc(double *p, double *s, double *a, int *n, double *c)
{
	double ee;
	double ps,ps1,thres,d;
	/*int*/long long m,n1;
	m = dchlsky(p,n,c,&d);
	dlwrtrn(p,n,c,s);
	thres = 1.0e-31;
	n1 = *n + 1;
	ps = *(a + *n);
	ps1 = 1.e-8*ps;
	this->static_ppl = p + *n * m;
	m = 0;
	for(this->static_pp=p;this->static_pp<this->static_ppl;this->static_pp+=n1){
		if(*this->static_pp<thres)break;
		m++;
	}
	ee = ps;
	this->static_ppl = c + m; pa = a;
	for(this->static_pp=c;this->static_pp<this->static_ppl;this->static_pp++){
		ee = ee - *this->static_pp * *this->static_pp;
		if(ee<thres)break;
		if(ee<ps1)fprintf(stderr,"*w* covlpc is losing accuracy\n");
		*(this->static_pa++) = sqrt(ee);
	}
	m = this->static_pa - a;
	*c = - *c/sqrt(ps);
	this->static_ppl = c + m; this->static_pa = a;
	for(this->static_pp=c+1;this->static_pp<this->static_ppl;this->static_pp++)
		*this->static_pp = - *this->static_pp / *(this->static_pa++);
	dreflpc(c,a,&m);
	this->static_ppl = a + *n;
	for(this->static_pp=a+m+1;this->static_pp<=this->static_ppl;this->static_pp++)*this->static_pp=0.;
	return(m);
}

//convert ref to lpc
//c - ref
//a - polyn
//n - no of coef
void SC_Feature_Formant::dreflpc(double *c, double *a, /*int*/long long *n)
{
double ta1;
*a = 1.;
*(a+1) = *c;
this->static_pc = c; this->static_pa2 = a+ *n;
for(this->static_pa1=a+2;this->static_pa1<=this->static_pa2;this->static_pa1++)
	{
	this->pc++;
	*this->static_pa1 = *this->static_pc;
	this->static_pa5 = a + (this->static_pa1-a)/2;
	this->static_pa4 = this->static_pa1 - 1;
	for(this->static_pa3=a+1;this->static_pa3<=this->static_pa5;this->static_pa3++,this->static_pa4--)
		{
		ta1 = *this->static_pa3 + *this->static_pc * *this->static_pa4;
		*this->static_pa4 = *this->static_pa4 + *this->static_pa3 * *this->static_pc;
		*this->static_pa3 = ta1;
		}
	}
}

int SC_Feature_Formant::lpc(int lpc_ord, double lpc_stabl, int wsize, short *data, double *lpca, double *ar, double *lpck, double *normerr, double *rms, double preemp, int type)
{
  //static double *dwind=NULL;
  //static int nwind=0;
  double rho[MAXORDER+1], k[MAXORDER], a[MAXORDER+1],*r,*kp,*ap,en,er;
  double wfact = 1.0;

  if((wsize <= 0) || (!data) || (lpc_ord > MAXORDER)) return(FALSE);
  
  if(this->static_nwind != wsize) {
    if(this->static_dwind) this->static_dwind = (double*)realloc((void *)this->static_dwind,wsize*sizeof(double));
    else this->static_dwind = (double*)malloc(wsize*sizeof(double));
    if(!this->static_dwind) {
      printf("Can't allocate scratch memory in lpc()\n");
      return(FALSE);
    }
    this->static_nwind = wsize;
  }
  
  w_window(data, this->static_dwind, wsize, preemp, type);
  if(!(r = ar)) r = rho;
  if(!(kp = lpck)) kp = k;
  if(!(ap = lpca)) ap = a;
  autoc( wsize, this->static_dwind, lpc_ord, r, &en );
  if(lpc_stabl > 1.0) { /* add a little to the diagonal for stability */
    int i;
    double ffact;
    ffact =1.0/(1.0 + exp((-lpc_stabl/20.0) * log(10.0)));
    for(i=1; i <= lpc_ord; i++) rho[i] = ffact * r[i];
    *rho = *r;
    r = rho;
    if(ar)
      for(i=0;i<=lpc_ord; i++) ar[i] = r[i]; /* copy out for possible use later */
  }
  durbin ( r, kp, &ap[1], lpc_ord, &er);

/* After talkin with David T., I decided to take the following correction
out -- if people window, the resulting output (spectrum and power) should
correspond to the windowed data, and we shouldn't try to back-correct
to un-windowed data. */

/*  switch(type) {		/ rms correction for window */
/*   case 0:
    wfact = 1.0;		/ rectangular */
/*    break;
  case 1:
    wfact = .630397;		/ Hamming */
/*    break;
  case 2:
    wfact = .443149;		/ (.5 - .5*cos)^4 */
/*    break;
  case 3:
    wfact = .612372;		/ Hanning */
/*    break;
  }
*/

  *ap = 1.0;
  if(rms) *rms = en/wfact;
  if(normerr) *normerr = er;
  return(TRUE);
}

//Compute the pp+1 autocorrelation lags of the windowsize samples in s.
//Return the normalized autocorrelation coefficients in r.
//The rms is returned in e.
void SC_Feature_Formant::autoc(register int windowsize, register double *s, register int p, register double *r, register double *e)
{
  register int i, j;
  register double *q, *t, sum, sum0;

  for ( i=0, q=s, sum0=0.; i< windowsize; q++, i++){
	sum0 += (*q) * (*q);
  }
  *r = 1.;  /* r[0] will always =1. */
  if ( sum0 == 0.){   /* No energy: fake low-energy white noise. */
  	*e = 1.;   /* Arbitrarily assign 1 to rms. */
		   /* Now fake autocorrelation of white noise. */
	for ( i=1; i<=p; i++){
		r[i] = 0.;
	}
	return;
  }
  for( i=1; i <= p; i++){
	for( sum=0., j=0, q=s, t=s+i; j < (windowsize)-i; j++, q++, t++){
		sum += (*q) * (*t);
	}
	*(++r) = sum/sum0;
  }
  if(sum0 < 0.0) printf("lpcfloat:autoc(): sum0 = %f\n",sum0);
  *e = sqrt(sum0/windowsize);
}

//Compute the AR and PARCOR coefficients using Durbin's recursion. 
//Note: Durbin returns the coefficients in normal sign format.
//(i.e. a[0] is assumed to be = +1.)
void SC_Feature_Formant::durbin(register double *r, register double *k, register double *a, register int p, register double *ex)
{
  double b[MAXORDER];
  register int i, j;
  register double e, s;

    e = *r;
    *k = -r[1]/e;
    *a = *k;
    e *= (1. - (*k) * (*k));
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
	e *= ( 1. - (k[i] * k[i]) );
    }
    *ex = e;
}

#endif
