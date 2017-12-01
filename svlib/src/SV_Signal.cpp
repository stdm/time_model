//************************************************************************
//    Load RAW speech signal into memory
//
//
//    Author  : Jialong HE
//    Date    : March 11, 1999
//************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include "SV_Signal.h"
#include "SV_Error.h"

using namespace std;

static char SV_LibID[] = "Copyright (c) by Jialong He";

//==========================================
// default constructor
//==========================================
SV_Signal::SV_Signal() {
	Buf_L            = NULL;
	Buf_R            = NULL;
	Len              = 0;
	SigPar.SRate     = 16000;
	SigPar.Encode    = PCM_2;
	SigPar.NChannel  = 1;
	SigPar.StByte    = 0;
	ClassSig         = 1964121;
}


//==========================================
// destructor
//==========================================
SV_Signal::~SV_Signal() {

	if (Buf_L != NULL) {
		MFree_1D (Buf_L);
	}

	if (Buf_R != NULL) {
		MFree_1D (Buf_R);
	}

	ClassSig         = 0;

}

//=======================================================*/
//  Test if current class is valid        
//  return 1: yes, it is valid             
//  return 0: no, it is not valid, maybe deleted.
//=======================================================*/
int SV_Signal::Valid (void) {
	if (ClassSig == 1964121) {return (1);}
	else {return(0);}
}


//=======================================================*/
//  Test machine's byte order                            */
//  return 1: BigEndian (e.g., Motorola or SUN)          */
//  return 0: LittleEndian (Intel)                       */
//=======================================================*/
int SV_Signal::IsBigEndian (void) {

    short Sample = 1;
    unsigned char *BytePnt;

    BytePnt = (unsigned char *)&Sample;

    if (*BytePnt == 0) {
		return(1);
	}
    else {return(0);}
};

//=======================================================*/
//  Swap byte order in Buf (short integer)               */
//=======================================================*/
void SV_Signal::SwapByteOrder (void) {

	if ( (Len != 0) && (Buf_L != NULL)) {
       //swab((char*)Buf_L, (char*)Buf_L, Len*2);
		this->io.swapBytes((char*)this->Buf_L, 2, this->Len*2); //by thilo: os-independant
    }

	if ( (Len != 0) && (Buf_R != NULL)) {
       //swab((char*)Buf_R, (char*)Buf_R, Len*2);
		this->io.swapBytes((char*)this->Buf_L, 2, this->Len*2); //by thilo: os-independant
    }
};

//=======================================================*/
//  Find out diskfile size in bytes                      */
//=======================================================*/
long SV_Signal::FileSize (char *FName) {

	ifstream InFile(FName);
    streampos EndPos;
	InFile.seekg(0, ios::end);
	EndPos = InFile.tellg();
    InFile.close();
	return(EndPos);
}

//=======================================================*/
//  Load signal into memory                              */
//=======================================================*/
long SV_Signal::LoadSignal (char *FName) {

	ifstream InFile;
	long NBytes, NSamples;
	char *pByte;
	int  Cnt;

	//--------------------------------------
	//  Currently, only support one channel 
	//--------------------------------------
	if (SigPar.NChannel != 1) {
 	   REPORT_ERROR(SVLIB_BadData, "Only one channel supported!");   
	}

	NBytes   = FileSize(FName) - SigPar.StByte;
	NSamples = NBytes;    // default PCM_ulaw

	if (SigPar.Encode == PCM_2) {
		NSamples = NBytes >> 1;   
	}

	//------------------------------------
	// Allocate memory for signal
	//------------------------------------
	if (Len != NSamples) {
		MFree_1D(Buf_L);
		MArray_1D(Buf_L, NSamples, short, "Buf_L"); 
		Len = NSamples;
	}

	pByte = (char*)Buf_L;
	//------------------------------------
	// Read signal into memory
	//------------------------------------
	InFile.open(FName, ios::in|ios::binary); //InFile.open(FName, ios::in|ios::nocreate|ios::binary);
	InFile.seekg(SigPar.StByte, ios::beg);
	InFile.read(pByte, NBytes);
	if (InFile.fail()) {
		REPORT_ERROR(SVLIB_FileErr, "LoadSignal"); 	
	}
	InFile.close();
    
	//------------------------------------
	// for ulaw, convert samples to short
	// from last sample --> begining
	//------------------------------------
	if (SigPar.Encode == ULAW) {
		for (Cnt=NSamples-1; Cnt>=0; Cnt--) {
			Buf_L[Cnt] = (short)ulaw2linear(pByte[Cnt]);
		}
		
	}

	return(Len);
}

//=======================================================*/
//  Save signal to diskfile                              */
//=======================================================*/
long SV_Signal::SaveSignal (char *FName) {

	ofstream OutFile;

	OutFile.open(FName, ios::out|ios::binary);
//	OutFile.seekg(SigPar.StByte, ios::beg);
	OutFile.write((char*)Buf_L, Len*sizeof(short));
	if (OutFile.fail()) {
		REPORT_ERROR(SVLIB_FileErr, "SaveSignal"); 	
	}

	OutFile.close();
	return (Len*2);
}



/*=====================================================================*/
/* This procedure implements endpoint detection algorithm proposed by  */
/* Rabiner & Sambur (See AT&T Technical Journal, Vol. 54, No. 2,       */
/* Feb. 1975, pp. 297-31)                                              */
/*                                                                     */
/*  Example:                                                           */
/*      short Buf[8192];                                               */
/*      int Start, End;                                                */
/*       (load speech samples)                                         */
/*      endpoint(Buf, 8192, &Start, &End);                             */
/*---------------------------------------------------------------------*/
/* Author : Jialong He                                                 */
/*  Date  : April 25, 1995                                             */
/*                                                                     */
/*=====================================================================*/


/*=========================================================*/
/* find endpoint from given Energy and Zero contours       */
/*=========================================================*/
int search_endpoint(float *Energy, int *Zero, int FrameNum,
	     double EnergyThresL, double EnergyThresH, int ZeroThres) {

    int FrameInd, AboveInd, CurrEndPoint;
    int AboveLow=0, AboveZero=0;

    /*-----------------------------------------------*/
    /* decide endpoint from energy contour           */
    /*-----------------------------------------------*/
    for (FrameInd=0; FrameInd<FrameNum; FrameInd++) {

      if (Energy[FrameInd] > EnergyThresL)
	if (!AboveLow) {            /* exceed low threshold */
	  AboveLow = 1;
	  AboveInd = 0;
	  CurrEndPoint = FrameInd;
	}

	else
	  if ((Energy[FrameInd] > EnergyThresH) && (AboveInd < 10))
	      break;        /* find the end point */
	  else  if (AboveInd < 10) AboveInd++;
		else AboveLow = 0;   /*state between low and high too long*/

      else AboveLow = 0;

    };      /* for FrameCnt */

    /*-----------------------------------------------*/
    /* fine tune CurrEndPoint based on Zero contour  */
    /*-----------------------------------------------*/
    for (FrameInd=CurrEndPoint; FrameInd>0; FrameInd--)
      if (Zero[FrameInd] > ZeroThres) AboveZero++;
      else break;

   if (AboveZero>3) CurrEndPoint -=AboveZero;   /* move endpoint ahead */
   return(CurrEndPoint);
};

/*=========================================================*/
/* call this procedure for start and end point             */
/*=========================================================*/
void SV_Signal::EndPoint(short *SpeechBuf, int SampleNum, int *Start, int *End) {

   short ThisSample, LastSample;
   int   StartFrame, EndFrame;
   int   FrameNum, WinSize, StepSize;
   int   FrameCnt, SampleCnt;
   float CurrE, *Energy;                    /* energy contour */
   float MaxEnergy, MinEnergy;              /* high and low peak of energy */
   double EnergyThresH, EnergyThresL;        /* two thresholds for energy */
   int   CurrZ, *Zero, ZeroThres;           /* zero crossing contour */

   WinSize  = SigPar.SRate / 100;           /* 10 ms */
   StepSize = WinSize / 2;                  /* overlay half */
   FrameNum = (SampleNum - WinSize) / StepSize + 1;
   //----------------------------------
   // Too short, can not do it.
   //----------------------------------
   if (FrameNum<10) {
	 *Start = 0;
	 *End   = SampleNum;
	 return;	   
   }

   MArray_1D(Energy, FrameNum, float, "Energy");      /* allo memory */
   MArray_1D(Zero,   FrameNum, int,   "Zero");

   /*------------------------------------------*/
   /* Calculate Energy and Zero Crossing Rate  */
   /*------------------------------------------*/
   for (FrameCnt=0; FrameCnt<FrameNum; FrameCnt++) {
     CurrE = 0.0;
     CurrZ = 0;
     for (SampleCnt=1; SampleCnt<WinSize; SampleCnt++) {
	ThisSample = SpeechBuf[FrameCnt * StepSize + SampleCnt];
	LastSample = SpeechBuf[FrameCnt * StepSize + SampleCnt - 1];
	CurrE +=(float) abs(ThisSample);  /* Sum ABS as energy */
	if (((ThisSample>0) && (LastSample<0)) || ((ThisSample<0) && (LastSample>0)))
	   CurrZ++;           /* Sign changed, a zero crossing found */
     }
     Energy[FrameCnt] = CurrE;
     Zero[FrameCnt]   = CurrZ;
   }  /* for FrameCnt */

   /*----------------------------------------------*/
   /* caculate ZeroThres from the first 10 frame   */
   /*----------------------------------------------*/
   ZeroThres = 0;
   for (FrameCnt=0; FrameCnt<10; FrameCnt++)
     if (FrameCnt<FrameNum)    ZeroThres += Zero[FrameCnt];
   ZeroThres = ZeroThres / 10;         /* Note: only mean, deviation is ignored */
   if (ZeroThres > 25) ZeroThres = 25;   /*max 25 zero crossing per 10 ms*/

   /*----------------------------------------------*/
   /* caculate low and high energy thresholds      */
   /*----------------------------------------------*/
   MaxEnergy = 0;
   MinEnergy = LARGE;

   for (FrameCnt=0; FrameCnt<FrameNum; FrameCnt++) {
     if (Energy[FrameCnt] > MaxEnergy) MaxEnergy = Energy[FrameCnt];
     if (Energy[FrameCnt] < MinEnergy) MinEnergy = Energy[FrameCnt];
   }

   EnergyThresL = 0.03 * (MaxEnergy - MinEnergy) + MinEnergy;
   if (EnergyThresL > 4.0 * MinEnergy)  EnergyThresL = 4.0 * MinEnergy;
   EnergyThresH  = 5.0 * EnergyThresL;

   /*------------------------------------*/
   /* if very low SNR                    */
   /*------------------------------------*/
   if (EnergyThresH > MaxEnergy/3.0) EnergyThresH = MaxEnergy / 3.0;

   StartFrame = search_endpoint(Energy, Zero, FrameNum, EnergyThresL, EnergyThresH, ZeroThres);
   /*--------------------------------------------*/
   /* reverse frame sequence for end point       */
   /*--------------------------------------------*/
   for (FrameCnt = 0; FrameCnt<FrameNum/2; FrameCnt++) {
      CurrE = Energy[FrameCnt];
      CurrZ = Zero[FrameCnt];
      Energy[FrameCnt] = Energy[FrameNum - FrameCnt - 1];
      Zero[FrameCnt]   = Zero[FrameNum - FrameCnt - 1];
      Energy[FrameNum - FrameCnt -1] = CurrE;
      Zero[FrameNum - FrameCnt - 1]  = CurrZ;
   }

   EndFrame = search_endpoint(Energy, Zero, FrameNum, EnergyThresL, EnergyThresH, ZeroThres);
   /*--------------------------------------------*/
   /* convert frame number to sample number      */
   /*--------------------------------------------*/
   *Start = StartFrame * StepSize;
   *End   = (FrameNum - EndFrame) * StepSize;


   MFree_1D(Energy);
   MFree_1D(Zero);
}




/*
** This routine converts from linear to ulaw.
**
** Craig Reese: IDA/Supercomputing Research Center
** Joe Campbell: Department of Defense
** 29 September 1989
**
** References:
** 1) CCITT Recommendation G.711  (very difficult to follow)
** 2) "A New Digital Technique for Implementation of Any
**     Continuous PCM Companding Law," Villeret, Michel,
**     et al. 1973 IEEE Int. Conf. on Communications, Vol 1,
**     1973, pg. 11.12-11.17
** 3) MIL-STD-188-113,"Interoperability and Performance Standards
**     for Analog-to_Digital Conversion Techniques,"
**     17 February 1987
**
** Input: Signed 16 bit linear sample
** Output: 8 bit ulaw sample
*/

#define ZEROTRAP    /* turn on the trap as per the MIL-STD */
#undef ZEROTRAP
#define BIAS 0x84   /* define the add-in bias for 16 bit samples */
#define CLIP 32635

unsigned char SV_Signal::linear2ulaw(int sample) {

  static int exp_lut[256] = {0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,
                             4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                             5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                             5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                             6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                             6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                             6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                             6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                             7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
  int sign, exponent, mantissa;
  unsigned char ulawbyte;

  /* Get the sample into sign-magnitude. */
  sign = (sample >> 8) & 0x80;		/* set aside the sign */
  if(sign != 0) sample = -sample;		/* get magnitude */
  if(sample > CLIP) sample = CLIP;		/* clip the magnitude */

  /* Convert from 16 bit linear to ulaw. */
  sample = sample + BIAS;
  exponent = exp_lut[( sample >> 7 ) & 0xFF];
  mantissa = (sample >> (exponent + 3)) & 0x0F;
  ulawbyte = ~(sign | (exponent << 4) | mantissa);
#ifdef ZEROTRAP
  if (ulawbyte == 0) ulawbyte = 0x02;	/* optional CCITT trap */
#endif

  return(ulawbyte);
}

/*
** This routine converts from ulaw to 16 bit linear.
**
** Craig Reese: IDA/Supercomputing Research Center
** 29 September 1989
**
** References:
** 1) CCITT Recommendation G.711  (very difficult to follow)
** 2) MIL-STD-188-113,"Interoperability and Performance Standards
**     for Analog-to_Digital Conversion Techniques,"
**     17 February 1987
**
** Input: 8 bit ulaw sample
** Output: signed 16 bit linear sample
*/

int SV_Signal::ulaw2linear(unsigned char ulawbyte){

  static int exp_lut[8] = { 0, 132, 396, 924, 1980, 4092, 8316, 16764 };
  int sign, exponent, mantissa, sample;

  ulawbyte = ~ulawbyte;
  sign = (ulawbyte & 0x80);
  exponent = (ulawbyte >> 4) & 0x07;
  mantissa = ulawbyte & 0x0F;
  sample = exp_lut[exponent] + (mantissa << (exponent + 3));
  if(sign != 0) sample = -sample;

  return(sample);
}

/*========================================================
 * This source code is a product of Sun Microsystems, Inc., modified
 * by CMU,  and is provided  for unrestricted use.  Users may copy
 * or modify this source code without  charge.
 *
==========================================================*/

/*
 * g711.c
 *
 * u-law, A-law and linear PCM conversions.
 */
#define	SIGN_BIT	(0x80)		/* Sign bit for a A-law byte. */
#define	QUANT_MASK	(0xf)		/* Quantization field mask. */
#define	NSEGS		(8)		/* Number of A-law segments. */
#define	SEG_SHIFT	(4)		/* Left shift for segment number. */
#define	SEG_MASK	(0x70)		/* Segment field mask. */

static short seg_end[8] = {0xFF, 0x1FF, 0x3FF, 0x7FF,
			    0xFFF, 0x1FFF, 0x3FFF, 0x7FFF};

/* copy from CCITT G.711 specifications */
unsigned char _u2a[128] = {			/* u- to A-law conversions */
	1,	1,	2,	2,	3,	3,	4,	4,
	5,	5,	6,	6,	7,	7,	8,	8,
	9,	10,	11,	12,	13,	14,	15,	16,
	17,	18,	19,	20,	21,	22,	23,	24,
	25,	27,	29,	31,	33,	34,	35,	36,
	37,	38,	39,	40,	41,	42,	43,	44,
	46,	48,	49,	50,	51,	52,	53,	54,
	55,	56,	57,	58,	59,	60,	61,	62,
	64,	65,	66,	67,	68,	69,	70,	71,
	72,	73,	74,	75,	76,	77,	78,	79,
	81,	82,	83,	84,	85,	86,	87,	88,
	89,	90,	91,	92,	93,	94,	95,	96,
	97,	98,	99,	100,	101,	102,	103,	104,
	105,	106,	107,	108,	109,	110,	111,	112,
	113,	114,	115,	116,	117,	118,	119,	120,
	121,	122,	123,	124,	125,	126,	127,	128};

unsigned char _a2u[128] = {			/* A- to u-law conversions */
	1,	3,	5,	7,	9,	11,	13,	15,
	16,	17,	18,	19,	20,	21,	22,	23,
	24,	25,	26,	27,	28,	29,	30,	31,
	32,	32,	33,	33,	34,	34,	35,	35,
	36,	37,	38,	39,	40,	41,	42,	43,
	44,	45,	46,	47,	48,	48,	49,	49,
	50,	51,	52,	53,	54,	55,	56,	57,
	58,	59,	60,	61,	62,	63,	64,	64,
	65,	66,	67,	68,	69,	70,	71,	72,
	73,	74,	75,	76,	77,	78,	79,	79,
	80,	81,	82,	83,	84,	85,	86,	87,
	88,	89,	90,	91,	92,	93,	94,	95,
	96,	97,	98,	99,	100,	101,	102,	103,
	104,	105,	106,	107,	108,	109,	110,	111,
	112,	113,	114,	115,	116,	117,	118,	119,
	120,	121,	122,	123,	124,	125,	126,	127};

static int
search(
	int		val,
	short		*table,
	int		size)
{
	int		i;

	for (i = 0; i < size; i++) {
		if (val <= *table++)
			return (i);
	}
	return (size);
}

/*
 * linear2alaw() - Convert a 16-bit linear PCM value to 8-bit A-law
 *
 * linear2alaw() accepts an 16-bit integer and encodes it as A-law data.
 *
 *		Linear Input Code	Compressed Code
 *	------------------------	---------------
 *	0000000wxyza			000wxyz
 *	0000001wxyza			001wxyz
 *	000001wxyzab			010wxyz
 *	00001wxyzabc			011wxyz
 *	0001wxyzabcd			100wxyz
 *	001wxyzabcde			101wxyz
 *	01wxyzabcdef			110wxyz
 *	1wxyzabcdefg			111wxyz
 *
 * For further information see John C. Bellamy's Digital Telephony, 1982,
 * John Wiley & Sons, pps 98-111 and 472-476.
 */
unsigned char SV_Signal::linear2alaw(int pcm_val) {
	/* 2's complement (16-bit range) */

	int		mask;
	int		seg;
	unsigned char	aval;

	if (pcm_val >= 0) {
		mask = 0xD5;		/* sign (7th) bit = 1 */
	} else {
		mask = 0x55;		/* sign bit = 0 */
		pcm_val = -pcm_val - 8;
	}

	/* Convert the scaled magnitude to segment number. */
	seg = search(pcm_val, seg_end, 8);

	/* Combine the sign, segment, and quantization bits. */

	if (seg >= 8)		/* out of range, return maximum value. */
		return (0x7F ^ mask);
	else {
		aval = seg << SEG_SHIFT;
		if (seg < 2)
			aval |= (pcm_val >> 4) & QUANT_MASK;
		else
			aval |= (pcm_val >> (seg + 3)) & QUANT_MASK;
		return (aval ^ mask);
	}
}

/*
 * alaw2linear() - Convert an A-law value to 16-bit linear PCM
 *
 */
int SV_Signal::alaw2linear(unsigned char a_val) {

	int		t;
	int		seg;

	a_val ^= 0x55;

	t = (a_val & QUANT_MASK) << 4;
	seg = ((unsigned)a_val & SEG_MASK) >> SEG_SHIFT;
	switch (seg) {
	case 0:
		t += 8;
		break;
	case 1:
		t += 0x108;
		break;
	default:
		t += 0x108;
		t <<= seg - 1;
	}
	return ((a_val & SIGN_BIT) ? t : -t);
}

#ifdef  AnotherIMP //===============================

#define	BIAS		(0x84)		/* Bias for linear code. */

/*
 * linear2ulaw() - Convert a linear PCM value to u-law
 *
 * In order to simplify the encoding process, the original linear magnitude
 * is biased by adding 33 which shifts the encoding range from (0 - 8158) to
 * (33 - 8191). The result can be seen in the following encoding table:
 *
 *	Biased Linear Input Code	Compressed Code
 *	------------------------	---------------
 *	00000001wxyza			000wxyz
 *	0000001wxyzab			001wxyz
 *	000001wxyzabc			010wxyz
 *	00001wxyzabcd			011wxyz
 *	0001wxyzabcde			100wxyz
 *	001wxyzabcdef			101wxyz
 *	01wxyzabcdefg			110wxyz
 *	1wxyzabcdefgh			111wxyz
 *
 * Each biased linear code has a leading 1 which identifies the segment
 * number. The value of the segment number is equal to 7 minus the number
 * of leading 0's. The quantization interval is directly available as the
 * four bits wxyz.  * The trailing bits (a - h) are ignored.
 *
 * Ordinarily the complement of the resulting code word is used for
 * transmission, and so the code word is complemented before it is returned.
 *
 * For further information see John C. Bellamy's Digital Telephony, 1982,
 * John Wiley & Sons, pps 98-111 and 472-476.
 */
unsigned char linear2ulaw( int	pcm_val) {
	/* 2's complement (16-bit range) */
	int		mask;
	int		seg;
	unsigned char	uval;

	/* Get the sign and the magnitude of the value. */
	if (pcm_val < 0) {
		pcm_val = BIAS - pcm_val;
		mask = 0x7F;
	} else {
		pcm_val += BIAS;
		mask = 0xFF;
	}

	/* Convert the scaled magnitude to segment number. */
	seg = search(pcm_val, seg_end, 8);

	/*
	 * Combine the sign, segment, quantization bits;
	 * and complement the code word.
	 */
	if (seg >= 8)		/* out of range, return maximum value. */
		return (0x7F ^ mask);
	else {
		uval = (seg << 4) | ((pcm_val >> (seg + 3)) & 0xF);
		return (uval ^ mask);
	}

}

/*
 * ulaw2linear() - Convert a u-law value to 16-bit linear PCM
 *
 * First, a biased linear code is derived from the code word. An unbiased
 * output can then be obtained by subtracting 33 from the biased code.
 *
 * Note that this function expects to be passed the complement of the
 * original code word. This is in keeping with ISDN conventions.
 */
int
ulaw2linear(
	unsigned char	u_val)
{
	int		t;

	/* Complement to obtain normal u-law value. */
	u_val = ~u_val;

	/*
	 * Extract and bias the quantization bits. Then
	 * shift up by the segment number and subtract out the bias.
	 */
	t = ((u_val & QUANT_MASK) << 3) + BIAS;
	t <<= ((unsigned)u_val & SEG_MASK) >> SEG_SHIFT;

	return ((u_val & SIGN_BIT) ? (BIAS - t) : (t - BIAS));
}

#endif //  AnotherIMP ===============================


/* A-law to u-law conversion */
unsigned char SV_Signal::alaw2ulaw(unsigned char	aval) {

	aval &= 0xff;
	return ((aval & 0x80) ? (0xFF ^ _a2u[aval ^ 0xD5]) :
	    (0x7F ^ _a2u[aval ^ 0x55]));
}

/* u-law to A-law conversion */
unsigned char SV_Signal::ulaw2alaw(unsigned char	uval) {
	uval &= 0xff;
	return ((uval & 0x80) ? (0xD5 ^ (_u2a[0xFF ^ uval] - 1)) :
	    (0x55 ^ (_u2a[0x7F ^ uval] - 1)));
}

