//************************************************************************
//    Calculate LPC residual.
//
//
//    Author  : Jun Zhou
//    Date    : April 15, 2006
//************************************************************************

#include <stdlib.h>
#include <stdio.h>
#include "SC_Feature_LPCresidual.h"
#include "SC_Aux.h"
#include "SC_Signal_WAVE.h"
#include "SC_Conversion.h"
#include "SC_Transform.h"
#include <GN_LPC.h>
#include <SV_Error.h>

//=================================================================
//  Default Constructor
//=================================================================
SC_Feature_LPCresidual::SC_Feature_LPCresidual(int sampleRate, int frameSize, int frameStep, double preemphasis, int order): SV_Feature() {
  this->Para.StpSz = frameStep; //by thilo
	this->Para.WinSz = frameSize; //by thilo
	this->Para.Alpha = preemphasis; //by thilo
	this->Para.LPC_Order = order;	//by thilo: changed the parameter from extra protected variable to this one inside the Para-struct
	this->Para.SRate = sampleRate; //by thilo
}

//=================================================================
//  Default Destructor 
//=================================================================
SC_Feature_LPCresidual::~SC_Feature_LPCresidual() { 

}

//=================================================================
//  extract the residual
//=================================================================
SV_Data *SC_Feature_LPCresidual::ExtractFeature(void) {
	float	*Segment;
  float *Signal;
  long SigLen;
	int frameCount;
	SV_Data *DataSet = NULL;
		
	if (IsSigLoaded()) { 
		Signal = GetSig();
		SigLen = GetLen();
	}	else {
		return (NULL);
	}

	frameCount = (int)(sclib::getRowCount(SigLen, Para.WinSz, Para.StpSz)); //(SigLen / Para.StpSz) - (Para.WinSz / Para.StpSz) + 1;
	if (frameCount > 0) {
		// Preeamphasize if wished
		if (this->Para.Alpha != 0.0) {
			PreEmphasize(this->Para.Alpha);
		}

		DataSet = new SV_Data;
		if (DataSet==NULL) {
			REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
		} 

		DataSet->Row = frameCount;
		DataSet->Col = this->Para.WinSz;
		DataSet->Alloc();
		DataSet->Hdr.frameSize = this->Para.WinSz; //by thilo
		DataSet->Hdr.frameStep = this->Para.StpSz; //by thilo
		DataSet->Hdr.sampleRate = this->Para.SRate; //by thilo
		DataSet->Hdr.ID = sclib::featureLPCresidual; //by thilo
	 
		GN_LPC  LPC_Eng;	// LPC analysis engine
		double  *LpcBuf;
		float  *Rst;

		MArray_1D(LpcBuf, this->Para.LPC_Order, double, "LpcBuf");
		MArray_1D(Segment, this->Para.WinSz, float, "Segment");
		MArray_1D(Rst, this->Para.WinSz,float,"Rst");

		for (int FrmCnt=0; FrmCnt<DataSet->Row; FrmCnt++) {

			//---------------------------------------------------
			// content of Signal copy to Segment
			//---------------------------------------------------
			for (int Col=0; Col< this->Para.WinSz; Col++) {
				Segment[Col] = Signal[FrmCnt*Para.StpSz + Col];
 			}

			//--------------------------------------
			// Calculte LPC coef and residual for one segment
			//--------------------------------------
			LPC_Eng.CalcLPC (Segment, this->Para.WinSz, LpcBuf, this->Para.LPC_Order);
			LPC_Eng.Residual(Segment, Rst, this->Para.WinSz, LpcBuf, this->Para.LPC_Order);

			//copy to Matrix "DataSet"
			for (int Col=0; Col<DataSet->Col; Col++) {
						DataSet->Mat[FrmCnt][Col] = float(Rst[Col]);	
			}
		}	// end of for(FrmCnt)

		MFree_1D(Segment);
		MFree_1D(LpcBuf);
		MFree_1D(Rst);
	}

	return(DataSet);
}

//====================================================================================================================
//	constructs a residual signal (that e.g. can be saved and listened to) from the extracted framed (overlapping) 
//  residual samples
//  by thilo
//====================================================================================================================
SC_Signal* SC_Feature_LPCresidual::getResidualSignal(SV_Data *pResidual) {
	SC_Signal *pResidualSignal = NULL;
	short *pSamples = NULL;
	long int numSamples, start, end, count = 0;
	int halfStep = pResidual->Hdr.frameStep/2; 
	int halfSize = pResidual->Hdr.frameSize/2; 

	numSamples = (pResidual->Row-1)*pResidual->Hdr.frameStep + pResidual->Hdr.frameSize;
	MArray_1D(pSamples, numSamples, short, "SC_Feature_LPCresidual.getResidualSignal: pSamples");

	for (int t = 0; t < pResidual->Row; t++) {
		start = (t > 0) ? halfSize-halfStep : 0; //decide which samples to take from which frame in order to construct the nonoverlapping unframed residual signal
		end = (t < pResidual->Row-1) ? halfSize+halfStep : pResidual->Col;
		for (int s = start; s < end; s++) {
			pSamples[count++] = (short)(sclib::round)(pResidual->Mat[t][s]);
		}
	}	

	pResidualSignal = new SC_Signal_WAVE();
	pResidualSignal->SigPar.SRate = this->Para.SRate;
	pResidualSignal->SigPar.NChannel = 1;
	pResidualSignal->setBuf_L(pSamples, numSamples);

	return pResidualSignal;
}

//====================================================================================================================
//	construct a residual signal (see above; stored in a one column SV_Data container with 1 sample per row) from the 
//  given residual frames of voiced speech that conatins only those samples within a short window around the glottal
//  closure instants; see Dhananjaya, Yegnanarayana, "Speaker Change Detection in Casual Conversations Using 
//  Excitation Source Features", 2008; Naylor, Kounoudes, Gudnason, Brookes, "Estimation of Glottal Closure Instants 
//  in Voiced Speech Using the DYPSA Algorithm", 2007
//  phaseSlopeWindowSize is the size of the window used to compute the average energy weighted phase slope function in 
//  [ms].
//  windowSizeAroundGCI is the size of the window around the found GCIs from which the samples are taken for the final 
//  result (in [ms]).
//  by thilo
//====================================================================================================================
SV_Data* SC_Feature_LPCresidual::getResidualSignalAroundGCIs(SV_Data *pVoicedSpeechResidual, float phaseSlopeWindowSize, float windowSizeAroundGCI) {
	SV_Data *pGCIresidual = NULL;
	double *aewPhaseSlope = NULL; //the averaged energy weighted phase slope function
	long int numSamples, row, col, psWindowSize, gciWindowSize, halfWindowSize, numFinalSamples = 0;
	SC_Conversion converter(pVoicedSpeechResidual->Hdr.sampleRate);
	SC_Transform transformer;
	double *hammingWindow = NULL, energy, windowEnergy, frameMovementRatio = (double)(pVoicedSpeechResidual->Hdr.frameSize)/(double)(pVoicedSpeechResidual->Hdr.frameStep);
	int i, m, start, end, lastZeroCrossing = 0, offset; //, numZeroCrossings = 0;
	//double avgTimeBetweenZeroCrossings = 0.0;
	
	//compute the aew phase slope function for each sample (i.e. for a window starting at each sample with size windowSize)
	//as the center of gravity of the time-domain waveform of the windowed voiced speech residual signal
	psWindowSize = converter.ms2sample_f(phaseSlopeWindowSize);
	numSamples = (pVoicedSpeechResidual->Row-1)*pVoicedSpeechResidual->Hdr.frameStep + pVoicedSpeechResidual->Hdr.frameSize; //number of samples in the residual *signal*
	hammingWindow = transformer.hamming(psWindowSize);
	MArray_1D(aewPhaseSlope, numSamples-psWindowSize, double, "SC_Feature_LPCresidual.getResidualSignalAroundGCIs: aewPhaseSlope");
	for (i = 0; i < numSamples-psWindowSize; i++) {
		aewPhaseSlope[i] = 0.0;
		windowEnergy = 0.0;

		//row = sclib::min(pVoicedSpeechResidual->Row-1, sclib::round(sclib::max(0.0, (double)(i)-((double)(pVoicedSpeechResidual->Hdr.frameSize)/2.0)) / (double)(pVoicedSpeechResidual->Hdr.frameStep)));
		//col = sclib::min(pVoicedSpeechResidual->Col-1, (i) - (pVoicedSpeechResidual->Hdr.frameStep * row));
		//sclib::scalarOut("residualSignal.txt", pVoicedSpeechResidual->Mat[row][col], NULL, true);

		for (m = 0; m < psWindowSize; m++) {
			//compute the coordinates of the (i+m)th sample in the frame-set
			row = sclib::min(pVoicedSpeechResidual->Row-1, sclib::round(sclib::max(0.0, (double)(i+m)-((double)(pVoicedSpeechResidual->Hdr.frameSize)/2.0)) / (double)(pVoicedSpeechResidual->Hdr.frameStep))); //RUNDEN(MAX(0; A5-($B$1/2))/$B$2; 0)
			col = sclib::min(pVoicedSpeechResidual->Col-1, (i + m) - (pVoicedSpeechResidual->Hdr.frameStep * row)); //A5-($B$2*B5)

			//energy = square of signal*window
			energy = hammingWindow[m] * pVoicedSpeechResidual->Mat[row][col];
			energy *= energy;

			aewPhaseSlope[i] += m * energy;
			windowEnergy += energy;
		}
		if (windowEnergy > 0.0) {
			aewPhaseSlope[i] /= windowEnergy;
		} //no else necessary here; if windowEnergy is zero, so is the phase slope value
	}
	MFree_1D(hammingWindow);

	//shift the phase slope so that the window used to compute aewPhaseSlope[i] is centered around i
	offset = (psWindowSize - 1) / 2;
	for (i = numSamples-psWindowSize-2; i >= 0; i--) {
		if (i >= offset) {
			aewPhaseSlope[i] = aewPhaseSlope[i - offset] - offset;
		} else {
			aewPhaseSlope[i] = aewPhaseSlope[i+1]; //avoid introducing a spurious zero crossing in the first subwindow
		}
	}

	//sclib::vectorOut("centeredAewPhaseSlope.txt", aewPhaseSlope, numSamples-psWindowSize, true);

	//count nr of zero crossing in the phase slope and allocate the result set accordingly
	gciWindowSize = converter.ms2sample_f(windowSizeAroundGCI);
	halfWindowSize = gciWindowSize / 2;
	lastZeroCrossing = -halfWindowSize;
	for (i = 0; i < numSamples-psWindowSize-1; i++) {
		if (aewPhaseSlope[i] >= 0.0 && aewPhaseSlope[i+1] < 0.0) { //a zero crossing (form + to -) is found
			start = sclib::max(lastZeroCrossing + halfWindowSize + 1, i - halfWindowSize); //be sure not to count samples more than once
			end = sclib::min(i + halfWindowSize, numSamples - 1);
			numFinalSamples += end - start + 1;
			//numZeroCrossings++;
			//avgTimeBetweenZeroCrossings += 1.0/(double)(numZeroCrssings) * ((i-((i==0)?0:lastZeroCrossing)) - avgTimeBetweenZeroCrossings); //update mean time between GCIs
			lastZeroCrossing = i;
		}
	}
	pGCIresidual = new SV_Data(numFinalSamples, 1);

	//construct the final result by taking the samples in a window with the found GCIs (+ -> - zero crossings in the slope function) in the middle
	numFinalSamples = 0;
	lastZeroCrossing = -halfWindowSize;
	for (i = 0; i < numSamples-psWindowSize-1; i++) {
		if (aewPhaseSlope[i] >= 0.0 && aewPhaseSlope[i+1] < 0.0) { //a zero crossing (form + to -) is found
			start = sclib::max(lastZeroCrossing + halfWindowSize + 1, i - halfWindowSize);
			end = sclib::min(i + halfWindowSize, numSamples - 1);
			
			for (m = start; m <= end; m++) {
				row = sclib::min(pVoicedSpeechResidual->Row-1, sclib::round(sclib::max(0.0, (double)(m)-((double)(pVoicedSpeechResidual->Hdr.frameSize)/2.0)) / (double)(pVoicedSpeechResidual->Hdr.frameStep))); //RUNDEN(MAX(0; A5-($B$1/2))/$B$2; 0)
				col = sclib::min(pVoicedSpeechResidual->Col-1, (m) - (pVoicedSpeechResidual->Hdr.frameStep * row)); //A5-($B$2*B5)
				pGCIresidual->Mat[numFinalSamples++][0] = pVoicedSpeechResidual->Mat[row][col];
			}
			
			lastZeroCrossing = i;
		}
	}
	MFree_1D(aewPhaseSlope);

	return pGCIresidual;
}

//====================================================================================================================
//	converts the result from getResidualSignalAroundGCIs() back to a frame-based representation with 1 frame per row
//  in the result-set; 5ms frameSize should be used
//  by thilo
//====================================================================================================================
SV_Data* SC_Feature_LPCresidual::residualSignal2frames(SV_Data *pResidualSignal, int frameSize, int frameStep, bool doNormalization) {
	SV_Data *pFrames = NULL;
	int x, y, frameCount = sclib::getRowCount(pResidualSignal->Row, frameSize, frameStep);
	double rms;

	if (frameCount > 0 && frameSize > 0) {
		pFrames = new SV_Data(frameCount, frameSize);

		for (y = 0; y < frameCount; y++) {
			rms = 0.0;
			for (x = 0; x < frameSize; x++) {
				pFrames->Mat[y][x] = pResidualSignal->Mat[frameStep*y + x][0];
				if (doNormalization == true) {
					rms += pFrames->Mat[y][x]*pFrames->Mat[y][x];
				}
			}
			if (doNormalization == true) {
				rms = sqrt(rms / (double)(frameSize)); //the Root Mean Square (RMS) value of the samples in the frame
				for (x = 0; x < frameSize; x++) { //normalize to unit norm; see Dhananjaya & Yegnanarayana's paper
					pFrames->Mat[y][x] = (float)(pFrames->Mat[y][x] / rms);
				}
			}
		}
	}	

	return pFrames;
}

//====================================================================================================================
//	converts the result from getResidualSignal() back to a frame-based representation with 1 frame per row in the 
//  result-set; 5ms frameSize should be used
//  by thilo
//====================================================================================================================
SV_Data* SC_Feature_LPCresidual::residualSignal2frames(SC_Signal *pResidualSignal, int frameSize, int frameStep, bool doNormalization) {
	SV_Data *pFrames = NULL;
	int x, y, frameCount = sclib::getRowCount(pResidualSignal->GetLen(), frameSize, frameStep);
	short *signal = pResidualSignal->GetBuf_L();
	double rms;

	if (frameCount > 0 && frameSize > 0) {
		pFrames = new SV_Data(frameCount, frameSize);

		for (y = 0; y < frameCount; y++) {
			rms = 0.0;
			for (x = 0; x < frameSize; x++) {
				pFrames->Mat[y][x] = signal[frameStep*y + x];
				if (doNormalization == true) {
					rms += pFrames->Mat[y][x]*pFrames->Mat[y][x];
				}
			}
			if (doNormalization == true) {
				rms = sqrt(rms / (double)(frameSize)); //the Root Mean Square (RMS) value of the samples in the frame
				for (x = 0; x < frameSize; x++) { //normalize to unit norm; see Dhananjaya & Yegnanarayana's paper
					pFrames->Mat[y][x] = (float)(pFrames->Mat[y][x] / rms);
				}
			}
		}
	}	

	return pFrames;
}

//====================================================================================================================
//	(re)-set the class members
//  by thilo
//====================================================================================================================
void SC_Feature_LPCresidual::setParameters(int sampleRate, int frameSize, int frameStep, double preemphasis, int order) {
  this->Para.StpSz = frameStep;
	this->Para.WinSz = frameSize;
	this->Para.Alpha = preemphasis;
	this->Para.LPC_Order = order;
	this->Para.SRate = sampleRate;

	return;
}
