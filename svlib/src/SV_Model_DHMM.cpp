//************************************************************************
//    Implementation of DHMM
//
//
//    Author  : Jialong HE
//    Date    : May 13, 1999
//    
//    TODO : too short sequence (<StaNum) may have problem
//           When test has longer string than training?
//   
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "SV_Model_DHMM.h"
#include "SV_Error.h"
#include "GN_Matrix.h"
#include "GN_Rand.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//==========================================
// default constructor
//==========================================
SV_Model_DHMM::SV_Model_DHMM() {

	States = NULL; 
	StaNum = 0;

	NState = 3;       // default state number, used when AllocModel()
	SymNum = 0;       // codebook size for easy access 

	//-------------------------------------
	// Temp varible for training
	//-------------------------------------
	num_strings = 0;
	strings = NULL;
	string_len = NULL;

	alpha = NULL;                     // matrix for alpha trellis
	beta  = NULL;                    // matrix for beta trellis
	gamma_next   = NULL;               // matrix of gamma transition probs
	gamma_recur  = NULL;             // matrix of gamma recurrence probs
	scaling_factors  = NULL;         // array of alpha and beta scaling factors
	a_numer_sum_recur = NULL;        // array of numerators for a_ij's
	a_numer_sum_next  = NULL;         // array of numerators for a_ij's
	a_denom_sum_recur = NULL;        // array of denomonators for a_ij's
	a_denom_sum_next  = NULL;         // array of denomonators for a_ij's
	b_numer_sum_recur = NULL;       // array of numerators for b_ij's
	b_numer_sum_next  = NULL;        // array of numerators for b_ij's

}

//==========================================
// default destructor
//==========================================
SV_Model_DHMM::~SV_Model_DHMM() {
	
	FreeModel();
}

//==========================================
// Allocate memory for DHMM
//==========================================
void SV_Model_DHMM::AllocModel(void) {

	SymNum = GetCNum();
	StaNum = NState;        // NState is set from outside

	if (SymNum <=0 || StaNum <= 0) {
		REPORT_ERROR(SVLIB_BadArg, "GMM:: AllocModel");
	}

	
	//------------------------------------------
	// In order to call constructor of state
	// use column vector (one col matrix) rather
	// than an array of state.
	//------------------------------------------
	States = new state* [StaNum];
	for (int i=0; i < StaNum; i++) {
		States[i] = new state(SymNum);
	}

}

//==========================================
// Free memory for DHMM
//==========================================
void SV_Model_DHMM::FreeModel(void) {
	
	if (States != NULL) {
		for (int i=0; i < StaNum; i++) {
			delete States[i];
		}
		delete  States;
	}

	States = NULL;
}

//===============================================
// Test DHMM
//===============================================
SV_Data* SV_Model_DHMM::TestModel(SV_Data *pData) {

	SV_Data *Score, *DataCurr;
	int PatNum, PatCnt;

	if (pData->Col != Dim) {
		REPORT_ERROR(SVLIB_BadData, "Model and Data have different Dim!");
	}

	PatNum = CountPat(pData);
	//----------------------------
	// allocate score space
	//----------------------------
	Score = new SV_Data;
	Score->Row = PatNum;	
	Score->Col = 1;	         // one score for each pattern
	Score->Alloc();
	
	//----------------------------------------------
	// Convert vector sequences to discret symbols
	//----------------------------------------------
	trellis_width = GenSymbol(pData) + 1;
	MArray_2D(alpha, trellis_width, StaNum, double, "alpha");
	MArray_1D(scaling_factors, trellis_width, double, "scaling_factors");

	//----------------------------
	// test score for each pattern
	//----------------------------
	DataCurr = pData;
	for (PatCnt=0; PatCnt<PatNum; PatCnt++) {

		Score->Mat[PatCnt][0] = (float)test(strings[PatCnt], string_len[PatCnt]);

		DataCurr = DataCurr->Next;		   // point to next
	}

	MFree_2D(alpha);
	MFree_1D(scaling_factors);
	return (Score);
}

//================================================================
// Save the model to current opened model file
// if success, return total bytes writed, otherwise, REPORT_ERROR
// by thilo: chaged all calls to write()
//================================================================
int SV_Model_DHMM::SaveModel(void) {

	int RtCode, StaCnt, SymCnt;
	//int TotalByte; //by thilo
	double *Buf;
	int bytes; //by thilo
	SV_DataIO io; //by thilo


	Hdr.ModelType = MT_DHMM;
    RtCode = SaveHdr();
	if (RtCode == SVLIB_Fail) {
		REPORT_ERROR(SVLIB_Fail, "Save DHMM model Failed!");
	}
	
	//--------------------------
	// VQ codebook parameters
	//--------------------------
	//DFile.write((char*)(&Dim), sizeof(int));
	bytes = io.writeScalar(&(this->DFile), this->Dim);
	//DFile.write((char*)(&CNum), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->CNum);
	//DFile.write((char*)(CBook[0]), CNum*Dim*sizeof(float));
	bytes += io.writeMatrix(&(this->DFile), this->CBook, this->CNum, this->Dim);
	if (DFile.good() != TRUE) {
		REPORT_ERROR(SVLIB_Fail, "Save VQ Codebook Failed!");
	}

	//TotalByte = MHLen + 2*sizeof(int) +  CNum*Dim*sizeof(float);

	//DFile.write((char*)(&StaNum), sizeof(int));
	bytes += io.writeScalar(&(this->DFile), this->StaNum);

	MArray_1D(Buf, CNum+1, double, "Buf");
	//----------------------------------------------
	// Save State Prob
	//----------------------------------------------
	for (StaCnt=0; StaCnt < StaNum; StaCnt++) {
		//------------------------------------------
		// Copy recur prob into Buf and save Buf
		//------------------------------------------
		Buf[0] = States[StaCnt]->recur_trans;
		for (SymCnt=0;  SymCnt<States[StaCnt]->SymNum; SymCnt++) {
			Buf[SymCnt+1] =	States[StaCnt]->recur_out[SymCnt];
		}

		//DFile.write((char*)Buf, (CNum+1)*sizeof(double));
		bytes += io.writeArray(&(this->DFile), Buf, this->CNum+1);
		if (DFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Save recur prob. of DHMM Failed!");
		}

		//------------------------------------------
		// Copy next prob into Buf and save Buf
		//------------------------------------------
		Buf[0] = States[StaCnt]->next_trans;
		for (SymCnt=0;  SymCnt<States[StaCnt]->SymNum; SymCnt++) {
			Buf[SymCnt+1] =	States[StaCnt]->next_out[SymCnt];
		}

		//DFile.write((char*)Buf, (CNum+1)*sizeof(double));
		bytes += io.writeArray(&(this->DFile), Buf, this->CNum+1);
		if (DFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Save next prob. of DHMM Failed!");
		}

	}

	//TotalByte += sizeof(int) + StaNum * (1 + CNum) * 2 * sizeof(double);
	MFree_1D(Buf);
	//return(TotalByte);
	return bytes + MHLen; //MHLen is maybe not what was really written by SaveHdr()...
}

//===========================================================
// Load model's parameter from current opened model file
// if success, return (this) pointer, if fail, return (NULL)
//===========================================================
SV_Model * SV_Model_DHMM::LoadModel(void) {

	int RtCode, NewDim, NewSize, StaCnt, SymCnt;
	double *Buf;
	int bytes; //by thilo
	SV_DataIO io; //by thilo
	SV_DataIO::SV_DatatypeSizes codeSizes, fileSizes; //by thilo
	io.getCurrentDatatypeSizes(codeSizes); //by thilo

    RtCode = LoadHdr(fileSizes); //by thilo: changed to catch the machine dependant header
	if (RtCode == SVLIB_Fail) {
		return(NULL);
	}

	//--------------------------
	// check if valid header
	//--------------------------
	if (Hdr.ModelType != MT_DHMM) {
		return(NULL);
	}

	//--------------------------
	// VQ model's parameters
	//--------------------------
	//DFile.read((char*)(&NewDim), sizeof(int));
	bytes = io.readScalar(&(this->DFile), NewDim, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewDim == 0) {
		return(NULL);
	}

	
	//DFile.read((char*)(&NewSize), sizeof(int));
	bytes += io.readScalar(&(this->DFile), NewSize, codeSizes, fileSizes);
	if (DFile.good() != TRUE || NewSize == 0) {
		return(NULL);
	}

	//--------------------------------------------
	// Allocate CBook
	//--------------------------------------------
	if (NewDim != Dim || NewSize != CNum) {
		MFree_2D(CBook);
		CNum = NewSize;
		Dim  = NewDim;
		MArray_2D(CBook, CNum, Dim, float, "LoadVQ");

	}; 
	
	//DFile.read((char*)(CBook[0]), CNum*Dim*sizeof(float));
	bytes += io.readMatrix(&(this->DFile), this->CBook, this->CNum, this->Dim, codeSizes, fileSizes);
	
	if (DFile.good() != TRUE) {
		return(NULL);
	}

	//-------------------------------------
	// Load States
	//-------------------------------------
	//DFile.read((char*)(&NState), sizeof(int));
	bytes += io.readScalar(&(this->DFile), this->NState, codeSizes, fileSizes);
	FreeModel();	// free old memory of States
	AllocModel();	// use NState allocate new space for States
	MArray_1D(Buf, CNum+1, double, "Buf");

	//----------------------------------------------
	// Load State Prob
	//----------------------------------------------
	for (StaCnt=0; StaCnt < StaNum; StaCnt++) {
		//DFile.read(( char*)Buf, (CNum+1)*sizeof(double));
		bytes += io.readArray(&(this->DFile), Buf, this->CNum+1, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Load recur prob. of DHMM Failed!");
		}

		//------------------------------------------
		// Copy recur prob into Buf and save Buf
		//------------------------------------------
		States[StaCnt]->recur_trans = Buf[0];
		for (SymCnt=0;  SymCnt<States[StaCnt]->SymNum; SymCnt++) {
			States[StaCnt]->recur_out[SymCnt] = Buf[SymCnt+1];
		}


		//DFile.read((char*)Buf, (CNum+1)*sizeof(double));
		bytes += io.readArray(&(this->DFile), Buf, this->CNum+1, codeSizes, fileSizes);
		if (DFile.good() != TRUE) {
			REPORT_ERROR(SVLIB_Fail, "Load next prob. of DHMM Failed!");
		}

		//------------------------------------------
		// Copy next prob into Buf and save Buf
		//------------------------------------------
		States[StaCnt]->next_trans = Buf[0];
		for (SymCnt=0;  SymCnt<States[StaCnt]->SymNum; SymCnt++) {
			States[StaCnt]->next_out[SymCnt] = Buf[SymCnt+1];
		}

	}

	MFree_1D(Buf);
	return(this);
}

//=============================================
// Dump model's parameter in ASCII 
//=============================================
ostream& operator<< (ostream& OutS, SV_Model_DHMM& Data) {

	int Row, Col;
	int CNum, Dim;
	float **CBook;

	//--------------------------------
	// Output Codebook
	//--------------------------------
	CNum  = Data.GetCNum();
	Dim   = Data.GetCDim();
	CBook = Data.GetCBook();

	OutS << CNum <<" "<<Dim << endl;
	for (Row=0; Row<CNum; Row++) {
		for (Col=0; Col<Dim; Col++)
			OutS<<CBook[Row][Col]<< " ";
		OutS<<endl;
	}

	//------------------------------------
	// Output transition and output prob.
	//------------------------------------
	OutS << Data.StaNum << endl;
	for (Row=0; Row<Data.StaNum; Row++) {

		OutS <<Data.States[Row]->recur_trans << " " <<Data.States[Row]->next_trans << endl;  // tran. prob.

		for (Col=0; Col<CNum; Col++) {
			OutS<<Data.States[Row]->recur_out[Col]<< " ";
		}
		OutS<<endl;

		for (Col=0; Col<CNum; Col++) {
			OutS<<Data.States[Row]->next_out[Col]<< " ";
		}
		OutS<<endl;
	}


	
	return(OutS);
}




//============================================================
// Init prob. to random 
// Suppose "States" have been allocted by AllocModel()
//============================================================
void SV_Model_DHMM::InitModel(void) {

	GN_Rand REng;
	int SymCnt, StaCnt;
	double Sum, SumRecur, SumNext;

	REng.srandom(RandSeed);
	
	//----------------------------------------
	// For each state, Init Tran and OutProb 
	//----------------------------------------
	for (StaCnt=0; StaCnt< StaNum; StaCnt++) {
		//----------------------------
		// Init recur_out and next_out
		//----------------------------
		SumRecur  = 0.0;
		SumNext   = 0.0; 
		for (SymCnt=0; SymCnt<States[StaCnt]->SymNum; SymCnt++) {
			States[StaCnt]->recur_out[SymCnt] = fabs(double(REng.random()));
			SumRecur += States[StaCnt]->recur_out[SymCnt]; 

			States[StaCnt]->next_out[SymCnt] = fabs(double(REng.random()));
			SumNext += States[StaCnt]->next_out[SymCnt]; 
		}

		//---------------------------------
		// Normalize recur_out and next_out
		//---------------------------------
		for (SymCnt=0; SymCnt<States[StaCnt]->SymNum; SymCnt++) {
			States[StaCnt]->recur_out[SymCnt] /= SumRecur;
			States[StaCnt]->next_out[SymCnt] /= SumNext;
		}

		//---------------------------------
		// Init and normalize Tran Prob
		//---------------------------------
		States[StaCnt]->recur_trans = fabs(double(REng.random()));
		States[StaCnt]->next_trans = fabs(double(REng.random()));

		Sum = States[StaCnt]->recur_trans + States[StaCnt]->next_trans;
		States[StaCnt]->recur_trans /= Sum;
		States[StaCnt]->next_trans /= Sum;

	}  // for (StaCnt)

	//---------------------------------------
	// last state, only self loop
	//---------------------------------------
    States[StaNum-1]->recur_trans = 1; 
    States[StaNum-1]->next_trans  = 0.0;

};

//==============================================
// Convert vector sequences to discrete symbols
// using the current CBook.
//
//	short **strings;
//	int *string_len;
//	int num_strings;
//
//  return max length of strings
//==============================================
int SV_Model_DHMM::GenSymbol(SV_Data *pData) {

	SV_Data *DataCurr;

	float **CodeBook;
	int Row, Dim, CSize;
	int MaxLen, StrNum;
	
	//-------------------------------
	// check input 
	//-------------------------------
	if (pData == NULL) {
		REPORT_ERROR(SVLIB_BadArg, "GenSymbol");
	} 

	//-------------------------------
	// check codebook
	//-------------------------------
	CodeBook = GetCBook();
	CSize	 = GetCNum();   
	Dim      = GetCDim();
	if ( CodeBook == NULL || CSize <= 0) {
		REPORT_ERROR(SVLIB_BadArg, "GenSymbol:: no codebook");
	}

	//------------------------------------------------------
	// figure out number of patterns and max string length
	//------------------------------------------------------
	DataCurr = pData;
	MaxLen = -1;
	StrNum = 0;
	while (DataCurr != NULL) {
		if (Dim != DataCurr->Col) {
			REPORT_ERROR(SVLIB_BadData, "Different Dimension"); 	
		}
		
		StrNum++;
		if (DataCurr->Row > MaxLen) {
			MaxLen = DataCurr->Row;
		}
		DataCurr = DataCurr->Next;
	}  // while

	//------------------------------------------------------
	// allocate memory for strings, string_len;
	//------------------------------------------------------
    MFree_2D(strings);
	MFree_1D(string_len);
	MArray_2D(strings, StrNum, MaxLen, short, "GenSymbol::strings");
	MArray_1D(string_len, StrNum, int, "string_len");
	num_strings = StrNum;

	//------------------------------------------------------
	// find nearest code for each vector
	//------------------------------------------------------
	DataCurr = pData;
	StrNum   = 0;
	while (DataCurr != NULL) {
		for (Row=0; Row<DataCurr->Row; Row++) {		
			strings[StrNum][Row] = nearest_code (CodeBook, DataCurr->Mat[Row], CSize, Dim);
		}
		string_len[StrNum] = DataCurr->Row;

		StrNum++;
		DataCurr = DataCurr->Next;
	}  // while

	return(MaxLen);
}

//==========================================
// train DHMM model
//==========================================
int SV_Model_DHMM::TrainModel(SV_Data *pData) { //by thilo: changed return-value from void to int (needed in derived class to indicate error)
	
	
//	SV_Data *pTData;

//	int Row, Col, MixCnt, DimCnt, DataDim, CurrMix;

	//------------------------------------------
	// Generate Codebook by calling VQ's method
	//------------------------------------------
	SV_Model_VQ::TrainModel(pData);

	//------------------------------------------
	// Convert vector sequence to discrete symbol
	// return longest string length
	//------------------------------------------
	trellis_width = GenSymbol(pData) + 1;

	//------------------------------------------
	// allocate emission prob memory and init 
	//------------------------------------------
	AllocModel();			
	InitModel();			
	
	//------------------------------------------
	// Generate model 
	//------------------------------------------
    BatchTrain();

    return 0;
}

//========================================================//
// train hmm  until the sum of the probability            //
// changes is less than min_delta_psum.                   //
//========================================================//
double SV_Model_DHMM::BatchTrain(void) {

	double ThisProb = -LARGE / 2.0, LastProb = -LARGE;
	int Iter=0, ThresNum = 0;

	SymNum = GetCNum();
	//-------------------------------------
	// Allocate temp memory for training
	//-------------------------------------
	MArray_2D(alpha, trellis_width, StaNum, double, "alpha");
	MArray_2D(beta, trellis_width, StaNum, double, "beta");
	MArray_2D(gamma_next, trellis_width, StaNum, double, "gamma_next");
	MArray_2D(gamma_recur, trellis_width, StaNum, double, "gamma_recur");

	MArray_1D(scaling_factors, trellis_width, double, "scaling_factors");

	MArray_1D(a_numer_sum_next, StaNum, double, "a_numer_sum_next");
	MArray_1D(a_numer_sum_recur, StaNum, double, "a_numer_sum_recur");
	MArray_1D(a_denom_sum_next, StaNum, double, "a_denom_sum_next");
	MArray_1D(a_denom_sum_recur, StaNum, double, "a_denom_sum_recur");

	MArray_2D(b_numer_sum_next, StaNum, SymNum, double, "b_number_sum_next");
	MArray_2D(b_numer_sum_recur, StaNum, SymNum, double, "b_number_sum_recur");

	while ( ThresNum <= 3 && Iter <= MaxIter ) {

		/*---------------------------------------------------------*/
		/* count the iteration number coutinuously lower the thres */
		/*---------------------------------------------------------*/
		if ((ThisProb - LastProb) / fabs(ThisProb) < 0.001) {
			ThresNum++;
		}
		else {ThresNum = 0;}

		LastProb = ThisProb;
		ThisProb = set_cumulative_ab_counts();
		reestimate();

		Iter++;
		if (Verbose) {
			cerr << ThisProb << endl;
		}
	} 

	//-------------------------------------
	// Free temp memory for training
	//-------------------------------------
//	MFree_2D(alpha);
	MFree_2D(beta);
	MFree_2D(gamma_next);
	MFree_2D(gamma_recur);

	MFree_1D(scaling_factors);

	MFree_1D(a_numer_sum_next);
	MFree_1D(a_numer_sum_recur);
	MFree_1D(a_denom_sum_next);
	MFree_1D(a_denom_sum_recur);

	MFree_2D(b_numer_sum_next);
	MFree_2D(b_numer_sum_recur);


    return (ThisProb);

}  

//---------------------------------------------------
// reestimate parameters after calculating ab_counts
//---------------------------------------------------
double SV_Model_DHMM::reestimate(void) {

	double prob, diff_sum = 0.0;
	int k, i, j;
	
	// loop all the transition probs
	for ( i=0, j=1; j < StaNum; i++, j++) {
		//----------------------------------------------
		// do all the prob transitions except last one
		//----------------------------------------------
		prob = a_numer_sum_recur[i] / a_denom_sum_recur[i];
		diff_sum += fabs(prob - States[i]->set_recur_trans());
		States[i]->set_recur_trans(prob);

		for (k=0; k < SymNum; k++) {
			prob = b_numer_sum_recur[i][k] / a_numer_sum_recur[i];
			diff_sum += fabs(prob -	States[i]->set_recur_out(k));
			States[i]->set_recur_out(k,prob);
		}
		
		//----------------------------------
		// do all the next transitions
		//----------------------------------
		prob = a_numer_sum_next[i] / a_denom_sum_next[i];
		diff_sum += fabs(prob - States[i]->set_next_trans());
		States[i]->set_next_trans(prob);
		
		for (k=0; k < SymNum; k++) {
			prob = b_numer_sum_next[i][k] / a_numer_sum_next[i];
			diff_sum += fabs(prob - States[i]->set_next_out(k));
			States[i]->set_next_out(k,prob);
		}
	}

	//----------------------------------------------
	// calc the recurrent prob for last state
	//----------------------------------------------
	i = StaNum-1;
	prob = a_numer_sum_recur[i] / a_denom_sum_recur[i];
	diff_sum += fabs(prob - States[i]->set_recur_trans());
	States[i]->set_recur_trans(prob);
	
	for (k=0; k < SymNum; k++) {
		prob = b_numer_sum_recur[i][k] / a_numer_sum_recur[i];
		diff_sum += fabs(prob - States[i]->set_recur_out(k));
		States[i]->set_recur_out(k,prob);
	}
	
	return diff_sum;
}

//---------------------------------------------------------------
// set cumulative ab counts from an array of equal length strings.
//---------------------------------------------------------------
double SV_Model_DHMM::set_cumulative_ab_counts(void) {

	double alpha_tot = 0.0, log_alpha = 0.0;
	int i, j, k, s;

	//--------------------------------
	// clear cumulative sum arrays
	//--------------------------------
	for(i=0; i < StaNum; i++) {
		a_numer_sum_recur[i] = 0.0;
		a_denom_sum_recur[i] = 0.0;
		a_numer_sum_next[i] = 0.0;
		a_denom_sum_next[i] = 0.0;
		for (int j=0; j < SymNum; j++) {
			b_numer_sum_recur[i][j] = 0.0;
			b_numer_sum_next[i][j] = 0.0;
		}
	}
	
	//-----------------------------------------
	// loop all the strings calc a,b sums
	//-----------------------------------------
	for(s=0; s < num_strings; s++) {

		// fill the alpha matrix
		alpha_tot += log(alpha_F(strings[s], string_len[s]));
		
		// calc log probability from scaling factors
		if (scaling_factors[0]>0.0) { // if scaling_factors set
			for (int t=0; t < string_len[s]+1; t++) {
				log_alpha += log(scaling_factors[t]);
			}
		}

		// fill the beta matrix
		beta_I(strings[s], string_len[s]);

		// fill the gamma_next and gamma_recur matrices
		compute_gamma(strings[s], string_len[s]);

		for (i=0,j=1; j < StaNum; i++,j++) {
			a_numer_sum_recur[i] += a_numer(i,i, string_len[s]);
			a_numer_sum_next[i] += a_numer(i,j, string_len[s]);
			a_denom_sum_recur[i] += a_denom(i,i, string_len[s]);
			a_denom_sum_next[i] += a_denom(i,j, string_len[s]);
			for (k=0; k < SymNum; k++) {
				b_numer_sum_recur[i][k] += b_numer(i,i,k,strings[s], string_len[s]);
				b_numer_sum_next[i][k] += b_numer(i,j,k,strings[s], string_len[s]);
			}
		}

		// compute recurrent probs for last state
		i = StaNum-1;
		a_numer_sum_recur[i] += a_numer(i,i, string_len[s]);
		a_denom_sum_recur[i] += a_denom(i,i, string_len[s]);
		for (k=0; k < SymNum; k++) {
			b_numer_sum_recur[i][k] += b_numer(i,i,k,strings[s], string_len[s]);
		}
	}

	return(log_alpha/num_strings);
}


//-------------------------------------------
// rescale alpha values (from Rabiner).
//------------------------------------------- 
void SV_Model_DHMM::rescale_alphas(int col) {

	int i;
	scaling_factors[col] = 0.0;
	for (i=0; i < StaNum; i++) {
		scaling_factors[col] += alpha[col][i];
	}
	
	// rescale all the alpha's and smooth them
	for (i=0; i <  StaNum; i++) {
		alpha[col][i] /= scaling_factors[col];
	}
}

//----------------------------------------------
// rescale beta values after rescaling alphas.
//----------------------------------------------
void SV_Model_DHMM::rescale_betas(int col) {

	int i;	
	// rescale all the beta's w alpha's factors
	for (i=0; i < StaNum; i++) {
		beta[col][i] /= scaling_factors[col];
	}
}

//-------------------------------------------------------
// calculate alpha trellis and return final alpha value 
// (recognition prob).
//------------------------------------------------------- 
double SV_Model_DHMM::alpha_F(short* symbol_array, int symbol_count ) {

  double accum;
  int i,j;

  symbol_count=(symbol_count<0?trellis_width-1:symbol_count);

  // clear the trellis and set the first column
  for (i=0; i < trellis_width; i++)
    for (j=0; j < StaNum; j++) {
      alpha[i][j] = 0.0;
    }
  alpha[0][0] = 1.0;  // first state starts with alpha = 1.0
  rescale_alphas(0);  // first col rescale is div by one

  // calculate the alphas for rest of cols
  for (i = 0; i < symbol_count; i++) { 		
    for (j=StaNum-1; j > 0; j--) {
      accum = alpha[i][j] * 
	States[j]->set_recur_trans() *
	  States[j]->set_recur_out(symbol_array[i]);
      alpha[i+1][j] = alpha[i][j-1] * 
	States[j-1]->set_next_trans() * 
	  States[j-1]->set_next_out(symbol_array[i]) + 
	    accum;
    }
    alpha[i+1][0] = alpha[i][0] * States[0]->set_recur_trans() *
      States[0]->set_recur_out(symbol_array[i]);
    rescale_alphas(i+1);
  }

  return (alpha[symbol_count][SymNum-1]);

}


//----------------------------------------------------------
// calculate beta trellis and return initial beta value,
// instead of going bottom to top, go top to bottom; the bottom
// is the exception this time.  we look at transistions from current
// state (single) to the next states (plural).
//----------------------------------------------------------
double SV_Model_DHMM::beta_I(short* symbol_array, int symbol_count) {

     double accum;
     int i, j, t;

     symbol_count=(symbol_count<0?trellis_width-1:symbol_count);     

     // clear the trellis and set the last column
     for ( i=0; i < trellis_width; i++)
       for (j=0; j < StaNum; j++) {
	 beta[i][j] = 0.0;
       }
     beta[symbol_count][StaNum-1] = 1.0;  /* final beta is 1.0 */
     rescale_betas(symbol_count); // this is really div by 1 so no effect

     // begin setting all the cols except last one
     for ( t = symbol_count-1; t >= 0; t--) {
       for ( j=0; j < StaNum-1; j++) {
	 accum = beta[t+1][j] * 
	   States[j]->set_recur_trans() *
	     States[j]->set_recur_out(symbol_array[t]);
	 beta[t][j] = beta[t+1][j+1] * 
	   States[j]->set_next_trans() * 
	     States[j]->set_next_out(symbol_array[t]) + 
	       accum;
       }

       beta[t][StaNum-1] = beta[t+1][StaNum-1] *
	 States[StaNum-1]->set_recur_trans() *
	   States[StaNum-1]->set_recur_out(symbol_array[t]);
       rescale_betas(t);
     }

     return (beta[0][0]);
   }

//---------------------------------------------------------------
// compute gamma values and fill gamma tables for entire trellis.
//---------------------------------------------------------------
void SV_Model_DHMM::compute_gamma(short* symbol_array, int symbol_count ) {

  int i, j, t;
  symbol_count=(symbol_count<0?trellis_width-1:symbol_count);

  // clear gammas
  for (i=0; i < symbol_count; i++)
    for ( j=0; j < StaNum; j++) {
	gamma_recur[i][j] = 0.0;
	gamma_next[i][j] = 0.0;
      }
  // calc the gamma table
  for ( t = 0; t < symbol_count; t++) {
    for ( i=0,j=1; i < StaNum-1; i++,j++) {
      gamma_next[t][i] = (alpha[t][i]*
			  States[i]->set_next_trans() *
			  States[i]->set_next_out(symbol_array[t]) *
			  beta[t+1][j]) /
			    alpha[symbol_count][StaNum-1];
      gamma_recur[t][i] = (alpha[t][i]* 
			   States[i]->set_recur_trans() * 
			   States[i]->set_recur_out(symbol_array[t]) *
			   beta[t+1][i]) /
			     alpha[symbol_count][StaNum-1];  
    }
    gamma_recur[t][StaNum-1] = (alpha[t][StaNum-1] * 
				    States[StaNum-1]->set_recur_trans() * 
		       States[StaNum-1]->set_recur_out(symbol_array[t]) *
			 	    beta[t+1][StaNum-1]) /
				      alpha[symbol_count][StaNum-1];
    gamma_next[t][StaNum-1] = 0.0;
  }
}

//-----------------------------------------------------------
// this is the numerator of the a_ij reestimate and denom of b_ij reestimate,
// this func assumes gammas have been calc'd.
//------------------------------------------------------------
double SV_Model_DHMM::a_numer(int i, int j, int symbol_count) {

  double sum = 0.0;
  int t;

  symbol_count=(symbol_count<0?trellis_width-1:symbol_count);

  for(t=0; t < symbol_count; t++)
    {
      if (i==j) {
	sum += gamma_recur[t][i];
      }
      else if ( i < StaNum-1 ) {
	sum += gamma_next[t][i];
      } else {
	cerr << "WARNING: gamma_next[t]["<<i<<"] shouldn't be requested\n";
	break;
      }
    }
  return sum;
}

//-------------------------------------------------
// this is the denominator of the a_ij reestimate.
//-------------------------------------------------
double SV_Model_DHMM::a_denom(int i, int j, int symbol_count ) {

  int t;	
  symbol_count=(symbol_count<0?trellis_width-1:symbol_count);

  j = j;  // j not used in this func - pass for consistency

  double sum = 0.0;
  for( t=0; t < symbol_count; t++)
    if (i<StaNum-1)
      sum += gamma_recur[t][i] + gamma_next[t][i];
    else
      sum += gamma_recur[t][i];

  return sum;
}

//-----------------------------------------------------
// this is the numerator of the b_ij reestimate.
//----------------------------------------------------- 
double SV_Model_DHMM::b_numer(int i, int j, int sym, short *symbol_array, int symbol_count ) {

  int t;
  symbol_count=(symbol_count<0?trellis_width-1:symbol_count);

  double sum = 0.0;
  for( t=0; t < symbol_count; t++) {
    if ( symbol_array[t] == sym ) 
      {
	if (i==j)
	  sum += gamma_recur[t][i];
	else  {
	  if (i < StaNum-1)
	    sum += gamma_next[t][i];
	  else {
	    cerr << "WARNING: gamma_next[t]["<<i<<"] shouldn't be requested\n";
	    break;
	  }
	}
      }
  }
  return sum;
}

//===================================
// test model for a single string.
//===================================
double SV_Model_DHMM::test(short* string, int symbol_count ) {

  double log_alpha;

  // fill alpha trellis and scaling coeffs
  log_alpha=log(alpha_F(string,symbol_count));

  // calc log probability from coeffs
  if (scaling_factors[0]>0.0) { // if scaling_factors set
    log_alpha=0.0;
    for (int t=0; t < symbol_count+1; t++)
      log_alpha += log(scaling_factors[t]);
  }
  return log_alpha;
}

//****************** Implementation for State Class *****************
//------------------------------------
// Constructor for State class
//------------------------------------
state::state(int num_symbols) {

	MArray_1D(recur_out, num_symbols, double, "recur_out");
	MArray_1D(next_out, num_symbols, double, "next_out");
	SymNum = num_symbols;
}

//------------------------------------
// Destructor for State class
//------------------------------------
state::~state() {

	MFree_1D(recur_out);
	MFree_1D(next_out);
}

//--------------------------------------------------------------
// set/return the probability of generating a particluar symbol 
// during a recurrent transition. (smooth with MIN_PROB)
//--------------------------------------------------------------
double state::set_recur_out(int symbol, double prob) {

  if (prob != -1.0) {
    recur_out[symbol] = prob;
  }
  if (recur_out[symbol] < MIN_PROB)
    recur_out[symbol] = MIN_PROB;
  return recur_out[symbol];
}

//---------------------------------------------------------------
// set/return the probability of generating a particluar symbol 
// during a transition to the next state. (smooth with MIN_PROB)
//---------------------------------------------------------------
double state::set_next_out(int symbol, double prob) {

  if (prob != -1.0) {
    next_out[symbol] = prob;
  }
  if (next_out[symbol] < MIN_PROB)
    next_out[symbol] = MIN_PROB;
  return next_out[symbol];
}

//------------------------------------------------------------
// set/return the probability of making a recurrent transition.
// (smooth with MIN_PROB)
//------------------------------------------------------------
double state::set_recur_trans(double prob) {

  if (prob != -1.0) {
    recur_trans = prob;
  }
  if (recur_trans < MIN_PROB)
    recur_trans = MIN_PROB;
  return recur_trans;
}

//------------------------------------------------------------------
// set/return the probability of making transitioning to the next state.
// (smooth with MIN_PROB)
//------------------------------------------------------------------
double state::set_next_trans(double prob) {

  if (prob != -1.0) {
    next_trans = prob;
  }
  if (next_trans < MIN_PROB)
    next_trans = MIN_PROB;
  return next_trans;
}


