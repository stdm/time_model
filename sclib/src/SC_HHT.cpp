//simple method for Hilbert-Transform
//Author: Jun,Zhou
//Date:   March, 2007

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <vld.h> // Include Visual Leak Detector

#include "SC_HHT.h"
#include "SC_Aux.h"
#include "SC_TweakableParameters.h"
#include <SV_General.h>
#include <SV_Data.h>
#include <SV_Error.h>

//===========================================================
//  Constructor
//===========================================================
SC_HHT::SC_HHT(SC_TweakableParameters *pTweak) {
	this->pTweak=pTweak;
}


//===========================================================
//  Destructor
//===========================================================
SC_HHT::~SC_HHT() {

}

SV_Data* SC_HHT::emd(short *Signal, long int Len) 
{	
	int NBSYM = 2;
//	int nbsym;
	double length_d = (double) Len;
	long int nh_r,nl_r,nh_h,nl_h,n_zero,nh_pre,nl_pre;
	double *mean, *r,*h,*a_par,*delta;
	int max_imf = 51;

	long int n1;//number of samples,which are smaller than first thredhold
	double g1 = 0.05;//first thredhold 
	double g2 = 0.5;//second thredhold=10* first thredhold
	double y1;//rate of numbers of n1 in the whole samples
	bool stop_sifting;
	bool stop_EMD = false;
	bool overThred2 = false;
	//for example, aktiv_imf = 2, means we got the 3.IMF from Signal.
	int aktiv_imf = 0;
	//every signal will be run by "sifting process" at most 2000 times
	int MAX_SIFT = 2000; 
	int sift=0;
	SV_Data *max_result = new SV_Data(max_imf,Len);

	MArray_1D(mean,Len,double,"array for average from envelope");
	MArray_1D(r,Len,double,"array for residue");
	MArray_1D(a_par,Len,double,"array for parameter a = (up-down)/2");
	MArray_1D(delta,Len,double,"array for parameter delta = mean/a");
	MArray_1D(h,Len,double,"array for signal,which will be tested, Imf or not");

	//from beginning we do not konw how many IMFs we can find, 
	//so we at first define a large Matrix(temp) to store the IMFs and count the numbers of IMFs(aktiv_imf)
	//at the end we copy all IMFs to the final Matrix "result(aktiv_imf,Len)".

	//flag=-2 means, sifting prozess must repeat, we got no IMF
	//flag=-1 means, we got an IMF, but we should continue to find the next IMF
	//flag=0  means, we got all IMF, and residue r[] is monoton.
	int flag=-1;
	double max=0.0;
	double min=0.0;
	double mid=0.0;

	for(long int i=0; i<Len; i++)	{
	r[i] = 0.0;
	}

	//we define first residue r0 = Signal X(t)
	for(long int i=0; i<Len; i++)	
	{	
		r[i]=Signal[i];
		if(max<r[i]) max = r[i];
		else if (min>r[i]) min = r[i];
	}

//mid is the largest value between maxima(local) and minima(local).
	mid = fabs(max-min);

	double *up_spline, *down_spline;

//EMD process, to find all IMF
while((flag!=0)&&(aktiv_imf<(max_imf))&&(!stop_EMD))
{
  //
	//init the number of Max. und Min. points
	nh_r=0;
	nl_r=0;
	overThred2 = false;
	stop_EMD = false;

	for(long int i=1; i<Len-1; i++)	{
		if((r[i-1]<r[i])&&(r[i]>r[i+1]))	nh_r++;
		if((r[i-1]>r[i])&&(r[i]<r[i+1]))	nl_r++;
	}

	//when no (not enough, too many)Max. or Min. points found, or maximal IMF found,stop the algorithmus and exit.
	if((nh_r==0)||(nl_r==0)||((nh_r+nl_r)<=7)||(aktiv_imf>=max_imf-1)||(nh_r+nl_r)>70000) 
	{ 	
		flag=0;
		break;
	}

	nh_pre = nh_r;
	nl_pre = nl_r;

	for(long int i=0; i<Len; i++) 	h[i]=r[i];
	stop_sifting = false;

//sifting process
  while((sift<=MAX_SIFT)&&(!stop_sifting)&&(!stop_EMD)){
	this->up_j=0;
	this->down_j=0;
	nh_h=0;
	nl_h=0;
	double zero_h = 0.0;
	n_zero = 0;
	stop_sifting = false;
	overThred2 = false;

	//too many times sifting
	if(sift>=MAX_SIFT){
	for(long int i=0;i<Len;i++){	
		max_result->Mat[aktiv_imf][i]=(float)(h[i]); //by thilo: explicit cast to avoid warning
		r[i]=r[i]-h[i];
		}
	  aktiv_imf++;
	  sift=0;
	  flag = -1;
	  stop_sifting = true;
	  break;
	}

//calculate the numbers of all Extrema, and zero-crossing
	for(long int i=1; i<Len-1; i++)	{
		if((h[i-1]<h[i])&&(h[i]>h[i+1]))	nh_h++;
		if((h[i-1]>h[i])&&(h[i]<h[i+1]))	nl_h++;		
	}

	for(long int i=1; i<Len-1; i++)	{
	  if ((h[i-1]*h[i])<0) n_zero++;
	  if ((h[i-1]*h[i])==0) {
		if((h[i-1]*h[i+1])<0) n_zero++;
		}
	  if ((h[i-1]*h[i])>0) continue;
	}

	long int *up_xs,*down_xs;
	double *up_ys,*down_ys;

//when too many Extrem or more and more Extrema then last Sift, stop the sifting process	
	if(((nh_h+nl_h)>=70000)||((nh_pre+nl_pre)<(nh_h+nl_h))){
		for(long int i=0;i<Len;i++){	
		max_result->Mat[aktiv_imf][i]=(float)(h[i]); 
		r[i]=r[i]-h[i];
		}
	  aktiv_imf++;
	  sift=0;
	  flag = -1;
	  stop_sifting = true;
	  break;
	}
	else{
	MArray_1D(up_xs,nh_h+2,long int,"array for position of maxima(local)");
	MArray_1D(down_xs,nl_h+2,long int,"array for position of minima(local)");
	MArray_1D(up_ys,nh_h+2,double,"array for maxima(local)");
	MArray_1D(down_ys,nl_h+2,double,"array for minima(local)");
	MArray_1D(up_spline,Len,double,"spline of up-envlope ");
	MArray_1D(down_spline,Len,double,"spline of down-envlope ");
	}

  //store all Extrema
  //store all position of local Max to up_xs
  //store all value of local Max to up_ys
  //store all position of local Min to down_xs
  //store all value of local Min to down_ys
  max_min(h,up_xs,down_xs,up_ys,down_ys,Len);

  //interpolate up- and down-envlope
  //store up-envlope to up_spline
  //store down-envlope to down_spline
  interpolation(up_xs,up_ys,this->up_j+2,Len,up_spline);
  interpolation(down_xs,down_ys,this->down_j+2,Len,down_spline);
  
  MFree_1D(up_xs);
  MFree_1D(up_ys);
  MFree_1D(down_xs);
  MFree_1D(down_ys);

  	for(long int i=0; i<Len; i++)
	{
	//we calculate the average of 2 envelope	
		mean[i] = (up_spline[i]+down_spline[i])/(double)2.0;
	//we define an amplitude-function a[], and an envaluation-function delta[]
		a_par[i]= (up_spline[i]-down_spline[i])/(double)2.0;
		if(a_par[i]!=0.0)
			delta[i] = fabs(mean[i]/a_par[i]);
		else delta[i] = 10.0;
	}

	MFree_1D(up_spline);
	MFree_1D(down_spline);

	n1=0;//number of samples,which are smaller than first thredhold
	y1=0.0;//rate of numbers of n1 in the whole samples

	//when at least one time point with delta>thredhold(0.5), it is not IMF,and continue sifting prozess 
	for(long int i=0;i<Len;i++){
	  if(delta[i]>g2) {
		overThred2 = true; 
		break;
	  }
	}

	if(overThred2) {
	  for(long int i=0;i<Len;i++)	h[i]=h[i]-mean[i];
		sift++;
		nh_pre = nh_h;
		nl_pre = nl_h;
		flag=-2;
		stop_sifting = false;
		continue;
	  }

	//when the difference from Extrema and zero-crossing more than 1, it is not IMF and continue sifting prozess
	if(abs(nh_h+nl_h-n_zero)>1){
	  	for(long int i=0;i<Len;i++)	h[i]=h[i]-mean[i];
		sift++;
		nh_pre = nh_h;
		nl_pre = nl_h;
		flag=-2;
		stop_sifting = false;
		continue;
	}
	else{
  //calculate the numbers of samples (n1), which delta-function smaller than 0.05 
	for(long int i=0;i<Len;i++){
		if(delta[i]<=g1) n1++;
	}
	//calculate the percent of n1, which delta-function smaller than 0.05 , 
	double n1_d = (double)n1;  
	y1 = (n1_d/length_d);	
	//test, whether 95% smaller then 0.05, when yes, we found an IMF,stop criterium by Rilling
	if(y1>=0.95) {
	for(long int i=0;i<Len;i++){	
		max_result->Mat[aktiv_imf][i]=(float)(h[i]); //by thilo: explicit cast to avoid warning
		r[i]=r[i]-h[i];
		}
	  aktiv_imf++;
	  sift=0;
	  flag = -1;
	  stop_sifting = true;
	  break;
	  }	
	else {//wenn not, continue sifting
		for(long int i=0;i<Len;i++)	h[i]=h[i]-mean[i];
			sift++;
			nh_pre = nh_h;
			nl_pre = nl_h;
			flag=-2;
			stop_sifting = false;
		  }
		}
  }
  
 //test, when Oscillation of the Residue is in the Threshold(Experiment value 10^-10) 
 if(tooSmall_simple(r,Len,0.0000000001,mid)) {
		stop_EMD = true;
		break;
	}
 }

  //copy all IMFs and residue to a new SV_data
  if(aktiv_imf>=max_imf-1) {
	  for(long int j=0;j<Len;j++){
		  max_result->Mat[max_imf-1][j]=(float)(r[j]);
	  }
	  return max_result;
  }
  else {
  SV_Data *result=new SV_Data(aktiv_imf+1,Len);

  for(long int i=0;i<aktiv_imf;i++)
		{
			for( long int j=0;j<Len;j++)
			{
				result->Mat[i][j]=max_result->Mat[i][j];
			}
		}
  for( long int j=0;j<Len;j++){
	result->Mat[aktiv_imf][j]=(float)(r[j]);
  }
  MFree_1D(max_result);
  MFree_1D(mean);
  MFree_1D(delta);
  MFree_1D(h);
  MFree_1D(a_par);
  MFree_1D(r);
  return result;
  }
 }

SV_Data *SC_HHT::insten_frequency(SV_Data* imf,long int row, long int col,double SampleRate,SV_Data *F,SV_Data *A){
  float tmp=0.0;
  SV_Data *F_tmp = new SV_Data(row,col);
  SV_Data *H_tmp = new SV_Data(row,col);

  //wenn imf[i][j]=0, change to a small value, to avoid null-division
  for(long int i=0;i<imf->Row;i++)
	for(long int j=0;j<imf->Col;j++)
	  if(imf->Mat[i][j]==0.0) imf->Mat[i][j] = (float)(0.001);

  //calculate Y[t] = (1/PI)*Ingration(C(t2)/(t1-t2))(dt2)
  for(long int i=0;i<F->Row;i++){
	for(long int j=0;j<F->Col;j++){
	  for(long int k=0;k<F->Col;k++){
		if(k!=j) tmp= tmp+(imf->Mat[i][k]/(j-k));
		}
	  H_tmp->Mat[i][j]=(float)(tmp/sclib::pi);
	  if((H_tmp->Mat[i][j]<0.001)&&(H_tmp->Mat[i][j]>-0.001)) F_tmp->Mat[i][j] = 0.0;
  //calculate Seta(t) = acrtan(Y(t)/C(t))
	  else F_tmp->Mat[i][j]=fabs(atan(H_tmp->Mat[i][j]/imf->Mat[i][j]));
	  tmp = 0.0;
	}
  }

  //calculate A(t) = sqrt(Y(t)^2+C(t)^2)
  //calculate w(t) = d(Seta(t))/dt
  //calculate f(t) = w(t)/(2*PI)
  for(long int i=0;i<F->Row;i++){
	for(long int j=0;j<F->Col;j++){
	  A->Mat[i][j] = (float) sqrt(imf->Mat[i][j]*imf->Mat[i][j]+H_tmp->Mat[i][j]*H_tmp->Mat[i][j]);
	  F->Mat[i][j] =(float)((int)((F_tmp->Mat[i][j]*SampleRate*SampleRate)/(2*col*sclib::pi)));
	}
  }

  MFree_0D(F_tmp);
  MFree_0D(H_tmp);

  //calculate the averange value of je 256 Samples(a Frame with 16000 SampleRate)
  int c = (int)(imf->Col/256);
  SV_Data *F_frame=new SV_Data(imf->Row,c);

  for(long int i=0;i<imf->Row;i++){
	for(long int k=0;k<c;k++){
	  for(long int j=k*256;j<(k+1)*256;j++){
		F_frame->Mat[i][k]+=F->Mat[i][j];
	  }
	F_frame->Mat[i][k] = (float)((int)((F_frame->Mat[i][k]/256)));
	}
  }

  //all instantaneous frequency values pro Sample are stored in Matrix F
  //all instantaneous frequency values pro Frame are stored in Matrix F_frame and return to the super class
  return F_frame;
}


SV_Data* SC_HHT::insten_frequency_op(SV_Data* imf,long int row, long int col,double SampleRate,SV_Data *F){
  //SV_Data *F = new SV_Data(row,col);//the Matrix of instantaneous frequency of all IMFs with time interval 0 to col

  double *abs_imf,*n_imf,*dq,*devi,*up_ys,*down_ys,*up_spline, *down_spline,*omcos;
  int Nnormal;
  double rangetop,rangebott;
  long int nh,nl;
  long int *up_xs,*down_xs;

  MArray_1D(abs_imf,col,double,"array for positive value of an imf");
  MArray_1D(n_imf,col,double,"array for positive value of an imf");
  MArray_1D(dq,col,double,"array for value of 1-imf[i]*imf[i]");
  MArray_1D(devi,col,double,"");
  MArray_1D(down_spline,col,double,"spline of down-envlope ");
  MArray_1D(omcos,col,double,"");

  for(long int i=0;i<row;i++){
	Nnormal=5;
	rangetop=0.90;

	for(long int j=0;j<col;j++){
		if(imf->Mat[i][j]<0)
		  abs_imf[j]=-imf->Mat[i][j];
		else abs_imf[j]=imf->Mat[i][j];
	  }

	for(int k=0;k<Nnormal;k++){
	  nh=0;
	  nl=0;
	  this->up_j=0;
	  this->down_j=0;

	  for(long int j=1;j<col-1;j++){
	  if((abs_imf[j-1]<abs_imf[j])&&(abs_imf[j]>abs_imf[j+1])) nh++;
	  else if ((abs_imf[j-1]>abs_imf[j])&&(abs_imf[j]<abs_imf[j+1])) nl++;
	}

	  MArray_1D(up_xs,nh+2,long int,"array for position of maxima(local)");
	  MArray_1D(down_xs,nl+2,long int,"array for position of minima(local)");
	  MArray_1D(up_ys,nh+2,double,"array for maxima(local)");
	  MArray_1D(down_ys,nl+2,double,"array for minima(local)");
	  MArray_1D(up_spline,col,double,"spline of up-envlope ");

	  max_min(abs_imf,up_xs,down_xs,up_ys,down_ys,col);

	  interpolation(up_xs,up_ys,this->up_j+2,col,up_spline);

	  for(long int j=0;j<col;j++){
	    if(up_spline[j]<0.0001&&up_spline[j]>-0.0001) abs_imf[j]=0;
	    else abs_imf[j] = abs_imf[j]/up_spline[j];
	  }
	MFree_1D(up_xs);
	MFree_1D(down_xs);
	MFree_1D(up_ys);
	MFree_1D(down_ys);
	MFree_1D(up_spline);
	}
	this->up_j=0;
	this->down_j=0;

	for(long int j=0;j<col;j++){
	  n_imf[j]=abs_imf[j];
	  if(imf->Mat[i][j]<0)
		n_imf[j]=-abs_imf[j];
	}
	
	for(long int j=0;j<col;j++)
	  dq[j]=sqrt(1-n_imf[j]*n_imf[j]);  
	
	for(long int j=1;j<col-1;j++){
	  devi[j]=n_imf[j+1]-n_imf[j-1];
	  if((devi[j]>0)&&(n_imf[j]<1))
		dq[j]=-dq[j];
	}
 
	rangebott=-rangetop;  

	for(long int j=1;j<col-1;j++){
	  if(n_imf[j]>rangetop) omcos[j]=-9999;
	  else if(n_imf[j]<rangebott) omcos[j]=-9999; 
	  else 
		omcos[j]=abs(n_imf[j+1]-n_imf[j-1])*0.5/sqrt(1-n_imf[j]*n_imf[j]);
	}
	omcos[0]=-9999;
	omcos[col-1]=-9999;


//	for(long int j=1;j<col-1;j++)
//	  omcos[j]=abs(n_imf[j+1]-n_imf[j-1])*0.5/sqrt(1-n_imf[j]*n_imf[j]);



	for(long int j=0;j<col;j++)
	  F->Mat[i][j] = (float)(omcos[j]*SampleRate*SampleRate/(2*sclib::pi*col));
  }
	
  for(long int i=0;i<F->Row;i++)
	for(long int j=0;j<F->Col;j++)
	  if((F->Mat[i][j]<=0.00001)&&(F->Mat[i][j]>-0.00001)) F->Mat[i][j]=0.0;

  for(long int i=0;i<F->Row;i++)
	for(long int j=0;j<F->Col;j++)
	  F->Mat[i][j]=(float)((int)(F->Mat[i][j]));

  int c = (int)(imf->Col/256);
  SV_Data *Ftime=new SV_Data(imf->Row,c);

  for(long int i=0;i<imf->Row;i++)
	for(long int k=0;k<c;k++)
	  Ftime->Mat[i][k]=0.0;

  for(long int i=0;i<imf->Row;i++){
	for(long int k=0;k<c;k++){
	  for(long int j=k*256;j<(k+1)*256;j++){
		//if(F->Mat[i][j]>Ftime->Mat[i][k]) Ftime->Mat[i][k]=F->Mat[i][j];
		Ftime->Mat[i][k]+=F->Mat[i][j];
	  }
	Ftime->Mat[i][k] = (float)((int)(Ftime->Mat[i][k]/256));
	}
  }


  MFree_1D(down_spline);
  MFree_1D(abs_imf);
  MFree_1D(n_imf);
  MFree_1D(dq);
  MFree_1D(devi);
  MFree_1D(omcos);

  return Ftime;
}

SV_Data* SC_HHT::hilbert_spectrum(SV_Data *frequency, SV_Data *A, SV_Data *imf,long int row, long int col){
	long int max_fre = 0;
	long int **fre_int;
	long int temp_fre = 0;
	MArray_2D(fre_int,row,col,long int,"array for frequency with type int");

//change the type of Frequency, "double"->"int" 
	for(long int i=0;i<row;i++){
		for(long int j=0;j<col;j++){ 
		  fre_int[i][j] = (long int)fabs(frequency->Mat[i][j]/10);
			//we got the largest Frequency value
		  if (fre_int[i][j]>max_fre) max_fre = fre_int[i][j];
		}
	}

//we create a new Matrix with max_fre+1	rows (a row with 10 Hz changed)
	if(max_fre>500) max_fre=500;

	SV_Data *result = new SV_Data(max_fre+1,col);
	
//initilation	
	for(long int i=0;i<=max_fre;i++)
		for(long int j=0;j<col;j++)
			result->Mat[i][j] = 0.0;
//copy the amplitude to the suitable Frequency-row in every time-point
	for(long int i=0;i<row;i++){
		for(long int j=0;j<col;j++)
		{	temp_fre = fre_int[i][j];
			if(temp_fre>max_fre) temp_fre=max_fre;
			result->Mat[temp_fre][j] += fabs(A->Mat[i][j]);
		}
	}
	MFree_2D(fre_int);
	return result;
}

//use this method to interpolate the cubic spline with Extrema positions *xs und values *ys
void SC_HHT::interpolation(long int *xs,double *ys,long int Length,long int Len,double *spline)
{
	double *ys2,*temp;
	double p,sig,a,b,c,d,e,f,g,a0,a1,a2,a3,xc,aaa,bbb,ccc;
	long int prev,cur,j,jfin;

	MArray_1D(ys2,Length,double,"");
	MArray_1D(temp,Length,double,"");

  ys2[0]=0.0;
  temp[0]=0.0;
  for (long int i=1;i<Length-1;i++) {
    sig=(xs[i]-xs[i-1])/(xs[i+1]-xs[i-1]);
    p=sig*ys2[i-1]+2.0;
    ys2[i]=(sig-1.0)/p;
	aaa = (ys[i+1] - ys[i]) / (double)(xs[i+1] - xs[i]);
    bbb = (ys[i] - ys[i-1]) / (double)(xs[i] - xs[i-1]);
    ccc = (double)xs[i+1] - xs[i-1];
    temp[i] = (6.0 * (aaa - bbb) / ccc - sig * temp[i-1])/p;
  }

  ys2[Length-1] = 0.0;  

  //get 2-th diff value
  for (long int k=Length-2;k>=0;k--) {
	ys2[k]=ys2[k]*ys2[k+1]+temp[k];
  }

  MFree_1D(temp);

  /* Compute the spline coefficients */
  cur=0;
  j=0;
  jfin=Length-1;
  while (xs[j+1]<0) j++;
  while (xs[jfin-1]>Len) jfin--;


  for (;j<=jfin;j++) {
    /* Compute the coefficients of the polynomial between two knots */
    a=xs[j];
    b=xs[j+1];
    c=b-a;
    d=ys[j];
    e=ys[j+1];
    f=ys2[j];
    g=ys2[j+1];
    a0=(b*d-a*e+(b*b*b)*f/6-(a*a*a)*g/6)/c+c*(a*g-b*f)/6;
    a1=(e-d-(b*b)*f/2+(a*a)*g/2)/c+c*(f-g)/6;
    a2=(b*f-a*g)/(2*c);
    a3=(g-f)/(6*c);

    prev=cur;
    while ((cur<Len-1) && ((j==jfin) || (cur<xs[j+1]))) cur++;

    /* Compute the value of the spline at the sampling times between up_xs[j] and up_xs[j+1] */
    for (long int i=prev;i<cur;i++) {
      xc=i;
      spline[i]=a0+a1*xc+a2*(xc*xc)+a3*(xc*xc*xc);
	}
  }
  MFree_1D(ys2);
}

//use this method to test, whether the Amplitude of Residue too small is
bool SC_HHT::tooSmall_simple(double *res,long int length,double alfa,double Thredhold){
	bool tooSmall = true;
	long int nh=0;
	long int nl=0;
	long int up_k=0;
	long int down_k=0;
	long int *up_x,*down_x;
	double *up_y,*down_y;
	double z=0.0;
	double a=0.0;
	double s=0.0;
	for(long int i=1; i<length-1; i++)	{
		if((res[i-1]<res[i])&&(res[i]>res[i+1]))	nh++;
		if((res[i-1]>res[i])&&(res[i]<res[i+1]))	nl++;
	}

	if(nh==0||nl==0||((nh+nl)<=7)||((nh+nl)>=40000)) return true;
	MArray_1D(up_x,nh,long int,"array for position of maxima(local)");
	MArray_1D(down_x, nl,long int,"array for position of minima(local)");
	MArray_1D(up_y,nh,double,"array for maxima(local)");
	MArray_1D(down_y, nl,double,"array for minima(local)");

	for(long int i=1; i<length-1; i++)	{
		if((res[i-1]<res[i])&&(res[i]>res[i+1]))	{
			up_x[up_k]=i;
			up_y[up_k]=res[i];
			up_k++;
		}
		else if((res[i-1]>res[i])&&(res[i]<res[i+1]))	{
			down_x[down_k]=i;
			down_y[down_k]=res[i];
			down_k++;
		}
	}
	up_k = 0;
	down_k = 0;

	for(long int i=0;i<nh;i++) 
	{	a+=up_y[i];
		z+=down_y[i];
	}

	s=fabs((a-z)/(double)(nh+nl));

	if (s<=(alfa*Thredhold)) tooSmall = true;
	else tooSmall = false;
	
  MFree_1D(up_x);
  MFree_1D(down_x);
  MFree_1D(up_y);
  MFree_1D(down_y);

  return tooSmall;	
}

//use this method to find all Maxima and Minima, and in suitable Matrix stored
void SC_HHT::max_min(double *h,long int *up_xs,long int *down_xs,double *up_ys,double *down_ys,long int Len){
  
	double up_s1,up_s2,down_s1,down_s2;
	double up_temp1,up_temp2,down_temp1,down_temp2;

	up_xs[0] = 0;
	down_xs[0] = 0;
	up_ys[0] = h[0];
	down_ys[0] = h[0];

	for(long int i=1; i<Len-1; i++)	{
		if((h[i-1]<h[i])&&(h[i]>h[i+1]))	{
			up_xs[this->up_j+1]=i;
			up_ys[this->up_j+1]=h[i];
			this->up_j++;
		}
		else if((h[i-1]>h[i])&&(h[i]<h[i+1]))	{
			down_xs[this->down_j+1]=i;
			down_ys[this->down_j+1]=h[i];
			this->down_j++;
		}
	}
	up_xs[this->up_j+1] = Len-1;
	down_xs[this->down_j+1] = Len-1;
	up_ys[this->up_j+1] = h[Len-1];
	down_ys[this->down_j+1] = h[Len-1];

	up_s1 = (up_ys[1]-up_ys[2])/(up_xs[1]-up_xs[2]);
	up_temp1 = up_s1*(up_xs[0]-up_xs[1])+up_ys[1];
	if(up_temp1>up_ys[0]) 
	  up_ys[0] = up_temp1;

	up_s2 = (up_ys[this->up_j]-up_ys[this->up_j-1])/(up_xs[this->up_j]-up_xs[this->up_j-1]);
	up_temp2 = up_s2*(up_xs[this->up_j+1]-up_xs[this->up_j])+up_ys[this->up_j];
	if(up_temp2>up_ys[this->up_j+1])
	  up_ys[this->up_j+1] = up_temp2;

	down_s1 = (down_ys[1]-down_ys[2])/(down_xs[1]-down_xs[2]);
	down_temp1 = down_s1*(down_xs[0]-down_xs[1])+down_ys[1];
	if(down_temp1<down_ys[0])
	  down_ys[0] = down_temp1;

	down_s2 = (down_ys[this->down_j]-down_ys[this->down_j-1])/(down_xs[this->down_j]-down_xs[this->down_j-1]);
	down_temp2 = down_s2*(down_xs[this->down_j+1]-down_xs[this->down_j])+down_ys[this->down_j];
	if(down_temp2<down_ys[this->down_j+1])
	  down_ys[this->down_j+1] = down_temp2;
}
