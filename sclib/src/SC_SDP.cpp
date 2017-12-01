/**************************************************************************/
/*    This class creates Symmetrized Dot Patterns from input signals as   */
/*    proposed by Clifford Pickover in "Computes, Pattern, Chaos and      */
/*    Beauty"                                                             */
/*                                                                        */
/*    Author  : Bing Shi        																					*/
/*    Date    : 15.12.2006																								*/
/**************************************************************************/

#include "SC_SDP.h"
#include "SC_Aux.h" //for sclib::bufferSize define
#include "SC_MatrixFunctions.h"
#include "SC_Centroid_Point.h"

#include "SC_Api.h"
#ifdef SC_USE_EASYBMP
	#include <EasyBMP.h>
#else 
	#define ebmpBYTE unsigned long
#endif

#include <SV_Error.h>
#include <GN_Filter.h>

//====================================================================================================================
//	The constructor
//====================================================================================================================
SC_SDP::SC_SDP(int m, int lag, int color, int n, int pictureSize, int tau) {
	this->m = m;
	this->lag = lag; //lag in [samples] in the correlation analysis
	this->color = color;
	this->n = n;
	this->pictureSize = pictureSize;
	this->tau = tau; //waveform is scaled to [0..tau] 
}

//====================================================================================================================
//	The destructor
//====================================================================================================================
SC_SDP::~SC_SDP() {

}

//====================================================================================================================
// create the SDP with the data
// Get signal pro frame and the length of the signal
// return the SDP in BMP
//====================================================================================================================
unsigned long int ** SC_SDP::createSDP(double* signal, unsigned long int length, bool withColor, bool bitMapLike){
	unsigned long int ** sdp;
	unsigned long int  j;
	int x, y, i;
	double r, zeta, fie, zeta1,highValue, lowValue, t = 0.0;
	double positionX1, positionY1, positionX2, positionY2;
	int x1_bmp, y1_bmp, x2_bmp, y2_bmp;
	SC_MatrixFunctions mFunc;

	MArray_2D(sdp, this->pictureSize, this->pictureSize, unsigned long int, "SC_SDP.createSDP: sdp");

	// init. the sdp
	for(x = 0; x<this->pictureSize; x++){
		for(y=0; y<this->pictureSize; y++){
			sdp[x][y] = 0;
		}
	}

	highValue = mFunc.max(signal, length);
	lowValue = mFunc.min(signal, length);

	//TODO by thilo: what is this??? -> only stored for the last loop!
	for (j = 0; j<length - this->lag; j++) {
		r = ((signal[j] - lowValue) / (highValue - lowValue)) * this->tau;
		if (withColor == true) {
			t = signal[j+lag] - (int) signal[j];
		}// for Color

		for(i=1; i<= this->m; i++){
			zeta1 = (360/m) * i;
			zeta = zeta1 + (((signal[j+this->lag] - lowValue) / (highValue - lowValue)) * tau);
			fie = zeta1 - (((signal[j+this->lag] - lowValue) / (highValue - lowValue)) * tau);
			positionX1 = cos(( zeta / 180) * sclib::pi) * r;
			positionY1= sin(( zeta / 180) * sclib::pi) * r;
			positionX2 =cos(( fie / 180) * sclib::pi) * r;
			positionY2 = sin(( fie / 180) * sclib::pi) * r;
			// calculate the x und y for bmp
			x1_bmp = approximate(positionX1);
			y1_bmp = approximate(positionY1);
			x2_bmp = approximate(positionX2);
			y2_bmp = approximate(positionY2);	
			//save the positions in bmp
			if(withColor){
				if(t >0){
					unsigned short R = 0;
					unsigned short G = (unsigned short)(256 - (t / 2.0));
					unsigned short B =(unsigned short) t;
					unsigned long int c = rgb2long(R, G, B);
					sdp[x1_bmp][y1_bmp] = c;
					sdp[x2_bmp][y2_bmp] = c;
				}else{
					unsigned short R = (unsigned short)(256 + (t / 2.0));
					unsigned short G = 0;
					unsigned short B = (unsigned short) (t * -1);
					unsigned long int c = rgb2long(R, G, B);
					sdp[x1_bmp][y1_bmp] = c;
					sdp[x2_bmp][y2_bmp] = c;
				}
			}else{
				sdp[x1_bmp][y1_bmp]++;
				sdp[x2_bmp][y2_bmp]++;
			}
		}
	}

	//block by thilo
	if (bitMapLike==true && withColor==false) {
		float w = 255 / 64.0;
		for(x=0; x<this->pictureSize; x++){
			for(y=0; y<this->pictureSize; y++){
				if(sdp[x][y] > 0){ //by thilo: here, white is 0, not 255, hence the "255-..." part
					sdp[x][y] = 255 - (ebmpBYTE)(100 + (sdp[x][y]-1)*w); //50;
				}else{
					sdp[x][y] = 255 - 255;
				}
			}
		}
	}

	return sdp;
}

//====================================================================================================================
//	save the 2-d Arrary in a bmp-file;
//====================================================================================================================
bool  SC_SDP::saveBitMap(char* fileName, unsigned long int **data, bool withColor) {
#ifdef SC_USE_EASYBMP
	BMP bitMap;
	unsigned short int r, g, b; //by thilo

	bitMap.SetSize(this->pictureSize, this->pictureSize);
	bitMap.SetBitDepth(this->color);
  int x, y;
	if(withColor == true){
		for(x=0; x<this->pictureSize; x++){
			for(y=0; y<this->pictureSize; y++){
				if(data[x][y] > 0){
					//const unsigned short* rgb = long2rgb(data[x][y]);
					long2rgb(data[x][y], r, g, b); //by thilo
					bitMap(x, y)->Blue = (ebmpBYTE)r;
					bitMap(x,y)->Green =(ebmpBYTE)g;
					bitMap(x,y)->Red = (ebmpBYTE)b;
				}else{
					bitMap(x, y)->Blue =0;
					bitMap(x,y)->Green =0;
					bitMap(x,y)->Red =0;
				}
			}
		}
	}else{
		float w = 255 / 64.0;
		for(x=0; x<this->pictureSize; x++){
			for(y=0; y<this->pictureSize; y++){
				if(data[x][y] > 0){
					bitMap(x, y)->Blue = (ebmpBYTE) (100 + (data[x][y]-1)*w); //50;
					bitMap(x,y)->Green = (ebmpBYTE) (100 + (data[x][y]-1)*w); //50;
					bitMap(x,y)->Red =(ebmpBYTE) (100 + (data[x][y]-1)*w);
				}else{
					bitMap(x, y)->Blue = 255;
					bitMap(x,y)->Green = 255;
					bitMap(x,y)->Red = 255;
				}
			}
		}
	}
	bitMap.WriteToFile(fileName);
  return true;
#else
	return false;
#endif
}

bool SC_SDP::saveBitMap(char* fileName, unsigned long int **data, int X, int Y) {
#ifdef SC_USE_EASYBMP
	BMP bitMap;
	bitMap.SetSize(X, Y);
	bitMap.SetBitDepth(this->color);
  int x, y;
	float w = 255 / 64.0;
	for(x = 0; x < X; x++){
		for(y = 0; y < Y; y++){
			if(data[x][y] > 0){
				bitMap(x, y)->Blue = 50;
				bitMap(x,y)->Green =50;
				bitMap(x,y)->Red =(ebmpBYTE) (100 + (data[x][y]-1)*w);
			}else{
				bitMap(x, y)->Blue = 255;
				bitMap(x,y)->Green = 255;
				bitMap(x,y)->Red = 255;
			}
		}
	}
	bitMap.WriteToFile(fileName);
  return true;
#else
	return false;
#endif
}

unsigned long int ** SC_SDP::loadBitMap(char* fileName) {
#ifdef SC_USE_EASYBMP
	BMP bitmap;
	unsigned long int ** pSDP = NULL;
	int X = this->pictureSize, Y = this->pictureSize, x = 0, y = 0;
	unsigned short int r, g, b;

	bitmap.ReadFromFile(fileName);
	MArray_2D(pSDP, X, Y, unsigned long int, "SC_SDP::loadBitMap.pSDP");

	for(y = 0;  y < Y; y++){
		for(x = 0; x < X; x++){
			r = bitmap(x, y)->Red;
			g = bitmap(x, y)->Green;
			b = bitmap(x, y)->Blue;

			if(r == 255 && g == 255 && b == 255)
				pSDP[x][y] = 0;
			else
				pSDP[x][y] = 1;
		}
	}
	return pSDP;
#else
	return NULL;
#endif
}

SV_Data * SC_SDP::sdp2sv_data(unsigned long int** pSDP, int frameSize, int frameStep) {
	SV_Data * pData = NULL;
	pData = new SV_Data;
	if (pData==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}

	pData->Row = this->pictureSize;
	pData->Col = this->pictureSize;
	pData->Alloc();
	pData->Hdr.frameSize = frameSize;
	pData->Hdr.frameStep = frameStep;
  unsigned long x, y;
	
	// copy the pSDP to pData
	for (x = 0; x < (unsigned long)pData->Row; x++){
		for(y = 0; y < (unsigned long)pData->Col; y++){
				pData->Mat[x][y] = (float)pSDP[x][y];
		}
	}
 return pData;
}

void SC_SDP::sdp2vector(unsigned long int ** pSDP, float * vector){
	int length = this->pictureSize * this->pictureSize;

	for (int x = 0; x < this->pictureSize; x++){
		for(int y = 0; y < this->pictureSize; y++){
			vector[x*this->pictureSize + y] =(float)pSDP[x][y];
		}
	}

	return;
}

unsigned long int ** SC_SDP::sv_data2sdp( SV_Data* pData){
	if (this->pictureSize != pData->Row) {
		REPORT_ERROR(SVLIB_BadArg, "the length must be same");
	}
	unsigned long int ** sdp;
	MArray_2D(sdp,this->pictureSize, this->pictureSize, unsigned long int, "SC_SDP::sv_data2sdp");
	int x, y;
	// copy the pSDP to pData
	for (x = 0; x < this->pictureSize; x++){
		for(y = 0; y < this->pictureSize; y++){
				sdp[x][y] =(unsigned long int)pData->Mat[x][y];
		}
	}
 return sdp;
}

unsigned long int** SC_SDP::vector2sdp(float *vector){
	unsigned long int ** sdp;
	MArray_2D(sdp,this->pictureSize, this->pictureSize, unsigned long int, "SC_SDP::sv_data2sdp");
	int x, y;
	// copy the pSDP to pData
	for (x = 0; x < this->pictureSize; x++){
		for(y = 0; y < this->pictureSize; y++){
				sdp[x][y] =(unsigned long int)vector[x*this->pictureSize + y];
		}
	}
 return sdp;
}

SC_Signature *  SC_SDP::sv_data2signature(SV_Data* pData){
	SC_Centroid **centroids;
 // SC_Signature  pSignature = NULL; 
  double *weights, point[2], sum = 0.0;
	double smallestUnnormalizedWeight = std::numeric_limits<double>::max() ;
  long int x, X, y ,Y, n; 
  long int number = 0;
	X = (long int)pData->Row;
  Y =( long int)pData ->Col ;
	
	//Number of the elements in the pData
	for(x = 0; x < X; x++){
		for(y = 0; y < Y; y++){
			if(pData->Mat[x][y] > 0.0){
				number++;
				sum += pData->Mat[x][y];
			}
		}
	}

	MArray_1D(centroids, number, SC_Centroid*, "SC_SDP:sv_data2signature: centroids");
  MArray_1D(weights, number, double, "SC_SDP:sv_data2signature: weights");
  
	n = 0; 
	for(x = 0; x < X; x++){
		for(y = 0; y < Y; y++){
			if(pData->Mat[x][y] > 0.0){
				point[0] = x;
				point[1] = y;
				centroids[n] = new SC_Centroid_Point(NULL, 2, point, false); //by thilo: no tweak necessary here because distance metrix can only be the euclidean distance
				weights[n] = pData->Mat[x][y] / sum;
				if(smallestUnnormalizedWeight >= pData->Mat[x][y]){
					smallestUnnormalizedWeight = pData->Mat[x][y];
				}
				n++;
			}
			
		}
	}

	SC_Signature* pSignature = new SC_Signature(centroids, weights, number, true, smallestUnnormalizedWeight);

	return pSignature;
}

SC_Signature* SC_SDP::vector2signature(float *vector, int length){
	SC_Centroid **centroids;
  double *weights, point[2], sum = 0.0;
	double smallestUnnormalizedWeight = std::numeric_limits<double>::max() ;
  long int x, y, n; 
  long int number = 0;
	int pictureSize = (int)(std::sqrt((double)length));
	
	//Number of the elements in the vector
	for(int i = 0; i < length; i++){
		if(vector[i] > 0.0){
			number++;
			sum += vector[i];
		}
	}

	MArray_1D(centroids, number, SC_Centroid*, "SC_SDP:sv_data2signature: centroids");
  MArray_1D(weights, number, double, "SC_SDP:sv_data2signature: weights");
  
	n = 0;
	for(int i = 0; i < length; i++){
		if(vector[i] > 0.0){
			x = i / pictureSize;
			y = i % pictureSize;
			point[0] = x;
				point[1] = y;
				centroids[n] = new SC_Centroid_Point(NULL, 2, point, false); //by thilo: no tweak necessary here because distance metrix can only be the euclidean distance
				weights[n] = vector[i] / sum;
				if(smallestUnnormalizedWeight >= vector[i]){
					smallestUnnormalizedWeight = vector[i];
				}
				n++;
		}
	}

	SC_Signature* pSignature = new SC_Signature(centroids, weights, number, true, smallestUnnormalizedWeight);
	pSignature->setJustLink(false, false);

	return pSignature;
}

SC_Signature** SC_SDP::sv_data2signatures(SV_Data *pData){
	int count = pData->Row;
	int length = pData->Col;
	SC_Signature **pSignatures;
	MArray_1D(pSignatures, count, SC_Signature*, "SC_SDP:sv_data2signature: pSignatures");
	
	for(int i = 0; i < count; i++){
		pSignatures[i] = vector2signature(pData->Mat[i], length);
	}
	return pSignatures;
}

SV_Data* SC_SDP::signature2sv_data(SC_Signature *pSignature, int frameSize, int frameStep) {
	SV_Data * pData = NULL;
	pData = new SV_Data;
	if (pData==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}

	pData->Row = this->pictureSize; //(sigLen - this->Para.WinSz + this->Para.StpSz) / this->Para.StpSz;
	pData->Col = this->pictureSize;
	pData->Hdr.frameSize = frameSize;
	pData->Hdr.frameStep = frameStep;
	pData->Alloc();
	
  unsigned long x, y, t;
	
	//init
	for (x = 0; x < (unsigned long)pData->Row; x++){
		for(y = 0; y < (unsigned long)pData->Col; y++){
				pData->Mat[x][y] =0.0;
		}
	}
	// read centroids from signature into sv_data
	SC_Centroid_Point** temp = 	(SC_Centroid_Point**)pSignature->getCentroids();
	double* coord, *weight;
	double coordX, coordY, factor, minW;
	unsigned int w;
	weight = pSignature->getWeight();
	minW =  std::numeric_limits<double>::max() ;

	//seach the min weight for calculating the sum of the unnormalized Weight
	for (t = 0; t < (unsigned int)(pSignature->getN()); t++) {
		if (minW > weight[t] ) {
			minW = weight[t];
		}
	}
	factor = pSignature->getSmallestUnnormalizedWeight() / minW;
	for (t = 0; t < (unsigned int)(pSignature->getN()); t++) { //by thilo: added explicit typecast to avoid warning
		coord = temp[t]->getCoordinate();
		coordX = coord[0];
		coordY = coord[1];
		w = (unsigned int)(floor( weight[t] * factor + 0.5));
		coordX =( coordX > (this->pictureSize-1)) ? (this->pictureSize-1): coordX;
		coordY = (coordY > (this->pictureSize-1)) ? (this->pictureSize-1): coordY;
		pData->Mat[(unsigned int)coordX][(unsigned int)coordY] = (float)(w); //by thilo: explicit cast to float to avoid warning
	}

 return pData;
}

float* SC_SDP::signature2vector(SC_Signature *pSignature) { //by thilo
	float *vector = NULL;
	int length = this->pictureSize * this->pictureSize;

	MArray_1D(vector, length, float, "SC_SDP:sv_data2signature: vector");
	signature2vector(pSignature, vector);	

	return vector;
}

void SC_SDP::signature2vector(SC_Signature *pSignature, float *vector) { //by thilo: altered the function that it now accepts an already allocated vector and created the above stump that fullfills the old interface that returns a new vector
  unsigned t;
	int length = this->pictureSize * this->pictureSize;
	
	//init
	for(int i = 0; i < length; i++){
		vector[i] = 0.0;
	}

	// read centroids from signature into sv_data
	SC_Centroid_Point** temp = 	(SC_Centroid_Point**)pSignature->getCentroids();
	double* coord, *weight;
	double  factor, minW;
	int coordX, coordY;
	unsigned int w;
	weight = pSignature->getWeight();
	minW =  std::numeric_limits<double>::max() ;

	//seach the min weight for calculating the sum of the unnormalized Weight
	for (t = 0; t < (unsigned int)(pSignature->getN()); t++) {
		if (minW > weight[t] ) {
			minW = weight[t];
		}
	}

	factor = pSignature->getSmallestUnnormalizedWeight() / minW;

	for (t = 0; t < (unsigned int)(pSignature->getN()); t++) { //by thilo: added explicit typecast to avoid warning
		coord = temp[t]->getCoordinate();
		coordX = (int)coord[0];
		coordY = (int)coord[1];
		w = (unsigned int)(floor( weight[t] * factor + 0.5));
		coordX =( coordX > (this->pictureSize-1)) ? (this->pictureSize-1): coordX;
		coordY = (coordY > (this->pictureSize-1)) ? (this->pictureSize-1): coordY;
		vector[int(coordX * this->pictureSize + coordY)] = (float)(w); //by thilo: explicit cast to float to avoid warning
	}

	return;
}

SV_Data* SC_SDP::signatures2sv_data(SC_Signature **pSignatures, int numOfSignatures, int frameSize, int frameStep){
	SV_Data * pData = NULL;
	pData = new SV_Data;
	if (pData==NULL) {
		REPORT_ERROR(SVLIB_NoMem, "No memory for DataSet");
	}

	pData->Row = numOfSignatures;
	pData->Col = this->pictureSize * this->pictureSize;
	pData->Alloc();
	pData->Hdr.frameSize = frameSize;
	pData->Hdr.frameStep = frameStep;

	for(int i = 0; i < numOfSignatures; i++){
		//pData->Mat[i] = signature2vector(pSignatures[i]); //by thilo: TODO: this will give error while freeing because Mat[i] is already an allocated vector, but the pointer to it will be lost after this operation!!!
		signature2vector(pSignatures[i], pData->Mat[i]); //by thilo to overcome the above error
	}
	
	return pData;
}

//approximate the coordinate system
int SC_SDP::approximate(double data){
	int pos = sclib::round(data * (this->pictureSize/2/this->tau)) + (this->pictureSize/2) -1;
	
	if (pos < 0) {
		pos = 0;
	} else if (pos > this->pictureSize -1) {
		pos = this->pictureSize-1;
	}

	return pos;
}

unsigned long int SC_SDP::rgb2long(unsigned short R, unsigned short G, unsigned short B){
 unsigned long int data =( unsigned long int) (R*255*255 + G*255 + B);
 return data;
}

void SC_SDP::long2rgb(unsigned long int data, unsigned short &r, unsigned short &g, unsigned short &b){ //by thilo: changed way of returning the rgb values so that no pointer to a local variable has to be transmitted
	/*static unsigned short rgb[3];//rgb[0] = red, rgb[1] = green, rgb[2] = blue
	//MArray_1D(rgb,3,unsigned short,"SC_SDP::long2rgb.rgb");
	rgb[0] =(unsigned short) (data / ( 255*255));
	rgb[1] = (unsigned short)((data - rgb[0]*255*255)/255);
	rgb[2] = (unsigned short)(data- ( rgb[0]*255*255 + rgb[1]*255));
	return rgb; //TODO by thilo: this array is local, it can't be returned (allocate dynamic memory via MArray_1D() to do this...)
	*/
	
	r =(unsigned short) (data / ( 255*255));
	g = (unsigned short)((data - r*255*255)/255);
	b = (unsigned short)(data- ( r*255*255 + g*255));

	return;
}
