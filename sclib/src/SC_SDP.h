/**************************************************************************/
/*    This class creates Symmetrized Dot Patterns from input signals as   */
/*    proposed by Clifford Pickover in "Computes, Pattern, Chaos and      */
/*    Beauty"                                                             */
/*                                                                        */
/*    Author  : Bing Shi        																					*/
/*    Date    : 15.12.2006																								*/
/**************************************************************************/

#ifndef __SC_SDP_H__
#define __SC_SDP_H__

#include "SC_Api.h"
#include "SC_Aux.h"
#include "SC_Signature.h"
#include <SV_Error.h>
#include <SV_Data.h>

class SCLIB_API SC_SDP {

	private:
		int m;
		int lag;
		int color;
		int n;
		int pictureSize;
		int tau;
		
		int approximate(double data);

	protected:
        
  public:

    SC_SDP(int m, int lag, int color, int n, int pictureSize, int tau);
    ~SC_SDP();

		//====================================================================================================================
		//	save and load the 2-d Arrary in a bmp-file;
		//====================================================================================================================
		bool saveBitMap(char* fileName, unsigned long int **data, bool withColor);
		bool saveBitMap(char* fileName, unsigned long int **data, int X, int Y);
		unsigned long int** loadBitMap(char* fileName);
   
		//====================================================================================================================
		// create the SDP with the data
		// Get signal pro frame and the length of the signal
		// return the SDP in BMP
		//====================================================================================================================
		unsigned long int ** createSDP(double* signal, unsigned long int length, bool withColor, bool bitMapLike = false);

		SV_Data * sdp2sv_data(unsigned long int ** pSDP, int frameSize, int frameStep); //by thilo: frameSize and frameStep are only here to fill the SV_Data's header
		void sdp2vector(unsigned long int ** pSDP, float * vector);
		unsigned long int ** vector2sdp(float *vectors);
		SC_Signature* sv_data2signature(SV_Data * pData);
		SC_Signature** sv_data2signatures(SV_Data *pData);
		SC_Signature* vector2signature(float *vector, int length);
		float* signature2vector(SC_Signature *pSingature);
		void signature2vector(SC_Signature *pSingature, float *vector); //by thilo: same as above, but relies on an already alocated vector
		SV_Data* signature2sv_data(SC_Signature *pSignature, int frameSize, int frameStep); //by thilo: frameSize and frameStep are only here to fill the SV_Data's header
		SV_Data* signatures2sv_data(SC_Signature **pSignatures, int numOfSignatures, int frameSize, int frameStep); //by thilo: frameSize and frameStep are only here to fill the SV_Data's header
		unsigned long int ** sv_data2sdp(SV_Data *pData);

		unsigned long int rgb2long(unsigned short R, unsigned short G, unsigned short B);
		void long2rgb(unsigned long data, unsigned short &r, unsigned short &g, unsigned short &b);

		double getHighValue (double * signal,  unsigned long int length);
		double getLowValue(double* signal, unsigned long int length);
};

#endif
