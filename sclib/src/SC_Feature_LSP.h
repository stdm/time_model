/**************************************************************************/
/*    This class implements the Line Spectra Pairs feature, that is       */
/*    derived from LPCs (or the raw signal in this case, too). It         */
/*    includes the same information as LPCs but inherits more robustness. */
/*    The implementation is based on the corresponding speex code         */
/*    (http://www.speex.org/) or ephone MELP Proposed Federal Standard    */
/*    speech coder code and not much more than a wrapper around both.     */
/*    (see copyright etc. in the corresponding .cpp file)                 */
/*																																				*/
/*    Author  : Thilo Stadelmann            															*/
/*    Date    : 12.02.2006																								*/
/**************************************************************************/

#ifndef __SC_Feature_LSP_H__
#define __SC_Feature_LSP_H__

#include "SC_Api.h"
#include <SV_Feature.h>

class SC_Feature_LSP : public SV_Feature {

	private :

		//from speex:
#ifdef SC_USE_SPEEXLSP
		float cheb_poly_eva(float *coef,float x,int m,char *stack);
		int lpc_to_lsp (float *a,int lpcrdr,float *freq,int nb,float delta, char *stack);
		void lsp_to_lpc(float *freq,float *ak,int lpcrdr, char *stack);
		void lsp_enforce_margin(float *lsp, int len, float margin);
#endif

		//from ephone/melp:
#ifdef SC_USE_MELPLSP
		int lpc_pred2lsp(float *a,float *w,int p,float lsp_delta = 0.0,float root_delta = 0.00781250,int root_bisections = 4,int clmp_max_loops = 10);
		int lsp_roots(float *w,float **c,int p2,float delta = 0.00781250,int bisections = 4);
		float lsp_g(float x,float *c,int p2);
		int lpc_clmp(float *w, float delta, int p, int max_loops = 10);
#endif

	protected :

		int method; //choice of algorithm to derive LSPs from LPCs: SPEEX or MELP (SPEEX recommended)
		double delta; //delta for root-searching
		int bisections; //number of bisections during root-search
		double minSeparation; //minimum separation delta for LSPs and MELP method 
		int maxLoops; //maximum number of loops of separation check and MELP method

	public :

    //====================================================================================================================
    // constructor / destructor
    //====================================================================================================================
		SC_Feature_LSP(int sampleRate, int frameLength, int frameStep, unsigned int window, double preemphasize, int LPCorder, int method, double delta = 0.00781250, int bisections = 4, double minSeparation = 0.0, int maxLoops = 10);
		virtual ~SC_Feature_LSP();

    //====================================================================================================================
		// override base class method, return log-spectral vector sequence
    //====================================================================================================================
		virtual SV_Data *ExtractFeature(void);
    
    //====================================================================================================================
		// if LPCs are already computed, just use them
    //====================================================================================================================
    SV_Data *ExtractFeature(SV_Data *pLPCs);
};

#endif
