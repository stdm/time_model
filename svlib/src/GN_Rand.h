//########################################################################
//  
// A C++ class library for automatic speech recognition and 
// speaker recognition (identification and verification). 
//
// This class library is provided "as is" without any express 
// or implied warranty of any kind with respect to this software. 
// In particular the author shall not be liable for any direct, 
// indirect, special, incidental or consequential damages arising 
// in any way from use of the software.
//
//
// Author   : Jialong He,  Copyright (C), all rights reserved. 
// Date     : May, 1999. 
//
// Contact  : Jialong_He@bigfoot.com, Jialong_He@homemail.com
// Web Page : www.bigfoot.com/~Jialong_He
//########################################################################

//************************************************************************
//    Generate uniform or normal distribution random numbers.
//
//
//    Author  : Jialong HE
//    Date    : April 27, 1999
//************************************************************************
#ifndef __GN_Rand_H__
#define __GN_Rand_H__

#include "SV_General.h"

class GN_Rand {

private :

	//by thilo: made all this previously non-class-members to classmembers to not flood the global namespace with new symbnols...

	/*
	* For each of the currently supported random number generators, we have a
	* break value on the amount of state information (you need at least this
	* many bytes of state info to support this random number generator), a degree
	* for the polynomial (actually a trinomial) that the R.N.G. is based on, and
	* the separation between the two lower order coefficients of the trinomial.
	*/
	#define	TYPE_0		0		/* linear congruential */
	#define	BREAK_0		8
	#define	DEG_0		0
	#define	SEP_0		0

	#define	TYPE_1		1		/* x**7 + x**3 + 1 */
	#define	BREAK_1		32
	#define	DEG_1		7
	#define	SEP_1		3

	#define	TYPE_2		2		/* x**15 + x + 1 */
	#define	BREAK_2		64
	#define	DEG_2		15
	#define	SEP_2		1

	#define	TYPE_3		3		/* x**31 + x**3 + 1 */
	#define	BREAK_3		128
	#define	DEG_3		31
	#define	SEP_3		3

	#define	TYPE_4		4		/* x**63 + x + 1 */
	#define	BREAK_4		256
	#define	DEG_4		63
	#define	SEP_4		1

	/*
	* Array versions of the above information to make code run faster --
	* relies on fact that TYPE_i == i.
	*/
	#define	MAX_TYPES	5		/* max number of types above */

	static char sccsid[];
	static int degrees[MAX_TYPES];
	static int seps [MAX_TYPES];
	static long randtbl[DEG_3 + 1];
	static long *fptr;
	static long *rptr;
	static long *state;
	static int rand_type;
	static int rand_deg;
	static int rand_sep;
	static long *end_ptr;

	//end by thilo

	char * initstate(int seed, char *arg_state, int n);
	char * setstate(char *arg_state); //by thilo: moved delcaration into the class to avoid error on machines having differing declarations of this method in gobal scope

public :

	//------------------------------- 
	// constructor/destructor
	//------------------------------- 
	GN_Rand();
	virtual ~GN_Rand();

	//------------------------------- 
	// public methods
	//------------------------------- 
	long getmax(void) {return (0x7fffffff);}    // maximum of uniform rand
	void srandom(int Seed);                     // set random seed
	long random(void);                          // uniform distribution
    double rand_gaus(double Mean, double Std); // gaus distribution

};   // class GN_Rand

#endif   // GN_Rand
