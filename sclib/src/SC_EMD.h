/**************************************************************************/
/*    emd.h                                                               */
/*                                                                        */
/*    Last update: 3/24/98                                                */
/*                                                                        */
/*    An implementation of the Earth Movers Distance.                     */
/*    Based of the solution for the Transportation problem as described   */
/*    in "Introduction to Mathematical Programming" by F. S. Hillier and  */
/*    G. J. Lieberman, McGraw-Hill, 1990.                                 */
/*                                                                        */
/*    Copyright (C) 1998 Yossi Rubner                                     */
/*    Computer Science Department, Stanford University                    */
/*    E-Mail: rubner@cs.stanford.edu                                      */
/*    URL: http://vision.stanford.edu/~rubner                             */
/*                                                                        */
/*    This class is based on the above-mentioned C code by Yossi Rubner.  */
/*    It is nothing more than an object-orientaded wrapper around it with */
/*    the following changes:																							*/
/*     - The fixed static arrays of MAX_SIG_SIZE size have been	changed		*/
/*       to be allocated dynamically to reduce the memory-sage for large  */
/*			 signatures                                                       */
/*     - SC_Signature and SC_Centroid are now classes (no more structs)		*/
/*       and therefore have a method for ground-distance computation, so	*/
/*       the need for a function pointer to the distance measure has			*/
/*       vanished																													*/
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 13.06.2005																								*/
/**************************************************************************/

#ifndef __SC_EMD_H__
#define __SC_EMD_H__

#include "SC_Signature.h"
#include "SC_TweakableParameters.h"
#include "SC_Api.h"
#include <SV_Error.h>

/* DEFINITIONS */
//#define SCLIB_EMD_MAX_SIG_SIZE   150
//#define SCLIB_EMD_MAX_SIG_SIZE   100000
//#define SCLIB_EMD_MAX_ITERATIONS 1000
//#define SCLIB_EMD_INFINITY       1e20
//#define SCLIB_EMD_EPSILON        1e-6
//#define SCLIB_EMD_MAX_SIG_SIZE1  (SCLIB_EMD_MAX_SIG_SIZE+1)  /* FOR THE POSIBLE DUMMY FEATURE */
//#define SCLIB_EMD_DEBUG_LEVEL    0 
																	/* SCLIB_EMD_DEBUG_LEVEL:
                                      0 = NO MESSAGES
                                      1 = PRINT THE NUMBER OF ITERATIONS AND THE FINAL RESULT
                                      2 = PRINT THE RESULT AFTER EVERY ITERATION
                                      3 = PRINT ALSO THE FLOW AFTER EVERY ITERATION
                                      4 = PRINT A LOT OF INFORMATION (PROBABLY USEFUL ONLY FOR THE AUTHOR)
                                  */

/*****************************************************************************/
/* SC_Centroid SHOULD BE MODIFIED BY THE USER TO REFLECT THE FEATURE TYPE    */
//typedef int SC_Centroid; => has been done in SC_Signature/SC_Centroid class
/*****************************************************************************/

typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  float amount;         /* Amount of flow from "from" to "to" */
} SC_Flow;

class SC_EMD {
	private :

#ifdef SC_USE_EARTHMOVERSDISTANCE
    /* NEW TYPES DEFINITION */
    typedef struct node1_t {    /* node1_t IS USED FOR SINGLE-LINKED LISTS */
      int i;
      double val;
      struct node1_t *Next;
    } node1_t;
    typedef struct node2_t {    /* node2_t IS USED FOR DOUBLE-LINKED LISTS */
      int i, j;
      double val;
      struct node2_t *NextC;               /* NEXT COLUMN */
      struct node2_t *NextR;               /* NEXT ROW */
    } node2_t;

    /* GLOBAL VARIABLE DECLARATION */
    int _n1, _n2;                          /* SIGNATURES SIZES */
    //float _C[SCLIB_EMD_MAX_SIG_SIZE1][SCLIB_EMD_MAX_SIG_SIZE1];/* THE COST MATRIX */
    float **_C; /* THE COST MATRIX */
    //node2_t _X[SCLIB_EMD_MAX_SIG_SIZE1*2];            /* THE BASIC VARIABLES VECTOR */
    node2_t *_X;            /* THE BASIC VARIABLES VECTOR */
    
    /* VARIABLES TO HANDLE _X EFFICIENTLY */
    node2_t *_EndX, *_EnterX;
    //char _IsX[SCLIB_EMD_MAX_SIG_SIZE1][SCLIB_EMD_MAX_SIG_SIZE1];
    char **_IsX;
    //node2_t *_RowsX[SCLIB_EMD_MAX_SIG_SIZE1], *_ColsX[SCLIB_EMD_MAX_SIG_SIZE1];
    node2_t **_RowsX, **_ColsX;
    double _maxW;
    float _maxC;

    /* DECLARATION OF FUNCTIONS */
    float init(SC_Signature *Signature1, SC_Signature *Signature2);
    bool findBasicVariables(node1_t *U, node1_t *V); //by thilo: changed return value from void to bool to indicate error
    int isOptimal(node1_t *U, node1_t *V);
    int findLoop(node2_t **Loop);
    void newSol();
    void russel(double *S, double *D);
    void addBasicVariable(int minI, int minJ, double *S, double *D, node1_t *PrevUMinI, node1_t *PrevVMinJ, node1_t *UHead);
    void printSolution();
#endif

	protected :

		const double infinity; //instead of defined macro
		const double epsilon;
		SC_TweakableParameters *pTweak; //to hold some of the parameters that have previously been defined macros

  public :
		
		SC_EMD(SC_TweakableParameters *pTweak);
    virtual ~SC_EMD();

    float emd(SC_Signature *Signature1, SC_Signature *Signature2, SC_Flow *Flow, int *FlowSize);
};

#endif

