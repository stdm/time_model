#ifndef _VECTOR_H
#define _VECTOR_H
/**************************************************************************

VECTOR.H - GENERIC VECTOR MACROS INCLUDE FILE

Performs operations on vectors (arrays) of data.  Must give
pointers and type of data string to each macro.  A string
indicating the type of each vector is used to allocate a temp
pointer for use within each macro.  This allows the macros to
work with any legal c data type (int, float, double, unsigned
int, char, unsigned char, ect.)  More than one pointer can point
to the same vector, if desired.

More than one pointer can point to the same vector, if desired.
Thus, a = a*b can be obtained with MULT_VEC(a,b,a,10,int,int,int)
where a and b are 10 element int arrays.

Pointer arguments (a,b,c) can be expressions since they are used
only once.  Thus, ADD_VEC(z+2,x+5,y,10,int,int,int)
will work as expected on the 10 elements begining at z+2 and x+5
(note that x, y, and z are all int arrays of data).

CONTAINS THE FOLLOWING MACROS:

ADD_VEC(a,b,c,len,typea,typeb,typec)       add vectors, c = a + b
SUB_VEC(a,b,c,len,typea,typeb,typec)       subtract vectors, c = a - b
MULT_VEC(a,b,c,len,typea,typeb,typec)      multiply vectors, c = a * b
COPY_VEC(a,b,len,typea,typeb)              copy vector a to b,  b = a
DOTP_VEC(a,b,s,len,typea,typeb)            dot product s=sum(a[i]*b[i])
DIST_VEC(a,b,s,len,typea,typeb)            Euclidean distance s=sum((a-b) * (a-b))
SCALE_VEC(a,b,s,len,typea,typeb)           scale a vector b[i]=a[i] * s
ASSIGN_VEC(a,s,len,typea)                  assignment a[*] = s
SUM_VEC(a,s,len,typea)                     sum of vector, s=sum(a)
NORM_VEC(a,s,len,typea)                    sum of square vector, s=sum(a*a)
MAX_IND(a, s, len, typea)                  find index of maximum value
MIN_IND(a, s, len, typea)                  find index of minimum value

*************************************************************************/

/*
ADD_VEC macro:

ADDS TWO VECTORS (a,b) POINT BY POINT (PROMOTING AS REQUIRED) AND 
PUTS THE RESULT IN THE c VECTOR (DEMOTING IF REQUIRED).

ADD_VEC(a,b,c,len,typea,typeb,typec)

    a       pointer to first vector.
    b       pointer to second vector.
    c       pointer to result vector.
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.
    typec   legal C type describing the type of c data.
*/

#define ADD_VEC(a,b,c,len,typea,typeb,typec) {  \
                  typea *_PTA = a;  \
                  typeb *_PTB = b;  \
                  typec *_PTC = c;  \
                  int _IX;  \
                      for(_IX = 0 ; _IX < (len) ; _IX++)  \
                          *_PTC++ = (typec)((*_PTA++) + (*_PTB++));  \
                  }

/*
SUB_VEC macro:

SUBTRACTS TWO VECTORS (a,b) POINT BY POINT (PROMOTING AS REQUIRED) AND
PUTS THE RESULT IN THE c VECTOR (DEMOTING IF REQUIRED).

SUB_VEC(a,b,c,len,typea,typeb,typec)

    a       pointer to first vector.
    b       pointer to second vector.
    c       pointer to result vector.
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.
    typec   legal C type describing the type of c data.

*/

#define SUB_VEC(a,b,c,len,typea,typeb,typec) {  \
                  typea *_PTA = a;  \
                  typeb *_PTB = b;  \
                  typec *_PTC = c;  \
                  int _IX;  \
                      for(_IX = 0 ; _IX < (len) ; _IX++)  \
                          *_PTC++ = (typec)((*_PTA++) - (*_PTB++));  \
                  }

/*
MULT_VEC macro:

MULTIPLIES TWO VECTORS (a,b) POINT BY POINT (PROMOTING AS REQUIRED) AND
PUTS THE RESULT IN THE c VECTOR (DEMOTING IF REQUIRED).

MULT_VEC(a,b,c,len,typea,typeb,typec)

    a       pointer to first vector.
    b       pointer to second vector.
    c       pointer to result vector.
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.
    typec   legal C type describing the type of c data.

WARNING: The input data vectors are not cast to the type of c.
         This means that at least one of the input types must
         be able to represent the individual products without
         overflow.

*/

#define MULT_VEC(a,b,c,len,typea,typeb,typec) {  \
                   typea *_PTA = a;  \
                   typeb *_PTB = b;  \
                   typec *_PTC = c;  \
                   int _IX;  \
                       for(_IX = 0 ; _IX < (len) ; _IX++)  \
                           *_PTC++ = (typec)((*_PTA++) * (*_PTB++));  \
                   }


/*

COPY_VEC(a,b,len,typea,typeb)

    a       pointer to input vector.
    b       pointer to output vector.
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.

*/
#define COPY_VEC(a,b,len,typea,typeb) {  \
		       typea *_PTA = (typea *)a;  \
		       typeb *_PTB = (typeb *)b;  \
		       int _IX;  \
		       for(_IX = 0 ; _IX < (len) ; _IX++)  \
			   *_PTB++ = (typeb) (*_PTA++);  \
		    }


/*
DOTP_VEC macro:

FORMS THE SUM OF PRODUCTS OF TWO VECTORS (a,b) AND
PUTS THE RESULT IN THE PREVIOUSLY DEFINED VARIABLE s.

DOTP_VEC(a,b,s,len,typea,typeb)

    a       pointer to first vector.
    b       pointer to second vector.
    s       variable used to store result (not a pointer).
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.

WARNING: The input data vectors are not cast to the type of s.
         This means that at least one of the input types must
         be able to represent the individual products without
         overflow.
*/

#define DOTP_VEC(a,b,s,len,typea,typeb) {  \
                       typea *_PTA = a;  \
                       typeb *_PTB = b;  \
                       int _IX;  \
                       s = (*_PTA++) * (*_PTB++);  \
                       for(_IX = 1 ; _IX < (len) ; _IX++)  \
                           s += (*_PTA++) * (*_PTB++);  \
                   }

/*
SUM_VEC macro:

FORMS THE SUM THE VECTOR a AND PUT THE RESULT IN THE
PREVIOUSLY DEFINED VARIABLE s.

SUM_VEC(a,s,len,typea)

    a       pointer to first vector.
    s       variable used to store result (not a pointer).
    len     length of vector (integer).
    typea   legal C type describing the type of a data.

*/

#define SUM_VEC(a,s,len,typea) {  \
                       typea *_PTA = a;  \
                       int _IX;  \
                       s = (*_PTA++);  \
                       for(_IX = 1 ; _IX < (len) ; _IX++)  \
                           s += (*_PTA++);  \
		   }

/*

NORM_VEC macro:
FORMS THE SUM THE VECTOR a AND PUT THE RESULT IN THE
PREVIOUSLY DEFINED VARIABLE s.

NORM_VEC(a,s,len,typea)

    a       pointer to first vector.
    s       variable used to store result (not a pointer).
    len     length of vector (integer).
    typea   legal C type describing the type of a data.

*/

#define NORM_VEC(a,s,len,typea) {  \
		       typea *_PTA = a;  \
		       typea *_PTB = a;  \
		       int _IX;  \
		       s = (*_PTA++) * (*_PTB++);  \
		       for(_IX = 1 ; _IX < (len) ; _IX++) \
			   s += (*_PTA++) * (*_PTB++); \
		   }

/*
ASSIGN_VEC macro:

ASSIGN EACH ELEMENT OF VECTOR TO A CONSTANT VALUE s.

ASSIGN_VEC(a,s,len,typea)

    a       pointer to first vector.
    s       constant variable
    len     length of vector (integer).
    typea   legal C type describing the type of a data.

*/

#define ASSIGN_VEC(a,s,len,typea) {  \
		       typea *_PTA = a;  \
		       int _IX;  \
		       for(_IX = 0 ; _IX < (len) ; _IX++)  *_PTA++ = s;  \
		   }

/*

SCALE_VEC macro:

SCALES AND/OR CONVERTS (PROMOTES OR DEMOTES) THE INPUT VECTOR a
(of typea) AND COPIES THE SCALED VECTOR INTO ANOTHER VECTOR b
(of typeb).

SCALE_VEC(a,b,s,len,typea,typeb)

    a       pointer to input vector.
    b       pointer to output vector.
    s       variable used to scale output vector (not a pointer).
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.

*/

#define SCALE_VEC(a,b,s,len,typea,typeb) {  \
		       typea *_PTA = (typea *)a;  \
		       typeb *_PTB = (typeb *)b;  \
		       int _IX;  \
		       for(_IX = 0 ; _IX < (len) ; _IX++)  \
			   *_PTB++ = (typeb)(s * (*_PTA++));  \
		    }




/*

DIST_VEC macro:

Euclidean DISTANCE OF TWO VECTORS (a,b) AND
PUTS THE RESULT IN THE PREVIOUSLY DEFINED VARIABLE s.

DIST_VEC(a,b,s,len,typea,typeb)

    a       pointer to first vector.
    b       pointer to second vector.
    s       variable used to store result (not a pointer).
    len     length of vectors (integer).
    typea   legal C type describing the type of a data.
    typeb   legal C type describing the type of b data.

WARNING: The input data vectors are not cast to the type of s.
         This means that at least one of the input types must
         be able to represent the individual products without
         overflow.

*/

#define DIST_VEC(a,b,s,len,typea,typeb) {  \
		       typea *_PTA = a;  \
		       typeb *_PTB = b;  \
		       typea buf_a;  \
		       typeb buf_b;  \
		       int _IX;  \
		       buf_a = (*_PTA++) ;\
		       buf_b = (*_PTB++) ;\
		       s =  (buf_a - buf_b) * (buf_a - buf_b);  \
		       for (_IX = 1 ; _IX < (len) ; _IX++) {  \
			  buf_a = (*_PTA++) ;\
			  buf_b = (*_PTB++) ;\
			  s += (buf_a - buf_b) * (buf_a - buf_b);  \
		       } \
		   }



/*
MAX_IND macro:
find maximum value from array a, return index in s

MAX_IND(a,s,len,typea)

    a       pointer to first vector.
    s       integer variable for retruned index.
    len     length of vector (integer).
    typea   legal C type describing the type of a data.

*/

#define MAX_IND(a,s,len,typea) {  \
		       typea *_PTA = a, a_MAX;  \
		       int _IX;  \
		       a_MAX = (*_PTA++); s = 0; \
                       for(_IX = 1 ; _IX < (len) ; _IX++, _PTA++)  \
			 if ((*_PTA) > a_MAX) {a_MAX = (*_PTA); s = _IX;}; \
		   }



/*
MIN_IND macro:
find minmum value from array a, return index in s

MIN_IND(a,s,len,typea)

    a       pointer to first vector.
    s       integer variable for retruned index.
    len     length of vector (integer).
    typea   legal C type describing the type of a data.

*/

#define MIN_IND(a,s,len,typea) {  \
		       typea *_PTA = a, a_MIN;  \
		       int _IX;  \
		       a_MIN = (*_PTA++); s = 0; \
                       for(_IX = 1 ; _IX < (len) ; _IX++, _PTA++)  \
                           if ((*_PTA) < a_MIN) { a_MIN = (*_PTA); s = _IX;}; \
		   }

#endif

