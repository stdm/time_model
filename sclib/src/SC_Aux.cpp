/**************************************************************************/
/*    Some auxiliary functions needed by SC_*															*/
/*																																				*/
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 28.02.2004																								*/
/**************************************************************************/

#ifndef _MSC_VER
	#define _nextafter(x, y) nextafter((x), (y))
#endif

#include "SC_Aux.h"
#include <stdio.h>
#include <cstdlib> //for rand()
#include <string.h>
#include <assert.h>
#include <math.h>
#include <SV_Error.h>
#include <SV_DataIO.h> //for SV_Int* data types

//====================================================================================================================
//  returns a once for all initialized and seeded random number generator (for e.g. gaussian random deviates)
//====================================================================================================================
GN_Rand* sclib::getRandomizer(void) {
	static GN_Rand randomizer;

	return &randomizer;
}

//====================================================================================================================
//  the sigmoid function in its basic form has an s-shape, centered around 0 (where it has the value of 0.5), below
//  this it drops fastly towards zero, above it aproaching 1 (p is a scaling parameter that makes the curve steeper as 
//  it is increased):
//
//                  y ^
//                    |          
//                  1 |  __ _ _----------
//                    | /                                        1
//                    ||                           sig(x) = ------------    (inversion: 1-sig(x))
//                 .5 /                                      1 + e^(-p*x)
//                   ||
//  __________- - --/ |
//  -------------------------------------> 
//  -inf              0              inf x
//
//  this version lets one define the x-position of the infliction-point (pointFiveAt), the point where it shall reach 
//  1 (in fact, where it is >=1-eps) and where it shall reach zero (i.e. <=0+eps); if zeroAt>oneAt, the function is 
//  inversed such that it comes from 1 and approaches zero with rising x
//====================================================================================================================
double sclib::sigmoid(double x, double zeroAt, double pointFiveAt, double oneAt, double eps, double p) {
	const double default_e_inf = -13.81550955796377, default_e_0 = 13.81551055796427, defaultEps = 0.000001;
	double scaled_x, e_inf, e_0;
	double res;

	if (eps == defaultEps) {
		e_inf = default_e_inf; //e^(-e_inf) ~ infinity
		e_0 = default_e_0; //e^(-e_0) ~ 0
	} else {
		e_inf = -1.0 * log(1.0/eps);
		e_0 = -1.0 * log(eps);
	}

	if (zeroAt < oneAt) { //standard form as depicted above: the function comes from x/y=-inf/0 and aproaches x/y=inf/1
		if (x < zeroAt) {
			res = 0.0;
		} else if (x >= zeroAt && x < pointFiveAt) {
			scaled_x = scaleToInterval(x, zeroAt, pointFiveAt, e_inf, 0.0);
      res = sigmoid(scaled_x, p);
		} else if (x >= pointFiveAt && x < oneAt) {
			scaled_x = scaleToInterval(x, pointFiveAt, oneAt, 0.0, e_0);
			res = sigmoid(scaled_x, p);
		} else { // x >= oneAt
			res = 1.0;
		}
	} else { //inverted form: the function comes from x/y=-inf/1 and aproaches x/y=inf/0
		if (x < oneAt) {
			res = 1.0;
		} else if (x >= oneAt && x < pointFiveAt) {
			scaled_x = scaleToInterval(x, oneAt, pointFiveAt, e_inf, 0.0);
      res = 1.0 - sigmoid(scaled_x, p);
		} else if (x >= pointFiveAt && x < zeroAt) {
			scaled_x = scaleToInterval(x, pointFiveAt, zeroAt, 0.0, e_0);
			res = 1.0 - sigmoid(scaled_x, p);
		} else { // x >= zeroAt
			res = 0.0;
		}
	}

	return res;
}

//====================================================================================================================
//  this is the inverted version of the sigmoidal function: give a sigmoid value [0..1] and get the x value according
//  to the parameters: minAt describes the x-value for y=0, midAt the x-value for y=0.5, maxAt the value for
//  y=1
//====================================================================================================================
double sclib::invSigmoid(double y, double eps, double p) {
	double res;

	if (y <= 0.0+eps) {
		y = eps;
	} 
	if (y >= 1.0-eps) {
		y = 1.0-eps; 
	}

	res = sclib::ln((1/y) - 1);

	return res / (-1.0 * p);
}
double sclib::invSigmoid(double y, double minAt, double midAt, double maxAt, double eps, double p) {
	double res, minX, maxX, x, m, b;

	if (minAt > maxAt) { //inverted form: the function comes from x/y=-inf/1 and aproaches x/y=inf/0
		y = 1.0 - y; //then invert the sigmoid value
	}

	if (y <= 0.5) {
		if (p == 0.0) { //if p==0, a linar function is used
			m = 0.5 / (midAt-minAt);
			b = (0.5*minAt) / (midAt-minAt);
			res = sclib::getBetween(minAt, (y-b)/m, midAt);
		} else {
			minX = sclib::invSigmoid(0+eps, eps, 1.0);
			x = sclib::getBetween(minX, sclib::invSigmoid(y, eps, p), 0.0);
			res = sclib::scaleToInterval(x, minX, 0.0, minAt, midAt); 
		}
	} else {
		if (p == 0.0) { //if p==0, a linar function is useds
			m = 0.5 / (maxAt-midAt);
			b = 0.5 - ((0.5*midAt) / (maxAt-midAt));
			res = sclib::getBetween(midAt, (y-b)/m, maxAt);
		} else {
			maxX = sclib::invSigmoid(1-eps, eps, 1.0);
			x = sclib::getBetween(0.0, sclib::invSigmoid(y, eps, p), maxX);
			res = sclib::scaleToInterval(x, 0.0, maxX, midAt, maxAt);
		}
	}

	return res;
}

//====================================================================================================================
//  calculates the binomial coefficient, i.e.  /n\, with an iterative programm based on one by B.R.Preiss, see 
//                                             \k/
//  http://www.brpreiss.com/books/opus4/html/page467.html; the return value is a double to avoid overflow because the
//  numbers grow big quick here, the downside is: rounding off errors!
//====================================================================================================================
double sclib::binomi(unsigned long int n, unsigned long int k) {
	double *b, result;

	MArray_1D(b, n+1, double, "SC_Aux.binomi: b");
	b[0] = 1;

  //"The implementation shown uses an array of length n to represent a row of Pascal's triangle. Consequently, 
	//instead of a table of size O(n^2), the algorithm gets by with O(n) space. The implementation has been coded 
	//carefully so that the computation can be done in place. I.e., the elements of S_i+1 are computed in reverse 
	//so that they can be written over the elements of S_i that are no longer needed".
	//(S_i is the ith row of Pascal's triangle)
	for (unsigned long int i = 1; i <= n; i++) {
		b[i] = 1.0;
		for (unsigned long int j = i-1UL; j > 0; j--) {
			b[j] += b[j-1UL];
		}
	}

	result = b[k];
	MFree_1D(b);

	return result;
}

//====================================================================================================================
//  generates a random number between min and max
//====================================================================================================================
double sclib::rand(double min, double max) {
  double r;

	r = (double)(std::rand()) / RAND_MAX; //random number in the range 0..1
  r *= max - min; //random number in the range 0..(max-min)
  r += min; //random number in the range min..max
  
  return r;
}

//-------------------------------------------------------------------------------------------------------------------
//  generates a random index between 0 and maxIdx
//-------------------------------------------------------------------------------------------------------------------
unsigned int sclib::rand(unsigned int maxIdx) {
	double random;
	unsigned int r;

	random = (double)(std::rand()) / (double)(RAND_MAX); //random number in [0..1]
	random *= maxIdx; //random number in [0..maxIdx]
	r = sclib::round(random);

	return r;
}

//====================================================================================================================
//  generates n random numbers according to the given multivariate gaussian distribution (mean and covariance matrix)
//====================================================================================================================
double** sclib::randN(unsigned int n, unsigned int dim, double *mean, double **covar) {
	unsigned int t, d, dd;
	double **tmp, *sampleMean, **cholesky, **randomNumbers = NULL; 
	SC_MatrixFunctions matFunc;

	cholesky = matFunc.choleskyDecomposition(covar, dim); //sort of sqrt(covar)
	if (cholesky != NULL) {
		MArray_2D(tmp, (int)(n), (int)(dim), double, "SC_Aux.randN: tmp");
		sampleMean = matFunc.zeros(dim);
		for (t = 0; t < n; t++) {
			for (d = 0; d < dim; d++) {
				tmp[t][d] = sclib::getRandomizer()->rand_gaus(0.0, 1.0); //normally, i.e. N(0,1), distributed random numbers
				sampleMean[d] += tmp[t][d] / (double)(n); //gradually build up the mean of the pseudo-random numbers
			}
		}
		for (t = 0; t < n; t++) {
			for (d = 0; d < dim; d++) {
				tmp[t][d] -= sampleMean[d]; //subtract the mean => now they are really zero-mean!
			}
		}

		MArray_2D(randomNumbers, (int)(n), (int)(dim), double, "SC_Aux.randN: randomNumbers");
		for (t = 0; t < n; t++) {
			for (d = 0; d < dim; d++) {
				randomNumbers[t][d] = mean[d]; //this adds the desired mean => we have pseudo random numbers according to N(mean, covar)!
				for (dd = 0; dd < dim; dd++) {
					randomNumbers[t][d] += tmp[t][dd] * cholesky[d][dd]; //this is tmp * cholesky', which models the desired correlation (note the transpose on cholesky, because it is a lower triangular matrix and needs to be upper triangular here, see gaussian-faq by John D'Errico)
				}
			}
		}

		MFree_2D(tmp);
	} //cholesky != NULL

	return randomNumbers;
}

//====================================================================================================================
//  given an array with probabilities (i.e. single values between 0..1 that altogether sum up to 1), an index into 
//  this array is drawn according to the given distribution. this is useful e.g. for randomly choosing a mixture 
//  component of a gmm (or a state in a hmm) according to the mixture weights (state transition probabilities).
//====================================================================================================================
unsigned int sclib::drawIndexFromDistribution(double *weight, unsigned int weightCount) {
	unsigned int idx = 0;
	double randomValue, weightSum = 0.0;
  
	//randomly draw a index respecting the weights
	randomValue = sclib::rand(0.0, 1.0 - 0.000001);
  while (weightSum < randomValue) {
    weightSum += weight[idx++];
  }
  if (idx > 0) {
    idx--;
  }

	return idx;
}

//-------------------------------------------------------------------------------------------------------------------
//	checks (and returns an altered, definitively valid version) of the given fftSize; it must be a power of 2 and
//  greater or equal to the given signal length
//-------------------------------------------------------------------------------------------------------------------
int sclib::checkFftSize(int givenFftSize, int givenSignalLength) {
	double tmp = sclib::log2(givenFftSize);
	int validSize = sclib::round(pow(2.0, ceil(tmp)));

	while (validSize < givenSignalLength) {
		validSize *= 2;
	}

	return validSize;
}

//====================================================================================================================
//	returns the "next" (next representable) double that is greater than x
//====================================================================================================================
double sclib::incrementDouble(double x) {
	double next = x;
	
	next = _nextafter(next, std::numeric_limits<double>::max());

	/*
	if (fabs(next) <= 1.0) {
		next += std::numeric_limits<double>::epsilon()
	} else {
		//see: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
		//in short: floats are ordered lexicographically, that means, if you interpret their bit pattern as an equally sized integer
		//and increment it, and reinterpret it as a float, you get the "next" reachable float
		if (sizeof(double) == sizeof(SV_Int32)) {
			(*(SV_Int32*)&next) += 1;
		} else if (sizeof(double) == sizeof(SV_Int64)) {
			(*(SV_Int64*)&next) += 1;
		} else {
			REPORT_ERROR(SVLIB_Fail, "couldn't find integer of same size than double.");
		}
	}
	*/

	return next;
}

//====================================================================================================================
//	returns the "next" (next representable) double that is smaller than x
//====================================================================================================================
double sclib::decrementDouble(double x) {
	double previous = x;

	previous = _nextafter(previous, std::numeric_limits<double>::max()*-1.0);

	/*
	if (fabs(previous <= 1.0) {
		previous -= std::numeric_limits<double>::epsilon();
	} else {
		//see: http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
		//in short: floats are ordered lexicographically, that means, if you interpret their bit pattern as an equally sized integer
		//and decrement it, and reinterpret it as a float, you get the "previous" reachable float
		if (sizeof(double) == sizeof(SV_Int32)) {
			(*(SV_Int32*)&previous) -= 1;
		} else if (sizeof(double) == sizeof(SV_Int64)) {
			(*(SV_Int64*)&previous) -= 1;
		} else {
			REPORT_ERROR(SVLIB_Fail, "couldn't find integer of same size than double.");
		}
	}
	*/

	return previous;
}

//====================================================================================================================
//  given are the bounaries of two segments; returns amount of intersection if the two segments intersect
//====================================================================================================================
unsigned long int sclib::intersect(long int start1, long int end1, long int start2, long int end2) {
	unsigned long int res = 0;
	
	//find the minimum start-point of the two segments
	//if the two segments intersect, the greater start must lie between the smaller start and it's corresponding end
	if (start1 < start2) {
		if (start2 <= end1) {
			res = sclib::min(end1, end2) - start2 + 1;
		}
	} else {
		if (start1 <= end2) {
			res = sclib::min(end2, end1) - start1 + 1;
		}
	}

	return res;
}

//====================================================================================================================
// Returns the number of bits equal to "1" in the given integer; optimal for sparse ones; 
// See http://www-db.stanford.edu/~manku/bitcount/bitcount.html (2006-10-03) for more algorithms and explanations
//====================================================================================================================
unsigned int sclib::bitCount(unsigned long int n) {  
	int count = 0;

	while (n) {
		count++;
		n &= (n - 1);
	}

	return count;
}

//====================================================================================================================
//	read a complete, single line from an ascii-file
//====================================================================================================================
unsigned int sclib::readline(FILE* file, char* buffer, unsigned int maxlen) {
  unsigned int tmp;
  unsigned int counter = 0;

  do {
    tmp = fgetc(file);
    if (tmp != '\n' && tmp != EOF) {
      *(buffer+counter++) = (char) tmp;
    }
  } while (tmp != '\n' && tmp != EOF && counter<maxlen-1); //-1 so that the \0 below still fits into the maxLen-sized buffer

  *(buffer+counter++) = 0;
  return counter;
}

//====================================================================================================================
//	return the number of lines in the given file; 0 in case of error or empty file
//====================================================================================================================
unsigned long int sclib::countLines(const char* fileName) {
  unsigned long int res = 0;
  int tmp;
  FILE *file = NULL; 

  if ((file = fopen(fileName, "r")) != NULL) {
    do {
      tmp = fgetc(file);
      if (tmp == '\n' || tmp == EOF) {
        res++;
      }
    } while (tmp != EOF);
    fclose(file);
  }

  return res;
}

//====================================================================================================================
//	check if a character contains a number (didn't find this in the stdlib :-)
//====================================================================================================================
bool sclib::isNum(char test) {
  switch(test) {
		case '0' :	case '1' :	case '2' :	case '3' :	case '4' :	case '5' :	case '6' :	case '7' :	case '8' :
		case '9' : {return true; break;}
		default  : {return false; break;}
	}
}

//====================================================================================================================
//	cut the frontmost integer out of a string containing numbers separated by whitespaces
//====================================================================================================================
int sclib::getNextIntFromString(char* buffer, int size, const char *separators) {
  char* temp = new char[size+1];
  int i = SVLIB_Fail, posStart, posEnd, len = (int)(strnlen(buffer, size));

  posStart = (unsigned long)strcspn(buffer, "0123456789"); //strcspn() return the length of the initial segment that doesn't contain any numeral
	if (posStart == len) { //no integer found => empty the buffer
		sprintf(buffer, "");
	} else { //integer found, proceed
		//extract the integer
		posEnd = (unsigned long)strcspn((buffer+posStart), separators);
		strncpy(temp, (buffer + posStart), posEnd);
		temp[posEnd] = '\0';
		i = atoi(temp);

		//crop the integer out of the buffer
		sprintf(buffer, "%s", (buffer+posStart+posEnd));
	}

	MFree_1D(temp);

  return i;
}

//====================================================================================================================
//	extract the frontmost integer (containing numbers separated by whitespaces) out of a string starting at startPos, 
//  return last read position to provide a new starting-pos
//====================================================================================================================
int sclib::getNextIntFromString(const char* buffer, int size, int &res, int startPos, const char *separators) {
  char* temp = new char[size+1], *movedBuf = const_cast<char*>(buffer)+startPos;
  int i = SVLIB_Fail, posStart, posEnd, movedSize = size-startPos;
	int movedLen = (int)(strnlen(movedBuf, movedSize)), lastReadPos;

  posStart = (unsigned long)strcspn(movedBuf, "0123456789"); //strcspn() return the length of the initial segment that doesn't contain any numeral
	if (posStart == movedLen) { //no integer found
		lastReadPos = (int)(strnlen(buffer, size));
	} else { //integer found, proceed
		//extract the integer
		posEnd = (unsigned long)strcspn((movedBuf+posStart), separators);
		strncpy(temp, (movedBuf + posStart), posEnd);
		temp[posEnd] = '\0';
		res = atoi(temp);

		//crop the integer out of the buffer
		lastReadPos = startPos + posStart + posEnd + 1;
	}

	MFree_1D(temp);

  return lastReadPos;
}

//====================================================================================================================
//	cut the frontmost double out of a string containing numbers separated by whitespaces
//====================================================================================================================
double sclib::getNextDoubleFromString(char* buffer, int size, const char *separators) {
  char* temp = new char[size+1];
  int posStart, posEnd, len = (int)(strnlen(buffer, size));
	double d = -1.0;

  posStart = (unsigned long)strcspn(buffer, "0123456789."); //strcspn() return the length of the initial segment that doesn't contain any numeral
	if (posStart == len) { //no integer found => empty the buffer
		sprintf(buffer, "");
	} else { //integer found, proceed
		//extract the integer
		posEnd = (unsigned long)strcspn((buffer+posStart), separators);
		strncpy(temp, (buffer + posStart), posEnd);
		temp[posEnd] = '\0';
		d = atof(temp);

		//crop the integer out of the buffer
		sprintf(buffer, "%s", (buffer+posStart+posEnd));
	}

	MFree_1D(temp);

  return d;
}

//====================================================================================================================
//	extract the frontmost double (containing numbers separated by whitespaces) out of a string starting at startPos, 
//  return last read position to provide a new starting-pos
//====================================================================================================================
int sclib::getNextDoubleFromString(const char* buffer, int size, double &res, int startPos, const char * separators) {
  char* temp = new char[size+1], *movedBuf = const_cast<char*>(buffer)+startPos;
  int posStart, posEnd, movedSize = size-startPos;
	int movedLen = (int)(strnlen(movedBuf, movedSize)), lastReadPos;

  posStart = (unsigned long)strcspn(movedBuf, "0123456789."); //strcspn() return the length of the initial segment that doesn't contain any numeral
	if (posStart == movedLen) { //no integer found
		lastReadPos = (int)(strnlen(buffer, size));
	} else { //integer found, proceed
		//extract the integer
		posEnd = (unsigned long)strcspn((movedBuf+posStart), separators);
		strncpy(temp, (movedBuf + posStart), posEnd);
		temp[posEnd] = '\0';
		res = atof(temp);

		//crop the integer out of the buffer
		lastReadPos = startPos + posStart + posEnd + 1;
	}

	MFree_1D(temp);

  return lastReadPos;
}

//====================================================================================================================
//	cut the frontmost string out of a string containing strings separated by whitespaces
//====================================================================================================================
char* sclib::getNextStringFromString(char *buffer, int size, const char *separators) {
  char *temp = new char[size+1];
  int posStart = 0, posEnd;

	while ((buffer[posStart] <= 32) && (posStart < size)) {
	  posStart++;
	}
  posEnd = (unsigned long)strcspn((buffer + posStart), separators);

	//extract the first string
  strncpy(temp, (buffer + posStart), posEnd);
	temp[posEnd] = '\0';

	//crop it out
	sprintf(buffer, "%s", (buffer+posStart+posEnd));

  return temp;
}

//====================================================================================================================
//	returns the remainder of inLine after a '#'-sign, if any; otherwise returns NULL
//====================================================================================================================
char* sclib::extractComment(const char* inLine, int size, const char commentPrefix) {
  char* temp = new char[sclib::bufferSize];
  int posStart = (int)strcspn(inLine, &commentPrefix);

	if (inLine[posStart] == commentPrefix) {
		sprintf(temp, "%s", (inLine+posStart+1));
	} else {
		sprintf(temp, "");
	}

  return temp;
}

//====================================================================================================================
//  returns strIn without heading whitespaces (ASCII <= 32); a pointer to a new string is returned
//====================================================================================================================
char* sclib::lTrim(const char* strIn) {
	unsigned int x = 0;
	char* strOut;
	
	while ((strIn[x] <= 32) && (x < strlen(strIn))) {
	  x++;
	}

	if (x < (strlen(strIn)-1)) {
  	strOut = new char[strlen(strIn) - x + 1];
		sprintf(strOut, "%s\0", strIn + (x * sizeof(char)));

  } else {
		strOut = new char[1];
		sprintf(strOut, "\0");
	}

	return strOut;
}

//====================================================================================================================
//  returns strIn without trailing whitespaces (ASCII <= 32); a pointer to a new string is returned
//====================================================================================================================
char * sclib::rTrim(const char* strIn) {
  char *x = const_cast<char*>(strIn) + strlen(strIn);
  char *strOut;

  while ((*x <= 32) && (x >= strIn) && (*x > 0)) {
    x--;
  }
  
  strOut = new char[x - strIn + 1];
  strncpy(strOut, strIn, x-strIn);
  strOut[x - strIn] = '\0';

  return strOut;
}

//====================================================================================================================
//  returns strIn without heading and trailing whitespaces (ASCII <= 32); a pointer to a new string is returned
//====================================================================================================================
char* sclib::trim(const char* strIn) {
  char* temp = sclib::lTrim(strIn);
  char* strOut = sclib::rTrim(temp);

  MFree_1D(temp);
  
  return strOut;
}

//====================================================================================================================
//  returns strIn without heading and trailing whitespaces (ASCII <= 32); the original StrIn's content is changed and 
//  returned
//====================================================================================================================
char* trimInPlace(char* strIn, bool left, bool right) {
  char *tmp1, *tmp2;

	if (left == true) {
		tmp1 = sclib::lTrim(strIn);
		if (right == true) {
			tmp2 = sclib::rTrim(tmp1);
			sprintf(strIn, "%s", tmp2);
			MFree_1D(tmp1);
			MFree_1D(tmp2);
		} else {
			sprintf(strIn, "%s", tmp1);
			MFree_1D(tmp1);
		}
	} else if (right == true) {
		tmp1 = sclib::rTrim(strIn);
		sprintf(strIn, "%s", tmp1);
		MFree_1D(tmp1);
	}

	return strIn;
}

//====================================================================================================================
// Converts all lower-case lettes in 'text' to upper-case
//====================================================================================================================
char* sclib::uCase(char *text) {
	char *c = text;
	
	while (*c != '\0') {
		*c = toupper(*c);
		c++;
	}

	return text;
}

//====================================================================================================================
//	A simple 'like'-function written originally by Stefan Hogedal. It is supposed to mimic the behaviour of MS SQL 
//  Server 1.1
// 
//  The code comes from the newsgroup-thread "C++ has not a Like Operator?", 9 Apr. 1997 09:00, from 
//  microsoft.public.win32.programmer.gdi 
//  http://groups.google.de/group/microsoft.public.win32.programmer.gdi/browse_frm/thread/c1285fe694e5c900?page=end&q=like+operator+%22c%2B%2B%22&hl=de&
//  24.09.2005
//
//  Additions by thilo: Character-escaping within the likeExpr with the '\'-character, const parameters
//====================================================================================================================
bool sclib::like(const char *text, const char *likeExpr) { 
  char *start = const_cast<char*>(text); //by thilo
	char *textPtr = const_cast<char*>(text); //by thilo
	char *likeExprPtr = const_cast<char*>(likeExpr); //by thilo
  char *temp1, *temp2; //by thilo

  while (*likeExprPtr) { 
    switch (*likeExprPtr) { 
      case '\\': //by thilo
        if ((likeExprPtr > start) && (*(likeExprPtr-1) == '\\')) {
          if (*likeExprPtr++ != *textPtr++)
            return false;
          break; 
        }
        break;

      case '[': 
        if ((likeExprPtr > start) && (*(likeExprPtr-1) == '\\')) { //by thilo
          if (*likeExprPtr++ != *textPtr++)
            return false;
          break; 
        }
        if (strchr(likeExprPtr, ']') == NULL) 
          return false; 
        if (!*textPtr) 
          return false; 
        if ((*(++likeExprPtr) == '^') && ((likeExprPtr > start && *(likeExprPtr-1) != '\\') || likeExprPtr == start)) { //by thilo
          if (sclib::inSet(++likeExprPtr, *textPtr++)) 
            return false; 
        } else { 
          if (!sclib::inSet(likeExprPtr, *textPtr++)) 
            return false; 
        }

        //instead of doing this (which doesn't compile in the gcc)...
        //likeExpr = strchr(likeExpr, ']') + 1; 
        //if (likeExpr == (char*)(NULL+1)) //by thilo: maybe there was a closing bracket, but it was escaped and the real closing bracket for the set is missing, then the above line evaluates to NULL+1
        //  return false;

        //...we use this equivalent code, avoiding the expression "NULL + 1"
        likeExprPtr = strchr(likeExprPtr, ']');
        if (likeExprPtr == NULL) {
          return false;
        } else {
          likeExprPtr++;
        }
                
        temp1 = strchr(likeExprPtr, ']'); //by thilo: maybe there was an escaped ']' in the set, so we are now just 1 character behind this, not behind the closing bracket
        temp2 = strchr(likeExprPtr, '[');
        if ((temp1 != NULL && temp2 != NULL && temp1 < temp2) || (temp1 != NULL && temp2 == NULL))
          likeExprPtr = temp1 + 1;

        break; 
      
      case '%': 
        if ((likeExprPtr > start) && (*(likeExprPtr-1) == '\\')) { //by thilo
          if (*likeExprPtr++ != *textPtr++)
            return false;
          break; 
        }
        if (*(++likeExprPtr) == '\0') 
          return true; 
        while (*textPtr) { 
          if (sclib::like(textPtr, likeExprPtr)) {
            return true; 
          } else { 
            textPtr++; 
          } 
        } 
        return false; 
        /* 941124 break; */ 
      
      case '_':
        if ((likeExprPtr > start) && (*(likeExprPtr-1) == '\\')) { //by thilo
          if (*likeExprPtr++ != *textPtr++)
            return false;
          break; 
        }
        likeExprPtr++; 
        if (!*textPtr++) 
          return false; 
        break; 
      
      default: 
        if (*likeExprPtr++ != *textPtr++) 
          return false; 
        break; 
    } 
  }

  return *textPtr ? false : true; 
} 

//====================================================================================================================
//	An auxiliary function to the above simple 'like'-function written originally by Stefan Hogedal. 
// 
//  The code comes from the newsgroup-thread "C++ has not a Like Operator?", 9 Apr. 1997 09:00, from 
//  microsoft.public.win32.programmer.gdi 
//  http://groups.google.de/group/microsoft.public.win32.programmer.gdi/browse_frm/thread/c1285fe694e5c900?page=end&q=like+operator+%22c%2B%2B%22&hl=de&
//  24.09.2005
//====================================================================================================================
bool sclib::inSet(const char *likeExpr, char c) { 
  char* start = const_cast<char*>(likeExpr); //by thilo

  while (((*likeExpr != ']') || (likeExpr > start && *(likeExpr-1) == '\\')) && (*likeExpr != '\0')) { //by thilo
    if (*likeExpr == '\\') //by thilo
      *likeExpr++;

    if (*likeExpr == '-' && *(likeExpr-1) != '[' && *(likeExpr-1) != '^' && *(likeExpr+1) != ']' && *(likeExpr-1) != '\\') { //by thilo
      if (c >= *(likeExpr-1) && c <= *(likeExpr+1)) 
        return true; 
      else { 
        likeExpr += 2; 
        continue; 
      } 
    } 
    if( c == *likeExpr++ ) 
      return true; 
  } 

  return false;
} 

//====================================================================================================================
// Converts all upper-case lettes in 'text' to lower-case
//====================================================================================================================
char* sclib::lCase(char *text) {
	char *c = text;
	
	while (*c != '\0') {
		*c = tolower(*c);
		c++;
	}

	return text;
}

//====================================================================================================================
// Compares the two strings on max. len characters ignoring case; return values are the same as for strncmp()
//====================================================================================================================
int sclib::strincmp(const char *str1, const char *str2, int len) {
	int res;
	char *tmp1 = new char[len], *tmp2 = new char[len];

	strncpy(tmp1, str1, len);
	strncpy(tmp2, str2, len);

	res = strncmp(sclib::lCase(tmp1), sclib::lCase(tmp2), len);

	MFree_1D(tmp1);
	MFree_1D(tmp2);

	return res;
}

//====================================================================================================================
//  searches 'strIn' and replaces all occurences of 'in' with 'out'; changes the value of the argument, doesn't 
//  allocate new space or touch the pointer
//====================================================================================================================
void sclib::strReplace(const char *strIn, char in, char out) {
  char *x = const_cast<char*>(strIn) + strlen(strIn); //pointer to last character in strIn
  
  while (x >= strIn) { //move pointer from last to first character in strIn
    if (*x == in) {
      *x = out;
    }
    x--;
  }

  return;
}

//====================================================================================================================
//	this function takes the fileName and changes it's extension (if anyone, else: appends the new one) to the new one.
//  the space for the new fileName is allocated by this function and the pointer is returned, the 2 parameters remain
//  untouched. The new extension is meant to be given with the heading '.'
//====================================================================================================================
char* sclib::exchangeFileExtension(const char* fileName, const char* newExtension) {
  char* newFileName = NULL;
  char* lastDot = const_cast<char*>(strrchr(fileName, '.'));
  int length = 0;

  if (lastDot == NULL || *(lastDot+1) == '/') {
    length = (int)strlen(fileName);
  } else {
    length = (int)(lastDot - fileName);
  }

  newFileName = new char[length + strlen(newExtension) + 1];
  newFileName = strncpy(newFileName, fileName, length);
  strcpy(newFileName+length, newExtension);
  newFileName[length + strlen(newExtension)] = '\0'; 

  return newFileName;
}

//====================================================================================================================
//	returns a new string (space is allocated within, original pointer remains untouched) containing only a filename
//  without path information
//====================================================================================================================
char* sclib::extractFileName(const char* fullPath) {
  char* fileName = NULL;
  char* lastSlash = const_cast<char*>(strrchr(fullPath, '/'));
  int length = 0;

  if (lastSlash == NULL) {
    fileName = new char[strlen(fullPath) + 1];
    fileName = strcpy(fileName, fullPath);
  } else {
    length = (int)(fullPath + strlen(fullPath) - lastSlash);
    fileName = new char[length + 1];
    fileName = strcpy(fileName, lastSlash + 1);
  }

  return fileName;
}

//====================================================================================================================
//	returns a new string (space is allocated within, original pointer remains untouched) containg only a path
//  without trailing fileName; if giveLastSlash is true (the default), the last charcter will be the '/', so that a 
//  filename can immediately be added to this path and no information is lost by splitting a full path into its 
//  ingredients, filename and path.
//====================================================================================================================
char* sclib::extractPath(const char *fullPath, bool giveLastSlash) {
  char* path = NULL;
  char* lastSlash = const_cast<char*>(strrchr(fullPath, '/'));
  int length = 0;

  if (lastSlash == NULL) {
    path = new char;
    path = '\0';
  } else {
		length = (int)(lastSlash - fullPath + ((giveLastSlash==true)?1:0));
    path = new char[length+1];
    path = strncpy(path, fullPath, length);
    path[length] = '\0';
  }

  return path;
}

//====================================================================================================================
//	extract and return the extension of the given filename in a new string
//====================================================================================================================
char* sclib::extractExtension(const char *fileName) {
  char* extension = NULL;
  char* lastDot = const_cast<char*>(strrchr(fileName, '.'));
  unsigned int length;

  if (lastDot != NULL) {
    length = (unsigned int)(strlen(lastDot+1) + 1); //extension + \0
    extension = new char[length];
    extension = strncpy(extension, lastDot+1, length-1);
    extension[length-1] = '\0';
  }

  return extension;
}

//====================================================================================================================
//	add a postfix to a given filename by inserting it directly before the dot introducing the extension (or directly 
//  at the end, if there is no extension); a new string is returned
//====================================================================================================================
char* sclib::addPostfixToFilename(const char *fileName, const char *postfix) {
	char *newFilename = new char[strlen(fileName) + strlen(postfix) + 1];
	char *lastDot = const_cast<char*>(strrchr(fileName, '.'));
	unsigned int prefixLength, extensionLength, postfixLength;

	if (lastDot != NULL) {
		prefixLength = (unsigned int)(lastDot - fileName);
		extensionLength = (unsigned int)(strlen(lastDot));
		postfixLength = (unsigned int)(strlen(postfix));

		strncpy(newFilename, fileName, prefixLength);
		strncpy(newFilename+prefixLength, postfix, postfixLength);
		strncpy(newFilename+prefixLength+postfixLength, lastDot, extensionLength);
		newFilename[prefixLength+postfixLength+extensionLength] = '\0';
	} else { //fileName has no extension, so just concat the two strings
		sprintf(newFilename, "%s%s\0", fileName, postfix);
	}

	return newFilename;
}

//====================================================================================================================
//	takes something that should be a path and returns it in a from that can be used as a path throughout this 
//  library; this means: convert '\' and '\\' to '/' and adda trailing '/'
//====================================================================================================================
char* sclib::makePath(const char *pathCandidate) {
	char *path = new char[sclib::bufferSize];
	char *pp = path;
	const char *pc = pathCandidate;
	
	while (*pc != '\0') {
		if (*pc=='\\') { //replace backshlash with slash
			if (*(pc+1)!='\\') { //make double- to single slashes
				*pp = '/';
				pp++;
			}
		} else if (!(*pc=='/' && *(pc+1)=='/')) { //just copy everything else instead of double-slashes (make them single)
			*pp = *pc;
			pp++;
		}

		pc++;
	}
	
	if (*(pp-1) != '/') { //add trailing slash if absent
		*pp = '/';
		pp++;
	}
	*pp = '\0'; //add terminating null

	return path;
}

//====================================================================================================================
//	takes a string of the form "key=value" and return the key and value in newly allocated buffers
//====================================================================================================================
bool sclib::extractKeyValue(const char *pair, char* &key, char* &value) {
	const char *pos;

	MFree_1D(key);
	MFree_1D(value);
	
	pos = strchr(pair, '=');
	if (pos != NULL) { //pair has the form "key=value"
		key = new char[pos - pair + 1];
		strncpy(key, pair, pos-pair);
		key[pos-pair] = '\0';

		value = new char[pair + strlen(pair) - pos];
		sprintf(value, "%s\0", pos+1);
	} //form "key=value"

	return ((pos!=NULL) ? true : false);
}

//====================================================================================================================
//	checks if a given filename corresponds with a real existing file
//====================================================================================================================
bool sclib::fileExists(const char* fileName) {
	FILE*	testFile = NULL;
	bool	res = false;

	if (fileName != NULL && strncmp(fileName, "", sclib::bufferSize) != 0 && strncmp(fileName, "\0", sclib::bufferSize) != 0 && strlen(fileName) > 0) {
		testFile = fopen(fileName, "r");
		if (testFile != NULL) {
			res = true;
			fclose(testFile);
		}
	}

	return res;
}

//====================================================================================================================
//	checks if a given pathname corresponds with a real existing path which is accessible for writing
//====================================================================================================================
bool sclib::pathExists(const char* pathName) {
	FILE*	testFile = NULL;
	char* fileName = NULL;
	bool	res = false;

	if (strncmp(pathName, "", sclib::bufferSize) != 0 && strncmp(pathName, "\0", sclib::bufferSize) != 0 && strlen(pathName) > 0) {
		if (((((strrchr(pathName, '/') - pathName) / sizeof(char)) + sizeof(char)) < strlen(pathName)) != 0) {
			fileName = new char[strlen(pathName) + 6];
			sprintf(fileName, "%s/test\0", pathName);
		} else {
			fileName = new char[strlen(pathName) + 5];
			sprintf(fileName, "%stest\0", pathName);
		}

		testFile = fopen(fileName, "a");
		if (testFile != NULL) {
			res = true;
			fclose(testFile);
			remove(fileName);
		}
	}

  MFree_1D(fileName);
	return res;
}

//====================================================================================================================
//	return an empty string if it is a NULL pointer
//====================================================================================================================
const char* sclib::nullFilter(const char *text) {
	return (text!=NULL) ? text : "";
}

//====================================================================================================================
//	print a string to a file
//====================================================================================================================
void sclib::stringOut(const char *fileName, const char *string, SC_TweakableParameters *pTweak, const char* separator) {
	fstream	fileOut;
	char* fName;
	
  if (pTweak != NULL) {
	  fName = new char[strlen(pTweak->debug.debugDir) + strlen(fileName) + 1];
	  sprintf(fName, "%s%s\0", pTweak->debug.debugDir, fileName);
  } else {
	  fName = new char[strlen(fileName) + 1];
    sprintf(fName, "%s\0", fileName);
  }
	
  fileOut.open(fName, ios_base::out|ios_base::app);
	fileOut << string << separator;
	fileOut.close();

	MFree_1D(fName);

	return;
}

//====================================================================================================================
//	this function sets or retrieves the status of what the error handler does in release mode: throwing an exception
//  (as needed when the library is called via JNI) or exiting with error code
//====================================================================================================================
bool sclib::errorHandlerThrows(bool justGet, bool newValue) {
	static bool throws = false;

	if (justGet == false) {
		throws = newValue;
	}

	return throws;
}

//====================================================================================================================
//	new svlib-like error-handler, which allows resuming 
//====================================================================================================================
void sclib::errorHandler(int ErrorCode, const char* ErrorMsg, const char* FName, int LNum) {
  //ignore some strange errors
	if (ErrorCode == -3) {
		if (strcmp(ErrorMsg, "erfc underflow!") == 0) { //ignore erf underflows, they are handled otherwise (meaningful return-value of '0')
  		return;
		}
  } else if (ErrorCode == -1) {
		if (strcmp(ErrorMsg, "SweepCount > slimit") == 0) { //ignore this strange error which may arise during SVD computation in GN_Matrix
			return;
		}
  }

#ifdef NDEBUG
  fprintf(stderr, "\n%s (%d)\n[%s:%d]\nWill exit...", ErrorMsg, ErrorCode, FName, LNum);
	if (sclib::errorHandlerThrows(true) == true) {
		throw ErrorMsg;
	} else {
		exit(ErrorCode);
	}
#else
  //raise an assertion error so that debugging is possible in the IDE
  assert(ErrorCode == 0);

  char *c = new char[sclib::bufferSize];
  fprintf(stderr, "\n%s (%d)\n[%s:%d]\nPress '0' to exit, any key to resume.", ErrorMsg, ErrorCode, FName, LNum);
  fscanf(stdin, "%s", c);

	if (strcmp(c, "0") == 0) {
		MFree_1D(c);
		exit(ErrorCode);
  } else {
  	MFree_1D(c);
  }
#endif

	return;
}

//====================================================================================================================
// Generate synthetic data to test the models
//====================================================================================================================
void sclib::generateTestData(SV_Data* &pNoise, SV_Data* &pSpeech, unsigned long int D, unsigned long int J, unsigned long int I) {
	unsigned long int d, t, T;
  unsigned long int i, j;
	double mean[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8};
	double variance = 0.01;

	T = J * 500;
	pNoise = new SV_Data(T, D);
	for (t = 0; t < T; t++) {
		for (d = 0; d < D; d++) {
			//pNoise->Mat[t][d] = randomizer.rand_gaus(0.5, 0.01);
			j = (unsigned long)(floor(((float)(sclib::getRandomizer()->random()) / (float)(sclib::getRandomizer()->getmax())) * J));

			if (d == 0) {
				pNoise->Mat[t][d]	= (float)(sclib::getRandomizer()->rand_gaus(mean[j] * 10, sqrt(variance) * 10));
			} else {
				pNoise->Mat[t][d]	= (float)(sclib::getRandomizer()->rand_gaus(mean[j] * 10, sqrt(variance) * 10));
			}
		}
	}
	
	T = I * 500;
	pSpeech = new SV_Data(T, D);
	for (t = 0; t < T; t++) {
		for (d = 0; d < D; d++) {
			//pSpeech->Mat[t][d] = randomizer.rand_gaus(2.0, 2.0);
			i = (unsigned long)(floor(((float)(sclib::getRandomizer()->random()) / (float)(sclib::getRandomizer()->getmax())) * I));

			if (d == 0) {
				pSpeech->Mat[t][d]	= (float)(sclib::getRandomizer()->rand_gaus(mean[i] * 10, sqrt(variance) * 10));
			} else {
				pSpeech->Mat[t][d]	= (float)(sclib::getRandomizer()->rand_gaus(mean[i] * 10, sqrt(variance) * 10));
			}
		}
	}

	return;
}

//====================================================================================================================
// assuming that the currentValue is between 1 and maxValue, the percentage of how much of the maxValue is already 
// reached is printed on stdOut if there is at least minIncrease percentage points increase as compared to the 
// lastPercentage; set firstPrint=true if this is the first percentage to be printed in this position so that 
// indentation is handled correctly; the last printed percentage is returned (for the next lastPercentage parameter)
//====================================================================================================================
double sclib::printPercentage(unsigned long int maxValue, unsigned long int currentValue, double lastPercentage, double minIncrease, bool firstPrint) {
	double currentPercent, res;

	currentPercent = ((double)(currentValue) / (double)(maxValue)) * 100.0;

	if (currentPercent >= lastPercentage+minIncrease || firstPrint == true) {
		if (firstPrint == true) {
			printf("       ");
		}

		printf("\b\b\b\b\b\b\b%#5.1f%% ", currentPercent);

		res = currentPercent;
	} else {
		res = lastPercentage;
	}
	
	return res;
}
