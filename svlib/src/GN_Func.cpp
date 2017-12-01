//************************************************************************
//    
//
//    Author  : Jialong HE
//    Date    : May 18, 1999
//************************************************************************
#include <stdlib.h>
#include <math.h>
#include "SV_Error.h"
#include "GN_Func.h"

static char SV_LibID[] = "Copyright (c) by Jialong He";
//===========================================================
//  Constructor
//===========================================================
GN_Func::GN_Func (){



}


//===========================================================
//  Destructor
//===========================================================
GN_Func::~GN_Func (){



}

//===========================================================
//	  
//   Obtain the machine EPSILON
//  i.e. the smallest positive number which, been added to 1., 
//  yields the result other than 1.
//
//===========================================================
double GN_Func::Epsilon(void) {

	double eps;
	eps = 1.0;
	while( eps + 1.0 != 1.0 ) {
		eps /= 2.0;
	}
	
	eps *= 2.0;

	return (eps);
}


//===========================================================
//  find a root for func(x) within the range (x1, x2)
//===========================================================
#define ITMAX 100
double GN_Func::Fzero(double (*func)(double), double x1, double x2, double tol) {

 
	int iter;
	double a=x1,b=x2,c,d,e,min1,min2;
	double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;
	
	double EPS = Epsilon();
	
	if (fb*fa > 0.0) { 
		REPORT_ERROR(SVLIB_BadArg, "Root must be bracketed!");
    }

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if (fb*fc > 0.0) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb=(*func)(b);
	}

	REPORT_ERROR(SVLIB_Fail, "Maximum number of iterations exceeded");
	return (-1);  // not reached
}


//===========================================================
//  find a min value for func(x) within range (x1, x2)
//===========================================================
/*
 ************************************************************************
 *	    		    C math library
 * function FMINBR - one-dimensional search for a function minimum
 *			  over the given range
 *
 * Input
 *	double fmin(a,b,f,tol)
 *	double a; 			Minimum will be seeked for over
 *	double b;  			a range [a,b], a being < b.
 *	double (*f)(double x);		Name of the function whose minimum
 *					will be seeked for
 *	double tol;			Acceptable tolerance for the minimum
 *					location. It have to be positive
 *					(e.g. may be specified as EPSILON)
 *
 * Output
 *	Fminbr returns an estimate for the minimum location with accuracy
 *	3*SQRT_EPSILON*abs(x) + tol.
 *	The function always obtains a local minimum which coincides with
 *	the global one only if a function under investigation being
 *	unimodular.
 *	If a function being examined possesses no local minimum within
 *	the given range, Fminbr returns 'a' (if f(a) < f(b)), otherwise
 *	it returns the right range boundary value b.
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.202 of the Russian edition
 *
 *	The function makes use of the "gold section" procedure combined with
 *	the parabolic interpolation.
 *	At every step program operates three abscissae - x,v, and w.
 *	x - the last and the best approximation to the minimum location,
 *	    i.e. f(x) <= f(a) or/and f(x) <= f(b)
 *	    (if the function f has a local minimum in (a,b), then the both
 *	    conditions are fulfiled after one or two steps).
 *	v,w are previous approximations to the minimum location. They may
 *	coincide with a, b, or x (although the algorithm tries to make all
 *	u, v, and w distinct). Points x, v, and w are used to construct
 *	interpolating parabola whose minimum will be treated as a new
 *	approximation to the minimum location if the former falls within
 *	[a,b] and reduces the range enveloping minimum more efficient than
 *	the gold section procedure. 
 *	When f(x) has a second derivative positive at the minimum location
 *	(not coinciding with a or b) the procedure converges superlinearly
 *	at a rate order about 1.324
 *
 ************************************************************************
 */
#define SQRT_EPSILON	1.49012e-08
#define EPSILON	        2.22045e-16

double GN_Func::Fmin(double (*f)(double x), double a, double b, double tol) {       	        /* An estimate to the min location*/ 

  double x,v,w;				/* Abscissae, descr. see above	*/
  double fx;				/* f(x)				*/
  double fv;				/* f(v)				*/
  double fw;				/* f(w)				*/
  const double r = (3.-sqrt(5.0))/2;	/* Gold section ratio		*/

  v = a + r*(b-a);  fv = (*f)(v);       /* First step - always gold section*/
  x = v;  w = v;
  fx=fv;  fw=fv;

  for(;;)		/* Main iteration loop	*/
  {
    double range = b-a;			/* Range over which the minimum */
					/* is seeked for		*/
    double middle_range = (a+b)/2;
    double tol_act =			/* Actual tolerance		*/
		SQRT_EPSILON*fabs(x) + tol/3;
    double new_step;      		/* Step at this iteration       */

       

    if( fabs(x-middle_range) + range/2 <= 2*tol_act )
      return x;				/* Acceptable approx. is found	*/

					/* Obtain the gold section step	*/
    new_step = r * ( x<middle_range ? b-x : a-x );


    			/* Decide if the interpolation can be tried	*/
    if( fabs(x-w) >= tol_act  )		/* If x and w are distinct      */
    {					/* interpolatiom may be tried	*/
	double p; 		/* Interpolation step is calcula-*/
	double q;              /* ted as p/q; division operation*/
                                        /* is delayed until last moment	*/
	double t;

	t = (x-w) * (fx-fv);
	q = (x-v) * (fx-fw);
	p = (x-v)*q - (x-w)*t;
	q = 2*(q-t);

	if( q>(double)0 )		/* q was calculated with the op-*/
	  p = -p;			/* posite sign; make q positive	*/
	else				/* and assign possible minus to	*/
	  q = -q;			/* p				*/

	if( fabs(p) < fabs(new_step*q) &&	/* If x+p/q falls in [a,b]*/
	    p > q*(a-x+2*tol_act) &&		/* not too close to a and */
	    p < q*(b-x-2*tol_act)  )            /* b, and isn't too large */
	  new_step = p/q;			/* it is accepted         */
					/* If p/q is too large then the	*/
					/* gold section procedure can 	*/
					/* reduce [a,b] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act )	/* Adjust the step to be not less*/
      if( new_step > (double)0 )	/* than tolerance		*/
	new_step = tol_act;
      else
	new_step = -tol_act;

				/* Obtain the next approximation to min	*/
    {				/* and reduce the enveloping range	*/
      double t = x + new_step;	/* Tentative point for the min	*/
      double ft = (*f)(t);
      if( ft <= fx )
      {                                 /* t is a better approximation	*/
	if( t < x )			/* Reduce the range so that	*/
	  b = x;                        /* t would fall within it	*/
	else
	  a = x;
      
	v = w;  w = x;  x = t;		/* Assign the best approx to x	*/
	fv=fw;  fw=fx;  fx=ft;
      }
      else                              /* x remains the better approx  */
      {        		             
	if( t < x )			/* Reduce the range enclosing x	*/
	  a = t;                   
	else
	  b = t;
      
        if( ft <= fw || w==x )
        {
           v = w;  w = t;
	   fv=fw;  fw=ft;
        }
        else if( ft<=fv || v==x || v==w )
        {
           v = t;
	   fv=ft;
        }
      }
      
    }			/* ----- end-of-block ----- */
  }		/* ===== End of loop ===== */

}

//=======================================================
//
//
//=======================================================
#define EPS 1.0e-5
#define JMAX 20
#define FUNC(x) ((*func)(x))
#define JMAXP JMAX+1
#define K 5

double trapzd( double (*func)(double), double a, double b, int n) {

	double x,tnm,sum,del;
	static double s;
	int it, j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1, j=1; j<n-1; j++) {it<<=1;}
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {sum += FUNC(x);}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}


double qtrap(double (*func)(double), double a, double b) {

	int j;
	double s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}

	REPORT_ERROR(SVLIB_Fail, "Too many steps");
    return(0.0);

}

/*---------------------------------------------------------*/

double qsimp(double (*func)(double), double a, double b) {

	int j;
	double s,st,ost,os;

	ost = os =  -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}

	REPORT_ERROR(SVLIB_Fail, "Too many steps");
    return(s);
}

/*---------------------------------------------------------*/

void polint(double *xa, double *ya, int n, double x, double *y, double *dy) {

	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	MArray_1D(c, n+1, double, "poling:c");
	MArray_1D(d, n+1, double, "poling:d");

	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
				REPORT_ERROR(SVLIB_Fail, "Error in routine POLINT");
			}

			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}

	MFree_1D(d);
	MFree_1D(c);
}


//=======================================================
//
//  Numerical integration for func(x) from a to b 
//=======================================================
double GN_Func::Finte(double (*func) (double) , double a, double b) {

	double ss, dss;
	double s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
			if (fabs(dss) < EPS*fabs(ss)) {
				return (ss);
			}
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}

	REPORT_ERROR(SVLIB_Fail, "Too many steps");
    return(ss);

}

/*---------------------------------------------------------*/
/*  Calculating error function erf(x), erfc(x) and inverse */
/*  error function x=erfinv(y)                             */
/*                                                         */
/*  error function is defined as                           */
/*                                                         */
/*                                                         */
/*                           x                             */
/*                           -                             */
/*                 2         | |          2                */
/*   erf(x)  =  --------     |    exp( - t  ) dt.          */
/*              sqrt(pi)   | |                             */
/*                          -                              */
/*                           0                             */
/*                                                         */
/*                                                         */
/*                           inf.                          */
/*                             -                           */
/*                  2         | |          2               */
/*   erfc(x)  =  --------     |    exp( - t  ) dt          */
/*               sqrt(pi)   | |                            */
/*                           -                             */
/*                            x                            */
/*                                                         */
/*   erfc(x) = 1 - erf(x)                                  */
/*   when  -Inf < x < Inf,  -1 < erf(x) < 1                */
/*                                                         */
/*                                                         */
/*   Author : Jialong He                                   */
/*   Date   : 10.11.97                                     */
/*                                                         */
/*---------------------------------------------------------*/
#define MAXLOG 7.08396418532264106224E2     /* log 2**1022 */
#define Pi     3.141592653589793

static double P[] = {
 2.46196981473530512524E-10,
 5.64189564831068821977E-1,
 7.46321056442269912687E0,
 4.86371970985681366614E1,
 1.96520832956077098242E2,
 5.26445194995477358631E2,
 9.34528527171957607540E2,
 1.02755188689515710272E3,
 5.57535335369399327526E2
};

static double Q[] = {
/* 1.00000000000000000000E0,*/
 1.32281951154744992508E1,
 8.67072140885989742329E1,
 3.54937778887819891062E2,
 9.75708501743205489753E2,
 1.82390916687909736289E3,
 2.24633760818710981792E3,
 1.65666309194161350182E3,
 5.57535340817727675546E2
};
static double R[] = {
 5.64189583547755073984E-1,
 1.27536670759978104416E0,
 5.01905042251180477414E0,
 6.16021097993053585195E0,
 7.40974269950448939160E0,
 2.97886665372100240670E0
};
static double S[] = {
/* 1.00000000000000000000E0,*/
 2.26052863220117276590E0,
 9.39603524938001434673E0,
 1.20489539808096656605E1,
 1.70814450747565897222E1,
 9.60896809063285878198E0,
 3.36907645100081516050E0
};
static double T[] = {
 9.60497373987051638749E0,
 9.00260197203842689217E1,
 2.23200534594684319226E3,
 7.00332514112805075473E3,
 5.55923013010394962768E4
};
static double U[] = {
/* 1.00000000000000000000E0,*/
 3.35617141647503099647E1,
 5.21357949780152679795E2,
 4.59432382970980127987E3,
 2.26290000613890934246E4,
 4.92673942608635921086E4
};

static double A[] = {0.886226899, -1.645349621,  0.914624893, -0.140543331};
static double B[] = {-2.118377725, 1.442710462, -0.329097515, 0.012229801};
static double C[] = {-1.970840454, -1.624906493, 3.429567803, 1.641345311};
static double D[] = {3.543889200, 1.637067800};

/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */

double polevl( double x, double coef[], int N ) {

double ans;
int i;
double *p;

p = coef;
ans = *p++;
i = N;

do
	ans = ans * x  +  *p++;
while( --i );

return( ans );
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evl(double x, double coef[], int N ) {

double ans;
double *p;
int i;

p = coef;
ans = x + *p++;
i = N-1;

do
	ans = ans * x  + *p++;
while( --i );

return( ans );
}


/*-----------------------------------*/
/* complementary error function      */
/*-----------------------------------*/
double GN_Func::erfc( double a) {

double p,q,x,y,z;

if( a < 0.0 )
	x = -a;
else
	x = a;

if( x < 1.0 )
	return( 1.0 - erf(a) );

z = -a * a;

if( z < -MAXLOG )
	{
under:
	REPORT_ERROR(SVLIB_BadArg, "erfc underflow!");
	if( a < 0 )
		return( 2.0 );
	else
		return( 0.0 );
	}

z = exp(z);

if( x < 8.0 )
	{
	p = polevl( x, P, 8 );
	q = p1evl( x, Q, 8 );
	}
else
	{
	p = polevl( x, R, 5 );
	q = p1evl( x, S, 6 );
	}
y = (z * p)/q;

if( a < 0 )
	y = 2.0 - y;

if( y == 0.0 )
	goto under;

return(y);
}

/*---------------------------------*/
/* Gaussian error function         */
/*---------------------------------*/
double GN_Func::erf(double x) {

double y, z;

if( fabs(x) > 1.0 )
	return( 1.0 - erfc(x) );
z = x * x;
y = x * polevl( z, T, 4 ) / p1evl( z, U, 5 );
return( y );

}

/*---------------------------------*/
/* Inverse Gaussian error function */
/*---------------------------------*/
double GN_Func::inverf(double Y) {

  double X, Z, Y0=0.7;

  if (fabs(Y)>1.0) {
	REPORT_ERROR(SVLIB_BadArg, "No inverse erf(x)");
  }

  if (fabs(Y)<Y0) {

    Z = Y * Y;
    X = Y * (((A[3]*Z+A[2])*Z+A[1])*Z+A[0]) / ((((B[3]*Z+B[2])*Z+B[1])*Z+B[0])*Z+1);

  }

  if (Y>Y0 && Y<1)  {

    Z = sqrt(-log((1.0-Y)/2.0));
    X = (((C[3]*Z+C[2])*Z+C[1])*Z+C[0]) / ((D[1]*Z+D[0])*Z+1);

  }

  if (Y<-Y0 && Y>-1)  {

    Z = sqrt(-log((1.0+Y)/2.0));
    X = -(((C[3]*Z+C[2])*Z+C[1])*Z+C[0]) / ((D[1]*Z+D[0])*Z+1);

  }
                                                             
  X = X - (erf(X) - Y) / (2/sqrt(Pi) * exp(-X*X));

  return(X);
}


//###########################################################
static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//------------------------------
// SVD function
//------------------------------
void SVD_NEW(double **a, int m, int n, double *w, double **v) {

	int flag,i,its,j,jj,k,l,nm;
	double c,f,h,s,x,y,z;
	double anorm=0.0,g=0.0,scale=0.0;
	double *rv1, *rv1_orig;

	if (m < n) {
		REPORT_ERROR(SVLIB_Fail, "SVD: You must augment A with extra zero rows");
	}

	MArray_1D(rv1, n+1, double, "rv1");
	rv1_orig = rv1;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=PYTHAG(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) {
				REPORT_ERROR(SVLIB_Fail, "No convergence in 30 SVD iterations");			
			};

			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

	MFree_1D(rv1_orig);
}

#undef SIGN
#undef MAX
#undef PYTHAG

//-----------------------------------------------
//
//-----------------------------------------------
void svbksb(double **u, double *w, double **v, int m, int n, double *b, double*x) {

	int jj,j,i;
	double s, *tmp;

	MArray_1D(tmp, n+1, double, "tmp");

	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	MFree_1D(tmp);
}

//-----------------------------------------------
//
//-----------------------------------------------
void svdfit(double *x, double *y, double *sig, int ndata, double *a, int ma, 
			double **u, double **v, double *w, double *chisq, void (*funcs)(double,double *,int)) {

	int j,i;
	double wmax,tmp,thresh,sum,*b,*afunc;

	MArray_1D(b, ndata+1, double, "b");
	MArray_1D(afunc, ma+1, double, "afunc");
	
	
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
		b[i]=y[i]*tmp;
	}

	SVD_NEW(u,ndata,ma,w,v);
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh= 2.22e-16 * wmax;
	for (j=1;j<=ma;j++)
		if (w[j] < thresh) w[j]=0.0;
	svbksb(u,w,v,ndata,ma,b,a);
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
	MFree_1D(afunc);
	MFree_1D(b);
}


#ifdef NEW_SVD_ALGORITHM
//=========================================
// Test SVD_NEW 
//=========================================
void main (void) {

     double **a, *w, **v;
     int row, col, m = 3, n = 3;

     MArray_2D(a, m+1, n+1, double, "a");
     MArray_2D(v, m+1, n+1, double, "v");
     MArray_1D(w, m+1, double, "w");
     a[1][1]=1; a[1][2]=3; a[1][3] = 5;
     a[2][1]=-4; a[2][2]=3; a[2][3] = 8;
     a[3][1]=0; a[3][2]=0; a[3][3] = 0;

/*

     a[0][0]=1; a[0][1]=3;
     a[1][0]=-4; a[1][1]=3;
 */

     SVD(a,m,n,w,v);

     printf ("Singular values:\n");
     for (row = 1; row <= m; row++ ) {
	printf ("%f ", w[row]);
     }
     printf ("\n");


     printf ("Matrix U:\n");
     for (row = 1; row <= m; row++ ) {
       for (col = 1; col <= n; col++ )
	 printf ("%f ", a[row][col]);
       printf ("\n");
     }

     printf ("Matrix V:\n");
     for (row = 1; row <= m; row++ ) {
       for (col = 1; col <= n; col++ )
	 printf ("%f ", v[row][col]);
       printf ("\n");
     }

}

//###########################################################
#endif


//===========================================================
// Function fitting
//  
// Given a set of data points X[1..ndata], Y[1..ndata]  
// with individual standard deviations W[1..ndata].
// find fitting coefficients Coef for the approximation
// Y = f(X) with
//
//  Y= Sum( Coef(i) * afunc(i, X) )
//
//  Funcs(x, afunc, NCoef), user provided functions.
//  return NCoef base fucntion values in "afunc"
//   
//===========================================================
double GN_Func::FitFuncs (double *X, double *Y, double *W, int ndata,
							double *Coef, int NCoef,
 						    void (*Funcs)(double a, double* afunc, int ma)) {

	
	double chisq;
	double **U, **V, *S;
	
	MArray_2D(U, ndata+2, NCoef+2, double, "FitFuncs");
	MArray_2D(V, NCoef+2, NCoef+2, double, "FitFuncs");
	MArray_1D(S, NCoef+2, double, "FitFuncs");
	
	//--------------------------------------------------
	// call function from "Numerical Recipes in C"
	//--------------------------------------------------
	svdfit(X, Y, W, ndata, Coef, NCoef, U, V, S, &chisq, Funcs);

	MFree_2D(U);
	MFree_2D(V);
	MFree_1D(S);
	return (chisq);

}

//========================================================
// find ploynomial roots
//========================================================
#define MAXM 100

typedef struct FCOMPLEX {double r,i;} fcomplex;

fcomplex Cadd(const fcomplex &a, fcomplex &b) {

	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(const fcomplex &a,const fcomplex &b) {

	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}

fcomplex Cmul(const fcomplex &a, fcomplex &b) {

	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(double re, double im) {

	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex &z) {

	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Cdiv(const fcomplex &a, fcomplex &b) {

	fcomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double Cabs(fcomplex &z) {

	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(const fcomplex &z) {

	fcomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(double x,const fcomplex &a) {

	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

//=======================================================
//
//
//=======================================================
void laguer(fcomplex *a, int m, fcomplex *x, double eps, int polish) {

	int j,iter;
	double err,cdx,abx;
	fcomplex sq,h,gp,gm,g2,g,b,d,dx,f,x1;

	for (iter=1;iter<=MAXM;iter++) {
		b=a[m];
		err=Cabs(b);
		d=f=Complex(0.0,0.0);
		abx=Cabs(*x);
		for (j=m-1;j>=0;j--) {
			f=Cadd(Cmul(*x,f),d);
			d=Cadd(Cmul(*x,d),b);
			b=Cadd(Cmul(*x,b),a[j]);
			err=Cabs(b)+abx*err;
		}
		err *= EPS;
		if (Cabs(b) <= err) return;
		g=Cdiv(d,b);
		g2=Cmul(g,g);
		h=Csub(g2,RCmul(2.0,Cdiv(f,b)));
		sq=Csqrt(RCmul((double) (m-1),Csub(RCmul((double) m,h),g2)));
		gp=Cadd(g,sq);
		gm=Csub(g,sq);
		if (Cabs(gp) < Cabs(gm))gp=gm;
		dx=Cdiv(Complex((double) m,0.0),gp);
		x1=Csub(*x,dx);
		if (x->r == x1.r && x->i == x1.i) return;
		*x=x1;
		cdx=Cabs(dx);
		if (!polish)
			if (cdx <= eps*Cabs(*x)) return;
	}
	REPORT_ERROR(SVLIB_Fail, "Too many iterations in routine LAGUER");
}

//=======================================================
// Given the degree m and the m+1 complex coefficients
// a[0..m] of polynomial sum( a(i) * x^i ), finds all
// roots in roots[1..m]
//
//=======================================================
void zroots(fcomplex *a, int m, fcomplex *roots, int polish) {

	int jj,j,i;
	fcomplex x,b,c,ad[MAXM];

	for (j=0;j<=m;j++) ad[j]=a[j];

	for (j=m;j>=1;j--) {
		x=Complex(0.0,0.0);
		laguer(ad,j,&x,EPS,0);
		if (fabs(x.i) <= (2.0*EPS*fabs(x.r))) x.i=0.0;
		roots[j]=x;
		b=ad[j];
		for (jj=j-1;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=Cadd(Cmul(x,b),c);
		}
	}
	if (polish)
		for (j=1;j<=m;j++)
			laguer(a,m,&roots[j],EPS,1);
	for (j=2;j<=m;j++) {
		x=roots[j];


		for (i=j-1;i>=1;i--) {
			if (roots[i].r <= x.r) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
}

#ifdef TEST_POLY_ROOT
//========================================
//  test find polynomials roots
//========================================
void main(void) {

    fcomplex a[10], roots[10];
    int i;

    a[0].r = 2; a[0].i = 2;
    a[1].r = 1; a[1].i = -2;
    a[2].r = -5; a[2].i = 0;
    a[3].r = 10; a[3].i = -2;
    a[4].r = -7; a[4].i = -2;
    a[5].r = -1; a[5].i = -3;
    a[6].r = -45; a[6].i = -1;


   zroots(a, 6, roots, 1); // a[0..m], roots[1..m]

   for (i=1; i<=6; i++)
    printf ("%f %f\n", roots[i].r, roots[i].i);

   printf ("\n\n");

}
#endif

//===========================================================
//  Find all roots of a polynomial
//
//  Polynomial is given sum (C(i) * x^i), i = 0..M;
//  M roots in Roots[0..M-1].
//
//  NOTE *C has M+1 elements from 0..M   
//===========================================================
void GN_Func::PolyRoots (COMPLEX *C, COMPLEX *Roots, int M) {

	
	fcomplex *a, *r;
	int Cnt;

	MArray_1D(a, M+1, fcomplex, "PolyRoots");
	MArray_1D(r, M+1, fcomplex, "PolyRoots");

	//--------------------------------------
	// copy poly coefficients
	//--------------------------------------
	for (Cnt=0; Cnt<=M; Cnt++) {
		a[Cnt].r = C[Cnt].real;
		a[Cnt].i = C[Cnt].imag;
	}

	//-----------------------------------------------
	// call function from "Numerical Recipes in C"
	//-----------------------------------------------
	zroots(a, M, r, 1);
	
	//--------------------------------------
	// copy roots 
	//--------------------------------------
	for (Cnt=0; Cnt<M; Cnt++) {
		Roots[Cnt].real = r[Cnt+1].r;
		Roots[Cnt].imag = r[Cnt+1].i;
	}
	
	MFree_1D(a);
	MFree_1D(r);
}

