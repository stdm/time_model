/**************************************************************************/
/*    This class is based on the libsvm-2.81-implementation of Support    */
/*    Vector Machines.                                                    */
/*    It is nothing more than an object-orientated wrapper around it.     */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.03.2006																								*/
/**************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include "SC_SVM.h"
#include "SC_Aux.h"

#ifdef SC_USE_LIBSVM
/*
//these are already defined in SC_Aux.h or the C++ standard library...
#ifndef min
  template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
#endif

#ifndef max
  template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
#endif

template <class T> inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
*/

template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
	dst = new T[n];
	memcpy((void *)dst,(void *)src,sizeof(T)*n);
}

#define INF HUGE_VAL
#define TAU 1e-12
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

/*
#if 1
  void info(char *fmt,...) {
	  va_list ap;
	  va_start(ap,fmt);
	  vprintf(fmt,ap);
	  va_end(ap);
  }
  void info_flush() {
	  fflush(stdout);
  }
#else
  void info(char *fmt,...) {}
  void info_flush() {}
#endif
*/
//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size limit in bytes
//

void SC_SVM::info(char *fmt,...) {
  if (this->verbose == true) {
    va_list ap;
    va_start(ap,fmt);
    vprintf(fmt,ap);
    va_end(ap);
  }
}
void SC_SVM::info_flush() {
  if (this->verbose == true) {
    fflush(stdout);
  }
}

void SC_SVM::Solver::info(char *fmt,...) {
  if (this->verbose == true) {
    va_list ap;
    va_start(ap,fmt);
    vprintf(fmt,ap);
    va_end(ap);
  }
}
void SC_SVM::Solver::info_flush() {
  if (this->verbose == true) {
    fflush(stdout);
  }
}

SC_SVM::Cache::Cache(int l_,int size_):l(l_),size(size_)
{
	head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
	size /= sizeof(Qfloat);
	size -= l * sizeof(head_t) / sizeof(Qfloat);
	size = sclib::max(size, 2*l);	// cache must be large enough for two columns
	lru_head.next = lru_head.prev = &lru_head;
}

SC_SVM::Cache::~Cache()
{
	for(head_t *h = lru_head.next; h != &lru_head; h=h->next)
		free(h->data);
	free(head);
}

void SC_SVM::Cache::lru_delete(head_t *h)
{
	// delete from current location
	h->prev->next = h->next;
	h->next->prev = h->prev;
}

void SC_SVM::Cache::lru_insert(head_t *h)
{
	// insert to last position
	h->next = &lru_head;
	h->prev = lru_head.prev;
	h->prev->next = h;
	h->next->prev = h;
}

int SC_SVM::Cache::get_data(const int index, Qfloat **data, int len)
{
	head_t *h = &head[index];
	if(h->len) lru_delete(h);
	int more = len - h->len;

	if(more > 0)
	{
		// free old space
		while(size < more)
		{
			head_t *old = lru_head.next;
			lru_delete(old);
			free(old->data);
			size += old->len;
			old->data = 0;
			old->len = 0;
		}

		// allocate new space
		h->data = (Qfloat *)realloc(h->data,sizeof(Qfloat)*len);
		size -= more;
		swap(h->len,len);
	}

	lru_insert(h);
	*data = h->data;
	return len;
}

void SC_SVM::Cache::swap_index(int i, int j)
{
	if(i==j) return;

	if(head[i].len) lru_delete(&head[i]);
	if(head[j].len) lru_delete(&head[j]);
	swap(head[i].data,head[j].data);
	swap(head[i].len,head[j].len);
	if(head[i].len) lru_insert(&head[i]);
	if(head[j].len) lru_insert(&head[j]);

	if(i>j) swap(i,j);
	for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
	{
		if(h->len > i)
		{
			if(h->len > j)
				swap(h->data[i],h->data[j]);
			else
			{
				// give up
				lru_delete(h);
				free(h->data);
				size += h->len;
				h->data = 0;
				h->len = 0;
			}
		}
	}
}

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//

void SC_SVM::Kernel::swap_index(int i, int j) const	// no so const...
{
	swap(x[i],x[j]);
	if(x_square) swap(x_square[i],x_square[j]);
}

SC_SVM::Kernel::Kernel(int l, SC_SVMnode * const * x_, const SC_SVMparameter& param) :kernel_type(param.kernel_type), degree(param.degree), gamma(param.gamma), coef0(param.coef0)
{
	switch(kernel_type)
	{
		case SCLIB_SVM_KERNEL_LINEAR:
			kernel_function = &Kernel::kernel_linear;
			break;
		case SCLIB_SVM_KERNEL_POLY:
			kernel_function = &Kernel::kernel_poly;
			break;
		case SCLIB_SVM_KERNEL_RBF:
			kernel_function = &Kernel::kernel_rbf;
			break;
		case SCLIB_SVM_KERNEL_SIGMOID:
			kernel_function = &Kernel::kernel_sigmoid;
			break;
	}

	clone(x,x_,l);

	if(kernel_type == SCLIB_SVM_KERNEL_RBF)
	{
		x_square = new double[l];
		for(int i=0;i<l;i++)
			x_square[i] = dot(x[i],x[i]);
	}
	else
		x_square = 0;
}

SC_SVM::Kernel::~Kernel()
{
	delete[] x;
	delete[] x_square;
}

double SC_SVM::Kernel::dot(const SC_SVMnode *px, const SC_SVMnode *py)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			sum += px->value * py->value;
			++px;
			++py;
		}
		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}			
	}
	return sum;
}

double SC_SVM::Kernel::k_function(const SC_SVMnode *x, const SC_SVMnode *y, const SC_SVMparameter& param)
{
	switch(param.kernel_type)
	{
		case SCLIB_SVM_KERNEL_LINEAR:
			return dot(x,y);
		case SCLIB_SVM_KERNEL_POLY:
			return pow(param.gamma*dot(x,y)+param.coef0,param.degree);
		case SCLIB_SVM_KERNEL_RBF:
		{
			double sum = 0;
			while(x->index != -1 && y->index !=-1)
			{
				if(x->index == y->index)
				{
					double d = x->value - y->value;
					sum += d*d;
					++x;
					++y;
				}
				else
				{
					if(x->index > y->index)
					{	
						sum += y->value * y->value;
						++y;
					}
					else
					{
						sum += x->value * x->value;
						++x;
					}
				}
			}

			while(x->index != -1)
			{
				sum += x->value * x->value;
				++x;
			}

			while(y->index != -1)
			{
				sum += y->value * y->value;
				++y;
			}
			
			return exp(-param.gamma*sum);
		}
		case SCLIB_SVM_KERNEL_SIGMOID:
			return tanh(param.gamma*dot(x,y)+param.coef0);
		default:
			return 0;	/* Unreachable */
	}
}

double SC_SVM::Kernel::kernel_linear(int i, int j) const
{
	return dot(x[i],x[j]);
}
double SC_SVM::Kernel::kernel_poly(int i, int j) const
{
	return pow(gamma*dot(x[i],x[j])+coef0,degree);
}
double SC_SVM::Kernel::kernel_rbf(int i, int j) const
{
	return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
}
double SC_SVM::Kernel::kernel_sigmoid(int i, int j) const
{
	return tanh(gamma*dot(x[i],x[j])+coef0);
}

// Generalized SMO+SVMlight algorithm
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + b^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//	Q, b, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping criterion
//
// solution will be put in \alpha, objective value will be put in obj
//

double SC_SVM::Solver::get_C(int i)
{
  //return (this->y[i] > 0)? this->Cp : this->Cn;
	return C[i]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
}

//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu): the following modification is for maintaining the alpha_status
void SC_SVM::Solver::update_alpha_status(int i)
{
	//if(alpha[i] >= get_C(i))
	//	alpha_status[i] = UPPER_BOUND;
	//else if(alpha[i] <= 0)
	//	alpha_status[i] = LOWER_BOUND;
	//else alpha_status[i] = FREE;
		alpha_status[i] = FREE;
		if (alpha[i] <= 0)
			alpha_status[i] |= LOWER_BOUND;
		if (alpha[i] >= get_C(i))
			alpha_status[i] |= UPPER_BOUND;
}

void SC_SVM::Solver::swap_index(int i, int j)
{
	Q->swap_index(i,j);
	swap(y[i],y[j]);
	swap(G[i],G[j]);
	swap(alpha_status[i],alpha_status[j]);
	swap(alpha[i],alpha[j]);
	swap(b[i],b[j]);
	swap(active_set[i],active_set[j]);
	swap(G_bar[i],G_bar[j]);
	swap(C[i], C[j]); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
}

void SC_SVM::Solver::reconstruct_gradient()
{
	// reconstruct inactive elements of G from G_bar and free variables

	if(active_size == l) return;

	int i;
	for(i=active_size;i<l;i++)
		G[i] = G_bar[i] + b[i];
	
	for(i=0;i<active_size;i++)
		if(is_free(i))
		{
			const Qfloat *Q_i = Q->get_Q(i,l);
			double alpha_i = alpha[i];
			for(int j=active_size;j<l;j++)
				G[j] += alpha_i * Q_i[j];
		}
}

//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
//void SC_SVM::Solver::Solve(int l, const QMatrix& Q, const double *b_, const schar *y_, double *alpha_, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking)
void SC_SVM::Solver::Solve(int l, const QMatrix& Q, const double *b_, const schar *y_, double *alpha_, const double* C_, double eps, SolutionInfo* si, int shrinking)

{
	this->l = l;
	this->Q = &Q;
	QD=Q.get_QD();
	clone(b, b_,l);
	clone(y, y_,l);
	clone(alpha,alpha_,l);
	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	//this->Cp = Cp;
	//this->Cn = Cn;
	clone(C, C_, l);
	this->eps = eps;
	unshrinked = false;

	// initialize alpha_status
	{
		alpha_status = new char[l];
		for(int i=0;i<l;i++)
			update_alpha_status(i);
	}

	// initialize active set (for shrinking)
	{
		active_set = new int[l];
		for(int i=0;i<l;i++)
			active_set[i] = i;
		active_size = l;
	}

	// initialize gradient
	{
		G = new double[l];
		G_bar = new double[l];
		int i;
		for(i=0;i<l;i++)
		{
			G[i] = b[i];
			G_bar[i] = 0;
		}
		for(i=0;i<l;i++)
			if(!is_lower_bound(i))
			{
				const Qfloat *Q_i = Q.get_Q(i,l);
				double alpha_i = alpha[i];
				int j;
				for(j=0;j<l;j++)
					G[j] += alpha_i*Q_i[j];
				if(is_upper_bound(i))
					for(j=0;j<l;j++)
						G_bar[j] += get_C(i) * Q_i[j];
			}
	}

	// optimization step

	int iter = 0;
	int counter = sclib::min(l,1000)+1;

	while(1)
	{
		// show progress and do shrinking

		if(--counter == 0)
		{
			counter = sclib::min(l,1000);
			if(shrinking) do_shrinking();
			info("."); info_flush();
		}

		int i,j;
		if(select_working_set(i,j)!=0)
		{
			// reconstruct the whole gradient
			reconstruct_gradient();
			// reset active set size and check
			active_size = l;
			info("*"); info_flush();
			if(select_working_set(i,j)!=0)
				break;
			else
				counter = 1;	// do shrinking next iteration
		}
		
		++iter;

		// update alpha[i] and alpha[j], handle bounds carefully
		
		const Qfloat *Q_i = Q.get_Q(i,active_size);
		const Qfloat *Q_j = Q.get_Q(j,active_size);

		double C_i = get_C(i);
		double C_j = get_C(j);

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];

		if(y[i]!=y[j])
		{
			double quad_coef = Q_i[i]+Q_j[j]+2*Q_i[j];
			if (quad_coef <= 0)
				quad_coef = TAU;
			double delta = (-G[i]-G[j])/quad_coef;
			double diff = alpha[i] - alpha[j];
			alpha[i] += delta;
			alpha[j] += delta;
			
			if(diff > 0)
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = diff;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = -diff;
				}
			}
			if(diff > C_i - C_j)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = C_i - diff;
				}
			}
			else
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = C_j + diff;
				}
			}
		}
		else
		{
			double quad_coef = Q_i[i]+Q_j[j]-2*Q_i[j];
			if (quad_coef <= 0)
				quad_coef = TAU;
			double delta = (G[i]-G[j])/quad_coef;
			double sum = alpha[i] + alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;

			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum - C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum - C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
		}

		// update G

		double delta_alpha_i = alpha[i] - old_alpha_i;
		double delta_alpha_j = alpha[j] - old_alpha_j;
		
		for(int k=0;k<active_size;k++)
		{
			G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
		}

		// update alpha_status and G_bar

		{
			bool ui = is_upper_bound(i);
			bool uj = is_upper_bound(j);
			update_alpha_status(i);
			update_alpha_status(j);
			int k;
			if(ui != is_upper_bound(i))
			{
				Q_i = Q.get_Q(i,l);
				if(ui)
					for(k=0;k<l;k++)
						G_bar[k] -= C_i * Q_i[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_i * Q_i[k];
			}

			if(uj != is_upper_bound(j))
			{
				Q_j = Q.get_Q(j,l);
				if(uj)
					for(k=0;k<l;k++)
						G_bar[k] -= C_j * Q_j[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_j * Q_j[k];
			}
		}
	}

	// calculate rho

	si->rho = calculate_rho();

	// calculate objective value
	{
		double v = 0;
		int i;
		for(i=0;i<l;i++)
			v += alpha[i] * (G[i] + b[i]);

		si->obj = v/2;
	}

	// put back the solution
	{
		for(int i=0;i<l;i++)
			alpha_[active_set[i]] = alpha[i];
	}

	// juggle everything back
	/*{
		for(int i=0;i<l;i++)
			while(active_set[i] != i)
				swap_index(i,active_set[i]);
				// or Q.swap_index(i,active_set[i]);
	}*/

	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	//si->upper_bound_p = Cp;
	//si->upper_bound_n = Cn;

	info("\noptimization finished, #iter = %d\n",iter);

	delete[] b;
	delete[] y;
	delete[] C; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	delete[] alpha;
	delete[] alpha_status;
	delete[] active_set;
	delete[] G;
	delete[] G_bar;
}

//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
void SC_SVM::Solver::Solve(int l, const QMatrix& Q, const double *b_, const schar *y_, double *alpha_, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking)
{
	double* C_ = new double[l];
	for (int i = 0; i < l; ++i)
		C_[i] = (y_[i] > 0 ? Cp : Cn);
	Solve(l, Q, b_, y_, alpha_, C_, eps, si, shrinking);
	si->upper_bound_p = Cp;
	si->upper_bound_n = Cn;
	delete[] C_;
}

// return 1 if already optimal, return 0 otherwise
int SC_SVM::Solver::select_working_set(int &out_i, int &out_j)
{
	// return i,j such that
	// i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
	// j: minimizes the decrease of obj value
	//    (if quadratic coefficeint <= 0, replace it with tau)
	//    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
	
	double Gmax = -INF;
	double Gmax2 = -INF;
	int Gmax_idx = -1;
	int Gmin_idx = -1;
	double obj_diff_min = INF;

	for(int t=0;t<active_size;t++)
		if(y[t]==+1)	
		{
			if(!is_upper_bound(t))
				if(-G[t] >= Gmax)
				{
					Gmax = -G[t];
					Gmax_idx = t;
				}
		}
		else
		{
			if(!is_lower_bound(t))
				if(G[t] >= Gmax)
				{
					Gmax = G[t];
					Gmax_idx = t;
				}
		}

	int i = Gmax_idx;
	const Qfloat *Q_i = NULL;
	if(i != -1) // NULL Q_i not accessed: Gmax=-INF if i=-1
		Q_i = Q->get_Q(i,active_size);

	for(int j=0;j<active_size;j++)
	{
		if(y[j]==+1)
		{
			if (!is_lower_bound(j))
			{
				double grad_diff=Gmax+G[j];
				if (G[j] >= Gmax2)
					Gmax2 = G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef=Q_i[i]+QD[j]-2*y[i]*Q_i[j];
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						Gmin_idx=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
		else
		{
			if (!is_upper_bound(j))
			{
				double grad_diff= Gmax-G[j];
				if (-G[j] >= Gmax2)
					Gmax2 = -G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef=Q_i[i]+QD[j]+2*y[i]*Q_i[j];
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						Gmin_idx=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
	}

	if(Gmax+Gmax2 < eps)
		return 1;

	out_i = Gmax_idx;
	out_j = Gmin_idx;
	return 0;
}

// return 1 if already optimal, return 0 otherwise
int SC_SVM::Solver::max_violating_pair(int &out_i, int &out_j)
{
	// return i,j: maximal violating pair

	double Gmax1 = -INF;		// max { -y_i * grad(f)_i | i in I_up(\alpha) }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;		// max { y_i * grad(f)_i | i in I_low(\alpha) }
	int Gmax2_idx = -1;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)	// y = +1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] >= Gmax1)
				{
					Gmax1 = -G[i];
					Gmax1_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] >= Gmax2)
				{
					Gmax2 = G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y = -1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] >= Gmax2)
				{
					Gmax2 = -G[i];
					Gmax2_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] >= Gmax1)
				{
					Gmax1 = G[i];
					Gmax1_idx = i;
				}
			}
		}
	}

	if(Gmax1+Gmax2 < eps)
 		return 1;

	out_i = Gmax1_idx;
	out_j = Gmax2_idx;
	return 0;
}

void SC_SVM::Solver::do_shrinking()
{
	int i,j,k;
	if(max_violating_pair(i,j)!=0) return;
	double Gm1 = -y[j]*G[j];
	double Gm2 = y[i]*G[i];

	// shrink
	
	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm1) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || -(Gm1 + Gm2) > eps*10) return;
	
	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm1) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double SC_SVM::Solver::calculate_rho()
{
	double r;
	int nr_free = 0;
	double ub = INF, lb = -INF, sum_free = 0;
	for(int i=0;i<active_size;i++)
	{
		double yG = y[i]*G[i];

		if(is_lower_bound(i))
		{
			if(y[i] > 0)
				ub = sclib::min(ub,yG);
			else
				lb = sclib::max(lb,yG);
		}
		else if(is_upper_bound(i))
		{
			if(y[i] < 0)
				ub = sclib::min(ub,yG);
			else
				lb = sclib::max(lb,yG);
		}
		else
		{
			++nr_free;
			sum_free += yG;
		}
	}

	if(nr_free>0)
		r = sum_free/nr_free;
	else
		r = (ub+lb)/2;

	return r;
}

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//

void SC_SVM::Solver_NU::Solve(int l, const QMatrix& Q, const double *b, const schar *y, double *alpha, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking)
{
	this->si = si;
	Solver::Solve(l,Q,b,y,alpha,Cp,Cn,eps,si,shrinking);
}


// return 1 if already optimal, return 0 otherwise
int SC_SVM::Solver_NU::select_working_set(int &out_i, int &out_j)
{
	// return i,j such that y_i = y_j and
	// i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
	// j: minimizes the decrease of obj value
	//    (if quadratic coefficeint <= 0, replace it with tau)
	//    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

	double Gmaxp = -INF;
	double Gmaxp2 = -INF;
	int Gmaxp_idx = -1;

	double Gmaxn = -INF;
	double Gmaxn2 = -INF;
	int Gmaxn_idx = -1;

	int Gmin_idx = -1;
	double obj_diff_min = INF;

	for(int t=0;t<active_size;t++)
		if(y[t]==+1)
		{
			if(!is_upper_bound(t))
				if(-G[t] >= Gmaxp)
				{
					Gmaxp = -G[t];
					Gmaxp_idx = t;
				}
		}
		else
		{
			if(!is_lower_bound(t))
				if(G[t] >= Gmaxn)
				{
					Gmaxn = G[t];
					Gmaxn_idx = t;
				}
		}

	int ip = Gmaxp_idx;
	int in = Gmaxn_idx;
	const Qfloat *Q_ip = NULL;
	const Qfloat *Q_in = NULL;
	if(ip != -1) // NULL Q_ip not accessed: Gmaxp=-INF if ip=-1
		Q_ip = Q->get_Q(ip,active_size);
	if(in != -1)
		Q_in = Q->get_Q(in,active_size);

	for(int j=0;j<active_size;j++)
	{
		if(y[j]==+1)
		{
			if (!is_lower_bound(j))	
			{
				double grad_diff=Gmaxp+G[j];
				if (G[j] >= Gmaxp2)
					Gmaxp2 = G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef = Q_ip[ip]+QD[j]-2*Q_ip[j];
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						Gmin_idx=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
		else
		{
			if (!is_upper_bound(j))
			{
				double grad_diff=Gmaxn-G[j];
				if (-G[j] >= Gmaxn2)
					Gmaxn2 = -G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef = Q_in[in]+QD[j]-2*Q_in[j];
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						Gmin_idx=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
	}

	if(sclib::max(Gmaxp+Gmaxp2,Gmaxn+Gmaxn2) < eps)
 		return 1;

	if (y[Gmin_idx] == +1)
		out_i = Gmaxp_idx;
	else
		out_i = Gmaxn_idx;
	out_j = Gmin_idx;

	return 0;
}

void SC_SVM::Solver_NU::do_shrinking()
{
	double Gmax1 = -INF;	// max { -y_i * grad(f)_i | y_i = +1, i in I_up(\alpha) }
	double Gmax2 = -INF;	// max { y_i * grad(f)_i | y_i = +1, i in I_low(\alpha) }
	double Gmax3 = -INF;	// max { -y_i * grad(f)_i | y_i = -1, i in I_up(\alpha) }
	double Gmax4 = -INF;	// max { y_i * grad(f)_i | y_i = -1, i in I_low(\alpha) }

	// find maximal violating pair first
	int k;
	for(k=0;k<active_size;k++)
	{
		if(!is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] > Gmax1) Gmax1 = -G[k];
			}
			else	if(-G[k] > Gmax3) Gmax3 = -G[k];
		}
		if(!is_lower_bound(k))
		{
			if(y[k]==+1)
			{	
				if(G[k] > Gmax2) Gmax2 = G[k];
			}
			else	if(G[k] > Gmax4) Gmax4 = G[k];
		}
	}

	// shrinking

	double Gm1 = -Gmax2;
	double Gm2 = -Gmax1;
	double Gm3 = -Gmax4;
	double Gm4 = -Gmax3;

	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm4) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || sclib::max(-(Gm1+Gm2),-(Gm3+Gm4)) > eps*10) return;
	
	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm4) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double SC_SVM::Solver_NU::calculate_rho()
{
	int nr_free1 = 0,nr_free2 = 0;
	double ub1 = INF, ub2 = INF;
	double lb1 = -INF, lb2 = -INF;
	double sum_free1 = 0, sum_free2 = 0;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)
		{
			if(is_lower_bound(i))
				ub1 = sclib::min(ub1,G[i]);
			else if(is_upper_bound(i))
				lb1 = sclib::max(lb1,G[i]);
			else
			{
				++nr_free1;
				sum_free1 += G[i];
			}
		}
		else
		{
			if(is_lower_bound(i))
				ub2 = sclib::min(ub2,G[i]);
			else if(is_upper_bound(i))
				lb2 = sclib::max(lb2,G[i]);
			else
			{
				++nr_free2;
				sum_free2 += G[i];
			}
		}
	}

	double r1,r2;
	if(nr_free1 > 0)
		r1 = sum_free1/nr_free1;
	else
		r1 = (ub1+lb1)/2;
	
	if(nr_free2 > 0)
		r2 = sum_free2/nr_free2;
	else
		r2 = (ub2+lb2)/2;
	
	si->r = (r1+r2)/2;
	return (r1-r2)/2;
}

//
// Q matrices for various formulations
//
SC_SVM::SVC_Q::SVC_Q(const SC_SVMproblem& prob, const SC_SVMparameter& param, const schar *y_):Kernel(prob.l, prob.x, param)
{
	clone(y,y_,prob.l);
	cache = new Cache(prob.l,(int)(param.cache_size*(1<<20)));
	QD = new Qfloat[prob.l];
	for(int i=0;i<prob.l;i++)
		QD[i]= (Qfloat)(this->*kernel_function)(i,i);
}
	
Qfloat* SC_SVM::SVC_Q::get_Q(int i, int len) const
{
	Qfloat *data;
	int start;
	if((start = cache->get_data(i,&data,len)) < len)
	{
		for(int j=start;j<len;j++)
			data[j] = (Qfloat)(y[i]*y[j]*(this->*kernel_function)(i,j));
	}
	return data;
}

Qfloat* SC_SVM::SVC_Q::get_QD() const
{
	return QD;
}

void SC_SVM::SVC_Q::swap_index(int i, int j) const
{
	cache->swap_index(i,j);
	Kernel::swap_index(i,j);
	swap(y[i],y[j]);
	swap(QD[i],QD[j]);
}

SC_SVM::SVC_Q::~SVC_Q()
{
	delete[] y;
	delete cache;
	delete[] QD;
}

SC_SVM::ONE_CLASS_Q::ONE_CLASS_Q(const SC_SVMproblem& prob, const SC_SVMparameter& param):Kernel(prob.l, prob.x, param)
{
	cache = new Cache(prob.l,(int)(param.cache_size*(1<<20)));
	QD = new Qfloat[prob.l];
	for(int i=0;i<prob.l;i++)
		QD[i]= (Qfloat)(this->*kernel_function)(i,i);
}

Qfloat* SC_SVM::ONE_CLASS_Q::get_Q(int i, int len) const
{
	Qfloat *data;
	int start;
	if((start = cache->get_data(i,&data,len)) < len)
	{
		for(int j=start;j<len;j++)
			data[j] = (Qfloat)(this->*kernel_function)(i,j);
	}
	return data;
}

Qfloat* SC_SVM::ONE_CLASS_Q::get_QD() const
{
	return QD;
}

void SC_SVM::ONE_CLASS_Q::swap_index(int i, int j) const
{
	cache->swap_index(i,j);
	Kernel::swap_index(i,j);
	swap(QD[i],QD[j]);
}

SC_SVM::ONE_CLASS_Q::~ONE_CLASS_Q()
{
	delete cache;
	delete[] QD;
}

SC_SVM::SVR_Q::SVR_Q(const SC_SVMproblem& prob, const SC_SVMparameter& param):Kernel(prob.l, prob.x, param)
{
	l = prob.l;
	cache = new Cache(l,(int)(param.cache_size*(1<<20)));
	QD = new Qfloat[2*l];
	sign = new schar[2*l];
	index = new int[2*l];
	for(int k=0;k<l;k++)
	{
		sign[k] = 1;
		sign[k+l] = -1;
		index[k] = k;
		index[k+l] = k;
		QD[k]= (Qfloat)(this->*kernel_function)(k,k);
		QD[k+l]=QD[k];
	}
	buffer[0] = new Qfloat[2*l];
	buffer[1] = new Qfloat[2*l];
	next_buffer = 0;
}

void SC_SVM::SVR_Q::swap_index(int i, int j) const
{
	swap(sign[i],sign[j]);
	swap(index[i],index[j]);
	swap(QD[i],QD[j]);
}

Qfloat* SC_SVM::SVR_Q::get_Q(int i, int len) const
{
	Qfloat *data;
	int real_i = index[i];
	if(cache->get_data(real_i,&data,l) < l)
	{
		for(int j=0;j<l;j++)
			data[j] = (Qfloat)(this->*kernel_function)(real_i,j);
	}

	// reorder and copy
	Qfloat *buf = buffer[next_buffer];
	next_buffer = 1 - next_buffer;
	schar si = sign[i];
	for(int j=0;j<len;j++)
		buf[j] = si * sign[j] * data[index[j]];
	return buf;
}

Qfloat* SC_SVM::SVR_Q::get_QD() const
{
	return QD;
}

SC_SVM::SVR_Q::~SVR_Q()
{
	delete cache;
	delete[] sign;
	delete[] index;
	delete[] buffer[0];
	delete[] buffer[1];
	delete[] QD;
}

//
// construct and solve various formulations
//
//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
//void SC_SVM::solve_c_svc(const SC_SVMproblem *prob, const SC_SVMparameter* param, double *alpha, Solver::SolutionInfo* si, double Cp, double Cn)
void SC_SVM::solve_c_svc(const SC_SVMproblem *prob, const SC_SVMparameter* param, double *alpha, Solver::SolutionInfo* si, const double *C_)
{
	int l = prob->l;
	double *minus_ones = new double[l];
	schar *y = new schar[l];

	int i;

	for(i=0;i<l;i++)
	{
		alpha[i] = 0;
		minus_ones[i] = -1;
		if(prob->y[i] > 0) y[i] = +1; else y[i]=-1;
	}

  SC_SVM::Solver s(this->verbose);
	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	//s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y, alpha, Cp, Cn, param->eps, si, param->shrinking);
	s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y, alpha, C_, param->eps, si, param->shrinking);

	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	//double sum_alpha=0;
	//for(i=0;i<l;i++)
	//	sum_alpha += alpha[i];
  //
	//if (Cp==Cn)
	//	info("nu = %f\n", sum_alpha/(Cp*prob->l));

	for(i=0;i<l;i++)
		alpha[i] *= y[i];

	delete[] minus_ones;
	delete[] y;
}

void SC_SVM::solve_nu_svc(const SC_SVMproblem *prob, const SC_SVMparameter *param,	double *alpha, Solver::SolutionInfo* si)
{
	int i;
	int l = prob->l;
	double nu = param->nu;

	schar *y = new schar[l];

	for(i=0;i<l;i++)
		if(prob->y[i]>0)
			y[i] = +1;
		else
			y[i] = -1;

	double sum_pos = nu*l/2;
	double sum_neg = nu*l/2;

	for(i=0;i<l;i++)
		if(y[i] == +1)
		{
			alpha[i] = sclib::min(1.0,sum_pos);
			sum_pos -= alpha[i];
		}
		else
		{
			alpha[i] = sclib::min(1.0,sum_neg);
			sum_neg -= alpha[i];
		}

	double *zeros = new double[l];

	for(i=0;i<l;i++)
		zeros[i] = 0;

  SC_SVM::Solver_NU s(this->verbose);
	s.Solve(l, SVC_Q(*prob,*param,y), zeros, y,	alpha, 1.0, 1.0, param->eps, si,  param->shrinking);
	double r = si->r;

	info("C = %f\n",1/r);

	for(i=0;i<l;i++)
		alpha[i] *= y[i]/r;

	si->rho /= r;
	si->obj /= (r*r);
	si->upper_bound_p = 1/r;
	si->upper_bound_n = 1/r;

	delete[] y;
	delete[] zeros;
}

void SC_SVM::solve_one_class(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, SC_SVM::Solver::SolutionInfo* si)
{
	int l = prob->l;
	double *zeros = new double[l];
	schar *ones = new schar[l];
	int i;

	int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound

	for(i=0;i<n;i++)
		alpha[i] = 1;
	if(n<prob->l)
		alpha[n] = param->nu * prob->l - n;
	for(i=n+1;i<l;i++)
		alpha[i] = 0;

	for(i=0;i<l;i++)
	{
		zeros[i] = 0;
		ones[i] = 1;
	}

  SC_SVM::Solver s(this->verbose);
	s.Solve(l, SC_SVM::ONE_CLASS_Q(*prob,*param), zeros, ones, alpha, 1.0, 1.0, param->eps, si, param->shrinking);

	delete[] zeros;
	delete[] ones;
}


//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
//void SC_SVM::solve_epsilon_svr(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, Solver::SolutionInfo* si)
void SC_SVM::solve_epsilon_svr(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, Solver::SolutionInfo* si, const double *C_)
{
	int l = prob->l;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	double *C2 = new double[2*l]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	schar *y = new schar[2*l];
	int i;

	for(i=0;i<l;i++)
	{
		alpha2[i] = 0;
		linear_term[i] = param->p - prob->y[i];
		y[i] = 1;
		C2[i] = C_[i]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)

		alpha2[i+l] = 0;
		linear_term[i+l] = param->p + prob->y[i];
		y[i+l] = -1;
		C2[i+l] = C_[i]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	}

  SC_SVM::Solver s(this->verbose);
	//s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y, alpha2, param->C, param->C, param->eps, si, param->shrinking);
	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,	alpha2, C2, param->eps, si, param->shrinking);

	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	//double sum_alpha = 0;
	//for(i=0;i<l;i++)
	//{
	//	alpha[i] = alpha2[i] - alpha2[i+l];
	//	sum_alpha += fabs(alpha[i]);
	//}
	//info("nu = %f\n",sum_alpha/(param->C*l));
	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];

	delete[] alpha2;
	delete[] linear_term;
	delete[] C2; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	delete[] y;
}

void SC_SVM::solve_nu_svr(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, SC_SVM::Solver::SolutionInfo* si)
{
	int l = prob->l;
	double C = param->C;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	schar *y = new schar[2*l];
	int i;

	double sum = C * param->nu * l / 2;
	for(i=0;i<l;i++)
	{
		alpha2[i] = alpha2[i+l] = sclib::min(sum,C);
		sum -= alpha2[i];

		linear_term[i] = - prob->y[i];
		y[i] = 1;

		linear_term[i+l] = prob->y[i];
		y[i+l] = -1;
	}

	SC_SVM::Solver_NU s(this->verbose);
	s.Solve(2*l, SC_SVM::SVR_Q(*prob,*param), linear_term, y, alpha2, C, C, param->eps, si, param->shrinking);

	info("epsilon = %f\n",-si->r);

	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];

	delete[] alpha2;
	delete[] linear_term;
	delete[] y;
}


SC_SVM::SC_SVMdecisionFunction SC_SVM::svm_train_one(const SC_SVMproblem *prob, const SC_SVMparameter *param, double Cp, double Cn)
{
	double *alpha = Malloc(double,prob->l);
	Solver::SolutionInfo si;
	double *C_ = NULL; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)

	//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu): complete if-else-else-block is new
	if (prob->W){
		C_ = new double[prob->l];
		for(int i=0;i<prob->l;i++) C_[i] = prob->W[i] * param->C;
	}
	else if (param->svm_type == SCLIB_SVM_TYPE_CSVC){
		C_ = new double[prob->l];
		for(int i=0;i<prob->l;i++) C_[i] = (prob->y[i] > 0 ? Cp : Cn);
	}
	else if (param->svm_type == SCLIB_SVM_TYPE_EPSILONSVR){
		C_ = new double[prob->l];
		for(int i=0;i<prob->l;i++) C_[i] = param->C;
	}

	switch(param->svm_type)
	{
		case SCLIB_SVM_TYPE_CSVC:
			//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			//solve_c_svc(prob,param,alpha,&si,Cp,Cn);
			solve_c_svc(prob,param,alpha,&si,C_);
			break;
		case SCLIB_SVM_TYPE_NUSVC:
			info("Warning!! weighted mode not supported for NU_SVC"); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			solve_nu_svc(prob,param,alpha,&si);
			break;
		case SCLIB_SVM_TYPE_ONECLASS:
			info("Warning!! weighted mode not supported for one-class SVM"); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			solve_one_class(prob,param,alpha,&si);
			break;
		case SCLIB_SVM_TYPE_EPSILONSVR:
			//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu
			//solve_epsilon_svr(prob,param,alpha,&si);
			solve_epsilon_svr(prob,param,alpha,&si,C_);
			break;
		case SCLIB_SVM_TYPE_NUSVR:
			info("Warning!! weighted mode not supported for NU_SVR"); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			solve_nu_svr(prob,param,alpha,&si);
			break;
	}

	info("obj = %f, rho = %f\n",si.obj,si.rho);

	// output SVs

	int nSV = 0;
	int nBSV = 0;
	for(int i=0;i<prob->l;i++)
	{
		if(fabs(alpha[i]) > 0)
		{
			++nSV;
			if(C_){ //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
				if(fabs(alpha[i]) >= C_[i]) //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
					++nBSV; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			} //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			else //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			if(prob->y[i] > 0)
			{
				if(fabs(alpha[i]) >= si.upper_bound_p)
					++nBSV;
			}
			else
			{
				if(fabs(alpha[i]) >= si.upper_bound_n)
					++nBSV;
			}
		}
	}

	if (C_) delete[] C_; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	info("nSV = %d, nBSV = %d\n",nSV,nBSV);

	SC_SVMdecisionFunction f;
	f.alpha = alpha;
	f.rho = si.rho;
	return f;
}

// Platt's binary SVM Probablistic Output: an improvement from Lin et al.
void SC_SVM::sigmoid_train(int l, const double *dec_values, const double *labels, double& A, double& B)
{
	double prior1=0, prior0 = 0;
	int i;

	for (i=0;i<l;i++)
		if (labels[i] > 0) prior1+=1; //TODO by thilo: bring in individual weights here?
		else prior0+=1;
	
	int max_iter=100; 	// Maximal number of iterations
	double min_step=1e-10;	// Minimal step taken in line search
	double sigma=1e-3;	// For numerically strict PD of Hessian
	double eps=1e-5;
	double hiTarget=(prior1+1.0)/(prior1+2.0);
	double loTarget=1/(prior0+2.0);
	double *t=Malloc(double,l);
	double fApB,p,q,h11,h22,h21,g1,g2,det,dA,dB,gd,stepsize;
	double newA,newB,newf,d1,d2;
	int iter; 
	
	// Initial Point and Initial Fun Value
	A=0.0; B=log((prior0+1.0)/(prior1+1.0));
	double fval = 0.0;

	for (i=0;i<l;i++)
	{
		if (labels[i]>0) t[i]=hiTarget;
		else t[i]=loTarget;
		fApB = dec_values[i]*A+B;
		if (fApB>=0)
			fval += t[i]*fApB + log(1+exp(-fApB));
		else
			fval += (t[i] - 1)*fApB +log(1+exp(fApB));
	}
	for (iter=0;iter<max_iter;iter++)
	{
		// Update Gradient and Hessian (use H' = H + sigma I)
		h11=sigma; // numerically ensures strict PD
		h22=sigma;
		h21=0.0;g1=0.0;g2=0.0;
		for (i=0;i<l;i++)
		{
			fApB = dec_values[i]*A+B;
			if (fApB >= 0)
			{
				p=exp(-fApB)/(1.0+exp(-fApB));
				q=1.0/(1.0+exp(-fApB));
			}
			else
			{
				p=1.0/(1.0+exp(fApB));
				q=exp(fApB)/(1.0+exp(fApB));
			}
			d2=p*q;
			h11+=dec_values[i]*dec_values[i]*d2;
			h22+=d2;
			h21+=dec_values[i]*d2;
			d1=t[i]-p;
			g1+=dec_values[i]*d1;
			g2+=d1;
		}

		// Stopping Criteria
		if (fabs(g1)<eps && fabs(g2)<eps)
			break;

		// Finding Newton direction: -inv(H') * g
		det=h11*h22-h21*h21;
		dA=-(h22*g1 - h21 * g2) / det;
		dB=-(-h21*g1+ h11 * g2) / det;
		gd=g1*dA+g2*dB;


		stepsize = 1; 		// Line Search
		while (stepsize >= min_step)
		{
			newA = A + stepsize * dA;
			newB = B + stepsize * dB;

			// New function value
			newf = 0.0;
			for (i=0;i<l;i++)
			{
				fApB = dec_values[i]*newA+newB;
				if (fApB >= 0)
					newf += t[i]*fApB + log(1+exp(-fApB));
				else
					newf += (t[i] - 1)*fApB +log(1+exp(fApB));
			}
			// Check sufficient decrease
			if (newf<fval+0.0001*stepsize*gd)
			{
				A=newA;B=newB;fval=newf;
				break;
			}
			else
				stepsize = stepsize / 2.0;
		}

		if (stepsize < min_step)
		{
			info("Line search fails in two-class probability estimates\n");
			break;
		}
	}

	if (iter>=max_iter)
		info("Reaching maximal iterations in two-class probability estimates\n");
	free(t);
}

double SC_SVM::sigmoid_predict(double decision_value, double A, double B)
{
	double fApB = decision_value*A+B;
	if (fApB >= 0)
		return exp(-fApB)/(1.0+exp(-fApB));
	else
		return 1.0/(1+exp(fApB)) ;
}

// Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
void SC_SVM::multiclass_probability(int k, double **r, double *p)
{
	int t,j;
	int iter = 0, max_iter=sclib::max(100,k);
	double **Q=Malloc(double *,k);
	double *Qp=Malloc(double,k);
	double pQp, eps=0.005/k;
	
	for (t=0;t<k;t++)
	{
		p[t]=1.0/k;  // Valid if k = 1
		Q[t]=Malloc(double,k);
		Q[t][t]=0;
		for (j=0;j<t;j++)
		{
			Q[t][t]+=r[j][t]*r[j][t];
			Q[t][j]=Q[j][t];
		}
		for (j=t+1;j<k;j++)
		{
			Q[t][t]+=r[j][t]*r[j][t];
			Q[t][j]=-r[j][t]*r[t][j];
		}
	}
	for (iter=0;iter<max_iter;iter++)
	{
		// stopping condition, recalculate QP,pQP for numerical accuracy
		pQp=0;
		for (t=0;t<k;t++)
		{
			Qp[t]=0;
			for (j=0;j<k;j++)
				Qp[t]+=Q[t][j]*p[j];
			pQp+=p[t]*Qp[t];
		}
		double max_error=0;
		for (t=0;t<k;t++)
		{
			double error=fabs(Qp[t]-pQp);
			if (error>max_error)
				max_error=error;
		}
		if (max_error<eps) break;
		
		for (t=0;t<k;t++)
		{
			double diff=(-Qp[t]+pQp)/Q[t][t];
			p[t]+=diff;
			pQp=(pQp+diff*(diff*Q[t][t]+2*Qp[t]))/(1+diff)/(1+diff);
			for (j=0;j<k;j++)
			{
				Qp[j]=(Qp[j]+diff*Q[t][j])/(1+diff);
				p[j]/=(1+diff);
			}
		}
	}
	if (iter>=max_iter)
		info("Exceeds max_iter in multiclass_prob\n");
	for(t=0;t<k;t++) free(Q[t]);
	free(Q);
	free(Qp);
}

// Cross-validation decision values for probability estimates
void SC_SVM::svm_binary_svc_probability(const SC_SVMproblem *prob, const SC_SVMparameter *param, double Cp, double Cn, double& probA, double& probB)
{
	int i;
	int nr_fold = 3; //by thilo: changed from 5 to 3 to save time...
	int *perm = Malloc(int,prob->l);
	double *dec_values = Malloc(double,prob->l);

	// random shuffle
	for(i=0;i<prob->l;i++) perm[i]=i;
	for(i=0;i<prob->l;i++)
	{
		int j = i+rand()%(prob->l-i);
		swap(perm[i],perm[j]);
	}
	for(i=0;i<nr_fold;i++)
	{
		int begin = i*prob->l/nr_fold;
		int end = (i+1)*prob->l/nr_fold;
		int j,k;
		SC_SVMproblem subprob;

		subprob.l = prob->l-(end-begin);
		subprob.x = Malloc(SC_SVMnode*,subprob.l);
		subprob.y = Malloc(double,subprob.l);

		if (prob->W) subprob.W = Malloc(double, subprob.l); //by thilo: for Weighted Version (treatment of W here was somehow missing)
		else subprob.W = NULL; //by thilo: for Weighted Version (treatment of W here was somehow missing)

		k=0;
		for(j=0;j<begin;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			if (subprob.W) subprob.W[k] = prob->W[perm[j]]; //by thilo: for Weighted Version (treatment of W here was somehow missing)
			++k;
		}
		for(j=end;j<prob->l;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			if (subprob.W) subprob.W[k] = prob->W[perm[j]]; //by thilo: for Weighted Version (treatment of W here was somehow missing)
			++k;
		}
		int p_count=0,n_count=0;
		for(j=0;j<k;j++)
			if(subprob.y[j]>0)
				p_count++;
			else
				n_count++;

		if(p_count==0 && n_count==0)
			for(j=begin;j<end;j++)
				dec_values[perm[j]] = 0;
		else if(p_count > 0 && n_count == 0)
			for(j=begin;j<end;j++)
				dec_values[perm[j]] = 1;
		else if(p_count == 0 && n_count > 0)
			for(j=begin;j<end;j++)
				dec_values[perm[j]] = -1;
		else
		{
			SC_SVMparameter subparam = *param;
			subparam.probability=0;
			subparam.C=1.0;
			subparam.nr_weight=2;
			subparam.weight_label = Malloc(int,2);
			subparam.weight = Malloc(double,2);
			subparam.weight_label[0]=+1;
			subparam.weight_label[1]=-1;
			subparam.weight[0]=Cp;
			subparam.weight[1]=Cn;
			SC_SVMmodel *submodel = svm_train(&subprob,&subparam);
			for(j=begin;j<end;j++)
			{
				svm_predict_values(submodel,prob->x[perm[j]],&(dec_values[perm[j]])); 
				// ensure +1 -1 order; reason not using CV subroutine
				dec_values[perm[j]] *= submodel->label[0];
			}		
			svm_destroy_model(submodel);
			svm_destroy_param(&subparam);
			free(subprob.x);
			free(subprob.y);
			if (subprob.W) free(subprob.W); //by thilo: for Weighted Version (treatment of W here was somehow missing)
		}
	}		
	sigmoid_train(prob->l,dec_values,prob->y,probA,probB);
	free(dec_values);
	free(perm);
}

// Return parameter of a Laplace distribution 
double SC_SVM::svm_svr_probability(const SC_SVMproblem *prob, const SC_SVMparameter *param)
{
	int i;
	int nr_fold = 5;
	double *ymv = Malloc(double,prob->l);
	double mae = 0;

	SC_SVMparameter newparam = *param;
	newparam.probability = 0;
	svm_cross_validation(prob,&newparam,nr_fold,ymv);
	for(i=0;i<prob->l;i++)
	{
		ymv[i]=prob->y[i]-ymv[i];
		mae += fabs(ymv[i]);
	}		
	mae /= prob->l;
	double std=sqrt(2*mae*mae);
	int count=0;
	mae=0;
	for(i=0;i<prob->l;i++)
	        if (fabs(ymv[i]) > 5*std) 
                        count=count+1;
		else 
		        mae+=fabs(ymv[i]);
	mae /= (prob->l-count);
	info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma= %g\n",mae);
	free(ymv);
	return mae;
}


// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
void SC_SVM::svm_group_classes(const SC_SVMproblem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm)
{
	int l = prob->l;
	int max_nr_class = 16;
	int nr_class = 0;
	int *label = Malloc(int,max_nr_class);
	int *count = Malloc(int,max_nr_class);
	int *data_label = Malloc(int,l);	
	int i;

	for(i=0;i<l;i++)
	{
		int this_label = (int)prob->y[i];
		int j;
		for(j=0;j<nr_class;j++)
		{
			if(this_label == label[j])
			{
				++count[j];
				break;
			}
		}
		data_label[i] = j;
		if(j == nr_class)
		{
			if(nr_class == max_nr_class)
			{
				max_nr_class *= 2;
				label = (int *)realloc(label,max_nr_class*sizeof(int));
				count = (int *)realloc(count,max_nr_class*sizeof(int));
			}
			label[nr_class] = this_label;
			count[nr_class] = 1;
			++nr_class;
		}
	}

	int *start = Malloc(int,nr_class);
	start[0] = 0;
	for(i=1;i<nr_class;i++)
		start[i] = start[i-1]+count[i-1];
	for(i=0;i<l;i++)
	{
		perm[start[data_label[i]]] = i;
		++start[data_label[i]];
	}
	start[0] = 0;
	for(i=1;i<nr_class;i++)
		start[i] = start[i-1]+count[i-1];

	*nr_class_ret = nr_class;
	*label_ret = label;
	*start_ret = start;
	*count_ret = count;
	free(data_label);
}
#endif

SC_SVM::SC_SVM(bool verbose) {
  this->verbose = verbose;

#ifdef SC_USE_LIBSVM
  //for compatibility-issues with saved models from libsvm, retain the old values here
  MArray_1D(this->svm_type_table, 6, char*, "SC_SVM: svm_type_table");
  MArray_1D(this->svm_type_table[0], 6, char, "SC_SVM: svm_type_table[0]");
  MArray_1D(this->svm_type_table[1], 7, char, "SC_SVM: svm_type_table[1]");
  MArray_1D(this->svm_type_table[2], 10, char, "SC_SVM: svm_type_table[2]");
  MArray_1D(this->svm_type_table[3], 12, char, "SC_SVM: svm_type_table[3]");
  MArray_1D(this->svm_type_table[4], 7, char, "SC_SVM: svm_type_table[4]");
  sprintf(this->svm_type_table[0], "%s\0", "c_svc");
  sprintf(this->svm_type_table[1], "%s\0", "nu_svc");
  sprintf(this->svm_type_table[2], "%s\0", "one_class");
  sprintf(this->svm_type_table[3], "%s\0", "epsilon_svr");
  sprintf(this->svm_type_table[4], "%s\0", "nu_svr");
  this->svm_type_table[5] = NULL;

  MArray_1D(this->kernel_type_table, 5, char*, "SC_SVM: kernel_type_table");
  MArray_1D(this->kernel_type_table[0], 7, char, "SC_SVM: kernel_type_table[0]");
  MArray_1D(this->kernel_type_table[1], 11, char, "SC_SVM: kernel_type_table[1]");
  MArray_1D(this->kernel_type_table[2], 4, char, "SC_SVM: kernel_type_table[2]");
  MArray_1D(this->kernel_type_table[3], 8, char, "SC_SVM: kernel_type_table[3]");
  sprintf(this->kernel_type_table[0], "%s\0", "linear");
  sprintf(this->kernel_type_table[1], "%s\0", "polynomial");
  sprintf(this->kernel_type_table[2], "%s\0", "rbf");
  sprintf(this->kernel_type_table[3], "%s\0", "sigmoid");
  this->kernel_type_table[4] = NULL;
#else
	REPORT_ERROR(SVLIB_Fail, "libsvm code not available"
#endif
}

SC_SVM::~SC_SVM() {
#ifdef SC_USE_LIBSVM
	for (int i = 0; i < 5; i++) {
    if (i < 4) {
      MFree_1D(this->kernel_type_table[i]);
    }
    MFree_1D(this->svm_type_table[i]);
  }
  MFree_1D(this->kernel_type_table);
  MFree_1D(this->svm_type_table);
#endif
}

//
// Interface functions
//
SC_SVMmodel* SC_SVM::svm_train(const SC_SVMproblem *prob, const SC_SVMparameter *param)
{
#ifdef SC_USE_LIBSVM
	SC_SVMmodel *model = Malloc(SC_SVMmodel,1);
	model->param = *param;
	model->free_sv = 0;	// XXX

	if(param->svm_type == SCLIB_SVM_TYPE_ONECLASS ||
	   param->svm_type == SCLIB_SVM_TYPE_EPSILONSVR ||
	   param->svm_type == SCLIB_SVM_TYPE_NUSVR)
	{
		// regression or one-class-svm
		model->nr_class = 2;
		model->label = NULL;
		model->nSV = NULL;
		model->probA = NULL; model->probB = NULL;
		model->sv_coef = Malloc(double *,1);

		if(param->probability && 
		   (param->svm_type == SCLIB_SVM_TYPE_EPSILONSVR ||
		    param->svm_type == SCLIB_SVM_TYPE_NUSVR))
		{
			model->probA = Malloc(double,1);
			model->probA[0] = svm_svr_probability(prob,param);
		}

    SC_SVM::SC_SVMdecisionFunction f = svm_train_one(prob,param,0,0);
		model->rho = Malloc(double,1);
		model->rho[0] = f.rho;

		int nSV = 0;
		int i;
		for(i=0;i<prob->l;i++)
			if(fabs(f.alpha[i]) > 0) ++nSV;
		model->l = nSV;
		model->SV = Malloc(SC_SVMnode *,nSV);
		model->sv_coef[0] = Malloc(double,nSV);
		int j = 0;
		for(i=0;i<prob->l;i++)
			if(fabs(f.alpha[i]) > 0)
			{
				model->SV[j] = prob->x[i];
				model->sv_coef[0][j] = f.alpha[i];
				++j;
			}		

		free(f.alpha);
	}
	else
	{
		// classification
		int l = prob->l;
		int nr_class;
		int *label = NULL;
		int *start = NULL;
		int *count = NULL;
		int *perm = Malloc(int,l);

		// group training data of the same class
		svm_group_classes(prob,&nr_class,&label,&start,&count,perm);		
		SC_SVMnode **x = Malloc(SC_SVMnode *,l);
		double *W = NULL; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
		if (prob->W) W = Malloc(double, l); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)

		int i;
		for(i=0;i<l;i++)
		{ //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			x[i] = prob->x[perm[i]];
			if (prob->W) W[i] = prob->W[perm[i]]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
		} //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)

		// calculate weighted C

		double *weighted_C = Malloc(double, nr_class);
		for(i=0;i<nr_class;i++)
			weighted_C[i] = param->C;
		for(i=0;i<param->nr_weight;i++)
		{	
			int j;
			for(j=0;j<nr_class;j++)
				if(param->weight_label[i] == label[j])
					break;
      if(j == nr_class) {
        char *errStr = new char[sclib::bufferSize];
        sprintf(errStr, "warning: class label %d specified in weight is not found\n", param->weight_label[i]);
        REPORT_ERROR(SVLIB_BadArg, errStr);
        MFree_1D(errStr);
        //fprintf(stderr,"warning: class label %d specified in weight is not found\n", param->weight_label[i]);
      } else {
				weighted_C[j] *= param->weight[i];
      }
		}

		// train k*(k-1)/2 models
		
		bool *nonzero = Malloc(bool,l);
		for(i=0;i<l;i++)
			nonzero[i] = false;
		SC_SVMdecisionFunction *f = Malloc(SC_SVMdecisionFunction,nr_class*(nr_class-1)/2);

		double *probA=NULL,*probB=NULL;
		if (param->probability)
		{
			probA=Malloc(double,nr_class*(nr_class-1)/2);
			probB=Malloc(double,nr_class*(nr_class-1)/2);
		}

		int p = 0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				SC_SVMproblem sub_prob;
				int si = start[i], sj = start[j];
				int ci = count[i], cj = count[j];
				sub_prob.l = ci+cj;
				sub_prob.x = Malloc(SC_SVMnode *,sub_prob.l);
				sub_prob.y = Malloc(double,sub_prob.l);
				if (W) sub_prob.W = Malloc(double, sub_prob.l); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
				else sub_prob.W = NULL; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
				int k;
				for(k=0;k<ci;k++)
				{
					sub_prob.x[k] = x[si+k];
					sub_prob.y[k] = +1;
					if (sub_prob.W) sub_prob.W[k] = W[si+k]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
				}
				for(k=0;k<cj;k++)
				{
					sub_prob.x[ci+k] = x[sj+k];
					sub_prob.y[ci+k] = -1;
					if (sub_prob.W) sub_prob.W[ci+k] = W[sj+k]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
				}

				if(param->probability)
					svm_binary_svc_probability(&sub_prob,param,weighted_C[i],weighted_C[j],probA[p],probB[p]);

				f[p] = svm_train_one(&sub_prob,param,weighted_C[i],weighted_C[j]);
				for(k=0;k<ci;k++)
					if(!nonzero[si+k] && fabs(f[p].alpha[k]) > 0)
						nonzero[si+k] = true;
				for(k=0;k<cj;k++)
					if(!nonzero[sj+k] && fabs(f[p].alpha[ci+k]) > 0)
						nonzero[sj+k] = true;
				free(sub_prob.x);
				free(sub_prob.y);
				if (sub_prob.W) free(sub_prob.W); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
				++p;
			}

		// build output

		model->nr_class = nr_class;
		
		model->label = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
			model->label[i] = label[i];
		
		model->rho = Malloc(double,nr_class*(nr_class-1)/2);
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			model->rho[i] = f[i].rho;

		if(param->probability)
		{
			model->probA = Malloc(double,nr_class*(nr_class-1)/2);
			model->probB = Malloc(double,nr_class*(nr_class-1)/2);
			for(i=0;i<nr_class*(nr_class-1)/2;i++)
			{
				model->probA[i] = probA[i];
				model->probB[i] = probB[i];
			}
		}
		else
		{
			model->probA=NULL;
			model->probB=NULL;
		}

		int total_sv = 0;
		int *nz_count = Malloc(int,nr_class);
		model->nSV = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
		{
			int nSV = 0;
			for(int j=0;j<count[i];j++)
				if(nonzero[start[i]+j])
				{	
					++nSV;
					++total_sv;
				}
			model->nSV[i] = nSV;
			nz_count[i] = nSV;
		}
		
		info("Total nSV = %d\n",total_sv);

		model->l = total_sv;
		model->SV = Malloc(SC_SVMnode *,total_sv);
		p = 0;
		for(i=0;i<l;i++)
			if(nonzero[i]) model->SV[p++] = x[i];

		int *nz_start = Malloc(int,nr_class);
		nz_start[0] = 0;
		for(i=1;i<nr_class;i++)
			nz_start[i] = nz_start[i-1]+nz_count[i-1];

		model->sv_coef = Malloc(double *,nr_class-1);
		for(i=0;i<nr_class-1;i++)
			model->sv_coef[i] = Malloc(double,total_sv);

		p = 0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				// classifier (i,j): coefficients with
				// i are in sv_coef[j-1][nz_start[i]...],
				// j are in sv_coef[i][nz_start[j]...]

				int si = start[i];
				int sj = start[j];
				int ci = count[i];
				int cj = count[j];
				
				int q = nz_start[i];
				int k;
				for(k=0;k<ci;k++)
					if(nonzero[si+k])
						model->sv_coef[j-1][q++] = f[p].alpha[k];
				q = nz_start[j];
				for(k=0;k<cj;k++)
					if(nonzero[sj+k])
						model->sv_coef[i][q++] = f[p].alpha[ci+k];
				++p;
			}
		
		free(label);
		free(probA);
		free(probB);
		free(count);
		free(perm);
		free(start);
		if (W) free(W); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
		free(x);
		free(weighted_C);
		free(nonzero);
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			free(f[i].alpha);
		free(f);
		free(nz_count);
		free(nz_start);
	}
	return model;
#else
	return NULL;
#endif
}

// Stratified cross validation
void SC_SVM::svm_cross_validation(const SC_SVMproblem *prob, const SC_SVMparameter *param, int nr_fold, double *target, double &svFraction) //last parameter by thilo to return the average fraction of support vectors 
{
#ifdef SC_USE_LIBSVM
	int i;
	int *fold_start = Malloc(int,nr_fold+1);
	int l = prob->l;
	int *perm = Malloc(int,l);
	int nr_class;

	svFraction = 0.0; //by thilo

	// stratified cv may not give leave-one-out rate
	// Each class to l folds -> some folds may have zero elements
	if((param->svm_type == SCLIB_SVM_TYPE_CSVC ||
	    param->svm_type == SCLIB_SVM_TYPE_NUSVC) && nr_fold < l)
	{
		int *start = NULL;
		int *label = NULL;
		int *count = NULL;
		svm_group_classes(prob,&nr_class,&label,&start,&count,perm);

		// random shuffle and then data grouped by fold using the array perm
		int *fold_count = Malloc(int,nr_fold);
		int c;
		int *index = Malloc(int,l);
		for(i=0;i<l;i++)
			index[i]=perm[i];
		for (c=0; c<nr_class; c++) 
			for(i=0;i<count[c];i++)
			{
				int j = i+rand()%(count[c]-i);
				swap(index[start[c]+j],index[start[c]+i]);
			}
		for(i=0;i<nr_fold;i++)
		{
			fold_count[i] = 0;
			for (c=0; c<nr_class;c++)
				fold_count[i]+=(i+1)*count[c]/nr_fold-i*count[c]/nr_fold;
		}
		fold_start[0]=0;
		for (i=1;i<=nr_fold;i++)
			fold_start[i] = fold_start[i-1]+fold_count[i-1];
		for (c=0; c<nr_class;c++)
			for(i=0;i<nr_fold;i++)
			{
				int begin = start[c]+i*count[c]/nr_fold;
				int end = start[c]+(i+1)*count[c]/nr_fold;
				for(int j=begin;j<end;j++)
				{
					perm[fold_start[i]] = index[j];
					fold_start[i]++;
				}
			}
		fold_start[0]=0;
		for (i=1;i<=nr_fold;i++)
			fold_start[i] = fold_start[i-1]+fold_count[i-1];
		free(start);	
		free(label);
		free(count);	
		free(index);
		free(fold_count);
	}
	else
	{
		for(i=0;i<l;i++) perm[i]=i;
		for(i=0;i<l;i++)
		{
			int j = i+rand()%(l-i);
			swap(perm[i],perm[j]);
		}
		for(i=0;i<=nr_fold;i++)
			fold_start[i]=i*l/nr_fold;
	}

	for(i=0;i<nr_fold;i++)
	{
		int begin = fold_start[i];
		int end = fold_start[i+1];
		int j,k;
		SC_SVMproblem subprob;

		subprob.l = l-(end-begin);
		subprob.x = Malloc(SC_SVMnode*,subprob.l);
		subprob.y = Malloc(double,subprob.l);
		if (prob->W) subprob.W = Malloc(double,subprob.l); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
		else subprob.W = NULL; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			
		k=0;
		for(j=0;j<begin;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			if (subprob.W) subprob.W[k] = prob->W[perm[j]]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			++k;
		}
		for(j=end;j<l;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			if (subprob.W) subprob.W[k] = prob->W[perm[j]]; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
			++k;
		}
		SC_SVMmodel *submodel = svm_train(&subprob,param);
		if(param->probability && 
		   (param->svm_type == SCLIB_SVM_TYPE_CSVC || param->svm_type == SCLIB_SVM_TYPE_NUSVC))
		{
			double *prob_estimates=Malloc(double,svm_get_nr_class(submodel));
			for(j=begin;j<end;j++)
				target[perm[j]] = svm_predict_probability(submodel,prob->x[perm[j]],prob_estimates);
			free(prob_estimates);			
		}
		else
			for(j=begin;j<end;j++)
				target[perm[j]] = svm_predict(submodel,prob->x[perm[j]]);
		
		//block by thilo
		int nSV = 0;
		if (submodel->nSV == NULL) {
			nSV += submodel->l;
		} else {
			for (j = 0; j < nr_class; j++) {
				nSV += submodel->nSV[j];
			}
		}
		svFraction += (double)(nSV) / (double)(subprob.l);

		svm_destroy_model(submodel);
		free(subprob.x);
		free(subprob.y);
		if (subprob.W) free(subprob.W); //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
	}	
	svFraction /= (double)(nr_fold); //by thilo
	free(fold_start);
	free(perm);	
#endif
}


int SC_SVM::svm_get_svm_type(const SC_SVMmodel *model)
{
#ifdef SC_USE_LIBSVM
	return model->param.svm_type;
#else
	return SVLIB_Fail;
#endif
}

int SC_SVM::svm_get_nr_class(const SC_SVMmodel *model)
{
#ifdef SC_USE_LIBSVM
	return model->nr_class;
#else
	return SVLIB_Fail;
#endif
}

void SC_SVM::svm_get_labels(const SC_SVMmodel *model, int* label)
{
#ifdef SC_USE_LIBSVM
	if (model->label != NULL)
		for(int i=0;i<model->nr_class;i++)
			label[i] = model->label[i];
#endif
}

double SC_SVM::svm_get_svr_probability(const SC_SVMmodel *model)
{
#ifdef SC_USE_LIBSVM
	if ((model->param.svm_type == SCLIB_SVM_TYPE_EPSILONSVR || model->param.svm_type == SCLIB_SVM_TYPE_NUSVR) &&
	    model->probA!=NULL)
		return model->probA[0];
	else
	{
		info("Model doesn't contain information for SVR probability inference\n");
		return 0;
	}
#else
	return (double)(SVLIB_Fail);
#endif
}

void SC_SVM::svm_predict_values(const SC_SVMmodel *model, const SC_SVMnode *x, double* dec_values)
{
#ifdef SC_USE_LIBSVM
	if(model->param.svm_type == SCLIB_SVM_TYPE_ONECLASS ||
	   model->param.svm_type == SCLIB_SVM_TYPE_EPSILONSVR ||
	   model->param.svm_type == SCLIB_SVM_TYPE_NUSVR)
	{
		double *sv_coef = model->sv_coef[0];
		double sum = 0;
		for(int i=0;i<model->l;i++)
			sum += sv_coef[i] * Kernel::k_function(x,model->SV[i],model->param);
		sum -= model->rho[0];
		*dec_values = sum;
	}
	else
	{
		int i;
		int nr_class = model->nr_class;
		int l = model->l;
		
		double *kvalue = Malloc(double,l);
		for(i=0;i<l;i++)
			kvalue[i] = Kernel::k_function(x,model->SV[i],model->param);

		int *start = Malloc(int,nr_class);
		start[0] = 0;
		for(i=1;i<nr_class;i++)
			start[i] = start[i-1]+model->nSV[i-1];

		int p=0;
		int pos=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				double sum = 0;
				int si = start[i];
				int sj = start[j];
				int ci = model->nSV[i];
				int cj = model->nSV[j];
				
				int k;
				double *coef1 = model->sv_coef[j-1];
				double *coef2 = model->sv_coef[i];
				for(k=0;k<ci;k++)
					sum += coef1[si+k] * kvalue[si+k];
				for(k=0;k<cj;k++)
					sum += coef2[sj+k] * kvalue[sj+k];
				sum -= model->rho[p++];
				dec_values[pos++] = sum;
			}

		free(kvalue);
		free(start);
	}
#endif
}

double SC_SVM::svm_predict(const SC_SVMmodel *model, const SC_SVMnode *x)
{
#ifdef SC_USE_LIBSVM
	if(model->param.svm_type == SCLIB_SVM_TYPE_ONECLASS ||
	   model->param.svm_type == SCLIB_SVM_TYPE_EPSILONSVR ||
	   model->param.svm_type == SCLIB_SVM_TYPE_NUSVR)
	{
		double res;
		svm_predict_values(model, x, &res);
		
		if(model->param.svm_type == SCLIB_SVM_TYPE_ONECLASS)
			return (res>0)?1:-1;
		else
			return res;
	}
	else
	{
		int i;
		int nr_class = model->nr_class;
		double *dec_values = Malloc(double, nr_class*(nr_class-1)/2);
		svm_predict_values(model, x, dec_values);

		int *vote = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
			vote[i] = 0;
		int pos=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				if(dec_values[pos++] > 0)
					++vote[i];
				else
					++vote[j];
			}

		int vote_max_idx = 0;
		for(i=1;i<nr_class;i++)
			if(vote[i] > vote[vote_max_idx])
				vote_max_idx = i;
		free(vote);
		free(dec_values);
		return model->label[vote_max_idx];
	}
#else
	return (double)(SVLIB_Fail);
#endif
}

double SC_SVM::svm_predict_probability(const SC_SVMmodel *model, const SC_SVMnode *x, double *prob_estimates)
{
#ifdef SC_USE_LIBSVM
	if ((model->param.svm_type == SCLIB_SVM_TYPE_CSVC || model->param.svm_type == SCLIB_SVM_TYPE_NUSVC) &&
	    model->probA!=NULL && model->probB!=NULL)
	{
		int i;
		int nr_class = model->nr_class;
		double *dec_values = Malloc(double, nr_class*(nr_class-1)/2);
		svm_predict_values(model, x, dec_values);

		double min_prob=1e-7;
		double **pairwise_prob=Malloc(double *,nr_class);
		for(i=0;i<nr_class;i++)
			pairwise_prob[i]=Malloc(double,nr_class);
		int k=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				pairwise_prob[i][j]=sclib::min(sclib::max(sigmoid_predict(dec_values[k],model->probA[k],model->probB[k]),min_prob),1-min_prob);
				pairwise_prob[j][i]=1-pairwise_prob[i][j];
				k++;
			}
		multiclass_probability(nr_class,pairwise_prob,prob_estimates);

		int prob_max_idx = 0;
		for(i=1;i<nr_class;i++)
			if(prob_estimates[i] > prob_estimates[prob_max_idx])
				prob_max_idx = i;
		for(i=0;i<nr_class;i++)
			free(pairwise_prob[i]);
		free(dec_values);
                free(pairwise_prob);	     
		return model->label[prob_max_idx];
	}
	else 
		return svm_predict(model, x);
#else
	return (double)(SVLIB_Fail);
#endif
}

int SC_SVM::svm_save_model(const char *model_file_name, const SC_SVMmodel *model)
{
#ifdef SC_USE_LIBSVM
	FILE *fp = fopen(model_file_name,"w");
	if(fp==NULL) return SVLIB_Fail;

	const SC_SVMparameter& param = model->param;

  fprintf(fp,"svm_type %s\n", this->svm_type_table[param.svm_type]);
	fprintf(fp,"kernel_type %s\n", this->kernel_type_table[param.kernel_type]);

	if(param.kernel_type == SCLIB_SVM_KERNEL_POLY)
		fprintf(fp,"degree %g\n", param.degree);

	if(param.kernel_type == SCLIB_SVM_KERNEL_POLY || param.kernel_type == SCLIB_SVM_KERNEL_RBF || param.kernel_type == SCLIB_SVM_KERNEL_SIGMOID)
		fprintf(fp,"gamma %g\n", param.gamma);

	if(param.kernel_type == SCLIB_SVM_KERNEL_POLY || param.kernel_type == SCLIB_SVM_KERNEL_SIGMOID)
		fprintf(fp,"coef0 %g\n", param.coef0);

	int nr_class = model->nr_class;
	int l = model->l;
	fprintf(fp, "nr_class %d\n", nr_class);
	fprintf(fp, "total_sv %d\n",l);
	
	{
		fprintf(fp, "rho");
		for(int i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->rho[i]);
		fprintf(fp, "\n");
	}
	
	if(model->label)
	{
		fprintf(fp, "label");
		for(int i=0;i<nr_class;i++)
			fprintf(fp," %d",model->label[i]);
		fprintf(fp, "\n");
	}

	if(model->probA) // regression has probA only
	{
		fprintf(fp, "probA");
		for(int i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->probA[i]);
		fprintf(fp, "\n");
	}
	if(model->probB)
	{
		fprintf(fp, "probB");
		for(int i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->probB[i]);
		fprintf(fp, "\n");
	}

	if(model->nSV)
	{
		fprintf(fp, "nr_sv");
		for(int i=0;i<nr_class;i++)
			fprintf(fp," %d",model->nSV[i]);
		fprintf(fp, "\n");
	}

	fprintf(fp, "SV\n");
	const double * const *sv_coef = model->sv_coef;
	const SC_SVMnode * const *SV = model->SV;

	for(int i=0;i<l;i++)
	{
		for(int j=0;j<nr_class-1;j++)
			fprintf(fp, "%.16g ",sv_coef[j][i]);

		const SC_SVMnode *p = SV[i];
		while(p->index != -1)
		{
			fprintf(fp,"%d:%.8g ",p->index,p->value);
			p++;
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	return SVLIB_Ok;
#else
	return SVLIB_Fail;
#endif
}

SC_SVMmodel* SC_SVM::svm_load_model(const char *model_file_name)
{
#ifdef SC_USE_LIBSVM
	FILE *fp = fopen(model_file_name,"rb");
	if(fp==NULL) return NULL;
	
	// read parameters

	SC_SVMmodel *model = Malloc(SC_SVMmodel,1);
	SC_SVMparameter& param = model->param;
	model->rho = NULL;
	model->probA = NULL;
	model->probB = NULL;
	model->label = NULL;
	model->nSV = NULL;

	//by thilo: to allow freeing
	param.weight = NULL;
	param.weight_label = NULL;

	char cmd[81];
	while(1)
	{
		fscanf(fp,"%80s",cmd);

		if(strcmp(cmd,"svm_type")==0)
		{
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;svm_type_table[i];i++)
			{
				if(strcmp(svm_type_table[i],cmd)==0)
				{
					param.svm_type=i;
					break;
				}
			}
			if(svm_type_table[i] == NULL)
			{
        REPORT_ERROR(SVLIB_BadData, "unknown svm type.");
				free(model->rho);
				free(model->label);
				free(model->nSV);
				free(model);
				return NULL;
			}
		}
		else if(strcmp(cmd,"kernel_type")==0)
		{		
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;kernel_type_table[i];i++)
			{
				if(strcmp(kernel_type_table[i],cmd)==0)
				{
					param.kernel_type=i;
					break;
				}
			}
			if(kernel_type_table[i] == NULL)
			{
        REPORT_ERROR(SVLIB_BadData, "unknown kernel function.");
				free(model->rho);
				free(model->label);
				free(model->nSV);
				free(model);
				return NULL;
			}
		}
		else if(strcmp(cmd,"degree")==0)
			fscanf(fp,"%lf",&param.degree);
		else if(strcmp(cmd,"gamma")==0)
			fscanf(fp,"%lf",&param.gamma);
		else if(strcmp(cmd,"coef0")==0)
			fscanf(fp,"%lf",&param.coef0);
		else if(strcmp(cmd,"nr_class")==0)
			fscanf(fp,"%d",&model->nr_class);
		else if(strcmp(cmd,"total_sv")==0)
			fscanf(fp,"%d",&model->l);
		else if(strcmp(cmd,"rho")==0)
		{
			int n = model->nr_class * (model->nr_class-1)/2;
			model->rho = Malloc(double,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%lf",&model->rho[i]);
		}
		else if(strcmp(cmd,"label")==0)
		{
			int n = model->nr_class;
			model->label = Malloc(int,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%d",&model->label[i]);
		}
		else if(strcmp(cmd,"probA")==0)
		{
			int n = model->nr_class * (model->nr_class-1)/2;
			model->probA = Malloc(double,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%lf",&model->probA[i]);
		}
		else if(strcmp(cmd,"probB")==0)
		{
			int n = model->nr_class * (model->nr_class-1)/2;
			model->probB = Malloc(double,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%lf",&model->probB[i]);
		}
		else if(strcmp(cmd,"nr_sv")==0)
		{
			int n = model->nr_class;
			model->nSV = Malloc(int,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%d",&model->nSV[i]);
		}
		else if(strcmp(cmd,"SV")==0)
		{
			while(1)
			{
				int c = getc(fp);
				if(c==EOF || c=='\n') break;	
			}
			break;
		}
		else
		{
      REPORT_ERROR(SVLIB_BadData, "unknown text in model file.");
			free(model->rho);
			free(model->label);
			free(model->nSV);
			free(model);
			return NULL;
		}
	}

	// read sv_coef and SV

	int elements = 0;
	long pos = ftell(fp);

	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
			case '\n':
				// count the '-1' element
			case ':':
				++elements;
				break;
			case EOF:
				goto out;
			default:
				;
		}
	}
out:
	fseek(fp,pos,SEEK_SET);

	int m = model->nr_class - 1;
	int l = model->l;
	model->sv_coef = Malloc(double *,m);
	int i;
	for(i=0;i<m;i++)
		model->sv_coef[i] = Malloc(double,l);
	model->SV = Malloc(SC_SVMnode*,l);
	SC_SVMnode *x_space=NULL;
	if(l>0) x_space = Malloc(SC_SVMnode,elements);

	int j=0;
	for(i=0;i<l;i++)
	{
		model->SV[i] = &x_space[j];
		for(int k=0;k<m;k++)
			fscanf(fp,"%lf",&model->sv_coef[k][i]);
		while(1)
		{
			int c;
			do {
				c = getc(fp);
				if(c=='\n') goto out2;
			} while(isspace(c));
			ungetc(c,fp);
			fscanf(fp,"%d:%lf",&(x_space[j].index),&(x_space[j].value));
			++j;
		}	
out2:
		x_space[j++].index = -1;
	}

	fclose(fp);

	model->free_sv = 1;	// XXX
	return model;
#else
	return NULL;
#endif
}

void SC_SVM::svm_destroy_model(SC_SVMmodel* model)
{
#ifdef SC_USE_LIBSVM

	if(model->free_sv && model->l > 0)
		free((void *)(model->SV[0]));
	for(int i=0;i<model->nr_class-1;i++)
	 free(model->sv_coef[i]);
	free(model->SV);
	free(model->sv_coef);
	free(model->rho);
	free(model->label);
	free(model->probA);
	free(model->probB);
	free(model->nSV);
	free(model->param.weight); //by thilo
	free(model->param.weight_label); //by thilo
	free(model);
#endif
}

void SC_SVM::svm_destroy_param(SC_SVMparameter* param)
{
#ifdef SC_USE_LIBSVM
	if (param->weight_label != NULL) { //by thilo: handling of NULL
		free(param->weight_label);
		param->weight_label = NULL; //by thilo: handling of NULL
	} //by thilo: handling of NULL
	if (param->weight != NULL) { //by thilo: handling of NULL
		free(param->weight);
		param->weight = NULL; //by thilo: handling of NULL
	}
	param->nr_weight = 0; //by thilo
#endif
}

const char* SC_SVM::svm_check_parameter(const SC_SVMproblem *prob, const SC_SVMparameter *param)
{
#ifdef SC_USE_LIBSVM
	// svm_type

	int svm_type = param->svm_type;
	if(svm_type != SCLIB_SVM_TYPE_CSVC &&
	   svm_type != SCLIB_SVM_TYPE_NUSVC &&
	   svm_type != SCLIB_SVM_TYPE_ONECLASS &&
	   svm_type != SCLIB_SVM_TYPE_EPSILONSVR &&
	   svm_type != SCLIB_SVM_TYPE_NUSVR)
		return "unknown svm type";
	
	// kernel_type
	
	int kernel_type = param->kernel_type;
	if(kernel_type != SCLIB_SVM_KERNEL_LINEAR &&
	   kernel_type != SCLIB_SVM_KERNEL_POLY &&
	   kernel_type != SCLIB_SVM_KERNEL_RBF &&
	   kernel_type != SCLIB_SVM_KERNEL_SIGMOID)
		return "unknown kernel type";

	// cache_size,eps,C,nu,p,shrinking

	if(param->cache_size <= 0)
		return "cache_size <= 0";

	if(param->eps <= 0)
		return "eps <= 0";

	if(svm_type == SCLIB_SVM_TYPE_CSVC ||
	   svm_type == SCLIB_SVM_TYPE_EPSILONSVR ||
	   svm_type == SCLIB_SVM_TYPE_NUSVR)
		if(param->C <= 0)
			return "C <= 0";

	if(svm_type == SCLIB_SVM_TYPE_NUSVC ||
	   svm_type == SCLIB_SVM_TYPE_ONECLASS ||
	   svm_type == SCLIB_SVM_TYPE_NUSVR)
		if(param->nu <= 0 || param->nu > 1)
			return "nu <= 0 or nu > 1";

	if(svm_type == SCLIB_SVM_TYPE_EPSILONSVR)
		if(param->p < 0)
			return "p < 0";

	if(param->shrinking != 0 &&
	   param->shrinking != 1)
		return "shrinking != 0 and shrinking != 1";

	if(param->probability != 0 &&
	   param->probability != 1)
		return "probability != 0 and probability != 1";

	if(param->probability == 1 &&
	   svm_type == SCLIB_SVM_TYPE_ONECLASS)
		return "one-class SVM probability output not supported yet";


	// check whether nu-svc is feasible
	
	if(svm_type == SCLIB_SVM_TYPE_NUSVC)
	{
		int l = prob->l;
		int max_nr_class = 16;
		int nr_class = 0;
		int *label = Malloc(int,max_nr_class);
		int *count = Malloc(int,max_nr_class);

		int i;
		for(i=0;i<l;i++)
		{
			int this_label = (int)prob->y[i];
			int j;
			for(j=0;j<nr_class;j++)
				if(this_label == label[j])
				{
					++count[j];
					break;
				}
			if(j == nr_class)
			{
				if(nr_class == max_nr_class)
				{
					max_nr_class *= 2;
					label = (int *)realloc(label,max_nr_class*sizeof(int));
					count = (int *)realloc(count,max_nr_class*sizeof(int));
				}
				label[nr_class] = this_label;
				count[nr_class] = 1;
				++nr_class;
			}
		}
	
		for(i=0;i<nr_class;i++)
		{
			int n1 = count[i];
			for(int j=i+1;j<nr_class;j++)
			{
				int n2 = count[j];
				if(param->nu*(n1+n2)/2 > sclib::min(n1,n2))
				{
					free(label);
					free(count);
					return "specified nu is infeasible";
				}
			}
		}
		free(label);
		free(count);
	}
#endif

	return NULL;
}

int SC_SVM::svm_check_probability_model(const SC_SVMmodel *model)
{
#ifdef SC_USE_LIBSVM
	return ((model->param.svm_type == SCLIB_SVM_TYPE_CSVC || model->param.svm_type == SCLIB_SVM_TYPE_NUSVC) &&
		model->probA!=NULL && model->probB!=NULL) ||
		((model->param.svm_type == SCLIB_SVM_TYPE_EPSILONSVR || model->param.svm_type == SCLIB_SVM_TYPE_NUSVR) &&
		 model->probA!=NULL);
#else
	return SVLIB_Fail;
#endif
}
