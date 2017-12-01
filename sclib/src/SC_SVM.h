/**************************************************************************/
/*    This class is based on the libsvm-2.81-implementation of Support    */
/*    Vector Machines.                                                    */
/*    It is nothing more than an object-orientated wrapper around it.     */
/*                                                                        */
/*    Author  : Thilo Stadelmann																					*/
/*    Date    : 08.03.2006																								*/
/**************************************************************************/

#ifndef __SC_SVM_H__
#define __SC_SVM_H__

#include "SC_Api.h"
#include <SV_Error.h>
#include <SV_Data.h>

//used internally...
typedef float Qfloat;
typedef signed char schar;

//a SVM-node represents one entry in a vector
typedef struct {
	int index;
	double value;
} SC_SVMnode;

//a SVM-problem represents a feature-set including class labels
typedef struct {
	int l;
	double *y;
	SC_SVMnode **x;
	double *W; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
} SC_SVMproblem;

//a SVM-parameter represents all important prarameters to configure a SVM
struct SC_SVMparameter {
	int svm_type;
	int kernel_type;
	double degree;	    /* for SCLIB_SVM_KERNEL_POLY */
	double gamma;	      /* for SCLIB_SVM_KERNEL_POLY/SCLIB_SVM_KERNEL_RBF/SCLIB_SVM_KERNEL_SIGMOID */
	double coef0;	      /* for SCLIB_SVM_KERNEL_POLY/SCLIB_SVM_KERNEL_SIGMOID */

	/* these are for training only */
	double cache_size;  /* in MB */
	double eps;     	  /* stopping criteria */
	double C;	          /* for SCLIB_SVM_TYPE_CSVC, SCLIB_SVM_TYPE_EPSILONSVR and SCLIB_SVM_TYPE_NUSVR */
	int nr_weight;		  /* for SCLIB_SVM_TYPE_CSVC */
	int *weight_label;	/* for SCLIB_SVM_TYPE_CSVC */
	double* weight;		  /* for SCLIB_SVM_TYPE_CSVC */
	double nu;	        /* for SCLIB_SVM_TYPE_NUSVC, SCLIB_SVM_TYPE_ONECLASS, and SCLIB_SVM_TYPE_NUSVR */
	double p;	          /* for SCLIB_SVM_TYPE_EPSILONSVR */
	int shrinking;	    /* use the shrinking heuristics */
	int probability;    /* do probability estimates */

	SC_SVMparameter() { //by thilo: to support freeing
		this->weight = NULL;
		this->weight_label = NULL;
	};
};

//a SVM-model prepresents a trained SVM
struct SC_SVMmodel {
	SC_SVMparameter param;	// parameter
	int nr_class;       		// number of classes, = 2 in regression/one class svm
	int l;		            	// total #SV
	SC_SVMnode **SV;		    // SVs (SV[l])
	double **sv_coef;	      // coefficients for SVs in decision functions (sv_coef[n-1][l])
	double *rho;		        // constants in decision functions (rho[n*(n-1)/2])
	double *probA;          // pairwise probability information
	double *probB;
	
  // for classification only
	int *label;		          // label of each class (label[n])
	int *nSV;		            // number of SVs for each class (nSV[n])
  			                  // nSV[0] + nSV[1] + ... + nSV[n-1] = l
  // XXX
	int free_sv;		        // 1 if SC_SVMmodel is created by svm_load_model
				                  // 0 if SC_SVMmodel is created by svm_train
	SC_SVMmodel() { //by thilo: to support freeing
		this->param.weight = NULL;
		this->param.weight_label = NULL;
	}
};

//some types used inside
enum {SCLIB_SVM_TYPE_CSVC, SCLIB_SVM_TYPE_NUSVC, SCLIB_SVM_TYPE_ONECLASS, SCLIB_SVM_TYPE_EPSILONSVR, SCLIB_SVM_TYPE_NUSVR };	/* svm_type */
enum {SCLIB_SVM_KERNEL_LINEAR, SCLIB_SVM_KERNEL_POLY, SCLIB_SVM_KERNEL_RBF, SCLIB_SVM_KERNEL_SIGMOID };	/* kernel_type */

class SC_SVM {
  private:

  protected:

#ifdef SC_USE_LIBSVM
    class Cache {
    public:
	    Cache(int l,int size);
	    virtual ~Cache();
	    int get_data(const int index, Qfloat **data, int len); // request data [0,len), return some position p where [p,len) need to be filled (p >= len if nothing needs to be filled)
	    void swap_index(int i, int j);	// future_option
    private:
	    int l;
	    int size;
	    struct head_t {
		    head_t *prev, *next;	// a cicular list
		    Qfloat *data;
		    int len;		// data[0,len) is cached in this entry
	    };
	    head_t *head;
	    head_t lru_head;
	    void lru_delete(head_t *h);
	    void lru_insert(head_t *h);
    };

    class QMatrix {
      public:
	      virtual Qfloat *get_Q(int column, int len) const = 0;
	      virtual Qfloat *get_QD() const = 0;
	      virtual void swap_index(int i, int j) const = 0;
	      virtual ~QMatrix() {}
    };

    class Kernel: public SC_SVM::QMatrix {
      public:
	      Kernel(int l, SC_SVMnode * const * x, const SC_SVMparameter& param);
	      virtual ~Kernel();
	      static double k_function(const SC_SVMnode *x, const SC_SVMnode *y, const SC_SVMparameter& param);
	      virtual Qfloat *get_Q(int column, int len) const = 0;
	      virtual Qfloat *get_QD() const = 0;
	      virtual void swap_index(int i, int j) const;	// no so const...
      protected:
        double (SC_SVM::Kernel::*kernel_function)(int i, int j) const;
      private:
	      const SC_SVMnode **x;
	      double *x_square;
	      const int kernel_type; // SC_SVMparameter
	      const double degree;
	      const double gamma;
	      const double coef0;
	      static double dot(const SC_SVMnode *px, const SC_SVMnode *py);
	      double kernel_linear(int i, int j) const;
	      double kernel_poly(int i, int j) const;
	      double kernel_rbf(int i, int j) const;
	      double kernel_sigmoid(int i, int j) const;
    };

    class Solver {
      public:
	      Solver(bool verbose) {this->verbose = verbose;};
	      virtual ~Solver() {};
	      struct SolutionInfo {
		      double obj;
		      double rho;
		      double upper_bound_p;
		      double upper_bound_n;
		      double r;	// for Solver_NU
	      };
	      void Solve(int l, const QMatrix& Q, const double *b_, const schar *y_, double *alpha_, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking);
				//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu): the solver with different C_i
				void Solve(int l, const QMatrix& Q, const double *b_, const schar *y_, double *alpha_, const double *C_, double eps, SolutionInfo* si, int shrinking);
      protected:
	      int active_size;
	      schar *y;
	      double *G;		// gradient of objective function
				//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu): the enum is removed to allow examples to be both at the upper bound and at the lower bound (if C = 0) and replaced by static chars
				//enum { LOWER_BOUND, UPPER_BOUND, FREE };
				static const char FREE = 0;
				static const char LOWER_BOUND = 1; // bit 1
				static const char UPPER_BOUND = 2; // bit 2
				char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	      double *alpha;
	      const QMatrix *Q;
	      const Qfloat *QD;
	      double eps;
	      //double Cp,Cn;
				double *C; //Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu): the weighted C_i
	      double *b;
	      int *active_set;
	      double *G_bar;		// gradient, if we treat free variables as 0
	      int l;
	      bool unshrinked;	// XXX
	      double get_C(int i);
	      void update_alpha_status(int i);
				//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu): the following modification is for maintaining the alpha_status
				//bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
	      //bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
				bool is_upper_bound(int i) { return ((alpha_status[i] & UPPER_BOUND) != 0) ? true : false; }
				bool is_lower_bound(int i) { return ((alpha_status[i] & LOWER_BOUND) != 0) ? true : false; }
	      bool is_free(int i) { return alpha_status[i] == FREE; }
	      void swap_index(int i, int j);
	      void reconstruct_gradient();
	      virtual int select_working_set(int &i, int &j);
	      virtual int max_violating_pair(int &i, int &j);
	      virtual double calculate_rho();
	      virtual void do_shrinking();
        void info(char *fmt,...);
        void info_flush();
        bool verbose; //to decide whether to print debug informations...
    };

    class Solver_NU : public SC_SVM::Solver {
    public:
      Solver_NU(bool verbose) : Solver(verbose)  {}
	    void Solve(int l, const QMatrix& Q, const double *b, const schar *y, double *alpha, double Cp, double Cn, double eps, SolutionInfo* si, int shrinking);
    private:
	    SolutionInfo *si;
	    int select_working_set(int &i, int &j);
	    double calculate_rho();
	    void do_shrinking();
    };

    class SVC_Q: public SC_SVM::Kernel { 
    public:
	    SVC_Q(const SC_SVMproblem& prob, const SC_SVMparameter& param, const schar *y_);
	    Qfloat *get_Q(int i, int len) const;
	    Qfloat *get_QD() const;
      void swap_index(int i, int j) const;
	    ~SVC_Q();
    private:
	    schar *y;
	    Cache *cache;
	    Qfloat *QD;
    };

    class ONE_CLASS_Q: public SC_SVM::Kernel {
      public:
	      ONE_CLASS_Q(const SC_SVMproblem& prob, const SC_SVMparameter& param);
	      Qfloat *get_Q(int i, int len) const;
	      Qfloat *get_QD() const;
	      void swap_index(int i, int j) const;
	      ~ONE_CLASS_Q();
      private:
	      Cache *cache;
	      Qfloat *QD;
    };

    class SVR_Q: public SC_SVM::Kernel { 
      public:
	      SVR_Q(const SC_SVMproblem& prob, const SC_SVMparameter& param);
	      void swap_index(int i, int j) const;
	      Qfloat *get_Q(int i, int len) const;
	      Qfloat *get_QD() const;
	      ~SVR_Q();
      private:
	      int l;
	      Cache *cache;
	      schar *sign;
	      int *index;
	      mutable int next_buffer;
	      Qfloat *buffer[2];
	      Qfloat *QD;
    };

    // decision_function
    typedef struct {
	    double *alpha;
	    double rho;	
    } SC_SVMdecisionFunction;

    void svm_binary_svc_probability(const SC_SVMproblem *prob, const SC_SVMparameter *param, double Cp, double Cn, double& probA, double& probB);
    double svm_svr_probability(const SC_SVMproblem *prob, const SC_SVMparameter *param);
    SC_SVM::SC_SVMdecisionFunction svm_train_one(const SC_SVMproblem *prob, const SC_SVMparameter *param, double Cp, double Cn);
    void svm_group_classes(const SC_SVMproblem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm);
    void multiclass_probability(int k, double **r, double *p);
    void sigmoid_train(int l, const double *dec_values, const double *labels, double& A, double& B);
    double sigmoid_predict(double decision_value, double A, double B);
    
		//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
    //void solve_epsilon_svr(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, SC_SVM::Solver::SolutionInfo* si);
		void solve_epsilon_svr(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, Solver::SolutionInfo* si, const double *C_);
    void solve_one_class(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, SC_SVM::Solver::SolutionInfo* si);
    void solve_nu_svc(const SC_SVMproblem *prob, const SC_SVMparameter *param,	double *alpha, SC_SVM::Solver::SolutionInfo* si);
    void solve_nu_svr(const SC_SVMproblem *prob, const SC_SVMparameter *param, double *alpha, SC_SVM::Solver::SolutionInfo* si);
		//Weighted Version - Modified from LIBSVM 2.81 by Hsuan-Tien Lin (htlin at caltech.edu)
		//void solve_c_svc(const SC_SVMproblem *prob, const SC_SVMparameter* param, double *alpha, SC_SVM::Solver::SolutionInfo* si, double Cp, double Cn);
		void solve_c_svc(const SC_SVMproblem *prob, const SC_SVMparameter* param, double *alpha, Solver::SolutionInfo* si, const double *C_);

    void info(char *fmt,...);
    void info_flush();

    char **svm_type_table;
    char **kernel_type_table;
#endif

    bool verbose; //should debug-output be printed?

  public:

    SC_SVM(bool verbose = true);
    virtual ~SC_SVM();   

		void setVerbose(bool verbose) {this->verbose = verbose; return;} //by thilo

    SC_SVMmodel* svm_train(const SC_SVMproblem *prob, const SC_SVMparameter *param);
		void svm_cross_validation(const SC_SVMproblem *prob, const SC_SVMparameter *param, int nr_fold, double *target) { //by thilo: function to retain interface in the presence of newly introduced parameter below
			double tmp;
			return svm_cross_validation(prob, param ,nr_fold, target, tmp);
		}
		void svm_cross_validation(const SC_SVMproblem *prob, const SC_SVMparameter *param, int nr_fold, double *target, double &svFraction); //last parameter by thilo to return the average fraction of support vectors 

    int svm_save_model(const char *model_file_name, const SC_SVMmodel *model);
    SC_SVMmodel* svm_load_model(const char *model_file_name);

    int svm_get_svm_type(const SC_SVMmodel *model);
    int svm_get_nr_class(const SC_SVMmodel *model);
    void svm_get_labels(const SC_SVMmodel *model, int *label);
    double svm_get_svr_probability(const SC_SVMmodel *model);

    void svm_predict_values(const SC_SVMmodel *model, const SC_SVMnode *x, double* dec_values);
    double svm_predict(const SC_SVMmodel *model, const SC_SVMnode *x);
    double svm_predict_probability(const SC_SVMmodel *model, const SC_SVMnode *x, double* prob_estimates);

    void svm_destroy_model(SC_SVMmodel *model);
    void svm_destroy_param(SC_SVMparameter *param);

    const char *svm_check_parameter(const SC_SVMproblem *prob, const SC_SVMparameter *param);
    int svm_check_probability_model(const SC_SVMmodel *model);

		static double evaluateKernel(const SC_SVMnode *x, const SC_SVMnode *y, const SC_SVMparameter& param) {
#ifdef SC_USE_LIBSVM
			return Kernel::k_function(x, y, param);
#else
			return (double)(SVLIB_Fail);
#endif
		}
};

#endif
