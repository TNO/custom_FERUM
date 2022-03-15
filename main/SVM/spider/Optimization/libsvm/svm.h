#ifndef _LIBSVM_H
#define _LIBSVM_H


struct svm_node
{
	int index;
	double value;
};

struct svm_problem
{
	int l;
	double *y;
	svm_node **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID };	/* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	double degree;	// for poly
	double gamma;	// for poly/rbf/sigmoid
	double coef0;	// for poly/sigmoid

	// these are for training only
	double cache_size; // in MB
	double eps;	// stopping criteria
	double C;	// for C_SVC, EPSILON_SVR and NU_SVR
	int nr_weight;		// for C_SVC
	int *weight_label;	// for C_SVC
	double* weight;		// for C_SVC
	double nu;	// for NU_SVC, ONE_CLASS, and NU_SVR
	double p;	// for EPSILON_SVR
	int shrinking;	// use the shrinking heuristics
};

//
// svm_model
//
struct svm_model
{
	svm_parameter param;	// parameter
	int nr_class;		// number of classes, = 2 in regression/one class svm
	int l;			// total #SV
	svm_node **SV;		// SVs (SV[l])
	double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[n-1][l])
	double *rho;		// constants in decision functions (rho[n*(n-1)/2])

	// for classification only

	int *label;		// label of each class (label[n])
	int *nSV;		// number of SVs for each class (nSV[n])
				// nSV[0] + nSV[1] + ... + nSV[n-1] = l
	// XXX
	int free_sv;		// 1 if svm_model is created by svm_load_model
				// 0 if svm_model is created by svm_train
};


#endif /* _LIBSVM_H */
