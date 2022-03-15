//
// svm_model
//
public class svm_model
{
	svm_parameter param;	// parameter
	int nr_class;		// number of classes, = 2 in regression/one class svm
	int l;			// total #SV
	svm_node[][] SV;	// SVs (SV[l])
	double[][] sv_coef;	// coefficients for SVs in decision functions (sv_coef[n-1][l])
	double[] rho;		// constants in decision functions (rho[n*(n-1)/2])

	// for classification only

	int[] label;		// label of each class (label[n])
	int[] nSV;		// number of SVs for each class (nSV[n])
				// nSV[0] + nSV[1] + ... + nSV[n-1] = l
};
