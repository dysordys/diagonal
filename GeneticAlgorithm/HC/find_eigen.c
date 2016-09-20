#include "common.h"
// Find the largest eigenvalue of the nonsymmetric matrix M
double find_max_eigen(gsl_matrix * M,
					  gsl_matrix * tmpM,
					  gsl_vector_complex * eval,
					  gsl_eigen_nonsymm_workspace * w){
	//fprintf(stderr, "Copying \n");
	// Copy the matrix M, as it will be destroyed
	gsl_matrix_memcpy(tmpM, M);
	//fprintf(stderr, "Eigenvalues \n");
	// Find the eigenvalues
	gsl_eigen_nonsymm(tmpM, eval, w);
	//fprintf(stderr, "Finished computing \n");
	// return the maximum real part of the eigenvalues
	int n = M->size1;
	double Rel1 = -100000.0;
	int i;
	for (i = 0; i < n; i++){
    	if (GSL_REAL(gsl_vector_complex_get(eval, i)) > Rel1){
        	Rel1 = GSL_REAL(gsl_vector_complex_get(eval, i));
    	}
	} 
	return Rel1;
}

// Setup the environment for computing the eigenvalues
int eigenvalues_setup(int n,
					  gsl_matrix ** tmpM,
					  gsl_vector_complex ** eval,
					  gsl_eigen_nonsymm_workspace ** w){
	// Allocate temporary matrix for calculations
	*tmpM = gsl_matrix_calloc(n, n);
	// Allocate vector for storing eigenvalues
	*eval = gsl_vector_complex_alloc(n);
	// Allocate workspace for eigenvalue calculation
	*w = gsl_eigen_nonsymm_alloc(n);
	return 0;
}

// Free memory associated with eigenvalues
int eigenvalues_free(gsl_matrix ** tmpM,
					 gsl_vector_complex ** eval,
					 gsl_eigen_nonsymm_workspace ** w){
	gsl_matrix_free(*tmpM);
	gsl_vector_complex_free(*eval);
	gsl_eigen_nonsymm_free(*w);
	return 0;
}

double evaluate_arrangement(gsl_vector_int * mysol,
                      int k,
                      gsl_matrix * DD,
                      double d,
                      gsl_matrix * M,
					  gsl_matrix * tmpM,
					  gsl_vector_complex * eval,
					  gsl_eigen_nonsymm_workspace * w){
    gsl_matrix_set_zero(DD);
    int i;
	for(i = 0; i < k; i++){
		gsl_matrix_set(DD, gsl_vector_int_get(mysol, i), gsl_vector_int_get(mysol, i), d);
		//fprintf(stderr, "%d %f\n", gsl_vector_int_get(mysol, i), gsl_matrix_get(DD, gsl_vector_int_get(mysol, i), gsl_vector_int_get(mysol, i)));
	}
	gsl_matrix_add(DD, M);
    return find_max_eigen(DD, tmpM, eval, w);				  
}
