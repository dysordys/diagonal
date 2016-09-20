extern double find_max_eigen(gsl_matrix * M,
                             gsl_matrix * tmpM,
                             gsl_vector_complex * eval,
                             gsl_eigen_nonsymm_workspace * w);
extern int eigenvalues_setup(int n,
                             gsl_matrix ** tmpM,
                             gsl_vector_complex ** eval,
                             gsl_eigen_nonsymm_workspace ** w);
extern int eigenvalues_free(gsl_matrix ** tmpM,
                            gsl_vector_complex ** eval,
                            gsl_eigen_nonsymm_workspace ** w);
extern double evaluate_arrangement(gsl_vector_int * mysol, int k,
                      gsl_matrix * DD,
                      double d,
                      gsl_matrix * M,
					  gsl_matrix * tmpM,
					  gsl_vector_complex * eval,
					  gsl_eigen_nonsymm_workspace * w);
