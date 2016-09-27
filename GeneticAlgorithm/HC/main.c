#include "common.h"
#include "read_write_matrix.h"
#include "find_eigen.h"
#include "random_number.h"
#include "mutate.h"

int main (int argc, char *argv[]){
  int n = atoi(argv[1]); // number of species
  char * filename = argv[2]; // file storing the matrix
  int num_steps = atoi(argv[3]); // number of steps without improvement before quitting
  int seed = atoi(argv[4]); // random seed
  int k =  atoi(argv[5]); // number of nonzero terms
  double d =  atof(argv[6]); // value for nonzero elements, or ratio
  double ratio = 0.0;

  int i = 0;

  // read in the matrix
  gsl_matrix * M = NULL;
  M = read_matrix(n, filename);
  // divide by max element
  double scaling = gsl_matrix_max(M);
  gsl_matrix_scale(M, 1.0 / scaling);

  // initialize the random number generator
  gsl_rng * r = NULL;
  i = random_setup(&r, seed);

  // setup eigenvalues
  gsl_matrix * tmpM = NULL; // temporary matrix for eigenvalue calculation
  gsl_matrix * DD = gsl_matrix_calloc(n, n); // matrix for diagonal
  gsl_vector_complex * eval = NULL; // vector for storing eigenvalues
  gsl_eigen_nonsymm_workspace * w = NULL; // workspace for eigenvalues
  i = eigenvalues_setup(n, &tmpM, &eval, &w);
	
  // l1 original matrix
  double l1_original = 0.0;
  l1_original = find_max_eigen(M, tmpM, eval, w);
  fprintf(stderr, "Original Re(lambda_1) = %f\n", l1_original);

  // find the right diagonal element
  if (d < 0){
    ratio = -d / l1_original;
  }
  else {
    ratio = d;
    d = - ratio * l1_original;
  }
  fprintf(stderr, "ratio = %f, d = %f\n", ratio, d);
	
  // solutions
  gsl_vector_int * allcoeffs = gsl_vector_int_calloc(n);
  for (i = 0; i < n; i++) gsl_vector_int_set(allcoeffs, i, i);
  gsl_vector_int * cur_sol = gsl_vector_int_calloc(n);
  //gsl_vector_int * mut_sol = gsl_vector_int_calloc(n);
  gsl_vector_int * best_sol = gsl_vector_int_calloc(n);

  double cur_fit;
  //double mut_fit;
  double best_fit;

  // initialize current solution
  gsl_ran_choose (r, cur_sol->data, n, allcoeffs->data, n, sizeof(int));
  //for(i = 0; i < k; i++) fprintf(stderr, "%d ", gsl_vector_int_get(cur_sol, i));
	
  cur_fit = evaluate_arrangement(cur_sol, k, DD, d, M, tmpM, eval, w);
  best_fit = cur_fit;
  gsl_vector_int_memcpy(best_sol, cur_sol);
	
  fprintf(stderr, "\nInitial with diagonal Re(lambda_1) = %f\n", cur_fit);

  int how_many_steps = 0; 
  int counter = 0; 
  int success = 0;
  while(how_many_steps < num_steps){
    how_many_steps++;
    counter++;
    // copy the best solution
    gsl_vector_int_memcpy(cur_sol, best_sol);
    // mutate
    mutate_sol(n, k, cur_sol, r);
    cur_fit = evaluate_arrangement(cur_sol, k, DD, d, M, tmpM, eval, w);
    if (cur_fit < best_fit){
      fprintf(stderr, "Found %e \n", cur_fit);
      best_fit = cur_fit;
      gsl_vector_int_memcpy(best_sol, cur_sol);
      how_many_steps = 0;
      // We're only interested in finding a stable solution!
      if (best_fit < -0.000001) {
        how_many_steps = num_steps + 1;
        fprintf(stderr, "Terminating...\n");
        fprintf(stdout, "SUCCESS\n");
        success = 1;
      }
    }
  }
	
  fprintf(stderr, "\nFinal with diagonal Re(lambda_1) = %f\n", best_fit);
  for(i = 0; i < k; i++) fprintf(stderr, "%d ", gsl_vector_int_get(best_sol, i)); 
  fprintf(stderr, "\n");
	
  // free matrix M
  gsl_matrix_free(M);
  gsl_matrix_free(DD);
  gsl_vector_int_free(allcoeffs);
  gsl_vector_int_free(cur_sol);
  gsl_vector_int_free(best_sol);
  // free eigenvalue-related variables
  i = eigenvalues_free(&tmpM, &eval, &w);
  // free random number generator
  i = random_free(&r);
  if (success == 0) fprintf(stdout, "COULD NOT FIND STABLE MATRIX WITH k NONZERO\n");
  return 0;
}
