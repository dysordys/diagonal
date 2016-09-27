#include "common.h"
#include "read_write_matrix.h"
#include "find_eigen.h"
#include "random_number.h"
#include "mutate.h"

int main (int argc, char *argv[]){
  int n = atoi(argv[1]); // number of species
  char * filename = argv[2]; // file storing the matrix
  int num_steps = atoi(argv[3]); // number of steps without improvement before quitting
  int pop_size = atoi(argv[4]); // population size
  int seed = atoi(argv[5]); // random seed
  int k =  atoi(argv[6]); // number of nonzero terms
  double d =  atof(argv[7]); // value for nonzero elements, or ratio
  double ratio = 0.0;

  int i = 0;
  int j = 0;
  int l = 0;

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
  gsl_matrix * tmpM = NULL; // temporary matrix for eigenvalue
  // calculation
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
  gsl_vector_int * bestsol = gsl_vector_int_calloc(n);
  for (i = 0; i < n; i++) gsl_vector_int_set(allcoeffs, i, i);
  gsl_vector_int * pop1[pop_size];
  gsl_vector_int * pop2[pop_size];
  gsl_vector * fitness;
	
  for (i = 0; i < pop_size; i++){
    pop1[i] = gsl_vector_int_calloc(n);
    pop2[i] = gsl_vector_int_calloc(n);
    gsl_ran_choose (r, pop1[i]->data, n, allcoeffs->data, n, sizeof(int));
  }

  double cur_fit;
  double best_fit = l1_original;
  fitness = gsl_vector_calloc(pop_size);	
  int how_many_steps = 0; 
  int counter = 0; 
  int success = 0;
  while(how_many_steps < num_steps){
    how_many_steps++;
    counter++;
    // calculate solutions
    for (i = 0; i < pop_size; i++){
      cur_fit = evaluate_arrangement(pop1[i], k, DD, d, M, tmpM, eval, w);
      gsl_vector_set(fitness, i, cur_fit);
      if (cur_fit < best_fit){
        fprintf(stderr, "Found %e \n", cur_fit);
        best_fit = cur_fit;
        gsl_vector_int_memcpy(bestsol, pop1[i]);
        how_many_steps = 0;
        // We're only interested in finding a stable solution!
        if (best_fit < -0.000001) {
          how_many_steps = num_steps + 1;
          fprintf(stderr, "Terminating...\n");
          fprintf(stdout, "SUCCESS\n");
          success = 1;
          break;
        }
      }
    }
    if (success == 1) break;
    for (l = 0; l < pop_size; l++){
      // reproduction
      // choose parent
      i = gsl_rng_uniform_int(r, pop_size);
      j = gsl_rng_uniform_int(r, pop_size);
      if (gsl_vector_get(fitness, i) < gsl_vector_get(fitness, j)){
        gsl_vector_int_memcpy(pop2[l], pop1[i]);
      } 
      else {
        gsl_vector_int_memcpy(pop2[l], pop1[j]);
      }
      // mutate
      mutate_sol(n, k, pop2[l], r);
    }
    for (l = 0; l < pop_size; l++){
      gsl_vector_int_memcpy(pop1[l], pop2[l]);
    }
  }
	
  fprintf(stderr, "\nFinal with diagonal Re(lambda_1) = %f\n", best_fit);
  for(i = 0; i < k; i++) fprintf(stderr, "%d ", gsl_vector_int_get(bestsol, i)); 
  fprintf(stderr, "\n");
	
  // free matrix M
  gsl_matrix_free(M);
  gsl_matrix_free(DD);
  gsl_vector_int_free(allcoeffs);
  gsl_vector_free(fitness);	
  gsl_vector_int_free(bestsol);
  fprintf(stderr, "1");
  for (i = 0; i < pop_size; i++){
    gsl_vector_int_free(pop1[i]);
    gsl_vector_int_free(pop2[i]);
  }
  fprintf(stderr, "2");
  // free eigenvalue-related variables
  i = eigenvalues_free(&tmpM, &eval, &w);
  // free random number generator
  i = random_free(&r);
  if (success == 0) fprintf(stdout, "COULD NOT FIND STABLE MATRIX WITH k NONZERO\n");
  return 0;
}
