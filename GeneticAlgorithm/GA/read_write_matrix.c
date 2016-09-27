#include "common.h"

// read a matrix stored in filename and assign to M
// M is initialized in the process
gsl_matrix * read_matrix(int n, char * filename){
  gsl_matrix * M = gsl_matrix_calloc(n, n);
  FILE * F;
  F = fopen(filename, "rb");
  gsl_matrix_fscanf(F, M);
  fclose(F);
  return M;
}

// print a matrix either to a file or to a stream
int print_matrix(int n, FILE * F, gsl_matrix * M){
  int i, j;
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      fprintf(F, "%f ", gsl_matrix_get(M, i, j));
    }
    fprintf(F, "\n");
  }
  return 0;
}
