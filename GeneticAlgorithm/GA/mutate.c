#include "common.h"
int mutate_sol(int n, int k, gsl_vector_int * cur_sol, gsl_rng * r){
    int i,j,z;
    i = gsl_rng_uniform_int(r, k);
    j = k + gsl_rng_uniform_int(r, n - k);
    z = gsl_vector_int_get(cur_sol, j);
    gsl_vector_int_set(cur_sol, j, gsl_vector_int_get(cur_sol, i));
    gsl_vector_int_set(cur_sol, i, z);
    return 0;
}
