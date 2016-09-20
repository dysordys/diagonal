#include "common.h"

int random_setup(gsl_rng ** r, int seed){
    const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	*r = gsl_rng_alloc(T);
	gsl_rng_set (*r, seed);
	return 0;
}

int random_free(gsl_rng ** r){
    gsl_rng_free(*r);
    return 0;
}
