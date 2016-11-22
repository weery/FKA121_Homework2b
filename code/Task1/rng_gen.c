// this include is needed, remember to add the gsl library when you compile as in E1.
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>



static gsl_rng *q;

void Initialize_Generator()
{
	gsl_rng_env_setup();
    const gsl_rng_type *T;
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q,time(NULL));
}

double randq()
{
    return gsl_rng_uniform(q);
}

void Free_Generator()
{
    gsl_rng_free(q);
}
