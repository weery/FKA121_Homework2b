/*
 MD_main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rng_gen.h"

double trial_wave(double*, double*, int);


/* Main program */
int main()
{
    // Initialize the gsl random number generator
    Initialize_Generator();

    // Variables
    double h_bar, e, m_e, e4pi;

    // Initialize Variables
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;



    // Free the gsl random number generator
    Free_Generator();
    return 0;
}

double trial_wave(double* r_1, double* r_2, int nbr_of_particles)
{
    return 0;
}
