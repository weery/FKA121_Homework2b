/*
 MD_main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rng_gen.h"

double trial_wave(double*, double*, int, double);

double array_abs(double*, int);

void array_diff(double*, double*, int, double*);

double local_energy(double*, double*, int, double);

double array_mult(double*, double*, int);

void new_configuration(double *, double*, int);

/* Main program */
int main()
{
    // Initialize the gsl random number generator
    Initialize_Generator();

    // Variables
    double h_bar;
    double e;
    double m_e;
    double e4pi;
    double alpha;

    double[3] r_1;
    double[3] r_2;

    // Initialize Variables
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;
    alpha   = 0.1;

    new_configuration(r_1,r_2);




    // Free the gsl random number generator
    Free_Generator();
    return 0;
}

double trial_wave(double* r_1, double* r_2, int nbr_of_dimensions, double alpha)
{
    double r1 = array_abs(r_1,nbr_of_dimensions);
    double r2 = array_abs(r_2,nbr_of_dimensions);
    double[3] r_12;
    array_diff(r_1,r_2,nbr_of_dimensions,r12);
    double r_12 = array_abs(r_12,nbr_of_dimensions);

    double f_val = exp(-2*r_1)*exp(-2*r_2)*exp(r12/(2*(1+alpha*r12)));

    return f_val;
}

double array_abs(double* vals, int N)
{
    double sum  = 0;
    for (int i = 0; i < N; i++)
    {
        sum += vals[i]*vals[i];
    }
    sum = sqrt(sum);
    return sum;
}

void array_diff(double* arr_1, double* arr_2, int N, double* diff)
{
    for (int i = 0; i < N; i++)
        diff[i]=arr_1[i]-arr_2[i];
}

double local_energy(double* r_1, double* r_2, int nbr_of_dimensions, double alpha)
{
    double r1 = array_abs(r_1,nbr_of_dimensions);
    double r2 = array_abs(r_2,nbr_of_dimensions);
    double[3] r_12;
    array_diff(r_1,r_2,nbr_of_dimensions,r12);
    double r_12 = array_abs(r_12,nbr_of_dimensions);

    double* unit_1 = r_1/r1;
    double* unit_2 = r_2/r2;
    double[3] unit_diff;
    array_diff(unit_1,unit_2,nbr_of_dimensions,unit_diff);

    double term1 = array_mult(unit_diff,r_12)/(r12*(1+alpha*r12)*(1+alpha*r12));
    double term2 = - 1/(r12*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    double term3 = - 1/(4*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    double term4 = 1/r12;

    double energy = 0;
    energy = -4 +term1+term2+term3+term4;

    return energy;
}


double array_mult(double* arr_1, double* arr_2, int N)
{
    double sum = 0;
    for (int i = 0; i < N ; i++)
        sum += arr_1[i]*arr_2[i];
    return sum;
}

void new_configuration(double * r_1, double* r_2, int nbr_of_dimensions)
{
    for (int i = 0; i < nbr_of_dimensions; i++)
    {
        r_1[i]=randq();
        r_2[i]=randq();
    }
}
