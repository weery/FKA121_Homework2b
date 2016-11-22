/*
 MD_main.c

 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rng_gen.h"
#define nbr_of_dimensions 3

double trial_wave(double*, double*, int, double);

double array_abs(double*, int);

void array_diff(double*, double*, int, double*);

double local_energy(double*, double*, int, double);

double array_mult(double*, double*, int);

void new_configuration(double *, double*, int);

void array_scalar(double*,double*, int, double);

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

    double r_1[nbr_of_dimensions];
    double r_2[nbr_of_dimensions];

    // Initialize Variables
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;
    alpha   = 0.1;



    new_configuration(r_1,r_2, nbr_of_dimensions);

    double trial = trial_wave(r_1,r_2,nbr_of_dimensions,alpha);

    double energy = local_energy(r_1,r_2,nbr_of_dimensions,alpha);

    printf("%e \t %e \n", trial,energy);
    // Free the gsl random number generator
    Free_Generator();
    return 0;
}

double trial_wave(double* r_1, double* r_2, int dims, double alpha)
{
    double r1 = array_abs(r_1,dims);
    double r2 = array_abs(r_2,dims);
    double r_12[dims];
    array_diff(r_1,r_2,dims,r_12);
    double r12 = array_abs(r_12,dims);

    double f_val = exp(-2*r1)*exp(-2*r2)*exp(r12/(2*(1+alpha*r12)));

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

double local_energy(double* r_1, double* r_2, int dims, double alpha)
{
    double r1 = array_abs(r_1,dims);
    double r2 = array_abs(r_2,dims);
    double r_12[dims];
    array_diff(r_1,r_2,dims,r_12);
    double r12 = array_abs(r_12,dims);

    double* unit_1 = malloc(sizeof(double)*dims);
    double* unit_2 = malloc(sizeof(double)*dims);
    array_scalar(unit_1, r_1,dims,1/r1);
    array_scalar(unit_2, r_2,dims,1/r2);
    double unit_diff[dims] ;
    array_diff(unit_1,unit_2,dims,unit_diff);

    double term1 = array_mult(unit_diff,r_12,dims)/(r12*(1+alpha*r12)*(1+alpha*r12));
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

void new_configuration(double * r_1, double* r_2, int dims)
{
    for (int i = 0; i < dims; i++)
    {
        r_1[i]=randq();
        r_2[i]=randq();
    }
}


void array_scalar(double* arr_out, double* arr , int N , double scalar)
{
    for (int i = 0; i < N; i++)
        arr_out[i]=arr[i]* scalar;
}
