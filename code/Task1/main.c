#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rng_gen.h"

// DEFINES
#define nbr_of_dimensions 3


// HELPER FUNCTION DECLARATIONS
 // -----------------------------------------------------------------
double  trial_wave(double*, double*, int, double);
double  array_abs(double*, int);
void    array_diff(double*, double*, int, double*);
double  local_energy(double*, double*, int, double);
double  array_mult(double*, double*, int);
void    new_configuration(double *, double*, int);
void    array_scalar(double*,double*, int, double);
double  montecarlo(int N, double(), int, double);
double  relative_probability(double*, double* ,double* , double* , int , double , double ());
double  mean(double*, int);
double  density_probability(double,double);

// MAIN PROGRAM
// ------------------------------------------------------------------
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



    // Initialize Variables
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;
    alpha   = 0.1;

    double monte = montecarlo(1000000, trial_wave, nbr_of_dimensions, alpha);


    printf("%e \n", monte );

    // Free the gsl random number generator
    Free_Generator();
    return 0;
}

// HELPER FUNCTION DEFINITIONS
// ------------------------------------------------------------------
double relative_probability(double* r_1, double* r_2,double* R_1, double* R_2, int dims, double alpha, double (*f)(double*,double*,int,double))
{
    double trial_1 = f(r_1,r_2,dims,alpha);
    double trial_2 = f(R_1,R_2,dims,alpha);

    trial_1 = abs(trial_1);
    trial_2 = abs(trial_2);

    return trial_1/trial_2;
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

    free(unit_1);
    free(unit_2);

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


double  montecarlo(int N, double (*f)(double*,double*,int,double), int dims, double alpha)
{
    double r_1[dims];
    double r_2[dims];
    new_configuration(r_1,r_2, dims);

    double* energy = malloc(sizeof(double)*N);

    double delta = 0.5;

    for (int i = 0; i < N; i++)
    {
        double r_1_new[dims];
        double r_2_new[dims];
        double r;
        for (int d = 0; d < dims; d++)
        {
            r=randq();
            r_1_new[d] = r_1[d]+delta*(r-0.5);
            r=randq();
            r_2_new[d] = r_2[d]+delta*(r-0.5);
        }

        double relative_prob = relative_probability(r_1_new,r_2_new,r_1,r_2,dims,alpha,f);

        r = randq();
        if (relative_prob > r)
        {
            for (int d = 0; d < dims; d++)
            {
                r_1[d] = r_1_new[d];
                r_2[d] = r_2_new[d];
            }
        }

        energy[i] = local_energy(r_1,r_2,dims,alpha);
    }

    int equilibrium_time= N/10;

    double mean_energy =mean(&energy[equilibrium_time],N-equilibrium_time);

    free (energy);
    return mean_energy;
}

double mean(double* arr, int N)
{
    double sum =0;
    for (int i = 0; i < N; i++)
        sum += arr[i];
    return sum/(double)(N);
}


double density_probability(double r, double Z)
{
    double prob=pow(z,3)*4*pow(r,2)*exp(-2*Z*r);
    return prob;
}
