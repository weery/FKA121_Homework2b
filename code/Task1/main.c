#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "rng_gen.h"

// DEFINES
#define d_param 0.1
#define nbr_of_dimensions 3


// HELPER FUNCTION DECLARATIONS
 // -----------------------------------------------------------------
double  trial_wave(double*, double*, double);
double  array_abs(double*);
void    array_diff(double*, double*, double*);
double  local_energy(double*, double*, double);
double  array_mult(double*, double*);
void    new_configuration(double *, double*);
void    array_scalar(double*,double*, double);
double  montecarlo(int N, double(), double, double*);
double  relative_probability(double*, double* ,double* , double*, double , double ());
double  calc_mean(double*, int);
double  calc_var(double*, int);
double  density_probability(double,double);
double calc_auto_corr(double* , int , int );

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


    int nbr_of_trials=1000000;
    double* rads = (double*)malloc(nbr_of_trials*2*sizeof(double));



    // Initialize Variables
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;
    alpha   = 0.1;

    double monte = montecarlo(nbr_of_trials, trial_wave, alpha,rads);

    printf("%e \n", monte );

    FILE* file;
    file = fopen("rads.dat","w");
    for (int i = 0; i < nbr_of_trials; i++)
    {
        fprintf(file, "%e\n", rads[i] );
    }
    fclose(file);
    // Free the gsl random number generator
    Free_Generator();
    return 0;
}

// HELPER FUNCTION DEFINITIONS
// ------------------------------------------------------------------
double relative_probability(double* r_1, double* r_2,double* R_1, double* R_2, double alpha, double (*f)(double*,double*,int,double))
{
    double trial_1 = f(r_1,r_2,nbr_of_dimensions,alpha);
    double trial_2 = f(R_1,R_2,nbr_of_dimensions,alpha);

    return trial_1/trial_2;
}

double trial_wave(double* r_1, double* r_2, double alpha)
{
    double r1 = array_abs(r_1);
    double r2 = array_abs(r_2);
    double r_12[nbr_of_dimensions];
    array_diff(r_1,r_2,r_12);
    double r12 = array_abs(r_12);

    double f_val = exp(-2*r1)*exp(-2*r2)*exp(r12/(2*(1+alpha*r12)));

    return f_val;
}

double array_abs(double* vals)
{
    double sum  = 0;
    for (int i = 0; i < nbr_of_dimensions; i++)
    {
        sum += vals[i]*vals[i];
    }
    sum = sqrt(sum);
    return sum;
}

void array_diff(double* arr_1, double* arr_2, double* diff)
{
    for (int i = 0; i < nbr_of_dimensions; i++)
        diff[i]=arr_1[i]-arr_2[i];
}

double local_energy(double* r_1, double* r_2, double alpha)
{
    double r1 = array_abs(r_1);
    double r2 = array_abs(r_2);
    double r_12[nbr_of_dimensions];
    array_diff(r_1,r_2,r_12);
    double r12 = array_abs(r_12);

    double* unit_1 = malloc(sizeof(double)*nbr_of_dimensions);
    double* unit_2 = malloc(sizeof(double)*nbr_of_dimensions);
    array_scalar(unit_1, r_1, 1.0/r1);
    array_scalar(unit_2, r_2, 1.0/r2);
    double unit_diff[nbr_of_dimensions] ;
    array_diff(unit_1,unit_2,unit_diff);

    double term1 = array_mult(unit_diff,r_12)/(r12*(1+alpha*r12)*(1+alpha*r12));
    double term2 = - 1/(r12*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    double term3 = - 1/(4*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    double term4 = 1/r12;

    double energy = 0;
    energy = -4 +term1+term2+term3+term4;

    free(unit_1);
    free(unit_2);

    return energy;
}


double array_mult(double* arr_1, double* arr_2)
{
    double sum = 0;
    for (int i = 0; i < nbr_of_dimensions; i++)
        sum += arr_1[i]*arr_2[i];
    return sum;
}

void new_configuration(double * r_1, double* r_2)
{
    for (int i = 0; i < nbr_of_dimensions; i++)
    {
        double r1 = randq()-0.5;
        double r2 = randq()-0.5;
        r_1[i] += d_param * r1;
        r_2[i] += d_param * r2;
    }
}


void array_scalar(double* arr_out, double* arr, double scalar)
{
    for (int i = 0; i < nbr_of_dimensions; i++)
        arr_out[i]=arr[i]* scalar;
}


double  montecarlo(int N, double (*f)(double*,double*,int,double), double alpha, double* rads)
{
    double r_1[nbr_of_dimensions];
    double r_2[nbr_of_dimensions];

    r_1[1]=0.1;
    r_2[1]=-0.1;

    double* energy = malloc(sizeof(double)*N);

    double delta = 0.5;

    for (int i = 0; i < N; i++)
    {
        // Allocate memory for trial state
        double r_1_new[nbr_of_dimensions];
        double r_2_new[nbr_of_dimensions];

        // Copy values from previous arrays
        memcpy(r_1_new, r_1, nbr_of_dimensions*sizeof(double));
        memcpy(r_2_new, r_2, nbr_of_dimensions*sizeof(double));

        new_configuration(r_1_new, r_2_new);


        double relative_prob = relative_probability(r_1_new,r_2_new,r_1,r_2,alpha,f);

        double r = randq();
        if (relative_prob > r)
        {
            memcpy(r_1, r_1_new, nbr_of_dimensions*sizeof(double));
            memcpy(r_2, r_2_new, nbr_of_dimensions*sizeof(double));
        }
        rads[2*i]=array_abs(r_1);
        rads[2*i+1]=array_abs(r_2);
        energy[i] = local_energy(r_1,r_2,alpha);
    }

    // Remove the first tenth of the simulation as equilibrium state
    int equilibrium_time= N/10;


    double mean_energy =calc_mean(&energy[equilibrium_time],N-equilibrium_time);

    free (energy);
    return mean_energy;
}

double calc_mean(double* arr, int N)
{
    double sum =0;
    for (int i = 0; i < N; i++)
        sum += arr[i];
    return sum/(double)(N);
}

double density_probability(double r, double Z)
{
    double prob=pow(Z,3)*4*pow(r,2)*exp(-2*Z*r);
    return prob;
}


// Methods that are to be in  Task 2
double calc_var(double* arr, int N)
{
    double var = 0;
    double mean = calc_mean(arr,N);
    for (int i = 0; i < N; i++)
        var += (arr[i]-mean)*(arr[i]-mean);
    var /= (double) N;

    return var;
}




// Calcultes the correlation length of the data
double auto_correlation(double* data, int N)
{
    double decay_value = exp(-2);
    double current_decay = 1;
    int k = 0;
    while (1)
    {
        current_decay = calc_auto_corr(data, N, ++k);
        if (current_decay < decay_value)
            break;
    }
    return k;
}


// Calculates the auto correlation with a given k - correlation length
double calc_auto_corr(double* data, int N, int k)
{
    double mean = calc_mean(data,N);
    double var = calc_var(data,N);

    double corr = 0;
    for (int i = 0; i < N-k;i++)
    {
        corr += data[i]*data[i+k];
    }
    corr/= (double)(N-k);
    return (corr-mean*mean)/var;
}

double variance_block(double* data, int N, int B)
{
    int n_blocks=floor(N/B);
    double* block = malloc((n_blocks)*sizeof(double));

    for (int i = 0; i < n_blocks; i++)
        block[i] = calc_mean(data+(i*B),B);

    double var = 0;

    var = calc_var(block,n_blocks);

    free (block);

    return var;
}

double block_correlation(double* data, int N, int B)
{
    double block_var = variance_block(data,N,B);
    double var = calc_var(data,N);
    double s = (double)(B)*block_var/var;
    return s;
}
