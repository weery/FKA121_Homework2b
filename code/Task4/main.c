#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../rng_gen.h"
#include "../helper.h"

// DEFINES
#define d_param 0.8
#define nbr_of_dimensions 3


// HELPER FUNCTION DECLARATIONS
 // -----------------------------------------------------------------
double trial_wave(double*, double*, double);
double local_energy(double*, double*, double);
void   new_configuration(double *, double*);
double montecarlo(int, int, double(), double(), double);
double montecarlo_trial(int,int,double(),double(),double);
double montecarlo_E(int,int,double(),double(),double);
double density_probability(double,double);

void montecarlo_E_trial(int,int,double(),double(),double(),double,double*);

//void montecarlo_E_trial(int N ,int equilibrium_time,double (*local_e)(double*,double*,double),double(*gradient)(double*,double*,double,int),double(*prob)(double*,double*,double),double alpha,double* out)

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
    double alpha_0;



    // Task Specific parameters
    double A;
    double beta;

    int nbr_of_trials       =   100000;
    int nbr_of_trials_eq    =   30000;
    int max_p               =   1000;
    int nbr_of_runs         =   5;

    #define alpha_p(i,a,p) (alpha_p_arr[i*max_p*2+a*max_p+p])
    double* alpha_p_arr = (double*)malloc(2*max_p*nbr_of_runs*sizeof(double));

    // Initialize Variables
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;
    alpha_0 = 0.1;
    A       = 1;
    beta    = 0.75;

    for (int i = 0; i < nbr_of_runs; i++)
    {
        double current_alpha = alpha_0;
        printf("New Simulation\n");
        for (int p = 0; p < max_p; p++)
        {
            double* output = (double*)malloc(3*sizeof(double));
            montecarlo_E_trial(nbr_of_trials,nbr_of_trials_eq,local_energy,gradient_alpha, trial_wave, current_alpha,output);

            double current_energy = output[0];
            double current_gradie = output[1];
            double current_produc = output[2];

            double grad_energy = 2*(current_produc-current_energy*current_gradie);

            double step = step_length(A,p+1,beta);
            alpha_p(i,0,p)=current_alpha;
            alpha_p(i,1,p)=current_energy;
            //printf("%e \t %e \n", alpha_p(i,0,p), alpha_p(i,1,p));
            current_alpha -= step*grad_energy;
        }
    }

    FILE* file;

    file = fopen("alpha.dat","w");
    for (int p = 0; p < max_p; p++)
    {
        for (int i = 0; i < nbr_of_runs; i++)
        {
            fprintf(file,"%e \t %e \t", alpha_p(i,0,p), alpha_p(i,1,p));
        }
        fprintf(file, "\n");
    }
    fclose(file);

    free(alpha_p_arr);
    // Free the gsl random number generator
    Free_Generator();
    return 0;
}

// HELPER FUNCTION DEFINITIONS
// ------------------------------------------------------------------
double trial_wave(double* r_1, double* r_2, double alpha)
{
    double r1 = array_abs(r_1,nbr_of_dimensions);
    double r2 = array_abs(r_2,nbr_of_dimensions);
    double r_12[nbr_of_dimensions];
    array_diff(r_1,r_2,r_12, nbr_of_dimensions);
    double r12 = array_abs(r_12,nbr_of_dimensions);

    double f_val = exp(-2*r1) * exp(-2*r2) * exp(r12/(2*(1+alpha*r12)));

    return f_val;
}

double local_energy(double* r_1, double* r_2, double alpha)
{
    double r1 = array_abs(r_1,nbr_of_dimensions);
    double r_12[nbr_of_dimensions];
    double r2 = array_abs(r_2,nbr_of_dimensions);
    array_diff(r_1,r_2,r_12,nbr_of_dimensions);
    double r12 = array_abs(r_12,nbr_of_dimensions);

    double* unit_1 = malloc(sizeof(double)*nbr_of_dimensions);
    double* unit_2 = malloc(sizeof(double)*nbr_of_dimensions);
    array_scalar(unit_1, r_1, 1.0/r1,nbr_of_dimensions);
    array_scalar(unit_2, r_2, 1.0/r2,nbr_of_dimensions);
    double unit_diff[nbr_of_dimensions] ;
    array_diff(unit_1,unit_2,unit_diff,nbr_of_dimensions);

    double term1 = array_mult(unit_diff,r_12,nbr_of_dimensions)/(r12*(1+alpha*r12)*(1+alpha*r12));
    double term2 = - 1/(r12*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    double term3 = - 1/(4*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12)*(1+alpha*r12));
    double term4 = 1/r12;

    double energy = 0;
    energy = -4 +term1+term2+term3+term4;

    free(unit_1);
    free(unit_2);

    return energy;
}

void new_configuration(double * r_1, double* r_2)
{
    for (int i = 0; i < nbr_of_dimensions; i++)
    {
        r_1[i]+= d_param*(randq()-0.5);
        r_2[i]+= d_param*(randq()-0.5);
    }
}


void montecarlo_E_trial(int N ,int equilibrium_time,double (*local_e)(double*,double*,double),double(*gradient)(double*,double*,double,int),double(*prob)(double*,double*,double),double alpha,double* out)
{
    double r_1[nbr_of_dimensions] = { 0 };
    double r_2[nbr_of_dimensions] = { 0 };

    double* energy              = malloc(sizeof(double)*(N-equilibrium_time));
    double* grad_trial          = malloc(sizeof(double)*(N-equilibrium_time));
    double* grad_trial_energy   = malloc(sizeof(double)*(N-equilibrium_time));

    for (int i = 0; i < N; i++)
    {
        // Allocate memory for trial state
        double r_1_new[nbr_of_dimensions];
        double r_2_new[nbr_of_dimensions];

        // Copy values from previous arrays
        memcpy(r_1_new, r_1, nbr_of_dimensions*sizeof(double));
        memcpy(r_2_new, r_2, nbr_of_dimensions*sizeof(double));

        new_configuration(r_1_new, r_2_new);

        double relative_prob = relative_probability(r_1_new,r_2_new,r_1,r_2,alpha,prob,nbr_of_dimensions);

        double r = randq();
        if (relative_prob > r)
        {
            memcpy(r_1, r_1_new, nbr_of_dimensions*sizeof(double));
            memcpy(r_2, r_2_new, nbr_of_dimensions*sizeof(double));
        }
        if (i >= equilibrium_time)
        {
            energy[i-equilibrium_time]      = local_e(r_1,r_2,alpha);
            grad_trial[i-equilibrium_time]  = gradient(r_1,r_2,alpha,nbr_of_dimensions);
            grad_trial_energy[i-equilibrium_time]  = grad_trial[i-equilibrium_time]*energy[i-equilibrium_time];
        }
    }

    double mean_energy              = calc_mean(energy,N-equilibrium_time);
    double mean_grad_trial          = calc_mean(grad_trial,N-equilibrium_time);
    double mean_grad_trial_energy   = calc_mean(grad_trial_energy,N-equilibrium_time);

    free(energy);
    free(grad_trial);

    out[0]=mean_energy;
    out[1]=mean_grad_trial;
    out[2]=mean_grad_trial_energy;
}


double montecarlo_trial(int N ,int equilibrium_time,double(*gradient)(double*,double*,double,int),double(*prob)(double*,double*,double),double alpha)
{
    double r_1[nbr_of_dimensions] = { 0 };
    double r_2[nbr_of_dimensions] = { 0 };

    double* energy = malloc(sizeof(double)*(N-equilibrium_time));


    for (int i = 0; i < N; i++)
    {
        // Allocate memory for trial state
        double r_1_new[nbr_of_dimensions];
        double r_2_new[nbr_of_dimensions];

        // Copy values from previous arrays
        memcpy(r_1_new, r_1, nbr_of_dimensions*sizeof(double));
        memcpy(r_2_new, r_2, nbr_of_dimensions*sizeof(double));

        new_configuration(r_1_new, r_2_new);

        double relative_prob = relative_probability(r_1_new,r_2_new,r_1,r_2,alpha,prob,nbr_of_dimensions);

        double r = randq();
        if (relative_prob > r)
        {
            memcpy(r_1, r_1_new, nbr_of_dimensions*sizeof(double));
            memcpy(r_2, r_2_new, nbr_of_dimensions*sizeof(double));
        }
        if (i >= equilibrium_time)
            energy[i-equilibrium_time] = gradient(r_1,r_2,alpha,nbr_of_dimensions);
    }

    double mean_energy =calc_mean(energy,N-equilibrium_time);

    free (energy);

    return mean_energy;
}


double  montecarlo_E(int N, int equilibrium_time,double (*local_e)(double*,double*,double), double (*f)(double*,double*,double), double alpha)
{
    double r_1[nbr_of_dimensions] = { 0 };
    double r_2[nbr_of_dimensions] = { 0 };

    double* energy = malloc(sizeof(double)*(N-equilibrium_time));


    for (int i = 0; i < N; i++)
    {
        // Allocate memory for trial state
        double r_1_new[nbr_of_dimensions];
        double r_2_new[nbr_of_dimensions];

        // Copy values from previous arrays
        memcpy(r_1_new, r_1, nbr_of_dimensions*sizeof(double));
        memcpy(r_2_new, r_2, nbr_of_dimensions*sizeof(double));

        new_configuration(r_1_new, r_2_new);

        double relative_prob = relative_probability(r_1_new,r_2_new,r_1,r_2,alpha,f,nbr_of_dimensions);


        double r = randq();
        if (relative_prob > r)
        {
            memcpy(r_1, r_1_new, nbr_of_dimensions*sizeof(double));
            memcpy(r_2, r_2_new, nbr_of_dimensions*sizeof(double));
        }
        if (i >= equilibrium_time)
            energy[i-equilibrium_time] = local_e(r_1,r_2,alpha);
    }



    double mean_energy =calc_mean(energy,N-equilibrium_time);

    free (energy);

    return mean_energy;
}

double density_probability(double r, double Z)
{
    double prob=pow(Z,3)*4*pow(r,2)*exp(-2*Z*r);
    return prob;
}
