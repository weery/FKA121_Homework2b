#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../rng_gen.h"
#include "../helper.h"

// DEFINES
#define d_param 1.0
#define nbr_of_dimensions 3


// HELPER FUNCTION DECLARATIONS
 // -----------------------------------------------------------------
double  trial_wave(double*, double*, double);
double  local_energy(double*, double*, double);
void    new_configuration(double *, double*);
void    montecarlo(int, double(), double(), double, double*);
double  density_probability(double,double);

// MAIN PROGRAM
// ------------------------------------------------------------------
int main()
{
    // Initialize the gsl random number generator
    Initialize_Generator();

    // ---- Variable Declarations ----
    // -------------------------------
    double h_bar;
    double e;
    double m_e;
    double e4pi;
    double alpha;

    int nbr_of_trials       =   100000;
    int nbr_of_block_trials =   1000;

    double* energy                = (double*)malloc(nbr_of_trials*sizeof(double));
    double* block_error_estimate  = (double*)malloc(nbr_of_block_trials*sizeof(double));


    // ---- Initialize Variables ----
    // ------------------------------
    h_bar   = 1;
    e       = 1;
    m_e     = 1;
    e4pi    = 1;
    alpha   = 0.1;


    // ========== RUN MONTE CARLO SIMULATION ================
    printf("Start Monte Carlo simulation\n");
    montecarlo(nbr_of_trials,local_energy, trial_wave, alpha,energy);
    printf("Monte Carlo simulaton done\n\n");

    int nbr_of_trials_eq    =   15000;
    block_error_estimates(&energy[nbr_of_trials_eq], block_error_estimate, nbr_of_trials-nbr_of_trials_eq, nbr_of_block_trials);
    double s2= auto_correlation(&energy[nbr_of_trials_eq],nbr_of_trials-nbr_of_trials_eq);

    printf("%.f\n", s2 );
    

    // ------ Print results to file -------
    // ------------------------------------
    FILE* file;
    file = fopen("energy.dat","w");
    for (int i = 0; i < nbr_of_trials; i++)
    {
        fprintf(file, "%e\n", energy[i] );
    }
    fclose(file);

    file = fopen("block_length.dat","w");
    for (int i = 0; i < nbr_of_block_trials; i++)
    {
        fprintf(file, "%e\n", block_error_estimate[i] );
    }
    fclose(file);
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
    double r2 = array_abs(r_2,nbr_of_dimensions);
    double r_12[nbr_of_dimensions];
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


void  montecarlo(int N,double (*local_e)(double*,double*,double), double (*f)(double*,double*,double), double alpha, double* energy)
{
    double r_1[nbr_of_dimensions] = { 0 };
    double r_2[nbr_of_dimensions] = { 0 };
    int rejects = 0;

    r_1[1] = 0.1;
    r_2[1] = -0.1;

    for (int i = 0; i < N; i++)
    {
        // Allocate memory for trial state
        double r_1_new[nbr_of_dimensions];
        double r_2_new[nbr_of_dimensions];

        // Copy values from previous arrays
        memcpy(r_1_new, r_1, nbr_of_dimensions*sizeof(double));
        memcpy(r_2_new, r_2, nbr_of_dimensions*sizeof(double));

        // Generate new configuration
        new_configuration(r_1_new, r_2_new);


        double relative_prob = relative_probability(r_1_new,r_2_new,r_1,r_2,alpha,f,nbr_of_dimensions);

        double r = randq();
        if (relative_prob > r)
        {
            memcpy(r_1, r_1_new, nbr_of_dimensions*sizeof(double));
            memcpy(r_2, r_2_new, nbr_of_dimensions*sizeof(double));
        } else
            rejects++;

        energy[i] = local_e(r_1,r_2,alpha);
    }

    printf("Estimated rejection prob. : %f\n", (double)rejects/N);
}

double density_probability(double r, double Z)
{
    double prob=pow(Z,3)*4*pow(r,2)*exp(-2*Z*r);
    return prob;
}
