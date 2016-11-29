#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./helper.h"

// r phi theta
void to_spherical(double* r, double* s)
{
    s[0] = array_abs(r,3);
    s[1] = atan(r[1]/r[0]);
    s[2] = acos(r[2]/s[0]);
}

double array_abs(double* vals, int nbr_of_dimensions)
{
    double sum  = 0;
    for (int i = 0; i < nbr_of_dimensions; i++)
    {
        sum += vals[i]*vals[i];
    }
    sum = sqrt(sum);
    return sum;
}

double relative_probability(double* r_1, double* r_2,double* R_1, double* R_2, double alpha, double (*f)(double*,double*,double), int nbr_of_dimensions)
{
    double trial_1 = f(r_1,r_2,alpha);
    double trial_2 = f(R_1,R_2,alpha);

    trial_1 = trial_1*trial_1;

    trial_2 = trial_2*trial_2;

    return trial_1/trial_2;
}


void array_diff(double* arr_1, double* arr_2, double* diff, int d)
{
    for (int i = 0; i < d; i++)
        diff[i]=arr_1[i]-arr_2[i];
}

double array_mult(double* arr_1, double* arr_2, int d)
{
    double sum = 0;
    for (int i = 0; i < d; i++)
        sum += arr_1[i]*arr_2[i];
    return sum;
}

void array_scalar(double* arr_out, double* arr, double scalar, int d)
{
    for (int i = 0; i < d; i++)
    {
        arr_out[i]=arr[i]* scalar;
    }
}

double calc_mean(double* arr, int N)
{
    double sum =0;
    for (int i = 0; i < N; i++)
        sum += arr[i];
    return sum/(double)(N);
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

void block_error_estimates(double* data, double* block_error, int N, int b_max)
{
    for (int i = 1; i < b_max; i++)
    {
        double s = block_correlation(data,N,i);
        block_error[i]=s;
    }
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

void to_cartesian(double* s, double* r)
{
    r[0]=s[0]*cos(s[1])*sin(s[2]);
    r[1]=s[0]*sin(s[1])*sin(s[2]);
    r[2]=s[0]*cos(s[2]);
}

double step_length(double A,int p, double beta)
{
    //printf("A: %e\t p: %i\t beta: %e\n",A,p,beta );
    double gamma = A*pow(p,-beta);
    return gamma;
}


double gradient_alpha(double* r_1,double* r_2, double alpha, int nbr_of_dimensions)
{
    double r1 = array_abs(r_1,nbr_of_dimensions);
    double r_12[nbr_of_dimensions];
    double r2 = array_abs(r_2,nbr_of_dimensions);
    array_diff(r_1,r_2,r_12,nbr_of_dimensions);
    double r12 = array_abs(r_12,nbr_of_dimensions);

    double quotient = r12/(1+alpha*r12);
    quotient*=quotient;
    return -quotient/2;
}

// DEBUG TOOL
void print_list(double* data, int N)
{
    for (int i = 0; i < N ; i++)
    {
        printf("%e\t", data[i]);
        if (i%10==9)
            printf("\n");
    }
    printf("\n");
}

double calc_alpha_exp(double x, double min, double max)
{
    return min + (max-min)*(exp(pow(x,3))-1.0)/(exp(1.0) - 1.0);
}
