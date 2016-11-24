#include <math.h>
#include <stdlib.h> 
#include "./helper.h"

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

double relative_probability(double* r_1, double* r_2,double* R_1, double* R_2, double alpha, double (*f)(double*,double*,int,double), int nbr_of_dimensions)
{
    double trial_1 = f(r_1,r_2,nbr_of_dimensions,alpha);
    double trial_2 = f(R_1,R_2,nbr_of_dimensions,alpha);
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
        arr_out[i]=arr[i]* scalar;
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
