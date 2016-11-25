#ifndef _helper_h_
#define _helper_h_

void    to_spherical(double*,double*);
double  array_abs(double*, int);
double  relative_probability(double*, double* ,double* , double*, double , double (),int);
void    array_diff(double*, double*, double*, int);
double  array_mult(double*, double*, int);
void    array_scalar(double*,double*, double, int);
double  calc_mean(double*, int);
double  calc_var(double*, int);
double  calc_auto_corr(double* , int , int );
void    to_cartesian(double*,double*);
void block_error_estimates(double*, double*, int, int);
double block_correlation(double*, int, int);
double variance_block(double*, int, int);
double auto_correlation(double*, int);

#endif
