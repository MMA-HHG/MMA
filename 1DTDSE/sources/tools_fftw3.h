/**
 * @file tools_fftw3.h
 * @brief Header for FFTW3 routines.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TOOLS_FFTW_H
#define TOOLS_FFTW_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fftw3.h>

double * FourInterp(int , double *, int);
void calcFFTW3(int, double, double, double *, double **, double **, double **, int *);
double ** GaborTransform(double *signal, double dt, int N, int N_freq, int N_t, 
                         double t_min, double t_max, double a);



#endif