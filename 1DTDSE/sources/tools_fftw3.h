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

double* FourInterp(int , double * , int );
void printFFTW3(FILE *, FILE *, double *, int, double );
void print2FFTW3(FILE *, FILE *, double *, double *, int, double, double);
void print2FFTW3binary(FILE *, FILE *, FILE *, FILE *,FILE *, FILE *,FILE *, FILE *, FILE *, double *, double *, int, double, double);
void printGaborFFTW3(FILE *, FILE *, FILE *, FILE *, double *, int, double, double, double, double);
void printlimitedFFTW3(FILE *, double *, int, double, double, double);
void printGaborFFTW3binary(FILE *, FILE *, FILE *, FILE *, double *, int, double, double, double, double);
void calc2FFTW3(int, double, double, double *, double *, double **, double **, double **, double **, double **, double **, int *);
void calcFFTW3(int, double, double, double *, double **, double **, double **, double **, int *);



#endif