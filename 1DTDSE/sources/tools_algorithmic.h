/**
 * @file tools_algorithmic.h
 * @brief Header containing interpolating tools.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TOOLS_ALGORITHMIC_H
#define TOOLS_ALGORITHMIC_H

#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>

void nxtval_init(int, int *);
void nxtval_strided(int, int *);
void coarsen_grid_real(double *, int, double **, int *, int, int);
double interpolate( int , double , double* , double* );
void findinterval(int , double , double* , int* , int* );
double findnextinterpolatedzero(int, double, double* , double* );
double ** create_2Darray_accessor_real(int *, double *);


#endif