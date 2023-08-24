/**
 * @file tridiag.h
 * @brief Header containing routines for tridiagonal matrices.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TRIDIAG_H
#define TRIDIAG_H

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "structures.h"

int QL_Tridiagonal_Symmetric_Matrix( double *, double *,double *, int, int);
void Inv_Tridiagonal_Matrix( double *, double *, double *, double *, double *, int);
void Inv_Tridiagonal_Matrix_complex( double *, double *, double *, double *, double *, int );
void Inv_Tridiagonal_Matrix_complex_Numerov( double *, double *, double *, double *, double*, int );
double Einitialise( trg_def, double *,double *,double *,double *,double *,double,double ,int);
double E_calculation(double *,double *,double *,double *,int);
double E_calculation_numerov( trg_def, double *,double ,double *,int );
double Calculate_Shift(double, double, double); 
void Transform_Matrix(double*, double, double, int, int);


#endif