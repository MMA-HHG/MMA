/**
 * @file tools.h
 * @brief Header containing additional functions for the main TDSE code.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TOOLS_H
#define TOOLS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "structures.h"


double potential(double, trg_def);
double gradpot(double, trg_def);
double norme(double *,int);
void normalise(double *,int);
void Initialise_grid_and_D2(double, int, double **, double **, double **);
void Initialise_grid_and_ground_state(inputs_def *);
double * projection_analysis_EV(inputs_def inputs, double *psi, int num_E, double dE);
double * window_analysis(inputs_def inputs, double *psi, int num_E, double dE, double Estep, double E_start);
double* extend_grid(double *,int,int,int);
void free_arr(double *);
void get_filename(char *);


#endif


