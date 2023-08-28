/**
 * @file prop.h
 * @brief Header for the propagation functions of the TDSE solver.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef PROP_H
#define PROP_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include "structures.h"

double* propagation(inputs_def inputs, outputs_def outputs);
void compute_population(trg_def, Efield_var, int, double *, int, double *, double, double *, double, double, double, double, double, outputs_def);



#endif