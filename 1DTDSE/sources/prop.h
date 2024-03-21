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
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include "structures.h"

double* propagation(inputs_def * inputs, outputs_def * outputs, double * in_field);
void compute_expectation_values(inputs_def *, int, double *, outputs_def *);



#endif