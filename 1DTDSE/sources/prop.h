/**
 * @file prop.h
 * @brief Header for the propagation functions of the TDSE solver.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef PROP
#define PROP

double* propagation(trg_def, Efield_var, double, int, int, double, int, int, double, double *, double *, double *,
				    FILE *, FILE *, double, double, double*, double*, int, int, double, analy_def, outputs_def);
void compute_population(trg_def, Efield_var, int, double *, int, double *, double, double *, double, double, double, double, double, outputs_def);



#endif