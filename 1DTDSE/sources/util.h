#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "structures.h"


// This file contains the variables used by various part of the code.
// It should not contain any external constructions (fftw3,hdf5,mpi,...), so it can be used in any subtask using none of these.


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The field and atomic target:



// functions
double potential(double, trg_def);
double gradpot(double, trg_def);
double norme(double *,int);
void normalise(double *,int);
void Initialise_grid_and_D2(double, int, double **, double **, double **);
void Initialise_grid_and_ground_state(inputs_def *);


void window_analysis(trg_def,double,double,double,int,int,double,double*,double*,double*,double*,double*);
void dipole_analysis(double,double,double*,double*,int,int);
double* extend_grid(double *,int,int,int);



