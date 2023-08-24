/**
 * @file singleTDSE.h
 * @brief Header containing 1D TDSE wrapper.
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <time.h>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "structures.h"

outputs_def call1DTDSE(inputs_def inputs);

// Electric field structure
Efield_var Efield;
// Target information - type of gas
trg_def trg;
// ??
analy_def analy;

// Spatial computation grid 'x'
double *x;
// Spatial stepsize
double dx;

double *timet;

double *dipole;
double Eguess,phi,omega,x_int,textend;

int gauge,transformgauge,fieldinau,input0,Ntinterp,InterpByDTorNT;

double dt, tmax,tmin;
int Nt;

// Initial wavefunction \psi, 
double *psi0;
// Wavefunction in time t_n
double *psi;
double E_start,ton,toff,dw;
int num_E,num_exp,num_w,N_t,dumint;


double a_Gabor, omegaMaxGabor, dtGabor, tmin1window, tmax1window, tmin2window, tmax2window, IonFilterThreshold;
int PrintGaborAndSpectrum, IonisationFilterForTheSourceTerm;

int PrintOutputMethod;


int i,num_r,num_t,size,k1,k2;
clock_t start, finish;


FILE *timef,*timef2,*param;
