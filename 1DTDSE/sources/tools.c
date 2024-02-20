/**
 * @file tools.c
 * @brief Contains functions for the core TDSE code.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "constants.h"
#include "structures.h"
#include "tridiag.h"
#include "tools_algorithmic.h"
#include "tools.h"

/**
 * @brief Initialises grid, ground state energy and psi0 in input structure.
 * 
 * @details Comment to the choice of the CV criterion for Einitialise: 
 * This number has to be small enough to assure a good conergence of the wavefunction. 
 * If it is not the case, then the scalar product of the the ground state and the excited states 
 * is not quite 0 and those excited states appear in the energy analysis of the ground states. 
 * So the value CV = 1E-25 has been choosen to have a scalar product of 10^-31 with 
 * the third excited state for num_r = 5000 and dx=0.1
 * 
 * @param in Input data for the TDSE.
 */
void Initialise_grid_and_ground_state(inputs_def *in) {
	int k1;
	int size = 2*((*in).num_r+1);
	double *off_diagonal = NULL;
	double *diagonal = NULL;

	// Alloc psi0
	(*in).psi0 = calloc(size,sizeof(double));
	for(k1 = 0; k1 <= (*in).num_r; k1++) {
		(*in).psi0[2*k1] = 1.0; 
		(*in).psi0[2*k1+1] = 0.;
	}
	// Alloc and declare x_grid and diagonal elements for E computation.
	Initialise_grid_and_D2((*in).dx, (*in).num_r, &((*in).x), &diagonal, &off_diagonal); 
	(*in).Einit = Einitialise((*in).trg, (*in).psi0, off_diagonal, diagonal, 
							  off_diagonal, (*in).x, (*in).Eguess, (*in).CV, (*in).num_r); 
	free(diagonal); 
	diagonal = NULL;
	free(off_diagonal);
	off_diagonal = NULL;
}

/**
 * @brief Initialises grid and diagonals.
 * 
 * @param dx Spatial step.
 * @param num_r Number of points in the grid.
 * @param x Grid.
 * @param diagonal Diagonal array.
 * @param off_diagonal Off diagonal array – same for super-/sub-diagonal.
 */
void Initialise_grid_and_D2(double dx, int num_r, double **x, double **diagonal, double **off_diagonal) // Initialise ground-state
{
    int k1;
    double xmax = 0.5*num_r*dx;
	*x = calloc((num_r+1),sizeof(double));
	*off_diagonal = calloc(2*(num_r+1),sizeof(double));
	*diagonal = calloc(2*(num_r+1),sizeof(double));	

	// Declare x grid and diagonals (complex)
	for(k1=0;k1<=num_r;k1++)
	{
		(*x)[k1] = (double)k1 * dx - xmax;
		(*off_diagonal)[2*k1] = -0.5/(dx*dx); 
		(*off_diagonal)[2*k1 + 1] = 0.;
		(*diagonal)[2*k1] = 1./(dx*dx); 
		(*diagonal)[2*k1 + 1] = 0.;
	}	
}

/**
 * @brief Computes potential of the Hamiltonian.
 * 
 * @details See the choice of different soft-core potentials in 
 * https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.023401
 * 
 * @param x Position on the grid.
 * @param trg Target specification.
 * @return double 
 */
double potential(double x, trg_def trg)
{
	// Soft core potential V = -1/sqrt(x^2 + a^2)
	return -1.0/sqrt(trg.a*trg.a+x*x);
	//return -0.5*trg.a/sqrt(x*x + 1/(4*trg.a*trg.a));
}

/**
 * @brief Computes gradient of potential from Hamiltonian.
 * 
 * @details See the choice of different soft-core potentials in 
 * https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.023401
 * 
 * @param x Position on the grid.
 * @param trg Target specification.
 * @return double 
 */
double gradpot(double x,  trg_def trg)
{
	// Returns -grad V
  	return x*pow(trg.a*trg.a+x*x,-1.5);
	//return 0.5*trg.a*x*pow(x*x + 1/(4*trg.a*trg.a), -1.5);
}

/**
 * @brief Computes norm of wavefunction.
 * 
 * @param x Wavefunction.
 * @param Num_r Wavefunction size.
 * @return double 
 */
double norme(double *x,int Num_r) {
	int i;
	double sum = 0.;
	for(i = 0; i <= Num_r; i++) {
		sum = sum + x[2*i]*x[2*i] + x[2*i+1]*x[2*i+1];
	}
	return sum;
}

/**
 * @brief Normalises wavefunction.
 * 
 * @param x Wavefunction.
 * @param Num_r Wavefunction size.
 */
void normalise(double *x, int Num_r) {
	int i;
	double norm;
	norm = sqrt(norme(x, Num_r));
	for(i = 0; i <= Num_r; i++) {
		x[2*i] /= norm;
		x[2*i+1] /= norm;
	}
}

/**
 * @brief Extends wavefunction on the grid.
 * 
 * @param pold Old wavefunction.
 * @param size Additional number of points for extension.
 * @param oldsize Old grid size.
 * @param shift Shift of the wavefunction on the grid.
 * @return double* extended wavefunction
 * 
 * @warning Deprecated, needs to be adjusted to work properly with the grid
 * extension.
 */
double * extend_grid(double *pold, int size, int oldsize, int shift)
{
	int i;
	double *pnew;

	pnew = calloc(2*(size+oldsize+1), sizeof(double));

	if ((shift >= 0) && (size >= shift))
	{
		for (i = oldsize + 1; i <= oldsize + size; i++) {
			pnew[2*i] = 0.; 
			pnew[2*i+1] = 0.;
		}

		for (i = 0; i <= oldsize; i++) {
			pnew[2*shift+2*(oldsize-i)] = pold[2*(oldsize-i)]; 
			pnew[2*shift+2*(oldsize-i)+1] = pold[2*(oldsize-i)+1]; 
		}

		for(i = 0; i < shift; i++) {
			pnew[2*i] = 0.; 
			pnew[2*i+1] = 0.;
		}
	}
	else {
		printf("\nExpansion of the grid incorrect : shift <0 or size < shift \n\n");
	}
	
	free(pold);
	pold = NULL;
	return pnew;
}

/**
 * @brief Computes photoelectron spectrum (PES)
 * 
 * @warning Not implemented into the main TDSE code. Accessible only from Python
 * API.
 * 
 * @param inputs Input structure.
 * @param psi Wavefunction for PES computation.
 * @param num_E Number of energy points. 
 * @param dE Integration range.
 * @param Estep Energy step.
 * @param E_start Starting energy for the PES.
 * @return double* 
 */
double * window_analysis(inputs_def inputs, double *psi, int num_E, double dE, double Estep, double E_start)
{	
	// Diagonals and intermediate results
	double *dnew, *dnew2, *dinfnew, *dsupnew, *res, *res2, *psi2, *dinfnew2, *dsupnew2;
	// Photoelectron spectrum for 1 energy value
	double prob;
	// Iterables
	int i, j, k1;
	// Diagonals of the Hamiltonian matrix
	double *diagonal, *off_diagonal;
	// x grid
	double *x;
	// Photoelectron spectrum
	double *PES;

	// Grid step
	double dx = inputs.dx;
	// Number of x points
	int num_r = inputs.num_r;
	// Size of the x grid
	double xmax = 0.5*num_r*dx;

	// Allocation of arrays
	dnew = calloc(2*(num_r+1),sizeof(double));
	dnew2 = calloc(2*(num_r+1),sizeof(double)); 
	dinfnew = calloc(2*(num_r+1),sizeof(double)); 
	dsupnew = calloc(2*(num_r+1),sizeof(double));
	dinfnew2 = calloc(2*(num_r+1),sizeof(double));
	dsupnew2 = calloc(2*(num_r+1),sizeof(double));
	res = calloc(2*(num_r+1),sizeof(double));
	res2 = calloc(2*(num_r+1),sizeof(double));
	psi2 = calloc(2*(num_r+1),sizeof(double));
	off_diagonal = calloc(2*(num_r+1),sizeof(double));
	diagonal = calloc(2*(num_r+1),sizeof(double));	
	x = calloc((num_r+1),sizeof(double));
	PES = calloc(num_E, sizeof(double));

	// Declare diagonals (complex)
	for (k1 = 0; k1 <= num_r; k1++)
	{
		x[k1] = (double)k1 * dx - xmax;
		off_diagonal[2*k1] = -0.5/(dx*dx); 
		off_diagonal[2*k1 + 1] = 0.;
		diagonal[2*k1] = 1./(dx*dx); 
		diagonal[2*k1 + 1] = 0.;
	}

	// Main cycle for PES:
	// Computes PES using the resolvent operator method.
	for (i = 0; i < num_E; i++) {
		for (j = 0; j <= num_r; j++) 
		{	
			dnew[2*j] = 10*(E_start+Estep*i-dE-potential(x[j],inputs.trg))/12. - diagonal[2*j]; 
			dnew[2*j+1] = -diagonal[2*j+1]-10*dE/12;
			dnew2[2*j] = 10*(E_start+Estep*i+dE-potential(x[j],inputs.trg))/12. - diagonal[2*j]; 
			dnew2[2*j+1] = -diagonal[2*j+1]+10*dE/12;
			
			dinfnew[2*j] = (E_start+Estep*i-dE-potential(x[j],inputs.trg))/12. - off_diagonal[2*j]; 
			dinfnew[2*j+1] = -off_diagonal[2*j+1]-dE/12;
			dsupnew[2*j] = (E_start+Estep*i-dE-potential(x[j+1],inputs.trg))/12. - off_diagonal[2*j]; 
			dsupnew[2*j+1] = -off_diagonal[2*j+1]-dE/12;	 

			dinfnew2[2*j] = (E_start+Estep*i+dE-potential(x[j],inputs.trg))/12. - off_diagonal[2*j]; 
			dinfnew2[2*j+1] = -off_diagonal[2*j+1]+dE/12;
			dsupnew2[2*j] = (E_start+Estep*i+dE-potential(x[j+1],inputs.trg))/12. - off_diagonal[2*j]; 
			dsupnew2[2*j+1] = -off_diagonal[2*j+1]+dE/12;
		}

		Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew, dnew, dsupnew, psi, res, num_r);
		
		for(j = 0; j <= num_r; j++) {
			psi2[2*j] = res[2*j];
			psi2[2*j+1] = res[2*j+1];
		}

		Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew2, dnew2, dsupnew2, psi2, res2, num_r);

		prob = norme(res2,num_r);
		prob = prob*dx*pow(dE,4.);

		PES[i] = prob;
	}

	// Free arrays
	free(dnew); 
	dnew = NULL;
	free(dnew2); 
	dnew2 = NULL;
	free(dinfnew);
	dinfnew = NULL;
	free(res);
	res = NULL;
	free(dsupnew);
	dsupnew = NULL;
	free(dinfnew2);
	dinfnew2 = NULL;
	free(dsupnew2);
	dsupnew2 = NULL;
	free(res2);
	res2 = NULL;
	free(psi2); 
	psi2 = NULL;
	free(diagonal);
	diagonal = NULL;
	free(off_diagonal);
	off_diagonal = NULL;
	free(x);
	x = NULL;

	return PES;
}

/**
 * @brief Does projection of wavefunction onto the ground state
 * 
 * @warning Not implemented into the main TDSE code.
 * 
 * @param inputs Input structure
 * @param psi Wavefunction for the projection
 * @param num_E Number of energy points
 * @param dE Energy step
 * @return double* 
 */
double * projection_analysis_EV(inputs_def inputs, double *psi, int num_E, double dE)
{
	double prob, CV, Eguess, E, ps_re, ps_im;
	// Local wavefunction corresponding to energy E
	double *psi_EV;
	// Diagonals of the Hamiltonian matrix
	double *diagonal, *off_diagonal;
	// x grid
	double *x;
	//double num_E, dE;
	double dx = inputs.dx;
	int num_r = inputs.num_r;
	double xmax = 0.5*num_r*dx;
	double *projection;
	int i, j, k1;
	
	off_diagonal = calloc(2*(num_r+1),sizeof(double));
	diagonal = calloc(2*(num_r+1),sizeof(double));	
	x = calloc((num_r+1),sizeof(double));
	projection = calloc((num_E),sizeof(double));
	psi_EV = calloc(2*(num_r+1),sizeof(double));

	// Declare diagonals (complex)
	for (k1 = 0; k1 <= num_r; k1++)
	{
		x[k1] = (double)k1 * dx - xmax;
		off_diagonal[2*k1] = -0.5/(dx*dx); 
		off_diagonal[2*k1 + 1] = 0.;
		diagonal[2*k1] = 1./(dx*dx); 
		diagonal[2*k1 + 1] = 0.;
	}

	CV = 1E-10; // CV criteria  
		
	for (i = 0; i < num_E; i++)
	{
		// Initialise psi_EV for Einitialise
		for (j = 0; j <= num_r; j++) {
			psi_EV[2*j] = 1; 
			psi_EV[2*j+1] = 0.;
		}
		normalise(psi_EV, num_r); 

	 	Eguess = i*dE; 
		
		E = Einitialise(inputs.trg, psi_EV, off_diagonal, diagonal, 
						off_diagonal, x, Eguess, CV, num_r);

		ps_re = 0; 
		ps_im = 0;
		for (j = 0; j <= num_r; j++)
		{
			ps_re += (psi[2*j]*psi_EV[2*j]-psi[2*j+1]*psi_EV[2*j+1])*x[j];
			ps_im += (psi[2*j]*psi_EV[2*j+1]+psi[2*j+1]*psi_EV[2*j])*x[j];
		}

		projection[i] = ps_re*ps_re + ps_im*ps_im; 

	}

	free(psi_EV);
	psi_EV = NULL;
	free(diagonal);
	diagonal = NULL;
	free(off_diagonal);
	off_diagonal = NULL;
	free(x);
	x = NULL;

	return projection;
}

/**
 * @brief Deletes a 1D array.
 * 
 * @details Used for the Python wrapper - Python TDSE. 
 * 
 * @param buffer Pointer to be freed.
 */
void free_arr(double * buffer) {
	free(buffer);
	buffer = NULL;
}

/**
 * @brief Set the time and field of the input structure
 * 
 * @details Used for the Python wrapper - Python TDSE.
 * 
 * @param in Input structure pointer
 * @param time Time array pointer
 * @param field Field array pointer
 * @param N Size of the arrays
 */
void set_time_and_field(inputs_def * in, double * time, double * field, int N) {
	int i;
	(*in).Efield.tgrid = calloc(N, sizeof(double));
	(*in).Efield.Field = calloc(N, sizeof(double));
	for (i = 0; i < N; i++) {
		(*in).Efield.tgrid[i] = time[i];
		(*in).Efield.Field[i] = field[i];
	}
}

void get_filename(char * filename) {
    FILE *fileStream;

    // Open the file for reading
    if ((fileStream = fopen("msg.tmp", "r")) == NULL) {
        printf("Error! File cannot be opened.\n");
        exit(1);
    }

    // Read the first line until a newline character is encountered
    fscanf(fileStream, "%[^\n]", filename);

    // Close the file
    fclose(fileStream);
}