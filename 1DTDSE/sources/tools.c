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
#include <hdf5.h>
#include "constants.h"
#include "tools_hdf5.h"
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
	free(off_diagonal);
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
 * @param x Position on the grid.
 * @param trg Target specification.
 * @return double 
 */
double potential(double x, trg_def trg)
{
	// Soft core potential V = -1/sqrt(x^2 + a^2)
	return -1.0/sqrt(trg.a*trg.a+x*x);
}

/**
 * @brief Computes gradient of potential from Hamiltonian.
 * 
 * @param x Position on the grid.
 * @param trg Target specification.
 * @return double 
 */
double gradpot(double x,  trg_def trg)
{
	// Returns -grad V
  	return x*pow(trg.a*trg.a+x*x,-1.5);
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
	return pnew;
}

/**
 * @brief Computes photoelectron spectrum.
 * 
 * @param trg 
 * @param dE 
 * @param Estep 
 * @param E_start 
 * @param num_E 
 * @param num_r 
 * @param dx 
 * @param psi 
 * @param dinf 
 * @param d 
 * @param dsup 
 * @param x 
 * 
 * @warning Not implemented into the main TDSE code.
 */
/*void window_analysis(trg_def trg, double dE, double Estep, double E_start, 
					 int num_E, int num_r, double dx, double *psi, double *dinf,
					 double *d, double *dsup, double *x)*/
double * window_analysis(inputs_def inputs, double *psi, int num_E, double dE, double Estep, double E_start)
{	
	double *dnew,*dnew2,*dinfnew,*dsupnew,*res,*res2,*psi2;
	double *dinfnew2,*dsupnew2;
	double prob;
	int i, j;
	//FILE *fel; 
	// Diagonals of the Hamiltonian matrix
	double *diagonal, *off_diagonal;
	// x grid
	double *x;
	double *PES;
	double dx = inputs.dx;
	int num_r = inputs.num_r;
	double xmax = 0.5*num_r*dx;

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
	for(int k1 = 0; k1 <= num_r; k1++)
	{
		x[k1] = (double)k1 * dx - xmax;
		off_diagonal[2*k1] = -0.5/(dx*dx); 
		off_diagonal[2*k1 + 1] = 0.;
		diagonal[2*k1] = 1./(dx*dx); 
		diagonal[2*k1 + 1] = 0.;
	}


	/*fel = fopen("electron_spectrum.dat","w");
	if (fel == NULL) {
		printf("Cannot open electron_spectrum.dat"); 
		exit(1);
	} 
	printf("Working on the bin 00000");
	*/

	for (i = 0; i < num_E; i++) {
		for (j = 0; j <= num_r; j++) 
		{	
			dnew[2*j] = 10*(E_start+Estep*i-dE/sqrt(8)-potential(x[j],inputs.trg))/12. - diagonal[2*j]; 
			dnew[2*j+1] = -diagonal[2*j+1]-10*dE/(12*sqrt(8));
			dnew2[2*j] = 10*(E_start+Estep*i+dE/sqrt(8)-potential(x[j],inputs.trg))/12. - diagonal[2*j]; 
			dnew2[2*j+1] = -diagonal[2*j+1]+10*dE/(12*sqrt(8));
			
			dinfnew[2*j] = (E_start+Estep*i-dE/sqrt(8)-potential(x[j],inputs.trg))/12. - off_diagonal[2*j]; 
			dinfnew[2*j+1] = -off_diagonal[2*j+1]-dE/(12*sqrt(8));
			dsupnew[2*j] = (E_start+Estep*i-dE/sqrt(8)-potential(x[j+1],inputs.trg))/12. - off_diagonal[2*j]; 
			dsupnew[2*j+1] = -off_diagonal[2*j+1]-dE/(12*sqrt(8));	 

			dinfnew2[2*j] = (E_start+Estep*i+dE/sqrt(8)-potential(x[j],inputs.trg))/12. - off_diagonal[2*j]; 
			dinfnew2[2*j+1] = -off_diagonal[2*j+1]+dE/(12*sqrt(8));
			dsupnew2[2*j] = (E_start+Estep*i+dE/sqrt(8)-potential(x[j+1],inputs.trg))/12. - off_diagonal[2*j]; 
			dsupnew2[2*j+1] = -off_diagonal[2*j+1]+dE/(12*sqrt(8));
		}

		Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew, dnew, dsupnew, psi, res, num_r);
		
		for(j = 0; j <= num_r; j++) {
			psi2[2*j] = res[2*j];
			psi2[2*j+1] = res[2*j+1];
		}

		Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew2, dnew2, dsupnew2, psi2, res2, num_r);

		prob = norme(res2,num_r);
		prob = prob*dx*pow(dE,4.);

		//fprintf(fel,"%e\t%e\n",E_start+Estep*i,prob);

		PES[i] = prob;
		
		//printf("\b\b\b\b\b%5d", i); fflush(stdout); 

	}



	free(dnew); 
	free(dnew2); 
	free(dinfnew);
	free(res);
	free(dsupnew);
	free(dinfnew2);
	free(dsupnew2);
	free(res2);
	free(psi2); 
	free(diagonal);
	free(off_diagonal);
	free(x);

	return PES;
}
//Does projection of wavefunction onto the ground state
//@warning Not implemented into the main TDSE code.

double * projection_analysis_EV(inputs_def inputs, double *psi, int num_E, double dE, double E_start) // inti procedure incompatible
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
	//double E_start = -0.6;
	double *projection;
	//FILE *fel;
	int i, j;
	//num_E = 10000;
	//dE = 0.001;	

	off_diagonal = calloc(2*(num_r+1),sizeof(double));
	diagonal = calloc(2*(num_r+1),sizeof(double));	
	x = calloc((num_r+1),sizeof(double));
	projection = calloc((num_E),sizeof(double));
	psi_EV = calloc(2*(num_r+1),sizeof(double));

	// Declare diagonals (complex)
	for(int k1 = 0; k1 <= num_r; k1++)
	{
		x[k1] = (double)k1 * dx - xmax;
		off_diagonal[2*k1] = -0.5/(dx*dx); 
		off_diagonal[2*k1 + 1] = 0.;
		diagonal[2*k1] = 1./(dx*dx); 
		diagonal[2*k1 + 1] = 0.;
	}



	//fel = fopen("electron_spectrum.dat", "w" );
	//printf("Working on the bin 00000");


	CV = 1E-10; // CV criteria  
		
	for(i = 0; i < num_E; i++)
	{
		// Initialise psi_EV for Einitialise
		for(j = 0; j <= num_r; j++) {
			psi_EV[2*j] = 1; 
			psi_EV[2*j+1] = 0.;
		}
		normalise(psi_EV, num_r); 

	 	Eguess = i*dE; 
		
		E = Einitialise(inputs.trg, psi_EV, off_diagonal, diagonal, 
						off_diagonal, x, Eguess, CV, num_r);
		//printf("%e\t%e\n", Eguess, E);

		ps_re = 0; 
		ps_im = 0;
		for(j = 0; j <= num_r; j++)
		{
			ps_re += (psi[2*j]*psi_EV[2*j]-psi[2*j+1]*psi_EV[2*j+1])*x[j];
			ps_im += (psi[2*j]*psi_EV[2*j+1]+psi[2*j+1]*psi_EV[2*j])*x[j];
		}

		projection[i] = ps_re*ps_re + ps_im*ps_im; 
		
		//fprintf(fel,"%e\t%e\n", E, prob); 

	}

	//printf("\n");

	//fclose(fel); 
	free(psi_EV);
	free(diagonal);
	free(off_diagonal);
	free(x);

	return projection;
}


