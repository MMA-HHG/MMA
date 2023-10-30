/**
 * @file prop.c
 * @brief Contains propagator of the TDSE.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "constants.h"
#include "prop.h"
#include "tools_fftw3.h"
#include "tridiag.h"
#include "structures.h"
#include "tools.h"
#include "tools_algorithmic.h"

/**
 * @brief Propagates the wavefunction using split operator technique
 * 
 * @param inputs Input parameters.
 * @param outputs Computation results.
 * @return double* wavefunction in final time.
 */
double * propagation(inputs_def *inputs, outputs_def *outputs)
{
	// Intermediate and tridiagonal variables
	double *res1, *dnew1, *dinfnew1, *dsupnew1, *psi_inter1;
	// Intermediate and tridiagonal variables
	double *res2, *dnew2, *dinfnew2, *dsupnew2, *psi_inter2;
	double Field, tt = (*inputs).tmin, coef, Apot = 0;
	// Zeroeing the first part of the pulse
	bool add_zeros = true;
	// Iterables
	int j, k, k1, k2, k3, k4;	
	// Potential strength
	double cpot;
	// Electron probability
	double ion_prob2; 

	// Define local variables for the computation
	int num_r = (*inputs).num_r;
	// Spatial grid
	double *x = (*inputs).x;
	// Init wavefunction
	double *psi0 = (*inputs).psi0;
	// Wavefunction
	double *psi;
	// Spatial step
	double dx = (*inputs).dx;
	// Integration limit for ionization 
	double x_int = (*inputs).x_int;
	// Time step
	double dt = (*inputs).Efield.dt;
	// Number of time steps for the TDSE
	int Nt = (*outputs).Nt;
	// Target information
	trg_def trg = (*inputs).trg;
	// Steps for writing wavefunction
	int steps_per_dt = floor(inputs->analy.tprint/dt);
	// Wavefunction write iteration
	int i_wf = 0;
	// Number of stored wavefunctions
	int size = Nt/steps_per_dt;

	// Allocate arrays
	psi = calloc(2*(num_r+1),sizeof(double));
	psi_inter1 = calloc(2*(num_r+1),sizeof(double));
	res1 = calloc(2*(num_r+1),sizeof(double));
	dnew1 = calloc(2*(num_r+1),sizeof(double));
	dinfnew1 = calloc(2*(num_r+1),sizeof(double));
	dsupnew1 = calloc(2*(num_r+1),sizeof(double));
	psi_inter2 = calloc(2*(num_r+1),sizeof(double));
	res2 = calloc(2*(num_r+1),sizeof(double));
	dnew2 = calloc(2*(num_r+1),sizeof(double));
	dinfnew2 = calloc(2*(num_r+1),sizeof(double));
	dsupnew2 = calloc(2*(num_r+1),sizeof(double)); 
	// If write wavefunction, allocate the matrix of wavefunctions in t
	if (inputs->analy.writewft) {
		// Allocate for every timestep t after every 'tprint' time interval
		outputs->psi = malloc(sizeof(double *) * size);
		for (j = 0; j < size; j++) {
			outputs->psi[j] = calloc(2*(num_r+1),sizeof(double));
		}
	}

	// Gauge independent probability of the electron being between -x_int and x_int
	k1 = 0; k2 = 0; k3 = 0; k4 = 0;
	findinterval(num_r, -x_int, x, &k1, &k2);
	findinterval(num_r, x_int, x, &k3, &k4);
	ion_prob2 = 0;
	for(j = k1; j <= k4; j++) {
		ion_prob2 = ion_prob2 + psi0[2*j]*psi0[2*j] + psi0[2*j+1]*psi0[2*j+1];
	}

	Field = 0.;
	cpot = 1.;
	(*outputs).tgrid[0] = tt; 
	(*outputs).sourceterm[0] = 0.; 
	(*outputs).Efield[0] = Field; 
	(*outputs).PopInt[0] = ion_prob2; 
	(*outputs).PopTot[0] = 1.0; 
	(*outputs).expval[0] = 0.0;

	// Save the initial state to psi
	for(j = 0; j <= num_r; j++) {
		psi[2*j] = psi0[2*j]; 
		psi[2*j+1] = psi0[2*j+1];
	}
	


	/************
	 * MAIN LOOP
	*/
	for(k = 0; k < Nt; k++)
	{
		tt = tt + dt;		
		coef = 0.5*dt/(dx*dx);
		
		/* Introduce real zero into the field until the field crosses 0 for
		*  the first time, then the field becomes nonzero. This is useful to
		*  get rid of the initial nonzero pulse front.
		*/
		if (add_zeros) {
			if((*inputs).Efield.Field[k]*(*inputs).Efield.Field[k+1] <= 0.0) {
				add_zeros = false;
			}
			Field = 0.;
		} else {
			Field = (*inputs).Efield.Field[k]; 
		}


		// Numerov matrix M_2 product with M_2^-1 * (d^2/dx^2 + V)
		for(j = 0 ; j < num_r ; j++) 
		{	
			// Subdiagonal, real and imaginary
			dinfnew1[2*j] = 1/12.; 
			dinfnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[j],trg));
			// Diagonal, real and imaginary
			dnew1[2*j] = 10/12.; 
			dnew1[2*j+1] = 0.5*dt*( 1./(dx*dx) )+0.5*dt*10/12.*(cpot*potential(x[j],trg));
			// Superdiagonal, real and imaginary
			dsupnew1[2*j] = 1/12.; 
			dsupnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[j+1],trg));			
		}
		// Last elements of the tridiagonal matric, j = num_r
		// Subdiagonal, real and imaginary
		dinfnew1[2*num_r] = 1/12.; 
		dinfnew1[2*num_r+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[num_r],trg));
		// Diagonal, real and imaginary
		dnew1[2*num_r] = 10/12.; 
		dnew1[2*num_r+1] = 0.5*dt*( 1./(dx*dx) )+0.5*dt*10/12.*(cpot*potential(x[num_r],trg));
		// Superdiagonal, real and imaginary, x[num_r] is the final element of the array
		dsupnew1[2*num_r] = 1/12.; 
		dsupnew1[2*num_r +1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[num_r]+dx,trg));	

		// first part of the evolution (H0+V)
		psi_inter1[0] = (10/12.)*psi[0]+coef*psi[1]+1/12.*psi[2]-0.5*coef*psi[3];
		psi_inter1[0] = psi_inter1[0]+0.5*dt*((10/12.)*psi[1]*(cpot*potential(x[0],trg))
						+(1/12.)*psi[3]*(cpot*potential(x[1],trg)));

		psi_inter1[1] = (10/12.)*psi[1]-coef*psi[0]+1/12.*psi[3]+0.5*coef*psi[2];	
		psi_inter1[1] = psi_inter1[1]-0.5*dt*((10/12.)*psi[0]*(cpot*potential(x[0],trg))
						+(1/12.)*psi[2]*(cpot*potential(x[1],trg)));

		for(j = 1; j < num_r; j++) {
			psi_inter1[2*j] = (10/12.)*psi[2*j]+coef*psi[2*j+1]+1/12.*psi[2*(j+1)]
							  +1/12.*psi[2*(j-1)]-0.5*coef*(psi[2*(j-1)+1]+psi[2*(j+1)+1]);
			psi_inter1[2*j] = psi_inter1[2*j]+0.5*dt*((10/12.)*psi[2*j+1]*(cpot*potential(x[j],trg))
							  +(1/12.)*psi[2*(j-1)+1]*(cpot*potential(x[j-1],trg))
							  +(1/12.)*psi[2*(j+1)+1]*(cpot*potential(x[j+1],trg)));

			psi_inter1[2*j+1] = (10/12.)*psi[2*j+1]-coef*psi[2*j]+1/12.*psi[2*(j+1)+1]
							  +1/12.*psi[2*(j-1)+1]+0.5*coef*(psi[2*(j-1)]+psi[2*(j+1)]);
			psi_inter1[2*j+1] = psi_inter1[2*j+1]-0.5*dt*((10/12.)*psi[2*j]*(cpot*potential(x[j],trg))
							  +(1/12.)*psi[2*(j-1)]*(cpot*potential(x[j-1],trg))
							  +(1/12.)*psi[2*(j+1)]*(cpot*potential(x[j+1],trg)));
		}

		psi_inter1[2*num_r] = (10/12.)*psi[2*num_r]+coef*psi[2*num_r+1]+1/12.*psi[2*(num_r-1)]-0.5*coef*psi[2*(num_r-1)+1];
		psi_inter1[2*num_r] = psi_inter1[2*num_r]+0.5*dt*((10/12.)*psi[2*num_r+1]*(cpot*potential(x[num_r],trg))
							  +(1/12.)*psi[2*(num_r-1)+1]*(cpot*potential(x[num_r-1],trg)));

		psi_inter1[2*num_r+1] = (10/12.)*psi[2*num_r+1]-coef*psi[2*num_r]+1/12.*psi[2*(num_r-1)+1]+0.5*coef*psi[2*(num_r-1)];
		psi_inter1[2*num_r+1] = psi_inter1[2*num_r+1]-0.5*dt*((10/12.)*psi[2*num_r]*(cpot*potential(x[num_r],trg))
							    +(1/12.)*psi[2*(num_r-1)]*(cpot*potential(x[num_r-1],trg)));
		
		// Solve for psi, tridiagonal matrix system
		Inv_Tridiagonal_Matrix_complex(dinfnew1,dnew1,dsupnew1,psi_inter1,res1,num_r+1);

		// second part of the evolution (Hint)
		// Depending on gauge (velocity/length), we apply the corresponding propagator exp(-i V(t))
		if ((*inputs).gauge == 0)
		{
			// Length gauge: exp(-i x * E) = cos(...) - i sin(...)
			for (j = 0; j <= num_r ; j++) 
			{
				psi[2*j] = cos(Field*dt*x[j])*res1[2*j]-sin(Field*dt*x[j])*res1[2*j+1]; 
				psi[2*j+1] = cos(Field*dt*x[j])*res1[2*j+1]+sin(Field*dt*x[j])*res1[2*j];
			}
		}
		else // velocity gauge (Apot has to be available): exp(A * d/dx), derivative approximated using Numerov
		{
			// Tridiagonal matrix init
			for (j = 0; j <= num_r; j++) 
			{			
				dinfnew2[2*j] = 1/6.+0.5*dt*Apot*0.5/dx; 
				dinfnew2[2*j+1] = 0;
				dnew2[2*j] = 4/6.; 
				dnew2[2*j+1] = 0;
				dsupnew2[2*j] = 1/6.-0.5*dt*Apot*0.5/dx; 
				dsupnew2[2*j+1] = 0;
			}

			// RHS vector
			psi_inter2[0] = 4/6.*res1[0]+(1/6.+0.5*dt*Apot*0.5/dx)*res1[2];
			psi_inter2[1] = 4/6.*res1[1]+(1/6.+0.5*dt*Apot*0.5/dx)*res1[3];
			for (j = 1; j < num_r; j++)
			{
				psi_inter2[2*j] = 4/6.*res1[2*j] + (1/6. + 0.5*dt*Apot*0.5/dx)*res1[2*(j+1)];
				psi_inter2[2*j] = psi_inter2[2*j] + (1/6. - 0.5*dt*Apot*0.5/dx)*res1[2*(j-1)];
				psi_inter2[2*j+1] = 4/6.*res1[2*j+1] + (1/6. + 0.5*dt*Apot*0.5/dx)*res1[2*(j+1)+1];
				psi_inter2[2*j+1] = psi_inter2[2*j+1] + (1/6. - 0.5*dt*Apot*0.5/dx)*res1[2*(j-1)+1];

			}
			psi_inter2[2*num_r] = 4/6.*res1[2*num_r]+(1/6.-0.5*dt*Apot*0.5/dx)*res1[2*(num_r-1)];
			psi_inter2[2*num_r+1] = 4/6.*res1[2*num_r+1]+(1/6.-0.5*dt*Apot*0.5/dx)*res1[2*(num_r-1)+1];	

			// Find psi by solving a tridiagonal system
			Inv_Tridiagonal_Matrix_complex(dinfnew2, dnew2, dsupnew2, psi_inter2, res2, num_r+1);
			for(j = 0 ; j<= num_r ; j++) 
			{
				psi[2*j] = res2[2*j]; 
				psi[2*j+1] = res2[2*j+1];	
			}
		}
		(*outputs).tgrid[k+1] = tt;
		(*outputs).Efield[k+1] = Field;
		
		// Compute expectation values: position, current, grad V, population
		compute_expectation_values(inputs, k, psi, outputs);
		
		// Save wavefunction to outputs
		if (inputs->analy.writewft) {
			if (k%steps_per_dt == 0 && i_wf < size-1) {
				for (j = 0; j <= num_r; j++) {
					outputs->psi[i_wf][2*j] = psi[2*j];
					outputs->psi[i_wf][2*j+1] = psi[2*j+1];
				}
				i_wf++;
			}
			// Write the final wavefunction
			if (k == Nt-1 && i_wf == size-1) {
				for (j = 0; j <= num_r; j++) {
					outputs->psi[i_wf][2*j] = psi[2*j];
					outputs->psi[i_wf][2*j+1] = psi[2*j+1];
				}
			}
		}

	} // end of the main loop

	
	free(psi_inter1);
	free(res1);
	free(dnew1);
	free(dinfnew1);
	free(dsupnew1);
    free(psi_inter2);
	free(res2);
	free(dnew2);
	free(dinfnew2);
	free(dsupnew2);

	return psi;
}

/**
 * @brief Computes expectation values of position, current, grad V, electron probability 
 * density and population.
 * 
 * @details Remark that population and current are gauge dependent.
 * 
 * @param inputs Input parameters of the TDSE.
 * @param k Iteration of the main propagation loop.
 * @param psi Wavefunction in time t[k+1].
 * @param outputs Storage of expectation values.
 */
void compute_expectation_values(inputs_def * inputs, int k, double *psi, outputs_def * outputs)
{
	// Average value
	double pop_re, pop_im, pop_tot, current, position, ion_prob2, grad_pot;
	// Iterables
	int j, k1, k2, k3, k4;
	// Grid size
	int num_r = (*inputs).num_r;

	// the population in the ground state (gauge dependent)
	pop_re = 0.; 
	pop_im = 0.;
	for (j = 0 ; j <= num_r ; j++) {
		pop_re = pop_re + psi[2*j]*(*inputs).psi0[2*j] + psi[2*j+1]*(*inputs).psi0[2*j+1]; 
		pop_im = pop_im + psi[2*j]*(*inputs).psi0[2*j+1] - psi[2*j+1]*(*inputs).psi0[2*j];
	}
	pop_tot = pop_re*pop_re + pop_im*pop_im;

	// the gauge independent probability of the electron being between -x_int and x_int
	k1 = 0; k2 = 0; k3 = 0; k4 = 0;
	findinterval(num_r, -(*inputs).x_int, (*inputs).x, &k1, &k2);
	findinterval(num_r, (*inputs).x_int, (*inputs).x, &k3, &k4);
	ion_prob2 = 0;
	for(j = k1; j <= k4; j++) {
		ion_prob2 = ion_prob2 + psi[2*j]*psi[2*j] + psi[2*j+1]*psi[2*j+1];
	}

	// calculation of <x> (gauge independent)
	position = 0.;
	for(j = 0 ; j <= num_r ; j++)
	{
		position = position + (psi[2*j]*psi[2*j] + psi[2*j+1]*psi[2*j+1])*(*inputs).x[j]; 
	}

	// calculation of <grad V> (gauge independent)
	grad_pot = 0.; 
	for (k1 = 0 ; k1 <= num_r ; k1++) {
		grad_pot = grad_pot + (psi[2*k1]*psi[2*k1] + psi[2*k1+1]*psi[2*k1+1])*
				   gradpot((*inputs).x[k1], (*inputs).trg);
	}	
	grad_pot = -grad_pot + (*outputs).Efield[k+1];

	// calculation of current (gauge dependent, gauge independent current is j + A (vector potential))
	current = 0.; 
	for(j = 1; j <= num_r-1; j++) {
		current = current + psi[2*j]*(psi[2*(j+1)+1]-psi[2*(j-1)+1]) - psi[2*j+1]*(psi[2*(j+1)]-psi[2*(j-1)]);                 
	}
	current = current*0.5/(*inputs).dx;

	// save to outputs
	(*outputs).sourceterm[k+1] = grad_pot; 
	(*outputs).PopTot[k+1]=pop_tot;
	(*outputs).expval[k+1]=position;
	(*outputs).PopInt[k+1]=ion_prob2;	
}

