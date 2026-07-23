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

/* Private functions for the absorber */
static double smoothstep_absorber(double, double, double, int);

/**
 * @brief Propagates the wavefunction using split operator technique
 * 
 * @param inputs Input parameters.
 * @param outputs Computation results.
 * @param in_field Input field.
 * @return double* wavefunction in final time.
 */
double * propagation(inputs_def *inputs, outputs_def *outputs, double * in_field)
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
	int size;
	if (inputs->analy.writewft != 0) size = Nt/steps_per_dt;
	int num_absorber;
	// arrays for absorbers
	double *absorber_realwf, *absorber_imaginarywf;
	// helpers to speed up calculations
	double *potx; // evaluated potential
	

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

	potx = calloc(num_r+2,sizeof(double)); // one more ghost-point is needed in one matrix

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

	// prepare evaluated potential
	// currently, the code uses the same grid and timesteps for all the calculation, so we can precompute it
	for(j = 0; j <= num_r; j++) {
		potx[j] = potential(x[j],trg); 
	}
	potx[num_r+1] = potential(x[num_r]+dx,trg);


	// prepare absorbing boundaries
	// similarly to potx
	if( (inputs->absorber.type == 1) || (inputs->absorber.type == 2)){

		// find the index corresponding to the absorbing region, shared for absorbers 1 and 2
		findinterval(num_r, x[0]+inputs->absorber.x_cap, x, &num_absorber, &k1);

		// allocate the array for the absorber
		absorber_realwf = calloc(num_absorber+1,sizeof(double));

		if(inputs->absorber.type == 1){ // complex absorber (derived from the complex potential, see the main paper)
			double alpha_vcap = inputs->absorber.alpha;
			double x_rel_distance;
			for(j = 0; j <= num_absorber; j++){
				x_rel_distance = x[j]-inputs->absorber.x_cap-x[0]; // the relative distance of the actual x from the boundary				
				absorber_realwf[j] = exp(-alpha_vcap*x_rel_distance*x_rel_distance*dt); // the absorber
			}
		}else if(inputs->absorber.type == 2){ // smoothstep
			double x_min = x[0];
			double x_max = x[0]+inputs->absorber.x_cap;
			for(j = 0; j <= num_absorber; j++){
				absorber_realwf[j] = smoothstep_absorber(x[j],x_min,x_max,0); // the absorber
			}
		}

		// in both cases, the absorption is the same for the real and the imaginary part of the wavefunction
		absorber_imaginarywf = absorber_realwf;
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
			//if((*inputs).Efield.Field[k]*(*inputs).Efield.Field[k+1] <= 0.0) {
			if(in_field[k]*in_field[k+1] <= 0.0) {
				add_zeros = false;
			}
			Field = 0.;
		} else {
			//Field = (*inputs).Efield.Field[k]; 
			Field = in_field[k]; 
		}


		// Numerov matrix M_2 product with M_2^-1 * (d^2/dx^2 + V)
		for(j = 0 ; j < num_r ; j++) 
		{	
			// Subdiagonal, real and imaginary
			dinfnew1[2*j] = 1/12.; 
			dinfnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potx[j]);
			// Diagonal, real and imaginary
			dnew1[2*j] = 10/12.; 
			dnew1[2*j+1] = 0.5*dt*( 1./(dx*dx) )+0.5*dt*10/12.*(cpot*potx[j]);
			// Superdiagonal, real and imaginary
			dsupnew1[2*j] = 1/12.; 
			dsupnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potx[j+1]);			
		}
		// Last elements of the tridiagonal matric, j = num_r
		// Subdiagonal, real and imaginary
		dinfnew1[2*num_r] = 1/12.; 
		dinfnew1[2*num_r+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potx[num_r]);
		// Diagonal, real and imaginary
		dnew1[2*num_r] = 10/12.; 
		dnew1[2*num_r+1] = 0.5*dt*( 1./(dx*dx) )+0.5*dt*10/12.*(cpot*potx[num_r]);
		// Superdiagonal, real and imaginary, x[num_r] is the final element of the array
		dsupnew1[2*num_r] = 1/12.; 
		dsupnew1[2*num_r +1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potx[num_r+1]);	

		// first part of the evolution (H0+V)
		psi_inter1[0] = (10/12.)*psi[0]+coef*psi[1]+1/12.*psi[2]-0.5*coef*psi[3];
		psi_inter1[0] = psi_inter1[0]+0.5*dt*((10/12.)*psi[1]*(cpot*potx[0])
						+(1/12.)*psi[3]*(cpot*potx[1]));

		psi_inter1[1] = (10/12.)*psi[1]-coef*psi[0]+1/12.*psi[3]+0.5*coef*psi[2];	
		psi_inter1[1] = psi_inter1[1]-0.5*dt*((10/12.)*psi[0]*(cpot*potx[0])
						+(1/12.)*psi[2]*(cpot*potx[1]));

		for(j = 1; j < num_r; j++) {
			psi_inter1[2*j] = (10/12.)*psi[2*j]+coef*psi[2*j+1]+1/12.*psi[2*(j+1)]
							  +1/12.*psi[2*(j-1)]-0.5*coef*(psi[2*(j-1)+1]+psi[2*(j+1)+1]);
			psi_inter1[2*j] = psi_inter1[2*j]+0.5*dt*((10/12.)*psi[2*j+1]*(cpot*potx[j])
							  +(1/12.)*psi[2*(j-1)+1]*(cpot*potx[j-1])
							  +(1/12.)*psi[2*(j+1)+1]*(cpot*potx[j+1]));

			psi_inter1[2*j+1] = (10/12.)*psi[2*j+1]-coef*psi[2*j]+1/12.*psi[2*(j+1)+1]
							  +1/12.*psi[2*(j-1)+1]+0.5*coef*(psi[2*(j-1)]+psi[2*(j+1)]);
			psi_inter1[2*j+1] = psi_inter1[2*j+1]-0.5*dt*((10/12.)*psi[2*j]*(cpot*potx[j])
							  +(1/12.)*psi[2*(j-1)]*(cpot*potx[j-1])
							  +(1/12.)*psi[2*(j+1)]*(cpot*potx[j+1]));
		}

		psi_inter1[2*num_r] = (10/12.)*psi[2*num_r]+coef*psi[2*num_r+1]+1/12.*psi[2*(num_r-1)]-0.5*coef*psi[2*(num_r-1)+1];
		psi_inter1[2*num_r] = psi_inter1[2*num_r]+0.5*dt*((10/12.)*psi[2*num_r+1]*(cpot*potx[num_r])
							  +(1/12.)*psi[2*(num_r-1)+1]*(cpot*potx[num_r-1]));

		psi_inter1[2*num_r+1] = (10/12.)*psi[2*num_r+1]-coef*psi[2*num_r]+1/12.*psi[2*(num_r-1)+1]+0.5*coef*psi[2*(num_r-1)];
		psi_inter1[2*num_r+1] = psi_inter1[2*num_r+1]-0.5*dt*((10/12.)*psi[2*num_r]*(cpot*potx[num_r])
							    +(1/12.)*psi[2*(num_r-1)]*(cpot*potx[num_r-1]));
		
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

		// apply absorption for 1 and 2 options, the factors are already computed, we only multiply here
		if( (inputs->absorber.type == 1) || (inputs->absorber.type == 2)){
			for(j = 0; j <= num_absorber; j++){
				// left side of the interval
				psi[2*j]   *= absorber_realwf[j];
				psi[2*j+1] *= absorber_imaginarywf[j];

				// right side of the interval processed from the right to the left
				// the grid is symmetric, we can use the same factors
				psi[2*(num_r-j)]   *= absorber_realwf[j];
				psi[2*(num_r-j)+1] *= absorber_imaginarywf[j];
			}
		}

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
	free(potx);
	if( (inputs->absorber.type == 1) || (inputs->absorber.type == 2)){
		free(absorber_realwf);
		// absorber_imaginarywf: case is symmentric, the pointers point to the same array, free only once
	}

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

static double smoothstep_absorber(double x, double x_min, double x_max, int direction){
	// rescale to unit interval considering also the direction
	if (direction == 0){
			x = (x-x_min)/(x_max-x_min);
	} else if (direction == 1) {
			x = (x_max-x)/(x_max-x_min);
	} else {
		printf("The smoothstep direction can be only 0 or 1.");
        exit(EXIT_FAILURE);
	}
	return 3.*x*x-2.*x*x*x;	
}

