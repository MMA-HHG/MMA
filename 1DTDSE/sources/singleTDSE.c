/**
 * @file singleTDSE.c
 * @brief Contains wrapper around the 1D TDSE solver.
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "constants.h"
#include "tools_hdf5.h"
#include "singleTDSE.h"
#include "structures.h"
#include "prop.h"
#include "tools_fftw3.h"


/**
 * @brief Wrapper for running TDSE tasks from the input structure. 
 * 
 * @details This wrapper is called within many binaries compiled in the 1D TDSE
 * bundle. It is a function that takes input structure containing all the necessary 
 * variables for running the TDSE code. The output then contains the computed 
 * values. Which values will be written to the output structure depends on the 
 * requirements set within the inputs file.
 *  
 * @param inputs Input structure.
 * @return outputs_def 
 */
outputs_def call1DTDSE(inputs_def inputs)
{
	// Output structure
	outputs_def outputs;	

	// Wavefunction in time t_n
	double * psi;
	// Timestep
	double dt;
	// Maximum time
	double tmax;
	// Temporal grid size CUPRAD
	int Nt;
	// Temporal grid size TDSE
	int num_t;
	// Radial grid size CUPRAD
	int num_r;
	// Switch
	int input0 = 1;
	// Dummy integer
	int dumint;
	// Iterable
	int k1;
	// Time
	clock_t start, finish;
	// Dummy pointer
	double * dumptrs_real[2];

	
	// local copies of variables given by inputs
	num_r = inputs.num_r; 
	dt = inputs.dt;
	Efield_var Efield = inputs.Efield;

	// PREPARATIONAL COMPUTATIONS 
	
	// find dt from the grid around 0	
	switch (input0) {
	case 0: 
		dumint = 0; 
		break; 
	case 1: 
		dumint = round(Efield.Nt/2.); /* field centered around 0 */ 
		break;
	} // choosing the best resolution	

	Efield.dt = Efield.tgrid[dumint+1]-Efield.tgrid[dumint]; 
	
	// total length of the grid
	tmax = Efield.tgrid[Efield.Nt-1]-Efield.tgrid[0]; 
   
	// refine dt either by given number of points or by required dt
	if (inputs.InterpByDTorNT == 1) {
		k1 = inputs.Ntinterp + 1;
	} else { 
		k1 = floor(Efield.dt/dt); 
		inputs.Ntinterp = k1; 
		k1++; 
	}
	dt = Efield.dt/((double)k1); // redefine dt properly

	// make the interpolation, note: tgrid does not correspond any more
	Efield.Field = FourInterp(k1, Efield.Field, Efield.Nt); 

	Nt = k1*Efield.Nt + 1;
	
	// the length of one cycle for 800 nm (i.e. omega=0.057) 
	num_t = floor((2*Pi)/(0.057*dt)); 
	num_t++;  


	// ALLOCATE MEMORY, COPY INITIAL ARRAYS AND PREPARE THEM FOR THE PROPAGATOR
	// Inputs
	inputs.x = malloc((num_r+1)*sizeof(double)); 
	inputs.psi0 = malloc(2*(num_r+1)*sizeof(double));	
	// Outputs
	outputs.tgrid = calloc((Nt+1),sizeof(double));
	outputs.Efield = calloc((Nt+1),sizeof(double));
	outputs.sourceterm = calloc((Nt+1),sizeof(double));
	outputs.PopTot = calloc((Nt+1),sizeof(double));
	outputs.PopInt = calloc((Nt+1),sizeof(double));
	outputs.expval = calloc((Nt+1),sizeof(double));
	outputs.Nt = (Nt+1);

	// do the calculation
	start = clock();

	// Propagate the solution
	psi = propagation(inputs, outputs);
	
	finish = clock();

	// Compute FFT
	calcFFTW3(outputs.Nt, dt, tmax, outputs.Efield, &(dumptrs_real[0]), &(dumptrs_real[1]), 
			  &outputs.FEfield_data, &outputs.FEfieldM2, &outputs.Nomega); 
	free(dumptrs_real[0]); 
	free(dumptrs_real[1]);
	calcFFTW3(outputs.Nt, dt, tmax, outputs.sourceterm, &outputs.tgrid_fftw, 
			  &outputs.omegagrid, &outputs.Fsourceterm_data, &outputs.FsourcetermM2, 
			  &outputs.Nomega);

	free(inputs.x); 
	free(psi); 
	free(inputs.psi0);
	free(Efield.Field);
	
	return outputs;
}

