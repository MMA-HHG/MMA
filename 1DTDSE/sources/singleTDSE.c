/**
 * @file singleTDSE.c
 * @brief Contains wrapper around the 1D TDSE solver.
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "constants.h"
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
void call1DTDSE(inputs_def * inputs, outputs_def * outputs)
{
	// Wavefunction in time t_n
	double * psi;
	// Timestep
	double dt;
	// Maximum time
	double tmax;
	// Temporal grid size TDSE
	int Nt;
	// Points per cycle of 800nm field
	int num_t;
	// Dummy integer
	int dumint;
	// Iterable
	int k1;
	// Dummy pointer
	double * dum_ptr;
	// Local field
	double * field;

	
	// local copies of variables given by inputs
	dt = (*inputs).dt;

	// PREPARATIONAL COMPUTATIONS 
	(*inputs).Efield.dt = (*inputs).Efield.tgrid[1]-(*inputs).Efield.tgrid[0]; 
	
	// total length of the grid
	tmax = (*inputs).Efield.tgrid[(*inputs).Efield.Nt-1]-(*inputs).Efield.tgrid[0]; 
   
	// refine dt either by given number of points or by required dt
	if ((*inputs).InterpByDTorNT == 1) {
		k1 = (*inputs).Ntinterp + 1;
	} else { 
		k1 = floor((*inputs).Efield.dt/dt); 
		(*inputs).Ntinterp = k1; 
		if (k1 != 1) k1++; 
	}
	dt = (*inputs).Efield.dt/((double)k1); // redefine dt properly

	// Free field and allocate new field array
	// make the interpolation, note: tgrid does not correspond any more
	field = FourInterp(k1, (*inputs).Efield.Field, (*inputs).Efield.Nt); 
	free((*inputs).Efield.Field);
	Nt = k1*(*inputs).Efield.Nt + 1;
	(*inputs).Efield.Field = field;
	
	// the length of one cycle for 800 nm (i.e. omega=0.057) 
	num_t = floor((2*Pi)/(0.057*dt)); 
	num_t++;  

	// ALLOCATE MEMORY, COPY INITIAL ARRAYS AND PREPARE THEM FOR THE PROPAGATOR
	// Inputs	
	(*inputs).Efield.dt = dt;
	(*inputs).num_t = num_t;
	// Outputs
	(*outputs).tgrid = calloc((Nt+1),sizeof(double));
	(*outputs).Efield = calloc((Nt+1),sizeof(double));
	(*outputs).sourceterm = calloc((Nt+1),sizeof(double));
	(*outputs).PopTot = calloc((Nt+1),sizeof(double));
	(*outputs).PopInt = calloc((Nt+1),sizeof(double));
	(*outputs).expval = calloc((Nt+1),sizeof(double));
	(*outputs).Nt = (Nt+1);

	// do the calculation	
	// Propagate the solution
	psi = propagation(inputs, outputs);

	// Compute FFT
	calcFFTW3(outputs->Nt, dt, tmax, outputs->Efield, &dum_ptr, 
			  &(outputs->FEfield_data), &(outputs->FEfieldM2), &(outputs->Nomega)); 
	free(dum_ptr);
	calcFFTW3(outputs->Nt, dt, tmax, outputs->sourceterm, 
			  &(outputs->omegagrid), &(outputs->Fsourceterm_data), &(outputs->FsourcetermM2), 
			  &(outputs->Nomega));

	free(psi); 
	
}

