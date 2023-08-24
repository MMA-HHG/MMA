#include "numerical_constants.h"
#include "util_hdf5.h"
#include "singleTDSE.h"
#include "structures.h"
//#include "util.h"
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
	// declarations
	outputs_def outputs;	
	double * dumptrs_real[2];

	///////////////////////////////////////////////
	// local copies of variables given by inputs //
	///////////////////////////////////////////////

	// input in atomic units

	Eguess = inputs.Eguess; // Energy of the initial state
	num_r = inputs.num_r; // Number of points of the initial spatial grid 16000
	num_exp = inputs.num_exp; // Number of points of the spatial grid for the expansion
	dx = inputs.dx; // resolution for the grid
	InterpByDTorNT = inputs.InterpByDTorNT; // Number of points of the spatial grid for the expansion
	dt = inputs.dt; // resolution in time
	Ntinterp = inputs.Ntinterp; // Number of points of the spatial grid for the expansion
	textend = inputs.textend; // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
	analy.writewft = inputs.analy.writewft; // writewavefunction (1-writting every tprint)
	analy.tprint = inputs.analy.tprint; // time spacing for writing the wavefunction	
	x_int = inputs.x_int; // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
	PrintGaborAndSpectrum = inputs.PrintGaborAndSpectrum; // print Gabor and partial spectra (1-yes)
	a_Gabor = inputs.a_Gabor; // the parameter of the gabor window [a.u.]
	omegaMaxGabor = inputs.omegaMaxGabor; // maximal frequency in Gabor [a.u.]
	dtGabor = inputs.dtGabor; // spacing in Gabor
	tmin1window = inputs.tmin1window; // analyse 1st part of the dipole
	tmax1window = inputs.tmax1window; // analyse 1st part of the dipole
	tmin2window = inputs.tmin2window; // analyse 2nd part of the dipole
	tmax2window = inputs.tmax2window; // analyse 2nd part of the dipole
	PrintOutputMethod = inputs.PrintOutputMethod; // (0 - only text, 1 - only binaries, 2 - both)
	trg.a = inputs.trg.a; // the limit of the integral for the ionisation //2 2 works fine with the
 
	Efield = inputs.Efield;


	// gauge = 0;
	gauge = inputs.gauge;
	transformgauge = 0;
	input0 = 1;
	// printf("tgrid,  %e, %e \n",Efield.tgrid[0],Efield.tgrid[1]);


	////////////////////////////////
	// PREPARATIONAL COMPUTATIONS //
	////////////////////////////////

	// find dt from the grid around 0	
	switch ( input0 ){case 0: dumint = 0; break; case 1: dumint = round(Efield.Nt/2.); /* field centered around 0 */ break;} // choosing the best resolution	
	Efield.dt = Efield.tgrid[dumint+1]-Efield.tgrid[dumint]; // 
	
	tmax = Efield.tgrid[Efield.Nt-1]-Efield.tgrid[0]; // total length of the grid
   
	// refine dt either by given number of points or by required dt
	if (InterpByDTorNT == 1){k1 = Ntinterp + 1;} else { k1 = floor(Efield.dt/dt); Ntinterp = k1; k1++; }
	dt = Efield.dt/((double)k1); // redefine dt properly


	Efield.Field = FourInterp(k1, Efield.Field, Efield.Nt); // make the interpolation !!!!!! tgrid does not correspond any more
	// Efield.Field = dumptrs_real[0]; dumptrs_real[0] = NULL;

	Nt = k1*Efield.Nt + 1; // not very nice to compute it, FourInterp should do it
	
	num_t = floor((2*Pi)/(0.057*dt)); num_t++;  // the length of one cycle for 800 nm (i.e. omega=0.057) 

	size = 2*(num_r+1);// for complex number


	// ALLOCATE MEMORY, COPY INITIAL ARRAYS AND PREPARE THEM FOR THE PROPAGATOR

	x = malloc((num_r+1)*sizeof(double)); // we keep this construction and not use directly initial x due to the extensibility of the grid
	memcpy(x,inputs.x,(num_r+1)*sizeof(double));

	psi0 = malloc(size*sizeof(double));	
	memcpy(psi0,inputs.psi0,size*sizeof(double));

	psi = calloc(size,sizeof(double));
	// t = calloc(Nt,sizeof(double));
	// timet = calloc(Nt,sizeof(double));
	// dipole = calloc(2*Nt,sizeof(double));



	// prepare outputs (there should be written a constructor depending on required values and a destructor on allocated memory should be called in the main code)
	// ineficient allocate every time... Should be done by reference as well.
	// outputs_constructor(outputs,Nt);

	outputs.tgrid = calloc((Nt+1),sizeof(double));
	outputs.Efield = calloc((Nt+1),sizeof(double));
	outputs.sourceterm = calloc((Nt+1),sizeof(double));
	outputs.PopTot = calloc((Nt+1),sizeof(double));

	outputs.PopInt = calloc((Nt+1),sizeof(double));
	outputs.expval = calloc((Nt+1),sizeof(double));

	outputs.Nt = (Nt+1);

	// do the calculation
	start = clock();

	// printf("bprop \n"); fflush(NULL);	
	psi = propagation(trg,Efield,tmin,Nt,num_t,dt,num_r,num_exp,dx,psi0,psi,x,timef,timef2,ton,toff,timet,dipole,gauge,transformgauge,x_int,analy,outputs);
	// printf("TDSE done \n"); fflush(NULL);	
	
	// SAVE THE RESULTS
	calcFFTW3(outputs.Nt, dt, tmax, outputs.Efield, &(dumptrs_real[0]), &(dumptrs_real[1]), &outputs.FEfield_data, &outputs.FEfieldM2, &outputs.Nomega); free(dumptrs_real[0]); free(dumptrs_real[1]);
	calcFFTW3(outputs.Nt, dt, tmax, outputs.sourceterm, &outputs.tgrid_fftw, &outputs.omegagrid, &outputs.Fsourceterm_data, &outputs.FsourcetermM2, &outputs.Nomega);

	free(x); 
	free(psi); 
	free(psi0);
	free(Efield.Field); // this is a tricky free. We pass the inputs "by value", this means that the memory in the original code is unaffected by reallocation. ! BUT it can be freed if coded badly...
	
	return outputs;
}

