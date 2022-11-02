#include<time.h> 
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<hdf5.h>
#include<fftw3.h> // for finalising the plans

#include "numerical_constants.h"
#include "util.h"
#include "util_hdf5.h"

// hdf5 operation:
herr_t  h5error;
hid_t file_id; // file pointer

struct inputs_def inputs;
struct outputs_def outputs;

int k1, k2, k3;


clock_t start_clock, end_clock;


int main() 
{
	///////////////////////
	// PREPARATION PHASE //
	///////////////////////

	// initialise
	printf("program started\n"); fflush(NULL);
	start_clock = clock(); // the clock	
	Init_constants();
	inputs.Print = Initialise_Printing_struct(); // create printing driver


	// read inputs
	/* to check if exists use printf("link exists 1: %i\n",H5Lexists(file_id, "IRProp/lambda", H5P_DEFAULT)); */
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.
	ReadInputs(file_id, "TDSE_inputs/", &h5error, &inputs);
	Read_1_field_and_grid(file_id, "IRField/", &h5error, &inputs);
	h5error = H5Fclose(file_id); 

	// convert units
	// for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]*1e-15/TIMEau; inputs.Efield.Field[k1] = inputs.Efield.Field[k1]*1e9/EFIELDau;} // convert to atomic units (fs->a.u.), (GV/m->a.u.)

	printf("bconversion\n"); fflush(NULL);
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; /*inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;*/} // convert to atomic units (fs->a.u.), (GV/m->a.u.)
	printf("aconversion\n"); fflush(NULL);

	// Prepare the ground state (it's the state of the atom before the interaction)
	Initialise_grid_and_ground_state(&inputs);


	/////////////////////////
	// COMPUTATIONAL PHASE //
	/////////////////////////

	// calling the TDSE
	printf("call TDSE\n"); fflush(NULL);
	outputs = call1DTDSE(inputs); // THE TDSE
	printf("TDSE done, in the caller\n"); fflush(NULL);


	/////////////////////
	// PROCESS OUTPUTS //
	/////////////////////

	inputs.Print = Set_all_prints(); // set all possible printings for testing

	// store in the hdf5 archive
	file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // we use a different output file to testing, can be changed to have only one file
	hid_t g_id = H5Gcreate2(file_id, "/TDSEsingle", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	h5error = H5Gclose(g_id);	
	PrintOutputs(file_id, "/TDSEsingle/", &h5error, &inputs, &outputs);

	///////////////////////////////////
	// ANALYSES OF THE RUN + TESTING //
	///////////////////////////////////

	// time domain from fftw3 (to check consistency)
	hsize_t output_dims[2];
	output_dims[0] = outputs.Nt; output_dims[1] = 0;
	print_nd_array_h5(file_id, "/TDSEsingle/tgrid_fftw", &h5error, 1, output_dims, outputs.tgrid_fftw, H5T_NATIVE_DOUBLE);
	outputs_destructor(&outputs);
	inputs_destructor(&inputs);

	// durations of calculation
	output_dims[0] = 1; output_dims[1] = 0;	
	end_clock = clock();
	double elapsed_time = (double)((end_clock - start_clock)/CLOCKS_PER_SEC);
	print_nd_array_h5(file_id, "/TDSEsingle/full_runtime", &h5error, 1, output_dims, &elapsed_time, H5T_NATIVE_DOUBLE);
	add_units_1D_h5(file_id, "/TDSEsingle/full_runtime", &h5error, "[s]");


	h5error = H5Fclose(file_id);
	fftw_cleanup();
    printf("Done \n");
}
