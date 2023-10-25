/**
 * @file single_stride.c
 * @brief Executes single TDSE computing on the output field from the CUPRAD code.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <time.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "tools_hdf5.h"
#include "singleTDSE.h"
#include "structures.h"
#include "tools_algorithmic.h"
#include "tools.h"

// hdf5 operation
herr_t  h5error;
// file pointer
hid_t file_id; 
// dataspace pointers
hid_t filespace, dataspace_id, dataset_id, dset_id, dspace_id;

// Input structure
inputs_def inputs;
// Output structure
outputs_def * outputs;

// Iterable
int k1;

// Clock variables
clock_t start_main, finish_main;


int main() 
{
	// vars:
	const char filename_stub[] = "hdf5_temp_";
	char local_filename[50];

	// dummy
	int dum3int[3];
	hsize_t * dims, * dims_input; 
    int ndims; 
    hid_t datatype; // ! hot-fixed to have input dimension different
	char dumchar1[50];
	// Processing the queue
	// counter of simulations, indices in the Field array
    int kr, kz; 

	int comment_operation = 1;

	////////////////////////
	// PREPARATION PHASE  //
	////////////////////////

	// Initialise the code values
	Init_constants();
	inputs.Print = Initialise_Printing_struct();
	inputs.Print = Set_all_prints();

    // Start the clock	
    start_main = clock(); 

	// Create parameters & load initial data
	file_id = H5Fopen("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); 
	ReadInputs(file_id, "TDSE_inputs/", &h5error, &inputs);
	inputs.Print = Set_prints_from_HDF5(file_id, "TDSE_inputs/", &h5error);
	dims = get_dimensions_h5(file_id, "outputs/output_field", &h5error, &ndims, &datatype);
	dims_input = get_dimensions_h5(file_id, "outputs/output_field", &h5error, &ndims, &datatype);

    // Load dimensions of the output field
	hsize_t dim_t = *get_dimensions_h5(file_id, "outputs/tgrid", &h5error, &ndims, &datatype);
    hsize_t dim_r = *get_dimensions_h5(file_id, "outputs/rgrid", &h5error, &ndims, &datatype);
    hsize_t dim_z = *get_dimensions_h5(file_id, "outputs/zgrid", &h5error, &ndims, &datatype); // label the dims by physical axes	

    dims[0] = dim_t; 
    dims[1] = dim_r; 
    dims[2] = dim_z;
    dims_input[0] = dim_z; 
    dims_input[1] = dim_t; 
    dims_input[2] = dim_r;

    printf("Radial (r) dimension: %llu \n", dims[1]);
    printf("Propagation (z) dimension: %llu \n", dims[2]);

    kr = dims[1];
    kz = dims[2];

    // Prompt the particular field from the output field array
    while (kr >= (int)dims[1] || kr < 0) {
        printf("Set the index in the radial dimension: ");
        scanf("%d", &kr);
    }
    while (kz >= (int)dims[2] || kz < 0) {
        printf("Set the index in the propagation dimension: ");
        scanf("%d", &kz);
    }

	// Allocate space for the fields & load the tgrid
	inputs.Efield.Field = malloc(((int)dims[0])*sizeof(double));
	inputs.Efield.tgrid = readreal1Darray_fort(file_id, "outputs/tgrid",&h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs
	
    // coarsing procedure
    int kz_step, Nz_max, kr_step, Nr_max;
    readint(file_id, "TDSE_inputs/kz_step", &h5error, &kz_step);
    readint(file_id, "TDSE_inputs/Nz_max", &h5error, &Nz_max);
    readint(file_id, "TDSE_inputs/kr_step", &h5error, &kr_step);
    readint(file_id, "TDSE_inputs/Nr_max", &h5error, &Nr_max);

    // redefine dimensions, t-not affected
    dim_z = Nz_max/kz_step; 
    dim_r = Nr_max/kr_step;
    dims[0] = dim_z; 
    dims[1] = dim_t; 
    dims[2] = dim_r;
    
    // Close file access
    h5error = H5Fclose(file_id);

    // convert units to atomic units
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++) {
        inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; 
    } 

    // Print information output
    if (comment_operation == 1) {
        printf("dx = %e \n",inputs.dx);
    }
    if (comment_operation == 1) {
        printf("Fields dimensions (t,r,z) = (%llu,%llu,%llu)\n",dims[0],dims[1],dims[2]);
        printf("Fields dimensions (z,t,r) = (%llu,%llu,%llu)\n",dims_input[0],dims_input[1],dims_input[2]);
    }

	// Prepare the ground state (it's the state of the atom before the interaction)
	Initialise_grid_and_ground_state(&inputs);
    printf("Ground state found.\n");

	//////////////////////////
	// COMPUTATIONAL PHASE  //
	//////////////////////////

    // Computational grids
    double *rgrid_coarse, *zgrid_coarse;
    double *rgrid_CUPRAD, *zgrid_CUPRAD;

    // find proper simulation & load the field
    file_id = H5Fopen("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    
    dum3int[0] = kz; 
    dum3int[1] = -1; 
    dum3int[2] = kr;

    // Read the electric field from the hdf5 file with the indices
    rw_real_fullhyperslab_nd_h5(file_id,"outputs/output_field",&h5error,3,dims_input,dum3int,inputs.Efield.Field,"r");

    // Number of points in the grid
    int Nz_CUPRAD, Nr_CUPRAD;
    rgrid_CUPRAD = readreal1Darray_fort(file_id, "outputs/rgrid", &h5error, &Nr_CUPRAD);
    zgrid_CUPRAD = readreal1Darray_fort(file_id, "outputs/zgrid", &h5error, &Nz_CUPRAD);

    // Close HDF5 file
    h5error = H5Fclose(file_id);

    // convert units
    for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++) {
        inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;
    }

    // do the TDSE calculation
    printf("Starting the computation.\n");
    call1DTDSE(&inputs, outputs); // THE TDSE

    // resize grids
    int Nr_coarse, Nz_coarse;
    coarsen_grid_real(rgrid_CUPRAD, &rgrid_coarse, &Nr_coarse, kr_step, Nr_max);
    coarsen_grid_real(zgrid_CUPRAD, &zgrid_coarse, &Nz_coarse, kz_step, Nz_max);

    // create local output .h5 file
    local_filename[0] = '\0'; 
    dumchar1[0] = '\0'; 
    sprintf(dumchar1, "%07d", 0);
    strcat(local_filename,filename_stub); 
    strcat(local_filename,dumchar1); 
    strcat(local_filename,".h5");		

    // Create a new temporary HDF5 file	
    file_id = H5Fcreate (local_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    prepare_local_output_fixed_print_grids_h5(file_id, "", &h5error, &inputs, outputs, 1, dims);
    print_local_output_fixed_h5(file_id,"", &h5error, &inputs, outputs, 1, 0, 0);

    // print coarser grids
    hsize_t output_dims[2];
    output_dims[0] = Nr_coarse;
    print_nd_array_h5(file_id, "rgrid_coarse", &h5error, 1, output_dims, rgrid_coarse, H5T_NATIVE_DOUBLE);
    output_dims[0] = Nz_coarse;
    print_nd_array_h5(file_id, "zgrid_coarse", &h5error, 1, output_dims, zgrid_coarse, H5T_NATIVE_DOUBLE);

    // print GS etc.
    output_dims[0] = inputs.num_r + 1; 
    output_dims[1] = 2;
    print_nd_array_h5(file_id, "xgrid_micro", &h5error, 1, output_dims, inputs.x, H5T_NATIVE_DOUBLE);
    print_nd_array_h5(file_id, "ground_state", &h5error, 2, output_dims, inputs.psi0, H5T_NATIVE_DOUBLE);

    // Close .h5 file
    h5error = H5Fclose(file_id); 
    // clean outputs
    outputs_destructor(outputs); 
    // clean inputs
    inputs_destructor(&inputs);

    // Free memory
    free(rgrid_coarse); 
    free(zgrid_coarse);
    free(rgrid_CUPRAD); 
    free(zgrid_CUPRAD);
	free(dims);
 
    finish_main = clock();
    printf("Computation finished.\n");
    printf("Program took %f s to complete. \n", (double)(finish_main - start_main)/CLOCKS_PER_SEC);
	return 0;	
}
