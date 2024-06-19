/**
 * @file strided_scheduler_separated_outs_coarsen_test.c
 * @brief Executes MPI processes, each computing on a part of the field 
 * computed from the CUPRAD code.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <time.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include "constants.h"
#include "tools.h"
#include "tools_hdf5.h"
#include "tools_MPI-RMA.h"
#include "singleTDSE.h"
#include "tools_algorithmic.h"
#include "h5namelist.h"

// hdf5 operation
herr_t h5error;
// file pointer
hid_t file_id; 
// dataspace pointers
hid_t filespace, dataspace_id, dataset_id, dset_id, dspace_id; 

// Input structure
inputs_def inputs;

// Iterable
int k1;

// Clock variables
clock_t start_main, finish1_main;

int main(int argc, char *argv[]) 
{
	// Root name of the temporary HDF5 file, global
	const char filename_stub[] = "hdf5_temp_";
	// Name of the temporary HDF5 file, local
	char local_filename[50];
	// Main HDF5 file name
	char h5_filename[1000];
	get_filename(h5_filename);

	// Dummy variables
	int dum3int[3];
	hsize_t * dims, * dims_input; int ndims; hid_t datatype; // ! hot-fixed to have input dimension different
	char dumchar1[50];
	char dum_h5_paths[3][100];
	// Processing the queue
	// counter of simulations
	int Nsim, Nsim_loc = -1;
	// indices in the Field array
	int kr, kz;

	int comment_operation = 1;
	double t_mpi[10]; 

	////////////////////////
	// PREPARATION PHASE  //
	////////////////////////

	// Initialise MPI
	int myrank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	// Initialise Nsim corresponding to a simulation number for each process:
	// Nsim = - Number_of_processors + process_rank (negative number initially)
	nxtval_init(-nprocs + myrank, &Nsim);

	// Initialise the code values
	Init_constants();
	inputs.Print = Initialise_Printing_struct();
	inputs.Print = Set_all_prints();

	if (comment_operation == 1 ){printf("Proc %i of %i started the program\n",myrank, nprocs);}
	t_mpi[0] = MPI_Wtime();	start_main = clock(); // the clock	

	// create parameters & load initial data
	file_id = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.
	ReadInputs(file_id, CTDSE_INPUTS, GLOBAL_INPUTS, &h5error, &inputs);
	inputs.Print = Set_prints_from_HDF5(file_id, CTDSE_INPUTS, &h5error);
	dims = get_dimensions_h5(file_id, CUPRAD_OUTPUTS_EFIELD, &h5error, &ndims, &datatype);
	dims_input = get_dimensions_h5(file_id, CUPRAD_OUTPUTS_EFIELD, &h5error, &ndims, &datatype);

	// Get dims from the arrays
	hsize_t *dim_t = get_dimensions_h5(file_id, CUPRAD_OUTPUTS_TGRID, &h5error, &ndims, &datatype), \
            *dim_r = get_dimensions_h5(file_id, CUPRAD_OUTPUTS_RGRID, &h5error, &ndims, &datatype), \
            *dim_z = get_dimensions_h5(file_id, CUPRAD_OUTPUTS_ZGRID, &h5error, &ndims, &datatype); 
	// label the dims by physical axes	
    dims[0] = *dim_t; 
	dims[1] = *dim_r; 
	dims[2] = *dim_z;
	dims_input[0] = *dim_z; 
	dims_input[1] = *dim_t; 
	dims_input[2] = *dim_r;

	// Allocate space for the fields & load the tgrid
	inputs.Efield.Field = malloc(((int)dims[0])*sizeof(double));
	inputs.Efield.tgrid =  readreal1Darray_fort(file_id, CUPRAD_OUTPUTS_TGRID, &h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs
	
    /*
	Coarsing procedure:
	Divide the original dimensions of the output laser field by number of
	steps in each dimension. More computationally effective.
	*/
    int kz_step, Nz_max, kr_step, Nr_max;
    readint(file_id, CTDSE_INPUTS_KZ_STEP, &h5error, &kz_step);
    readint(file_id, CTDSE_INPUTS_NZ_MAX,  &h5error, &Nz_max);
    readint(file_id, CTDSE_INPUTS_KR_STEP, &h5error, &kr_step);
    readint(file_id, CTDSE_INPUTS_NR_MAX,  &h5error, &Nr_max);

    // Coarsening: redefine dimensions, t-not affected
    *dim_z = Nz_max/kz_step; 
	*dim_r = Nr_max/kr_step;
    dims[0] = *dim_z; 
	dims[1] = *dim_t; 
	dims[2] = *dim_r;
    
	// Close file access (dataset with the field).
    h5error = H5Fclose(file_id);
    // convert units to atomic units
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++) {
		inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; 
	} 

	// Print information output
	if ((comment_operation == 1) && (myrank == 0)) {
		printf("Proc %i uses dx = %e \n",myrank,inputs.dx);
	}
	if ((comment_operation == 1) && (myrank == 0)) {
		printf("Fields dimensions (z,t,r) = (%llu,%llu,%llu)\n",
			   dims[0],dims[1],dims[2]);
		printf("Fields dimensions (z,t,r) = (%llu,%llu,%llu)\n",
			   dims_input[0],dims_input[1],dims_input[2]);
	}


	// Prepare the ground state (it's the state of the atom before the interaction)
	Initialise_grid_and_ground_state(&inputs);
	
	// Counter - queue length
	int Ntot = (*dim_r)*(*dim_z); 	

	/*
	//////////////////////////
	// COMPUTATIONAL PHASE  //
	//////////////////////////

	First simulation prepares the outputfile (we keep it for the purpose of 
	possible generalisations for parallel output).
	*/
	// Output structure
	outputs_def outputs;

	// Adjusts local simulation number by addition of the number of processors
	// to the initial simulation number Nsim (now non-negative value).
	nxtval_strided(nprocs, &Nsim); 
	// First parallel batch of jobs -> Nsim_loc = 0
	Nsim_loc++;
	t_mpi[1] = MPI_Wtime(); 

	if (Nsim < Ntot){
		double *rgrid_coarse, *zgrid_coarse;
		double *rgrid_CUPRAD, *zgrid_CUPRAD;


		// find proper simulation & load the field
		file_id = H5Fopen (h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		// compute offsets in each dimension
		kr = Nsim % (*dim_r); 
		kz = Nsim - kr;  
		kz = kz/(*dim_r); 
		// coarsen the access	
		dum3int[0] = kz_step*kz; 
		dum3int[1] = -1; 
		dum3int[2] = kr_step*kr;	

		// Read the electric field from the hdf5 file with the indices
		rw_real_fullhyperslab_nd_h5(file_id,CUPRAD_OUTPUTS_EFIELD,&h5error,3,dims_input,dum3int,inputs.Efield.Field,"r");

		int Nz_CUPRAD, Nr_CUPRAD;
		rgrid_CUPRAD = readreal1Darray_fort(file_id, CUPRAD_OUTPUTS_RGRID, &h5error, &Nr_CUPRAD);
		zgrid_CUPRAD = readreal1Darray_fort(file_id, CUPRAD_OUTPUTS_ZGRID, &h5error, &Nz_CUPRAD);

		// Kill the program if Nr_max > Nr_CUPRAD
		if (Nr_max > Nr_CUPRAD) {
			printf("'Nr_max' must be smaller than 'numerics_number_of_points_in_r in the fields array input!");
			free(dims);
			free(dims_input);
			free(rgrid_CUPRAD);
			free(zgrid_CUPRAD);
			inputs_destructor(&inputs);
			exit(1);
		}

		// Close HDF5 file
		h5error = H5Fclose(file_id);
		
		// convert units
		for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++) {
			inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;
		}

		// do the TDSE calculation
		call1DTDSE(&inputs, &outputs); 

		// resize grids
		int Nr_coarse, Nz_coarse;
		coarsen_grid_real(rgrid_CUPRAD, &rgrid_coarse, &Nr_coarse, kr_step, Nr_max);
		coarsen_grid_real(zgrid_CUPRAD, &zgrid_coarse, &Nz_coarse, kz_step, Nz_max);

		// create local output file
		// set name
		local_filename[0] = '\0'; 
		dumchar1[0] = '\0'; 
		sprintf(dumchar1, "%07d", myrank);
		strcat(local_filename,filename_stub); 
		strcat(local_filename,dumchar1); 
		strcat(local_filename,".h5");	

		// Create a new temporary HDF5 file	
		file_id = H5Fcreate (local_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		// Allocate space within temporary HDF5 file
		prepare_local_output_fixed_print_grids_h5(file_id, "", &h5error, &inputs, &outputs, Ntot/nprocs + 1, dims);
		// Write output into the temporary HDF5 file
		print_local_output_fixed_h5(file_id,"", &h5error, &inputs, &outputs, Ntot/nprocs + 1, Nsim, Nsim_loc);

		// print coarser grids
		hsize_t output_dims[2];
		output_dims[0] = Nr_coarse;
		print_nd_array_h5(file_id, "rgrid_coarse", &h5error, 1, output_dims, rgrid_coarse, H5T_NATIVE_DOUBLE);
		output_dims[0] = Nz_coarse;
		print_nd_array_h5(file_id, "zgrid_coarse", &h5error, 1, output_dims, zgrid_coarse, H5T_NATIVE_DOUBLE);

		// print GS etc.
		output_dims[0] = inputs.num_r + 1; output_dims[1] = 2;
		print_nd_array_h5(file_id, "xgrid_micro", &h5error, 1, output_dims, inputs.x, H5T_NATIVE_DOUBLE);
		print_nd_array_h5(file_id, "ground_state", &h5error, 2, output_dims, inputs.psi0, H5T_NATIVE_DOUBLE);

		// print soft-Coulomb parameter for reference
		output_dims[0] = 1;
    	print_nd_array_h5(file_id, "trg_a", &h5error, 1, output_dims, &(inputs.trg.a), H5T_NATIVE_DOUBLE);

		// Close .h5 file
		h5error = H5Fclose(file_id); 
		// clean outputs
		outputs_destructor(&outputs); 

		// Free memory
		free(rgrid_coarse); 
		free(zgrid_coarse);
		free(rgrid_CUPRAD); 
		free(zgrid_CUPRAD);
		
	}
	t_mpi[2] = MPI_Wtime(); 


	// process the MPI queue, update Nsim number for each process
	nxtval_strided(nprocs, &Nsim); 
	// update the number of simulations
	Nsim_loc++;

	if ((comment_operation == 1) && (myrank == 0)) { 
		printf("Proc %i c %i (%f %%); time %f sec \n", myrank, Nsim, 100.0*Nsim/Ntot, t_mpi[2]-t_mpi[1]); 
		fflush(NULL);
	}
	
	// run till queue is not treated
	while (Nsim < Ntot) { 
		t_mpi[3] = MPI_Wtime(); 

		// compute offsets in each dimension
		kr = Nsim % (*dim_r); 
		kz = Nsim - kr;  
		kz = kz/(*dim_r); 

		// prepare the part in the arrray to r/w
		// coarsen the access
 		dum3int[0] = kz_step*kz; 
		dum3int[1] = -1; 
		dum3int[2] = kr_step*kr;	

		//dims_input[0] = dim_t;

		// Alloc input field array again – it has been reallocated within call1DTDSE()
		//inputs.Efield.Field = malloc(((int)(*dim_t))*sizeof(double));

		// read the HDF5 file
		file_id = H5Fopen(h5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		rw_real_fullhyperslab_nd_h5(file_id,CUPRAD_OUTPUTS_EFIELD, &h5error, 3,
									dims_input, dum3int, inputs.Efield.Field, "r");
		h5error = H5Fclose(file_id);

		// convert units
		for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++) {
			inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;
		}

		// do the calculation
		call1DTDSE(&inputs, &outputs); // THE TDSE  

		// open file again
		file_id = H5Fopen(local_filename, H5F_ACC_RDWR, H5P_DEFAULT); 
		print_local_output_fixed_h5(file_id, "", &h5error, &inputs, &outputs, 
									Ntot/nprocs + 1, Nsim, Nsim_loc);
		// close file
		h5error = H5Fclose(file_id); 
		
		// free memory
		outputs_destructor(&outputs);
		nxtval_strided(nprocs,&Nsim); 
		Nsim_loc++;
		
		t_mpi[4] = MPI_Wtime();
		if ((comment_operation == 1) && (myrank == 0)) { 
			printf("Proc %i c %i; time %f sec, from start %f sec \n", myrank, 
				   Nsim, t_mpi[4]-t_mpi[3], t_mpi[4]-t_mpi[1]); 
			fflush(NULL);
		}
		
		t_mpi[5] = MPI_Wtime();
	}

	free(dims);
	free(dims_input);
	free(dim_r);
	free(dim_t);
	free(dim_z);
	// Free inputs
	inputs_destructor(&inputs);
 
	t_mpi[6] = MPI_Wtime();
	if (comment_operation == 1) {
		printf("Proc %i is going to finish., Total time %f sec.\n",myrank, t_mpi[6]-t_mpi[0]); 
		fflush(NULL);
	}

	// End MPI processes
	MPI_Finalize();
	return 0;	
}
