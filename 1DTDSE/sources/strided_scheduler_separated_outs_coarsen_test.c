#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>

#include<math.h>
#include "hdf5.h"

#include "numerical_constants.h"
#include "util.h"
#include "util_hdf5.h"
#include "util_mpi.h"

// hdf5 operation:
herr_t  h5error;
hid_t file_id; // file pointer
hid_t filespace, dataspace_id, dataset_id, dset_id, dspace_id; // dataspace pointers

struct inputs_def inputs;
struct outputs_def outputs;

int k1, k2, k3;
int one = 1;


clock_t start_main, finish2_main, finish1_main, finish3_main, finish4_main;


int main(int argc, char *argv[]) 
{
	// vars:

	const char filename_stub[] = "hdf5_temp_";
	char local_filename[50];

	// dummy
	int dum3int[3];
	hsize_t * dims, * dims_input; int ndims; hid_t datatype; // ! hot-fixed to have input dimension different
	char dumchar1[50], dumchar2[50];
	// Processing the queue
	int Nsim, Nsim_loc = -1, kr, kz; // counter of simulations, indices in the Field array

	int comment_operation = 1;
	double t_mpi[10]; 

	// hsize_t output_dims[4];

	////////////////////////
	// PREPARATION PAHASE //
	////////////////////////

	// Initialise MPI
	int myrank, nprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	nxtval_init(-nprocs+myrank,&Nsim);

	// Initialise the code values
	Init_constants();
	inputs.Print = Initialise_Printing_struct();
	inputs.Print = Set_all_prints();

	if (comment_operation == 1 ){printf("Proc %i of %i started the program\n",myrank, nprocs);}
	t_mpi[0] = MPI_Wtime();	start_main = clock(); // the clock	

	// create parameters & load initial data
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.
	ReadInputs(file_id, "TDSE_inputs/", &h5error, &inputs);
	inputs.Print = Set_prints_from_HDF5(file_id, "TDSE_inputs/", &h5error);
	dims = get_dimensions_h5(file_id, "outputs/output_field", &h5error, &ndims, &datatype);
	dims_input = get_dimensions_h5(file_id, "outputs/output_field", &h5error, &ndims, &datatype);

    // *dims = malloc((*ndims)*sizeof(hsize_t))

	hsize_t dim_t = *get_dimensions_h5(file_id, "outputs/tgrid", &h5error, &ndims, &datatype), \
            dim_r = *get_dimensions_h5(file_id, "outputs/rgrid", &h5error, &ndims, &datatype), \
            dim_z = *get_dimensions_h5(file_id, "outputs/zgrid", &h5error, &ndims, &datatype); // label the dims by physical axes	

    dims[0] = dim_t; dims[1] = dim_r; dims[2] = dim_z;
	dims_input[0] = dim_z; dims_input[1] = dim_t; dims_input[2] = dim_r;

	// create space for the fields & load the tgrid
	inputs.Efield.Field = malloc(((int)dims[0])*sizeof(double));
	inputs.Efield.tgrid =  readreal1Darray_fort(file_id, "outputs/tgrid",&h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs
	
    // coarsing procedure
    int kz_step, Nz_max, kr_step, Nr_max;
    readint(file_id, "TDSE_inputs/kz_step", &h5error, &kz_step);
    readint(file_id, "TDSE_inputs/Nz_max", &h5error, &Nz_max);
    readint(file_id, "TDSE_inputs/kr_step", &h5error, &kr_step);
    readint(file_id, "TDSE_inputs/Nr_max", &h5error, &Nr_max);

    // redefine dimensions, t-not affected
    dim_z = Nz_max/kz_step; dim_r = Nr_max/kr_step;
    dims[0] = dim_z; dims[1] = dim_t; dims[2] = dim_r;
    
    
    h5error = H5Fclose(file_id);
    // convert units
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; /*inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;*/} // convert to atomic units (fs->a.u.), (GV/m->a.u.)

	if (( comment_operation == 1 ) && ( myrank == 0 ) ){printf("Proc %i uses dx = %e \n",myrank,inputs.dx);}
	if ( ( comment_operation == 1 ) && ( myrank == 0 ) ){printf("Fields dimensions (t,r,z) = (%i,%i,%i)\n",dims[0],dims[1],dims[2]);
														 printf("Fields dimensions (z,t,r) = (%i,%i,%i)\n",dims_input[0],dims_input[1],dims_input[2]);}


	// Prepare the ground state (it's the state of the atom before the interaction)
	Initialise_grid_and_ground_state(&inputs);

	int Ntot = dim_r*dim_z; // counter (queue length)	


	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////

	// first simulation prepares the outputfile (we keep it for the purpose of possible generalisations for parallel output)
	nxtval_strided(nprocs,&Nsim); Nsim_loc++;
	t_mpi[1] = MPI_Wtime(); 

	if (Nsim < Ntot){

		double *rgrid_coarse, *zgrid_coarse;
		double *rgrid_CUPRAD, *zgrid_CUPRAD;

		// find proper simulation & load the field
		file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dim_r; // compute offsets in each dimension
		dum3int[0]=kz_step*kz; dum3int[1]=-1; dum3int[2]=kr_step*kr;	// coarsen the access	
		rw_real_fullhyperslab_nd_h5(file_id,"outputs/output_field",&h5error,3,dims_input,dum3int,inputs.Efield.Field,"r");

		int Nz_CUPRAD, Nr_CUPRAD;
		// double *rgrid_CUPRAD, *zgrid_CUPRAD;
		rgrid_CUPRAD = readreal1Darray_fort(file_id, "outputs/rgrid", &h5error, &Nr_CUPRAD);
		zgrid_CUPRAD = readreal1Darray_fort(file_id, "outputs/zgrid", &h5error, &Nz_CUPRAD);


		h5error = H5Fclose(file_id);

		// convert units
		for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;}

		// do the calculation
		outputs = call1DTDSE(inputs); // THE TDSE

		// resize grids
		// double *rgrid_coarse, *zgrid_coarse;

		int Nr_coarse, Nz_coarse;
		coarsen_grid_real(rgrid_CUPRAD, Nr_CUPRAD, &rgrid_coarse, &Nr_coarse, kr_step, Nr_max);
		coarsen_grid_real(zgrid_CUPRAD, Nz_CUPRAD, &zgrid_coarse, &Nz_coarse, kz_step, Nz_max);


		// create local output file
		local_filename[0] = '\0'; dumchar1[0] = '\0'; sprintf(dumchar1, "%07d", myrank);
		strcat(local_filename,filename_stub); strcat(local_filename,dumchar1); strcat(local_filename,".h5");		
		file_id = H5Fcreate (local_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		prepare_local_output_fixed_print_grids_h5(file_id, "", &h5error, &inputs, &outputs, Ntot/nprocs + 1, dims);
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


		h5error = H5Fclose(file_id); // file
		outputs_destructor(&outputs); // clean ouputs

		free(rgrid_coarse); free(zgrid_coarse);
		free(rgrid_CUPRAD); free(zgrid_CUPRAD);
		
	}
	t_mpi[2] = MPI_Wtime(); 


	// process the MPI queue
	nxtval_strided(nprocs,&Nsim); Nsim_loc++;
	if ( ( comment_operation == 1 ) && ( myrank == 0 ) ){ printf("Proc %i c %i; time %f sec \n",myrank,Nsim, t_mpi[2]- t_mpi[1]); fflush(NULL);}
	//t_mpi[7] = MPI_Wtime();
	//printf("Proc %i, reached the point 2  : %f sec\n",myrank,t_mpi[7]-t_mpi[0]);
	while (Nsim < Ntot){ // run till queue is not treated
		t_mpi[3] = MPI_Wtime(); 
		kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dim_r; // compute offsets in each dimension

		// prepare the part in the arrray to r/w
		// dum3int[0]=-1; dum3int[1]=kr; dum3int[2]=kz; // set offset as inputs for hdf5-procedures
 		dum3int[0]=kz_step*kz; dum3int[1]=-1; dum3int[2]=kr_step*kr;	// coarsen the access

		dims_input[0] = dim_t;

		// read the HDF5 file
		file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		rw_real_fullhyperslab_nd_h5(file_id,"outputs/output_field",&h5error,3,dims_input,dum3int,inputs.Efield.Field,"r");
		h5error = H5Fclose(file_id);

		// convert units
		for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;}

		// do the calculation
		// t_mpi[3] = MPI_Wtime(); finish3_main = clock();
		outputs = call1DTDSE(inputs); // THE TDSE  
		// t_mpi[1] = MPI_Wtime(); finish1_main = clock();


    

		file_id = H5Fopen (local_filename, H5F_ACC_RDWR, H5P_DEFAULT); // open file
		print_local_output_fixed_h5(file_id,"", &h5error, &inputs, &outputs, Ntot/nprocs + 1, Nsim, Nsim_loc);

		h5error = H5Fclose(file_id); // file

		
		// outputs_destructor(outputs); // free memory
		outputs_destructor(&outputs);
		nxtval_strided(nprocs,&Nsim); Nsim_loc++;
		// printf("Proc %i c %i\n",myrank,Nsim); fflush(NULL);
		t_mpi[4] = MPI_Wtime();
		if ( ( comment_operation == 1 ) && ( myrank == 0 ) ){ printf("Proc %i c %i; time %f sec, from start %f sec \n",myrank,Nsim, t_mpi[4]- t_mpi[3], t_mpi[4]- t_mpi[1]); fflush(NULL);}
		
		t_mpi[5] = MPI_Wtime();
	}
	//h5error = H5Sclose(memspace_id);

	free(dims);
 

	t_mpi[6] = MPI_Wtime();
	if ( ( comment_operation == 1 )){ printf("Proc %i is going to finish., Total time %f sec.\n",myrank,t_mpi[6]-t_mpi[0]); fflush(NULL);}
	MPI_Finalize();
	return 0;	
}



/* to check if exists use printf("link exists 1: %i\n",H5Lexists(file_id, "outputs/lambda", H5P_DEFAULT)); */
