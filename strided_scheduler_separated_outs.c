#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
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
	hsize_t * dims; int ndims; hid_t datatype;
	char dumchar1[50], dumchar2[50];
	// Processing the queue
	int Nsim, Nsim_loc = -1, kr, kz; // counter of simulations, indices in the Field array

	int comment_operation = 1;
	double t_mpi[10]; 

	hsize_t output_dims[4];

	////////////////////////
	// PREPARATION PAHASE //
	////////////////////////

	// Initialise MPI
	int myrank, nprocs;
	int MPE_MC_KEYVAL, MPE_C_KEYVAL, MPE_M_KEYVAL; // this is used to address the mutex and counter
	MPI_Win mc_win, c_win, m_win; // this is the shared window, it is used both  for mutices and counter
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	nxtval_init(-nprocs+myrank,&Nsim);

	Init_constants();
	inputs.Print = Initialise_Printing_struct();
	inputs.Print = Set_all_prints();

	if (comment_operation == 1 ){printf("Proc %i started the program\n",myrank);}

	t_mpi[0] = MPI_Wtime();	start_main = clock(); // the clock	

	// create parameters & load initial data
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.
	ReadInputs(file_id, "TDSE_inputs/", &h5error, &inputs);
	dims = get_dimensions_h5(file_id, "IRProp/Fields_rzt", &h5error, &ndims, &datatype);
	hsize_t dim_t = dims[0], dim_r = dims[1], dim_z = dims[2]; // label the dims by physical axes	

	// create space for the fields & load the tgrid
	inputs.Efield.Field = malloc(((int)dims[0])*sizeof(double));
	inputs.Efield.tgrid =  readreal1Darray_fort(file_id, "IRProp/tgrid",&h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs
	// convert units
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; /*inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;*/} // convert to atomic units (fs->a.u.), (GV/m->a.u.)


	if (comment_operation == 1 ){printf("Proc %i uses dx = %e \n",myrank,inputs.dx);}
	if ( ( comment_operation == 1 ) && ( myrank == 0 ) ){printf("Fields dimensions (t,r,z) = (%i,%i,%i)\n",dims[0],dims[1],dims[2]);}


	// Prepare the ground state (it's the state of the atom before the interaction)
	Initialise_grid_and_ground_state(&inputs);

	int Ntot = dim_r*dim_z; // counter (queue length)	



	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////


	// first simulation prepares the outputfile
	nxtval_strided(nprocs,&Nsim); Nsim_loc++;
	if (Nsim < Ntot){

		// find proper simulation & load the field
		file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
		kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dim_r; // compute offsets in each dimension
		dum3int[0]=-1; dum3int[1]=kr; dum3int[2]=kz;		
		rw_real_fullhyperslab_nd_h5(file_id,"IRProp/Fields_rzt",&h5error,3,dims,dum3int,inputs.Efield.Field,"r");
		h5error = H5Fclose(file_id);

		// convert units
		for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.Field[k1] = inputs.Efield.Field[k1]*1e-15;}

		// do the calculation
		outputs = call1DTDSE(inputs); // THE TDSE

		// prepare the output file
		dims[0] = outputs.Nt; // length defined by outputs



		local_filename[0] = '\0'; dumchar1[0] = '\0'; sprintf(dumchar1, "%07d", myrank);
		strcat(local_filename,filename_stub); strcat(local_filename,dumchar1); strcat(local_filename,".h5");
		
		
		//strcat(strcat(path,inpath),"Eguess");
		file_id = H5Fcreate (local_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		prepare_local_output_fixed_print_grids_h5(file_id, "", &h5error, &inputs, &outputs, Ntot/nprocs + 1);
		print_local_output_fixed_h5(file_id,"", &h5error, &inputs, &outputs, Ntot/nprocs + 1, Nsim, Nsim_loc);

		// output_dims[0] = Ntot/nprocs + 1;
		// //create_nd_array_h5(file_id, "/keys", &h5error, 1, output_dims, H5T_NATIVE_INT);
		// rw_hyperslab_nd_h5(file_id, "/keys", &h5error, one, &one, &Nsim_loc, &one, &Nsim, "w");

		// output_dims[0] = outputs.Nt; output_dims[1] = Ntot/nprocs + 1;
		// //create_nd_array_h5(file_id, "/Efield", &h5error, 2, output_dims, H5T_NATIVE_DOUBLE);
		// dum3int[0] = Nsim; dum3int[0] = -1;
		// //rw_real_fullhyperslab_nd_h5(file_id, "/Efield", &h5error, 2, output_dims, dum3int, outputs.Efield, "w");
		// // rw_real_fullhyperslab_nd_h5(file_id,"/SourceTerms",&h5error,3,dims,dum3int,outputs.Efield,"w");

		// output_dims[0] = outputs.Nomega; output_dims[1] = 2; output_dims[2] = Ntot/nprocs + 1;
		// create_nd_array_h5(file_id, "/FEfield", &h5error, 3, output_dims, H5T_NATIVE_DOUBLE);

		// int hcount[3] = {outputs.Nomega,2,1};
		// int hoffset[3] = {0,0,Nsim_loc};
		// int dimsloc[2] = {outputs.Nomega,2};
		// rw_hyperslab_nd_h5(file_id, "/FEfield", &h5error, 2, dimsloc, hoffset, hcount, outputs.FEfield_data, "w");

		h5error = H5Fclose(file_id); // file
	}

	MPI_Finalize();
	return 0;	




	//t_mpi[6] = MPI_Wtime();
	//printf("Proc %i, reached the point 1  : %f sec\n",myrank,t_mpi[6]-t_mpi[0]);
	
	// first process prepare file based on the first simulation
	// first process release mutex
 
	//t_mpi[4] = MPI_Wtime();	finish4_main = clock();
		

	// we now process the MPI queue
	nxtval_strided(nprocs,&Nsim); Nsim_loc++;
	printf("Proc %i c %i\n",myrank,Nsim); fflush(NULL);

	//t_mpi[7] = MPI_Wtime();
	//printf("Proc %i, reached the point 2  : %f sec\n",myrank,t_mpi[7]-t_mpi[0]);

	//t_mpi[5] = MPI_Wtime(); 
	while (Nsim < Ntot){ // run till queue is not treated
		kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dim_r; // compute offsets in each dimension

		// prepare the part in the arrray to r/w
		dum3int[0]=-1; dum3int[1]=kr; dum3int[2]=kz; // set offset as inputs for hdf5-procedures
		dims[0] = dim_t;

		// read the HDF5 file

		// MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL); // We now use different input and output file, input is for read-only, this mutex is here in the case we have only one file for I/O.
		file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // same as shown
		rw_real_fullhyperslab_nd_h5(file_id,"IRProp/Fields_rzt",&h5error,3,dims,dum3int,inputs.Efield.Field,"r");



		h5error = H5Fclose(file_id); // file
		// MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);

    		t_mpi[3] = MPI_Wtime(); finish3_main = clock();
		outputs = call1DTDSE(inputs); // THE TDSE  
		t_mpi[1] = MPI_Wtime(); finish1_main = clock();


    

		file_id = H5Fopen (dumchar2, H5F_ACC_RDWR, H5P_DEFAULT); // open file

		//dims[0] = outputs.Nt;
		//rw_real_fullhyperslab_nd_h5(file_id,"/SourceTerms",&h5error,3,dims,dum3int,outputs.Efield,"w");

		rw_hyperslab_nd_h5(file_id, "/keys", &h5error, one, &one, &Nsim_loc, &one, &Nsim, "w");

		output_dims[0] = outputs.Nt; output_dims[1] = Ntot/nprocs + 1;
		// create_nd_array_h5(file_id, "/Efield", &h5error, 2, output_dims, H5T_NATIVE_DOUBLE);
		dum3int[0] = Nsim; dum3int[0] = -1;
		rw_real_fullhyperslab_nd_h5(file_id, "/Efield", &h5error, 2, output_dims, dum3int, outputs.Efield, "w");


		h5error = H5Fclose(file_id); // file

		
		// outputs_destructor(outputs); // free memory
		nxtval_strided(nprocs,&Nsim); Nsim_loc++;
		printf("Proc %i c %i\n",myrank,Nsim); fflush(NULL);
		t_mpi[5] = MPI_Wtime();
	}
	//h5error = H5Sclose(memspace_id);

	free(dims);
 
//  if ( myrank == 0 )  // time
//	{
//  finish2 = clock();
//  printf("\nFirst processor measured time: %f sec\n\n",(double)(finish2 - start2) / CLOCKS_PER_SEC);
//  }
  // MPI_Win_free(&mc_win);
	MPI_Finalize();
	return 0;	
}



/* to check if exists use printf("link exists 1: %i\n",H5Lexists(file_id, "IRProp/lambda", H5P_DEFAULT)); */
