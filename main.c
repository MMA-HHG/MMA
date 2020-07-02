/*
The plot of the code:

1) See the comment at the end of the file for deatils about the MPI-scheduler.

2) The main program (MP) reads all the parameters from an input HDF5-archive (each process independently, read-only).

3) MP decides based on parameters the type of input field. Fist implementation: stick to numerical fields from CUPRAD stored in the archive.

4) MPI-scheduler executes a simulation in every point. THe data are directly written into the output using the mutex.

5) The code finishes.

DEVELOPMENT: wrap the call of singleTDSE, where propagation is called to test, erase this extra step after

-----------------------------------------------------------------------------------------------------------------------
Extensions/features already presented in 1DTDSE and needed to implement in a test mode. We wil add them as optional outputs.

The features already presented in TDSE:
	1) Print wavefunction
	2) Print Gabor transformation
	3) Photoelectron spectrum (Fabrice)

2) is computationally demanding; 1) 2) are both data-storage demandig. It should be used then only in few prescribed points.


There is possibility of various inputs:
	1) Numerical/analytic field
	2) Computation on velocity/length gauge

We have already an analytic model of a beam (Python/MATLAB), we will then rewrite it and construct the parameters on-the-fly.
The versatility in numeric-or-analytic field length-or-velocity gauge is ideal for testing of numerical vector potential that we can use after in SFA.


------------------------------------------------------------------------------------------------------------------------
Development notes:

1) We can use checks whether parameters exist in the input HDF5-archive. If not, we can create them with default values.
Implementation: since this is I/O operation with one file, we need r/w. Maybe read paramc only by proc 0 and the broadcast structure (see the MPI book for transfering structs).

For reading, it should be easy. ** R/W may occur simultaneously in in the MPI loop. Separate I/O at the instant or ensure it will work (R/W from independent datasets may be fine???).
https://support.hdfgroup.org/HDF5/Tutor/selectsimple.html


2) we get rid of mutexes and use rather parallel acces to files all the time.
2.develop) it seems that many-readers many-writers would be possible by HDF5 parallel since we will not modify the file much. However, we may also try stick with independent files and eventually 
https://stackoverflow.com/questions/49851046/merge-all-h5-files-using-h5py
https://portal.hdfgroup.org/display/HDF5/Collective+Calling+Requirements+in+Parallel+HDF5+Applications

*/

#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "util.h"

// hdf5 operation:
herr_t  h5error;
hid_t file_id; // file pointer
hid_t filespace, dataspace_id, dataset_id, dset_id, dspace_id; // dataspace pointers

struct inputs_def inputs;
struct outputs_def outputs;

int k1, k2, k3;


clock_t start_main, finish2_main, finish1_main, finish3_main, finish4_main;


int main(int argc, char *argv[]) 
{

	int comment_operation = 1;

	// Initialise MPI
	int myrank, nprocs;
	int MPE_MC_KEYVAL; // this is used to address the mutex and counter
	MPI_Win mc_win; // this is the shared window, it is used both  for mutices and counter
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	if (comment_operation == 1 ){printf("Proc %i started the program\n",myrank);}

	////////////////////////
	// PREPARATION PAHASE //
	////////////////////////	

	// READ DATA
	/* to check if exists use printf("link exists 1: %i\n",H5Lexists(file_id, "IRProp/lambda", H5P_DEFAULT)); */
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.

	readreal(file_id, "TDSE_inputs/Eguess"					,&h5error,&inputs.Eguess); // Energy of the initial state
	readint(file_id, "TDSE_inputs/N_r_grid"					,&h5error,&inputs.num_r); // Number of points of the initial spatial grid 16000
	readint(file_id, "TDSE_inputs/N_r_grid_exp"				,&h5error,&inputs.num_exp); // Number of points of the spatial grid for the expansion
	readreal(file_id, "TDSE_inputs/dx"						,&h5error,&inputs.dx); // resolution for the grid
	readint(file_id, "TDSE_inputs/InterpByDTorNT"			,&h5error,&inputs.InterpByDTorNT); 
	readreal(file_id, "TDSE_inputs/dt"						,&h5error,&inputs.dt); // resolution in time
	readint(file_id, "TDSE_inputs/Ntinterp"					,&h5error,&inputs.Ntinterp); // Number of points of the spatial grid for the expansion
	readreal(file_id, "TDSE_inputs/textend"					,&h5error,&inputs.textend); // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
	readint(file_id, "TDSE_inputs/analy_writewft"			,&h5error,&inputs.analy.writewft); // writewavefunction (1-writting every tprint)
	readreal(file_id, "TDSE_inputs/analy_tprint"			,&h5error,&inputs.analy.tprint); // time spacing for writing the wavefunction	
	readreal(file_id, "TDSE_inputs/x_int"					,&h5error,&inputs.x_int); // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
	readint(file_id, "TDSE_inputs/PrintGaborAndSpectrum"	,&h5error,&inputs.PrintGaborAndSpectrum); // print Gabor and partial spectra (1-yes)
	readreal(file_id, "TDSE_inputs/a_Gabor"					,&h5error,&inputs.a_Gabor); // the parameter of the gabor window [a.u.]
	readreal(file_id, "TDSE_inputs/omegaMaxGabor"			,&h5error,&inputs.omegaMaxGabor); // maximal frequency in Gabor [a.u.]
	readreal(file_id, "TDSE_inputs/dtGabor"					,&h5error,&inputs.dtGabor); // spacing in Gabor
	readreal(file_id, "TDSE_inputs/tmin1window"				,&h5error,&inputs.tmin1window); // analyse 1st part of the dipole
	readreal(file_id, "TDSE_inputs/tmax1window"				,&h5error,&inputs.tmax1window); // analyse 1st part of the dipole
	readreal(file_id, "TDSE_inputs/tmin2window"				,&h5error,&inputs.tmin2window); // analyse 2nd part of the dipole
	readreal(file_id, "TDSE_inputs/tmax2window"				,&h5error,&inputs.tmax2window); // analyse 2nd part of the dipole
	readint(file_id, "TDSE_inputs/PrintOutputMethod"		,&h5error,&inputs.PrintOutputMethod); // (0 - only text, 1 - only binaries, 2 - both)
//	readint(file_id, "TDSE_inputs/IonisationFilterForTheSourceTerm"	,&h5error,&inputs.IonisationFilterForTheSourceTerm); // filter source term by high-ionisation components (1-yes)
//	readreal(file_id, "TDSE_inputs/IonFilterThreshold"		,&h5error,&inputs.IonFilterThreshold); // threshold for the ionisation [-]

	if (comment_operation == 1 ){printf("Proc %i uses dx = %e \n",myrank,inputs.dx);}
	

	// load the tgrid
	inputs.Efield.tgrid =  readreal1Darray_fort(file_id, "IRProp/tgrid",&h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs
	// convert to atomic units
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]*1e15*41.34144728;}

	// dimension of the 3D array containing all the inputs
	hsize_t * dims; int ndims; hid_t datatype;
	dims = get_dimensions_h5(file_id, "IRProp/Fields_rzt", &h5error, &ndims, &datatype);
	hsize_t dim_t = dims[0], dim_r = dims[1], dim_z = dims[2]; // label the dims by physical axes
	if ( ( comment_operation == 1 ) && ( myrank == 0 ) ){printf("Fields dimensions (t,r,z) = (%i,%i,%i)\n",dims[0],dims[1],dims[2]);}


	// PREPARE DATA

	// We prepare the storage for data and set up the run	
	int Ntot = dim_r*dim_z; // counter (queue length)	
	hsize_t  offset[ndims], stride[ndims], count[ndims], block[ndims]; // selections (hyperslabs) are needed	
	hsize_t field_dims[1]; field_dims[0] = dims[0]; // a way to specify the length of the array for HDF5	
	hid_t memspace_id = H5Screate_simple(1,field_dims,NULL); // this memspace correspond to one Field/SourceTerm hyperslab, we will keep it accross the code
	double Fields[dims[0]], SourceTerms[dims[0]]; // Here we store the field and computed Source Term for every case




	// Prepare the ground state
	if ( myrank == 0 )
	{

	// Initialise vectors and Matrix 
	// Initialise_GS(inputs.num_r);
	
	//normalise(psi0,inputs.num_r); // Initialise psi0 for Einitialise
	//normalise(psiexc,inputs.num_r);

	double CV = 1E-20; // CV criteria

	/* This number has to be small enough to assure a good convregence of the wavefunction
	if it is not the case, then the saclar product of the the ground state and the excited states 
	is not quite 0 and those excited appears in the energy analysis of the gorund states, so the propagation !!
	CV = 1E-25 has been choosen to have a scalar product of 10^-31 with the third excited state for num_r = 5000 and dx=0.1
	*/
	
	printf("Calculation of the energy of the ground sate ; Eguess : %f\n",inputs.Eguess);
	int size = 2*(inputs.num_r+1);
	double *psi0, *off_diagonal, *diagonal, *x, *psiexc;
	double Einit = 0.0, Einit2 = 0.0;
	psi0 = calloc(size,sizeof(double));
	psiexc = calloc(size,sizeof(double));
	for(k1=0;k1<=inputs.num_r;k1++){psi0[2*k1] = 1.0; psi0[2*k1+1] = 0.; psiexc[2*k1] = 1; psiexc[2*k1+1] = 0.;}
	printf("binit\n");
	Initialise_grid_and_D2(inputs.dx, inputs.num_r, &x, &diagonal, &off_diagonal);
	printf("bEinit\n");
	Einit = Einitialise(inputs.trg,psi0,off_diagonal,diagonal,off_diagonal,x,inputs.Eguess,CV,inputs.num_r);
	//for(i=0;i<=inputs.num_r;i++) {fprintf(eingenvectorf,"%f\t%e\t%e\n",x[i],psi0[2*i],psi0[2*i+1]); fprintf(pot,"%f\t%e\n",x[i],potential(x[i],trg));}

	printf("Initial energy is : %1.12f\n",Einit);
	printf("first excited energy is : %1.12f\n",Einit2);
	exit(0);
	}

	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////

	// create counter and mutex in one pointer
	MPE_MC_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 2, &mc_win); // first is counter, second mutex


	if ( myrank == 0 )  // first process is preparing the file and the rest may do their own work (the file is locked in the case they want to write); it creates the resulting dataset
	{
	MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL);

	file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // we use a different output file to testing, can be changed to have only one file
	dataspace_id = H5Screate_simple(ndims, dims, NULL); // create dataspace for outputs
	dataset_id = H5Dcreate2(file_id, "/SourceTerms", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // create dataset

	// close it
	h5error = H5Dclose(dataset_id); // dataset
	h5error = H5Sclose(dataspace_id); // dataspace
	h5error = H5Fclose(file_id); // file

	MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);
 
	}
	// an empty dataset is prepared to be filled with the data
 
	start_main = clock(); // the clock
	finish4_main = clock();
		

	// we now process the MPI queue
	int Nsim, kr, kz; // counter of simulations, indices in the FIeld array
	MPE_Counter_nxtval(mc_win, 0, &Nsim, MPE_MC_KEYVAL); // get my first task (every worker calls)

	do { // run till queue is not treated
		kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dim_r; // compute offsets in each dimension

		// prepare the part in the arrray to r/w
		offset[0] = 0; offset[1] = kr; offset[2] = kz; 
		stride[0] = 1; stride[1] = 1; stride[2] = 1;
		count[0] = dims[0]; count[1] = 1; count[2] = 1; // takes all t
		block[0] = 1; block[1] = 1; block[2] = 1;

		// read the HDF5 file

		// MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL); // We now use different input and output file, input is for read-only, this mutex is here in the case we have only one file for I/O.
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i will read from (kr,kz)=(%i,%i), job %i \n",myrank,kr,kz,Nsim);}

		file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // same as shown
		dset_id = H5Dopen2 (file_id, "IRProp/Fields_rzt", H5P_DEFAULT); 
		dspace_id = H5Dget_space (dset_id);

		h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block); // operation with only a part of the array = hyperslab
		h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, Fields); // read only the hyperslab

		h5error = H5Dclose(dset_id); // dataset
		h5error = H5Sclose(dspace_id); // dataspace
		h5error = H5Fclose(file_id); // file

		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i finished read of job %i \n",myrank, Nsim);}
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i, job %i, field[0] = %e \n",myrank, Nsim,Fields[0]);}
		// MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);

		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i doing the job %i \n",myrank,Nsim);}

		// THE TASK IS DONE HERE, we can call 1D/3D TDSE, etc. here
		for (k1 = 0; k1 < inputs.Efield.Nt; k1++){SourceTerms[k1]=2.0*Fields[k1];}; // just 2-multiplication
   
		inputs.Efield.Field = Fields;
   
    		finish3_main = clock();
		outputs = call1DTDSE(inputs); // THE TDSE
   
   
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i finished TDSE job %i \n",myrank,Nsim);}
		printf("address2 %p \n",outputs.Efield);
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("%e \n",outputs.Efield[0]);}
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i, job %i some outputs are: %e, %e, %e \n",myrank,Nsim, outputs.tgrid[0], outputs.Efield[0], outputs.sourceterm[0]);}
		// for (k1 = 0; k1 < inputs.Efield.Nt; k1++){SourceTerms[k1]=outputs.sourceterm[k1];}; // assign results
    		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Proc %i finished TDSE job %i \n",myrank,Nsim);}
		
		// print the output in the file
		finish1_main = clock();
    
		MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL); // mutex is acquired

		finish2_main = clock();
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){
			printf("Proc %i will write in the hyperslab (kr,kz)=(%i,%i), job %i \n",myrank,kr,kz,Nsim);
			printf("Proc %i, returned mutex last time  : %f sec\n",myrank,(double)(finish4_main - start_main) / CLOCKS_PER_SEC);
			printf("Proc %i, before job started        : %f sec\n",myrank,(double)(finish3_main - start_main) / CLOCKS_PER_SEC);
			printf("Proc %i, clock the umnutexed value : %f sec\n",myrank,(double)(finish1_main - start_main) / CLOCKS_PER_SEC);
			printf("Proc %i, clock in the mutex block  : %f sec\n",myrank,(double)(finish2_main - start_main) / CLOCKS_PER_SEC);
			printf("first element to write: %e \n",SourceTerms[0]);
			fflush(NULL); // force write
    	}
    

		file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // open file
		dset_id = H5Dopen2 (file_id, "/SourceTerms", H5P_DEFAULT); // open dataset
		filespace = H5Dget_space (dset_id); // Get the dataspace ID   
		h5error = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block); // again the same hyperslab as for reading

		h5error = H5Dwrite(dset_id,datatype,memspace_id,filespace,H5P_DEFAULT,SourceTerms); // write the data

		// close
		h5error = H5Dclose(dset_id); // dataset
		h5error = H5Sclose(filespace); // dataspace
		h5error = H5Fclose(file_id); // file

		MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);
    	finish4_main = clock();
		MPE_Counter_nxtval(mc_win, 0, &Nsim, MPE_MC_KEYVAL); // get my next task
	} while (Nsim < Ntot);
	h5error = H5Sclose(memspace_id);
 
//  if ( myrank == 0 )  // time
//	{
//  finish2 = clock();
//  printf("\nFirst processor measured time: %f sec\n\n",(double)(finish2 - start2) / CLOCKS_PER_SEC);
//  }
  // MPI_Win_free(&mc_win);
	MPI_Finalize();
	return 0;	
}


/*
Here we show a simple scheduler of tasks from a queue. In this case, each task is specified by kr and kz that correpond to two indices in a table (:,kr,kz).
The elementary program process data from that cut. The processors takes tasks until the queue is empty. There is no synchoronisation at this point, so tasks
can be of various lenghts, etc. THe code than could be very easily adapted to any task.

The program oprates with one input file and one output file. It could be desirable to use only one file, this would require to either use SWMR technique or
mutex the accesses ( https://portal.hdfgroup.org/pages/viewpage.action?pageId=48812567 ). 

The limitation is that the mutex is used for writing, the writing should take then small amount of time compared to the atomic task.

Concrete decription of this tutorial: it takes the fields from results.h5[/IRProp/Fields_rzt] , multiplies the array by 2 using the forementioned procedure
and prints the result in results.h5[/SourceTerms]. Next, it also loads results.h5[/IRProp/tgrid] before the calculation.

The code may explain its work by setting comment_operation = 1; and muted by comment_operation = 0; Only first 20 tasks are shown.

Possible extensions are discussed at the end of the file.
*/
