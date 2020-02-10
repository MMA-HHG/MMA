/*
The HDF5 version will implement the following idea:

There is one parameter file shared by all simulations, and also only one I/O hdf5 file.
The parameters given to the code by slurm are two integers defining the indices in r and z.
There should be a logfile for noting succesful/failed simulations.
The *.output files should be stacked somewhere.

For reading, it should be easy. ** R/W may occur simultaneously in in the MPI loop. Separate I/O at the instant or ensure it will work (R/W from independent datasets may be fine???).
https://support.hdfgroup.org/HDF5/Tutor/selectsimple.html


After discussions, we try to test

a)mutex 
https://www.thegeekstuff.com/2012/05/c-mutex-examples/ https://www.geeksforgeeks.org/mutex-lock-for-linux-thread-synchronization/ .
https://computing.llnl.gov/tutorials/pthreads/#Mutexes
https://www.mcs.anl.gov/~robl/papers/ross_atomic-mpiio.pdf
https://stackoverflow.com/questions/37236499/mpi-ensure-an-exclusive-access-to-a-shared-memory-rma

b) temporary files
Each process writes in its own hdf5 file. There should be a way to do a "virtual" merging procedure: by using virtual datasets
https://portal.hdfgroup.org/display/HDF5/Introduction+to+the+Virtual+Dataset++-+VDS
For the instant, we may use a more direct method-store data in binary files etc. It may be easier for testing & debugging.



All the code will be encapsulated in an MPI-loop.

The plot of the code development:
1) we leave the original parametric file, the only difference will be omitting the filenames. Istead of this there gonna be two indices (r and z). Matrix size will be leaded from the hfd5 archive.
1.develop) first do only hdf5 stuff single run with fixed indices
2) we use strided MPI simulations.
2.develop) there should be an MPI-paradigm that allows to create a pool of jobs , we let a free process to take a job from the top of the buffer, implement it.
3) we use mutex to write into the hdf5 archive
3.test) we include a direct printing in separated files in the testing mode

*/
#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
// #include "hdf5.h"
#include "util.h"


//#pragma warning( disable : 4996 ) // warning for fopen in visual 2005



// vars
herr_t  h5error;
int k1;

int MPE_MUTEX_KEYVAL;


int main(int argc, char *argv[]) 
{	

    // MPI according to thread support for mutexes, see https://stackoverflow.com/questions/14836560/thread-safety-of-mpi-send-using-threads-created-with-stdasync and the Hristo Iliev's answer      

	// standard operation	 
	int myrank, nprocs;
	// MPI_Init(&argc, &argv); for non-threaded mpi

	// threaded operation
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided < MPI_THREAD_MULTIPLE)
	{
		printf("ERROR: The MPI library does not have full thread support\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


		printf("1 mutex  %d of %d\n", myrank, nprocs);
		printf("2 mutex  %d of %d\n", myrank, nprocs);
		printf("3 mutex  %d of %d\n", myrank, nprocs);
		printf("4 mutex  %d of %d\n", myrank, nprocs);
		printf("5 mutex  %d of %d\n", myrank, nprocs);
		printf("6 mutex  %d of %d\n", myrank, nprocs);
		printf("7 mutex  %d of %d\n", myrank, nprocs);


	// printf("I am node %d of %d\n", myrank, nprocs);

	// printf("program's running\n");

	// TESTING VARIABLES
	// int kz = 0;
	// int kr = 0;

		// OPEN HDF5 file	
        hid_t file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // open file // H5F_ACC_RDWR

		// find dimensions	
		hid_t dset_id = H5Dopen2 (file_id, "IRProp/tgrid", H5P_DEFAULT); // open dataset	     
        hid_t dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
		const int ndims = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions in the grid
		hsize_t dims[ndims]; // define dims variable
		// printf("ndim is: %i \n",ndims);
		H5Sget_simple_extent_dims(dspace_id, dims, NULL); // get dimensions


		// printf("Size 1 is: %i \n",dims[0]);
		// printf("Size 2 is: %i \n",dims[1]);
		// printf("Size is: %i \n",dims[ndims]);

		// read data
		hid_t datatype  = H5Dget_type(dset_id);     /* datatype handle */

		double tgrid[dims[0]][dims[1]];

		h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tgrid);

		// printf("test1: %lf \n",tgrid[5][0]);
		// printf("test2: %e \n",tgrid[5][0]); 

		/* Close the dataset. */
		h5error = H5Dclose(dset_id);
			

		// test hyperslab
		// find dimensions	
		dset_id = H5Dopen2 (file_id, "IRProp/Fields_rzt", H5P_DEFAULT); // open dataset	     
        dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
		const int ndims2 = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions in the grid
		hsize_t dims2[ndims2]; // define dims variable
		// printf("ndim is: %i \n",ndims2);
		H5Sget_simple_extent_dims(dspace_id, dims2, NULL); // get dimensions
		// printf("Size 1 is: %i \n",dims2[0]);	printf("Size 2 is: %i \n",dims2[1]); printf("Size 3 is: %i \n",dims2[2]);
		datatype  = H5Dget_type(dset_id);     /* datatype handle */

        int kz = 1;
		int kr = myrank;
		//make selectoin
		hsize_t  offset[ndims2];
        hsize_t  stride[ndims2];
        hsize_t  count[ndims2];
        hsize_t  block[ndims2];
		offset[0] = 0; offset[1] = kr; offset[2] = kz; 
		stride[0] = 1; stride[1] = 1; stride[2] = 1;
		count[0] = dims2[0]; count[1] = 1; count[2] = 1;
		block[0] = 1; block[1] = 1; block[2] = 1;

		double Fields[dims2[0]]; // offset adds these extra 1-dimensions... is there a way to remove them?

		hsize_t field_dims[1];
		field_dims[0] = dims2[0];
		hid_t memspace_id = H5Screate_simple(1,field_dims,NULL);

		h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block);
    	h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, Fields);

		// h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fields); // used for reading all

		// printf("test1: \n");
		// for(k1 = 0 ; k1 <= 5 ; k1++){
		// 	printf("%lf \n",Fields[k1]);
		// }
		

		// HERE WE TEST MUTEX FOR PRINTING, I added myrank as the offset for each worker

		// synchro according to https://www.mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf

		// https://cvw.cac.cornell.edu/MPIoneSided/lul

		int assert;
		// MPI_Win win;

		// MPI_Win_lock(MPI_LOCK_EXCLUSIVE, myrank, assert, win); // not sure with the type of the lock Indicates whether other processes may access the target window at the same time (if MPI_LOCK_SHARED) or not (MPI_LOCK_EXCLUSIVE)
	
		printf("1 I am node %d of %d\n", myrank, nprocs);
		printf("2 I am node %d of %d\n", myrank, nprocs);
		printf("3 I am node %d of %d\n", myrank, nprocs);
		printf("4 I am node %d of %d\n", myrank, nprocs);
		printf("5 I am node %d of %d\n", myrank, nprocs);
		printf("6 I am node %d of %d\n", myrank, nprocs);
		printf("7 I am node %d of %d\n", myrank, nprocs);
		printf("8 I am node %d of %d\n", myrank, nprocs);
		printf("9 I am node %d of %d\n", myrank, nprocs);
		printf("10 I am node %d of %d\n", myrank, nprocs);
		printf("11 I am node %d of %d\n", myrank, nprocs);
		printf("12 I am node %d of %d\n", myrank, nprocs);
		printf("13 I am node %d of %d\n", myrank, nprocs);
		printf("14 I am node %d of %d\n", myrank, nprocs);  


		// printf("process %d gives data %lf \n",myrank,Fields[2]);
		printf("process %d gives data \n",myrank);


		// MPI_Win_unlock(myrank, win);

		

		// printf("test1: %lf \n",Fields[1][1][1]);
		// printf("test2: %e \n",Fields[2][2][1]); 







        // h5sget_simple_extent_dims(dspace_id, dims, maxdims)  //Getting dims from dataspace

		// allocate fields

		// set hyperslab
		// dataset_id = H5Dopen2 (file_id, "IRprop/Fields_rzt", H5P_DEFAULT); // open dataset

		// load fields

		// close field
		
	
}