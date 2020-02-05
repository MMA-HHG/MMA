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
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "util.h"

//#pragma warning( disable : 4996 ) // warning for fopen in visual 2005



int main(void)
{	

	printf("program's running\n");

	// TESTING VARIABLES
	int kz = 0;
	int kr = 0;

		// OPEN HDF5 file	
        hid_t file_id = H5Fopen ("results.h5", H5F_ACC_RDWR, H5P_DEFAULT); // open file

		// find dimensions	
		hid_t dset_id = H5Dopen2 (file_id, "IRProp/tgrid", H5P_DEFAULT); // open dataset	     
        hid_t dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     

		const int ndims = H5Sget_simple_extent_ndims(dspace_id);
		hsize_t dims[ndims];
		printf("ndim is: %i \n",ndims);
		H5Sget_simple_extent_dims(dspace_id, dims, NULL);

		printf("Size 1 is: %i \n",dims[0]);
		printf("Size 2 is: %i \n",dims[1]);
		// printf("Size is: %i \n",dims[ndims]);

        // h5sget_simple_extent_dims(dspace_id, dims, maxdims)  //Getting dims from dataspace

		// allocate fields

		// set hyperslab
		// dataset_id = H5Dopen2 (file_id, "IRprop/Fields_rzt", H5P_DEFAULT); // open dataset

		// load fields

		// close field
		
	
}