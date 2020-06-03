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
The main idea comes from the MPI3 book.
(
https://www.thegeekstuff.com/2012/05/c-mutex-examples/ https://www.geeksforgeeks.org/mutex-lock-for-linux-thread-synchronization/ .
https://computing.llnl.gov/tutorials/pthreads/#Mutexes
https://www.mcs.anl.gov/~robl/papers/ross_atomic-mpiio.pdf
https://stackoverflow.com/questions/37236499/mpi-ensure-an-exclusive-access-to-a-shared-memory-rma
)

b) temporary files
Each process writes in its own hdf5 file. There should be a way to do a "virtual" merging procedure: by using virtual datasets
https://portal.hdfgroup.org/display/HDF5/Introduction+to+the+Virtual+Dataset++-+VDS
For the instant, we may use a more direct method-store data in binary files etc. It may be easier for testing & debugging.



All the code will be encapsulated in an MPI-loop.

The plot of the code development:
1) we leave the original parametric file, the only difference will be omitting the filenames. Istead of this there gonna be two indices (r and z). Matrix size will be leaded from the hfd5 archive.
1.develop) first do only hdf5 stuff single run with fixed indices
2) Pool of processes should be easy with NXTVAL from MPI3. We implement it directly. The RMA window for mutex and queue is shared.

note: mutexes will be there for writing, the counter will be used for assigning simulations to workers at the moment they finish their work. I think it ensures maximal fair-share of the load.

3) we use mutex to write into the hdf5 archive
3.test) we include a direct printing in separated files in the testing mode

4) we get rid of mutexes and use rather parallel acces to files all the time.
4.develop) it seems that many-readers many-writers would be possible by HDF5 parallel since we will not modify the file much. However, we may also try stick with independent files and eventually 
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



// vars
herr_t  h5error;
hid_t file_id; // file pointer
hid_t filespace
int k1;

int MPE_MC_KEYVAL;


int main(int argc, char *argv[]) 
{	

    // MPI according to thread support for mutexes, see https://stackoverflow.com/questions/14836560/thread-safety-of-mpi-send-using-threads-created-with-stdasync and the Hristo Iliev's answer      

	// standard operation	 
	int myrank, nprocs;
	// MPI_Init(&argc, &argv); for non-threaded mpi

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);


		// OPEN HDF5 file	- this will be only file for reading: we open it independently by all processes
		// https://portal.hdfgroup.org/pages/viewpage.action?pageId=48812567
        file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // open file // H5F_ACC_RDWR

		// we first start with the t-grid
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
			


		// now we proceed with the fields
		// test hyperslab
		// find dimensions	
		dset_id = H5Dopen2 (file_id, "IRProp/Fields_rzt", H5P_DEFAULT); // open dataset	     
        dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
		const int ndims2 = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions in the grid
		hsize_t dims2[ndims2]; // define dims variable
		// printf("ndim is: %i \n",ndims2);
		H5Sget_simple_extent_dims(dspace_id, dims2, NULL); // get dimensions
		printf("Size 1 is: %i \n",dims2[0]);	printf("Size 2 is: %i \n",dims2[1]); printf("Size 3 is: %i \n",dims2[2]);
		datatype  = H5Dget_type(dset_id);     /* datatype handle */

		hsize_t dim_t = dims2[0], dim_r = dims2[1], dim_z = dims2[2]; // label the dims by physical axes

		// based on dimendions, we establish a counter
		int Ntot = dim_r*dim_z;

		hsize_t  offset[ndims2];
        hsize_t  stride[ndims2];
        hsize_t  count[ndims2];
        hsize_t  block[ndims2];
		double Fields[dims2[0]]; // offset adds these extra 1-dimensions... is there a way to remove them?
		double SourceTerms[dims2[0]]; 
		hsize_t field_dims[1];
		field_dims[0] = dims2[0];
		hid_t memspace_id = H5Screate_simple(1,field_dims,NULL); // reusing memspace: try it

		// create common counter and mutex
		MPE_MC_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 2, &mc_win); // first is counter, second mutex


		if ( myrank == 0 )  // first process is preparing the file and the rest may do their own work; it creates the resulting dataset
		{
		MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL);
		file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // we use a differenti file to testing
		dataspace_id = H5Screate_simple(ndims2, dims2, NULL); // create dataspace

		dataset_id = H5Dcreate2(file_id, "/SourceTerms", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // create dataset

		// close it
		h5error = H5Dclose(dataset_id); // dataset
		h5error = H5Sclose(dataspace_id); // dataspace
		h5error = H5Fclose(file_id); // file

		MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);
		}
		

		// we now process the MPI queue
		MPE_Counter_nxtval(mc_win, 0, &Nsim, MPE_MC_KEYVAL); // get my first simulation

		do { // run till queue is not treated
			kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dimr; // compute offsets in each dimension

			// prepare the part in the file to r/w
			offset[0] = 0; offset[1] = kr; offset[2] = kz; 
			stride[0] = 1; stride[1] = 1; stride[2] = 1;
			count[0] = dims2[0]; count[1] = 1; count[2] = 1; // takes all t
			block[0] = 1; block[1] = 1; block[2] = 1;

			// read the HDF5 file
			/* !!!!!!!!!!!! We use only reading from separate datasets, need to test, since datasets are exclusive, we probably don't need a lock. 
			Lege artis would be use SWMR approach. */
			file_id = H5Fopen ("results2.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

			dset_id = H5Dopen2 (file_id, "IRProp/Fields_rzt", H5P_DEFAULT); // open dataset	     
			dspace_id = H5Dget_space (dset_id); // Get the dataspace ID 

			h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block);
			h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, Fields);

			h5error = H5Dclose(dataset_id); // dataset
			h5error = H5Sclose(dspace_id); // dataspace
			h5error = H5Fclose(file_id); // file

			// do the job here
			for (k1 = 0; k1 < dims2[0]; k1++){SourceTerms[k1]=2.0*Fields[k1];}; // just 2-multiplication
			
			// print the output in the file
			MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL); // mutex is acquired
			// write in the file
			// open dataset
			file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // open file // H5F_ACC_RDWR
			dset_id = H5Dopen2 (file_id, "/SourceTerms", H5P_DEFAULT); // open dataset	

			filespace = H5Dget_space (dset_id); // Get the dataspace ID   
			h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block);

			h5error = H5Dwrite(dset_id, datatype, memspace_id, filespace,H5P_DEFAULT, SourceTerms);

			// close
			h5error = H5Dclose(dset_id); // dataset
			h5error = H5Sclose(memspace_id); // dataspace
			h5error = H5Sclose(filespace); // dataspace
			h5error = H5Fclose(file_id); // file

			MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);
			MPE_Counter_nxtval(mc_win, 0, &Nsim, MPE_MC_KEYVAL);
		} while (Nsim < Ntot);











		// h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Fields); // used for reading all

		// printf("test1: \n");
		// for(k1 = 0 ; k1 <= 5 ; k1++){
		// 	printf("%lf \n",Fields[k1]);
		// }
		

		// HERE WE TEST MUTEX FOR PRINTING, I added myrank as the offset for each worker

		// synchro according to https://www.mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf

		// https://cvw.cac.cornell.edu/MPIoneSided/lul

		// int assert;
		// MPI_Win win;

		// MPI_Win_lock(MPI_LOCK_EXCLUSIVE, myrank, assert, win); // not sure with the type of the lock Indicates whether other processes may access the target window at the same time (if MPI_LOCK_SHARED) or not (MPI_LOCK_EXCLUSIVE)
	

		printf("process %d gives data %lf \n",myrank,Fields[2]);
		// printf("process %d gives data \n",myrank);


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