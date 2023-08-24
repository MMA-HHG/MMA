/*
This is for the locally the reading procedure.
*/

#include<time.h> 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
// include "/usr/local/hdf5/include/hdf5.h"
#include "hdf5.h"
#include "util.h"

// hdf5 operation:
herr_t  h5error;
hid_t file_id; // file pointer
hid_t filespace, dataspace_id, dataset_id; // dataspace pointers

int k1;




int main(int argc, char *argv[]) 
{

	int comment_operation = 1;



	// the file is opened for read only by all the processes independently, every process then has its own copy of variables.
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT);  

	// we first start with the t-grid
	hid_t dset_id = H5Dopen2 (file_id, "IRProp/tgrid", H5P_DEFAULT); // open dataset	     
	hid_t dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
	const int ndims = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions in the tgrid
	if ( comment_operation == 1 ){printf("dimensionality tgrid is: %i \n",ndims);}
	hsize_t dims[ndims]; // we need the size to allocate tgrid for us
	H5Sget_simple_extent_dims(dspace_id, dims, NULL); // get dimensions
	if ( comment_operation == 1 ){printf("Size 1 is: %i \nSize 2 is: %i \nGrid is from Fortran as a column, it gives the extra 1-dimension\n",dims[0],dims[1]);}
	hid_t datatype  = H5Dget_type(dset_id);     // we gat the type of data (SINGLE, DOUBLE, etc. from HDF5)
	double tgrid[dims[0]]; // allocate the grid
	h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tgrid); // read the grid
	if ( comment_operation == 1 ){printf("(t_init,t_end) = (%e,%e) \n",tgrid[0],tgrid[dims[0]-1]);}
	h5error = H5Dclose(dset_id);

	// we move to the Fields
	dset_id = H5Dopen2 (file_id, "IRProp/Fields_rzt", H5P_DEFAULT); // open dataset	     
	dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
	const int ndims2 = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions for the fields
	hsize_t dims2[ndims2]; // variable to access
	H5Sget_simple_extent_dims(dspace_id, dims2, NULL); // get dimensions
	if ( comment_operation == 1 ){printf("Fields dimensions (t,r,z) = (%i,%i,%i)\n",dims2[0],dims2[1],dims2[2]);}
	datatype  = H5Dget_type(dset_id);     // get datatype
	hsize_t dim_t = dims2[0], dim_r = dims2[1], dim_z = dims2[2]; // label the dims by physical axes

	// based on dimensions, we set a counter (queue length)
	int Ntot = dim_r*dim_z;
	// selections (hyperslabs) are needed
	hsize_t  offset[ndims2], stride[ndims2], count[ndims2], block[ndims2];

	double Fields[dims2[0]], SourceTerms[dims2[0]]; // Here we store the field and computed SOurce Term for every case	
	
	hsize_t field_dims[1]; // we need to specify the length of the array this way for HDF5
	field_dims[0] = dims2[0];

	hid_t memspace_id = H5Screate_simple(1,field_dims,NULL); // this memspace correspond to one Field/SourceTerm hyperslab, we will keep it accross the code


	// THE MAIN OPERATION LEADING TO OUTPUT STARTS HERE

	// prepare file
	printf("bprep: %i \n");
	
	file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // we use a different output file to testing, can be changed to have only one file
	dataspace_id = H5Screate_simple(ndims2, dims2, NULL); // create dataspace for outputs
	dataset_id = H5Dcreate2(file_id, "/SourceTerms", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // create dataset

	// close it
	h5error = H5Dclose(dataset_id); // dataset
	h5error = H5Sclose(dataspace_id); // dataspace
	h5error = H5Fclose(file_id); // file
        printf("aprep: %i \n");


        int Nsim, kr, kz; // counter of simulations, indices in the FIeld array

	Nsim = 0;




	// Hyperslabs wrongly working
/*	do { // run till queue is not treated*/
/*		kr = Nsim % dim_r; kz = Nsim - kr;  kz = kz / dim_r; // compute offsets in each dimension*/

/*		// prepare the part in the arrray to r/w*/
/*		offset[0] = 0; offset[1] = kr; offset[2] = kz; */
/*		stride[0] = 1; stride[1] = 1; stride[2] = 1;*/
/*		count[0] = dims2[0]; count[1] = 1; count[2] = 1; // takes all t*/
/*		block[0] = 1; block[1] = 1; block[2] = 1;*/

/*		// read the HDF5 file*/

/*		// MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL); // We now use different input and output file, input is for read-only, this mutex is here in the case we have only one file for I/O.*/
/*		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Will read from (kr,kz)=(%i,%i), job %i \n",kr,kz,Nsim);}*/

/*		file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // same as shown*/
/*		dset_id = H5Dopen2 (file_id, "IRProp/Fields_rzt", H5P_DEFAULT); */
/*		dspace_id = H5Dget_space (dset_id);*/

/*		h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block); // operation with only a part of the array = hyperslab*/
/*		h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, Fields); // read only the hyperslab*/

/*		h5error = H5Dclose(dset_id); // dataset*/
/*		h5error = H5Sclose(dspace_id); // dataspace*/
/*		h5error = H5Fclose(file_id); // file*/

/*		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Finished read of job %i \n", Nsim);}*/
/*		// MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);*/

/*		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Doing the job %i \n",Nsim);}*/

/*		// THE TASK IS DONE HERE, we can call 1D/3D TDSE, etc. here*/
/*		for (k1 = 0; k1 < dims2[0]; k1++){SourceTerms[k1]=2.0*Fields[k1];}; // just 2-multiplication*/
/*		*/
/*		// print the output in the file*/

/*		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Will write in the hyperslab (kr,kz)=(%i,%i), job %i \n",kr,kz,Nsim);}*/

/*		file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // open file*/
/*		dset_id = H5Dopen2 (file_id, "/SourceTerms", H5P_DEFAULT); // open dataset*/
/*		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Dset opened");}*/
/*		filespace = H5Dget_space (dset_id); // Get the dataspace ID   */
/*		h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block); // again the same hyperslab as for reading*/

/*		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){printf("Bwrite");}*/

/*		h5error = H5Dwrite(dset_id,datatype,memspace_id,filespace,H5P_DEFAULT,SourceTerms); // write the data*/

/*		// close*/
/*		h5error = H5Dclose(dset_id); // dataset*/
/*		h5error = H5Sclose(filespace); // dataspace*/
/*		h5error = H5Fclose(file_id); // file*/

/*		Nsim++;*/
/*	} while (Nsim < Ntot);*/
	h5error = H5Sclose(memspace_id);
	return 0;	
}



