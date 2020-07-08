#include<time.h> 
#include<stdio.h>
#include<string.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "mpi.h"


void rw_real_fullhyperslab_nd_h5(hid_t file_id, char *dset_name, herr_t *h5error, int ndims, hsize_t *dimensions, int *selection, double *array1D, char *rw) // This function reads full line from an n-D array, the selected dimension is given by (-1), the rest of selection is the offset
{ 
  int k1;
  hid_t memspace_id;
  hsize_t field_dims[1];  
  hsize_t  offset[ndims], stride[ndims], count[ndims], block[ndims];

  for(k1 = 0; k1 < ndims; k1++){
    stride[k1] = 1; block[k1] = 1;
    if (selection[k1] < 0){
      offset[k1] = 0;
      count[k1] = dimensions[k1];
      field_dims[0] = dimensions[k1]; // a way to specify the length of the array for HDF5	
      memspace_id = H5Screate_simple(1,field_dims,NULL);
    }else{
      offset[k1] = selection[k1];
      count[k1] = 1;
    }
  }
  hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT);
  hid_t dspace_id = H5Dget_space (dset_id);
  hid_t datatype  = H5Dget_type(dset_id);
  *h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block); // operation with only a part of the array = hyperslab	
  if (strcmp(rw,"r")==0){
    *h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, array1D); // read only the hyperslab
  } else if (strcmp(rw,"w")==0){
    *h5error = H5Dwrite (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, array1D); // write the data
  } else {
    printf("wrongly sepcified r/w: nothing done\n"); 
  }
  
  *h5error = H5Dclose(dset_id); // dataset
  *h5error = H5Sclose(dspace_id); // dataspace
  *h5error = H5Sclose(memspace_id); // dataspace
}


void readreal(hid_t file_id, char *dset_name, herr_t *h5error, double *value)
{
  hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset
  hid_t datatype  = H5Dget_type(dset_id);
  *h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  *h5error = H5Dclose(dset_id);
}


void readint(hid_t file_id, char *dset_name, herr_t *h5error, int *value)
{
  hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset
  hid_t datatype  = H5Dget_type(dset_id);
  *h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  *h5error = H5Dclose(dset_id);
}


double * readreal1Darray_fort(hid_t file_id, char *dset_name, herr_t *h5error, int *N_points) // fort is for the extra diemnsion due to fortran
{
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset	     
	hid_t dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
	const int ndims = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions in the tgrid
	hsize_t dims[ndims]; // we need the size to allocate tgrid for us
	H5Sget_simple_extent_dims(dspace_id, dims, NULL); // get dimensions
	hid_t datatype  = H5Dget_type(dset_id);     // we get the type of data (SINGLE, DOUBLE, etc. from HDF5)
  double *array = malloc((int)dims[0]*sizeof(double)); 
	/*see https://stackoverflow.com/questions/10575544/difference-between-array-type-and-array-allocated-with-malloc
	      https://stackoverflow.com/questions/216259/is-there-a-max-array-length-limit-in-c/216731#216731  */
	*h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array); // read the grid
  *h5error = H5Sclose(dspace_id);
  *h5error = H5Dclose(dset_id);	

  *N_points = (int)dims[0];
  return array;
}


hsize_t * get_dimensions_h5(hid_t file_id, char *dset_name, herr_t *h5error, int *ndims, hid_t *datatype)
{
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset	     
	hid_t dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
	*ndims = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions for the fields
  hsize_t *dims = malloc((*ndims)*sizeof(hsize_t));
	H5Sget_simple_extent_dims(dspace_id, dims, NULL); // get dimensions
 
	*datatype  = H5Dget_type(dset_id); // get datatype

  *h5error = H5Sclose(dspace_id);
	*h5error = H5Dclose(dset_id);
 
  return dims;
}

//int linkexists(hid_t file_id, char *link_name, herr_t *h5error)
//{
//  hid_t dset_id = H5Lexists (file_id, link_name, H5P_DEFAULT); // open dataset
//  hid_t datatype  = H5Dget_type(dset_id);
//  *h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
//  *h5error = H5Dclose(dset_id);
//}

void addone(int *val ){*val=*val+1;} // to test pointers
