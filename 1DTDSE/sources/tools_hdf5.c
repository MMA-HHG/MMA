/**
 * @file tools_hdf5.c
 * @brief Contains functions for read/write, printing and other 
 * operations with HDF5 files.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "tools_hdf5.h"
#include "constants.h"
#include "structures.h"

int one = 1;

/**
 * @brief Adds units to HDF5 varibles.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset to write.
 * @param h5error Status.
 * @param units_value Name of the unit.
 */
void add_units_1D_h5(hid_t file_id, char *dset_name, herr_t *h5error, char *units_value) {
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT);
	hsize_t dumh51D[1] = {1};
	hid_t aspace_id = H5Screate_simple(1, dumh51D, dumh51D);
	hid_t atype_id = H5Tcopy(H5T_C_S1);
	*h5error = H5Tset_size(atype_id, (size_t) 10);
	hid_t attr_id = H5Acreate2(dset_id, "units", atype_id, aspace_id, H5P_DEFAULT, H5P_DEFAULT);
	*h5error = H5Awrite(attr_id, atype_id, units_value);
	*h5error = H5Aclose(attr_id);
	*h5error = H5Tclose(atype_id);
	*h5error = H5Sclose(aspace_id);
	*h5error = H5Dclose(dset_id);
}

/**
 * @brief Shortcut to HDF5 datatypes.
 * 
 * @param foo Specifies the datatype, "d" = double, "s" = float, "i" = integer.
 * @return hid_t Datatype
 */
hid_t dtype_h5(char *foo)
{
	if (strcmp(foo, "d") == 0){
		return H5T_NATIVE_DOUBLE;
	} else if (strcmp(foo, "s") == 0){
		return H5T_NATIVE_FLOAT;
	} else if (strcmp(foo, "i") == 0){
		return H5T_NATIVE_INT;
	} else {
		fprintf(stderr, "Wrongly specified r/w: nothing done. \n"); 
		return H5T_ORDER_ERROR; // should choose better error, it's just a random valid error I've found.
	}
}

/**
 * @brief Writes data from the HDF5 input into the ```inputs_def``` structure.
 * 
 * @param file_id HDF5 file.
 * @param inpath Path to the input group.
 * @param h5error Status.
 * @param in Input structure.
 */
void ReadInputs(hid_t file_id, char *inpath, char *inpath_glob, herr_t *h5error, inputs_def *in)
{
	// Dummy string with path to the input value
	char path[50];
	char *dumstring;
	int gas_preset_exists=1;
	// Energy of the initial state
	path[0] = '\0';	strcat(strcat(path,inpath),"Eguess");
	readreal(file_id, path, h5error,&(*in).Eguess); 
	// Number of points of the initial spatial grid
	path[0] = '\0';	strcat(strcat(path,inpath),"Nx_max");
	readint(file_id, path, h5error,&(*in).num_r); 
	// resolution for the grid
	path[0] = '\0';	strcat(strcat(path,inpath),"dx");
	readreal(file_id, path, h5error,&(*in).dx); 
	// ??? Interpolation variables
	path[0] = '\0';	strcat(strcat(path,inpath),"InterpByDTorNT");
	readint(file_id, path, h5error,&(*in).InterpByDTorNT);
	// Resolution in time
	path[0] = '\0';	strcat(strcat(path,inpath),"dt");
	readreal(file_id, path, h5error,&(*in).dt); 
	// Number of points of the spatial grid for the expansion
	path[0] = '\0';	strcat(strcat(path,inpath),"Ntinterp");
	readint(file_id, path, h5error,&(*in).Ntinterp); 
	// the limit of the integral for the ionisation (works fine with the lenth gauge and strong fields)
	path[0] = '\0';	strcat(strcat(path,inpath),"x_int");
	readreal(file_id, path, h5error,&(*in).x_int);  

	// Target parameter
	path[0] = '\0';	strcat(strcat(path,inpath_glob),"gas_preset");

	printf("Path to preset %s.\n", path);
	printf("Path to preset exists %d.\n", H5Lexists(file_id, path, H5P_DEFAULT));

	if (H5Lexists(file_id, path, H5P_DEFAULT)>0){
		readstring(file_id, path, h5error, &dumstring);
		Soft_Coulomb_parameters(dumstring, &(*in).trg.a);
		free(dumstring);
		gas_preset_exists=0;
	}

	path[0] = '\0';	strcat(strcat(path,inpath),"trg_a");

	printf("Path to trg_a %s.\n", path);
	printf("Path to trg_a exists %d.\n", H5Lexists(file_id, path, H5P_DEFAULT));

	if (H5Lexists(file_id, path, H5P_DEFAULT)>0){
		readreal(file_id, path, h5error,&(*in).trg.a);
		if (gas_preset_exists==0){
			printf("WARNING: Using soft-Coulomb from CTDSE inputs, ignoring gas_preset.\n");
		}

	}else if (gas_preset_exists==1){
        printf("The soft-Coulomb parameter in 1D-TDSE not specified.");
        exit(EXIT_FAILURE);
	}
	 

	// *NOTE* CV criterion is added as an input 
	// Load CV criterion
	path[0] = '\0';	strcat(strcat(path,inpath),"CV_criterion_of_GS");
	readreal(file_id, path, h5error,&(*in).CV);
	// Gauge selection
	path[0] = '\0';	strcat(strcat(path,inpath),"gauge_type");
	readint(file_id, path, h5error,&(*in).gauge); // analyse 2nd part of the dipole

	strcpy((*in).precision,"d");
}

/**
 * @brief Reads electric field and time grid
 * 
 * @param file_id HDF5 file.
 * @param inpath Path to input.
 * @param h5error Status.
 * @param in Input structure.
 */
void Read_1_field_and_grid(hid_t file_id, char *inpath, herr_t *h5error,  inputs_def *in)
{
	char path[50];
	printf("read 1d\n"); fflush(NULL);
	path[0] = '\0';	strcat(strcat(path,inpath),"tgrid");
	(*in).Efield.tgrid =  readreal1Darray_fort(file_id, path, h5error, &(*in).Efield.Nt);
	path[0] = '\0';	strcat(strcat(path,inpath),"Field");
	(*in).Efield.Field =  readreal1Darray_fort(file_id, path, h5error, &(*in).Efield.Nt); 
}

/**
 * @brief Creates dataset and writes an nd array into HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset.
 * @param h5error Status.
 * @param ndims Number of dimensions of the dataset.
 * @param dimensions Dimensions of the dataset.
 * @param array Array to write.
 * @param datatype Datatype.
 */
void print_nd_array_h5(hid_t file_id, char *dset_name, herr_t *h5error, int ndims, hsize_t *dimensions, void * array, hid_t datatype)
{
	hid_t dspace_id = H5Screate_simple(ndims, dimensions, NULL);
	hid_t dset_id = H5Dcreate2(file_id, dset_name, datatype, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	*h5error = H5Dwrite(dset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
	*h5error = H5Sclose(dspace_id);
	*h5error = H5Dclose(dset_id);	
}

/**
 * @brief Creates dataset for an nd array in HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset.
 * @param h5error Status.
 * @param ndims Number of dimensions of the dataset.
 * @param dimensions Dimensions of the dataset.
 * @param datatype Datatype.
 */
void create_nd_array_h5(hid_t file_id, char *dset_name, herr_t *h5error, 
						int ndims, hsize_t *dimensions, hid_t datatype)
{
	hid_t dspace_id = H5Screate_simple(ndims, dimensions, NULL);
	hid_t dset_id = H5Dcreate2(file_id, dset_name, datatype, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	*h5error = H5Sclose(dspace_id);
	*h5error = H5Dclose(dset_id);	
}

/**
 * @brief Writes output of the TDSE computation into HDF5 file
 * 
 * @param file_id HDF5 file.
 * @param inpath Path to outputs in the HDF5 file.
 * @param h5error Status.
 * @param in Input structure.
 * @param out Output structure.
 */
void PrintOutputs(hid_t file_id, char *inpath, herr_t *h5error,  inputs_def *in,  
				  outputs_def *out)
{
	// Dimensions never exceeds 2 in this case, can be longer
	hsize_t output_dims[2]; 
	char path[50];

  	// time domain
	output_dims[0] = (*out).Nt; output_dims[1] = 0;
	if ( (*in).Print.Efield == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"Efield");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).Efield, H5T_NATIVE_DOUBLE);
	}

  	if ( (*in).Print.sourceterm == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"SourceTerm");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).sourceterm, H5T_NATIVE_DOUBLE);
	}

	if ( (*in).Print.PopTot == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"PopTot");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).PopTot, H5T_NATIVE_DOUBLE);
	}

	// the grid
	if ( (*in).Print.Efield == 1 || (*in).Print.sourceterm == 1 || (*in).Print.PopTot == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"tgrid");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).tgrid, H5T_NATIVE_DOUBLE);
	}

	// omega domain - complex
	output_dims[0] = (*out).Nomega; output_dims[1] = 2;
	if ( (*in).Print.FEfield == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FEfield");
		print_nd_array_h5(file_id, path, h5error, 2, output_dims, (*out).FEfield, H5T_NATIVE_DOUBLE);
	}

	if ( (*in).Print.Fsourceterm == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTerm");
		print_nd_array_h5(file_id, path, h5error, 2, output_dims, (*out).Fsourceterm, H5T_NATIVE_DOUBLE);
	}

	// omega domain - real
	output_dims[0] = (*out).Nomega; output_dims[1] = 0;

	if ( (*in).Print.FEfieldM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FEfieldM2");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).FEfieldM2, H5T_NATIVE_DOUBLE);
	}

	if ( (*in).Print.FsourceTermM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTermM2");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).FsourcetermM2, H5T_NATIVE_DOUBLE);
	}

	// the grid
	if ( (*in).Print.FEfield == 1 || (*in).Print.Fsourceterm == 1 || (*in).Print.FEfieldM2 == 1 || (*in).Print.FsourceTermM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"omegagrid");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).omegagrid, H5T_NATIVE_DOUBLE);
	}

	// various scalars
	output_dims[0] = 1; output_dims[1] = 0;
	path[0] = '\0'; strcat(strcat(path,inpath),"Energy_of_the_ground_state");
	print_nd_array_h5(file_id, path, h5error, 1, output_dims, &(*in).Einit, H5T_NATIVE_DOUBLE);
}

/**
 * @brief Reads/writes full line from an nd array, the selected dimension is given by (-1), 
 * the rest of selection is the offset.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset of the nd array.
 * @param h5error Status.
 * @param nhyperslab_dimensions Dimension of the nd array.
 * @param hyperslab_dimensions_int Dimensions of the nd array.
 * @param offset_int Offset in the chosen dimension.
 * @param count_int Total number of offsets.
 * @param array Array for reading/writing.
 * @param rw Read or write operation, {'r', 'w'}.
 */
void rw_hyperslab_nd_h5(hid_t file_id, char *dset_name, herr_t *h5error, 
						int nhyperslab_dimensions, int *hyperslab_dimensions_int, 
						int *offset_int, int *count_int, void *array, char *rw) 
{ 
	int ndims, k1;
	hsize_t *dum_dims;
	hid_t datatype;
	dum_dims = get_dimensions_h5(file_id, dset_name, h5error, &ndims, &datatype);
	hsize_t hyperslab_dimensions[nhyperslab_dimensions];
	for(k1 = 0; k1 < nhyperslab_dimensions; k1++){hyperslab_dimensions[k1] = hyperslab_dimensions_int[k1];}
	hid_t memspace_id = H5Screate_simple(nhyperslab_dimensions,hyperslab_dimensions,NULL);
	hsize_t  offset[ndims], count[ndims];
	for(k1 = 0; k1 < ndims; k1++) {
		offset[k1] = offset_int[k1]; 
		count[k1] = count_int[k1];
	}

	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT);
	hid_t dspace_id = H5Dget_space (dset_id);
	// operation with only a part of the array = hyperslab	
	*h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, NULL, count, NULL); 
	// read only the hyperslab
	if (strcmp(rw,"r")==0){
		*h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, array); 
	} else if (strcmp(rw,"w")==0){
		// write the data
		*h5error = H5Dwrite (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, array); 
	} else {
		printf("wrongly sepcified r/w: nothing done\n"); 
	}
	
	*h5error = H5Dclose(dset_id); // dataset
	*h5error = H5Sclose(dspace_id); // dataspace
	*h5error = H5Sclose(memspace_id); // dataspace
	free(dum_dims);
}

/**
 * @brief Reads/writes full line from an nd array, the selected dimension is given 
 * by (-1), the rest of selection is the offset.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset of the nd array.
 * @param h5error Status.
 * @param ndims Dimension of the nd array.
 * @param dimensions Dimensions of the nd array.
 * @param selection Indeces for array selection.
 * @param array1D Array for read/write.
 * @param rw Read or write operation, {'r', 'w'}.
 */
void rw_real_fullhyperslab_nd_h5(hid_t file_id, char *dset_name, herr_t *h5error, int ndims, hsize_t *dimensions, int *selection, double *array1D, char *rw) 
{ 
	int k1;
	hid_t memspace_id;
	hsize_t field_dims[1];  
	hsize_t  offset[ndims], stride[ndims], count[ndims], block[ndims];

	for(k1 = 0; k1 < ndims; k1++){
		stride[k1] = 1; 
		block[k1] = 1;
		if (selection[k1] < 0){
			offset[k1] = 0;
			count[k1] = dimensions[k1];
			// a way to specify the length of the array for HDF5	
			field_dims[0] = dimensions[k1]; 
		} else {
			offset[k1] = selection[k1];
			count[k1] = 1;
		}
	}
	memspace_id = H5Screate_simple(1,field_dims,NULL);
	
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT);
	hid_t dspace_id = H5Dget_space (dset_id);
	// operation with only a part of the array = hyperslab	
	*h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, stride, count, block); 
	if (strcmp(rw,"r")==0){
		// read only the hyperslab
		*h5error = H5Dread (dset_id, H5T_NATIVE_DOUBLE, memspace_id, dspace_id, H5P_DEFAULT, array1D); 
	} else if (strcmp(rw,"w")==0){
		// write the data
		*h5error = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memspace_id, dspace_id, H5P_DEFAULT, array1D); 
	} else {
		printf("wrongly sepcified r/w: nothing done\n"); 
	}
	
	*h5error = H5Dclose(dset_id); // dataset
	*h5error = H5Sclose(dspace_id); // dataspace
	*h5error = H5Sclose(memspace_id); // dataspace
}

/**
 * @brief Reads real variable from HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset.
 * @param h5error Status.
 * @param value Variable for storing the value.
 */
void readreal(hid_t file_id, char *dset_name, herr_t *h5error, double *value)
{
	// open dataset
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); 
	*h5error = H5Dread(dset_id,  H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
	*h5error = H5Dclose(dset_id);
}

/**
 * @brief Reads integer variable from HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset.
 * @param h5error Status.
 * @param value Variable for storing the value.
 */
void readint(hid_t file_id, char *dset_name, herr_t *h5error, int *value)
{
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); 
	*h5error = H5Dread(dset_id,  H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
	*h5error = H5Dclose(dset_id);
}

/**
 * @brief Reads string from HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Name of the dataset.
 * @param h5error Status.
 * @param value Variable for storing the value.
 */
void readstring(hid_t file_id, char *dset_name, herr_t *h5error, char **value)
{
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); 
	hid_t dtype_id = H5Dget_type(dset_id);
    // hid_t dspace_id = H5Dget_space(dataset_id);

	size_t str_size = H5Tget_size(dtype_id); // Get the size of the datatype (string length)
	*value = (char *)malloc((str_size + 1) * sizeof(char)); // +1 for null terminator

	H5Dread(dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, *value);


	(*value)[str_size] = '\0'; // make it a valid C string

	
	*h5error = H5Tclose(dtype_id);
	*h5error = H5Dclose(dset_id);
}

/**
 * @brief Returns real 1D array from the dataset created using Fortran code. 
 * 
 * @details Fort is for the extra diemnsion due to fortran.
 * Further notes: https://stackoverflow.com/questions/10575544/difference-between-array-type-and-array-allocated-with-malloc
	      		  https://stackoverflow.com/questions/216259/is-there-a-max-array-length-limit-in-c/216731#216731
 * 
 * @param file_id HDF5 file.
 * @param dset_name Dataset name.
 * @param h5error Status.
 * @param N_points Size of the dataset.
 * @return double*
 */
double * readreal1Darray_fort(hid_t file_id, char *dset_name, herr_t *h5error, int *N_points) 
{
	// open dataset	     
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); 
	// Get the dataspace ID     
	hid_t dspace_id = H5Dget_space (dset_id); 
	// number of dimensions in the tgrid
	const int ndims = H5Sget_simple_extent_ndims(dspace_id); 
	// we need the size to allocate tgrid for us
	hsize_t dims[ndims]; 
	// get dimensions
	H5Sget_simple_extent_dims(dspace_id, dims, NULL); 
  	double *array = malloc((int)dims[0]*sizeof(double)); 
	// read the grid
	*h5error = H5Dread(dset_id,  H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array); 
	*h5error = H5Sclose(dspace_id);
	*h5error = H5Dclose(dset_id);	

	*N_points = (int)dims[0];
	return array;
}

/**
 * @brief Returns dimensions of HDF5 dataset.
 * 
 * @param file_id HDF5 file.
 * @param dset_name Dataset name.
 * @param h5error Status.
 * @param ndims Dimension of the nd array.
 * @param datatype Datatype of the array.
 * @return hsize_t* 
 */
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

/**
 * @brief Allocates output arrays for the output of a single TDSE computation 
 * into a temporary HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param inpath Path to outputs in the HDF5 file.
 * @param h5error Status.
 * @param in Input structure.
 * @param out Output structure.
 * @param nsimulations Total simulations.
 * @param dims Field dimensions.
 */
void prepare_local_output_fixed_print_grids_h5(hid_t file_id, char *inpath, 
	herr_t *h5error,  inputs_def *in,  outputs_def *out, int nsimulations, 
	hsize_t *dims)
{
	hsize_t output_dims[3];
	char path[50];
	int k1;

	// keys for post-processing
	output_dims[0] = nsimulations;
	path[0] = '\0';	strcat(strcat(path,inpath),"keys");
	int initial_keys[nsimulations]; 
	for(k1 = 0 ; k1 < nsimulations; k1++) {
		initial_keys[k1] = -1;
	}  
	print_nd_array_h5(file_id, path, h5error, 1, output_dims, initial_keys, H5T_NATIVE_INT);

	// time domain
	output_dims[0] = (*out).Nt; output_dims[1] = nsimulations;
	if ( (*in).Print.Efield == 1 ) 
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"Efield");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	if ( (*in).Print.sourceterm == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"SourceTerm");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	if ( (*in).Print.PopTot == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"PopTot");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	if ( (*in).Print.PopInt == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"PopInt");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	if ( (*in).Print.expval_x == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"expval_x");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	// the grid
	output_dims[0] = (*out).Nt; output_dims[1] = 0;
	if ((*in).Print.Efield == 1 || (*in).Print.sourceterm == 1 || 
		(*in).Print.PopTot == 1  || (*in).Print.PopInt == 1 || 
		(*in).Print.expval_x == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"tgrid");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).tgrid, H5T_NATIVE_DOUBLE);
	}

	// omega domain - complex
	output_dims[0] = (*out).Nomega; output_dims[1] = 2; 
	output_dims[2] = nsimulations;
	if ( (*in).Print.FEfield == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FEfield");
		create_nd_array_h5(file_id, path, h5error, 3, output_dims, dtype_h5((*in).precision));
	}

	if ( (*in).Print.Fsourceterm == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTerm");
		create_nd_array_h5(file_id, path, h5error, 3, output_dims, dtype_h5((*in).precision));
	}

	// omega domain - real
	output_dims[0] = (*out).Nomega; output_dims[1] = nsimulations;

	if ( (*in).Print.FEfieldM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FEfieldM2");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	if ( (*in).Print.FsourceTermM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTermM2");
		create_nd_array_h5(file_id, path, h5error, 2, output_dims, dtype_h5((*in).precision));
	}

	// the grid
	output_dims[0] = (*out).Nomega; output_dims[1] = 0;
	if ((*in).Print.FEfield == 1 || (*in).Print.Fsourceterm == 1 || 
		(*in).Print.FEfieldM2 == 1 || (*in).Print.FsourceTermM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"omegagrid");
		print_nd_array_h5(file_id, path, h5error, 1, output_dims, (*out).omegagrid, H5T_NATIVE_DOUBLE);
	}

	// various scalars
	output_dims[0] = 1; output_dims[1] = 0;
	path[0] = '\0'; strcat(strcat(path,inpath),"Energy_of_the_ground_state");
	print_nd_array_h5(file_id, path, h5error, 1, output_dims, &(*in).Einit, H5T_NATIVE_DOUBLE);

	path[0] = '\0'; strcat(strcat(path,inpath),"Nr_orig"); int foo = dims[2];
	print_nd_array_h5(file_id, path, h5error, 1, output_dims, &foo, H5T_NATIVE_INT);

	path[0] = '\0'; strcat(strcat(path,inpath),"Nz_orig"); foo = dims[0];
	print_nd_array_h5(file_id, path, h5error, 1, output_dims, &foo, H5T_NATIVE_INT);

	path[0] = '\0'; strcat(strcat(path,inpath),"Nt_orig"); foo = dims[1];
	print_nd_array_h5(file_id, path, h5error, 1, output_dims, &foo, H5T_NATIVE_INT);

	path[0] = '\0'; strcat(strcat(path,inpath),"number_of_local_simulations");
	create_nd_array_h5(file_id, path, h5error, 1, output_dims, H5T_NATIVE_INT);
}


/**
 * @brief Writes output of a single TDSE computation into a temporary HDF5 file.
 * 
 * @param file_id HDF5 file.
 * @param inpath Path to outputs in the HDF5 file.
 * @param h5error Status.
 * @param in Input structure.
 * @param out Output structure.
 * @param nsimulations Total simulations.
 * @param Nsim Number of parallel jobs.
 * @param Nsim_loc Local job.
 */
void print_local_output_fixed_h5(hid_t file_id, char *inpath, herr_t *h5error,  
	inputs_def *in,  outputs_def *out, int nsimulations, int Nsim, int Nsim_loc)
{
	hsize_t output_dims[3]; 
  	int offsets[3];
	char path[50];

	// keys for post-processing
	output_dims[0] = nsimulations;
	path[0] = '\0';	strcat(strcat(path,inpath),"keys");
	rw_hyperslab_nd_h5(file_id, path, h5error, one, &one, &Nsim_loc, &one, &Nsim, "w");

	// time domain
	output_dims[0] = (*out).Nt; output_dims[1] = nsimulations;
	offsets[0] = -1; offsets[1] = Nsim_loc;
	if ( (*in).Print.Efield == 1 ) 
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"Efield");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).Efield, "w");
	}

	if ( (*in).Print.sourceterm == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"SourceTerm");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).sourceterm, "w");
	}

	if ( (*in).Print.PopTot == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"PopTot");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).PopTot, "w");
	}

	if ( (*in).Print.PopTot == 1 ) 
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"PopInt");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).PopInt, "w");
	}

	if ( (*in).Print.PopTot == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"expval_x");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).expval, "w");
	}

	// omega domain - complex
	int hcount[3] = {(*out).Nomega,2,1};
	int hoffset[3] = {0,0,Nsim_loc};
	int dimsloc[2] = {(*out).Nomega,2};
	output_dims[0] = (*out).Nomega; 
	output_dims[1] = 2; 
	output_dims[2] = nsimulations;
	
	if ( (*in).Print.FEfield == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FEfield");	
		rw_hyperslab_nd_h5(file_id, path, h5error, 2, dimsloc, hoffset, hcount, (*out).FEfield, "w");
	}

	if ( (*in).Print.Fsourceterm == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTerm");
		rw_hyperslab_nd_h5(file_id, path, h5error, 2, dimsloc, hoffset, hcount, (*out).Fsourceterm, "w");
	}

	// omega domain - real
	output_dims[0] = (*out).Nomega; 
	output_dims[1] = nsimulations;
	offsets[0] = -1; 
	offsets[1] = Nsim_loc;

	if ( (*in).Print.FEfieldM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FEfieldM2");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).FEfieldM2, "w");
	}

	if ( (*in).Print.FsourceTermM2 == 1 )
	{
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTermM2");
		rw_real_fullhyperslab_nd_h5(file_id, path, h5error, 2, output_dims, offsets, (*out).FsourcetermM2, "w");
	}

	// various scalars
	output_dims[0] = 1; output_dims[1] = 0;
	
	path[0] = '\0'; 
	strcat(strcat(path,inpath),"number_of_local_simulations");
	int foo = Nsim_loc + 1;
	hid_t dset_id = H5Dopen2 (file_id, path, H5P_DEFAULT);
	*h5error = H5Dwrite (dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &foo);
	*h5error = H5Dclose(dset_id); // dataset
}


/**
 * @brief Returns printing structure for HDF5 writing and sets all prints according
 * to the input HDF5 file.
 * 
 * @return output_print_def 
 */
output_print_def Set_prints_from_HDF5(hid_t file_id, char *inpath, herr_t *h5error)
{
	output_print_def res;
	char path[50];
	int dum_int;

	res = Initialise_Printing_struct();
	path[0] = '\0';	strcat(strcat(path,inpath),"print_Efield");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.Efield = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_Source_Term");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.sourceterm = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_F_Source_Term");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.Fsourceterm = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_F_Efield_M2");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.FEfieldM2 = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_F_Efield");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.FEfield = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_F_Source_Term_M2");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.FsourceTermM2 = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_GS_population");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.PopTot = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_integrated_population");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.PopInt = 1;}

	path[0] = '\0';	strcat(strcat(path,inpath),"print_x_expectation_value");
	readint(file_id, path, h5error,&dum_int);
	if(dum_int==1){res.expval_x = 1;}

	res.tgrid = 1; // not memory consuming
	res.omegagrid = 1; // not memory consuming

	return res;
}
