#include<time.h> 
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>

#include "hdf5.h"
#include "numerical_constants.h"
#include "util.h"
#include "util_hdf5.h"

hid_t dtype_h5(char *foo)
{
  if (strcmp(foo,"d")==0){
    return H5T_NATIVE_DOUBLE;
  } else if (strcmp(foo,"s")==0){
    return H5T_NATIVE_FLOAT;
  } else if (strcmp(foo,"i")==0){
    return H5T_NATIVE_INT;
  } else {
    fprintf(stderr,"wrongly sepcified r/w: nothing done\n"); 
    return H5T_ORDER_ERROR; // should choose better error, it's just a random valid error I've found.
  }
}

void ReadInputs(hid_t file_id, char *inpath, herr_t *h5error, struct inputs_def *in)
{
	char path[50];
	printf("t1 \n"); fflush(NULL);

  path[0] = '\0';	strcat(strcat(path,inpath),"Eguess");
  readreal(file_id, path, h5error,&(*in).Eguess); // Energy of the initial state

  path[0] = '\0';	strcat(strcat(path,inpath),"N_r_grid");
  readint(file_id, path, h5error,&(*in).num_r); // Number of points of the initial spatial grid 16000

  path[0] = '\0';	strcat(strcat(path,inpath),"N_r_grid_exp");
  readint(file_id, path, h5error,&(*in).num_exp); // Number of points of the spatial grid for the expansion

  path[0] = '\0';	strcat(strcat(path,inpath),"dx");
  readreal(file_id, path, h5error,&(*in).dx); // resolution for the grid

  path[0] = '\0';	strcat(strcat(path,inpath),"InterpByDTorNT");
  readint(file_id, path, h5error,&(*in).InterpByDTorNT);

  path[0] = '\0';	strcat(strcat(path,inpath),"dt");
  readreal(file_id, path, h5error,&(*in).dt); // resolution in time

  path[0] = '\0';	strcat(strcat(path,inpath),"Ntinterp");
  readint(file_id, path, h5error,&(*in).Ntinterp); // Number of points of the spatial grid for the expansion

  path[0] = '\0';	strcat(strcat(path,inpath),"textend");
  readreal(file_id, path, h5error,&(*in).textend); // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700

  path[0] = '\0';	strcat(strcat(path,inpath),"analy_writewft");
  readint(file_id, path, h5error,&(*in).analy.writewft); // writewavefunction (1-writting every tprint)

  path[0] = '\0';	strcat(strcat(path,inpath),"analy_tprint");
  readreal(file_id, path, h5error,&(*in).analy.tprint); // time spacing for writing the wavefunction	

  path[0] = '\0';	strcat(strcat(path,inpath),"x_int");
  readreal(file_id, path, h5error,&(*in).x_int); // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields

  path[0] = '\0';	strcat(strcat(path,inpath),"PrintGaborAndSpectrum");
  readint(file_id, path, h5error,&(*in).PrintGaborAndSpectrum); // print Gabor and partial spectra (1-yes)

  path[0] = '\0';	strcat(strcat(path,inpath),"a_Gabor");
  readreal(file_id, path, h5error,&(*in).a_Gabor); // the parameter of the gabor window [a.u.]

  path[0] = '\0';	strcat(strcat(path,inpath),"omegaMaxGabor");
  readreal(file_id, path, h5error,&(*in).omegaMaxGabor); // maximal frequency in Gabor [a.u.]

  path[0] = '\0';	strcat(strcat(path,inpath),"dtGabor");
  readreal(file_id, path, h5error,&(*in).dtGabor); // spacing in Gabor

  path[0] = '\0';	strcat(strcat(path,inpath),"tmin1window");
  readreal(file_id, path, h5error,&(*in).tmin1window); // analyse 1st part of the dipole

  path[0] = '\0';	strcat(strcat(path,inpath),"tmax1window");
  readreal(file_id, path, h5error,&(*in).tmax1window); // analyse 1st part of the dipole

  path[0] = '\0';	strcat(strcat(path,inpath),"tmin2window");
  readreal(file_id, path, h5error,&(*in).tmin2window); // analyse 2nd part of the dipole

  path[0] = '\0';	strcat(strcat(path,inpath),"tmax2window");
  readreal(file_id, path, h5error,&(*in).tmax2window); // analyse 2nd part of the dipole

  path[0] = '\0';	strcat(strcat(path,inpath),"PrintOutputMethod");
  readint(file_id, path, h5error,&(*in).PrintOutputMethod); // (0 - only text, 1 - only binaries, 2 - both)

  path[0] = '\0';	strcat(strcat(path,inpath),"trg_a");
  readreal(file_id, path, h5error,&(*in).trg.a); // analyse 2nd part of the dipole

  // CV criterion will be added as an input
  (*in).CV = 1E-20; 

	//(*in).Efield.tgrid =  readreal1Darray_fort(file_id, "IRField/tgrid",h5error,&(*in).Efield.Nt); // tgrid is not changed when program runs
	//(*in).Efield.Field =  readreal1Darray_fort(file_id, "IRField/Field",h5error,&(*in).Efield.Nt); // tgrid is not changed when program runs  


// these two aren't in this version waiting to reintroduce
//	readint(file_id, "TDSE_inputs/IonisationFilterForTheSourceTerm"	,&h5error,&inputs.IonisationFilterForTheSourceTerm); // filter source term by high-ionisation components (1-yes)
//	readreal(file_id, "TDSE_inputs/IonFilterThreshold"		,&h5error,&inputs.IonFilterThreshold); // threshold for the ionisation [-]


// these are old ones
	// readreal(file_id, "TDSE_inputs/Eguess"					,&h5error,&inputs.Eguess); // Energy of the initial state
	// readint(file_id, "TDSE_inputs/N_r_grid"					,&h5error,&inputs.num_r); // Number of points of the initial spatial grid 16000
	// readint(file_id, "TDSE_inputs/N_r_grid_exp"				,&h5error,&inputs.num_exp); // Number of points of the spatial grid for the expansion
	// readreal(file_id, "TDSE_inputs/dx"						,&h5error,&inputs.dx); // resolution for the grid
	// readint(file_id, "TDSE_inputs/InterpByDTorNT"			,&h5error,&inputs.InterpByDTorNT); 
	// readreal(file_id, "TDSE_inputs/dt"						,&h5error,&inputs.dt); // resolution in time
	// readint(file_id, "TDSE_inputs/Ntinterp"					,&h5error,&inputs.Ntinterp); // Number of points of the spatial grid for the expansion
	// readreal(file_id, "TDSE_inputs/textend"					,&h5error,&inputs.textend); // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
	// readint(file_id, "TDSE_inputs/analy_writewft"			,&h5error,&inputs.analy.writewft); // writewavefunction (1-writting every tprint)
	// readreal(file_id, "TDSE_inputs/analy_tprint"			,&h5error,&inputs.analy.tprint); // time spacing for writing the wavefunction	
	// readreal(file_id, "TDSE_inputs/x_int"					,&h5error,&inputs.x_int); // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
	// readint(file_id, "TDSE_inputs/PrintGaborAndSpectrum"	,&h5error,&inputs.PrintGaborAndSpectrum); // print Gabor and partial spectra (1-yes)
	// readreal(file_id, "TDSE_inputs/a_Gabor"					,&h5error,&inputs.a_Gabor); // the parameter of the gabor window [a.u.]
	// readreal(file_id, "TDSE_inputs/omegaMaxGabor"			,&h5error,&inputs.omegaMaxGabor); // maximal frequency in Gabor [a.u.]
	// readreal(file_id, "TDSE_inputs/dtGabor"					,&h5error,&inputs.dtGabor); // spacing in Gabor
	// readreal(file_id, "TDSE_inputs/tmin1window"				,&h5error,&inputs.tmin1window); // analyse 1st part of the dipole
	// readreal(file_id, "TDSE_inputs/tmax1window"				,&h5error,&inputs.tmax1window); // analyse 1st part of the dipole
	// readreal(file_id, "TDSE_inputs/tmin2window"				,&h5error,&inputs.tmin2window); // analyse 2nd part of the dipole
	// readreal(file_id, "TDSE_inputs/tmax2window"				,&h5error,&inputs.tmax2window); // analyse 2nd part of the dipole
	// readint(file_id, "TDSE_inputs/PrintOutputMethod"		,&h5error,&inputs.PrintOutputMethod); // (0 - only text, 1 - only binaries, 2 - both)

	// readreal(file_id, "TDSE_inputs/trg_a"		,&h5error,&inputs.trg.a);
}

void Read_1_field_and_grid(hid_t file_id, char *inpath, herr_t *h5error, struct inputs_def *in)
{
	char path[50];
	printf("read 1d\n"); fflush(NULL);

  path[0] = '\0';	strcat(strcat(path,inpath),"tgrid");
	(*in).Efield.tgrid =  readreal1Darray_fort(file_id, path, h5error, &(*in).Efield.Nt); // tgrid is not changed when program runs
  path[0] = '\0';	strcat(strcat(path,inpath),"Field");
	(*in).Efield.Field =  readreal1Darray_fort(file_id, path, h5error, &(*in).Efield.Nt); // tgrid is not changed when program runs  
}


void print_nd_array_h5(hid_t file_id, char *dset_name, herr_t *h5error, int ndims, hsize_t *dimensions, void * array, hid_t datatype)
{
  hid_t dspace_id = H5Screate_simple(ndims, dimensions, NULL);
  hid_t dset_id = H5Dcreate2(file_id, dset_name, datatype, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  *h5error = H5Dwrite(dset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
  *h5error = H5Sclose(dspace_id);
  *h5error = H5Dclose(dset_id);	
}

void create_nd_array_h5(hid_t file_id, char *dset_name, herr_t *h5error, int ndims, hsize_t *dimensions, hid_t datatype)
{
  hid_t dspace_id = H5Screate_simple(ndims, dimensions, NULL);
  hid_t dset_id = H5Dcreate2(file_id, dset_name, datatype, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  *h5error = H5Sclose(dspace_id);
  *h5error = H5Dclose(dset_id);	
}


void PrintOutputs(hid_t file_id, char *inpath, herr_t *h5error, struct inputs_def *in, struct outputs_def *out)
{
	hsize_t output_dims[2]; // never exceeds 2 in this case, can be longer
	char path[50];
	printf("t1 \n"); fflush(NULL);

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
		print_nd_array_h5(file_id, path, h5error, 2, output_dims, (*out).FEfield_data, H5T_NATIVE_DOUBLE);
  }

  if ( (*in).Print.Fsourceterm == 1 )
  {
		path[0] = '\0';	strcat(strcat(path,inpath),"FSourceTerm");
		print_nd_array_h5(file_id, path, h5error, 2, output_dims, (*out).Fsourceterm_data, H5T_NATIVE_DOUBLE);
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


void rw_hyperslab_nd_h5(hid_t file_id, char *dset_name, herr_t *h5error, int nhyperslab_dimensions, int *hyperslab_dimensions_int, int *offset_int, int *count_int, void *array, char *rw) // This function reads full line from an n-D array, the selected dimension is given by (-1), the rest of selection is the offset
{ 
  int ndims, k1;
  hid_t datatype;
  hsize_t *dims; dims = get_dimensions_h5(file_id, dset_name, h5error, &ndims, &datatype);
  hsize_t hyperslab_dimensions[nhyperslab_dimensions];
  for(k1 = 0; k1 < nhyperslab_dimensions; k1++){hyperslab_dimensions[k1] = hyperslab_dimensions_int[k1];}
  hid_t memspace_id = H5Screate_simple(nhyperslab_dimensions,hyperslab_dimensions,NULL);
  hsize_t  offset[ndims], count[ndims];
  for(k1 = 0; k1 < ndims; k1++){offset[k1] = offset_int[k1]; count[k1] = count_int[k1]; //printf("offset: %i \n", offset[k1]); fflush(NULL); printf("count: %i \n", count[k1]); fflush(NULL);
	}

  hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT);
  hid_t dspace_id = H5Dget_space (dset_id);
  *h5error = H5Sselect_hyperslab (dspace_id, H5S_SELECT_SET, offset, NULL, count, NULL); // operation with only a part of the array = hyperslab	
  if (strcmp(rw,"r")==0){
    *h5error = H5Dread (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, array); // read only the hyperslab
  } else if (strcmp(rw,"w")==0){
    *h5error = H5Dwrite (dset_id, datatype, memspace_id, dspace_id, H5P_DEFAULT, array); // write the data
  } else {
    printf("wrongly sepcified r/w: nothing done\n"); 
  }
  
  *h5error = H5Dclose(dset_id); // dataset
  *h5error = H5Sclose(dspace_id); // dataspace
  *h5error = H5Sclose(memspace_id); // dataspace
}


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

void rw_real_full2Dhyperslab_nd_h5(hid_t file_id, char *dset_name, herr_t *h5error, int ndims, hsize_t *dimensions, int *selection, double *array1D, char *rw) // This function reads full line from an n-D array, the selected dimension is given by (-1), the rest of selection is the offset
{ 
  int k1, k2 = 0;
  hid_t memspace_id;
  hsize_t field_dims[2];  
  hsize_t  offset[ndims], stride[ndims], count[ndims], block[ndims];

  for(k1 = 0; k1 < ndims; k1++){
    stride[k1] = 1; block[k1] = 1;
    if (selection[k1] < 0){
      offset[k1] = 0;
      count[k1] = dimensions[k1];
      field_dims[k2] = dimensions[k1]; // a way to specify the length of the array for HDF5
      k2++;	
    }else{
      offset[k1] = selection[k1];
      count[k1] = 1;
    }
  }

  memspace_id = H5Screate_simple(2,field_dims,NULL);

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
