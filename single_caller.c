#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<hdf5.h>

#include "numerical_constants.h"
#include "util.h"
#include "util_hdf5.h"

// hdf5 operation:
herr_t  h5error;
hid_t file_id; // file pointer
hid_t filespace, dataspace_id, dataset_id, dset_id, dspace_id; // dataspace pointers

struct inputs_def inputs;
struct outputs_def outputs;

int k1, k2, k3;


clock_t start_clock, end_clock;


int main() 
{
	// vars:
	// dummy
	int dum3int[3];
	hsize_t * dims; int ndims; hid_t datatype;
	// Processing the queue


	int comment_operation = 1;

	///////////////////////
	// PREPARATION PHASE //
	///////////////////////

	printf("program started\n"); fflush(NULL);
	start_clock = clock(); // the clock	
	Init_constants();
	inputs.Print = Initialise_Printing_struct(); // crete printing driver


	// read inputs

	/* to check if exists use printf("link exists 1: %i\n",H5Lexists(file_id, "IRProp/lambda", H5P_DEFAULT)); */
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.
	ReadInputs(file_id, "TDSE_inputs/", &h5error, &inputs);
	h5error = H5Fclose(file_id); 

	// convert units
	// for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]*1e-15/TIMEau; inputs.Efield.Field[k1] = inputs.Efield.Field[k1]*1e9/EFIELDau;} // convert to atomic units (fs->a.u.), (GV/m->a.u.)
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; /*inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;*/} // convert to atomic units (fs->a.u.), (GV/m->a.u.)

	// Prepare the ground state

	inputs.CV = 1E-20; 
	Initialise_grid_and_ground_state(&inputs);
	printf("Initial energy is : %1.12f\n",inputs.Einit); fflush(NULL);
	printf("xgrid, psi0 : %e %e %e %e %e %e\n", inputs.x[0],inputs.x[1],inputs.x[2],inputs.psi0[0],inputs.psi0[1],inputs.psi0[2]); fflush(NULL);

	/////////////////////////
	// COMPUTATIONAL PHASE //
	/////////////////////////


	outputs = call1DTDSE(inputs); // THE TDSE
	printf("TDSE done, in the caller\n"); fflush(NULL);


	/////////////////////
	// PROCESS OUTPUTS //
	/////////////////////


	printf("Printing the outputs \n"); fflush(NULL);
	//outputs.FEfield = create_2Darray_accessor_real(mydims, outputs.FEfield_data); // in the case we would like to access the array usual way
	//outputs.Fsourceterm = create_2Darray_accessor_real(mydims, outputs.Fsourceterm_data);
	//printf("efield out    : %e, %e, %e \n",outputs.Efield[0],outputs.Efield[1],outputs.Efield[2]); fflush(NULL);


	hsize_t output_dims[2]; // never exceeds 2 in this case, can be longer

	file_id = H5Fopen ("results2.h5", H5F_ACC_RDWR, H5P_DEFAULT); // we use a different output file to testing, can be changed to have only one file
	hid_t g_id = H5Gcreate2(file_id, "/TDSEsingle", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	h5error = H5Gclose(g_id);

	inputs.Print = Set_all_prints();
	PrintOutputs(file_id, "/TDSEsingle/", &h5error, &inputs, &outputs);

	// time domain
	output_dims[0] = outputs.Nt; output_dims[1] = 0;
	print_nd_array_h5(file_id, "/TDSEsingle/tgrid_fftw", &h5error, 1, output_dims, outputs.tgrid_fftw, H5T_NATIVE_DOUBLE);
	outputs_destructor(&outputs);

	// durations of calculation
	output_dims[0] = 1; output_dims[1] = 0;	
	end_clock = clock();
	double elapsed_time = (double)((end_clock - start_clock)/CLOCKS_PER_SEC);
	print_nd_array_h5(file_id, "/TDSEsingle/full_runtime", &h5error, 1, output_dims, &elapsed_time, H5T_NATIVE_DOUBLE);


	h5error = H5Fclose(file_id);
	





	// int mydims[2] = {outputs.Nomega,2};

	

	// //

	
    //     //printf("sourceterm out: %e, %e, %e \n",outputs.sourceterm[0],outputs.sourceterm[1],outputs.sourceterm[2]);
    //     printf("efield out    : %e, %e, %e \n",outputs.Efield[0],outputs.Efield[1],outputs.Efield[2]);
	// 	printf("\nFefield out    : \n%e, %e \n%e, %e \n%e, %e \n",outputs.FEfield[0][0],outputs.FEfield[0][1],outputs.FEfield[1][0],outputs.FEfield[1][1],outputs.FEfield[2][0],outputs.FEfield[2][1]);

	// 	printf("\nFSourceTerm: \n%e, %e,\n%e, %e \n%e, %e \n",
	// 	outputs.Fsourceterm[0][0],outputs.Fsourceterm[0][1],
	// 	outputs.Fsourceterm[1][0],outputs.Fsourceterm[1][1],
	// 	outputs.Fsourceterm[2][0],outputs.Fsourceterm[2][1]);

	// 	printf("FfieldM2      : %e, %e, %e \n",outputs.FEfieldM2[0],outputs.FEfieldM2[1],outputs.FEfieldM2[2]);
	// 	printf("FSourceTermM2 : %e, %e, %e \n",outputs.FsourcetermM2[0],outputs.FsourcetermM2[1],outputs.FsourcetermM2[2]);
	// 	printf("Nomega        : %i \n",outputs.Nomega);

    //     // prepare the dataset(s) for outputs

    //     // dims[0] = outputs.Nt; // length defined by outputs
        
	// hsize_t dims3[2]; dims3[0] = outputs.Nomega; dims3[1] = 2;

	// print_nd_array_h5(file_id, "/test4", &h5error, 2, dims3, outputs.FEfield_data, H5T_NATIVE_DOUBLE); // https://support.hdfgroup.org/HDF5/doc1.6/PredefDTypes.html
	
	// hsize_t dims2[1]; dims2[0] = outputs.Nt;
	// print_nd_array_h5(file_id, "/test", &h5error, 1, dims2, outputs.Efield, H5T_NATIVE_DOUBLE); // https://support.hdfgroup.org/HDF5/doc1.6/PredefDTypes.html

	
	// double myarray[outputs.Nomega][2];
	// for(k1 = 0; k1 < outputs.Nomega;k1++){myarray[k1][0] = outputs.FEfield[k1][0];myarray[k1][1] = outputs.FEfield[k1][1];}
	// //print_nd_array_h5(file_id, "/test2", &h5error, 2, dims3, (double*)((*outputs.FEfield) + outputs.Nomega), H5T_NATIVE_DOUBLE); // https://support.hdfgroup.org/HDF5/doc1.6/PredefDTypes.html
	// print_nd_array_h5(file_id, "/test2", &h5error, 2, dims3, myarray, H5T_NATIVE_DOUBLE); // https://support.hdfgroup.org/HDF5/doc1.6/PredefDTypes.html

	// double *myarray2;
	// myarray2 = (double*) malloc(2*outputs.Nomega*sizeof(double));
	// for(k1 = 0; k1 < outputs.Nomega;k1++){myarray2[2*k1] = outputs.FEfield[k1][0]; myarray2[2*k1+1] = outputs.FEfield[k1][1];}
	// print_nd_array_h5(file_id, "/test3", &h5error, 2, dims3, myarray2, H5T_NATIVE_DOUBLE); // https://support.hdfgroup.org/HDF5/doc1.6/PredefDTypes.html

	

	// double **array_accessor;
	// array_accessor = (double**) malloc(outputs.Nomega*sizeof(double));
	// for(k1 = 0; k1 < outputs.Nomega;k1++){array_accessor[k1] = &myarray2[2*k1];}
	// printf("\narray accessor    : \n%e, %e \n%e, %e \n%e, %e \n",array_accessor[0][0],array_accessor[0][1],array_accessor[1][0],array_accessor[1][1],array_accessor[2][0],array_accessor[2][1]);


	// array_accessor[0][0] = 1.6;
	// printf("changed array : %e \n",array_accessor[0][0]);

    //     // dataspace_id = H5Screate_simple(ndims, dims, NULL); // create dataspace for outputs
    //     // dataset_id = H5Dcreate2(file_id, "/SourceTerms", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // create dataset
    //     // h5error = H5Sclose(dataspace_id);
    //     // h5error = H5Dclose(dataset_id);
    //     // rw_real_fullhyperslab_nd_h5(file_id,"/SourceTerms",&h5error,3,dims,dum3int,outputs.Efield,"w");



   	
	



	// hsize_t dims4[1]; dims4[0] = outputs.Nomega;



   
    
    printf("Done \n");

}
