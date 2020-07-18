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

	// READ DATA
	/* to check if exists use printf("link exists 1: %i\n",H5Lexists(file_id, "IRProp/lambda", H5P_DEFAULT)); */
	file_id = H5Fopen ("results.h5", H5F_ACC_RDONLY, H5P_DEFAULT); // the file is opened for read only by all the processes independently, every process then has its own copy of variables.

	ReadInputs(file_id, "TDSE_inputs/", &h5error, &inputs);
	printf("Initial guess: %e\n", inputs.Eguess); fflush(NULL);	
	
	printf("Read tested\n"); fflush(NULL);
	
// 	// return 0;
// 	readreal(file_id, "TDSE_inputs/Eguess"					,&h5error,&inputs.Eguess); // Energy of the initial state
// 	readint(file_id, "TDSE_inputs/N_r_grid"					,&h5error,&inputs.num_r); // Number of points of the initial spatial grid 16000
// 	readint(file_id, "TDSE_inputs/N_r_grid_exp"				,&h5error,&inputs.num_exp); // Number of points of the spatial grid for the expansion
// 	readreal(file_id, "TDSE_inputs/dx"						,&h5error,&inputs.dx); // resolution for the grid
// 	readint(file_id, "TDSE_inputs/InterpByDTorNT"			,&h5error,&inputs.InterpByDTorNT); 
// 	readreal(file_id, "TDSE_inputs/dt"						,&h5error,&inputs.dt); // resolution in time
// 	readint(file_id, "TDSE_inputs/Ntinterp"					,&h5error,&inputs.Ntinterp); // Number of points of the spatial grid for the expansion
// 	readreal(file_id, "TDSE_inputs/textend"					,&h5error,&inputs.textend); // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
// 	readint(file_id, "TDSE_inputs/analy_writewft"			,&h5error,&inputs.analy.writewft); // writewavefunction (1-writting every tprint)
// 	readreal(file_id, "TDSE_inputs/analy_tprint"			,&h5error,&inputs.analy.tprint); // time spacing for writing the wavefunction	
// 	readreal(file_id, "TDSE_inputs/x_int"					,&h5error,&inputs.x_int); // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
// 	readint(file_id, "TDSE_inputs/PrintGaborAndSpectrum"	,&h5error,&inputs.PrintGaborAndSpectrum); // print Gabor and partial spectra (1-yes)
// 	readreal(file_id, "TDSE_inputs/a_Gabor"					,&h5error,&inputs.a_Gabor); // the parameter of the gabor window [a.u.]
// 	readreal(file_id, "TDSE_inputs/omegaMaxGabor"			,&h5error,&inputs.omegaMaxGabor); // maximal frequency in Gabor [a.u.]
// 	readreal(file_id, "TDSE_inputs/dtGabor"					,&h5error,&inputs.dtGabor); // spacing in Gabor
// 	readreal(file_id, "TDSE_inputs/tmin1window"				,&h5error,&inputs.tmin1window); // analyse 1st part of the dipole
// 	readreal(file_id, "TDSE_inputs/tmax1window"				,&h5error,&inputs.tmax1window); // analyse 1st part of the dipole
// 	readreal(file_id, "TDSE_inputs/tmin2window"				,&h5error,&inputs.tmin2window); // analyse 2nd part of the dipole
// 	readreal(file_id, "TDSE_inputs/tmax2window"				,&h5error,&inputs.tmax2window); // analyse 2nd part of the dipole
// 	readint(file_id, "TDSE_inputs/PrintOutputMethod"		,&h5error,&inputs.PrintOutputMethod); // (0 - only text, 1 - only binaries, 2 - both)
// //	readint(file_id, "TDSE_inputs/IonisationFilterForTheSourceTerm"	,&h5error,&inputs.IonisationFilterForTheSourceTerm); // filter source term by high-ionisation components (1-yes)
// //	readreal(file_id, "TDSE_inputs/IonFilterThreshold"		,&h5error,&inputs.IonFilterThreshold); // threshold for the ionisation [-]
// 	readreal(file_id, "TDSE_inputs/trg_a"		,&h5error,&inputs.trg.a);


	// load the tgrid
	inputs.Efield.tgrid =  readreal1Darray_fort(file_id, "IRField/tgrid",&h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs
	inputs.Efield.Field =  readreal1Darray_fort(file_id, "IRField/Field",&h5error,&inputs.Efield.Nt); // tgrid is not changed when program runs

	h5error = H5Fclose(file_id); // file

	// for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]*1e-15/TIMEau; inputs.Efield.Field[k1] = inputs.Efield.Field[k1]*1e9/EFIELDau;} // convert to atomic units (fs->a.u.), (GV/m->a.u.)
	for(k1 = 0 ; k1 < inputs.Efield.Nt; k1++){inputs.Efield.tgrid[k1] = inputs.Efield.tgrid[k1]/TIMEau; inputs.Efield.Field[k1] = inputs.Efield.Field[k1]/EFIELDau;} // convert to atomic units (fs->a.u.), (GV/m->a.u.)

	// crete printing driver
	inputs.Print = Initialise_Printing_struct();

	// Prepare the ground state

	// Initialise vectors and Matrix 
	// Initialise_GS(inputs.num_r);
	
	//normalise(psi0,inputs.num_r); // Initialise psi0 for Einitialise
	// printf("test\n"); fflush(NULL);

	// double CV = 1E-20; // CV criteria

	// /* This number has to be small enough to assure a good convregence of the wavefunction
	// if it is not the case, then the saclar product of the the ground state and the excited states 
	// is not quite 0 and those excited appears in the energy analysis of the gorund states, so the propagation !!
	// CV = 1E-25 has been choosen to have a scalar product of 10^-31 with the third excited state for num_r = 5000 and dx=0.1
	// */
	
	// printf("Calculation of the energy of the ground sate ; Eguess : %f\n",inputs.Eguess); fflush(NULL);
	// int size = 2*(inputs.num_r+1);
	// double *off_diagonal, *diagonal;
	// double Einit = 0.0;
	// inputs.psi0 = calloc(size,sizeof(double));
	// for(k1=0;k1<=inputs.num_r;k1++){inputs.psi0[2*k1] = 1.0; inputs.psi0[2*k1+1] = 0.;}
	// Initialise_grid_and_D2(inputs.dx, inputs.num_r, &inputs.x, &diagonal, &off_diagonal); // !!!! dx has to be small enough, it doesn't converge otherwise
	// Einit = Einitialise(inputs.trg,inputs.psi0,off_diagonal,diagonal,off_diagonal,inputs.x,inputs.Eguess,CV,inputs.num_r); // originally, some possibility to have also excited state

	// printf("Initial energy is : %1.12f\n",Einit); fflush(NULL);
	// printf("xgrid, psi0 : %e %e %e %e %e %e\n", inputs.x[0],inputs.x[1],inputs.x[2],inputs.psi0[0],inputs.psi0[1],inputs.psi0[2]); fflush(NULL);

	// free(inputs.x); free(inputs.psi0);

	inputs.CV = 1E-20; 	
	Initialise_grid_and_ground_state(&inputs, &Einit);
	printf("Initial energy is : %1.12f\n",Einit); fflush(NULL);
	printf("xgrid, psi0 : %e %e %e %e %e %e\n", inputs.x[0],inputs.x[1],inputs.x[2],inputs.psi0[0],inputs.psi0[1],inputs.psi0[2]); fflush(NULL);

	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////


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
