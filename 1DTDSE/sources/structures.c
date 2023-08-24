/**
 * @file structures.c
 * @brief Contains methods operating over defined structures.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "structures.h"
#include "util_hdf5.h"

void outputs_destructor(outputs_def *outputs) // frees memory allocated for outputs
{
	free((*outputs).tgrid);
	free((*outputs).omegagrid);
	free((*outputs).tgrid_fftw);

	free((*outputs).Efield);
	free((*outputs).sourceterm);
	free((*outputs).PopTot);

	free((*outputs).FEfield_data);
	free((*outputs).Fsourceterm_data);

	free((*outputs).FEfieldM2);
	free((*outputs).FsourcetermM2);

	free((*outputs).PopInt);
	free((*outputs).expval);

	(*outputs).Nt = 0;
	(*outputs).Nomega = 0;
}

void inputs_destructor(inputs_def *in) // frees memory allocated for inputs
{
	free((*in).psi0);
	free((*in).x);
	free((*in).timet);
	free((*in).dipole);
	free((*in).Efield.tgrid);
	free((*in).Efield.Field);
}

output_print_def Set_all_prints(void)
{
	output_print_def res;

	res.Efield = 1;
	res.FEfield = 1;
	res.sourceterm = 1;
	res.Fsourceterm = 1;
	res.FEfieldM2 = 1;
	res.FsourceTermM2 = 1;
	res.PopTot = 1;
	res.tgrid = 1;
	res.omegagrid = 1;
	res.PopInt = 1;
	res.expval_x = 1;

	return res;
}

output_print_def Initialise_Printing_struct(void)
{
	 output_print_def res;

	res.Efield = 0;
	res.FEfield = 0;
	res.sourceterm = 0;
	res.Fsourceterm = 0;
	res.FEfieldM2 = 0;
	res.FsourceTermM2 = 0;
	res.PopTot = 0;
	res.tgrid = 0;
	res.omegagrid = 0;
	res.PopInt = 0;
	res.expval_x = 0;

	return res;
}

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

	
	// res.FEfield = 1;
	// res.sourceterm = 1;
	// res.Fsourceterm = 1;
	// res.FEfieldM2 = 1;
	// res.FsourceTermM2 = 1;
	// res.PopTot = 1;

	res.tgrid = 1; // not memory consuming
	res.omegagrid = 1; // not memory consuming

	// res.PopInt = 1;
	// res.expval_x = 1;

	return res;
}