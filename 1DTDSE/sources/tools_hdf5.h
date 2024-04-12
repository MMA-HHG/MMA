/**
 * @file tools_hdf5.h
 * @brief Header containing HDF5 routines.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TOOLS_HDF5_H
#define TOOLS_HDF5_H

#include <time.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <math.h>
#include "structures.h"


void readreal(hid_t, char *, herr_t *, double *);
void readint(hid_t, char *, herr_t *, int *);
void rw_hyperslab_nd_h5(hid_t, char *, herr_t *, int, int*, int *, int *, void *, char *);
void rw_real_fullhyperslab_nd_h5(hid_t, char *, herr_t *, int, hsize_t *, int *, double *, char *);
void print_nd_array_h5(hid_t, char *, herr_t *, int, hsize_t *, void *, hid_t);
void create_nd_array_h5(hid_t, char *, herr_t *, int, hsize_t *, hid_t);
double * readreal1Darray_fort(hid_t, char *, herr_t *, int *);
hsize_t * get_dimensions_h5(hid_t, char *, herr_t *, int *, hid_t *);
void prepare_local_output_fixed_print_grids_h5(hid_t, char *, herr_t *, inputs_def *, outputs_def *, int, hsize_t *);
void ReadInputs(hid_t, char *, char *, herr_t *,  inputs_def *);
void Read_1_field_and_grid(hid_t, char *, herr_t *,  inputs_def *);
void PrintOutputs(hid_t, char *, herr_t *, inputs_def *, outputs_def *);
hid_t dtype_h5(char *);
output_print_def Set_prints_from_HDF5(hid_t, char *, herr_t *); // sequence of *.h-files is that it cannot be in util.h now
void print_local_output_fixed_h5(hid_t, char *, herr_t *, inputs_def *, outputs_def *, int, int, int);

#endif