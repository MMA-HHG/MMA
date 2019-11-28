/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 *   This example reads hyperslab from the SDS.h5 file
 *   created by h5_write.c program into two-dimensional
 *   plane of the three-dimensional array.
 *   Information about dataset in the SDS.h5 file is obtained.
 */

#include "hdf5.h"

#define H5FILE_NAME        
#define DATASETNAME "IntArray"
#define NX_SUB  3           /* hyperslab dimensions */
#define NY_SUB  4
#define NX 7           /* output buffer dimensions */
#define NY 7
#define NZ  3
#define RANK         2
#define RANK_OUT     3

int
main (void)
{
    hid_t       file, dataset;         /* handles */
    hid_t       datatype, dataspace;
    hid_t       memspace;
    H5T_class_t t_class;                 /* data type class */
    H5T_order_t order;                 /* data order */
    size_t      size;                  /*
				        * size of the data element
				        * stored in file
				        */
    hsize_t     dimsm[3];              /* memory space dimensions */
    hsize_t     dims_out[2];           /* dataset dimensions */
    herr_t      status;

    int         data_out[NX][NY][NZ ]; /* output buffer */

    hsize_t      count[2];              /* size of the hyperslab in the file */
    hsize_t      offset[2];             /* hyperslab offset in the file */
    hsize_t      count_out[3];          /* size of the hyperslab in memory */
    hsize_t      offset_out[3];         /* hyperslab offset in memory */
    int          i, j, k, status_n, rank;



    for (j = 0; j < NX; j++) {
	for (i = 0; i < NY; i++) {
	    for (k = 0; k < NZ ; k++)
		data_out[j][i][k] = 0;
	}
    }

    /*
     * Open the file and the dataset.
     */
    file = H5Fopen("data.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen2(file, "micro/params/dt", H5P_DEFAULT);

    /*
     * Get datatype and dataspace handles and then query
     * dataset class, order, size, rank and dimensions.
     */


    datatype  = H5Dget_type(dataset);     /* datatype handle */
    t_class     = H5Tget_class(datatype);
    if (t_class == H5T_INTEGER) printf("Data set has INTEGER type \n");
    order     = H5Tget_order(datatype);
    if (order == H5T_ORDER_LE) printf("Little endian order \n");

    
    size  = H5Tget_size(datatype);
    printf(" Data size is %d \n", (int)size);

/*    memspace = H5Screate_simple(1,1,NULL);*/


/*    status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace,*/
/*		     H5P_DEFAULT, dt);*/

    int pdumint[0];
    double dt[0];
/*    status = H5Dread(dataset,  H5T_STD_I64LE, H5S_ALL, H5S_ALL,*/
/*		     H5P_DEFAULT, pdumint);*/

    status = H5Dread(dataset,  datatype, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, dt);

    printf(" Value is %lf \n", dt[0]);


    H5Tclose(datatype);
    H5Dclose(dataset);
/*    H5Sclose(dataspace);*/
/*    H5Sclose(memspace);*/
    H5Fclose(file);

    return 0;
}
