#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "tools_hdf5.h"
#include "constants.h"

#define FILENAME "test_dset.h5"
#define DSETNAME "aaa"



int main() {
    hid_t file_id, dataset_id, dataspace_id, datatype_id;
    herr_t h5error;
    hsize_t dims[1];
    char *data, *data2;

    file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    double a = 1.2;
    hsize_t output_dims[1] = {1};
    print_nd_array_h5(file_id, DSETNAME, &h5error, 1, output_dims, &a, H5T_NATIVE_DOUBLE);

    H5Fclose(file_id);



    return 0;
}