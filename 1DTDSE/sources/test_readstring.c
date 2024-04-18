#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "tools_hdf5.h"
#include "constants.h"

#define FILENAME "results.h5"
#define PATH1 "global_inputs/"
#define DATASET_NAME PATH1 "xxx"
#define DATASET_NAME2 PATH1 "gas_preset"


int main() {
    hid_t file_id, dataset_id, dataspace_id, datatype_id;
    herr_t h5error;
    hsize_t dims[1];
    char *data, *data2;

    file_id = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

    readstring(file_id, DATASET_NAME, &h5error, &data);
    readstring(file_id, DATASET_NAME2, &h5error, &data2);
    H5Fclose(file_id);

    // Print the string
    printf("String read from dataset: %s\n", data);
    printf("String2 read from dataset: %s\n", data2);

    printf("------ TEST tabulated ------\n");

    double a;

    Soft_Coulomb_parameters(data2, &a);

    printf("a = %f\n", a);


    return 0;
}