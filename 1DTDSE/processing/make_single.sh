#!/bin/bash


python3 $TDSE_1D_HOME/prepare_fields/prepare_field_single.py -nodisplay

python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $TDSE_1D_HOME/processing/FreeFormInputs_TDSE_single.inp -ihdf5 results.h5 -ohdf5results.h5 -g TDSE_inputs 