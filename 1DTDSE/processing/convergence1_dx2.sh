#!/bin/bash


python3 $TDSE_1D_HOME/prepare_fields/prepare_fields.py -nodisplay -g fields_list -ohdf5 fields_list.h5 -g_params grids_for_scans -i $TDSE_1D_HOME/prepare_fields/TDSE_convergence1.inp

python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $TDSE_1D_HOME/processing/TDSE_convergence1_dx2.inp -ihdf5 fields_list.h5 -ohdf5 results.h5 -g TDSE_inputs 

cp results.h5 results2.h5