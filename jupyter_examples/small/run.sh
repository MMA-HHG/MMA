#!/bin/bash

# run CUPRAD pre-processor
var=$(echo *.h5) # use the h5 file within the directory
$CUPRAD_HOME/build/make_start.e <<INPUTS # run the pre-processor with inputs
$var
0
0
0
INPUTS                               
mpirun -n 4 --allow-run-as-root $CUPRAD_BUILD/cuprad.e      # run CUPRAD
python3 $TDSE_1D_PYTHON/prepare_TDSE_Nz.py                  # add correct number of planes in the medium
mpirun -n 4 --allow-run-as-root $TDSE_1D_BUILD/TDSE.e       # run TDSE
python3 $TDSE_1D_PYTHON/merge.py                            # merge TDSE results
python3 $HANKEL_HOME/Hankel_long_medium_parallel_cluster.py         # run Hankel
python3 $HANKEL_HOME/copy_results_tom_main.py                       # copy Hankel results to main file
