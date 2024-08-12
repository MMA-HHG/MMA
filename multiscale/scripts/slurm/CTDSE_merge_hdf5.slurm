#!/bin/bash

#SBATCH --job-name="HDF5_MRG"
#SBATCH --ntasks=1
#SBATCH -t 0-23:59 # time (D-HH:MM)
#SBATCH -o HDF5_MERGE.%j.%N.out # STDOUT
#SBATCH -e HDF5_MERGE.%j.%N.err # STDERR
#SBATCH --profile=All


module purge
load_python_modules


python3 $TDSE_1D_POST_PROCESSING/merge.py -keep-files