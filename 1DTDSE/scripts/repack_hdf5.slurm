#!/bin/bash

#SBATCH --job-name="REPACK_HDF5"
#SBATCH --ntasks=1
#SBATCH -t 0-02:00 # time (D-HH:MM)
#SBATCH -o HDF5_MERGE.%j.%N.out # STDOUT
#SBATCH -e HDF5_MERGE.%j.%N.err # STDERR
#SBATCH --profile=All


h5repack -i results.h5 -o results_repacked.h5
