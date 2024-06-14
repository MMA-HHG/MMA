#!/bin/bash
#SBATCH --job-name="MPI_TDSE"
#SBATCH -t 1-00:00 # time (D-HH:MM)
#SBATCH --partition=T24H
#SBATCH -o MPI_TDSE.%j.%N.out # STDOUT
#SBATCH -e MPI_TDSE.%j.%N.err # STDERR
#SBATCH --profile=All

### Purge modules
module purge

load_modules

mpirun $TDSE_1D_HOME/build/TDSE.e
