#!/bin/bash
#SBATCH --job-name="MPI_TDSE"
#SBATCH -t 0-02:00 # time (D-HH:MM)
#SBATCH --partition=MIX
#SBATCH -o MPI_TDSE.%j.%N.out # STDOUT
#SBATCH -e MPI_TDSE.%j.%N.err # STDERR
#SBATCH --profile=All

### Purge modules
module purge

### Load modules
load_modules

mpirun $TDSE_1D_HOME/build/TDSE.e
