#!/bin/bash
#SBATCH --job-name="MPI_TDSE"
#SBATCH --time=24:00:00
#SBATCH --account FTA-25-17
#SBATCH --partition qcpu
#SBATCH -o MPI_TDSE.%j.%N.out # STDOUT
#SBATCH -e MPI_TDSE.%j.%N.err # STDERR
#SBATCH --profile=All

### Purge modules
module purge

load_modules

mpirun $TDSE_1D_HOME/build/TDSE.e
