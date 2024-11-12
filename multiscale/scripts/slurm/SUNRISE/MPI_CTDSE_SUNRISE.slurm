#!/bin/bash
#SBATCH --job-name="MPI_TDSE"
#SBATCH -t 1-23:59 # time (D-HH:MM)
#SBATCH --partition=T48H
#SBATCH -o MPI_TDSE.%j.%N.out # STDOUT
#SBATCH -e MPI_TDSE.%j.%N.err # STDERR
#SBATCH --profile=All

### Purge modules
module purge

load_modules

mpirun socket $TDSE_1D_HOME/build/TDSE.e # numa

# possibilities to test different options for making MPI run more efficient https://www.open-mpi.org/doc/v3.0/man1/mpirun.1.php
# mpirun --bind-to socket $TDSE_1D_HOME/build/TDSE.e # numa