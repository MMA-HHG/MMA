#!/bin/bash
#SBATCH --job-name="MPI_TDSE"
#SBATCH --time=20:00:00
# #SBATCH --qos=qos_cpu-t4 # name Subscribe to a Quality of service
# #SBATCH --partition=T48H
#SBATCH -o MPI_TDSE.%j.%N.out # STDOUT
#SBATCH -e MPI_TDSE.%j.%N.err # STDERR
#SBATCH --profile=All

### Purge modules
module purge

load_modules

mpirun $TDSE_1D_HOME/build/TDSE.e
