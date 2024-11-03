#!/bin/bash

#SBATCH --job-name="MPI-TDSE-prepare"
#SBATCH --ntasks=1
#SBATCH --time=00:02:00
#SBATCH -o MPI-TDSE-prepare.%j.%N.out # STDOUT
#SBATCH -e MPI-TDSE-prepare.%j.%N.err # STDERR
#SBATCH --profile=All



module purge
load_python_modules


python3 $TDSE_1D_PYTHON/prepare_TDSE_Nz.py