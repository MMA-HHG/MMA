#!/bin/bash
#SBATCH -J MMA-CUPRAD-base
#SBATCH --nodes=1
#SBATCH --ntasks=32
# #SBATCH --mem=16G
# #SBATCH --cpus-per-task=4 # to allocate memory
#SBATCH --hint=nomultithread       # 1 MPI process per physical core (no hyperthreading)
#SBATCH --time=02:00:00
#SBATCH --output MMA-CUPRAD-base-%J.output

# module purge
# module load python
#module load hdf5
# module load intel-mkl
load_modules

# echo of launched commands
set -x

srun $CUPRAD_HOME/build/cuprad.e