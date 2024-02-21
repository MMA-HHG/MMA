#!/bin/bash
# Your job name
#SBATCH --job-name=CUPRAD-TDSE-base
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH -t 0-03:00
#SBATCH --partition=T12H
#SBATCH --profile=All
#SBATCH --output MMA-CUPRAD-base-%J.output

# #SBATCH -J CUPRAD-density-mod
# #SBATCH --nodes=1
# #SBATCH --ntasks=32
# # #SBATCH --mem=16G
# # #SBATCH --cpus-per-task=4 # to allocate memory
# #SBATCH --hint=nomultithread       # 1 MPI process per physical core (no hyperthreading)
# #SBATCH --time=2:00:00
# #SBATCH --output CUPRAD_density_mod-%J.output


# module purge
# module load python
# module load hdf5
# module load intel-mkl

load_modules

# echo of launched commands
set -x

srun $CUPRAD_HOME/build/cuprad.e