#!/bin/bash  
#SBATCH --account=vabek-acc1
#SBATCH --partition=queue48h
#SBATCH -N 10 # number of nodes 
#SBATCH --clusters=eclipse
#SBATCH -t 2-0:00 # time (D-HH:MM) 
#SBATCH -o slurm.%N.%j.out # STDOUT 
#SBATCH -e slurm.%N.%j.err # STDERR 
#SBATCH -J TDSE_Maker2



# module load OpenMPI/1.8.4-GCC-4.8.4
# mpirun -npernode 2  ./stream_mpi_test

module purge
module load GCC/6.1.0 OpenMPI/1.10.2-GCC-6.1.0-psm-dbg HDF5/1.8.13-OpenMPI-1.10.2-GCC-6.1.0-psm FFTW/3.3.5-OpenMPI-1.10.2-GCC-6.1.0-psm-dbg


mpirun -n 160 $TDSE_1D_HOME/processing/TDSEs_list_compilation1.e


echo "Done"


