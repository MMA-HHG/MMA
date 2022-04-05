#!/bin/bash
#SBATCH -J CUPRAD
#SBATCH --nodes=2
#SBATCH --constraint=BDW28
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --time=06:00:00
#SBATCH --output CUPRAD-%J.output

module purge
module load intel intelmpi hdf5

#module load GCC/6.1.0 OpenMPI/1.10.2-GCC-6.1.0-psm-dbg HDF5/1.8.13-OpenMPI-1.10.2-GCC-6.1.0-psm FFTW/3.3.5-OpenMPI-1.10.2-GCC-6.1.0-psm-dbg
#module load intel/19.3 intelmpi/2019.3.199 hdf5/1.10.5-intel-19.0.3-intelmpi-2019.3.199 
#module load intel/18.1
#module load intelmpi/2018.1.163

echo 'CUPRAD code executing'

srun --mpi=pmi2 -K1 -n $SLURM_NTASKS $CUPRAD_BINARY/cuprad_occigen.e

echo 'CUPRAD code done'
