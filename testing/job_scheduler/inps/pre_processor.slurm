#!/bin/bash
#SBATCH -J pre_processor
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output pre_processor-%J.output

# module purge
# module load intel intelmpi hdf5

#module load GCC/6.1.0 OpenMPI/1.10.2-GCC-6.1.0-psm-dbg HDF5/1.8.13-OpenMPI-1.10.2-GCC-6.1.0-psm FFTW/3.3.5-OpenMPI-1.10.2-GCC-6.1.0-psm-dbg
#module load intel/19.3 intelmpi/2019.3.199 hdf5/1.10.5-intel-19.0.3-intelmpi-2019.3.199 
#module load intel/18.1
#module load intelmpi/2018.1.163

$TESTPATH/pre_processor.sh

# srun --mpi=pmi2 -K1 -n $SLURM_NTASKS /home/vabekjan/git/CUPRAD_DEVELOP/binary/cuprad_occigen.e
