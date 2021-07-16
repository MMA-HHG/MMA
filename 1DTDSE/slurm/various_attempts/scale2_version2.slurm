#!/bin/bash
#SBATCH -J MPI-s2-v2
#SBATCH  --constraint=HSW24  
#SBATCH  --exclusive
#SBATCH  --nodes=20
#SBATCH  --ntasks=480
#SBATCH  --ntasks-per-node=24 
#SBATCH  --threads-per-core=1
#SBATCH  --time=24:00:00 
#SBATCH  --output=MPI_s2_v2-%J.out
module purge
module load intel/18.1 intelmpi/2018.1.163 hdf5/1.10.4 
module load fftw/3.3.8
# module load intel/19.4 intelmpi/2019.4.243
# module load intel/19.4 openmpi/intel/2.0.4
# module load gcc/8.3.0 openmpi/gnu/2.0.4

unset I_MPI_DAPL_DIRECT_COPY_THRESHOLD
export I_MPI_PIN_DOMAIN=auto:scatter
export SLURM_CPU_BIND=none
export I_MPI_PIN_DOMAIN=auto
export I_MPI_PIN_RESPECT_CPUSET=0
export I_MPI_DYNAMIC_CONNECTION=0 #Turn on/off the dynamic connection establishment.
export I_MPI_DAPL_UD=enable                   #( connectionless UD datagrams)
export I_MPI_FALLBACK=0

export I_MPI_ADJUST_ALLTOALLV=2
export I_MPI_ADJUST_ALLTOALL=3
export I_MPI_DEBUG=+6

ldd ./test_strided3.e

mpirun -n $SLURM_NTASKS ./test_strided3.e > logfile.log-$SLURM_JOBID