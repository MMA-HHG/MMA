#!/bin/bash
#SBATCH -J 1DTDSE_MPI
#SBATCH  --constraint=HSW24  
#SBATCH  --exclusive
#SBATCH  --nodes=5
#SBATCH  --ntasks=120
#SBATCH  --ntasks-per-node=24 
#SBATCH  --threads-per-core=1
#SBATCH  --time=20:00:00 
#SBATCH  --output=cines-TDSEScheduler-%J.out
module purge
module load intel/18.1 intelmpi/2018.1.163 hdf5/1.10.4 
module load fftw/3.3.8

 export I_MPI_PIN_PROCESSOR_LIST=allcores:map=spread
 export SLURM_CPU_BIND=none
#
  export I_MPI_DAPL_PROVIDER=ofa-v2-mlx5_0-1s
 export I_MPI_DAPL_DIRECT_COPY_THRESHOLD=655360
# 
export  I_MPI_HBW_POLICY=hbw_bind,hbw_preferred,hbw_bind
export I_MPI_ADJUST_ALLTOALLV=2
export I_MPI_ADJUST_ALLTOALL=3

 export I_MPI_PIN=1
 export I_MPI_DEBUG=+6

ldd ./test_strided3.e

mpirun -n $SLURM_NTASKS ./test_strided3.e > logfile.log-$SLURM_JOBID
