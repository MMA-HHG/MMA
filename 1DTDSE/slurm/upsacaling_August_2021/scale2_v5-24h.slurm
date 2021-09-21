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
fi_info
 export I_MPI_PIN_PROCESSOR_LIST=allcores:map=spread
 export SLURM_CPU_BIND=none

# CINES N.A
set -x
#export I_MPI_OFI_PROVIDER=ofi_rxm
#export I_MPI_FABRICS=shm:ofi
export I_MPI_FABRICS_LIST=ofa,dapl,tcp
unset I_MPI_DAPL_UD
export I_MPI_FALLBACK=enable
set +x
# fin CINES N.A

export  I_MPI_HBW_POLICY=hbw_bind,hbw_preferred,hbw_bind
export I_MPI_ADJUST_ALLTOALLV=2
export I_MPI_ADJUST_ALLTOALL=3

 export I_MPI_PIN=1
 export I_MPI_DEBUG=+6


ldd ./TDSEs_compilation1.e

cksum TDSEs_compilation1.e

mpirun -n $SLURM_NTASKS ./TDSEs_compilation1.e > logfile.log-$SLURM_JOBID
