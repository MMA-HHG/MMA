#!/bin/bash
#SBATCH -J Hankel_parallel
#SBATCH --nodes=3
#SBATCH --constraint=compute  # HSW24|BDW28
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=32
#SBATCH --hint=multithread
#SBATCH --ntasks-per-core=2
#SBATCH --time=24:00:00
#SBATCH --output Hankel.output
module load python
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS path

