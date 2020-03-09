#!/bin/bash
#SBATCH -J Hankel_parallel
#SBATCH --ntasks=1 
#SBATCH --constraint=BDW28      # SKL224
#SBATCH --cpus-per-task=28
#SBATCH --hint=multithread
#SBATCH --time=2:00:00
#SBATCH --output Hankel.output
module purge
module load intelpython3/2019.3
srun python Hankel_Ilist_and_phenom.py

echo "Job finished" 
