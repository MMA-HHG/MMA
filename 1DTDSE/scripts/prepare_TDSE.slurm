#!/bin/bash
#SBATCH -J TDSE_prepare
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output TDSE_prepare-%J.output

module purge
module load intelpython3/2019.3

srun python $TDSE_1D_POST_PROCESSING/prepare_TDSE.py

echo "Job finished" 
