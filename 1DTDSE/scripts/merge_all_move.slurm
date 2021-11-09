#!/bin/bash
#SBATCH -J TDSE_merge_and_move
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output TDSE_merge-%J.output

module purge
module load intelpython3/2019.3

srun python $TDSE_1D_POST_PROCESSING/merge_move.py

echo "Job finished" 