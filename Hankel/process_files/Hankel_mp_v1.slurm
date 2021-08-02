#!/bin/bash
#SBATCH -J Hankel-long
#SBATCH --nodes=1
#SBATCH --constraint=BDW28
#SBATCH --ntasks=28
#SBATCH --time=10:00:00
#SBATCH --output Hankel_mp-%J.output

module purge
module load intelpython3/2019.3

srun python $HANKEL_HOME/Hankel_long_medium.py -here -nodisplay

echo "Job finished" 