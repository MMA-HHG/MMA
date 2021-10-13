#!/bin/bash
#SBATCH -J Hankel-long-mp
#SBATCH --ntasks=1 
#SBATCH --constraint=BDW28
#SBATCH --cpus-per-task=28
#SBATCH --hint=multithread
#SBATCH --time=10:00:00
#SBATCH --output Hankel_mp-%J.output

module purge
module load intelpython3/2019.3

srun python $HANKEL_HOME/Hankel_long_medium_mp.py -here -nodisplay -input inputs_Hankel_all_cummulative_dr2.h5

echo "Job finished" 