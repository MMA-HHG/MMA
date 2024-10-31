#!/bin/bash
#SBATCH -J Hankel-long-mp
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24
#SBATCH --hint=multithread
#SBATCH -t 0-00:59
#SBATCH --partition=T1H
#SBATCH --output Hankel_mp-%J.output

module purge
load_python_modules

# srun python $HANKEL_HOME/Hankel_long_medium_mp.py -here -nodisplay

# echo "No Hankel at the instant" 

srun python3 $HANKEL_HOME/Hankel_long_medium_parallel_cluster.py

echo "Job finished" 