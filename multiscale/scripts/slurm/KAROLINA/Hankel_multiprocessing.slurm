#!/bin/bash
#SBATCH -J Hankel-long-mp
#SBATCH --account FTA-25-17
#SBATCH --partition qcpu
#SBATCH --nodes 1
#SBATCH --ntasks=1 
#SBATCH --hint=multithread
#SBATCH --time=72:00:00
#SBATCH --output Hankel_mp-%J.output

module purge
load_python_modules

# srun python $HANKEL_HOME/Hankel_long_medium_mp.py -here -nodisplay

# echo "No Hankel at the instant" 

srun python3 $HANKEL_HOME/Hankel_long_medium_parallel_cluster.py

echo "Job finished" 