#!/bin/bash
#SBATCH -J Hankel_testing
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --output Hankel_testing-%J.output

module purge
module load intelpython3/2019.3

python $HANKEL_HOME/Hankel_from_TDSE_z_planes_sum_diffraction2.py -here -nodisplay