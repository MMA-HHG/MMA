#!/bin/bash
#SBATCH -Jspectral-analyses
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --output spectral_analyses-%J.output

module purge
module load intelpython3/2019.3

python $CUPRAD_PYTHON/spectral_analyses.py -here -nodisplay
