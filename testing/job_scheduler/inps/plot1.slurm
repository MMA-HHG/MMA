#!/bin/bash
#SBATCH -J plot1
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output plot1-%J.output

module purge
module load intelpython3/2019.3

python $PYTHON_CUPRAD/compare_n_results_XUVshift.py -here -nodisplay