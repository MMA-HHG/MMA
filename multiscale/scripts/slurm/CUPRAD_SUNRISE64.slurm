#!/bin/bash
# Your job name
#SBATCH --job-name=CUPRAD-Bessel
#SBATCH --ntasks=128
#SBATCH -t 0-23:59
#SBATCH --partition=T24H
#SBATCH --profile=All
#SBATCH --output MMA-CUPRAD-base-%J.output

load_modules

# echo of launched commands
set -x

srun $CUPRAD_HOME/build/cuprad.e