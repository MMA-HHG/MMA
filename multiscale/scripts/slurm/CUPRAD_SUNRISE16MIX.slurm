#!/bin/bash
# Your job name
#SBATCH --job-name=CUPRAD-TDSE-base
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH -t 0-03:00
#SBATCH --partition=MIX
#SBATCH --profile=All
#SBATCH --output MMA-CUPRAD-base-%J.output

load_modules

# echo of launched commands
set -x

srun $CUPRAD_HOME/build/cuprad.e