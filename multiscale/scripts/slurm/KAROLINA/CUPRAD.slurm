#!/bin/bash
#SBATCH -J "MMA-CUPRAD-base"
#SBATCH --account FTA-25-17
#SBATCH --partition qcpu
#SBATCH --nodes 1
#SBATCH --ntasks=64
#SBATCH --hint=nomultithread       # 1 MPI process per physical core (no hyperthreading)
#SBATCH --time=02:00:00
#SBATCH --output MMA-CUPRAD-base-%J.output
#SBATCH --profile=All


load_modules

# echo of launched commands
set -x

srun $CUPRAD_HOME/build/cuprad.e