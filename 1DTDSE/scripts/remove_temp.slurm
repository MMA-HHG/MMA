#!/bin/bash
#SBATCH -J remove_temp
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=00:01:30
#SBATCH --output remove_temp-%J.output

module purge
module load intelpython3/2019.3

rm -r temp

echo "Job finished" 