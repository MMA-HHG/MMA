#!/bin/bash
#SBATCH -J tar_full
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output tar-%J.output

module purge
module load intelpython3/2019.3

./tar_scan3.sh

echo "Job finished" 