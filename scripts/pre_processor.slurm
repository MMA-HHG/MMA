#!/bin/bash
#SBATCH -J pre_processor
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output pre_processor-%J.output

$CUPRAD_SCRIPTS/pre_processor.sh