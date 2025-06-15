#!/usr/bin/bash

#SBATCH --job-name "CUPRAD_pre_processor"
#SBATCH --account FTA-25-17
#SBATCH --partition qcpu
#SBATCH --time 00:05:00
#SBATCH --nodes 1
#SBATCH -o CUPRAD_pre_processor.%j.%N.out # STDOUT
#SBATCH -e CUPRAD_pre_processor.%j.%N.err # STDERR
#SBATCH --profile=All


$CUPRAD_SCRIPTS/pre_processor.sh