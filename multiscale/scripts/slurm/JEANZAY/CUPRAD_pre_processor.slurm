#!/bin/bash

#SBATCH --job-name="CUPRAD_pre_processor"
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH -o CUPRAD_pre_processor.%j.%N.out # STDOUT
#SBATCH -e CUPRAD_pre_processor.%j.%N.err # STDERR
#SBATCH --profile=All


$CUPRAD_SCRIPTS/pre_processor.sh