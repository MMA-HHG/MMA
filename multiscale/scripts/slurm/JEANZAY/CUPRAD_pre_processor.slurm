#!/bin/bash

#SBATCH --job-name="CUPRAD_pre_processor"
#SBATCH --ntasks=1
#SBATCH -t 05:00 # time (HH:MM)
#SBATCH -o CUPRAD_pre_processor.%j.%N.out # STDOUT
#SBATCH -e CUPRAD_pre_processor.%j.%N.err # STDERR
#SBATCH --profile=All


$CUPRAD_SCRIPTS/pre_processor.sh