#!/bin/bash

#SBATCH --job-name="h5_repack_gzip"
#SBATCH --ntasks=1
#SBATCH --partition=T24H
#SBATCH -t 0-04:00 # time (D-HH:MM)
#SBATCH -o h5_repack+gzip.%j.%N.out # STDOUT
#SBATCH -e h5_repack+gzip.%j.%N.err # STDERR
#SBATCH --profile=All


h5file=$(cat msg.tmp)
h5file_new=${h5file%.h5}_gzip.h5

# echo $h5file_new

h5repack -f SHUF -f GZIP=9 $h5file $h5file_new

echo "Job finished" 