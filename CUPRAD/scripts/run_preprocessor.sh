#!/bin/bash

#SBATCH --job-name="PREPROCESSOR"
#SBATCH --ntasks=1
#SBATCH -t 0-00:01 # time (D-HH:MM)
#SBATCH -o PREPROCESSOR.%j.%N.out # STDOUT
#SBATCH -e PREPROCESSOR.%j.%N.err # STDERR
#SBATCH --profile=All

### Usage information
usage()
{
    echo "A script for preprocessor execution.
usage: [[[-i (--ihdf5) input HDF5 file]
         [-s (--slurm) execute script within slurm task manager]] | 
         [-h help]]"
}

### Slurm switch
slurm()
{
    module purge
    load_modules
}

### Main
while [ "$1" != "" ]; do
    case $1 in
        -i | --ihdf5 )          shift
                                h5_filename="$1"
                                ;;
        -s | --slurm )          slurm
                                ;;
        -h | --help )           usage
                                return 0
                                ;;
        * )                     usage
                                return 1
    esac
    shift
done

echo "Processing ${h5_filename}"

$CUPRAD_HOME/build/make_start.e << INPUTS
${h5_filename}
0
0
0
INPUTS
