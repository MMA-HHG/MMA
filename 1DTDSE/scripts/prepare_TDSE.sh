#!/bin/bash

#SBATCH --job-name="PREP_HDF5"
#SBATCH --ntasks=1
#SBATCH -t 0-00:01 # time (D-HH:MM)
#SBATCH -o PREPARE_TDSE.%j.%N.out # STDOUT
#SBATCH -e PREPARE_TDSE.%j.%N.err # STDERR
#SBATCH --profile=All

### Usage information
usage()
{
    echo "A script for converting parameter file into an HDF5 file.
usage: [[[-i (--inp) input parameter file (.inp)] 
         [-o (--ohdf5) output HDF5 file] 
         [-s (--slurm) execute script within slurm task manager]] | 
         [-h help]]"
}

### Slurm switch
slurm()
{
    module purge
    load_python_modules
}


### Main

while [ "$1" != "" ]; do
    case $1 in
        -i | --inp )            shift
                                inp_filename="$1"
                                ;;
        -o | --ohdf5 )          shift
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

if [ "$inp_filename" = "" ]
then
    echo "No input file (.inp) added. Use option [-i (--inp)]."
    return 1
fi

if [ "$h5_filename" = "" ]
then
    echo "No output HDF5 file added. Use option [-o (--ohdf5)]."
    return 1
fi


echo "Output hdf5 file: ${h5_filename}"
echo "Input file: ${inp_filename}"

python3 $TDSE_1D_POST_PROCESSING/prepare_TDSE.py -i $inp_filename -o $h5_filename

