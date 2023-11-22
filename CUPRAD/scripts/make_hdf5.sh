#!/bin/bash

#SBATCH --job-name="MAKE_H5"
#SBATCH --ntasks=1
#SBATCH -t 0-00:01 # time (D-HH:MM)
#SBATCH -o MAKE_H5.%j.%N.out # STDOUT
#SBATCH -e MAKE_H5.%j.%N.err # STDERR
#SBATCH --profile=All

### Usage information
usage()
{
    echo "A script for converting parameter file into an HDF5 file.
usage: [[[-i (--inp) input file (.inp)] 
         [-o (--ohdf5) output HDF5 file] 
         [-s (--slurm) execute script within slurm task manager]] | 
         [-h help]]"
}

### Slurm switch
slurm()
{
    module purge
    source $MULTISCALE_HOME/Modules/load_python_modules.sh
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
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

if [ "$inp_filename" = "" ]
then
    echo "No input file (.inp) added. Use option [-i (--inp)]."
    exit 1
fi

if [ "$h5_filename" = "" ]
then
    echo "No output HDF5 file added. Use option [-o (--ohdf5)]."
    exit 1
fi


echo "Output hdf5 file: ${h5_filename}"
echo "Input file: ${inp_filename}"

python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $inp_filename -ohdf5 $h5_filename -g inputs

