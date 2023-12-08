#!/bin/bash

#SBATCH --job-name="HDF5_MRG"
#SBATCH --ntasks=1
#SBATCH -t 0-02:00 # time (D-HH:MM)
#SBATCH -o HDF5_MERGE.%j.%N.out # STDOUT
#SBATCH -e HDF5_MERGE.%j.%N.err # STDERR
#SBATCH --profile=All

### Usage information
usage()
{
    echo "A script for merging the temporary HDF5 files.
usage: [[[-p (--printdata) select data to print] 
         [-s (--slurm) execute script within slurm task manager]] | 
         [-h help]]

--------------------------------------------
Run example (for Slurm): source merge_hdf5.sh -p ['Efield SourceTerm ...'] -s

Print data available options: 'Efield', 'FEfield', 'SourceTerm', 'FSourceTerm', 
                              'FEfieldM2', 'FSourceTermM2', 'PopTot', 'PopInt', 
                              'expval_x'"
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
        -p | --printdata )      shift
                                printdata="$1"
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

if [ "$printdata" = "" ]
then
    printdata="Efield FEfield SourceTerm FSourceTerm FEfieldM2 FSourceTermM2 PopTot PopInt expval_x"
fi

python3 $TDSE_1D_HOME/post_processing/merge.py --printdata $printdata