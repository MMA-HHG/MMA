#!/bin/sh

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
    module load cmake/3.26.3/gcc@11.2.0-rthuom3
    module load compiler/intel/2020.4.304
    module load mpi/intel/2020.4.304
    module load mkl/2023.0.0
    module load hdf5/1.10.5/openmpi_4.1.2/intel_2022.0.1
    module load hdf5/1.10.5/impi_2019.1.144/intel_2019.1.144
}

### Main
h5_filename=""

while [ "$1" != "" ]; do
    case $1 in
        -i | --ihdf5 )           shift
                                h5_filename="$1"
                                ;;
        -s | --slurm )          slurm
                                exit
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
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
