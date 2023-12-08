#!/bin/bash

### Usage information
usage()
{
    echo "A script for running the whole multiscale model.
----------------------------------------------------
usage: [[[-i (--inp) input file with CUPRAD parameters (.inp)] 
         [-o (--ohdf5) output HDF5 file] 
         [-n (--ncuprad) number of processes for CUPRAD execution (must fit input file!)] 
         [-m (--ntdse) number of processes for MPI-TDSE execution]  
         [-p (--printdata) data to print into the merged HDF5 file]] | 
         [-h help]]

----------------------------------------------------
In order to run the script properly, environment variables must be set using 'set_env_vars.sh' 
script in the current terminal shell as we execute the 'run_multiscale.sh' script.

Example run: source /path/to/script/../run_multiscale.sh -i [input_file.inp] -o [result.h5] -n [N_tasks_cuprad] -m [M_tasks_tdse] -p ['Efield FSourceTerm SourceTerm expval_x PopInt PopTot ...']


Print data available options: 'Efield', 'FEfield', 'SourceTerm', 'FSourceTerm', 
                              'FEfieldM2', 'FSourceTermM2', 'PopTot', 'PopInt', 
                              'expval_x'
Note: to use this option, the keywords must be in parentheses, e.g. "Efield SourceTerm".
"
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
        -n | --ncuprad)         shift
                                ntasks_cuprad="$1"
                                ;;
        -m | --ntdse)           shift
                                ntasks_tdse="$1"
                                ;;
        -p | --printdata)       shift
                                printdata="$1"
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit
    esac
    shift
done

if [ "$inp_filename" = "" ]
then
    echo "No input file (.inp) added. Use option [-i (--inp)]."
    exit
fi

if [ "$h5_filename" = "" ]
then
    echo "No output HDF5 file added. Use option [-o (--ohdf5)]."
    exit
fi

if [ "$ntasks_cuprad" = "" ]
then
    echo "Unspecified number of processes for CUPRAD. Use option [-n (--ncuprad)]."
    exit
fi

if [ "$ntasks_tdse" = "" ]
then
    echo "Unspecified number of processes for MPI-TDSE. Use option [-m (--ntdse)]."
    exit
fi

echo "Input file: ${inp_filename}"
echo "Output hdf5 file: ${h5_filename}"
echo "Number of processes for CUPRAD: ${ntasks_cuprad}"
echo "Number of processes for MPI-TDSE: ${ntasks_tdse}"

### ENVIRONMENT VARIABLES MUST BE SET!
#source set_env_vars.sh

### Create HDF5 file from parameters
JOB0=$(sbatch --parsable --export=ALL $CUPRAD_SCRIPTS/make_hdf5.sh \
        --inp $inp_filename --ohdf5 $h5_filename -s);

### Submit Preprocessor
JOB1=$(sbatch --parsable --export=ALL --dependency=afterok:$JOB0 \
        $CUPRAD_SCRIPTS/run_preprocessor.sh --ihdf5 $h5_filename -s);

### Submit CUPRAD
JOB2=$(sbatch --ntasks=$ntasks_cuprad --parsable --export=ALL --dependency=afterok:$JOB1 \
        $CUPRAD_SCRIPTS/CUPRAD.sh)

### Submit TDSE
JOB3=$(sbatch --ntasks=$ntasks_tdse --parsable --export=ALL --dependency=afterok:$JOB2 \
        $TDSE_1D_SCRIPTS/MPI_TDSE.sh)

### Merge TDSEs
JOB4=$(sbatch --parsable --export=ALL --dependency=afterok:$JOB3 \
        $TDSE_1D_SCRIPTS/merge_hdf5.sh -p $printdata -s)
