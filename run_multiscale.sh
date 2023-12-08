#!/bin/bash


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

if [ "$ntasks_cuprad" = "" ]
then
    echo "Unspecified number of processes for CUPRAD. Use option [-n (--ncuprad)]."
    exit 1
fi

echo "Input file: ${inp_filename}"
echo "Output hdf5 file: ${h5_filename}"
echo "Number of processes for CUPRAD: ${ntasks_cuprad}"

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
