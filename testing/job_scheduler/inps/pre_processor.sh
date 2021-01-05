#!/bin/bash


# echo 'test environment'
# env


# echo 'test environment'
h5filename=results_*.h5
echo $h5filename

var=$(echo $h5filename)

echo $var

module purge
module load intel intelmpi hdf5

/home/vabekjan/git/CUPRAD_DEVELOP/binary/make_start_occigen.e <<INPUTS
$var
0
0
0
INPUTS


echo "pre-processor" 
