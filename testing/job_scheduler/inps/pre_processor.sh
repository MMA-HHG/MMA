#!/bin/bash


echo 'test environment'
env

h5filename=results_*.h5
# echo $h5filename

module purge
module load intel intelmpi hdf5

# home/vabekjan/git/CUPRAD_DEVELOP/binary/make_start_occigen.e << INPUTS
# $h5filename
# 0
# 0
# 0
# INPUTS


echo "pre-processor" 
