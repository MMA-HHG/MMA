#!/bin/bash
h5filename=results_*.h5
var=$(echo $h5filename)

module purge
module load intel intelmpi hdf5

/home/vabekjan/git/CUPRAD_DEVELOP/binary/make_start_occigen.e <<INPUTS
$var
0
0
0
INPUTS