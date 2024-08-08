#!/bin/bash

var=$(echo $H5FILE)


echo $h5filename

load_modules

$CUPRAD_HOME/build/make_start.e <<INPUTS
$var
0
0
0
INPUTS