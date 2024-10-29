#!/bin/bash

# var=$(echo $H5FILE)
var=$(echo *.h5)


echo $h5filename

load_modules

$CUPRAD_HOME/build/make_start.e <<INPUTS
$var
0
0
0
INPUTS