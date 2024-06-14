#!/bin/bash
h5filename=results.h5
var=$(echo $h5filename)

load_modules

$CUPRAD_HOME/build/make_start.e <<INPUTS
$var
0
0
0
INPUTS