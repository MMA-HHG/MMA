#!/bin/bash

var=$(echo *.h5)
echo $h5filename

load_modules

$CUPRAD_HOME/build/make_start.e <<INPUTS
$var
INPUTS