#!/bin/bash

##################
# COMPILE SCRIPT #
##################

# This script compiles all the available codes and prepares the directories.

# READ THE INSTRUCTIONS in README.md on how to run the script!

# Important prerequisites: ...
#   - properly loaded all the important modules/set environment variables

CMAKE_Fort_Comp=mpifort
CMAKE_C_Comp=mpicc

# Creating the build directories for CMake
mkdir CUPRAD/build
mkdir 1DTDSE/build

# CUPRAD build
cd CUPRAD/build 
cmake -D CMAKE_Fortran_COMPILER=${CMAKE_Fort_Comp} ..

# TDSE build
cd CUPRAD/build 
cmake -D CMAKE_Fortran_COMPILER=${CMAKE_C_Comp} ..


