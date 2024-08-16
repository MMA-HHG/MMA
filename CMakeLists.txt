### CMake file for the whole CUPRAD TDSE Hankel project

cmake_minimum_required(VERSION 3.13)
project(CUPRAD_TDSE_Hankel)

### If there are problems with compilers, export CC=mpicc/mpiicc and FC=mpifort/mpiifort

add_subdirectory(1DTDSE 1DTDSE/build)
add_subdirectory(CUPRAD CUPRAD/build)