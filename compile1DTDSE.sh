#!/bin/bash

gcc -O2 -c rwhdf5.c
gcc *.o -lm -o rwhdf5.e
rm *.o