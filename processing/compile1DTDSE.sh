#!/bin/bash

icc -O2 -c $HOME/git/1DTDSE/rwhdf5.c
icc *.o -lfftw3 -lm -o rwhdf5.out
rm *.o