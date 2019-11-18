#!/bin/bash

gcc -O2 -c $HOME/git/1DTDSE/*.c -I/gpfs/softs/contrib/apps/fftw3/3.3.8/include/ -L/gpfs/softs/contrib/apps/fftw3/3.3.8/lib/D
#gcc *.o -lfftw3 -lm -o TDSE1D.out
gcc *.o -I/gpfs/softs/contrib/apps/fftw3/3.3.8/include/ -L/gpfs/softs/contrib/apps/fftw3/3.3.8/lib/ -lfftw3 -lm -o TDSE1D.out
rm *.o

