#!/bin/bash
#$ -cwd
#$ -pe smp 4
#$ -S /bin/bash

module purge 
module load hdf5-parallel-threadunsafe/1.14.1 fftw

mpirun -np 4 cuprad.e
