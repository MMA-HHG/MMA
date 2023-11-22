#!/bin/bash

### Load cmake, MPI, HDF5, FFTW/MKL
module load cmake/3.26.3/gcc@11.2.0-rthuom3
module load compiler/intel/2020.4.304
module load mpi/intel/2020.4.304
module load mkl/2023.0.0
module load hdf5/1.10.5/openmpi_4.1.2/intel_2022.0.1
module load hdf5/1.10.5/impi_2019.1.144/intel_2019.1.144