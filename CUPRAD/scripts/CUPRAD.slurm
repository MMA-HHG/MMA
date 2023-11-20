#!/bin/bash
#SBATCH --job-name="CUPRAD"
#SBATCH -t 0-01:00 # time (D-HH:MM)
#SBATCH -o CUPRAD.%j.%N.out # STDOUT
#SBATCH -e CUPRAD.%j.%N.err # STDERR
#SBATCH --profile=All

### Purge modules
module purge

### Load modules
module load cmake/3.26.3/gcc@11.2.0-rthuom3
module load compiler/intel/2020.4.304
module load mpi/intel/2020.4.304
module load mkl/2023.0.0
module load hdf5/1.10.5/openmpi_4.1.2/intel_2022.0.1
module load hdf5/1.10.5/impi_2019.1.144/intel_2019.1.144

mpirun $CUPRAD_HOME/build/cuprad.e