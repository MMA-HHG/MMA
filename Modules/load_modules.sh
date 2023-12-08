#!/bin/bash

### Define your modules and HPCs here
### Define separately modules for compiling and modules for Python
Curta_modules() {
    module load cmake/3.26.3/gcc@11.2.0-rthuom3
    module load compiler/intel/2020.4.304
    module load mpi/intel/2020.4.304
    module load mkl/2023.0.0
    module load hdf5/1.10.5/openmpi_4.1.2/intel_2022.0.1
    module load hdf5/1.10.5/impi_2019.1.144/intel_2019.1.144
}
export -f Curta_modules

Curta_python_modules() {
    module load python/3.9
}
export -f Curta_python_modules

### Load cmake, MPI, HDF5, FFTW/MKL
load_modules() {
    if [ "$HPC" == "Curta" ] 
    then
        Curta_modules
    fi

}
export -f load_modules

load_python_modules() {
    if ["$HPC" = "Curta"] 
    then
        Curta_python_modules
    fi

}
export -f load_python_modules

echo "Functions loaded into env."
echo "HPC = $HPC"
