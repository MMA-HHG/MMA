#!/bin/bash

### Define your modules and HPCs here
### Define separately modules for compiling and modules for Python
### We use 'module purge' policy for loading the module to prevent conflict with already modified environment.
Curta_modules() {
    module purge
    module load cmake/3.26.3/gcc@11.2.0-qh75jnq
    module load compiler/intel/2020.4.304
    module load mpi/intel/2020.4.304
    module load mkl/2023.0.0
    module load hdf5/1.10.5/openmpi_4.1.2/intel_2022.0.1
    module load hdf5/1.10.5/impi_2019.1.144/intel_2019.1.144
    export CPATH=$CPATH:$MKLROOT/include/fftw
}
export -f Curta_modules

Sunrise_modules() {
    module purge
    module load GCC OpenMPI HDF5 FFTW3 CMake 
}
export -f Sunrise_modules

JeanZay_intel_modules() {
    module purge
    module load hdf5 
    module load intel-mkl
    module load cmake
}
export -f JeanZay_intel_modules

Metacentrum_modules() {
    module purge
    module add fftw/3.3.10-gcc-10.2.1-yxsjm6z hdf5/1.12.2-gcc-10.2.1-gfdwqr3 
    module add openblas/0.3.20-gcc-10.2.1-p4skjks
    module load cmake
}
export -f Metacentrum_modules


### Python modules
Curta_python_modules() {
    module purge
    module load python/3.9
}
export -f Curta_python_modules

Sunrise_python_modules() {
    module purge
    module load GCC Python
}
export -f Sunrise_python_modules


JeanZay_python_modules() {
    module purge
    module load python/3.9.12
}
export -f JeanZay_python_modules

### Load cmake, MPI, HDF5, FFTW/MKL
load_modules() {
    if [ "$HPC" == "Curta" ] 
    then
        Curta_modules
    elif [ "$HPC" == "Sunrise" ]
    then
        Sunrise_modules
    elif [ "$HPC" == "JeanZay-intel" ]
    then
        JeanZay_intel_modules
    elif [ "$HPC" == "Metacentrum" ]
    then
        Metacentrum_modules
    fi

}
export -f load_modules

load_python_modules() {
    if [ "$HPC" == "Curta" ] 
    then
        Curta_python_modules
    elif [ "$HPC" == "Sunrise" ]
    then
        Sunrise_python_modules
    elif [ "$HPC" == "JeanZay-intel" ]
    then
        JeanZay_python_modules
    fi

}
export -f load_python_modules

echo "Functions loaded into env."
echo "HPC = $HPC"
