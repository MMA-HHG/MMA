### CMake file for cuprad project
cmake_minimum_required(VERSION 3.5)

project(cuprad)
enable_language(Fortran)      

message("Compiler: ${CMAKE_Fortran_COMPILER}")
message("Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")

if(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
    set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
    set(MKL_LIB_DIRS ${MKLROOT}/lib)
    set(FFTW_LIBS mkl_gf_lp64 mkl_sequential mkl_core)
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIB_DIRS})
#    set(CMAKE_Fortran_FLAGS "-O3")
    set(CMAKE_Fortran_FLAGS "-check all -warn all -traceback -fpe0 -fp-stack-check -diag-disable=8889") # disabel warning caused by the direct uses of mkl- & mpi-lib procedures
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
    find_library(fftw3 fftw3)
    set(CMAKE_Fortran_FLAGS "-O3 -Wall -std=gnu -fcheck=all -Wno-compare-reals -fbacktrace -ffixed-line-length-none -ffree-line-length-none -g")
    set(FFTW_LIBS fftw3 fftw3_mpi)
else()
    message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER}")
endif()

find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
find_package(MPI REQUIRED Fortran)

### Define variables
### ================
### Main propagation
set(PROPAGATION cuprad.e)
set(S_PROPAGATION sources/constants.f90 sources/linked_list_module.f90 
    sources/hdf5_helper_serial.f90 sources/hdf5_helper.f90 sources/array_helper.f90 
    sources/modules.f90 sources/pre_ionised.f90 sources/density_mod.f90 sources/code_ionisation.f90 
    sources/fft_rk.f90 sources/finalize.f90 sources/firststep.f90 sources/longstep_rk.f90 
    sources/output.f90 sources/cuprad.f90)

### Preprocessor - preparation of 'results.h5' file
set(MAKE_START make_start.e)
set(S_MAKE_START sources/constants.f90 sources/hdf5_helper_serial.f90 
    sources/hdf5_helper.f90 sources/write_start.f90 sources/calc_start.f90 
    sources/normalisation.f90 sources/write_listing.f90 sources/default_inputs.f90 
    sources/make_start.f90)

### Executable targets
### ==================
### cuprad.e
add_executable(${PROPAGATION} ${S_PROPAGATION})
target_link_libraries(${PROPAGATION} ${HDF5_LIBRARIES} ${FFTW_LIBS})
target_include_directories(${PROPAGATION} PRIVATE ${HDF5_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS})

### make_start.e
add_executable(${MAKE_START} ${S_MAKE_START})
target_include_directories(${MAKE_START} PRIVATE ${HDF5_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS})
target_link_libraries(${MAKE_START} ${HDF5_LIBRARIES} ${FFTW_LIBS})

### Add custom commands
### ===================
add_custom_target(clean_preprocessor COMMAND rm -f listing STOP *.DAT *.LOG)
add_custom_target(clean_modules COMMAND rm -f *.mod)
add_custom_target(code DEPENDS clean_preprocessor ${MAKE_START} ${PROPAGATION})
add_custom_target(preprocessor DEPENDS clean_preprocessor ${MAKE_START})
