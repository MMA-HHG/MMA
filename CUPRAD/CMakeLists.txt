### CMake file for cuprad project
cmake_minimum_required(VERSION 3.5)

#cmake_policy(SET CMP0074 NEW)

project(cuprad)
enable_language(Fortran)      

### Set compiler flags
set(CMAKE_Fortran_FLAGS "-std=gnu -O0 -fcheck=all -W -Wno-compare-reals -fbacktrace -ffixed-line-length-none -ffree-line-length-none -g")

### Find FFTW and HDF5 libraries and packages
find_library(fftw3 fftw3)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
find_package(MPI REQUIRED Fortran)

### Link FFTW3 library
#link_directories(/usr/local/lib)

### Define variables
### ================
### Main propagation
set(PROPAGATION cuprad.e)
set(S_PROPAGATION sources/constants.f90 sources/linked_list_module.f90 
    sources/hdf5_helper_serial.f90 sources/hdf5_helper.f90 sources/array_helper.f90 
    sources/modules.f90 sources/pre_ionised.f90 sources/code_ionisation.f90 
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
target_link_libraries(${PROPAGATION} fftw3 fftw3_mpi ${HDF5_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_include_directories(${PROPAGATION} PRIVATE ${HDF5_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_PATH})# /usr/local/include)

### make_start.e
add_executable(${MAKE_START} ${S_MAKE_START})
target_include_directories(${MAKE_START} PRIVATE ${HDF5_INCLUDE_DIRS})# /usr/local/include)
target_link_libraries(${MAKE_START} ${HDF5_LIBRARIES} fftw3 fftw3_mpi)# ${MPI_Fortran_LIBRARIES})


### Add custom commands
### ===================
add_custom_target(clean_preprocessor COMMAND rm -f listing STOP *.DAT *.LOG)
add_custom_target(clean_modules COMMAND rm -f *.mod)
add_custom_target(code DEPENDS clean_preprocessor ${MAKE_START} ${PROPAGATION})
add_custom_target(preprocessor DEPENDS clean_preprocessor ${MAKE_START})
