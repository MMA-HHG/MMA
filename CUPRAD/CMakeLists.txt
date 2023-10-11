### CMake file for cuprad project
cmake_minimum_required(VERSION 3.5)

project(cuprad)
enable_language(Fortran)      

message("Compiler: ${CMAKE_Fortran_COMPILER}")
message("Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")


if(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
    set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
    set(MKL_LIB_DIRS ${MKLROOT}/lib)
    set(MKL_LIBS mkl_gf_lp64 mkl_sequential mkl_core)
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIB_DIRS})
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
    find_library(fftw3 fftw3)
    set(MKL_LIBS fftw3 fftw3_mpi)
else()
    message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER}")
endif()
#execute_process(COMMAND which mpifort OUTPUT_VARIABLE MPIFORT)
#execute_process(COMMAND which mpiifort OUTPUT_VARIABLE MPIIFORT)
#execute_process(COMMAND MPIIFORT=${MPIIFORT} | tr -d '\n')
#execute_process(COMMAND MPIFORT=${MPIFORT} | tr -d '\n')
#message(${MPIFORT})
#message(${MPIIFORT})
#execute_process(COMMAND echo ${MPIIFORT} | tr -d '\n')

#if(${CMAKE_Fortran_COMPILER} STREQUAL ${MPIFORT})
    ### Set compiler flags
#    set(CMAKE_Fortran_FLAGS "-std=gnu -O0 -fcheck=all -W -Wno-compare-reals -fbacktrace -ffixed-line-length-none -ffree-line-length-none -g")

    ### Find FFTW and HDF5 libraries and packages
#    find_library(fftw3 fftw3)
#    find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
#    find_package(MPI REQUIRED Fortran)
#elseif(${CMAKE_Fortran_COMPILER} STREQUAL ${MPIIFORT})

#find_package(MKL CONFIG REQUIRED)
#include_directories(${MKL_INCLUDE_DIRS})
#link_directories(${MKL_LIB_DIRS})
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
find_package(MPI REQUIRED Fortran)

#set(CMAKE_Fortran_FLAGS "-O0")

#else()
#    message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER}")
#endif()

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
#target_link_libraries(${PROPAGATION} fftw3 fftw3_mpi ${HDF5_LIBRARIES} ${MPI_Fortran_LIBRARIES})
#target_compile_options(${PROPAGATION} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
#target_link_libraries(${PROPAGATION} PUBLIC $<LINK_ONLY:MKL::MKL> ${HDF5_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_link_libraries(${PROPAGATION} ${HDF5_LIBRARIES} ${MKL_LIBS})#fftw3 fftw3_mpi)
target_include_directories(${PROPAGATION} PRIVATE ${HDF5_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS})
#target_include_directories(${PROPAGATION} PRIVATE ${HDF5_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_PATH} $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)# /usr/local/include)
#target_compile_options(${PROPAGATION} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)


### make_start.e
add_executable(${MAKE_START} ${S_MAKE_START})
target_include_directories(${MAKE_START} PRIVATE ${HDF5_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS})# /usr/local/include)
target_link_libraries(${MAKE_START} ${HDF5_LIBRARIES} ${MKL_LIBS}) # fftw3 fftw3_mpi)# ${MPI_Fortran_LIBRARIES})
#target_compile_options(${MAKE_START} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
#target_link_libraries(${MAKE_START} PUBLIC $<LINK_ONLY:MKL::MKL> ${HDF5_LIBRARIES})
#target_include_directories(${MAKE_START} PRIVATE ${HDF5_INCLUDE_DIRS} $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)# /usr/local/include)
#target_compile_options(${MAKE_START} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)

### Add custom commands
### ===================
add_custom_target(clean_preprocessor COMMAND rm -f listing STOP *.DAT *.LOG)
add_custom_target(clean_modules COMMAND rm -f *.mod)
add_custom_target(code DEPENDS clean_preprocessor ${MAKE_START} ${PROPAGATION})
add_custom_target(preprocessor DEPENDS clean_preprocessor ${MAKE_START})
