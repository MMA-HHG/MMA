### CMake file for cuprad project
cmake_minimum_required(VERSION 3.5)

project(TDSE)
enable_language(C)      

message("Compiler: ${CMAKE_C_COMPILER}")
message("Compiler ID: ${CMAKE_C_COMPILER_ID}")

if(${CMAKE_C_COMPILER_ID} MATCHES Intel)
    set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
    set(MKL_LIB_DIRS ${MKLROOT}/lib)
    set(FFTW_LIBS mkl_gf_lp64 mkl_sequential mkl_core)
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIB_DIRS})
    #set(CMAKE_Fortran_FLAGS "-O0")
elseif(${CMAKE_C_COMPILER_ID} MATCHES GNU)
    find_library(fftw3 fftw3)
    #set(CMAKE_Fortran_FLAGS "-O0")
    set(FFTW_LIBS fftw3 fftw3_mpi)
elseif(${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    find_library(fftw3 fftw3)
    #set(CMAKE_Fortran_FLAGS "-O0")
    set(FFTW_LIBS fftw3 fftw3_mpi)
else()
    message(FATAL_ERROR "Unsupported C compiler: ${CMAKE_C_COMPILER}")
endif()

### Set compiler flags
set(CMAKE_C_FLAGS "-Wall -std=c99 -O3")

### Find HDF5 libraries and packages
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(MPI REQUIRED)

set(SOURCE sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c 
	sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c 
	sources/structures.c sources/tools.c sources/prop.c sources/singleTDSE.c 
	sources/strided_scheduler_separated_outs_coarsen.c)

set(SOURCE_STRIDE sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c 
	sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c 
	sources/structures.c sources/tools.c sources/prop.c sources/singleTDSE.c 
	sources/single_stride.c)

set(SOURCE_DLL sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c 
	sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c sources/structures.c 
	sources/tools.c sources/prop.c sources/singleTDSE.c)

set(TDSE TDSE.e)
set(TDSE_STRIDE TDSE_stride.e)
set(TDSE_LIB singleTDSE)

### Parallel TDSE
add_executable(${TDSE} ${SOURCE})
target_link_libraries(${TDSE} mpi ${FFTW_LIBS} ${HDF5_LIBRARIES})
target_include_directories(${TDSE} PRIVATE ${HDF5_INCLUDE_DIRS})

### Single TDSE
add_executable(${TDSE_STRIDE} ${SOURCE_STRIDE})
target_link_libraries(${TDSE_STRIDE} mpi ${FFTW_LIBS} ${HDF5_LIBRARIES})
target_include_directories(${TDSE_STRIDE} PRIVATE ${HDF5_INCLUDE_DIRS})

### Dynamic library
add_library(${TDSE_LIB} SHARED ${SOURCE_DLL})
target_link_libraries(${TDSE_LIB} PRIVATE mpi ${FFTW_LIBS} ${HDF5_LIBRARIES})
target_include_directories(${TDSE_LIB} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_options(${TDSE_LIB} PRIVATE -fPIC)