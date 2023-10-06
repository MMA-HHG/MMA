### CMake file for cuprad project
cmake_minimum_required(VERSION 3.5)

#cmake_policy(SET CMP0074 NEW)

project(TDSE)
enable_language(C)      

### Set compiler flags
set(CMAKE_C_FLAGS "-Wall -std=c99 -O3")

### Find FFTW and HDF5 libraries and packages
find_library(fftw3 fftw3)
find_package(HDF5 REQUIRED COMPONENTS C HL)

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
target_link_libraries(${TDSE} fftw3 hdf5 ${HDF5_LIBRARIES})
target_include_directories(${TDSE} PRIVATE ${HDF5_INCLUDE_DIRS})

### Single TDSE
add_executable(${TDSE_STRIDE} ${SOURCE_STRIDE})
target_link_libraries(${TDSE_STRIDE} fftw3 hdf5 ${HDF5_LIBRARIES})
target_include_directories(${TDSE_STRIDE} PRIVATE ${HDF5_INCLUDE_DIRS})

### Dynamic library
add_library(${TDSE_LIB} SHARED ${SOURCE_DLL})
target_link_libraries(${TDSE_LIB} PRIVATE fftw3 hdf5 ${HDF5_LIBRARIES})
target_include_directories(${TDSE_LIB} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_options(${TDSE_LIB} PRIVATE -fPIC)