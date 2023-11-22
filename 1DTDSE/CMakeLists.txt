### CMake file for the CTDSE
cmake_minimum_required(VERSION 3.13)

project(TDSE)
enable_language(C)      

message("Compiler: ${CMAKE_C_COMPILER}")
message("Compiler ID: ${CMAKE_C_COMPILER_ID}")

if(${CMAKE_C_COMPILER_ID} MATCHES Intel)
    set(MKL_INCLUDE_DIRS ${MKLROOT}/include/fftw)
    set(MKL_LIB_DIRS ${MKLROOT}/lib)
    set(FFTW_LIBS mkl_gf_lp64 mkl_sequential mkl_core)
elseif(${CMAKE_C_COMPILER_ID} MATCHES GNU)
    find_library(fftw3 fftw3)
    set(FFTW_LIBS fftw3) 
elseif(${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    find_library(fftw3 fftw3)
    set(FFTW_LIBS fftw3 fftw3_mpi)
else()
    message(FATAL_ERROR "Unsupported C compiler: ${CMAKE_C_COMPILER}")
endif()

### Find HDF5 libraries and packages
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(MPI)


### Dynamic library & single TDSE
set(SOURCE_STRIDE sources/constants.c sources/tools_algorithmic.c
	sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c 
	sources/structures.c sources/tools.c sources/prop.c sources/singleTDSE.c 
	sources/single_stride.c)

set(SOURCE_DLL sources/constants.c sources/tools_algorithmic.c 
	sources/tridiag.c sources/tools_fftw3.c sources/structures.c 
	sources/tools.c sources/prop.c sources/singleTDSE.c)

set(TDSE_STRIDE TDSE_stride.e)
set(TDSE_LIB singleTDSE)

### Single TDSE
add_executable(${TDSE_STRIDE} ${SOURCE_STRIDE})
target_link_libraries(${TDSE_STRIDE} -lm -ldl ${FFTW_LIBS} ${HDF5_LIBRARIES})
target_include_directories(${TDSE_STRIDE} PRIVATE ${HDF5_INCLUDE_DIRS})

### Dynamic library
add_library(${TDSE_LIB} SHARED ${SOURCE_DLL})
target_link_libraries(${TDSE_LIB} PRIVATE ${FFTW_LIBS} ${HDF5_LIBRARIES})
target_include_directories(${TDSE_LIB} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_options(${TDSE_LIB} PRIVATE -fPIC)


### MPI-TDSE scheduler (only if MPI is available)
if(MPI_FOUND)

    set(SOURCE sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c 
        sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c 
        sources/structures.c sources/tools.c sources/prop.c sources/singleTDSE.c 
        sources/MPI-TDSE.c)

    set(TDSE TDSE.e)

    ### Parallel TDSE
    add_executable(${TDSE} ${SOURCE})
    target_link_libraries(${TDSE} mpi ${FFTW_LIBS} ${HDF5_LIBRARIES})
    target_include_directories(${TDSE} PRIVATE ${HDF5_INCLUDE_DIRS})

else()
   message(WARNING "MPI not found. Skipping MPI-scheduler won't be built.")

endif()