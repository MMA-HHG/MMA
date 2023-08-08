### CMake file for cuprad project
cmake_minimum_required(VERSION 3.27)

#cmake_policy(SET CMP0074 NEW)


project(cuprad)
enable_language(Fortran)      

#list (APPEND CMAKE_PREFIX_PATH "/usr/local/Cellar/gcc/13.1.0/lib/gcc/current")
#set(CMAKE_Fortran_MODULE_DIRECTORY ".")


### Set compiler flags
#if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(dialect "-ffree-form -std=f2008 -fimplicit-none -fcheck=all -W 
#        -Wno-compare-reals -fbacktrace -ffixed-line-length-none -ffree-line-length-none -g")
#    set(bounds "-fbounds-check")
#endif()
#if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    set(dialect "-stand f08 -free -implicitnone")
#    set(bounds "-check bounds")
#endif()

#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_FLAGS "-std=gnu -O0 -fcheck=all -W -Wno-compare-reals -fbacktrace -ffixed-line-length-none -ffree-line-length-none -g")

### Set compiler
#SET(CMAKE_Fortran_COMPILER gfortran)

### Find FFTW and HDF5 libraries and packages
find_library(fftw3 fftw3)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
#find_package(MPI REQUIRED Fortran)
find_package(MPI REQUIRED Fortran)

### Link FFTW3 library
link_directories(/usr/local/lib)

### Add MPI
#include_directories(${MPI_Fortran_INCLUDE_PATH})
#add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

### Define variables
### ================
### Main propagation
set(PROPAGATION cuprad.e)
set(S_PROPAGATION sources/constants.f90 sources/linked_list_module.f90 
    sources/hdf5_helper_serial.f90 sources/hdf5_helper.f90 sources/array_helper.f90 
    sources/modules.f90 sources/pre_ionised.f90 sources/code_ionisation.f90 
    sources/fft_rk.f90 sources/finalize.f90 sources/firststep.f90 sources/longstep_rk.f90 
    sources/output.f90 sources/cuprad.f90)

set(MAKE_START make_start.e)
set(S_MAKE_START sources/constants.f90 sources/hdf5_helper_serial.f90 
    sources/hdf5_helper.f90 sources/write_start.f90 sources/calc_start.f90 
    sources/normalisation.f90 sources/write_listing.f90 sources/default_inputs.f90 
    sources/make_start.f90)

set(CLEAN clean.e)
set(S_CLEAN sources/clean.f90)

set(MERGE merge.e)
set(S_MERGE sources/merge.f90)

set(CALC_FT calc_ft.e)
set(S_CALC_FT sources/constants.f90 sources/linked_list_module.f90 
    sources/modules.f90 sources/fft.f90 sources/calc_ft.f90)

### Executable targets
### ==================
### cuprad.e
add_executable(${PROPAGATION} ${S_PROPAGATION})
target_link_libraries(${PROPAGATION} fftw3 fftw3_mpi ${HDF5_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_include_directories(${PROPAGATION} PRIVATE ${HDF5_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_PATH})

### make_start.e
add_executable(${MAKE_START} ${S_MAKE_START})
target_include_directories(${MAKE_START} PRIVATE ${HDF5_INCLUDE_DIRS})# /usr/local/include)
target_link_libraries(${MAKE_START} ${HDF5_LIBRARIES} fftw3 fftw3_mpi)# ${MPI_Fortran_LIBRARIES})

### clean.e
add_executable(${CLEAN} ${S_CLEAN})

### merge.e
add_executable(${MERGE} ${S_MERGE})

### calc_ft.e
add_executable(${CALC_FT} ${S_CALC_FT})
target_link_libraries(${CALC_FT} fftw3 fftw3_mpi ${MPI_Fortran_LIBRARIES})
target_include_directories(${CALC_FT} PRIVATE ${MPI_Fortran_INCLUDE_PATH})

### Linked list module
#add_executable(linked_list.mod sources/linked_list_module.f90)


### Add custom commands
#add_custom_target(mv_binary 
#    COMMAND if [ ! -d "../binary" ]; then mkdir ../binary; fi
#    COMMAND /bin/mv -f *.e ../binary/
#)
#add_custom_target(linked_list DEPENDS linked_list.mod)
add_custom_target(clean_preprocessor COMMAND rm -f listing STOP *.DAT *.LOG)
add_custom_target(code DEPENDS clean_preprocessor ${CALC_FT} ${MAKE_START} ${PROPAGATION} ${CLEAN} ${MERGE})
add_custom_target(preprocessor DEPENDS clean_preprocessor ${MAKE_START})
add_custom_target(util DEPENDS ${CALC_FT} ${MAKE_START} ${CLEAN} ${MERGE})

### Source: Conversation with Bing, 25. 7. 2023
### (1) cmake - Convert Makefile into CMakeLists, where to start - Stack Overflow. https://stackoverflow.com/questions/40860769/convert-makefile-into-cmakelists-where-to-start.
### (2) Converting Existing Systems To CMake — Mastering CMake. https://cmake.org/cmake/help/book/mastering-cmake/chapter/Converting%20Existing%20Systems%20To%20CMake.html.
### (3) c++ - Converting a Makefile project to Cmake - Stack Overflow. https://stackoverflow.com/questions/18680977/converting-a-makefile-project-to-cmake.