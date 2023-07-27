### CMake file for cuprad project
cmake_minimum_required(VERSION 3.11)

project(cuprad)
enable_language(Fortran)      

### Set compiler flags
set(CMAKE_Fortran_FLAGS "-fcheck=all -W -Wno-compare-reals -fbacktrace -ffpe-trap=zero -ffree-line-length-250 -g")

### Find FFTW and HDF5 libraries and packages
find_library(fftw3 fftw3)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)

### Link FFTW3 library
link_directories(/usr/local/lib)

### Define variables
### ================
### Main propagation
set(PROPAGATION cuprad_occigen.e)
set(S_PROPAGATION sources/constants.f90 sources/linked_list_module.f90 
    sources/hdf5_helper_serial.f90 sources/hdf5_helper.f90 sources/array_helper.f90 
    sources/modules.f90 sources/pre_ionised.f90 sources/code_ionisation.f90 
    sources/fft_rk.f90 sources/finalize.f90 sources/firststep.f90 sources/longstep_rk.f90 
    sources/output.f90 sources/cuprad.f90)

set(MAKE_START make_start_occigen.e)
set(S_MAKE_START sources/constants.f90 sources/hdf5_helper_serial.f90 
    sources/hdf5_helper.f90 sources/write_start.f90 sources/calc_start.f90 
    sources/normalisation.f90 sources/write_listing.f90 sources/default_inputs.f90 sources/make_start.f90)

set(CLEAN clean_occigen.e)
set(S_CLEAN sources/clean.f90)

set(MERGE merge_occigen.e)
set(S_MERGE sources/merge.f90)

set(CALC_FT calc_ft_occigen.e)
set(S_CALC_FT sources/constants.f90 sources/modules.f90 sources/fft.f90 sources/calc_ft.f90)

### Executable targets
### ==================
### cuprad_occigen.e
add_executable(${PROPAGATION} ${S_PROPAGATION})
target_link_libraries(${PROPAGATION} fftw3 fftw3_mpi ${HDF5_LIBRARIES})
target_include_directories(${PROPAGATION} PRIVATE ${HDF5_INCLUDE_DIRS} /usr/local/include)

### make_start_occigen.e
add_executable(${MAKE_START} ${S_MAKE_START})
target_include_directories(${MAKE_START} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${MAKE_START} ${HDF5_LIBRARIES} fftw3 fftw3_mpi)

### clean_occigen.e
add_executable(${CLEAN} ${S_CLEAN})

### merge_occigen.e
add_executable(${MERGE} ${S_MERGE})

### calc_ft_occigen.e
add_executable(${CALC_FT} ${S_CALC_FT})
target_link_libraries(${CALC_FT} fftw3 fftw3_mpi)


### Add custom commands
#add_custom_target(mv_binary 
#    COMMAND if [ ! -d "../binary" ]; then mkdir ../binary; fi
#    COMMAND /bin/mv -f *.e ../binary/
#)
add_custom_target(clean_preprocessor COMMAND rm -f listing STOP *.DAT *.LOG)
add_custom_target(code DEPENDS clean_preprocessor ${CALC_FT} ${MAKE_START} ${PROPAGATION} ${CLEAN} ${MERGE})
add_custom_target(preprocessor DEPENDS clean_preprocessor ${MAKE_START})
add_custom_target(util DEPENDS ${CALC_FT} ${MAKE_START} ${CLEAN} ${MERGE})

### Source: Conversation with Bing, 25. 7. 2023
### (1) cmake - Convert Makefile into CMakeLists, where to start - Stack Overflow. https://stackoverflow.com/questions/40860769/convert-makefile-into-cmakelists-where-to-start.
### (2) Converting Existing Systems To CMake — Mastering CMake. https://cmake.org/cmake/help/book/mastering-cmake/chapter/Converting%20Existing%20Systems%20To%20CMake.html.
### (3) c++ - Converting a Makefile project to Cmake - Stack Overflow. https://stackoverflow.com/questions/18680977/converting-a-makefile-project-to-cmake.