# CMake file for cuprad project
#cmake_minimum_required(VERSION 3.10)
cmake_minimum_required(VERSION 3.10)
project(cuprad)
enable_language(Fortran)      

# Set compiler flags
set(CMAKE_Fortran_FLAGS "-fbacktrace -ffpe-trap=zero -ffree-line-length-250")

# Set compiler
#set(MPI_Fortran_COMPILER "mpif90")

# Find FFTW and HDF5 libraries
#find_package(FFTW3 REQUIRED)
find_library(FFTW3 FFTW3 REQUIRED)
#find_library(HDF5 HDF5 REQUIRED COMPONENTS Fortran HL)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
include_directories(${HDF5_INCLUDE_DIRS})
#find_library(LIBRARY_HDF5 HDF5 REQUIRED)
#find_package(MPI REQUIRED)
#find_library(MPI REQUIRED)

# Add executable targets
add_executable(cuprad_occigen.e sources/constants.f90 sources/linked_list_module.f90 sources/hdf5_helper_serial.f90 sources/hdf5_helper.f90 sources/array_helper.f90 sources/modules.f90 sources/pre_ionised.f90 sources/code_ionisation.f90 sources/fft_rk.f90 sources/finalize.f90 sources/firststep.f90 sources/longstep_rk.f90 sources/output.f90 sources/cuprad.f90)
#target_link_libraries(cuprad_occigen.e PRIVATE LIBRARY_FFTW3 LIBRARY_HDF5)
#include_directories(sources)
#target_link_libraries(cuprad_occigen.e LIBRARY_FFTW3 HDF5)
target_link_libraries(cuprad_occigen.e FFTW3 ${HDF5_LIBRARIES})
#target_link_libraries(make_start_occigen.e PRIVATE LIBRARY_FFTW3::FFTW3 HDF5::HDF5_Fortran)

add_executable(make_start_occigen.e sources/constants.f90 sources/hdf5_helper_serial.f90 sources/hdf5_helper.f90 sources/write_start.f90 sources/calc_start.f90 sources/normalisation.f90 sources/write_listing.f90 sources/default_inputs.f90 sources/make_start.f90)
#target_link_libraries(make_start_occigen.e FFTW3 HDF5)
target_link_libraries(make_start_occigen.e FFTW3 ${HDF5_LIBRARIES})
#target_link_libraries(make_start_occigen.e PRIVATE LIBRARY_FFTW3 HDF5)

add_executable(clean_occigen.e sources/clean.f90)

add_executable(merge_occigen.e sources/merge.f90)

add_executable(calc_ft_occigen.e sources/constants.f90 sources/modules.f90 sources/fft.f90 sources/calc_ft.f90)
target_link_libraries(calc_ft_occigen.e FFTW3)
#target_link_libraries(calc_ft_occigen.e PRIVATE LIBRARY_FFTW3)

# Add custom commands
#add_custom_target(clean COMMAND rm -f *.e *.o *.mod)
add_custom_target(binary COMMAND if [ ! -d "../binary" ]; then mkdir ../binary; fi; mv *.e ../binary/)
add_custom_target(clean_preprocessor COMMAND rm -f listing STOP *.DAT *.LOG)
add_custom_target(code DEPENDS clean_preprocessor sources/linked_list_module.o sources/calc_ft_occigen.e sources/make_start_occigen.e sources/cuprad_occigen.e sources/clean_occigen.e sources/merge_occigen.e binary clean)
add_custom_target(preprocessor DEPENDS clean clean_preprocessor make_start_occigen.e)
add_custom_target(util DEPENDS calc_ft_occigen.e make_start_occigen.e clean_occigen.e merge_occigen.e binary clean)

# Source: Conversation with Bing, 25. 7. 2023
# (1) cmake - Convert Makefile into CMakeLists, where to start - Stack Overflow. https://stackoverflow.com/questions/40860769/convert-makefile-into-cmakelists-where-to-start.
# (2) Converting Existing Systems To CMake — Mastering CMake. https://cmake.org/cmake/help/book/mastering-cmake/chapter/Converting%20Existing%20Systems%20To%20CMake.html.
# (3) c++ - Converting a Makefile project to Cmake - Stack Overflow. https://stackoverflow.com/questions/18680977/converting-a-makefile-project-to-cmake.