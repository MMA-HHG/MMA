General remarks & problems occured along the way
================================================

* PREREQUISITE: MPI version of HDF5 has to be installed (should be on clusters anyway)
    <--- Fortran HDF5 is inherently not thread-safe!!!

* Compilation flags?? (which to select, i.e. -Wall etc.)

* Purpose of 'MKL_LIBS_MT' variable?
    <--- Likely a local variable on Occigen cluster

* CMAKE needs to have packages already installed – may be a drawback!!!

* In ```CMakeLists.txt``` the line ```include_directories(/usr/local/include)```
    has to be there else the system grabs incorrect version of hdf5.mod file 
    incompatible with the compiling language (older version of gfortran vs 
    mpif90 wrapper). This line had to be essentially hardcoded to compile with HDF5
    libraries. 
    <--- Fixed by correctly loading the HDF5 package and including via HDF5_INCLUDE_DIRS

* HDF5_ROOT has to be set to the proper installation of the HDF5. My local previous
    installation was set to conda installation, which didn't include the HDF5 with 
    MPI. The workaround was to set the environment variable as follows:
        ```export HDF5_ROOT=/usr/local/```
    Then it was able to find the proper installation and compile. I just had to 
    set:
        ```find_package(HDF5 REQUIRED COMPONENTS Fortran HL)```
    and then for the target library linking:
        ```target_link_libraries(cuprad_occigen.e fftw3 fftw3_mpi ${HDF5_LIBRARIES})```
    <--- This fix is no longer necessary in the non-Conda env!

* Need to find a way how to link FFTW3 library without giving the explicit path (maybe unavoidable)
    <--- Setting the proper environment variables LIBRARY_PATH

* Problem on Quantum cluster:
    ```Fatal Error: Can't open module file 'linked_list.mod' for reading at (1)```
    Has to be addressed! After a few attempts it compiles OK.
    <--- Addressed by including the file during the compilation

* The makefile is built using the following command:
    ```cmake -D CMAKE_Fortran_COMPILER=mpifort ..```
    <--- Not necessary with proper installation of the packages (non-Conda env).

* In order to compile the code using non-GNU standard (flag -std=gnu, -std=2003 etc.),
    several modifications had to be done. First, non-gnu standards do not support
    the following lines:
    ``` REAL (IMAG (Z)) ```, instead they support ```REAL (AIMAG (Z))```.

    Next, the variable buffer is defined as 32-bit float, however it was compared
    with 64-bit float ```0.D0```. Rewriting to ```0.``` fixed the issue.

    Last, the X descriptor was incorrectly written in ```write_listing.f90``` file.
    Instead of ```WRITE(100,'(a,t50,es12.4, x,a)')``` we need to add '1' in front of 
    'x' so the correct form is: ```WRITE(100,'(a,t50,es12.4, 1x,a)')```

* The code must be compiled in the default Mac environment, NOT Conda!!! With 
    the change of the environment the code compiles okay without specifying 
    the compiler explicitly during cmake call, now it is sufficient to compile 
    the code using ```cmake ..```.

* Predefining CMAKE_Fortran_COMPILER is probably a good approach

* Invoking CMake with the following commands yields more information about linking
    cmake -D CMAKE_Fortran_COMPILER=mpifort -D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE ..

* Local environment variables have to be set correctly otherwise it needs to be
    set manually in CMakeLists. This is the case of the FFTW3 library. One way is
    to set explicitly 
        ```link_directories(/usr/local/lib)```
    or check that the environment variable LIBRARY_PATH is set to a location 
    of FFTW libs. This can be set as follows:
        ```export LIBRARY_PATH=/usr/local/lib```