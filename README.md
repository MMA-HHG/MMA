# MMA for HHG modelling

Table of content

## Installations
Both `CUPRAD` and `CTDSE` are compiled from the source. Hankel is implemented in Python and together with the other Pythonic procedures and scripting around the model rely on setting up the environment variables (one might consider setting up [a virtual environment](https://stackoverflow.com/questions/9554087/setting-an-environment-variable-in-virtualenv) dedicated to this whole project).
<!--- We head for RC1, use virtual environment, pip, etc. in future releases by default --->  

Here is the list of requirements
* **CUPRAD**: MPI, FFTW3, parallel HDF5, CMake
* **CTDSE**: FFTW3, MPI, serial HDF5, CMake (MPI & HDF5 are not needed for the dynamic library), Python with h5py
    * Python ctypes library for the interactive CTDSE
* **Hankel**: numpy, scipy, h5py, multiprocessing + some usual Python libraries

We have not found any particular requirements for the versions of the libraries, and the code was successfully build using intel, GNU and AppleClang compilers. Alas some specific flags and settings are required for different compilers as discussed below.

Note that intel encapsulates FFTW3 into the [Math Kernel Library](https://en.wikipedia.org/wiki/Math_Kernel_Library).



<a id="setting-the-paths"></a>
### Setting the paths, installing libraries and modules
Here is the list of paths for running the model. The only customised path is the `$GIT_PATH`, which points to the parent directory of different git repositories used in the model. All the other paths are set relatively ot this path. Additioanlly to the paths, there are also bash functions for loading the necessary modules.


#### Paths
``` bash
export GIT_PATH=/users/k2255939/git # This is the only path that needs to be customised

# These paths are relative to the GIT_PATH
export UNIV_INPUT_PATH=$GIT_PATH/universal_input
export MSM_PATH=$GIT_PATH/CUPRAD_TDSE_Hankel
export PYTHONPATH=$PYTHONPATH:$GIT_PATH/python_modules:$MSM_PATH/shared_python
export CUPRAD_HOME=$MSM_PATH/CUPRAD
export CUPRAD_BINARY=$MSM_PATH/CUPRAD/binary
export CUPRAD_SCRIPTS=$MSM_PATH/CUPRAD/scripts
export CUPRAD_INPUTS=$MSM_PATH/CUPRAD/testing
export CUPRAD_PYTHON=$MSM_PATH/CUPRAD/python
export TDSE_1D_SOURCE=$MSM_PATH/1DTDSE/
export TDSE_1D_POST_PROCESSING=$MSM_PATH/1DTDSE/post_processing
export TDSE_1D_SCRIPTS=$MSM_PATH/1DTDSE/scripts
export TDSE_1D_SLURM=$MSM_PATH/1DTDSE/slurm
export TDSE_1D_BINARY=$MSM_PATH/1DTDSE/binary
export HANKEL_HOME=$MSM_PATH/Hankel
export TDSE_1D_HOME=$MSM_PATH/1DTDSE
export MULTISCALE_HOME=$MSM_PATH
export MULTISCALE_SCRIPTS=$MSM_PATH/multiscale/scripts
export PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON
```


#### Modules and libraries
When used locally on a personal computer, the libraries (FFTW3, CMake, …) are typically installed by a user. If using the code or a part of it, the corresponding libraries must be installed before. (Sharing needs to be enebled for FFTW3, see details below for the dyamic-library CTDSE.)

[The modules](https://hpc-wiki.info/hpc/Modules) provide all the necessary libraries for the code when using a computational cluster. The script `Modules/load_modules.sh` is used to load all the modules. There is a list of modules for various computational clusters specified by the variable `$HPC`. Another supercomputer (or compilation option intel/GNU/...) should be added there.

There are two `bash` functions `load_modules` and `load_python_modules`. The former is activated when running *CUPRAD* and *CTDSE*, while the latter is used for all Pythonic operations around the code. (The reason for this duality is that Python might need to load a compiler itself for some libraries, typically on intel.) 



### CUPRAD
All the source files are located in `CUPRAD/sources`. The CMake recipe is in `CUPRAD/CMakeLists.txt`. The code is supposedly built in `CUPRAD/build`.
There is the recipe for compilation (each point contains several notes about possible difficulties):

1) Run `load_modules`. [This can be verified by](https://hpc-wiki.info/hpc/Modules#:~:text=%24-,module%20list,-Currently%20Loaded%20Modulefiles) `module list`.
    * If the machine does not using modules, this step is replaced by installing the necessary libraries and setting up the environment.
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory.
    * We encountered CMake struggling to identify the proper MPI-Fortran compiler on several machines. CMake can be hinted to use the desired compiler by `cmake -D CMAKE_Fortran_COMPILER=mpifort ..` (GNU) or `cmake -D CMAKE_Fortran_COMPILER=mpiifort ..` (intel). This resolved the issue when we encountered it.
    * The CMake configuration can be manually adjusted using `ccmake`, see [link 1](https://cmake.org/cmake/help/latest/manual/ccmake.1.html) and [link 2](https://stackoverflow.com/a/1224652).
3) Compile the code from the CMake-generated `Makefile` by running `make code` in the `build` directory.


### CTDSE (+ the dynamic library)
All the source files are located in `1DTDSE/sources`. The CMake recipe is in `1DTDSE/CMakeLists.txt`. This recipe installs both the code and the interactive CTDSE library; if MPI is not present, only the dynamic library is installed. The code is supposedly built in `1DTDSE/build`.
There is the recipe for compilation:

1) Run `load_modules` or install the libraries if used on a personal computer.
    * The extension of the library may depend on your platform: `libsingleTDSE.so`, `libsingleTDSE.dynlib`, `libsingleTDSE.dll`, …
    * The fftw3 library needs to be installed as a shared library! See [link 1](https://www.fftw.org/fftw2_doc/fftw_6.html#:~:text=Note%20especially%20%2D%2Dhelp%20to%20list%20all%20flags%20and%20%2D%2Denable%2Dshared%20to%20create%20shared%2C%20rather%20than%20static%2C%20libraries.%20configure%20also%20accepts%20a%20few%20FFTW%2Dspecific%20flags%2C%20particularly), [link 2](https://stackoverflow.com/a/45327358).
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory.
    * The specification of the compiler might be neded similarly to CUPRAD: `cmake -D CMAKE_C_COMPILER=mpicc ..` or `cmake -D CMAKE_C_COMPILER=mpiicc ..`
    * (Intel:) The FFTW3 library might not be found within the MKL and might be needed to link manually by adding it into the environment: ```export CPATH=${CPATH}:${MKLROOT}/include/fftw``` (the location of fftw is not consistent across MKL versions and `fftw` needs to located within `$MKLROOT$`).
3) Compile the code by running `make` in the `build` directory.
4) Check that the `$PYTHONPATH` includes `1DTDSE/post_processing`, so the Pythonic scripting around the module works.

#### Local installation of the dynamic library on Ubuntu 22.04 (using WSL)
Here is an example of installing the interactive CTDSE library on Ubuntu 22.04 as a part of [The Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install). This option is convenient also for Windows users because WSL is also accessible directly from the Windows environment (for example by using VSCode).

* The following prerequisities are needed: ``sudo apt-get install build-essential``, ``sudo apt-get install libhdf5-dev``, `` sudo apt-get install hdf5-helpers``.
* The fftw3 library needs to be installed as a shared library! See [link 1](https://www.fftw.org/fftw2_doc/fftw_6.html#:~:text=Note%20especially%20%2D%2Dhelp%20to%20list%20all%20flags%20and%20%2D%2Denable%2Dshared%20to%20create%20shared%2C%20rather%20than%20static%2C%20libraries.%20configure%20also%20accepts%20a%20few%20FFTW%2Dspecific%20flags%2C%20particularly), [link 2](https://stackoverflow.com/a/45327358).
* `cmake` may have problems with `h5cc` wrapper (it worked with v3.15.7).
* Install using the previous recipe.

The dynamic library is then available, see [this example of a jupyter notebook integrating CTDSE](xxx).

### Hankel
This module becomes available by [including it into the `$PYTHONPATH`](#setting-the-paths).


## Inputs


### Global

### CUPRAD

### CTDSE

### Hankel

## Execution pipeline


# Standalone operations




# To do list
* Code + documentation
    * ***I/O handling + HDF5 structure***
        * 3 groups for each module + one global group for shared parameters (gas type) 
    * *Make compilation smooth and easy on various machines (ongoing)*
    * ***CUPRAD*** (Jan + Stefan), **[Doxygen]((https://www.doxygen.nl/)) documentation**
        * Jan - `output.f90` is missing
        * Stefan - the core routines
        * **pre-processor**: potentially very tedious, do up to some extent
    * ***1D-TDSE*** (Tadeáš) **Doxygen documentation** ✔️
    * ***interactive TDSE*** (Tadeáš)
        * Ensure compatibility of the dll library
    * ***Hankel*** (Jan) **Python docstrings** (process by [sphynx](https://www.sphinx-doc.org/en/master/tutorial/index.html))
        * *refactoring*
            * Implement density modulation
                * *Pythonic ecosystem around the whole code*
            * hard-code preset gases with their scattering factors
            * data processing (either load all $\mathscr{F}$-sources or chunk reading from the drive) ✔️
    * ***Assembling the code***
        * Work with all the modules
        * Make the execution smooth (slurm dependecies etc.)
        * Fine-adjustments of all codes
            * HDF5-organisation
            * array orderings
    * ***Pre/post-processing, data plotting***
        * organise Pythonic libraries, integrare `mynumerics` int the model, ...
        * *Jupyter notebook* showing processing the outputs
    * ***CodeOcean capsule*** (Tadeáš + Jan)
    * Enhance code + executables (compilation with a 'pedantic' compiler, valgrind) (to some extent)
* Prepare test cases (integrate within the CodeOcean capsule)
    * *HHG in gas jet*
    * *HHG in gas cell* (either this or the previus with a Gaussian density-modulation.)
    * *Absorption-limited HHG wit pre-ionisation in gas cell*
    * *HHG with Bessel beams*
    * *Interactive python* ✔️ via Jupyter
        * ***Promote invariant ionisation***
    * *?* Some filamentation form Stefan's papers, ...
* Writing the paper (all)
    * CPC: no new physics, only the code and its usage + already published results
    * see below independent repo + Overlaf
* Test cases & tutorials
    * ***Create a YouTube video-tutorial***
        * Go through the whole code installation and execution, comment it.
    * (?) Jupyter notebook  
* Documentation
    * README.md's in the git-repo + webpage genrated by doxygen + sphynx
* Check conventions
    * The conventions in the code shall correspond the written ones in supplmentary materials. I did checks during the development, that the code works intrinsically. However, there might be a mistake in the sign of the $e^-$ charge in TDSE etc.
    * Photoelectron spectrum norm.

## Notes
Here are notes to the development of the whole model. The development of the particular modules - ***CUPRAD***, ***1D-TDSE*** and ***Hankel*** - are in the `README.md`'s in their respective directories.

The work on the publication is being done in an independet [repo](https://github.com/sskupin/XUVIR) morrored with its [Overleaf](https://www.overleaf.com/project/64f9ee5fb00b8a641fe54780).

# General notes
## PYTHONPATH organisation

Several Python procedures require modules from various places. The actual approach is to include all the locations in PYTHONPATH

- There is needed to inlcude the independent git repository `python_modules`
- The path should also contain `shared_python` from the root of this repository
- There is an example of the paths added into the environment. The environment variables are used in some Python scripts, so they shall be added into `.bashrc` to ensure compatibility:
``` bash
export GIT_PATH=/users/k2255939/git

# export GIT_PATH=/gpfs/home/jvabek/git

export UNIV_INPUT_PATH=$GIT_PATH/universal_input

export MSM_PATH=$GIT_PATH/CUPRAD_TDSE_Hankel

# # export TESTPATH=/home/vabekjan/git/CUPRAD_DEVELOP/testing/job_scheduler/inps

export PYTHONPATH=$PYTHONPATH:$GIT_PATH/python_modules:$MSM_PATH/shared_python

export CUPRAD_HOME=$MSM_PATH/CUPRAD

export CUPRAD_BINARY=$MSM_PATH/CUPRAD/binary

export CUPRAD_SCRIPTS=$MSM_PATH/CUPRAD/scripts

export CUPRAD_INPUTS=$MSM_PATH/CUPRAD/testing

export CUPRAD_PYTHON=$MSM_PATH/CUPRAD/python

export TDSE_1D_SOURCE=$MSM_PATH/1DTDSE/

export TDSE_1D_POST_PROCESSING=$MSM_PATH/1DTDSE/post_processing

export TDSE_1D_SCRIPTS=$MSM_PATH/1DTDSE/scripts

export TDSE_1D_SLURM=$MSM_PATH/1DTDSE/slurm

export TDSE_1D_BINARY=$MSM_PATH/1DTDSE/binary

export HANKEL_HOME=$MSM_PATH/Hankel

export TDSE_1D_HOME=$MSM_PATH/1DTDSE

export MULTISCALE_HOME=$MSM_PATH

export MULTISCALE_SCRIPTS=$MSM_PATH/multiscale/scripts

# # export PYTHON_CUPRAD=/home/vabekjan/git/CUPRAD_DEVELOP/python

export PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON
```

These aliases may be separated from `.bashrc` by including them in `.bash_aliases` via
``` bash
if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi
```




## Merging independent repos into the multi-scale model:
https://stackoverflow.com/a/52933095
