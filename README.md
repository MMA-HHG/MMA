
<img src="misc/logos/eli.png" alt="ELI logo" style="height:60px;"/> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
<img src="misc/logos/Logo-celia-bleu-sans-slogan-300x184.png" alt="CELIA logo" style="height:60px;"/> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
<img src="misc/logos/ilm.jpg" alt="iLM logo" style="height:60px;"/> &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;
<img src="misc/logos/CTU_lion-1_tr.png" alt="CTU logo" style="height:60px;"/>

# Modular multiscale approach (MMA) for HHG modelling
This repository contains the sources for the [**Modular Multiscale Approach (MMA) for High-Harmonic Generation (HHG) modelling**](https://mma-hhg.github.io/) and tutorials showing examples of usage. It is still under active development and there will be frequent updates.

The model provides a full solution to treat HHG in gaseous media: the laser pulse propagation, the single-atom response, and the XUV propagation. The design of the code is modular and each of the respective modules can be used independently for its own purposes without requiring the other modules. The assumtions in the model are the linear polarisation and cyllindrical symmetry in the macroscopic propagation. (We are working on extensions of the approach, see notes below.)

[**An installation guide is below**](#installations), followed by [the description of input parameters](#inputs), a general overview is available in [**the pre-print**](https://arxiv.org/pdf/2507.04115). [***YouTube tutorials***](https://youtube.com/playlist?list=PLhgkvMYpHD8gA1C3FM_GB24HNYNsBGxQe&feature=shared) for the main operations are available. If you have difficulties with code, please [***raise an Issue***](https://github.com/MMA-HHG/MMA/issues) or [***start a Disscusion***](https://github.com/MMA-HHG/MMA/discussions) or contact us. The compiled version of the code is accessible through [***CodeOcean capsule***](https://codeocean.com/capsule/6775529/tree). Different modes of the operation of the code are provided in [***the jupyter tutorials***](./jupyter_examples/). The reference datasets for these demos are available [here](https://elibeamlines-my.sharepoint.com/:f:/g/personal/jan_vabek_eli-beams_eu/EpgdhQS_bQROnTYpj6y0zgYBTvz_sSh3WMQDiWgy9YgXfw?e=WnnpgR). The respective directories provide more technical details about the codes: [CUPRAD](./CUPRAD/), [1D-TDSE](./1DTDSE/README.md) and [the Hankel transform](./Hankel/).

The last point to address is the execution of the code. If used locally, CUPRAD and TDSE are standard MPI applications and Hankel a Pythonic module. Once running on HPC clusters, a scheduler is involved. The code is then executed as [a pipeline of the different modules](#execution-pipeline).

**Authors and contacts**: \
Jan Vábek: [Jan.Vabek@eli-beams.eu](Jan.Vabek@eli-beams.eu) \
Fabrice Catoire: [fabrice.catoire@u-bordeaux.fr](fabrice.catoire@u-bordeaux.fr) \
Stefan Skupin: [stefan.skupin@univ-lyon1.fr](stefan.skupin@univ-lyon1.fr) \
Tadeáš Němec

**Notes and disclaimers:**
* [***Subscribe for news!***](https://bit.ly/HHG-code-updates
)
* The code is provided as it is with open source. We cannot provide guarantee or take responsibility for its usage.
* We are a small developer's group and our primary occupation is science. We would be grateful to discuss the usage of the code. However, we cannot provide a commerce-level full-scale support at the instant.
* We would be grateful for your feedback!
* We are working on more advanced siblings of the codes (3D-vectorial pulse propagation, 3D-TDSE) that are not a part of the first release. Please contact authors for possible collaborations if you need the more advanced features.

## Table of contents
- [Modular multiscale approach (MMA) for HHG modelling](#modular-multiscale-approach-mma-for-hhg-modelling)
  - [Table of contents](#table-of-contents)
  - [Installations](#installations)
    - [**Recommended way** - Reference Docker installation](#recommended-way---reference-docker-installation)
    - [**Advanced** - Custom installations](#advanced---custom-installations)
      - [Installation steps](#installation-steps)
    - [HPC installation using Modules environment](#hpc-installation-using-modules-environment)
    - [CUPRAD installation](#cuprad-installation)
      - [Note: building only the pre-processor](#note-building-only-the-pre-processor)
    - [CTDSE (+ the dynamic library) installation](#ctdse--the-dynamic-library-installation)
    - [Hankel](#hankel)
  - [Running the code](#running-the-code)
  - [HDF5 data organisation](#hdf5-data-organisation)
  - [Inputs](#inputs)
    - [Global](#global)
    - [CUPRAD](#cuprad)
    - [CTDSE](#ctdse)
    - [Hankel](#hankel-1)


<a id="installations"></a>
## Installations
Both `CUPRAD` and `CTDSE` are compiled from the source. Hankel is implemented in Python and together with the other Pythonic procedures and scripting around the model rely on setting up the environment variables.
<!--- We head for RC1, use virtual environment, pip, etc. in future releases by default --->

We consider two usages os the code. First, [we provide full local conteinerised insallation within a Docker image](#recommended-way---reference-docker-installation). This approach is recommended to get familiar with the code and some smaller tasks. To upsacale the work, it is necessary to use HPC, this installation is usually machine specific and [we provide some a general overview of the used libraries, which shall together with the reference local container facilitate the usage of the code](#advanced---custom-installations).

In both cases, this repository is cloned by `git clone git@github.com:MMA-HHG/MMA.git` or `git clone https://github.com/MMA-HHG/MMA.git`.

The first option uses the Docker image. It is a direct multiplatform user-oriented way to obtain the executable model. This can be used for running the model locally. Moreover, this can be used as a direct reference for compiling the code the second way: directly from the source. This option might be neccessary for deploying the code on HPC clusters, develpoment, …

<a id="recommended-way---reference-docker-installation"></a>
### **Recommended way** - Reference Docker installation

The code can be accessed through Docker. Docker prepares the software environment, while the MMA repository is mounted into the container. The native executables are compiled automatically when the container is started for the first time. You can [**go directly to the first-installation commands.**](#install-direct)

1) [Docker](https://www.docker.com/) needs to be installed.

2) Clone the repository and enter its root directory.

        git clone git@github.com:MMA-HHG/MMA.git
        cd MMA

    If you do not use SSH keys with GitHub, clone the repository using https by `git clone https://github.com/MMA-HHG/MMA.git` or using any convenient way as [GitHub Desktop](https://desktop.github.com/download/).

3) Build the Docker image.

        docker build -t mma .

    The option `-t mma` specifies the name of the Docker image, i.e. in this case the name of the image is `mma`.

    *It is recommended to use [WSL in Windows](https://learn.microsoft.com/en-us/windows/wsl/install). Without WSL, it might be necessary to replace [Dockerfile](./Dockerfile) by the [Dockerfile for Windows](./docker/Dockerfile_Windows) due to a different character set.*

4) Start the container. Choose `CONTAINER_NAME`, for example `mma-c1`. The same name is also used as the container hostname, so the terminal prompt is easier to read.

        docker run --name CONTAINER_NAME --hostname CONTAINER_NAME -v .:/MMA -w /MMA -p 8888:8888 -it mma

   The repository is mounted into `/MMA`, so the compiled executables remain available in the parent filesystem. If port `8888` is already used on the  parent machine, change only the first number, for example `-p 8889:8888`.

5) JupyterLab can be started inside the container by:

        mma-jupyter

6) Two prepared tutorial workspaces are available through:

        teach-me-mma
        teach-me-tdse

   These commands print the corresponding JupyterLab workspace links.

7) To return to the same container later, use:

        docker start -ai CONTAINER_NAME

In short, the installation is:
<a id="install-direct"></a>
```bash
git clone git@github.com:MMA-HHG/MMA.git
cd MMA

docker build -t mma .
docker run --name mma-c1 --hostname mma-c1 -v .:/MMA -w /MMA -p 8888:8888 -it mma
```

Inside the container, start JupyterLab by:

```bash
mma-jupyter
```

and follow the instructions shown by these tutorials
```bash
teach-me-mma
teach-me-tdse
```
to learn basics of the code.

To return to the container later:

```bash
docker start -ai mma-c1
```

<a id="advanced---custom-installations"></a>
### **Advanced** - Custom installations
Here we provide a more detailed guide for the installation, this is intended for running on HPC clusters and developers. From our experience, each HPC cluster can use a different approaches, here can be found [some examples for different machines](./Modules/load_modules.sh). However, some interaction with HPC admin team could be necessary for compiling the code properly. The following specifications toghether with the reference DOcker installation shall serve as a guide to facilitate these deployments.

This is the list of requirements:
* **CUPRAD**: MPI, FFTW3, parallel HDF5 (with Fortran support and HL libraries), CMake
* **CTDSE**: FFTW3, MPI, serial HDF5, CMake (MPI & HDF5 are not needed for the dynamic library), Python with h5py
* **Hankel**: numpy, scipy, h5py, multiprocessing + some usual Python libraries

We have not found any particular requirements for the versions of the libraries, and the code was successfully build using intel, GNU and AppleClang compilers. Alas some specific flags and settings are required for different compilers as discussed below.

Note that intel encapsulates FFTW3 into the [Math Kernel Library](https://en.wikipedia.org/wiki/Math_Kernel_Library).

#### Installation steps
This guide assumes the user has clean installation of Linux (e.g. Ubuntu with [gcc](https://ubuntu.com/developers/docs/howto/gcc-setup/)), MacOS (with [XCode Command Line Tools](https://developer.apple.com/documentation/xcode/installing-the-command-line-tools) and [Homebrew](https://brew.sh)), or WSL subsystem on Windows with Ubuntu. Due to specific requirements for the dependencies, we recommend following this guide, and use default system paths for installing the dependencies.

To build the full MMA bundle locally from scratch with dependencies (using autotools and CMake), follow these steps:

1. Install the dependencies (skip if already acquired):
   - **Fortran Lang**:

      Ubuntu:
      ```bash
      apt update && sudo apt install gfortran
      ```
      Mac:
      ```bash
      brew install gfortran
      ```
   - **MPI** ([Official guide](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html))
      ```bash
      # Choose default download path
      prefix=~/Downloads

      # Download MPI - choose distribution here: https://www.open-mpi.org/software/ompi/v5.0/
      curl -L -o ${prefix}/openmpi-5.0.10.tar.gz "https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.10.tar.gz"

      # unpack
      tar xvf ${prefix}/openmpi-5.0.10.tar.gz -C ${prefix}
      cd ${prefix}/openmpi-5.0.10

      # configure installation:
      #  - we keep the default system path for installation
      #  - default path installation requires superuser (sudo)
      #  - in case of missing admin rights, add flag `-prefix=/custom/path`
      ./configure 2>&1 | tee config.out

      # Compile using 4 cores
      make -j4
      sudo make install
      ```
   - **HDF5 Parallel** ([official guide](https://github.com/HDFGroup/hdf5/blob/develop/docs/README_HPC.md)):

      Universal guide:
      ```bash
      prefix=~/Downloads

      # Download HDF5 - choose distribution here: https://github.com/HDFGroup/hdf5/releases
      curl -L -o ${prefix}/hdf5-2.2.0.tar.gz "https://github.com/HDFGroup/hdf5/releases/download/2.2.0/hdf5-2.2.0.tar.gz"

      # unpack
      tar xvf ${prefix}/hdf5-2.2.0.tar.gz -C ${prefix}
      cd ${prefix}/hdf5-2.2.0

      # configure installation:
      cmake -DCMAKE_BUILD_TYPE=Release \
        -DHDF5_ENABLE_PARALLEL=ON \
        -DBUILD_SHARED_LIBS=ON \
        -DHDF5_BUILD_FORTRAN=ON \
        -DHDF5_BUILD_HL_LIB=ON \
        ..

      # Compile and install
      cmake --build . --config Release
      sudo cmake --install .
      ```

      Simplified Ubuntu 22.04 approach:
      ```bash
      sudo apt install libhdf5-dev hdf5-helpers
      ```
   - **FFTW3** (GNU version, [official guide](https://www.fftw.org/fftw2_doc/fftw_6.html)):
      ```bash
      # Choose default download path
      prefix=~/Downloads

      # Download MPI - choose distribution here: https://www.open-mpi.org/software/ompi/v5.0/
      curl -L -o ${prefix}/fftw-3.3.11.tar.gz "http://www.fftw.org/fftw-3.3.11.tar.gz"

      # unpack
      tar xvf ${prefix}/fftw-3.3.11.tar.gz -C ${prefix}
      cd ${prefix}/fftw-3.3.11

      # configure installation:
      #  - we need to enable shared flag for dynamic linking of the library
      ./configure --enable-shared 2>&1 | tee config.out

      # Compile using 4 cores
      make -j4
      sudo make install
      ```
   - **CMake**:

      Ubuntu:
      ```bash
      apt update && sudo apt install cmake
      ```
      Mac:
      ```
      brew install cmake
      ```

   - **Python3**:

      Ubuntu:
      ```bash
      # for best compatibility we recommend Python >= 3.12
      apt update && sudo apt install python@3.14
      ```

      Mac:
      ```bash
      brew install python@3.14
      ```

   - Git:

      Ubuntu:
      ```bash
      apt update && sudo apt install git
      ```
      Mac:
      ```bash
      brew install git
      ```

2. Build the MMA code:
    ```bash
    # set git repository path
    prefix=/path/to/repository
    cd "${prefix}"

    # clone the repository
    git clone https://github.com/MMA-HHG/MMA.git
    cd MMA

    # prepare build
    mkdir build
    cd build

    # CMake and build
    # - specify compilers with mpi wrappers
    # - specify HDF5 path - by default in /usr/local/HDF_Group
    cmake .. \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_Fortran_COMPILER=mpifort \
      -DCMAKE_PREFIX_PATH=/usr/local/HDF_Group/HDF5/2.2.0

    # build
    make
    ```

3. Build Python virtual environment
    ```bash
    python -m venv venv
    source venv/bin/activate

    # install required Python modules
    pip install -r requirements.txt
    ```

### HPC installation using Modules environment
When used locally on a personal computer, the libraries (FFTW3, CMake, …) are typically installed by a user, see [previous section](#advanced---custom-installations). Installing the code on the HPC infrastructure usually requires loading dependencies using pre-installed modules.

[The modules](https://hpc-wiki.info/hpc/Modules) provide all the necessary libraries for the code when using a computational cluster. The script [`Modules/load_modules.sh`](./Modules/load_modules.sh) is used to load all the modules. This function are supposed to be added into the environemnt (e.g. by sourcing the script in `.bash_aliases` as done above). There is a list of modules for various computational clusters specified by the variable `$HPC`. Another supercomputer (or compilation option intel/GNU/...) should be added there.

There are two `bash` functions `load_modules` and `load_python_modules`. The former is activated when running *CUPRAD* and *CTDSE*, while the latter is used for all Pythonic operations around the code. (The reason for this duality is that Python might need to load a compiler itself for some libraries, typically on intel.)

If everything is set well, the following CMakes are wrapped in the master `CMakeList.txt` (the compilers tested on some machines and currently supported are GNU and Intel). Here is the recipe to install the code from its root directory.

1) Run `load_modules`. [This can be verified by](https://hpc-wiki.info/hpc/Modules#:~:text=%24-,module%20list,-Currently%20Loaded%20Modulefiles) `module list`.
    * If the machine does not using modules, this step is replaced by installing the necessary libraries and setting up the environment.
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory. We encountered CMake sometimes struggling to identify the proper MPI-Fortran compiler on several machines (the error is raised in the next step). There are more ways to hint CMake to find the compilers:
    * By providing environment variables with the compilers: `export CC=mpicc` and `export FC=mpifort` (GNU); `export CC=mpiicc` and `export FC=mpiifort` (Intel).
    * Controlling CMake directly during its execution `cmake -D CMAKE_Fortran_COMPILER=mpifort ..` (GNU) or `cmake -D CMAKE_Fortran_COMPILER=mpiifort ..` (intel). This resolved the issue when we encountered it.
    * The CMake configuration can be manually adjusted using `ccmake`, see [link 1](https://cmake.org/cmake/help/latest/manual/ccmake.1.html) and [link 2](https://stackoverflow.com/a/1224652).
    * Consider to run the compilation of the respective codes separately in the case problems occur.
3) Compile the code from the CMake-generated `Makefile` by running `make code` in the root directory.


Below are recipes for compilling the codes separately.

### CUPRAD installation
All the source files are located in `CUPRAD/sources`. The CMake recipe is in `CUPRAD/CMakeLists.txt`. The code is supposedly built in `CUPRAD/build`.
There is the recipe for compilation (each point contains several notes about possible difficulties):

1) Run `load_modules`. [This can be verified by](https://hpc-wiki.info/hpc/Modules#:~:text=%24-,module%20list,-Currently%20Loaded%20Modulefiles) `module list`.
    * If the machine does not using modules, this step is replaced by installing the necessary libraries and setting up the environment.
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory.
    * We encountered CMake struggling to identify the proper MPI-Fortran compiler on several machines. CMake can be hinted to use the desired compiler by `cmake -D CMAKE_Fortran_COMPILER=mpifort ..` (GNU) or `cmake -D CMAKE_Fortran_COMPILER=mpiifort ..` (intel). This resolved the issue when we encountered it.
    * (Intel:) The FFTW3 library might not be found within the MKL and might be needed to link manually by adding it into the environment: ```export CPATH=${CPATH}:${MKLROOT}/include/fftw``` (the location of fftw is not consistent across MKL versions and `fftw` needs to be located within `$MKLROOT$`).
    * The CMake configuration can be manually adjusted using `ccmake`, see [link 1](https://cmake.org/cmake/help/latest/manual/ccmake.1.html) and [link 2](https://stackoverflow.com/a/1224652).
3) Compile the code from the CMake-generated `Makefile` by running `make` in the `build` directory.


#### Note: building only the pre-processor
The pre-processor can built without MPI and parallel HDF5 (serial HDF5 is needed). The cmake should enter this option if no MPI is found.

### CTDSE (+ the dynamic library) installation
All the source files are located in `1DTDSE/sources`. The CMake recipe is in `1DTDSE/CMakeLists.txt`. This recipe installs both the code and the interactive CTDSE library; if MPI is not present, only the dynamic library is installed. The code is supposedly built in `1DTDSE/build`.
There is the recipe for compilation:

1) Run `load_modules` or install the libraries if used on a personal computer.
    * The extension of the library may depend on your platform: `libsingleTDSE.so`, `libsingleTDSE.dynlib`, `libsingleTDSE.dll`, …
    * The fftw3 library needs to be installed as a shared library! See [link 1](https://www.fftw.org/fftw2_doc/fftw_6.html#:~:text=Note%20especially%20%2D%2Dhelp%20to%20list%20all%20flags%20and%20%2D%2Denable%2Dshared%20to%20create%20shared%2C%20rather%20than%20static%2C%20libraries.%20configure%20also%20accepts%20a%20few%20FFTW%2Dspecific%20flags%2C%20particularly), [link 2](https://stackoverflow.com/a/45327358).
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory.
    * The specification of the compiler might be neded similarly to CUPRAD: `cmake -D CMAKE_C_COMPILER=mpicc ..` or `cmake -D CMAKE_C_COMPILER=mpiicc ..`
    * (Intel:) The FFTW3 library might not be found within the MKL and might be needed to link manually by adding it into the environment: ```export CPATH=${CPATH}:${MKLROOT}/include/fftw``` (the location of fftw is not consistent across MKL versions and `fftw` needs to be located within `$MKLROOT$`).
3) Compile the code by running `make` in the `build` directory.
4) Check that the `$PYTHONPATH` includes `1DTDSE/post_processing`, so the Pythonic scripting around the module works.

### Hankel
This module becomes available by [including it into the `$PYTHONPATH`](#setting-the-paths).

## Running the code

The model consists of three main jobs: 1) CUPRAD for the laser pulse propagation; 2) TDSE for the microscopic response; and 3) the Hankel transform for the far-field XUV distribution. There are some further auxiliary steps in the pipeline. **[This jupyter tutorial is desinged to guide the first execution of the code.](./jupyter_examples/mma_basics/teach_me_mma.ipynb)** The guide to open this tutorial though a jupyter server is shown by the command `teach-me-mma` from the Docker terminal. You can [**go directly to the pipeline.**](#run-direct)

The pipeline consists.

0) Sourcing the installation and run paths:
    ```bash
    # Root git path of MMA
    source ./set_env_vars.sh
    ```

1) CUPRAD pre-processor (`$CUPRAD_BUILD/make_start.e`)

        "$CUPRAD_BUILD/make_start.e" INPUT.h5

    * The pre-processor requires the name of the HDF5 input file.
    * The name of the file is then stored in `msg.tmp`, which transfers it
      through the execution pipeline.

2) Main MPI CUPRAD job (`$CUPRAD_BUILD/cuprad.e`)

        mpirun -n "$NUM_PROC_DEFAULT_CUPRAD" "$CUPRAD_BUILD/cuprad.e"

    * The design of the code requires the number of MPI processes to be a power
      of two.
    * Inside the Docker container, the default value is provided by
      `$NUM_PROC_DEFAULT_CUPRAD`.

3) Adjusting the TDSE parameters to the real number of steps in `z`
   (`$TDSE_1D_PYTHON/prepare_TDSE_Nz.py`)

        python3 "$TDSE_1D_PYTHON/prepare_TDSE_Nz.py"

4) Main MPI TDSE job (`$TDSE_1D_BUILD/TDSE.e`)

        mpirun -n "$NUM_PROC_DEFAULT_TDSE_1D" "$TDSE_1D_BUILD/TDSE.e"

    * In contrast to CUPRAD, TDSE does not require the number of MPI processes
      to be a power of two.
    * Inside the Docker container, the default value is provided by
      `$NUM_PROC_DEFAULT_TDSE_1D`.

5) Merge and clean the temporary TDSE files (`$TDSE_1D_PYTHON/merge.py`)

        python3 "$TDSE_1D_PYTHON/merge.py"

6) Hankel transform (`$HANKEL_HOME/Hankel_long_medium_parallel_cluster.py`)

        python3 "$HANKEL_HOME/Hankel_long_medium_parallel_cluster.py"

    * ***It has to be executed as a single multithreaded program.***
    * It uses multithreading parallelisation using [the multiprocessing library](https://docs.python.org/3/library/multiprocessing.html). *The number of threads has to be defined in the input hdf5-archive. Be careful, especially on HPC's, that these numbers match with hardware.*


In short, the pipeline can be run inside the Docker container as:
<a id="run-direct"></a>
```bash
$CUPRAD_BUILD/make_start.e INPUT.h5
mpirun -n $NUM_PROC_DEFAULT_CUPRAD $CUPRAD_BUILD/cuprad.e

python3 $TDSE_1D_PYTHON/prepare_TDSE_Nz.py
mpirun -n $NUM_PROC_DEFAULT_TDSE_1D $TDSE_1D_BUILD/TDSE.e
python3 $TDSE_1D_PYTHON/merge.py

python3 $HANKEL_HOME/Hankel_long_medium_parallel_cluster.py
```

Here `INPUT.h5` should be replaced by the prepared input file.

It is possible to run the process manually. However, computational clusters use jobs and queues for scheduling them. [Here](./multiscale/scripts/README.md) we discuss an example of this pipeline.


## HDF5 data organisation
All the inputs and outputs are stored in an HDF5 archive grouping together I/O for the different modules. The structure is flexible and defined in the *namelists*: [Python](shared_python/MMA_administration.py), [C](1DTDSE/sources/h5namelist.h), [Fortran](CUPRAD/sources/global_variables.f90).

## Inputs
Here is the exhaustive list of all the parameters. The bold **`parameters`** are obligatory to run the whole model with sourcing the default material constants. The other `parameters` are optional. If an optional parameter is present, it has priority.

### Global
The global inputs are stored in the `global_inputs` groups, they might be used by more than one module.
* `gas_preset`: The main specifier to define the gas. It defines all the material constants then. Implemented gases: `He`, `Ne`, `Ar`, `Kr`, `Xe`.
* **`medium_pressure_in_bar`**: Pressure of the medium in bars. (Preferrably defined it here for the multiscale usage.)
* **group `density_mod`**: This group defines the density modulation. The relative modulation to `medium_pressure_in_bar` is stored in `table` with the dimension corresponding to the grids (it can be 1- or 2-dimensional). Therefore, `zgrid` and/or `rgrid` needs to be present. The grids should be larger than the interaction volume. Note: the radial modulation is not fully available for Hankel (see detailed documentation of the module for details).
  * `zgrid`: The grid corresponding to the longitudinal coordinate.
  * `rgrid`: The grid corresponding to the radial coordinate.
  * `table`: The modulation on the grid(s) relative to `medium_pressure_in_bar`.
* **group `pre_ionised`**: Arbitrary pre-ionisation profile can be specified by `initial_electrons_ratio`. The organisation is analogical to the density modulation. The difference is that none of the grids is required. If there is no grid, the pre-ionisation is applied globally.
  * `zgrid`: The grid corresponding to the longitudinal coordinate.
  * `rgrid`: The grid corresponding to the radial coordinate.
  * `initial_electrons_ratio`: The pre-ionisation relative to the *local gas density*.

### CUPRAD
The input parameters of CUPRAD are stored in `CUPRAD/inputs` group. The default input beam and pulse are Gaussian profiles (see the example in [this jupyter notebook](./jupyter_examples/Bessel-Gauss_beams/prepare_Bessel.ipynb) for a customised field profile). Some inputs might be alternated (there are more ways to specify the geometry of the beam or the duration of the pulse, ...). The alternative inputs are stored in the `calculated` subgroup created by the pre-processor.
* **laser**
  * **`laser_wavelength`**: The central wavelength, $\lambda$, of the driving field.
  * **beam geometries and the entry intensity**: (It is obligatory to use one of the sets.)
    * Reference Gaussian beam: One option to specify the geometry of the beam is *the reference Gaussian beam*.
      * `laser_focus_beamwaist_Gaussian`: The beamwaist in the focus.
      * `laser_focus_position_Gaussian`: The position of the focus relative to the entry of the medium.
      * `laser_focus_intensity_Gaussian`: The peak intensity in the focus.
    * Specify the beam directly at the entry plane:
      * `laser_beamwaist_entry`: The beam radius at the entry plane.
      * `laser_focus_position_Gaussian`: The focal point of a virtual lens placed at the entry plane. (It imprints the according curvature in the entry plane.)
      * `laser_intensity_entry` or `laser_energy` or `laser_ratio_pin_pcr`
        * `laser_intensity_entry`: Peak intensity at the entrry plane.
        * `laser_energy`: The total energy in the laser pulse.
        * `laser_ratio_pin_pcr`: The peak intensity is inferred from the critical power $P_{\text{cr}}=\lambda^2/(2\pi n_2(p))$, where $n_2(p)$ is the non-linear refractive index charactersing the Kerr effect, at a given pressure $p$. The relation with the peak intensity $I_0$ is: $P_{\text{in}}/P_{\text{cr}}=n_2(p)I_0 (\pi w(z)/\lambda)^2$, where $P_{\text{in}}=\pi I_0 w^2(z)/2$. (This value is related with the possible beam collapse due to Kerr self-focusing, [see Sec. 3.1 here](https://iopscience.iop.org/article/10.1088/0034-4885/70/10/R03).)
  * **pulse duration specifications**: The pulse duration is specified by **either of these variables**.
    * `laser_pulse_duration_in_1_e_Efield`: The lenght of the pulse measured as the interval where the electric field amplitude exceeds $\mathcal{E}_{\text{max}}/\mathrm{e}$.
    * `laser_pulse_duration_in_1_e_Intensity`: The lenght of the pulse measured as the interval where the intensity exceeds $I_{\text{max}}/\mathrm{e}$.
    * `laser_pulse_duration_in_FWHM_Efield`: The lenght of the pulse measured as the interval where the electric field amplitude exceeds $\mathcal{E}_{\text{max}}/2$.
    * `laser_pulse_duration_in_FWHM_Intensity`: The lenght of the pulse measured as the interval where the intensity exceeds $I_{\text{max}}/2$.
    * `laser_pulse_duration_in_rms_Efield`: The lenght of the pulse measured by $`\tau = \sqrt{\int_{-\infty}^{+\infty}t^2\mathcal{E}_{\text{envelope}}(t)\,\mathrm{d}t/\int_{-\infty}^{+\infty}\mathcal{E}_{\text{envelope}}(t)\,\mathrm{d}t}`$ ([Ref. this discussion about the analogical spatial beam measuremet](https://en.wikipedia.org/w/index.php?title=Beam_diameter&oldid=1226051288#ISO11146_beam_width_for_elliptic_beams).)
    * `laser_pulse_duration_in_rms_Intensity`: Analogical to the previous one, but using hte intensity: $\tau = \sqrt{\int_{-\infty}^{+\infty}t^2 I(t)\,\mathrm{d}t/\int_{-\infty}^{+\infty}I(t)\,\mathrm{d}t}$.
  * `laser_degree_of_supergaussian`: The degree $d$ of the supergaussian anvelope in space $\mathcal{E}(\rho)\propto \mathrm{e}^{-(\rho/\rho_0)^{2d}}$.
  * `laser_degree_of_supergaussian_in_time`: The degree $d$ of the superaguassian anvelope in time $\mathcal{E}_{\text{envelope}}(\rho)\propto \mathrm{e}^{-(t/t_0)^{2d}}$.
  * `laser_initial_chirp_phase`: Initial phase modulation of the laser pulse, known as chirp.
* **medium**
  * **`medium_physical_distance_of_propagation`**: Physical distance over which the laser propagates in the medium.
  * `medium_pressure_in_bar`: Pressure of the medium in bars. (Preferrably defined it in the **global inputs** for the multiscale usage.)
  * `medium_effective_atmospheric_density_of_neutral_molecules`: Effective density of neutral molecules in the medium under atmospheric conditions. See [this reference](https://en.wikipedia.org/wiki/Number_density#Units).
  * `Kerr_nonlinear_refractive_index_kerr_coefficient`: [Coefficient $n_2$ that quantifies the nonlinear change in the refractive index due to the Kerr effect.](https://ieeexplore.ieee.org/document/5412129)
  * `Kerr_ionised_atoms_relative_Kerr_response`: The response of the ions relative to the neutrals, it equals $n_2^{\text{(ions)}}/n_2^{\text{(neutrals)}}$.
  * `Kerr_chi5_coefficient`: The fifth-order nonlinearity coefficient for Kerr effect, indicating the strength of the nonlinear response in a medium.
  * `Kerr_type_of_delayed_kerr_response`: The option to include delayed Kerr effect according to [Section 2.2 here](https://iopscience.iop.org/article/10.1088/0034-4885/70/10/R03). 1 -- instantaneous Kerr only (default); 2 -- instantaneous Kerr + simplified Raman response; 3 -- instantaneous Kerr + full intrapulse Raman response.
  * `Kerr_ratio_of_delayed_kerr_xdk`: Ratio $0\le x_{dK}\le 1$ of delayed Kerr (Raman) contribution; only active for `Kerr_type_of_delayed_kerr_response` = 2 or 3.
  * `Kerr_time_of_delayed_kerr_tdk`: Dipole dephasing time of Raman response in fs; only active for `Kerr_type_of_delayed_kerr_response` = 2 or 3.
  * `Kerr_frequency_in_delayed_kerr_wr`: Raman frequency in 1/fs; only active for `Kerr_type_of_delayed_kerr_response` = 3.
  * `dispersion_type_of_dispersion_law`: (DEPRECATED): Switch to select the dispersion law for neutrals. It is recomended to use the `gas_preset` to select it.
  * `ionization_ionization_potential_of_neutral_molecules`: The ionization potential of neutral molecules, indicating the energy required to ionize a molecule.
  * `ionization_model`: Model used to describe the ionization process. There are two options *PPT* and *ext*. The former computes the ionisation table using the PPT model. The latter reads user-inputted ionisation table from the group `CUPRAD/ionisation_model`.
  * `ionization_effective_residue_charge_for_method_3_4_7`: (DEPRECATED) Effective residual charge left after ionization for methods 3, 4, and 7.
  * `ionization_angular_momentum_for_method_3_7`: (DEPRECATED) Angular momentum of ionization for ionization methods 3 and 7.
  * `ionization_type_of_ionization_method`: (DEPRECATED) Type of ionization method used in the model.
  * `plasma_electron_colision_time`: Average time between collisions for electrons in the plasma to model collisional recombination.
  * `plasma_density_of_absorbing_molecules`: Density of molecules in the plasma that absorb radiation.
  * `plasma_initial_electron_density`: (DEPRECATED) Use the pre-ionisation module (`global_inputs/pre_ionised`) instead.
  * `plasma_linear_recombination_coefficient`: Coefficient for linear recombination processes in the plasma.
  * `plasma_number_of_photons_involved_in_the_n-absorption`: Number of photons involved in the nth order absorption process in the plasma.
  * `plasma_quadratic_recombination_(gasses)`: Coefficient for quadratic recombination processes in gaseous plasma.
  * `plasma_the_n-photon_absorption_cross_section`: Cross-section for the n-photon absorption process in the plasma.
* **numerics**
  * **`numerics_run_time_in_hours`**: Total run time provided to the simulation. Should be slightly smaller than the SLURM limit to include the overhead to finalise the simulation.
  * **`numerics_length_of_window_for_r_normalized_to_beamwaist`***: Length of the numerical window for the radial coordinate normalized to the beam waist.
  * **`numerics_length_of_window_for_t_normalized_to_pulse_duration`**: Length of the numerical window for the time coordinate normalized to the pulse duration ($1/\mathrm{e}$ duration in the electric field is used as the reference).
  * **`numerics_number_of_absorber_points_in_time`**: Number of absorber points at the edges of the time grid.
  * **`numerics_number_of_points_in_r`**: Number of grid points in the radial grid. Required to be a power of 2.
  * **`numerics_number_of_points_in_t`**: Number of grid points in the time grid. Required to be a power of 2.
  * **`numerics_operators_t_t-1`**: 1 -- Standard SVEA; 2 -- Higher order SVEA corrections (operators $T$, $T^{-1}$) used (default).
  * **`numerics_output_distance_in_z-steps_for_fluence_and_power`**: The number of steps in $z$ for storring the fluence and the power of the beam.
  * **`numerics_phase_threshold_for_decreasing_delta_z`**: The maximal phase variations between two consecutive $z$-planes to decrease the stepsize in $z$.
  * **`numerics_physical_first_stepwidth`**: Initial step width in $z$.
  * **`numerics_physical_output_distance_for_plasma_and_Efield`**: Output $z$ distance for storing the electric field and plasma density.
  * **`numerics_radius_for_diagnostics`**: Radius used for diagnostic calculations.
  * `numerics_type_of_input_beam`: (DEPRECATED) Formerly used to manage the input fields.
  * `numerics_noise_on_the_input_shape`: Artificial noise level applied to the input field. (Might be use to test the robustness of the calculation. *It should not be used for fields used for HHG! TDSE input is sensitive.*)
  * `numerics_spatial_noise_on_the_input_shape`: Artificial noise level applied to the input field in the spatial domain only.
  * `numerics_temporal_noise_on_the_input_shape`: Artificial noise level applied to the input field in the time domain only.




### CTDSE
Flags `print_xxx` define whether a given output is stored.
* `CV_criterion_of_GS`: Stopping criterion for convergence of the ground state computation using the resolvent-operation. (The iterations are stopped if $|E_{i+1}-E_i| < CV$.)
* **`dt`**: Time step size for the propagation in time.
* **`dx`**: Spatial resolution of the microcopic grid.
* **`Nx_max`**: Maximum number of grid points in the positive $x$ direction. (The total number of points is 2`Nx_max`+1)
* **`x_int`**: The "size of the atom" to compute the volumetric population of the ground state ([see Section 3.1 here](https://theses.hal.science/tel-04192431v1/document)).
* `InterpByDTorNT`: (DEPRECATED) Option to either refine the discretisation of the numerical input $\mathcal{E}$ to obtain `dt` or use `Ninterp` intermediate points.
* `Ninterp`: (DEPRECATED) Number of intermediate points for hte interpolation of $\mathcal{E}$. (Applied iff `InterpByDTorNT`=1).
* **`Nr_max`**: Maximal index in the macroscopic grid for which TDSE is computed.
* **`kr_step`**: Radial stride of the macroscopic radial grid for computing TDSE.
* **`kz_step`**: Longitudinal stride of the macroscopic radial grid for computing TDSE.
* **`print_Efield`**: Flag to print the electric field.
* **`print_F_Efield`**: Flag to print the Fourier-transformed electric field.
* **`print_F_Efield_M2`**: Flag to print the squared magnitude of the Fourier-transformed electric field, $`|\mathscr{F}[\mathcal{E}](\omega)|^2`$.
* **`print_Source_Term`**: Flag to print the source term $\partial_t j (t)$.
* **`print_F_Source_Term`**: Flag to print the Fourier-transformed source term $`\mathscr{F}[\partial_t j](\omega)`$.
* **`print_F_Source_Term_M2`**: Flag to print the squared magnitude of the Fourier-transformed source term $`|\mathscr{F}[\partial_t j](\omega)|^2`$.
* **`print_GS_population`**: Flag to print the population of the ground state.
* **`print_integrated_population`**: Flag to print the integrated volumetric population over time ($\int_{-x_{\text{int}}}^{x_{\text{int}}} |\psi(x,t)|^2 \, \mathrm{d}x$).
* **`print_x_expectation_value`**: Flag to print the expectation value of position $\braket{x}$.
* **`print_GS`**: Flag to print the ground state wavefunction.
* `absorber_type`: different options for absorbing boundaries (0-none, 2-complex potential, 3-smoothstep)
* `absorber_x_xap`: the spatial extent of the absorber from the boundaries
* `absorber_alpha`: the parameter of the absorber for option 1 $\sim \mathrm{e}^{-\alpha (|x-x_{\text{boundary}}|-x_{\mathrm{CAP}})^2 \mathrm{d}t}$




### Hankel
* **`Harmonic_range`** (2-component array): The spectral range in the spectral domain, relative to the fundamental frequency defined by the `laser_wavelength`.
* **`ko_step`**: The stride in the provided omega grid (this grid is inherited from TDSE).
* **`Nr_max`**: The maximal index of the provided radial grid in the interaction volume used for the integration.
* **`kr_step`**: The stride along the radial grid for the radial intagration.
* **`Nr_FF`**: "The number of pixels on the far-field detector."
* **`rmax_FF`**: The maximal far-field radial coordinate.
* **`distance_FF`**: The distance of the XUV detector (measured from the entry of the interaction volume).
* **`XUV_table_type_dispersion`**: The tables in the XUV range used for the dispersion ([`NIST`](https://physics.nist.gov/PhysRefData/FFast/html/form.html) and [`Henke`](https://henke.lbl.gov/optical_constants/asf.html) are available in the code.)
* **`XUV_table_type_dispersion`**: The tables in the XUV range used for the absorption ([`NIST`](https://physics.nist.gov/PhysRefData/FFast/html/form.html) and [`Henke`](https://henke.lbl.gov/optical_constants/asf.html) are available in the code.)
* **`store_cumulative_result`**: Option to keep the cumulative integral along $z$.
* **`Nthreads`**: The number of threads used by the multiprocessing.


