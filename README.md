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
Here is the exhaustive list of all the parameters. The bold **`parameters`** are obligatory to run the whole model with sourcing the default material constants. The other `parameters` are optional. If an optional parameter is present, it has priority. Some inputs might be alternated (there are more ways to specify the geometry of the beam or the duration of the pulse, ...). The alternative inputs are stored in the `calculated` subgroup created by the pre-processor.
* **laser group**
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
        * `laser_ratio_pin_pcr`: The peak intensity is inferred from the critical power $P_{\text{cr}}=\lambda^2/(2\pi n_2(p))$, where $n_2(p)$ is the non-linear refractive index charactersing the Kerr effect, at a given pressure $p$. The relation with the peak intensity $I_0$ is: $P_{\text{in}}/P_{\text{cr}}=n_2(p)I_0 (\pi w(z)/\lambda)^2$, where $P_{\text{in}}=I_0w^2(z)\pi /2$. (This value is related with the possible beam collapse due to Kerr self-focusing, [see Sec. 3.1 here](https://iopscience.iop.org/article/10.1088/0034-4885/70/10/R03).)
  * **pulse duration specifications**: The pulse duration is specified by either of these variables.
    * `laser_pulse_duration_in_1_e_Efield`: Duration of the laser pulse measured at the 1/e point of the electric field amplitude.
  * `laser_degree_of_supergaussian`: The degree $d$ of the supergaussian anvelope in space $\mathcal{E}(\rho)\propto \mathrm{e}^{-(\rho/\rho_0)^{2d}}$.
  * `laser_degree_of_supergaussian_in_time`: The degree $d$ of the superaguassian anvelope in time $\mathcal{E}_{\text{envelope}}(\rho)\propto \mathrm{e}^{-(t/t_0)^{2d}}$.
  * `laser_initial_chirp_phase`: Initial phase modulation of the laser pulse, known as chirp.
* **medium group**
  * `medium_effective_atmospheric_density_of_neutral_molecules`: Effective density of neutral molecules in the medium under atmospheric conditions.
  * **`medium_physical_distance_of_propagation`**: Physical distance over which the laser propagates in the medium.
  * **`medium_pressure_in_bar`**: Pressure of the medium in bars.
  * `Kerr_chi5_coefficient`: The fifth-order nonlinearity coefficient for Kerr effect, indicating the strength of the nonlinear response in a medium.
  * `Kerr_ionised_atoms_relative_Kerr_response`: The change in the Kerr effect response due to the presence of ionized atoms.
  * `Kerr_nonlinear_refractive_index_kerr_coefficient`: Coefficient that quantifies the nonlinear change in the refractive index due to the Kerr effect.
  * `Kerr_type_of_delayed_kerr_response`: Type of temporal delay in the Kerr response of a material.
  * `dispersion_type_of_dispersion_law`: Type of dispersion law used to describe the frequency dependence of the refractive index.
  * `ionization_angular_momentum_for_method_3_7`: Angular momentum of ionization for ionization methods 3 and 7.
  * `ionization_effective_residue_charge_for_method_3_4_7`: Effective residual charge left after ionization for methods 3, 4, and 7.
  * `ionization_ionization_potential_of_neutral_molecules`: The ionization potential of neutral molecules, indicating the energy required to ionize a molecule.
  * `ionization_model`: Model used to describe the ionization process.
  * `ionization_type_of_ionization_method`: Type of ionization method used in the model.
  * `plasma_density_of_absorbing_molecules`: Density of molecules in the plasma that absorb radiation.
  * `plasma_electron_colision_time`: Average time between collisions for electrons in the plasma.
  * `plasma_initial_electron_density`: Initial density of electrons in the plasma.
  * `plasma_linear_recombination_coefficient`: Coefficient for linear recombination processes in the plasma.
  * `plasma_number_of_photons_involved_in_the_n-absorption`: Number of photons involved in the nth order absorption process in the plasma.
  * `plasma_quadratic_recombination_(gasses)`: Coefficient for quadratic recombination processes in gaseous plasma.
  * `plasma_the_n-photon_absorption_cross_section`: Cross-section for the n-photon absorption process in the plasma.
* **numerics group**
  * `numerics_length_of_window_for_r_normalized_to_beamwaist`: Length of the numerical window for the radial coordinate normalized to the beam waist.
  * `numerics_length_of_window_for_t_normalized_to_pulse_duration`: Length of the numerical window for the time coordinate normalized to the pulse duration.
  * `numerics_noise_on_the_input_shape`: Noise level applied to the input shape in numerical simulations.
  * `numerics_number_of_absorber_points_in_time`: Number of absorber points used in the time domain for numerical simulations.
  * `numerics_number_of_points_in_r`: Number of grid points in the radial coordinate for numerical simulations.
  * `numerics_number_of_points_in_t`: Number of grid points in the time coordinate for numerical simulations.
  * `numerics_operators_t_t-1`: Operators used for advancing the solution from one time step to the next in numerical simulations.
  * `numerics_output_distance_in_z-steps_for_fluence_and_power`: Output distance in the z direction for fluence and power calculations in numerical simulations.
  * `numerics_phase_threshold_for_decreasing_delta_z`: Threshold for phase change used to decrease the step size in the z direction in numerical simulations.
  * `numerics_physical_first_stepwidth`: Initial step width in physical units for numerical simulations.
  * `numerics_physical_output_distance_for_plasma_and_Efield`: Output distance in physical units for plasma and electric field calculations.
  * `numerics_radius_for_diagnostics`: Radius used for diagnostic calculations in numerical simulations.
  * `numerics_run_time_in_hours`: Total run time of the numerical simulation in hours.
  * `numerics_spatial_noise_on_the_input_shape`: Spatial noise level applied to the input shape in numerical simulations.
  * `numerics_temporal_noise_on_the_input_shape`: Temporal noise level applied to the input shape in numerical simulations.
  * `numerics_type_of_input_beam`: Type of input beam used in numerical simulations.




### CTDSE

### Hankel

## Execution pipeline


# Standalone operations




