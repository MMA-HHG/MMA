# MMA for HHG modelling

[link](1DTDSE/README.md)

[link2](jupyter_examples/density_profile/prepare_inputs.ipynb)

Table of content

## Installations
Both `CUPRAD` and `CTDSE` are compiled from the source. Hankel is implemented in Python and together with the other Pythonic procedures and scripting around the model rely on setting up the environment variables (one might consider setting up [a virtual environment](https://stackoverflow.com/questions/9554087/setting-an-environment-variable-in-virtualenv) dedicated to this whole project).
<!--- We head for RC1, use virtual environment, pip, etc. in future releases by default --->

We provide two ways to obtain the code. The first option uses the Docker image. It is a direct multiplatform user-oriented way to obtain the executable model. This can be used for running the model locally. Moreover, this can be used as a direct reference for compiling the code the second way: directly from the source. This option might be neccessary for deploying the code on HPC clusters, develpoment, …  

## Reference Docker installation
The code can be accessed through Docker. We provide the direct Docker image, CodeOcean capsule and here we show how to build the code using Docker.

## Custom installations
Here we provide a more detailed guide for the installation, this is the list of requirements:
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
export MSM_PATH=$GIT_PATH/CUPRAD_TDSE_Hankel
export PYTHONPATH=$PYTHONPATH:$MSM_PATH/shared_python

export CUPRAD_HOME=$MSM_PATH/CUPRAD
export CUPRAD_BUILD=$MSM_PATH/CUPRAD/binary
export CUPRAD_SCRIPTS=$MSM_PATH/CUPRAD/scripts
export CUPRAD_PYTHON=$MSM_PATH/CUPRAD/python
export PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON

export TDSE_1D_HOME=$MSM_PATH/1DTDSE
export TDSE_1D_POST_PROCESSING=$MSM_PATH/1DTDSE/post_processing
export TDSE_1D_SCRIPTS=$MSM_PATH/1DTDSE/scripts
export TDSE_1D_SLURM=$MSM_PATH/1DTDSE/slurm
export TDSE_1D_BUILD=$MSM_PATH/1DTDSE/build
export PYTHONPATH=$PYTHONPATH:$TDSE_1D_HOME

export HANKEL_HOME=$MSM_PATH/Hankel

export MULTISCALE_HOME=$MSM_PATH
export MULTISCALE_SCRIPTS=$MSM_PATH/multiscale/scripts

export FSPA_PATH=$MSM_PATH/FSPA

export MULTISCALE_WORK_DIR=/mnt/d/data/work_dir
```


#### Modules and libraries
When used locally on a personal computer, the libraries (FFTW3, CMake, …) are typically installed by a user. If using the code or a part of it, the corresponding libraries must be installed before. (Sharing needs to be enebled for FFTW3, see details below for the dyamic-library CTDSE.)

[The modules](https://hpc-wiki.info/hpc/Modules) provide all the necessary libraries for the code when using a computational cluster. The script `Modules/load_modules.sh` is used to load all the modules. There is a list of modules for various computational clusters specified by the variable `$HPC`. Another supercomputer (or compilation option intel/GNU/...) should be added there.

There are two `bash` functions `load_modules` and `load_python_modules`. The former is activated when running *CUPRAD* and *CTDSE*, while the latter is used for all Pythonic operations around the code. (The reason for this duality is that Python might need to load a compiler itself for some libraries, typically on intel.) 

### Installing both CUPRAD and CTDSE
If everything is set well, the following CMakes are wrapped in the master `CMakeList.txt`. Here is the recipe to install the code from its root directory.

1) Run `load_modules`. [This can be verified by](https://hpc-wiki.info/hpc/Modules#:~:text=%24-,module%20list,-Currently%20Loaded%20Modulefiles) `module list`.
    * If the machine does not using modules, this step is replaced by installing the necessary libraries and setting up the environment.
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory. We encountered CMake sometimes struggling to identify the proper MPI-Fortran compiler on several machines (the error is raised in the next step). There are more ways to hint CMake to find the compilers:
    * By providing environment variables with the compilers: `export CC=mpicc` and `export FC=mpifort` (GNU); `export CC=mpiicc` and `export FC=mpiifort` (Intel).
    * Controlling CMake directly during its execution `cmake -D CMAKE_Fortran_COMPILER=mpifort ..` (GNU) or `cmake -D CMAKE_Fortran_COMPILER=mpiifort ..` (intel). This resolved the issue when we encountered it.
    * The CMake configuration can be manually adjusted using `ccmake`, see [link 1](https://cmake.org/cmake/help/latest/manual/ccmake.1.html) and [link 2](https://stackoverflow.com/a/1224652).
    * Consider to run the compilation of the respective codes separately in the case problems occur.
3) Compile the code from the CMake-generated `Makefile` by running `make code` in the root directory.


Below are recipes for compilling just one of the codes.


### CUPRAD
All the source files are located in `CUPRAD/sources`. The CMake recipe is in `CUPRAD/CMakeLists.txt`. The code is supposedly built in `CUPRAD/build`.
There is the recipe for compilation (each point contains several notes about possible difficulties):

1) Run `load_modules`. [This can be verified by](https://hpc-wiki.info/hpc/Modules#:~:text=%24-,module%20list,-Currently%20Loaded%20Modulefiles) `module list`.
    * If the machine does not using modules, this step is replaced by installing the necessary libraries and setting up the environment.
2) Prepare Makefile using `cmake` by running `cmake ..` in the `build` directory.
    * We encountered CMake struggling to identify the proper MPI-Fortran compiler on several machines. CMake can be hinted to use the desired compiler by `cmake -D CMAKE_Fortran_COMPILER=mpifort ..` (GNU) or `cmake -D CMAKE_Fortran_COMPILER=mpiifort ..` (intel). This resolved the issue when we encountered it.
    * The CMake configuration can be manually adjusted using `ccmake`, see [link 1](https://cmake.org/cmake/help/latest/manual/ccmake.1.html) and [link 2](https://stackoverflow.com/a/1224652).
3) Compile the code from the CMake-generated `Makefile` by running `make code` in the `build` directory.


#### Note: building only the pre-processor
The pre-processor can built without MPI and parallel HDF5 (serial HDF5 is needed). The compilation is similar except using **`make preprocessor`** instead of `make code`.

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
Here is the exhaustive list of all the parameters. The bold **`parameters`** are obligatory to run the whole model with sourcing the default material constants. The other `parameters` are optional. If an optional parameter is present, it has priority. 

### Global
The global inputs are stored in the `global_inputs` groups, they might be used by more than one module.
* `gas_preset`: The main specifier to define the gas. It defines all the material constants then. Implemented gases: `He`, `Ne`, `Ar`, `Kr`, `Xe`.
* **`medium_pressure_in_bar`**: (DEPRECATED) Pressure of the medium in bars. (Preferrably defined it in the **global inputs** for the multiscale usage.)
* **group `density_mod`**: This group defines the density modulation. The relative modulation to `medium_pressure_in_bar` is stored in `table` with the dimension corresponding to the grids (it can be 1- or 2-dimensional). Therefore, `zgrid` and/or `rgrid` needs to be present. The grids should be larger than the interaction volume. Note: the radial modulation is not fully available for Hankel (see detailed documentation of the module for details).
  * `zgrid`: The grid corresponding to the longitudinal coordinate.
  * `rgrid`: The grid corresponding to the radial coordinate.
  * `table`: The modulation on the grid(s) relative to `medium_pressure_in_bar`.
* **group `pre_ionised`**: Arbitrary pre-ionisation profile can be specified by `initial_electrons_ratio`. The organisation is analogical to the density modulation. The difference is that none of the grids is required. If there is no grid, the pre-ionisation is applied globally.
  * `zgrid`: The grid corresponding to the longitudinal coordinate.
  * `rgrid`: The grid corresponding to the radial coordinate.
  * `initial_electrons_ratio`: The pre-ionisation relative to the *local gas density*.

### CUPRAD
The input parameters of CUPRAD are stored in `CUPRAD/inputs` group. The default input beam and pulse are Gaussian profiles (see the example in <span style="color:red">JUPYTER NOTEBOOK</span> for a customised field profile). Some inputs might be alternated (there are more ways to specify the geometry of the beam or the duration of the pulse, ...). The alternative inputs are stored in the `calculated` subgroup created by the pre-processor.
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
    * `laser_pulse_duration_in_rms_Efield`: The lenght of the pulse measured by $\tau = \sqrt{\int_{-\infty}^{+\infty}t^2\mathcal{E}_{\text{envelope}}(t)\,\mathrm{d}t/\int_{-\infty}^{+\infty}\mathcal{E}_{\text{envelope}}(t)\,\mathrm{d}t}$ ([Ref. this discussion about the analogical spatial beam measuremet](https://en.wikipedia.org/w/index.php?title=Beam_diameter&oldid=1226051288#ISO11146_beam_width_for_elliptic_beams).)
    * `laser_pulse_duration_in_rms_Intensity`: Analogical to the previous one, but using hte intensity: $\tau = \sqrt{\int_{-\infty}^{+\infty}t^2 I(t)\,\mathrm{d}t/\int_{-\infty}^{+\infty}I(t)\,\mathrm{d}t}$.
  * `laser_degree_of_supergaussian`: The degree $d$ of the supergaussian anvelope in space $\mathcal{E}(\rho)\propto \mathrm{e}^{-(\rho/\rho_0)^{2d}}$.
  * `laser_degree_of_supergaussian_in_time`: The degree $d$ of the superaguassian anvelope in time $\mathcal{E}_{\text{envelope}}(\rho)\propto \mathrm{e}^{-(t/t_0)^{2d}}$.
  * `laser_initial_chirp_phase`: Initial phase modulation of the laser pulse, known as chirp.
* **medium**
  * **`medium_physical_distance_of_propagation`**: Physical distance over which the laser propagates in the medium.
  * `medium_pressure_in_bar`: (DEPRECATED) Pressure of the medium in bars. (Preferrably defined it in the **global inputs** for the multiscale usage.)
  * `medium_effective_atmospheric_density_of_neutral_molecules`: Effective density of neutral molecules in the medium under atmospheric conditions. See [this reference](https://en.wikipedia.org/wiki/Number_density#Units).
  * `Kerr_nonlinear_refractive_index_kerr_coefficient`: [Coefficient $n_2$ that quantifies the nonlinear change in the refractive index due to the Kerr effect.](https://ieeexplore.ieee.org/document/5412129)
  * `Kerr_ionised_atoms_relative_Kerr_response`: The response of the ions relative to the neutrals, it equals $n_2^{\text{(ions)}}/n_2^{\text{(neutrals)}}$.
  * `Kerr_chi5_coefficient`: The fifth-order nonlinearity coefficient for Kerr effect, indicating the strength of the nonlinear response in a medium.
  * `Kerr_type_of_delayed_kerr_response`: The option to include delayed Kerr effect according to [Section 2.2 here](https://iopscience.iop.org/article/10.1088/0034-4885/70/10/R03). Turned of in the default mode.
  * `dispersion_type_of_dispersion_law`: (DEPRECATED): Switch to select the dispersion law for neutrals. It is recomended to use the `gas_preset` to select it.
  * `ionization_ionization_potential_of_neutral_molecules`: The ionization potential of neutral molecules, indicating the energy required to ionize a molecule.
  * `ionization_model`: Model used to describe the ionization process. There are two options *PPT* and *ext*. The former computes the ionisation table using the PPT model. The latter reads user-inputted ionisation table from the group `CUPRAD/ionisation_model`.
  * `ionization_effective_residue_charge_for_method_3_4_7`: (DEPRECATED) Effective residual charge left after ionization for methods 3, 4, and 7.
  * `ionization_angular_momentum_for_method_3_7`: (DEPRECATED) Angular momentum of ionization for ionization methods 3 and 7.
  * `ionization_type_of_ionization_method`: (DEPRECATED) Type of ionization method used in the model.
  * `plasma_electron_colision_time`: Average time between collisions for electrons in the plasma to model collisional recmbination.
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
  * **`numerics_operators_t_t-1`**: <span style="color:red">Operators used for advancing the solution from one time step to the next. Stefan: Refer to the description in the paper.</span>
  * **`numerics_output_distance_in_z-steps_for_fluence_and_power`**: The number of steps in $z$ for storring the fluence and the power of the beam.
  * **`numerics_phase_threshold_for_decreasing_delta_z`**: The maximal phase variations between two consecutive $z$-planes to decrease the stepsize in $z$.
  * **`numerics_physical_first_stepwidth`**: Initial step width in $z$.
  * **`numerics_physical_output_distance_for_plasma_and_Efield`**: Output $z$ distance for storing the electric field and plasam density.
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
* **`print_F_Efield_M2`**: Flag to print the squared magnitude of the Fourier-transformed electric field, $|\mathscr{F}[\mathcal{E}](\omega)|^2$.
* **`print_Source_Term`**: Flag to print the source term $\partial_t j (t)$.
* **`print_F_Source_Term`**: Flag to print the Fourier-transformed source term $\mathscr{F}[\partial_t j](\omega)$.
* **`print_F_Source_Term_M2`**: Flag to print the squared magnitude of the Fourier-transformed source term $|\mathscr{F}[\partial_t j](\omega)|^2$.
* **`print_GS_population`**: Flag to print the population of the ground state.
* **`print_integrated_population`**: Flag to print the integrated volumetric population over time ($\int_{-x_{\text{int}}}^{x_{\text{int}}} |\psi(x,t)|^2 \, \mathrm{d}x$).
* **`print_x_expectation_value`**: Flag to print the expectation value of position $\braket{x}$.
* **`print_GS`**: Flag to print the ground state wavefunction.



### Hankel

## Execution pipeline


# Standalone operations




