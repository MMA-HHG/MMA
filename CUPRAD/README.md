# CUPRAD

## TODO
* To be removed
    * $\chi_5$-processes
    * $T$-operator case 4
    * (?) Raman remove option 2
* doxygen4FOrtran (useful links https://fortran-lang.discourse.group/t/fortran-2018-keywords/2555, https://github.com/cdslaborg/paramonte, https://cdslaborg.github.io/paramonte-kernel-doc/html)

### Documentation
Doxygen done for
 * `array_helper.f90` `density_mod.f90` `pre_ionised.f90` `linked_list_module.f90` `constants.f90`


 Missing: `calc_start.f90` `code_ionisation.f90` `cuprad.f90` `default_inputs.f90` `fft_rk.f90` `finalize.f90` `firststep.f90` `hdf5_*.f90` `longstep_rf.f90` `make_start.f90` `modules.f90` `normalisation.f90` `output.f90` `write_*.f90`

 ? Contain pre-processor in the documentation, or make an independent documetnaiton for it ?

## Development notes

* The code is not written in a single style according to various contributors.
    * I propose to keep a single style within each module, but to not unify it globally.
    * The in-line commenting of the code shall be extended.
* The use of the code requires a pre-processor.
    * The pre-processor does a lot of tedious work: mainly the conversion from SI inputs to *computational units **[C.U.]***. It will require a lot of thesting and work to Pythonise the pre-processor/include conversion in the main code/... I propse to clean up the pre-processor for now.
* Data treatment:
    * The storing procedures need a better organisation. For example, the field is printed twice.
    * There is a lot of redundant files produced form CUPRAD. We need to ensure that everything is in the HDF5-ouput. Other files will be removed.
* Overall functionality
    * There are still some extensions of the code that shall be either commented or deprecated[^1] to provide a clear code.
* Code continuation: it is now disallowed. THe purpose of it is to separate a long simulation into more jobs. It shall be easy to reintroduce: we just need to properly use the endplane as the input for the next run.

[^1]: It might be reintroduced later.


<!-- `module pre_ionised`: it allows to compute the pre-ionisation. It also encapsulates most of the work in the module and minimal changes are in the main code. Only firs-step changed slightly and on-the-fly calculations.

switch dispersion is applied only in the pre-processor to create the table -->



<!-- ## Development
Transformation of $1/e$ to FWHM has to be checked. 

Should we keep all inputs in one hdf5-group, or use further hierarchisation?

I commented most of the code. I found that the less readable are firstep and longstep_rk. It would be nice to rewiritte using helping modules or at least organise.

`SUBROUTINE mult_phase`: it's not so clear to me what effects are included in the propagator and which are applied in time domain afterwards.

### printing procedure
There are two procedures now,it should be rewritten and driven by an if-construct to select prints on demand. Solve naming problem by passing struct containing strings. -->

# Stub of documentation
This is a nutshell guide of the propagation code, it explains the basics of its operation and main ideas implemented in the code.
## INPUT/OUTPUT

The I/O are realised through a single HDF5-archive, currently called `results.h5`. One archive with all the inputs (except pre-computed tables) is used and all the ouputs of the code are also stored in this archive.

The input archive can be easily prepared from a text list using a general-purpose procedure for preparing the input files: https://github.com/vabekjan/universal_input. The code then works in two stages:

1. The pre-processor prepares inputs for the main code, the reason is that the code uses initial electric field, that has to be prepared according to the dimensions of grids. A second thing in the pre-processor is the conversion in the computational units (see later) of all the parameters. The starting fields are, however, stored in SI units. The reason is that we would like to simplify portability and allow to use any fields of user's interest.

2. The main code that computes the propagates the electric pulse in the medium. Except the field, there are various characteristics stored in the archive. The outputs are converted always to SI units (DELETE THIS WHEN CHECKED), again to simplify portability and reduce the need of a post-processor.

### Pre-computed tables
There are some other possible inputs from `calculated_tables.h5`. There are actually possible pre-computed ionisation tables (ionisation rate as a function of the electric-field amplitude).
Next quantity to be probably here are `Indexes` (to be tested, maybe include them in usual inputs using H5Lexist). This is an array of the variable index of refraction in the medium.

### Character of outputs

The main output are the electric fields and plasma density in time domain at given points in the medium. There are also some further characteristics (as peak power) computed on-the-fly and available in the ouputs. The electric fields are allowed to be printed with two different spacings at the same time. The reason is that another numerical code linked with cuprad may require a finer grid. This is currently memory-consuming, we eventually keep the data twice, since we may repack them once the second code is finished. Finally, there are also some reference ouputs as the ionisation table used for the calculation (now available only for PPT and external table).

### Linked-list buffering
There are some data stored in every step and the total amount is unknown at the beginning due to the adaptive step-size. These data are then buffered in the memory until the code finishes and data are then saved in a block. This approach is general, the limitation is that the data kept in RAM are sufficiently small.


## Key ideas of the code
The "core" of the code is in `cuprad.f90` and `longstep_rk.f90`. First contains the main program and executes the main loop. Second contains the detailed implementation of advancement in single steps.

The other modules used in the main program contain supplementary procedures (preparatory calculations-`firststep.f90`-, outputs processing and ionisation calculation).

### Used units

The code uses "computational units" (C.U.). It ensures to adapt the grids and other parameters for given paramters of the code, as focusing, maximal expected power of the field etc. The main conversion is done in the pre-processor. The code then works with that units. Conversions in the code are done for printing the outputs and for additional inputs as the ionisation model.

### Pre-porocessor

### List of f90-files
