# CUPRAD

This is a nutshell guide of the propagation code, it explains the basics of its operation and main ideas implemented in the code.

## Development
I commented most of the code. I found that the less readeble are firstep and longstep_rk. It would be nice to rewiritte using helping modules or at least organise.

`SUBROUTINE mult_phase`: it's not so clear to me what effects are included in the propagator and which are applied in time domain afterwards.

## INPUT/OUTPUT

The I/O are realised through a single HDF5-archive, currently called `results.h5`. One archive with all the inputs (except pre-computed tables) is used and all the ouputs of the code are also stored in this archive.

The input archive can be easily prepared from a text list using a general-purpose procedure for preparing the input files: https://github.com/vabekjan/universal_input. The code then works in two stages:

1) The pre-processor prepares inputs for the main code, the reason is that the code uses initial electric field, that has to be prepared according to the dimensions of grids. A second thing in the pre-processor is the conversion in the computational units (see later) of all the parameters. The starting fields are, however, stored in SI units. The reason is that we would like to simplify portability and allow to use any fields of user's interest.

2) The main code that computes the propagates the electric pulse in the medium. Except the field, there are various characteristics stored in the archive. The outputs are converted always to SI units (DELETE THIS WHEN CHECKED), again to simplify portability and reduce the need of a post-processor.

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
