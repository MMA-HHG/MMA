# CUPRAD


## The execution of the code
The executions of the code is split into two steps:
* First, the input is pre-processed by `/build/make_start.e`, the pre-processor defines the computational units based on the input parameters and constructs an input Gaussian field profile on the computational grids.
* The main parallel program `/build/cuprad.e` computes the non-linear pulse propagation.

## Loading data to Python & visualisation
THe inputs and outputs are organised within the hdf5-archive
