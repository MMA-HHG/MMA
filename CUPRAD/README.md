# CUPRAD


## The execution of the code
The executions of the code is split into two steps:
* First, the input is pre-processed by `/build/make_start.e`, the pre-processor defines the computational units based on the input parameters and constructs an input Gaussian field profile on the computational grids.
* The main parallel program `/build/cuprad.e` computes the non-linear pulse propagation.

## Loading data to Python & visualisation
The inputs and outputs are organised within the hdf5-archive. These data can be fetched into a Python class using the module `python/dataformat_CUPRAD.py`. This class then encapsulates all the data. By default, the field in the whole medium together with scalar and small data are available. Optionally, plasma density and other quantities can be loaded. Additionally to the data, the class contains several methods such as adjustments of the reference frame, the field complexification, etc.
