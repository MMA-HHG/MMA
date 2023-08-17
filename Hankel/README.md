# Development
The main API is the `HFn2.py`-module. There shall be some good examples of calling it. THe cannonical operation is pipelined after 1D-TDSE, it processes `*.hdf5` dipoles and compute their diffraction integral. Beside this, there is a few other uses of this code (*Maker fringes*, $\omega + 2\omega$ HHG, ...) The preparation shall consider following ideas:

## `Hfn2.py`
* Actual inplementation assumes all the dipoles are in the input array. This might be too memory-demanding. Maybe we can create a *class* providing the dipoles. This class could internally treat storing them or providing them on-the-fly. It would need some thinking about the efficiency.

## Calling the `Hfn2.py`
* The actual example for large-scale applications uses `multiprocessing`, i.e. parallelisation limited to multithreading. THe performace is thus limited by the # of threads per core... We can think about using `mpi4py`, it would require to develop a data treatment as they cannot be easily managed in MPI.

## General remarks
* The output is now in *arbitrary units*. Principally, all physical constants are included and we shall be able to retrieve true XUV intensity.[^1]
* Clear junk files. (This directory contains many junk files for testing etc.)


[^1]: This would require some testing and verifying. So far, we've been using *arbitrary units* for all comparisons with experiments etc. 