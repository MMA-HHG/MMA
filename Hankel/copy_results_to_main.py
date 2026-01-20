"""
This module contains collection of classes and functions for computing Hankel transform
(the diffraction integral to obtain far-field signal) incorporating also the longitudinal
integration accounting for phase-matching. THe content is the following:
- FSource_provider: a class transforming heterogenous input-streams into the form suitable for Hankel_long
- get_propagation_pre_factor_function: this function obtains the prefactor for the longitudinal integration
- HankelTransform: The core routine performing the Hankel transform from a single plane
- Hankel_long: The main class of this module providing Hankel transform if the longitudinaly integrated signal

@author: Jan Vábek
"""

import h5py
import MMA_administration as MMA

with open('msg.tmp','r') as msg_file:
    results_file = msg_file.readline()[:-1] # need to strip the last character due to Fortran msg.tmp

 
with h5py.File('results_Hankel.h5', 'r') as f_src, h5py.File(results_file, 'a') as f_dest:
    f_src.copy(f_src[MMA.paths['Hankel_outputs']],f_dest[MMA.paths['Hankel']])

print("Done")
