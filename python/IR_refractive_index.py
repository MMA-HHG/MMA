####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021)
#
# This module creates the interpolation functions from tables loaded from the archive
# The functions are called by f1, f2 = XUV_refractive_index.getf(gas,Energy), Energy is in eV
# alternatively, the full handle is XUV_refractive_index.index_funct[gas][f1/f2](Energy), Energy is in eV


import numpy as np
from scipy import interpolate
import h5py
import sys
import os
import shutil


susc_Ar = lambda x: (5.547e-4)*(1.0 + (5.15e5)/(x**2) + (4.19e11)/(x**4) + (4.09e17)/(x**6) + (4.32e23)/(x**8))
# susc_Ar = lambda x: susc_Ar(1e10*x)

susc_Kr = lambda x: (8.377e-4)*(1.0 + (6.7e5)/(x**2) + (8.84e11)/(x**4) + (1.49e18)/(x**6) + (2.74e24)/(x**8) + (5.1e30)/(x**10))
# susc_Kr = lambda x: susc_Kr(1e10*x)

susc_funct = {
    'Ar': susc_Ar,
    'Kr': susc_Kr
}

def getsusc(g,lambd):
    return susc_funct[g](1e10*lambd)