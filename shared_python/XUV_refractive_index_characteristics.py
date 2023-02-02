####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2023)


import numpy as np
from scipy import interpolate
import h5py
import sys
import os
import shutil

import XUV_refractive_index as XUV_index
import units
import mynumerics as mn

def L_abs(gas,omegaSI,pressure, N_atm = 2.7e19*1e6):
    f2_value = XUV_index.getf2(gas, mn.ConvertPhoton(omegaSI, 'omegaSI', 'eV'))
    lambdaXUV    = mn.ConvertPhoton(omegaSI, 'omegaSI', 'lambdaSI')    
    return 1.0 / (2.0 * pressure * N_atm * units.r_electron_classical * lambdaXUV * f2_value) 