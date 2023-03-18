####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021)
#
# This module creates the interpolation functions from tables loaded from the archive
# The functions are called by f1, f2 = XUV_refractive_index.getf(gas,Energy), Energy is in eV
# alternatively, the full handle is XUV_refractive_index.index_funct[gas][f1/f2](Energy), Energy is in eV


import numpy as np
from scipy import interpolate
import h5py
import os
import units
import mynumerics as mn

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


source_archive = os.path.join(THIS_DIR, 'XUV_refractive_index_tables.h5')
index_funct = {}

with h5py.File(source_archive, 'r') as SourceFile: # access option http://docs.h5py.org/en/stable/high/file.html#file
    gases = list(SourceFile.keys())
    index_table = {}
    print(gases)
    for gas in gases:
        local_table = {
            'Energy_f1': SourceFile[gas]['Energy_f1'][:],
            'Energy_f2': SourceFile[gas]['Energy_f2'][:],
            'f1': SourceFile[gas]['f1'][:],
            'f2': SourceFile[gas]['f2'][:]
        }
        index_table.update({gas: local_table})
        local_table = {
            'f1': interpolate.interp1d(SourceFile[gas]['Energy_f1'][:], SourceFile[gas]['f1'][:]),
            'f2': interpolate.interp1d(SourceFile[gas]['Energy_f2'][:], SourceFile[gas]['f2'][:])
        }
        index_funct.update({gas: local_table})

def getf(g,E):
    return index_funct[g]['f1'](E)[()], index_funct[g]['f2'](E)[()]

def getf1(g,E):
    return index_funct[g]['f1'](E)[()]

def getf2(g,E):
    return index_funct[g]['f2'](E)[()]



N_atm_default = 2.7e19*1e6

def dispersion_function(omega, pressure, gas, n_IR=1., N_atm=N_atm_default):
    # f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
    f1_value = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    nXUV     = 1.0 - pressure*N_atm*units.r_electron_classical * ((lambdaSI**2)*f1_value/(2.0*np.pi))           
    phase_velocity_XUV  = units.c_light / nXUV    
    phase_velocity_IR = units.c_light / n_IR
    return ((1./phase_velocity_IR) - (1./phase_velocity_XUV))

def beta_factor_atm(omega, gas, N_atm=N_atm_default):
    # omegaXUV    = Horder_foo*omegaSI
    # f2_value    = f2_funct(mn.ConvertPhoton(omegaXUV, 'omegaSI', 'eV'))
    f2_value    = getf2(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaXUV    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    beta_factor = N_atm*units.r_electron_classical * \
                  ((lambdaXUV**2)*f2_value/(2.0*np.pi))
    return beta_factor

def L_abs(omega, pressure, gas, N_atm=N_atm_default):
    f2_value    = getf2(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaXUV   = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 1.0 / (2.0 * pressure * N_atm * units.r_electron_classical * lambdaXUV * f2_value) 

def susc_atm(omega, gas, N_atm=N_atm_default):
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    return nXUV_atm**2 - 1

def polarisability(omega, gas, N_atm=N_atm_default):
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    susc_XUV_atm = nXUV_atm**2 - 1
    pol_XUV = susc_XUV_atm/N_atm
    return pol_XUV