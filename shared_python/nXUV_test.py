import numpy as np
import os
import time
import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn


# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index



gas_type = 'Kr'
XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
XUV_table_type_absorption = 'Henke' # {Henke, NIST} 
apply_diffraction = ['dispersion', 'absorption']
    


    





omega0 = mn.ConvertPhoton(800e-9,'lambdaSI','omegaSI')
N_atm = 2.7e19 * 1e6
   


# Here are the fuction to obtain the phase factors in SI units: exp(i*omega*function(omega))
def f1_funct(E):
    return XUV_index.getf1(gas_type+'_' + XUV_table_type_diffraction, E)
def f2_funct(E):
    return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)


def dispersion_function_def(omega,N_dens=N_atm):
    f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return (-N_dens*units.r_electron_classical * ((lambdaSI**2)*f1_value/(2.0*np.pi)))      



def absorption_function_def(omega,N_dens=N_atm):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    beta_factor = N_dens*units.r_electron_classical * \
                  ((lambdaSI**2)*f2_value/(2.0*np.pi))
    return beta_factor / units.c_light

def L_abs(omega,N_dens=N_atm):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 1.0 / (2.0 * N_dens * units.r_electron_classical * lambdaSI * f2_value) 


def nXUV(omega,N_dens=N_atm):
    return 1.0 + dispersion_function_def(omega,N_dens) + 1j*absorption_function_def(omega,N_dens)



print(nXUV(17*omega0, 0.020*N_atm))
print(L_abs(17*omega0, 0.020*N_atm))

v_XUV = units.c_light/nXUV(17*omega0, 0.020*N_atm).real

print(0.015*((units.c_light/v_XUV)-1))
print(mn.ConvertPhoton(17*omega0, 'omegaSI', 'lambdaSI'))




