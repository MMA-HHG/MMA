import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn
# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index

  
omegaSI = mn.ConvertPhoton(800e-9, 'lambdaSI', 'omegaSI') 
Horder = 17


gas_type = 'Kr'
XUV_table_type_absorption = 'Henke' # {Henke, NIST}    
XUV_table_type_dispersion = 'NIST'
def f2_funct(E):
    return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)

def sigma(omega):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 2.0 * units.r_electron_classical * lambdaSI * f2_value


susc_IR = IR_index.getsusc('Ar', mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))

N_atm = 2.7e19 * 1e6
# N_atm = 2.6867774e25
f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(Horder*omegaSI, 'omegaSI', 'eV'))
nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(Horder*omegaSI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
susc_XUV = nXUV_atm**2 - 1

delta_susc = (susc_IR - susc_XUV)/N_atm

fact = units.elcharge**2 / (omegaSI**2 *units.eps0*units.elmass)
def IXUV(eta):    
    return 1.0/((sigma(Horder*omegaSI)**2) + ((Horder*omegaSI)**2/units.c_light**2)*(delta_susc - eta * fact)**2 )

etagrid = np.linspace(0, 0.2, 200)

eta_optimal = (omegaSI**2 * units.eps0*units.elmass/(units.elcharge**2)) * 3.9e-29#delta_susc

print(eta_optimal)

IXUV_res = IXUV(etagrid)

IXUV_res = np.asarray(IXUV_res)
    
fig, ax = plt.subplots()     
plt.plot(etagrid,IXUV_res)
plt.show()
   
