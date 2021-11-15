import numpy as np
import os
import time
import shutil
import h5py
import sys
import units
import mynumerics as mn
import matplotlib.pyplot as plt
import re
import glob
import copy

# from contextlib import ExitStack
# import dataformat_CUPRAD as dfC


# files management
preion_extensions = ['half','end']

N_extensions = len(preion_extensions)
results = glob.glob(os.path.join('*','results_*.h5'))

ionisation_init = {}
for extension in preion_extensions:
    ionisation_init[extension] = []
ionisation_tmax = copy.deepcopy(ionisation_init); ionisation_end_pulse = copy.deepcopy(ionisation_init)

p_grid = []

for fname in results:
    numbers = re.findall(r'\d+',  os.path.basename(fname))
    p_grid.append(float(numbers[0]))
    
p_grid = np.unique(p_grid); Np = len(p_grid)
ionisation_init = {}
for extension in preion_extensions:
    ionisation_init[extension] = np.zeros((Np,))
ionisation_tmax = copy.deepcopy(ionisation_init); ionisation_end_pulse = copy.deepcopy(ionisation_init)

# load data
for fname in results:
    for extension in preion_extensions:
        if ('_'+extension) in fname: break
    else:
        ValueError('wrong fname-extension match')
        
    numbers = re.findall(r'\d+',  os.path.basename(fname)); p_value = float(numbers[0])
    k_press = np.where(p_grid == p_value)[0][0]
    
    print(fname)
    with h5py.File(fname, 'r') as InputArchive:

        pressure_mbar = 1e3*InputArchive['/inputs/medium_pressure_in_bar'][()]
        omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
        rho0_init = 1e6 * mn.readscalardataset(InputArchive, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')
        
        Nz = len(InputArchive['/outputs/zgrid'][:]) 
        
        preion = 100.*InputArchive['/pre_ionised/initial_electrons_ratio'][()]
        
        # find maximal intensity and tmax
        E_slice = InputArchive['/outputs/output_field'][0,:,0] # Nz-1
        plasma_slice = InputArchive['/outputs/output_plasma'][0,:,0]
        tgrid = InputArchive['/outputs/tgrid'][:]; Nt = len(tgrid)
        
        rem_fast_oscillations = np.exp(-1j*omega0*tgrid)
                    
        
        E_slice_envel = rem_fast_oscillations*mn.complexify_fft(E_slice)
        Intens_slice = mn.FieldToIntensitySI(abs(E_slice_envel))
        index_of_max = np.argmax(Intens_slice)
        print(tgrid[index_of_max+1])
        plasma_tmax = 100.*plasma_slice[index_of_max]/rho0_init
        
        
        ionisation_init[extension][k_press] = preion
        ionisation_tmax[extension][k_press] = plasma_tmax
        ionisation_end_pulse[extension][k_press] = 100.*plasma_slice[-1]/rho0_init
        




# svae results
out_h5name = 'ionisations.h5'

with h5py.File(out_h5name,'w') as OutFile: # this file contains numerical analyses
    mn.adddataset(OutFile,'p_grid',p_grid,'[mbar]')
    
    for extension in preion_extensions:
        mn.adddataset(OutFile,'ionisation_'+extension+'_init',ionisation_init[extension],'[%]')
        mn.adddataset(OutFile,'ionisation_'+extension+'_tmax',ionisation_tmax[extension],'[%]')
        mn.adddataset(OutFile,'ionisation_'+extension+'_end_pulse',ionisation_end_pulse[extension],'[%]')


print('done')