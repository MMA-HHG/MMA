import numpy as np
import os
import time
import shutil
import h5py
import sys
import units
import mynumerics as mn
import HHG
import matplotlib.pyplot as plt
import re
import glob
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index
# from contextlib import ExitStack
# import dataformat_CUPRAD as dfC


# files management
preion_extensions = ['half','end']
results = glob.glob(os.path.join('*','results_*.h5'))

available_points = {}
for extension in preion_extensions:
    available_points[extension] = []

p_grid = []

for fname in results:
    numbers = re.findall(r'\d+',  os.path.basename(fname))
    p_grid.append(float(numbers[0]))
p_grid = np.unique(p_grid)

print(results)
print(p_grid)
print(available_points)

# load data

# pressure_mbar = 1e3*InputArchive['/inputs/medium_pressure_in_bar'][()]
# omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
# rho0_init = 1e6 * mn.readscalardataset(InputArchive, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')

# Nz = len(InputArchive['/outputs/zgrid'][:]) 

# preion = InputArchive['/pre_ionised/initial_electrons_ratio'][()]

# E_slice = InputArchive['/outputs/output_field'][Nz-1,:,0]
# plasma_slice = InputArchive['/outputs/output_plasma'][Nz-1,:,0]
# tgrid = InputArchive['/outputs/tgrid'][:]; Nt = len(tgrid)

# rem_fast_oscillations = np.exp(-1j*omega0*tgrid)
            

# E_slice_envel = rem_fast_oscillations*mn.complexify_fft(E_slice)
# Intens_slice = mn.FieldToIntensitySI(abs(E_slice_envel))
# plasma_tmax = 100.*plasma_slice[index_of_max]/rho0_init
                
                
# find maximal intensity and tmax

# ionisations at tmax

# svae results