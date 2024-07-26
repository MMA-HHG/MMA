import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import copy
import units
import mynumerics as mn


import Hankel_transform as HT

from scipy import integrate


# import Hankel_tools
import MMA_administration as MMA


# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index


import plot_presets as pp

import dataformat_CUPRAD as dfC




# inputs from hdf5-input


# gas_type = 'Ar'
XUV_table_type_diffraction = 'NIST' # {Henke, NIST} # NIST for SciRep #Henke for Ar
XUV_table_type_absorption = 'Henke' # {Henke, NIST}  # Henke for SciRep
# apply_diffraction = ['dispersion', 'absorption']

Nr_max = 235 #470; 235; 155-still fine    
Hrange = [16, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]

kr_step = 2 # descending order, the last is "the most accurate"
ko_step = 1

rmax_FF = 6*1e-4
Nr_FF = 25 # 10 # 200
distance_FF = 1.

# FF_orders_plot = 4    
# Nz_max_sum = 5 # 41

file_CUPRAD = 'results.h5'
file_TDSE = 'results_merged.h5'
out_h5name = 'test_Hankel.h5'


arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
    results_CUPRAD = os.getcwd()
    results_TDSE = os.getcwd()
else:   

    results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t1")
    
    results_path = os.path.join("/mnt", "d","sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t1")

    


filename = "results.h5"

file = os.path.join(results_path,filename)





rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)


# load data
print('processing:', file)             
with h5py.File(file, 'r') as InpArch:
    
    gas_type = mn.readscalardataset(InpArch, MMA.paths['global_inputs']+
                                           '/gas_preset','S')
    
    print('gas_preset from hdf5: ', gas_type)
    # print(InpArch.keys())
    # print('h5path',MMA.paths['CUPRAD_inputs']+'/laser_wavelength')
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InpArch,
                                                        MMA.paths['CUPRAD_inputs']+
                                                        '/laser_wavelength','N'),'lambdaSI','omegaau')
    inverse_GV_IR = InpArch[MMA.paths['CUPRAD_logs']+'/inverse_group_velocity_SI'][()]; group_velocity_IR = 1./inverse_GV_IR
    
    
    # pressure_mbar = 1e3*InputArchiveCUPRAD['/inputs/medium_pressure_in_bar'][()]
    rho0_init = 1e6 * mn.readscalardataset(InpArch, MMA.paths['CUPRAD_inputs']+
                                           '/calculated/medium_effective_density_of_neutral_molecules','N') # SI
    
    pressure = MMA.pressure_constructor(InpArch)


    preset_gas = mn.readscalardataset(InpArch,MMA.paths['global_inputs']+'/gas_preset','S')
    
    effective_IR_refrective_index = inverse_GV_IR*units.c_light

    

    ogrid = InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:]
    rgrid_macro = InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:]
    zgrid_macro = InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:]
    

    


    
    ko_min = mn.FindInterval(ogrid/omega0, 16.8)
    ko_max = mn.FindInterval(ogrid/omega0, 17.3)
    

    
    omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
    ogridSI = omega_au2SI * ogrid
    omega0SI = omega_au2SI * omega0
    
    
    
    target_dynamic = HT.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ko_min,
                                                    ko_max = ko_max)
    
    
    absorption = True
    dispersion = True
    
    
    
    HL_res = HT.Hankel_long(target_dynamic,
                            distance_FF,
                            rgrid_FF,
                            preset_gas = preset_gas,
                            pressure = pressure,
                            absorption_tables = XUV_table_type_absorption,
                            include_absorption = absorption,
                            dispersion_tables = XUV_table_type_diffraction,
                            include_dispersion = dispersion,
                            effective_IR_refrective_index = effective_IR_refrective_index,
                            integrator_Hankel = integrate.trapz,
                            integrator_longitudinal = 'trapezoidal',
                            near_field_factor = True,
                            store_cummulative_result = True,
                            store_non_normalised_cummulative_result = True)
    

    


    # signal build-up for H17
    ko_17 = mn.FindInterval(target_dynamic.ogrid/omega0SI, 17)

    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "H17, normalisation test"
    image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_res.cummulative_field_no_norm[:,ko_17,:]),axis=1)]
    image.sf[1].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_res.cummulative_field[:,ko_17,:]),axis=1)]
    pp.plot_preset(image)

    


    