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
import Hfn
import Hfn2

from scipy import integrate


import Hankel_tools

# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index


import plot_presets as pp




# inputs from hdf5-input


gas_type = 'Ar'
XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
XUV_table_type_absorption = 'Henke' # {Henke, NIST} 
apply_diffraction = ['dispersion', 'absorption']

Nr_max = 235 #470; 235; 155-still fine    
Hrange = [16, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]

kr_step = 2 # descending order, the last is "the most accurate"
ko_step = 1

rmax_FF = 8*1e-4
Nr_FF = 50 # 10 # 200
distance_FF = 1.

FF_orders_plot = 4    
Nz_max_sum = 5 # 41

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

    results_CUPRAD = os.path.join("C:\sharepoint", "OneDrive - ELI Beamlines", "TEMP", "Sunrise","test1")

    results_TDSE = os.path.join("C:\sharepoint", "OneDrive - ELI Beamlines", "TEMP", "Sunrise","test1")
    




file_CUPRAD = os.path.join(results_CUPRAD,file_CUPRAD)
file_TDSE = os.path.join(results_TDSE,file_TDSE)


rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)


# load data
print('processing:', file_CUPRAD, file_TDSE)             
with h5py.File(file_CUPRAD, 'r') as InputArchiveCUPRAD, h5py.File(file_TDSE, 'r') as InputArchiveTDSE:
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchiveCUPRAD,'/inputs/laser_wavelength','N'),'lambdaSI','omegaau')
    inverse_GV_IR = InputArchiveCUPRAD['/logs/inverse_group_velocity_SI'][()]; group_velocity_IR = 1./inverse_GV_IR
    # pressure_mbar = 1e3*InputArchiveCUPRAD['/inputs/medium_pressure_in_bar'][()]
    rho0_init = 1e6 * mn.readscalardataset(InputArchiveCUPRAD, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N') # SI
   
    FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,:,0] + \
                        1j*InputArchiveTDSE['FSourceTerm'][:,:,:,1]
    ogrid = InputArchiveTDSE['omegagrid'][:]
    rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
    zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]
    
    FSourceTerm_sparse = InputArchiveTDSE['FSourceTerm'][:,:,0:-1:2,0] + \
                        1j*InputArchiveTDSE['FSourceTerm'][:,:,0:-1:2,1]
  
                        
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.sf[0].args = [ogrid, np.abs(FSourceTerm[0,0,:])]
    image.sf[0].method = plt.semilogy
    pp.plot_preset(image)
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.sf[0].args = [ogrid/omega0, np.abs(FSourceTerm[0,0,:])]
    image.sf[0].method = plt.semilogy
    pp.plot_preset(image)
    
    ko_min = mn.FindInterval(ogrid/omega0, 16)
    ko_max = mn.FindInterval(ogrid/omega0, 20)
    
    # print(type(InputArchiveTDSE))
    # grp =  InputArchiveCUPRAD['/logs']
    # print(type(grp))
    
    omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
    ogridSI = omega_au2SI * ogrid
    omega0SI = omega_au2SI * omega0
    
    target_static = Hankel_tools.FSources_provider(InputArchiveTDSE['zgrid_coarse'][:],
                                                   InputArchiveTDSE['rgrid_coarse'][:],
                                                   omega_au2SI*InputArchiveTDSE['omegagrid'][:],
                                                   FSource = np.transpose(FSourceTerm,axes=(1,2,0)),
                                                   data_source = 'static',
                                                   ko_min = ko_min,
                                                   ko_max = ko_max)
    
    # target_dynamic = Hankel_tools.FSources_provider(InputArchiveTDSE['zgrid_coarse'][:],
    #                                                InputArchiveTDSE['rgrid_coarse'][:],
    #                                                omega_au2SI*InputArchiveTDSE['omegagrid'][:],
    #                                                h5_handle = InputArchiveTDSE,
    #                                                h5_path = 'FSourceTerm',
    #                                                data_source = 'dynamic',
    #                                                ko_min = ko_min,
    #                                                ko_max = ko_max)
    
    # target_static_Ar = Hankel_tools.FSources_provider(InputArchiveTDSE['zgrid_coarse'][:],
    #                                                InputArchiveTDSE['rgrid_coarse'][:],
    #                                                omega_au2SI*InputArchiveTDSE['omegagrid'][:],
    #                                                FSource = np.transpose(FSourceTerm,axes=(1,2,0)),
    #                                                data_source = 'static',
    #                                                ko_min = ko_min,
    #                                                ko_max = ko_max)
    
    
    # plane1_dyn = next(target_dynamic.Fsource_plane)
    # plane2_dyn = next(target_dynamic.Fsource_plane)
        
    
    HL_end, HL_cum = Hfn2.HankelTransform_long(target_static, # FSourceTerm(r,z,omega)
                              distance_FF, rgrid_FF,
                              preset_gas = 'vacuum',
                              pressure = 1.,
                              absorption_tables = 'Henke',
                              include_absorption = True,
                              dispersion_tables = 'Henke',
                              include_dispersion = True,
                              effective_IR_refrective_index = 1.,
                              integrator_Hankel = integrate.trapz,
                              integrator_longitudinal = 'trapezoidal',
                              near_field_factor = True,
                              store_cummulative_result = True,
                              frequencies_to_trace_maxima = None,
                              )
    


    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.sf[0].args = [target_static.ogrid/omega0SI, rgrid_FF, np.abs(HL_cum[0].T)]
    image.sf[0].method = plt.pcolormesh
    pp.plot_preset(image)
    
    
    image.sf[0].args[-1] = np.abs(HL_cum[1].T)
    pp.plot_preset(image)


    image.sf[0].args[-1] = np.abs(HL_cum[3].T)
    pp.plot_preset(image)
    

    image.sf[0].args[-1] = np.abs(HL_cum[6].T)
    pp.plot_preset(image)
    

    image.sf[0].args[-1] = np.abs(HL_cum[9].T)
    pp.plot_preset(image)    
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.sf[0].args = [target_static.ogrid/omega0SI, rgrid_FF, np.abs(HL_end.T)]
    # image.sf[0].method = plt.pcolormesh
    # pp.plot_preset(image)
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.sf[0].args = [target_static.ogrid/omega0SI, rgrid_FF, np.abs(Hankel_long_static_Ar.T)]
    # image.sf[0].method = plt.pcolormesh
    # pp.plot_preset(image)






    # Hankel_long_static_Ar = Hfn2.HankelTransform_long(target_static_Ar, # FSourceTerm(r,z,omega)
    #                           distance_FF, rgrid_FF,
    #                           preset_gas = 'Ar',
    #                           pressure = 1.,
    #                           absorption_tables = 'Henke',
    #                           include_absorption = True,
    #                           dispersion_tables = 'Henke',
    #                           include_dispersion = True,
    #                           effective_IR_refrective_index = 1.,
    #                           integrator_Hankel = integrate.trapz,
    #                           integrator_longitudinal = 'trapezoidal',
    #                           near_field_factor = True,
    #                           store_cummulative_result = False,
    #                           frequencies_to_trace_maxima = None
    #                           )

    # Hankel_long_dynamic = Hfn2.HankelTransform_long(target_dynamic, # FSourceTerm(r,z,omega)
    #                           distance_FF, rgrid_FF,
    #                           preset_gas = 'vacuum',
    #                           pressure = 1.,
    #                           absorption_tables = 'Henke',
    #                           include_absorption = True,
    #                           dispersion_tables = 'Henke',
    #                           include_dispersion = True,
    #                           effective_IR_refrective_index = 1.,
    #                           integrator_Hankel = integrate.trapz,
    #                           integrator_longitudinal = 'trapezoidal',
    #                           near_field_factor = True,
    #                           store_cummulative_result = False,
    #                           frequencies_to_trace_maxima = None
    #                           )