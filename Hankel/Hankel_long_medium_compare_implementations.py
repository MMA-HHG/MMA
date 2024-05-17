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
import Hfn2
import Hfn2_v1

from scipy import integrate


import Hankel_tools
import MMA_administration as MMA


# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index


import plot_presets as pp




# inputs from hdf5-input


# gas_type = 'Ar'
XUV_table_type_diffraction = 'Henke' # {Henke, NIST}
XUV_table_type_absorption = 'Henke' # {Henke, NIST} 
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

    
    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","t1")    

    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","t2") 

    results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","SciRep","t1")

    results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","SciRep","t2")
    
    results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","t3")
    
    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","SciRep","t1")
    


file = "results_TDSEM.h5"
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
    
    
    
    pressure = Hankel_tools.pressure_constructor(InpArch)

    
    
    print(MMA.paths['global_inputs']+'/gas_preset')
    print(InpArch[MMA.paths['global_inputs']].keys())
    # xxx = InpArch[MMA.paths['global_inputs']+'/gas_preset'][()]
    # yyy = InpArch[MMA.paths['global_inputs']+'/gas_preset'].decode()
    preset_gas = mn.readscalardataset(InpArch,MMA.paths['global_inputs']+'/gas_preset','S')
    
    effective_IR_refrective_index = inverse_GV_IR*units.c_light
    
    # try:
    #     preset_gas = mn.readscalardataset(InpArch,MMA.paths['global_inputs']+'/gas_preset','S')
    # except:
    #     preset_gas = 'vacuum' 
    
    # pressure = 1.
    # preset_gas = 'vacuum'
    

    
    # FSourceTerm =    InpArch[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][:,:,:,0] + \
    #               1j*InpArch[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][:,:,:,1]
    ogrid = InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:]
    rgrid_macro = InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:]
    zgrid_macro = InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:]
    


    
    ko_min = mn.FindInterval(ogrid/omega0, 16.8)
    ko_max = mn.FindInterval(ogrid/omega0, 17.3)
    
    
    FSourceTerm =    InpArch[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][:,:,ko_min:ko_max,0] + \
                  1j*InpArch[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][:,:,ko_min:ko_max,1]
    

    
    omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
    ogridSI = omega_au2SI * ogrid
    omega0SI = omega_au2SI * omega0
    
    # target_static = Hankel_tools.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
    #                                                InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
    #                                                omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
    #                                                FSource = np.transpose(FSourceTerm,axes=(0,2,1)),
    #                                                data_source = 'static',
    #                                                ko_min = ko_min,
    #                                                ko_max = ko_max)
    
    target_dynamic = Hankel_tools.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ko_min,
                                                    ko_max = ko_max)
    
    
    
    # Here are the fuction to obtain the phase factors in SI units: exp(i*omega*function(omega))
    def f1_funct(E):
        return XUV_index.getf1(gas_type+'_' + XUV_table_type_diffraction, E)
    def f2_funct(E):
        return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)


    def dispersion_function_def(omega):
        f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
        lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
        nXUV     = 1.0 - rho0_init*units.r_electron_classical * \
                   ((lambdaSI**2)*f1_value/(2.0*np.pi))           
        phase_velocity_XUV  = units.c_light / nXUV
        return ((1./group_velocity_IR) - (1./phase_velocity_XUV))

    def absorption_function_def(omega):
        f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
        lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
        beta_factor = rho0_init*units.r_electron_classical * \
                      ((lambdaSI**2)*f2_value/(2.0*np.pi))
        return beta_factor / units.c_light
    
    
    # test dispersion functions
    df_t1 = dispersion_function_def(target_dynamic.ogrid[0])
    df_t2 = XUV_index.dispersion_function(target_dynamic.ogrid[0], pressure, gas_type+'_'+XUV_table_type_diffraction, n_IR=effective_IR_refrective_index) 
    
    print('disp functions_test:' , df_t1,  df_t2)
    
    # sys.exit(0)
    
    absorption = True
    dispersion = False
    

    
    HL_end, HL_cum, pf, HL_cum_test, abs_factor_Htools =  Hfn2.HankelTransform_long(target_dynamic, # FSourceTerm(r,z,omega)
                              distance_FF, rgrid_FF,
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
                              frequencies_to_trace_maxima = None,
                              )
    
    
    
    # original implementation


    if dispersion:  dispersion_function = dispersion_function_def
    else:           dispersion_function = None   
    
    if absorption:  absorption_function = absorption_function_def
    else:           absorption_function = None  

    
    HL_end_ref, HL_cum_ref, factor_e_ref, abs_factor_Hfn_ref = Hfn2_v1.HankelTransform_long(
                                                   target_dynamic.ogrid,
                                                   target_dynamic.rgrid,
                                                   target_dynamic.zgrid,
                                                   np.transpose(FSourceTerm,axes=(1,0,2)),
                                                   distance_FF,
                                                   rgrid_FF,
                                                   dispersion_function = dispersion_function, # None, #dispersion_function,
                                                   absorption_function = absorption_function,
                                                   store_cummulative_result = True)
    
    HL_cum_ref *= pressure
    

    
    

    No = len(target_dynamic.ogrid)
    dispersion_factor = np.empty(No)
    dispersion_factor_new = np.empty(No)
    for k1 in range(No):
        dispersion_factor[k1] = target_dynamic.ogrid[k1]*dispersion_function_def(target_dynamic.ogrid[k1]) 
        dispersion_factor_new[k1] = target_dynamic.ogrid[k1]*XUV_index.dispersion_function(
                                    target_dynamic.ogrid[k1], 
                                    pressure,                
                                    preset_gas+'_'+XUV_table_type_diffraction,
                                    n_IR = effective_IR_refrective_index)
    
    
    absorption_factor = np.empty(No)
    absorption_factor_new = np.empty(No)
    for k1 in range(No):
        absorption_factor[k1] = target_dynamic.ogrid[k1]*absorption_function_def(target_dynamic.ogrid[k1])
        absorption_factor_new[k1] = target_dynamic.ogrid[k1]*(pressure/units.c_light) *\
                                    XUV_index.beta_factor_ref(
                                        target_dynamic.ogrid[k1],
                                        preset_gas+'_'+XUV_table_type_absorption)
                                
      
    # compute z-evolution of the factors        
    factor_e_ref_reconstructed = pressure*np.exp(
                          1j*np.outer(target_dynamic.zgrid,dispersion_factor) +
                          np.outer(target_dynamic.zgrid-target_dynamic.zgrid[-1] ,absorption_factor)
                          )
    
    factor_e_Htools = np.asarray([ pf(k1)[0,:] for k1 in range(len(target_dynamic.zgrid))])


    test1 = np.abs(factor_e_ref/factor_e_Htools)   
    


    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.sf[0].args = [target_dynamic.ogrid/omega0SI, rgrid_FF, np.abs(HL_cum[-1].T)]
    image.sf[0].method = plt.pcolormesh
    image.sf[0].colorbar.show = True
    pp.plot_preset(image)
    
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.sf[0].args = [target_dynamic.ogrid/omega0SI, rgrid_FF, np.abs(HL_cum_ref[-1].T)]
    image.sf[0].method = plt.pcolormesh
    image.sf[0].colorbar.show = True
    pp.plot_preset(image)
    
    
    # signal build-up for H19
    ko_17 = mn.FindInterval(target_dynamic.ogrid/omega0SI, 17)
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "H17"
    image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_cum[:,ko_17,:]),axis=1)]
    image.sf[1].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_cum_ref[:,ko_17,:]),axis=1)]
    pp.plot_preset(image)
    
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "ratio H17 [-]"
    image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_cum[:,ko_17,:]),axis=1)/np.max(np.abs(HL_cum_ref[:,ko_17,:]),axis=1)]


    pp.plot_preset(image)
    
    
