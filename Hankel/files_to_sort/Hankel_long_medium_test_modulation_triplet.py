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

    
    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","t1")    

    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","t2") 

    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","SciRep","t1")

    # results_path2 = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","SciRep","t2")
    
    
    # results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","t3")
    
    # results_path2 = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
    #                 "data", "Sunrise","tmp","h5debug","TDSEs","t3mod")
    

    results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t1")
    
    results_path2 = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t2")
    
    results_path3 = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
                    "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t3")

    


file = "results_TDSEM.h5"
filename = "results.h5"

file = os.path.join(results_path,filename)
file2 = os.path.join(results_path2,filename)
file3 = os.path.join(results_path3,filename)




rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)


# load data
print('processing:', file)             
with h5py.File(file, 'r') as InpArch, h5py.File(file2, 'r') as InpArch2, h5py.File(file3, 'r') as InpArch3:
    
    gas_type = mn.readscalardataset(InpArch, MMA.paths['global_inputs']+
                                           '/gas_preset','S')
    
    print('gas_preset from hdf5: ', gas_type)
    # print(InpArch.keys())
    # print('h5path',MMA.paths['CUPRAD_inputs']+'/laser_wavelength')
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InpArch,
                                                        MMA.paths['CUPRAD_inputs']+
                                                        '/laser_wavelength','N'),'lambdaSI','omegaau')
    
    inverse_GV_IR = InpArch[MMA.paths['CUPRAD_logs']+'/inverse_group_velocity_SI'][()]; group_velocity_IR = 1./inverse_GV_IR
    inverse_GV_IR2 = InpArch2[MMA.paths['CUPRAD_logs']+'/inverse_group_velocity_SI'][()]; group_velocity_IR2 = 1./inverse_GV_IR2
    inverse_GV_IR3 = InpArch3[MMA.paths['CUPRAD_logs']+'/inverse_group_velocity_SI'][()]; group_velocity_IR2 = 1./inverse_GV_IR3
    
    # pressure_mbar = 1e3*InputArchiveCUPRAD['/inputs/medium_pressure_in_bar'][()]
    rho0_init = 1e6 * mn.readscalardataset(InpArch, MMA.paths['CUPRAD_inputs']+
                                           '/calculated/medium_effective_density_of_neutral_molecules','N') # SI

    
    pressure = MMA.pressure_constructor(InpArch)
    pressure2 = MMA.pressure_constructor(InpArch2)
    pressure3 = MMA.pressure_constructor(InpArch3)

    
    
    print('used gas:', MMA.paths['global_inputs']+'/gas_preset')

    preset_gas = mn.readscalardataset(InpArch,MMA.paths['global_inputs']+'/gas_preset','S')
    
    effective_IR_refrective_index = inverse_GV_IR*units.c_light
    effective_IR_refrective_index2 = inverse_GV_IR2*units.c_light
    effective_IR_refrective_index3 = inverse_GV_IR3*units.c_light
    
    
    ogrid = InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:]  

    
    ko_min = mn.FindInterval(ogrid/omega0, 16.8)
    ko_max = mn.FindInterval(ogrid/omega0, 17.3)
    

    
    omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
    ogridSI = omega_au2SI * ogrid
    omega0SI = omega_au2SI * omega0
    
    
    CUPRAD_res = dfC.get_data(InpArch)
    CUPRAD_res2 = dfC.get_data(InpArch2)
    CUPRAD_res3 = dfC.get_data(InpArch3)
    

    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "End fields orig"
    image.sf[0].args = [CUPRAD_res.tgrid, CUPRAD_res.E_zrt[-1,0,:]]
    image.sf[1].args = [CUPRAD_res2.tgrid, CUPRAD_res2.E_zrt[-1,0,:]]
    image.sf[2].args = [CUPRAD_res3.tgrid, CUPRAD_res3.E_zrt[-1,0,:]]
    pp.plot_preset(image)
       

    
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "End fields grid shift"
    # image.sf[0].args = [CUPRAD_res.tgrid+delta_t_tot1, CUPRAD_res.E_zrt[-1,0,:]]
    # image.sf[1].args = [CUPRAD_res2.tgrid+delta_t_tot2, CUPRAD_res2.E_zrt[-1,0,:]]
    image.sf[0].args = [CUPRAD_res.co_moving_t_grid(-1) , CUPRAD_res.E_zrt[-1,0,:]]
    image.sf[1].args = [CUPRAD_res2.co_moving_t_grid(-1) , CUPRAD_res2.E_zrt[-1,0,:]]
    image.sf[2].args = [CUPRAD_res3.co_moving_t_grid(-1) , CUPRAD_res3.E_zrt[-1,0,:]]
    pp.plot_preset(image)
    
    # sys.exit(0)
    
    CUPRAD_res.vacuum_shift(output='add')
    CUPRAD_res2.vacuum_shift(output='add')
    CUPRAD_res3.vacuum_shift(output='add')
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "End fields vac"
    image.sf[0].args = [CUPRAD_res.tgrid, CUPRAD_res.E_zrt_vac[-1,0,:]]
    image.sf[1].args = [CUPRAD_res2.tgrid, CUPRAD_res2.E_zrt_vac[-1,0,:]]
    image.sf[2].args = [CUPRAD_res3.tgrid, CUPRAD_res2.E_zrt_vac[-1,0,:]]
    pp.plot_preset(image)
    
    # sys.exit(0)
    
    CUPRAD_res.get_plasma(InpArch)
    CUPRAD_res2.get_plasma(InpArch2)
    CUPRAD_res3.get_plasma(InpArch3)
    
    CUPRAD_res.compute_spectrum()
    CUPRAD_res2.compute_spectrum()
    CUPRAD_res3.compute_spectrum()

    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "Spectra"
    image.sf[0].args = [CUPRAD_res.ogrid/CUPRAD_res.omega0, np.abs(CUPRAD_res.FE_zrt[-1,0,:])]
    image.sf[0].method = plt.semilogy
    image.sf[1].args = [CUPRAD_res2.ogrid/CUPRAD_res2.omega0, np.abs(CUPRAD_res2.FE_zrt[-1,0,:])]
    image.sf[1].method = plt.semilogy
    image.sf[2].args = [CUPRAD_res3.ogrid/CUPRAD_res3.omega0, np.abs(CUPRAD_res3.FE_zrt[-1,0,:])]
    image.sf[2].method = plt.semilogy
    pp.plot_preset(image)
    
    
    # End spectra
    HHG_Spectrum = (InpArch[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][-1,0,:,0] +
                    1j*InpArch[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][-1,0,:,1])
    HHG_Spectrum2 = (InpArch2[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][-1,0,:,0] +
                    1j*InpArch2[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][-1,0,:,1])
    HHG_Spectrum3 = (InpArch3[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][-1,0,:,0] +
                    1j*InpArch3[MMA.paths['CTDSE_outputs']+'/FSourceTerm'][-1,0,:,1])
    HHG_ogrid = InpArch3[MMA.paths['CTDSE_outputs']+'/omegagrid'][:]
    

    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "Harmonic Spectra"
    image.sf[0].args = [HHG_ogrid/omega0, np.abs(HHG_Spectrum)]
    image.sf[0].method = plt.semilogy
    image.sf[1].args = [HHG_ogrid/omega0, np.abs(HHG_Spectrum2)]
    image.sf[1].method = plt.semilogy
    image.sf[2].args = [HHG_ogrid/omega0, np.abs(HHG_Spectrum3)]
    image.sf[2].method = plt.semilogy
    pp.plot_preset(image)
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "Spectra entry"
    image.sf[0].args = [CUPRAD_res.ogrid/CUPRAD_res.omega0, np.abs(CUPRAD_res.FE_zrt[0,0,:])]
    image.sf[0].method = plt.semilogy
    image.sf[1].args = [CUPRAD_res2.ogrid/CUPRAD_res2.omega0, np.abs(CUPRAD_res2.FE_zrt[0,0,:]),'--']
    image.sf[1].method = plt.semilogy
    image.sf[2].args = [CUPRAD_res3.ogrid/CUPRAD_res2.omega0, np.abs(CUPRAD_res3.FE_zrt[0,0,:]),'--']
    image.sf[2].method = plt.semilogy
    pp.plot_preset(image)
    
    sys.exit(0)
       
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "Plasma"
    image.sf[0].args = [CUPRAD_res.plasma.tgrid, CUPRAD_res.plasma.value_zrt[-1,0,:]]
    image.sf[1].args = [CUPRAD_res2.plasma.tgrid, CUPRAD_res2.plasma.value_zrt[-1,0,:]]
    image.sf[2].args = [CUPRAD_res3.plasma.tgrid, CUPRAD_res3.plasma.value_zrt[-1,0,:]]
    pp.plot_preset(image)
    
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "Plasma shifted"
    image.sf[0].args = [CUPRAD_res.co_moving_t_grid(CUPRAD_res.zgrid[-1]), CUPRAD_res.plasma.value_zrt[-1,0,:]]
    image.sf[1].args = [CUPRAD_res2.co_moving_t_grid(CUPRAD_res2.zgrid[-1]), CUPRAD_res2.plasma.value_zrt[-1,0,:]]
    image.sf[2].args = [CUPRAD_res3.co_moving_t_grid(CUPRAD_res3.zgrid[-1]), CUPRAD_res3.plasma.value_zrt[-1,0,:]]
    pp.plot_preset(image)
    
    
    

    
    target_dynamic = HT.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ko_min,
                                                    ko_max = ko_max)
    

    target_dynamic2 = HT.FSources_provider(InpArch2[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch2[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch2[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch2,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ko_min,
                                                    ko_max = ko_max)      
    
    target_dynamic3 = HT.FSources_provider(InpArch2[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch2[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch3[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch3,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ko_min,
                                                    ko_max = ko_max)    

    
    # sys.exit(0)
    
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
    
    HL_res2 = HT.Hankel_long(target_dynamic2,
                            distance_FF,
                            rgrid_FF,
                            preset_gas = preset_gas,
                            pressure = pressure2,
                            absorption_tables = XUV_table_type_absorption,
                            include_absorption = absorption,
                            dispersion_tables = XUV_table_type_diffraction,
                            include_dispersion = dispersion,
                            effective_IR_refrective_index = effective_IR_refrective_index2,
                            integrator_Hankel = integrate.trapz,
                            integrator_longitudinal = 'trapezoidal',
                            near_field_factor = True,
                            store_cummulative_result = True,
                            store_non_normalised_cummulative_result = True)
    
    
    HL_res3 = HT.Hankel_long(target_dynamic3,
                            distance_FF,
                            rgrid_FF,
                            preset_gas = preset_gas,
                            pressure = pressure3,
                            absorption_tables = XUV_table_type_absorption,
                            include_absorption = absorption,
                            dispersion_tables = XUV_table_type_diffraction,
                            include_dispersion = dispersion,
                            effective_IR_refrective_index = effective_IR_refrective_index3,
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

    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "H17"
    image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_res.cummulative_field[:,ko_17,:]),axis=1)]
    image.sf[1].args = [1e3*target_dynamic2.zgrid[1:], np.max(np.abs(HL_res2.cummulative_field[:,ko_17,:]),axis=1)]
    image.sf[2].args = [1e3*target_dynamic3.zgrid[1:], np.max(np.abs(HL_res3.cummulative_field[:,ko_17,:]),axis=1)]
    pp.plot_preset(image)
        
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(32)]
    image.title = "H17 nonorm"
    image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_res.cummulative_field_no_norm[:,ko_17,:]),axis=1)]
    image.sf[1].args = [1e3*target_dynamic2.zgrid[1:], np.max(np.abs(HL_res2.cummulative_field_no_norm[:,ko_17,:]),axis=1)]
    image.sf[2].args = [1e3*target_dynamic3.zgrid[1:], np.max(np.abs(HL_res3.cummulative_field_no_norm[:,ko_17,:]),axis=1)]
    pp.plot_preset(image)        
        

    # No = len(target_dynamic.ogrid)
    # dispersion_factor = np.empty(No)
    # dispersion_factor_new = np.empty(No)
    # for k1 in range(No):
    #     dispersion_factor[k1] = target_dynamic.ogrid[k1]*dispersion_function_def(target_dynamic.ogrid[k1]) 
    #     dispersion_factor_new[k1] = target_dynamic.ogrid[k1]*XUV_index.dispersion_function(
    #                                 target_dynamic.ogrid[k1], 
    #                                 pressure,                
    #                                 preset_gas+'_'+XUV_table_type_diffraction,
    #                                 n_IR = effective_IR_refrective_index)
    
    
    # absorption_factor = np.empty(No)
    # absorption_factor_new = np.empty(No)
    # for k1 in range(No):
    #     absorption_factor[k1] = target_dynamic.ogrid[k1]*absorption_function_def(target_dynamic.ogrid[k1])
    #     absorption_factor_new[k1] = target_dynamic.ogrid[k1]*(pressure/units.c_light) *\
    #                                 XUV_index.beta_factor_ref(
    #                                     target_dynamic.ogrid[k1],
    #                                     preset_gas+'_'+XUV_table_type_absorption)
                                
      
    # # compute z-evolution of the factors        
    # factor_e_ref_reconstructed = pressure*np.exp(
    #                       1j*np.outer(target_dynamic.zgrid,dispersion_factor) +
    #                       np.outer(target_dynamic.zgrid-target_dynamic.zgrid[-1] ,absorption_factor)
    #                       )
    
    # factor_arg_e = np.outer(target_dynamic.zgrid-target_dynamic.zgrid[-1] ,absorption_factor)
    
    # factor_e_Htools = np.asarray([ pf(k1)[0,:] for k1 in range(len(target_dynamic.zgrid))])


    # test1 = np.abs(factor_e_ref/factor_e_Htools)   
    


    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.sf[0].args = [target_dynamic.ogrid/omega0SI, rgrid_FF, np.abs(HL_cum[-1].T)]
    # image.sf[0].method = plt.pcolormesh
    # image.sf[0].colorbar.show = True
    # pp.plot_preset(image)
    
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.sf[0].args = [target_dynamic.ogrid/omega0SI, rgrid_FF, np.abs(HL_cum_ref[-1].T)]
    # image.sf[0].method = plt.pcolormesh
    # image.sf[0].colorbar.show = True
    # pp.plot_preset(image)
    
    
    # # signal build-up for H19
    # ko_17 = mn.FindInterval(target_dynamic.ogrid/omega0SI, 17)
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "H17"
    # image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_cum[:,ko_17,:]),axis=1)]
    # image.sf[1].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_cum_ref[:,ko_17,:]),axis=1)]
    # pp.plot_preset(image)
    
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "ratio H17 [-]"
    # image.sf[0].args = [1e3*target_dynamic.zgrid[1:], np.max(np.abs(HL_cum[:,ko_17,:]),axis=1)/np.max(np.abs(HL_cum_ref[:,ko_17,:]),axis=1)]


    # pp.plot_preset(image)
    
    



















    # CUPRAD_res = dfC.get_data(InpArch)
    # CUPRAD_res2 = dfC.get_data(InpArch2)
    

    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "End fields orig"
    # image.sf[0].args = [CUPRAD_res.tgrid, CUPRAD_res.E_zrt[-1,0,:]]
    # image.sf[1].args = [CUPRAD_res2.tgrid, CUPRAD_res2.E_zrt[-1,0,:]]
    # pp.plot_preset(image)
       
    # # tgrid1_end_shift = CUPRAD_res.tgrid 
    # # tgrid2_end_shift = CUPRAD_res2.tgrid 
    
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "End fields grid shift"
    # # image.sf[0].args = [CUPRAD_res.tgrid+delta_t_tot1, CUPRAD_res.E_zrt[-1,0,:]]
    # # image.sf[1].args = [CUPRAD_res2.tgrid+delta_t_tot2, CUPRAD_res2.E_zrt[-1,0,:]]
    # image.sf[0].args = [CUPRAD_res.co_moving_t_grid(CUPRAD_res.zgrid[-1]) , CUPRAD_res.E_zrt[-1,0,:]]
    # image.sf[1].args = [CUPRAD_res2.co_moving_t_grid(CUPRAD_res2.zgrid[-1]) , CUPRAD_res2.E_zrt[-1,0,:]]
    # pp.plot_preset(image)
    
    # # sys.exit(0)
    
    # CUPRAD_res.vacuum_shift(output='add')
    # CUPRAD_res2.vacuum_shift(output='add')
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "End fields vac"
    # image.sf[0].args = [CUPRAD_res.tgrid, CUPRAD_res.E_zrt_vac[-1,0,:]]
    # image.sf[1].args = [CUPRAD_res2.tgrid, CUPRAD_res2.E_zrt_vac[-1,0,:],'--']
    # pp.plot_preset(image)
    
    # # sys.exit(0)
    
    # CUPRAD_res.get_plasma(InpArch)
    # CUPRAD_res2.get_plasma(InpArch2)
    
    # CUPRAD_res.compute_spectrum()
    # CUPRAD_res2.compute_spectrum()

    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "Spectra"
    # image.sf[0].args = [CUPRAD_res.ogrid/CUPRAD_res.omega0, np.abs(CUPRAD_res.FE_zrt[-1,0,:])]
    # image.sf[0].method = plt.semilogy
    # image.sf[1].args = [CUPRAD_res2.ogrid/CUPRAD_res2.omega0, np.abs(CUPRAD_res2.FE_zrt[-1,0,:]),'--']
    # image.sf[1].method = plt.semilogy
    # pp.plot_preset(image)
    

    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "Spectra entry"
    # image.sf[0].args = [CUPRAD_res.ogrid/CUPRAD_res.omega0, np.abs(CUPRAD_res.FE_zrt[0,0,:])]
    # image.sf[0].method = plt.semilogy
    # image.sf[1].args = [CUPRAD_res2.ogrid/CUPRAD_res2.omega0, np.abs(CUPRAD_res2.FE_zrt[0,0,:]),'--']
    # image.sf[1].method = plt.semilogy
    # pp.plot_preset(image)
       
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "Plasma"
    # image.sf[0].args = [CUPRAD_res.plasma.tgrid, CUPRAD_res.plasma.value_zrt[-1,0,:]]
    # image.sf[1].args = [CUPRAD_res2.plasma.tgrid, CUPRAD_res2.plasma.value_zrt[-1,0,:],'--']
    # pp.plot_preset(image)
    
    
    # image = pp.figure_driver()
    # image.sf = [pp.plotter() for k1 in range(32)]
    # image.title = "Plasma shifted"
    # image.sf[0].args = [CUPRAD_res.co_moving_t_grid(CUPRAD_res.zgrid[-1]), CUPRAD_res.plasma.value_zrt[-1,0,:]]
    # image.sf[1].args = [CUPRAD_res2.co_moving_t_grid(CUPRAD_res2.zgrid[-1]), CUPRAD_res2.plasma.value_zrt[-1,0,:],'--']
    # pp.plot_preset(image)
    
    # # sys.exit(0)
    