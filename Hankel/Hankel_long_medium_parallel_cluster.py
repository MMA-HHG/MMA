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


import XUV_refractive_index as XUV_index




import multiprocessing as mp

# Nthreads = 4


with open('msg.tmp','r') as msg_file:
    file = msg_file.readline()[:-1] # need to strip the last character due to Fortran msg.tmp


# inputs from hdf5-input
with h5py.File(file, 'r') as InpArch:
    
    inp_group = InpArch[MMA.paths['Hankel_inputs']]
    
    Nthreads = mn.readscalardataset(inp_group, 'Nthreads','N')  
    XUV_table_type_diffraction = mn.readscalardataset(
        inp_group, 'XUV_table_type_dispersion','S')    
    XUV_table_type_absorption = mn.readscalardataset(
        inp_group, 'XUV_table_type_absorption','S') 
    
    Nr_max = mn.readscalardataset(inp_group, 'Nr_max','N') 
    
    Hrange = inp_group['Harmonic_range'][:] 
    kr_step = mn.readscalardataset(inp_group, 'kr_step','N') 
    ko_step = mn.readscalardataset(inp_group, 'ko_step','N') 
    
    rmax_FF = mn.readscalardataset(inp_group, 'rmax_FF','N') 
    Nr_FF = mn.readscalardataset(inp_group, 'Nr_FF','N') 
    distance_FF = mn.readscalardataset(inp_group, 'distance_FF','N') 







# rmax_FF = 6*1e-4
# Nr_FF = 50 # 25 # 50 # 10 # 200
# distance_FF = 1.

# # FF_orders_plot = 4    
# # Nz_max_sum = 5 # 41

# file_CUPRAD = 'results.h5'
# file_TDSE = 'results_merged.h5'
# out_h5name = 'test_Hankel.h5'


# arguments = sys.argv

# showplots = not('-nodisplay' in arguments)

# if ('-here' in arguments):
#     results_path = os.getcwd()
#     results_CUPRAD = os.getcwd()
#     results_TDSE = os.getcwd()
# else:   

#     results_path = os.path.join("D:\sharepoint", "OneDrive - ELI Beamlines",
#                     "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t1")
    
#     results_path = os.path.join("/mnt", "d","sharepoint", "OneDrive - ELI Beamlines",
#                     "data", "Sunrise","tmp","h5debug","TDSEs","densmod","t1")

    


# filename = "results.h5"

# file = os.path.join(results_path,filename)





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
    omega0SI = omega_au2SI * omega0
    
    
    absorption = True
    dispersion = True
    
    
    ogrid_sel = ogrid[ko_min:ko_max]
    
    No_sel = len(ogrid_sel)
    
    print('No', No_sel, 'Nr_FF', Nr_FF)
    
    if (Nr_FF >= No_sel):
        rgrid_FF_parts = np.array_split(rgrid_FF, Nthreads)
        ogrid_parts = [ogrid_sel for _ in range(Nthreads)]
        
        rgrid_FF_indices = [mn.FindInterval(rgrid_FF, rgrid_FF_part[0]) for rgrid_FF_part in rgrid_FF_parts]
        
        ogrid_indices_start = Nthreads*[ko_min]
        ogrid_indices_end = Nthreads*[ko_max]
        
        
    else:
        ogrid_parts = np.array_split(ogrid_sel, Nthreads)
        rgrid_FF_parts = [rgrid_FF for _ in range(Nthreads)]
        rgrid_FF_indices = Nthreads*[0]
        
        ogrid_indices = [mn.FindInterval(ogrid, ogrid_part[0]) for ogrid_part in ogrid_parts] # From the original file used to load the data !!!
        ogrid_indices.append(ko_max)
        ogrid_indices_start = ogrid_indices[:-1]
        ogrid_indices_end = ogrid_indices[1:]
        
    targets = []
    for k1 in range(Nthreads):
        targets.append(HT.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ogrid_indices_start[k1],
                                                    ko_max = ogrid_indices_end[k1])
                       )
        
    # Parallel calculations of Hankel
    task_queue = mp.Queue()    
    def mp_handle(position, *args, **kwargs):
        task_queue.put(
                       [position, HT.Hankel_long(*args,**kwargs)] # the outputs are not ordered, keep the order in the result
                      )
    
    # Fsources
    processes = [mp.Process(target=mp_handle, # define processes
                            args=(k1,
                                  targets[k1],
                                  distance_FF,
                                  rgrid_FF_parts[k1]),
                            kwargs={
                                    'preset_gas': preset_gas,
                                    'pressure' : pressure,
                                    'absorption_tables' : XUV_table_type_absorption,
                                    'include_absorption' : absorption,
                                    'dispersion_tables' : XUV_table_type_diffraction,
                                    'include_dispersion' : dispersion,
                                    'effective_IR_refrective_index' : effective_IR_refrective_index,
                                    'integrator_Hankel' : integrate.trapz,
                                    'integrator_longitudinal' : 'trapezoidal',
                                    'near_field_factor' : True,
                                    'store_cummulative_result' : True,
                                    'store_non_normalised_cummulative_result' : True
                                   }
                            
                            ) for k1 in range(Nthreads)]

    for p in processes: p.start(); # run processes

    results = [task_queue.get() for p in processes] # there is no ordering!
    
    
    # merge results
    HL_res_p = copy.deepcopy(results[0][1])
    
    HL_res_p.ogrid = omega_au2SI*ogrid_sel
    ogrid_indices_new = [mn.FindInterval(ogrid_sel, ogrid_part[0]) for ogrid_part in ogrid_parts]
    
    oldshape = np.shape(HL_res_p.FF_integrated)
    newshape = (Nr_FF, No_sel)
    HL_res_p.FF_integrated = np.empty(newshape,dtype=HL_res_p.FF_integrated.dtype)

    
    
    oldshape = np.shape(HL_res_p.entry_plane_transform)
    newshape = (Nr_FF, No_sel)
    HL_res_p.entry_plane_transform = np.empty(newshape,dtype=HL_res_p.cummulative_field.dtype)
    HL_res_p.exit_plane_transform = np.empty(newshape,dtype=HL_res_p.cummulative_field.dtype)
    
    if 'cummulative_field' in dir(HL_res_p):
        oldshape = np.shape(HL_res_p.cummulative_field)
        newshape = (oldshape[0],) + (Nr_FF, No_sel)
        HL_res_p.cummulative_field = np.empty(newshape,dtype=HL_res_p.cummulative_field.dtype)
        
    if 'cummulative_field_no_norm' in dir(HL_res_p):
        oldshape = np.shape(HL_res_p.cummulative_field_no_norm)
        newshape = (oldshape[0],) + (Nr_FF, No_sel)
        HL_res_p.cummulative_field_no_norm = np.empty(newshape,dtype=HL_res_p.cummulative_field_no_norm.dtype)
    
    for k1, k_worker in enumerate([result[0] for result in results]):
        
        HL_res_p.FF_integrated[
        rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
        ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
        ] = results[k1][1].FF_integrated
        
        HL_res_p.entry_plane_transform[
        rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
        ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
        ] = results[k1][1].entry_plane_transform

        HL_res_p.exit_plane_transform[
        rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
        ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
        ] = results[k1][1].exit_plane_transform            
        
        if 'cummulative_field' in dir(HL_res_p):
            HL_res_p.cummulative_field[:,
            rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
            ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
            ] = results[k1][1].cummulative_field
            
        if 'cummulative_field_no_norm' in dir(HL_res_p):
            HL_res_p.cummulative_field_no_norm[:,
            rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
            ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
            ] = results[k1][1].cummulative_field_no_norm
    
    
    
    # target_dynamic = HT.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
    #                                                 InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
    #                                                 omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
    #                                                 h5_handle = InpArch,
    #                                                 h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
    #                                                 data_source = 'dynamic',
    #                                                 ko_min = ko_min,
    #                                                 ko_max = ko_max)
    
    

    
    
    
    # HL_res = HT.Hankel_long(target_dynamic,
    #                         distance_FF,
    #                         rgrid_FF,
    #                         preset_gas = preset_gas,
    #                         pressure = pressure,
    #                         absorption_tables = XUV_table_type_absorption,
    #                         include_absorption = absorption,
    #                         dispersion_tables = XUV_table_type_diffraction,
    #                         include_dispersion = dispersion,
    #                         effective_IR_refrective_index = effective_IR_refrective_index,
    #                         integrator_Hankel = integrate.trapz,
    #                         integrator_longitudinal = 'trapezoidal',
    #                         near_field_factor = True,
    #                         store_cummulative_result = True,
    #                         store_non_normalised_cummulative_result = True)
    
    with h5py.File('results_Hankel', 'a') as Hres_file:
        out_group = Hres_file.create_group(MMA.paths['Hankel_outputs'])
        
        mn.adddataset(out_group,
                      'FF_integrated',
                      np.stack((HL_res_p.FF_integrated.real, 
                                HL_res_p.FF_integrated.imag),axis=-1),
                      '[arb. u.]')
        
        mn.adddataset(out_group,
                      'entry_plane_transform',
                      np.stack((HL_res_p.entry_plane_transform.real, 
                                HL_res_p.entry_plane_transform.imag),axis=-1),
                      '[arb. u.]')

        mn.adddataset(out_group,
                      'exit_plane_transform',
                      np.stack((HL_res_p.exit_plane_transform.real, 
                                HL_res_p.exit_plane_transform.imag),axis=-1),
                      '[arb. u.]')
        
        mn.adddataset(out_group,
                      'ogrid',
                      HL_res_p.ogrid,
                      '[SI]')

        mn.adddataset(out_group,
                      'rgrid',
                      HL_res_p.rgrid,
                      '[SI]')

        
        if 'cummulative_field' in dir(HL_res_p):
            mn.adddataset(out_group,
                          'cummulative_field',
                          np.stack((HL_res_p.cummulative_field.real, 
                                    HL_res_p.cummulative_field.imag),axis=-1),
                          '[arb. u.]')
            
        if 'cummulative_field_no_norm' in dir(HL_res_p):
            mn.adddataset(out_group,
                          'cummulative_field_no_norm',
                          np.stack((HL_res_p.cummulative_field_no_norm.real, 
                                    HL_res_p.cummulative_field_no_norm.imag),
                                    axis=-1),
                          '[arb. u.]')
        
        if 'zgrid' in dir(HL_res_p):
            mn.adddataset(out_group,
                          'zgrid',
                          HL_res_p.zgrid,
                          '[SI]')
            
    

        