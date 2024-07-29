import numpy as np
from scipy import integrate
import h5py
import copy
import multiprocessing as mp


import MMA_administration as MMA
import units
import mynumerics as mn
import Hankel_transform as HT

omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')

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







rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)


# load data
print('processing:', file)             
with h5py.File(file, 'r') as InpArch:
    
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InpArch,
                                                        MMA.paths['CUPRAD_inputs']+
                                                        '/laser_wavelength','N'),'lambdaSI','omegaau')
    inverse_GV_IR = InpArch[MMA.paths['CUPRAD_logs']+'/inverse_group_velocity_SI'][()]; group_velocity_IR = 1./inverse_GV_IR
    rho0_init = 1e6 * mn.readscalardataset(InpArch, MMA.paths['CUPRAD_inputs']+
                                           '/calculated/medium_effective_density_of_neutral_molecules','N') # SI
    pressure = MMA.pressure_constructor(InpArch)
    preset_gas = mn.readscalardataset(InpArch,MMA.paths['global_inputs']+'/gas_preset','S')
    effective_IR_refrective_index = inverse_GV_IR*units.c_light
    

    ogrid = InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:]          # a.u.
    rgrid_macro = InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:] # SI
    zgrid_macro = InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:] # SI
    
    # the inidces of the selection in the frequency (harmonic) grid
    ko_min = mn.FindInterval(ogrid/omega0, Hrange[0])
    ko_max = mn.FindInterval(ogrid/omega0, Hrange[-1])

    
    
    omega0SI = omega_au2SI * omega0
    
    
    absorption = True
    dispersion = True
    
    
    ogrid_sel = ogrid[ko_min:ko_max]    
    No_sel = len(ogrid_sel)
    
    # print('No', No_sel, 'Nr_FF', Nr_FF)
    
    ## Parallel computing:
    # Decide which of the dimension on the screen is bigger and then apply
    # multiprocessing over this dimension. The process is then unified and
    # the code executes 'Nthreads' simulations using mp.Queue(). The
    if (Nr_FF >= No_sel): # If there are more radial points
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
    
    # prepare instances of 'FSources_provider' class, each of the instances
    # describes a subarray according to the splitting above
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
    
    task_queue = mp.Queue() # que to store the results from multiprocessing 
    def mp_handle(position, *args, **kwargs): # handle to trace the position of the subarrays
        task_queue.put(
                       [position, HT.Hankel_long(*args,**kwargs)] # the outputs are not ordered, keep the order in the result
                      )
    
    # define processes
    processes = [mp.Process(target=mp_handle, 
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

    # run the processes in parallel
    for p in processes: p.start()
    results = [task_queue.get() for p in processes] # The results are not ordered
    # result = [[position_index, Hankel_long-class-instance], ... ]
    
    
    ## Merge results ##
    
    ## Prepare new shared output class
    
    HL_res = copy.deepcopy(results[0][1]) # ensure all shared data are in the output class
    
    HL_res.ogrid = omega_au2SI*ogrid_sel  # reassing by the whole ogrid
    # find indices for merging subarrays
    ogrid_indices_new = [mn.FindInterval(ogrid_sel, ogrid_part[0]) for ogrid_part in ogrid_parts]
    
    # outputs that are always present
    # total integrated signal
    newshape = (Nr_FF, No_sel)
    HL_res.FF_integrated = np.empty(newshape,dtype=HL_res.FF_integrated.dtype)
    # Hankel transforms of the first and last planes
    newshape = (Nr_FF, No_sel)
    HL_res.entry_plane_transform = np.empty(newshape,dtype=HL_res.cummulative_field.dtype)
    HL_res.exit_plane_transform  = np.empty(newshape,dtype=HL_res.cummulative_field.dtype)
    
    # outputs that are present only of required
    if 'cummulative_field' in dir(HL_res):
        oldshape = np.shape(HL_res.cummulative_field)
        newshape = (oldshape[0],) + (Nr_FF, No_sel)
        HL_res.cummulative_field = np.empty(newshape,dtype=HL_res.cummulative_field.dtype)
        
    if 'cummulative_field_no_norm' in dir(HL_res):
        oldshape = np.shape(HL_res.cummulative_field_no_norm)
        newshape = (oldshape[0],) + (Nr_FF, No_sel)
        HL_res.cummulative_field_no_norm = np.empty(newshape,dtype=HL_res.cummulative_field_no_norm.dtype)
    
    ## copy data into the newly allocated class
    for k1, k_worker in enumerate([result[0] for result in results]):
        
        HL_res.FF_integrated[ # choosing proper indices to accomodate the subarrays
        rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
        ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
        ] = results[k1][1].FF_integrated
        
        HL_res.entry_plane_transform[
        rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
        ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
        ] = results[k1][1].entry_plane_transform

        HL_res.exit_plane_transform[
        rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
        ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
        ] = results[k1][1].exit_plane_transform            
        
        if 'cummulative_field' in dir(HL_res):
            HL_res.cummulative_field[:,
            rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
            ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
            ] = results[k1][1].cummulative_field
            
        if 'cummulative_field_no_norm' in dir(HL_res):
            HL_res.cummulative_field_no_norm[:,
            rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
            ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
            ] = results[k1][1].cummulative_field_no_norm
    
    
    # Save the results
    with h5py.File('results_Hankel.h5', 'a') as Hres_file:
        out_group = Hres_file.create_group(MMA.paths['Hankel_outputs'])        
        mn.adddataset(out_group,
                      'FF_integrated',
                      np.stack((HL_res.FF_integrated.real, # real and imaginary parts are in separate dimensions
                                HL_res.FF_integrated.imag),axis=-1),
                      '[arb. u.]')        
        mn.adddataset(out_group,
                      'entry_plane_transform',
                      np.stack((HL_res.entry_plane_transform.real, 
                                HL_res.entry_plane_transform.imag),axis=-1),
                      '[arb. u.]')
        mn.adddataset(out_group,
                      'exit_plane_transform',
                      np.stack((HL_res.exit_plane_transform.real, 
                                HL_res.exit_plane_transform.imag),axis=-1),
                      '[arb. u.]')        
        mn.adddataset(out_group,
                      'ogrid',
                      HL_res.ogrid,
                      '[SI]')
        mn.adddataset(out_group,
                      'rgrid',
                      HL_res.rgrid,
                      '[SI]')

        # optional outputs
        if 'cummulative_field' in dir(HL_res):
            mn.adddataset(out_group,
                          'cummulative_field',
                          np.stack((HL_res.cummulative_field.real, 
                                    HL_res.cummulative_field.imag),axis=-1),
                          '[arb. u.]')            
        if 'cummulative_field_no_norm' in dir(HL_res):
            mn.adddataset(out_group,
                          'cummulative_field_no_norm',
                          np.stack((HL_res.cummulative_field_no_norm.real, 
                                    HL_res.cummulative_field_no_norm.imag),
                                    axis=-1),
                          '[arb. u.]')        
        if 'zgrid' in dir(HL_res):
            mn.adddataset(out_group,
                          'zgrid',
                          HL_res.zgrid,
                          '[SI]')

print('The parallel Hankel transform finishes.')