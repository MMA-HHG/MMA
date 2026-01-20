"""
This script provides a multiprocessing wrapper of the Hankel transform module, namely the Hankel_long
routine, it acts as a part of the main computational pipeline. It does the following steps:

1. It reads the inputs from the main archive.
2. It splits the workload (subarrays of the far-field camera screen) and define que of processes using
   the multiprocessing module. 
3. It merges the results.
4. It stores the results in a separate file*. The results can be merged into the main file by the
additional script copy_results_to_main.py.

*The reason for using a separate file is that it might be desirable for users to rerun this part of the 
pipeline repeatedly: with different resolutions and ranges of the camera, test the physics without 
absorption, etc.

@author: Jan Vábek
"""

# Import standard modules
import numpy as np
import h5py
import copy
import multiprocessing as mp

# Import custom modules from this project
import MMA_administration as MMA
import units
import mynumerics as mn
import Hankel_transform as HT

# prepare conversion factor for the frequencies from atomic units to SI
omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')

# specify the input archive transferred in the temporary file 'msg.tmp'
with open('msg.tmp','r') as msg_file:
    file = msg_file.readline()[:-1] # need to strip the last character due to Fortran msg.tmp


# read the simulation parameters from the hdf5 archive

# # load the data from the hdf5 archive
# # ! We use the dynamic access to the data: it means the data are loaded during the calculation and that
# # the calculation must be done inside the with block.
# print('processing:', file)             
# with h5py.File(file, 'r') as InpArch:

# Open the source hdf5 archive:
# 1 - read input parameters
# 2 - perform the calculation
#     ! The whole calculation must be done inside the with block, data are accessed dynamically
with h5py.File(file, 'r') as InpArch:
    

    # We load the inputs from the file and store them in the new file together with possible defaults for a reference.
    inp_group = InpArch[MMA.paths['Hankel_inputs']]
    with h5py.File('results_Hankel.h5', 'a') as Hres_file:
        inps_full_group = Hres_file.create_group(MMA.paths['Hankel_inputs_full'])

        # number of threads for multiprocessing

        Nthreads = mn.readscalardataset(inp_group, 'Nthreads','N')  

        # disperison and absorption

        XUV_table_type_diffraction = mn.readscalardataset(
            inp_group, 'XUV_table_type_dispersion','S')    
        XUV_table_type_absorption = mn.readscalardataset(
            inp_group, 'XUV_table_type_absorption','S') 
        
        if ('include_dispersion' in inp_group):
            dispersion = (mn.readscalardataset(inp_group, 'include_dispersion','N') == 1) # optional
            inp_group.copy('include_dispersion',inps_full_group)
        else:
            dispersion = True   
            mn.adddataset(inps_full_group,'include_dispersion',1,'[-]')

        if ('include_absorption' in inp_group):
            absorption = (mn.readscalardataset(inp_group, 'include_absorption','N') == 1) # optional
            inp_group.copy('include_absorption',inps_full_group)
        else:
            absorption = True   
            mn.adddataset(inps_full_group,'include_absorption',1,'[-]')
    
        # Definition of the far-field camera screen

        distance_FF = mn.readscalardataset(inp_group, 'distance_FF','N') 

        rmax_FF = mn.readscalardataset(inp_group, 'rmax_FF','N') 
        Nr_FF = mn.readscalardataset(inp_group, 'Nr_FF','N')

        Hrange = inp_group['Harmonic_range'][:] 
        if ('ko_step' in inp_group):
            ko_step = mn.readscalardataset(inp_group, 'ko_step','N')  
            inp_group.copy('ko_step',inps_full_group)
        else:
            ko_step = 1
            mn.adddataset(inps_full_group,'ko_step',ko_step,'[-]')

        # Defintion of the integration variables
        kr_max = mn.readscalardataset(inp_group, 'Nr_max','N')

        if ('kr_step' in inp_group):
            kr_step = mn.readscalardataset(inp_group, 'kr_step','N')  
            inp_group.copy('kr_step',inps_full_group)
        else:
            kr_step = 1
            mn.adddataset(inps_full_group,'kr_step',kr_step,'[-]')

        if ('kz_step' in inp_group):
            kz_step = mn.readscalardataset(inp_group, 'kz_step','N')  
            inp_group.copy('kz_step',inps_full_group)
        else:
            kz_step = 1
            mn.adddataset(inps_full_group,'kz_step',kz_step,'[-]')

        # possible additional characteristics and advanced options 

        if ('store_cumulative_result' in inp_group):
            store_cumulative_result = (mn.readscalardataset(inp_group, 'store_cumulative_result','N') == 1) # optional
            inp_group.copy('store_cumulative_result',inps_full_group)
        else:
            store_cumulative_result = False
            mn.adddataset(inps_full_group,'store_cumulative_result',0,'[-]')

        if ('store_cumulative_result_nonorm' in inp_group):
            store_cumulative_result_nonorm = (mn.readscalardataset(inp_group, 'store_cumulative_result_nonorm','N') == 1) # optional
            inp_group.copy('store_cumulative_result_nonorm',inps_full_group)
        else:
            store_cumulative_result_nonorm = False
            mn.adddataset(inps_full_group,'store_cumulative_result_nonorm',0,'[-]')

        if ('near_field_factor' in inp_group):
            near_field_factor = (mn.readscalardataset(inp_group, 'near_field_factor','N') == 1) # optional
            inp_group.copy('near_field_factor',inps_full_group)
        else:
            near_field_factor = True
            mn.adddataset(inps_full_group,'near_field_factor',1,'[-]')

        # copy all the obligatory inputs for a reference:
        for dset in ['Nthreads', 'XUV_table_type_dispersion', 'XUV_table_type_absorption',
                     'Harmonic_range', 'Nr_max', 'rmax_FF', 'Nr_FF', 'distance_FF']:
            inp_group.copy(dset,inps_full_group)


    ## Prepare calculation

    rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)    
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InpArch,
                                                        MMA.paths['CUPRAD_inputs']+
                                                        '/laser_wavelength','N'),'lambdaSI','omegaau')
    omega0SI = omega_au2SI * omega0 
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
      
    # prepare the grid
    ogrid_sel = ogrid[ko_min:ko_max:ko_step]    
    No_sel = len(ogrid_sel)
    
    print('Calculation is running:', file)  
    print('The resolution (N_omega x N_r): (', No_sel, ')(', Nr_FF,')')
    print('------------------------------------------------')
    
    ## Parallel computing:
    # Decide which of the dimension on the screen is bigger and then apply
    # multiprocessing over this dimension. The process is then unified and
    # the code executes 'Nthreads' simulations using mp.Queue().
    if (Nr_FF >= No_sel): # If there are more radial points
        paralellised_in_r = True
        rgrid_FF_parts = np.array_split(rgrid_FF, Nthreads)
        ogrid_parts = [ogrid_sel for _ in range(Nthreads)]
        
        rgrid_FF_indices = [mn.FindInterval(rgrid_FF, rgrid_FF_part[0]) for rgrid_FF_part in rgrid_FF_parts]
        
        ogrid_indices_start = Nthreads*[ko_min]
        ogrid_indices_end = Nthreads*[ko_max]
        
        
    else:
        paralellised_in_r = False
        ogrid_parts = np.array_split(ogrid_sel, Nthreads)
        rgrid_FF_parts = [rgrid_FF for _ in range(Nthreads)]
        rgrid_FF_indices = Nthreads*[0]
        
        ogrid_indices = [mn.FindInterval(ogrid, ogrid_part[0]) for ogrid_part in ogrid_parts] # From the original file used to load the data !!!
        ogrid_indices.append(ko_max)
        ogrid_indices_start = ogrid_indices[:-1]
        ogrid_indices_end = ogrid_indices[1:]
    
    # prepare instances of 'FSources_provider' class, each of the instances
    # describes a subarray according to the splitting above, note the 'dynamic' option (reding from the file on-the-fly)
    targets = []
    for k1 in range(Nthreads):
        targets.append(HT.FSources_provider(InpArch[MMA.paths['CTDSE_outputs']+'/zgrid_coarse'][:],
                                                    InpArch[MMA.paths['CTDSE_outputs']+'/rgrid_coarse'][:],
                                                    omega_au2SI*InpArch[MMA.paths['CTDSE_outputs']+'/omegagrid'][:],
                                                    h5_handle = InpArch,
                                                    h5_path = MMA.paths['CTDSE_outputs']+'/FSourceTerm',
                                                    data_source = 'dynamic',
                                                    ko_min = ogrid_indices_start[k1],
                                                    ko_max = ogrid_indices_end[k1],
                                                    ko_step=ko_step,
                                                    kr_max=kr_max,
                                                    kr_step=kr_step,
                                                    kz_step=kz_step)

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
                                    'integrator_Hankel' : HT.trapezoidal_integrator,
                                    'integrator_longitudinal' : 'trapezoidal',
                                    'near_field_factor' : near_field_factor,
                                    'store_cumulative_result' : store_cumulative_result,
                                    'store_non_normalised_cumulative_result' : store_cumulative_result_nonorm
                                   }
                            
                            ) for k1 in range(Nthreads)]

    # run the processes in parallel
    for p in processes: p.start()
    results = [task_queue.get() for p in processes] # The results are not ordered
    # result = [[position_index, Hankel_long-class-instance], ... ]
    
    
    ## Merge results ##
    
    ## Prepare new shared output class
    # This part of the code copy the output from the 0th process, extends the array, and fill the arrays by data from the other processes
    
    HL_res = copy.deepcopy(results[0][1]) # ensure all shared data are in the output class
    
    HL_res.ogrid = omega_au2SI*ogrid_sel  # reassign by the main non-split ogrid
    # find indices for merging subarrays
    ogrid_indices_new = [mn.FindInterval(ogrid_sel, ogrid_part[0]) for ogrid_part in ogrid_parts]
    
    # outputs that are always present
    # total integrated signal
    newshape = (Nr_FF, No_sel)
    HL_res.FF_integrated = np.empty(newshape,dtype=HL_res.FF_integrated.dtype)
    # Hankel transforms of the first and last planes
    newshape = (Nr_FF, No_sel)
    HL_res.entry_plane_transform = np.empty(newshape,dtype=HL_res.entry_plane_transform.dtype)
    HL_res.exit_plane_transform  = np.empty(newshape,dtype=HL_res.exit_plane_transform.dtype)
    
    # outputs that are present only if required
    if 'cumulative_field' in dir(HL_res):
        oldshape = np.shape(HL_res.cumulative_field)
        newshape = (oldshape[0],) + (Nr_FF, No_sel)
        HL_res.cumulative_field = np.empty(newshape,dtype=HL_res.cumulative_field.dtype)
        
    if 'cumulative_field_no_norm' in dir(HL_res):
        oldshape = np.shape(HL_res.cumulative_field_no_norm)
        newshape = (oldshape[0],) + (Nr_FF, No_sel)
        HL_res.cumulative_field_no_norm = np.empty(newshape,dtype=HL_res.cumulative_field_no_norm.dtype)
    
    ## copy the original radial grid in the case it was split for parallelisation
    if paralellised_in_r:
        HL_res.rgrid = rgrid_FF

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
        
        if 'cumulative_field' in dir(HL_res):
            HL_res.cumulative_field[:,
            rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
            ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
            ] = results[k1][1].cumulative_field
            
        if 'cumulative_field_no_norm' in dir(HL_res):
            HL_res.cumulative_field_no_norm[:,
            rgrid_FF_indices[k_worker]:(rgrid_FF_indices[k_worker]+len(rgrid_FF_parts[k_worker])),
            ogrid_indices_new[k_worker]:(ogrid_indices_new[k_worker]+len(ogrid_parts[k_worker]))
            ] = results[k1][1].cumulative_field_no_norm
    
    
    ## Save the results
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
        if 'cumulative_field' in dir(HL_res):
            mn.adddataset(out_group,
                          'cumulative_field',
                          np.stack((HL_res.cumulative_field.real, 
                                    HL_res.cumulative_field.imag),axis=-1),
                          '[arb. u.]')            
        if 'cumulative_field_no_norm' in dir(HL_res):
            mn.adddataset(out_group,
                          'cumulative_field_no_norm',
                          np.stack((HL_res.cumulative_field_no_norm.real, 
                                    HL_res.cumulative_field_no_norm.imag),
                                    axis=-1),
                          '[arb. u.]')        
        if 'zgrid' in dir(HL_res):
            mn.adddataset(out_group,
                          'zgrid',
                          HL_res.zgrid,
                          '[SI]')

print('The parallel Hankel transform finishes.')