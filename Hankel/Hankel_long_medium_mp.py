import numpy as np
import os
import time
import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn
import Hfn
import Hfn2

# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index


arguments = sys.argv

if ("-input" in arguments):
    arg_index = arguments.index("-input")
    inputs_archive = arguments[arg_index+1]
else:
    inputs_archive = 'inputs_Hankel.h5'

# inputs from hdf5-input
try:
    with h5py.File(inputs_archive, 'r') as Parameters:
        inputs_group = Parameters['inputs']
        inputs_list = list(inputs_group.keys())
        
        file_CUPRAD = mn.readscalardataset(Parameters, 'inputs/file_CUPRAD', 'S')  
        file_TDSE = mn.readscalardataset(Parameters, 'inputs/file_TDSE', 'S')  
        out_h5name = mn.readscalardataset(Parameters, 'inputs/out_h5name', 'S')  
        
        gas_type = mn.readscalardataset(Parameters, 'inputs/gas_type', 'S')    
        XUV_table_type_diffraction = mn.readscalardataset(Parameters,
                                     'inputs/XUV_table_type_diffraction', 'S')
        XUV_table_type_absorption = mn.readscalardataset(Parameters,
                                     'inputs/XUV_table_type_absorption', 'S') 
        apply_diffraction = Parameters['inputs/apply_diffraction'][:].tolist()
        apply_diffraction = [k1.decode() for k1 in apply_diffraction]
        
        
        kr_step = mn.readscalardataset(Parameters, 'inputs/kr_step', 'N')
        ko_step = mn.readscalardataset(Parameters, 'inputs/ko_step', 'N')       
        Nr_max = mn.readscalardataset(Parameters, 'inputs/Nr_max', 'N')        
             
        Hrange = Parameters['inputs/Hrange'][:].tolist()   
        
        redefine_integration = False
        if ('Nz_first' in inputs_list):
            kz_start, kz_end = 0, mn.readscalardataset(Parameters, 'inputs/Nz_first', 'N')
            print('integrating over first Nz planes')
        elif ('Nz_last' in inputs_list):
            redefine_integration = True
            kz_end = 'end'
            Nz_last = mn.readscalardataset(Parameters, 'inputs/Nz_last', 'N')
            print('integrating over last Nz planes')
        elif ('kz_integrate' in inputs_list):
            kz_start, kz_end = Parameters['inputs/kz_integrate'][:].tolist()
            print('integrating over slected region')
        else:
            kz_start, kz_end = 0,'end'
            print('default integration: whole zgrid')

            
        # Nz_max_sum = mn.readscalardataset(Parameters, 'inputs/Nz_max_sum', 'N')  
        
        Nr_FF = mn.readscalardataset(Parameters, 'inputs/Nr_FF', 'N')
        rmax_FF = mn.readscalardataset(Parameters, 'inputs/rmax_FF', 'N')
        distance_FF = mn.readscalardataset(Parameters, 'inputs/distance_FF', 'N')
        
        
        FF_orders_plot = mn.readscalardataset(Parameters, 'inputs/FF_orders_plot', 'N')
        Workers = mn.readscalardataset(Parameters, 'inputs/N_threads', 'N')
        
        store_cummulative_result = False
        if ('store_cummulative_result' in inputs_list):
            store_cummulative_result = (mn.readscalardataset(Parameters, 'inputs/store_cummulative_result', 'N') == 1)

        print('reading from hdf5-input file')

        
except:
    print('error reading hdf5 input file, using defaults') 
    gas_type = 'Kr'
    XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
    XUV_table_type_absorption = 'Henke' # {Henke, NIST} 
    apply_diffraction = ['dispersion', 'absorption']
    
    Nr_max = 235 #470; 235; 155-still fine    
    Hrange = [17, 18] # [17, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]
    
    kr_step = 2 # descending order, the last is "the most accurate"
    ko_step = 2
    
    rmax_FF = 8*1e-4
    Nr_FF = 200 # 10 # 200
    distance_FF = 0.3
    
    FF_orders_plot = 4    
    # Nz_max_sum = 5 # 41
    kz_start, kz_end = 0,5
    
    file_CUPRAD = 'results_1.h5'
    file_TDSE = 'results_merged_Fsource.h5'
    out_h5name = 'Hankel.h5'
    
    Workers = 4
    
    store_cummulative_result = True


showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
    results_CUPRAD = os.getcwd()
    results_TDSE = os.getcwd()
else:
    results_CUPRAD = os.path.join("D:\data", "Discharges", "TDSE", "t6")
    results_CUPRAD = os.path.join("/mnt","d","data", "Discharges", "TDSE", "t6")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSEH1")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE40planes1")
    results_TDSE = os.path.join("C:\data", "Discharges", "TDSE", "TDSE50planes1") # fine up to 46  
    results_TDSE = os.path.join("/mnt","c","data", "Discharges", "TDSE", "TDSE50planes1")
    




file_CUPRAD = os.path.join(results_CUPRAD,file_CUPRAD)
file_TDSE = os.path.join(results_TDSE,file_TDSE)


# load data
print('processing:', file_CUPRAD, file_TDSE)             
with h5py.File(file_CUPRAD, 'r') as InputArchiveCUPRAD, h5py.File(file_TDSE, 'r') as InputArchiveTDSE:
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchiveCUPRAD,'/inputs/laser_wavelength','N'),'lambdaSI','omegaau')
    inverse_GV_IR = InputArchiveCUPRAD['/logs/inverse_group_velocity_SI'][()]; group_velocity_IR = 1./inverse_GV_IR
    # pressure_mbar = 1e3*InputArchiveCUPRAD['/inputs/medium_pressure_in_bar'][()]
    rho0_init = 1e6 * mn.readscalardataset(InputArchiveCUPRAD, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N') # SI
   
    # FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,:,0] + \
    #                     1j*InputArchiveTDSE['FSourceTerm'][:,:,:,1]
    ogrid = InputArchiveTDSE['omegagrid'][:]
    rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
    zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]; Nz_macro = len(zgrid_macro)    
    
    if (kz_end == 'end'): kz_end = Nz_macro  
    
    # Redefine integration, now only if we need last planes
    if redefine_integration:
        
        kz_start = Nz_macro - Nz_last
   


    omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
    ogridSI = omega_au2SI * ogrid
    Hgrid = ogrid/omega0
    H_indices = [mn.FindInterval(Hgrid,H_value) for H_value in Hrange]


    
    # load only the subarray off the data  # FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step]
    FSourceTerm_select = InputArchiveTDSE['FSourceTerm'][0:Nr_max:kr_step,kz_start:kz_end,H_indices[0]:H_indices[1]:ko_step,0] + \
                         1j*InputArchiveTDSE['FSourceTerm'][0:Nr_max:kr_step,kz_start:kz_end,H_indices[0]:H_indices[1]:ko_step,1]

print('data loaded:')
try:
    os.remove(out_h5name)
    print("previous results deleted")
except:
    print("no files deleted")  

rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)
ogrid_select_SI = ogridSI[H_indices[0]:H_indices[1]:ko_step]






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

def L_abs(omega):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 1.0 / (2.0 * rho0_init * units.r_electron_classical * lambdaSI * f2_value) 

if ('dispersion' in apply_diffraction): dispersion_function = dispersion_function_def
else: dispersion_function = None

if ('absorption' in apply_diffraction): absorption_function = absorption_function_def
else: absorption_function = None

## create subintervals to analyse the intensities and phases
Hgrid_I_study = mn.get_odd_interior_points(Hrange)
omega_I_study_intervals = []
omega0SI = mn.ConvertPhoton(omega0, 'omegaau', 'omegaSI')
for k1 in Hgrid_I_study:
    omega_I_study_intervals.append(
        [np.max([ogrid_select_SI[0],omega0SI*(k1-0.5)]), np.min([ogrid_select_SI[-1],omega0SI*(k1+0.5)])]
        )


# Workers = mp.cpu_count() # this is the number of workers


# split the work for the multiprocessing module
rgrid_parts = np.array_split(rgrid_FF,Workers) # r-grids for respective processes
rgrid_indices = [mn.FindInterval(rgrid_FF,r_part[0]) for r_part in rgrid_parts]

# define output
output = mp.Queue()    
def mp_handle(k_pos, *args, **kwargs):
    output.put(
               [k_pos, *Hfn2.HankelTransform_long(*args,**kwargs)] # the outputs are not ordered, keep the order in the result
              )

processes = [mp.Process(target=mp_handle, # define processes
                        args=(k1,
                              ogrid_select_SI,
                              rgrid_macro[0:Nr_max:kr_step],
                              zgrid_macro[kz_start:kz_end],
                              FSourceTerm_select,# FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step],
                              distance_FF,
                              rgrid_parts[k1]
                            ),
                        kwargs={
                              'dispersion_function': dispersion_function,
                              'absorption_function': absorption_function,
                              'frequencies_to_trace_maxima':  omega_I_study_intervals,
                              'store_cummulative_result': store_cummulative_result
                            }
                        
                        ) for k1 in range(Workers)]

for p in processes: p.start(); # run processes

results = [output.get() for p in processes] # there is no ordering!

# collect the results
# The results are unordered, we constructed the multiprocessing, that the first
# index of the result enumerates the process. We then remember the original 
# order and permute the results correctly.
print('collecting results')

FField_FF_integrated = np.empty((len(ogrid_select_SI),Nr_FF), dtype=np.cdouble)
source_maxima=[np.zeros(len(zgrid_macro[kz_start:kz_end])) for k1 in range(len(Hgrid_I_study))]
sim_ind = [result[0] for result in results]

if store_cummulative_result:
    cummulative_field = np.empty(
                        (len(zgrid_macro[kz_start:kz_end])-1, len(ogrid_select_SI),Nr_FF),
                        dtype=np.cdouble)

for k1 in range(Workers):   
   if (sim_ind[k1] == (Workers-1)):
       FField_FF_integrated[:,
                   rgrid_indices[sim_ind[k1]]:
                   ] = results[k1][1]
           
       if store_cummulative_result:
           cummulative_field[:,:,
                   rgrid_indices[sim_ind[k1]]:
                   ] = results[k1][3]
               
   else:
       FField_FF_integrated[:,
                   rgrid_indices[sim_ind[k1]]:rgrid_indices[sim_ind[k1]+1]
                   ] = results[k1][1]
       
       if store_cummulative_result:
           cummulative_field[:,:,
                   rgrid_indices[sim_ind[k1]]:rgrid_indices[sim_ind[k1]+1]
                   ] = results[k1][3]     
       
       
   for k2 in range(len(Hgrid_I_study)):
       for k3 in range(len(zgrid_macro[kz_start:kz_end])):
           source_maxima[k2][k3] = rho0_init * np.max([source_maxima[k2][k3],results[k1][2][k2][k3]])

# adjust density of emitters
cummulative_field = rho0_init * cummulative_field           
FField_FF_integrated = rho0_init * FField_FF_integrated
           
# Diagnostics on the phase           
Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]
Hgrid_phase_diagnostics_indices = [mn.FindInterval(Hgrid_select,H_value) for H_value in Hgrid_I_study]
phase_onax = []; phase_first_plane = []
for k1 in range(len(Hgrid_phase_diagnostics_indices)):
    phase_onax.append(
                    np.unwrap(np.angle(FSourceTerm_select[0,:,Hgrid_phase_diagnostics_indices[k1]]))
                    )
    phase_first_plane.append(
                    np.unwrap(np.angle(FSourceTerm_select[:,0,Hgrid_phase_diagnostics_indices[k1]]))
                    ) 
                     
# Diagnostics of the absorption length
L_abs_Hgrid = np.zeros( (len(ogrid_select_SI),) )
for k1 in range(len(ogrid_select_SI)):     
    L_abs_Hgrid[k1] = L_abs(ogrid_select_SI[k1])                   
           

print('before save')
# Save the data
with h5py.File(out_h5name,'w') as OutFile:
    grp = OutFile.create_group('XUV')
    grp.create_dataset('Spectrum_on_screen',
                                          data = np.stack((FField_FF_integrated.real, FField_FF_integrated.imag),axis=-1)
                                          )
    grp.create_dataset('Maxima_of_planes',
                                          data = np.asarray(source_maxima)
                                          )
    
    grp.create_dataset('Maxima_Hgrid',
                                          data = Hgrid_I_study
                                          )
    grp.create_dataset('rgrid_FF',
                                          data = rgrid_FF
                                          )    
    grp.create_dataset('Hgrid_select',
                                          data = Hgrid_select
                                          )
    grp.create_dataset('zgrid_integration',
                                          data = zgrid_macro[kz_start:kz_end]
                                          )
    grp.create_dataset('Phase_on_axis',
                                          data = np.asarray(phase_onax)
                                          )
    grp.create_dataset('Phase_first_plane',
                                          data = np.asarray(phase_first_plane)
                                          )
    grp.create_dataset('L_abs_Hgrid',
                                          data = L_abs_Hgrid
                                          )
    if store_cummulative_result:
            grp.create_dataset('Spectrum_on_screen_cummulative',
                                          data = np.stack((cummulative_field.real, cummulative_field.imag),axis=-1)
                                          )
    
# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(FField_FF_integrated.T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), integrated, log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# if showplots: plt.show()
# # plt.close(fig)
# # sys.exit()
