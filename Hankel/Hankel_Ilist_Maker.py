import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn
import Hfn
import Hfn2

# import mynumerics as mn
import matplotlib.pyplot as plt
import plot_presets as pp  

import XUV_refractive_index as XUV_index

from scipy import interpolate


# # inputs from hdf5-input
# try:
#     with h5py.File('inputs_Hankel.h5', 'r') as Parameters:
#         print('reading from hdf5-input file') 
#         file_CUPRAD = mn.readscalardataset(Parameters, 'inputs/file_CUPRAD', 'S')  
#         file_TDSE = mn.readscalardataset(Parameters, 'inputs/file_TDSE', 'S')  
#         out_h5name = mn.readscalardataset(Parameters, 'inputs/out_h5name', 'S')  
        
#         gas_type = mn.readscalardataset(Parameters, 'inputs/gas_type', 'S')    
#         XUV_table_type_diffraction = mn.readscalardataset(Parameters,
#                                      'inputs/XUV_table_type_diffraction', 'S')
#         XUV_table_type_absorption = mn.readscalardataset(Parameters,
#                                      'inputs/XUV_table_type_absorption', 'S') 
#         apply_diffraction = Parameters['inputs/apply_diffraction'][:].tolist()
#         apply_diffraction = [k1.decode() for k1 in apply_diffraction]
        
        
#         kr_step = mn.readscalardataset(Parameters, 'inputs/kr_step', 'N')
#         ko_step = mn.readscalardataset(Parameters, 'inputs/ko_step', 'N')       
#         Nr_max = mn.readscalardataset(Parameters, 'inputs/Nr_max', 'N')        
#         Nz_max_sum = mn.readscalardataset(Parameters, 'inputs/Nz_max_sum', 'N')        
#         Hrange = Parameters['inputs/Hrange'][:].tolist()     

#         Nr_FF = mn.readscalardataset(Parameters, 'inputs/Nr_FF', 'N')
#         rmax_FF = mn.readscalardataset(Parameters, 'inputs/rmax_FF', 'N')
#         distance_FF = mn.readscalardataset(Parameters, 'inputs/distance_FF', 'N')
        
        
#         FF_orders_plot = mn.readscalardataset(Parameters, 'inputs/FF_orders_plot', 'N')
        

        
# except:
#     print('error reading hdf5 file, using defaults') 
#     gas_type = 'Kr'
#     XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
#     XUV_table_type_absorption = 'Henke' # {Henke, NIST} 
#     apply_diffraction = ['dispersion', 'absorption']
    
#     Nr_max = 235 #470; 235; 155-still fine    
#     Hrange = [16, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]
    
#     kr_step = 2 # descending order, the last is "the most accurate"
#     ko_step = 2
    
#     rmax_FF = 8*1e-4
#     Nr_FF = 200 # 10 # 200
#     distance_FF = 0.3
    
#     FF_orders_plot = 4    
#     Nz_max_sum = 5 # 41
    
#     file_CUPRAD = 'results_1.h5'
#     file_TDSE = 'results_merged_Fsource.h5'
#     out_h5name = 'Hankel.h5'


arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
    results_CUPRAD = os.getcwd()
    results_TDSE = os.getcwd()
else:
    results_CUPRAD = os.path.join("D:\data", "Discharges", "TDSE", "t6")
    # results_CUPRAD = os.path.join("/mnt","d","data", "Discharges", "TDSE", "t6")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSEH1")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE40planes1")
    results_TDSE = os.path.join("C:\data", "Discharges", "TDSE", "TDSE50planes1") # fine up to 46  
    # results_TDSE = os.path.join("/mnt","c","data", "Discharges", "TDSE", "TDSE50planes1")
    



results_TDSE = os.path.join("D:\data", "TDSE_list", "Maker")
file_TDSE = 'results_merged.h5'

file_TDSE = os.path.join(results_TDSE,file_TDSE)


# load data
print('processing:', file_TDSE)             
with h5py.File(file_TDSE, 'r') as InputArchiveTDSE:

   
    FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,0] + \
                       1j*InputArchiveTDSE['FSourceTerm'][:,:,1]
   
    dum = InputArchiveTDSE['grids_for_scans/varying_params'][:]; Np = len(dum)
    varying_params = []
    for k1 in range(Np): varying_params.append(dum[k1].decode())
    
    E0_grid = InputArchiveTDSE[ 'grids_for_scans/param_'+str(varying_params.index('E0')) ][:]
             
    # E0grid =  InputArchiveTDSE['grids_for_scans'][:]
   
    ogrid = InputArchiveTDSE['omegagrid'][:]
    omega0 = InputArchiveTDSE['grids_for_scans/omega0'][()]
    Hgrid = ogrid/omega0
    # rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
    # zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]

print('data loaded:')

## Laser
NI0 = 3
I0_grid = np.linspace(0,E0_grid[-1]**2,NI0)
w0 = 25e-6

Gaussian_E_r = lambda r : np.exp(-(r/w0)**2)

Nr = 100
rgrid = np.linspace(0, 1.5*w0, Nr)


Hlimit = [10, 36]
# Hlimit = [24, 26]




## construct fields

domega = ogrid[1]-ogrid[0]
dE0 = E0_grid[1]-E0_grid[0]

k_omega_sel = list(mn.FindInterval(Hgrid,Hlimit))

ogrid_sel = ogrid[k_omega_sel[0]:k_omega_sel[1]]
No_sel = len(ogrid_sel)

FSourceTerm_sel = FSourceTerm[:,k_omega_sel[0]:k_omega_sel[1]]
FSourceTerm_interpE0 = interpolate.interp1d( E0_grid, FSourceTerm_sel ,axis=0)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k1 in range(16)]

image.sf[0].args = [ogrid_sel/omega0, np.real(FSourceTerm_sel[2500,:])]
image.sf[0].method = plt.semilogy

image.sf[1].args = [ogrid_sel/omega0, np.real(FSourceTerm_sel[2501,:])]
image.sf[1].method = plt.semilogy

image.sf[2].args = [ogrid_sel/omega0, np.real(FSourceTerm_interpE0(E0_grid[2500]+dE0/2.0))]
image.sf[2].method = plt.semilogy

# image.sf[1].args = [Hgrid[4], abs(FSourceTerm[7][15,:])]
# image.sf[1].method = plt.semilogy


pp.plot_preset(image)


# def Efield_r(r):
#     E0_r = Gaussian_E_r(r) # corresponding peak intensity




# omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
# ogridSI = omega_au2SI * ogrid
# Hgrid = ogrid/omega0
# H_indices = [mn.FindInterval(Hgrid,Hvalue) for Hvalue in Hrange]

# try:
#     os.remove(out_h5name)
#     print("previous results deleted")
# except:
#     print("no files deleted")  

# rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)
# ogrid_select_SI = ogridSI[H_indices[0]:H_indices[1]:ko_step]



# # Here are the fuction to obtain the phase factors in SI units: exp(i*omega*function(omega))
# def f1_funct(E):
#     return XUV_index.getf1(gas_type+'_' + XUV_table_type_diffraction, E)
# def f2_funct(E):
#     return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)


# def dispersion_function_def(omega):
#     f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
#     lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
#     nXUV     = 1.0 - rho0_init*units.r_electron_classical * \
#                ((lambdaSI**2)*f1_value/(2.0*np.pi))           
#     phase_velocity_XUV  = units.c_light / nXUV
#     return ((1./group_velocity_IR) - (1./phase_velocity_XUV))

# def absorption_function_def(omega):
#     f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
#     lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
#     beta_factor = rho0_init*units.r_electron_classical * \
#                   ((lambdaSI**2)*f2_value/(2.0*np.pi))
#     return beta_factor / units.c_light

# if ('dispersion' in apply_diffraction): dispersion_function = dispersion_function_def
# else: dispersion_function = None

# if ('absorption' in apply_diffraction): absorption_function = absorption_function_def
# else: dispersion_function = None

# ## create subintervals to analyse the intensities
# Hgrid_I_study = mn.get_odd_interior_points(Hrange)
# omega_I_study_intervals = []
# omega0SI = mn.ConvertPhoton(omega0, 'omegaau', 'omegaSI')
# for k1 in Hgrid_I_study:
#     omega_I_study_intervals.append(
#         [np.max([ogrid_select_SI[0],omega0SI*(k1-0.5)]), np.min([ogrid_select_SI[-1],omega0SI*(k1+0.5)])]
#         )

# # The main integration
# FField_FF_integrated, source_maxima = Hfn2.HankelTransform_long(
#                                                ogrid_select_SI,
#                                                rgrid_macro[0:Nr_max:kr_step],
#                                                zgrid_macro[:Nz_max_sum],
#                                                FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step],
#                                                distance_FF,
#                                                rgrid_FF,
#                                                dispersion_function = dispersion_function, # None, #dispersion_function,
#                                                absorption_function = absorption_function,
#                                                frequencies_to_trace_maxima = omega_I_study_intervals)


# # Save the data
# Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]
# with h5py.File(out_h5name,'w') as OutFile:
#     grp = OutFile.create_group('XUV')
#     grp.create_dataset('Spectrum_on_screen',
#                                           data = np.stack((FField_FF_integrated.real, FField_FF_integrated.imag),axis=-1)
#                                           )
#     grp.create_dataset('Maxima_of_planes',
#                                           data = np.asarray(source_maxima)
#                                           )
#     grp.create_dataset('Maxima_Hgrid',
#                                           data = Hgrid_I_study
#                                           )
#     grp.create_dataset('rgrid_FF',
#                                           data = rgrid_FF
#                                           )    
#     grp.create_dataset('Hgrid_select',
#                                           data = Hgrid_select
#                                           )
        


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
