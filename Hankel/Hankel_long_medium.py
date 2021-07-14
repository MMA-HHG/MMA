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

import XUV_refractive_index as XUV_index



try:
    with h5py.File('inputs_Hankel.h5', 'r') as Parameters:
        print('reading from hdf5-input file') 
        gas_type = mn.readscalardataset(Parameters, 'inputs/gas_type', 'S')
        #  XUV_table_type = mn.readscalardataset(Parameters, 'inputs/XUV_table_type', 'S') # 'NIST' # {Henke, NIST}        
        XUV_table_type_diffraction = mn.readscalardataset(Parameters,
                                     'inputs/XUV_table_type_diffraction', 'S')
        XUV_table_type_absorption = mn.readscalardataset(Parameters,
                                     'inputs/XUV_table_type_absorption', 'S')  
        
        Nr_max = mn.readscalardataset(Parameters, 'inputs/Nr_max', 'N') # 0.0
        
        Nz_max_sum = mn.readscalardataset(Parameters, 'inputs/Nz_max_sum', 'N')
        
        Hrange = Parameters['inputs/Hrange'][:].tolist() 

        
        kr_step = mn.readscalardataset(Parameters, 'inputs/kr_step', 'N')
        ko_step = mn.readscalardataset(Parameters, 'inputs/ko_step', 'N')
        
        rmax_FF = mn.readscalardataset(Parameters, 'inputs/rmax_FF', 'N')
        Nr_FF = mn.readscalardataset(Parameters, 'inputs/Nr_FF', 'N')
        
        FF_orders_plot = mn.readscalardataset(Parameters, 'inputs/FF_orders_plot', 'N')
        
except:
    print('error reading hdf5 file, using defaults') 
    gas_type = 'Kr'
    XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
    XUV_table_type_absorption = 'Henke' # {Henke, NIST} 
    Nr_max = 235 #470; 235; 155-still fine
    
    Hrange = [17, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]
    
    kr_step = 2 # descending order, tha last is "the most accurate"
    ko_step = 2
    
    rmax_FF = 8*1e-4
    Nr_FF = 200
    
    FF_orders_plot = 4
    
    Nz_max_sum = 5 # 41

    
    

# gas_type = 'Kr'
# XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
# XUV_table_type_absorption = 'Henke' # {Henke, NIST}


def f1_funct(E):
    return XUV_index.getf1(gas_type+'_' + XUV_table_type_diffraction, E)

def f2_funct(E):
    return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)

# f2 = XUV_index.getf2(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(50e-9, 'lambdaSI', 'eV'))
# f1f2 = XUV_index.getf(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(20e-9, 'lambdaSI', 'eV'))


arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
    results_CUPRAD = os.getcwd()
    results_TDSE = os.getcwd()
else:
    results_CUPRAD = os.path.join("D:\data", "Discharges", "TDSE", "t6")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSEH1")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE40planes1")
    results_TDSE = os.path.join("C:\data", "Discharges", "TDSE", "TDSE50planes1") # fine up to 46
    
    # results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE10planes1")


file_CUPRAD = 'results_1.h5'
file_TDSE = 'results_merged_Fsource.h5' # 'hdf5_temp_0000000.h5'


file_CUPRAD = os.path.join(results_CUPRAD,file_CUPRAD)
file_TDSE = os.path.join(results_TDSE,file_TDSE)


# sys.exit()

print('processing:', file_CUPRAD, file_TDSE)             
with h5py.File(file_CUPRAD, 'r') as InputArchiveCUPRAD, h5py.File(file_TDSE, 'r') as InputArchiveTDSE:
   # load data
   omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchiveCUPRAD,'/inputs/laser_wavelength','N'),'lambdaSI','omegaau')
   inverse_GV_IR = InputArchiveCUPRAD['/logs/inverse_group_velocity_SI'][()]; group_velocity_IR = 1./inverse_GV_IR
   # pressure_mbar = 1e3*InputArchiveCUPRAD['/inputs/medium_pressure_in_bar'][()]
   rho0_init = 1e6 * mn.readscalardataset(InputArchiveCUPRAD, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N') # SI
   
   FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,:,0] + \
                       1j*InputArchiveTDSE['FSourceTerm'][:,:,:,1]
   ogrid = InputArchiveTDSE['omegagrid'][:]
   rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
   zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]

   
   # GS_init = InputArchiveTDSE['ground_state'][:,0] + 1j*InputArchiveTDSE['ground_state'][:,1]





print('data loaded:')

out_h5name = 'Hankel.h5'

try:
    os.remove(out_h5name)
    print("previous results deleted")
except:
    print("no files deleted")
    
# sys.exit()

# Nr_max = 235 #470; 235; 155-still fine
# kr_step = 2 # descending order, tha last is "the most accurate"
# ko_step = 2

# rmax_FF = 8*1e-4
# Nr_FF = 200

# FF_orders_plot = 4

omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
ogridSI = omega_au2SI * ogrid


Hgrid = ogrid/omega0
# Hrange = [17, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]
H_indices = [mn.FindInterval(Hgrid,Hvalue) for Hvalue in Hrange]

# Nz_max_sum = 5 # 41



rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)
ogrid_select_SI = ogridSI[H_indices[0]:H_indices[1]:ko_step]



# Here are the fuction to obtain the phase factors in SI units: exp(i*omega*function(omega))
def dispersion_function(omega):
    f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    nXUV     = 1.0 - rho0_init*units.r_electron_classical * \
               ((lambdaSI**2)*f1_value/(2.0*np.pi))           
    phase_velocity_XUV  = units.c_light / nXUV
    return ((1./group_velocity_IR) - (1./phase_velocity_XUV))

def absorption_function(omega):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    beta_factor = rho0_init*units.r_electron_classical * \
                  ((lambdaSI**2)*f2_value/(2.0*np.pi))
    return beta_factor / units.c_light

FField_FF_int_test = Hfn2.HankelTransform_long(ogrid_select_SI,
                                               rgrid_macro[0:Nr_max:kr_step],
                                               zgrid_macro[:Nz_max_sum],
                                               FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step],
                                               0.3,
                                               rgrid_FF)

FField_FF_int_adj_test = Hfn2.HankelTransform_long(ogrid_select_SI,
                                               rgrid_macro[0:Nr_max:kr_step],
                                               zgrid_macro[:Nz_max_sum],
                                               FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step],
                                               0.3,
                                               rgrid_FF,
                                               dispersion_function = dispersion_function)

FField_FF_int_adj_abs_test = Hfn2.HankelTransform_long(ogrid_select_SI,
                                               rgrid_macro[0:Nr_max:kr_step],
                                               zgrid_macro[:Nz_max_sum],
                                               FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step],
                                               0.3,
                                               rgrid_FF,
                                               dispersion_function = dispersion_function,
                                               absorption_function = absorption_function)


 # FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_indices[0]:H_indices[1]:ko_step]).T
         
 #    # FSourceTerm_select[np.newaxis] # to enforce correct dimensions
 #    # FSourceTerm_select.reshape((1,) + FSourceTerm_select.shape) # to enforce correct dimensions
    
   
    
 #    FField_FF = Hfn2.HankelTransform(ogrid_select_SI,
 #                                      rgrid_macro[0:Nr_max:kr_step],
 #                                      FSourceTerm_select,
 #                                      0.3-zgrid_macro[k1],
 #                                      rgrid_FF)


with h5py.File(out_h5name,'w') as OutFile:
    grp = OutFile.create_group('XUV')
    grp.create_dataset('Spectrum_on_screen',
                                          data = np.stack((FField_FF_int_test.real, FField_FF_int_test.imag),axis=-1)
                                          )
    # grp.create_dataset('Spectrum_on_screen_abs',
    #                    data = np.stack((FField_FF_int.real, FField_FF_int.imag),axis=-1)
    #                    )
    
    grp.create_dataset('Spectrum_on_screen_disp',
                       data = np.stack((FField_FF_int_adj_test.real, FField_FF_int_adj_test.imag),axis=-1)
                       ) 

    grp.create_dataset('Spectrum_on_screen_disp+abs',
                       data = np.stack((FField_FF_int_adj_abs_test.real, FField_FF_int_adj_abs_test.imag),axis=-1)
                       )   
    grp.create_dataset('Spectrum_on_screen_test',
                                          data = np.stack((FField_FF_int_test.real, FField_FF_int_test.imag),axis=-1)
                                          )    

Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_int_test.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()


# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_int_adj_test.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated, dispersion, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()


# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_int_adj_abs_test.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated, disp + abs, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_int_adj_abs_test.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated, disp + abs, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()
