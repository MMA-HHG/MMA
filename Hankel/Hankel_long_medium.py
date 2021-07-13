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


gas_type = 'Kr'
XUV_table_type_diffraction = 'NIST' # {Henke, NIST}
XUV_table_type_absorption = 'Henke' # {Henke, NIST}


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
   
   
   # SourceTerm_TDSE = InputArchiveTDSE['SourceTerm'][0,0,:]
   FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,:,0] + \
                       1j*InputArchiveTDSE['FSourceTerm'][:,:,:,1]
   ogrid = InputArchiveTDSE['omegagrid'][:]
   rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
   zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]
   # PopTot_TDSE = InputArchiveTDSE['PopTot'][0,0,:]
   # PopInt_TDSE = InputArchiveTDSE['PopInt'][0,0,:]
   # expval_x_TDSE = InputArchiveTDSE['expval_x'][0,0,:]
   
   # GS_init = InputArchiveTDSE['ground_state'][:,0] + 1j*InputArchiveTDSE['ground_state'][:,1]





print('data loaded:')

out_h5name = 'Hankel.h5'

try:
    os.remove(out_h5name)
    print("previous results deleted")
except:
    print("no files deleted")
    
# sys.exit()

Nr_max = 235 #470; 235; 155-still fine
kr_step = 2 # descending order, tha last is "the most accurate"
ko_step = 2

rmax_FF = 8*1e-4
Nr_FF = 200

FF_orders_plot = 4

omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
ogridSI = omega_au2SI * ogrid


Hgrid = ogrid/omega0
Hrange = [17, 18] # [17, 18] # [14, 36] [17, 18] [16, 20] [14, 22]
H_indices = [mn.FindInterval(Hgrid,Hvalue) for Hvalue in Hrange]

Nz_max_sum = 5 # 41
kz_steps = [8,4,2,1] # descending order, tha last is "the most accurate"

# H_index = mn.FindInterval(Hgrid,Hvalue)

rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)
ogrid_select_SI = ogridSI[H_indices[0]:H_indices[1]:ko_step]


# experimenenting with numbers
# phase velocities
omega_SI2energy_eV = mn.ConvertPhoton(1.0, 'omegaSI', 'eV')
omega_SI2lambda_SI = mn.ConvertPhoton(1.0, 'omegaSI', 'lambdaSI')
ogrid_select_eV = ogrid_select_SI * omega_SI2energy_eV
ogrid_select_lambdaSI =  omega_SI2lambda_SI / ogrid_select_SI


f2_values = f2_funct(ogrid_select_eV)
f1_values = f1_funct(ogrid_select_eV)

nXUV = 1.0 - rho0_init*units.r_electron_classical * \
              ((ogrid_select_lambdaSI**2)*f1_values/(2.0*np.pi))
              
            
phase_velocities_XUV = units.c_light / nXUV

dephasing_factor = ogrid_select_SI * ((1./group_velocity_IR) - (1./phase_velocities_XUV))

dephase = np.outer(zgrid_macro,dephasing_factor)

dephase_e = np.exp(1j*dephase)



beta_factor = rho0_init*units.r_electron_classical * \
              ((ogrid_select_lambdaSI**2)*f2_values/(2.0*np.pi))
              
z_exit = zgrid_macro[Nz_max_sum-1]  # zgrid_macro[-1]

absorbing_factor = ogrid_select_SI * beta_factor / units.c_light

absorption = np.outer(zgrid_macro-z_exit ,absorbing_factor)

absorption_e = np.exp(absorption)

absorption_15mm = np.outer(zgrid_macro-15e-3 ,absorbing_factor)

absorption_e_15mm = np.exp(absorption_15mm)


# longitudinal phase

# on axis


Hphase = 17
Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]
k_Hphase = mn.FindInterval(Hgrid, Hphase)
k_Hphase_select = mn.FindInterval(Hgrid_select, Hphase)

fig = plt.figure()
plt.plot(1e3*zgrid_macro[:Nz_max_sum],np.unwrap(np.angle(FSourceTerm[0,:Nz_max_sum,k_Hphase]))) 
plt.title('on-axis phase, group-velocity frame, H17')
plt.xlabel('z [mm]')
plt.ylabel('phi [rad]')
if showplots: plt.show()


dum = np.squeeze(dephase_e[:Nz_max_sum,k_Hphase_select])*np.squeeze(FSourceTerm[0,:Nz_max_sum,k_Hphase])

fig = plt.figure()
plt.plot(1e3*zgrid_macro[:Nz_max_sum],np.unwrap(np.angle(dum[:Nz_max_sum]))) 
plt.title('on-axis phase, XUV frame, H17')
plt.xlabel('z [mm]')
plt.ylabel('phi [rad]')
if showplots: plt.show()

# sys.exit()
# sys.exit()



# include_dispersion = True   
# include_absorption  = True

# FField_FF = []
# FField_FF_scaled = []
for k1 in range(Nz_max_sum):
    print(k1)

    
    # FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
    
    # if (k1 == 0):  # to enforce correct dimensions
    #      dum = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_indices[0]:H_indices[1]:ko_step]).T
    #      FSourceTerm_select = np.zeros((1,)+dum.shape, dtype=np.cdouble)
    #      FSourceTerm_select = dum # np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
    # else:
    FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_indices[0]:H_indices[1]:ko_step]).T
         
    # FSourceTerm_select[np.newaxis] # to enforce correct dimensions
    # FSourceTerm_select.reshape((1,) + FSourceTerm_select.shape) # to enforce correct dimensions
    
   
    
    FField_FF = Hfn2.HankelTransform(ogrid_select_SI,
                                      rgrid_macro[0:Nr_max:kr_step],
                                      FSourceTerm_select,
                                      0.3-zgrid_macro[k1],
                                      rgrid_FF)
    
    if (k1 == 0):
        FField_FF_z = np.zeros( (Nz_max_sum,) + FField_FF.shape,dtype=np.cdouble) 
        FField_FF_z_adj = np.zeros( (Nz_max_sum,) + FField_FF.shape,dtype=np.cdouble)
        FField_FF_z_adj_abs = np.zeros( (Nz_max_sum,) + FField_FF.shape,dtype=np.cdouble)
    
        
    
    FField_FF_z[k1,:,:] = FField_FF    
    FField_FF_z_adj[k1,:,:] = np.outer(dephase_e[k1,:],np.ones(FField_FF.shape[1]))*FField_FF    
    FField_FF_z_adj_abs[k1,:,:] = np.outer(dephase_e[k1,:],np.ones(FField_FF.shape[1])) * \
                                  np.outer(absorption_e[k1,:],np.ones(FField_FF.shape[1])) * \
                                  FField_FF               


# def adjust_plane(k1,plane,include_dispersion,include_absorption,fact1,fact2):
#     pass


# trapezoid
for k1 in range(Nz_max_sum-1):    
    k_step = 1
    
    ## other versions
    if (k1 == 0):
        dum = 0.5*(zgrid_macro[(k1+1)*k_step]-zgrid_macro[k1*k_step]) * \
              (FField_FF_z[k1*k_step,:,:] + FField_FF_z[(k1+1)*k_step,:,:])
              
        dum2 = 0.5*(zgrid_macro[(k1+1)*k_step]-zgrid_macro[k1*k_step]) * \
              (FField_FF_z_adj[k1*k_step,:,:] + FField_FF_z_adj[(k1+1)*k_step,:,:])
              
        dum3 = 0.5*(zgrid_macro[(k1+1)*k_step]-zgrid_macro[k1*k_step]) * \
              (FField_FF_z_adj_abs[k1*k_step,:,:] + FField_FF_z_adj_abs[(k1+1)*k_step,:,:])
    else:
        dum = dum + \
          0.5*(zgrid_macro[(k1+1)*k_step]-zgrid_macro[k1*k_step]) * \
              (FField_FF_z[k1*k_step,:,:] + FField_FF_z[(k1+1)*k_step,:,:])

        dum2 = dum2 + \
          0.5*(zgrid_macro[(k1+1)*k_step]-zgrid_macro[k1*k_step]) * \
              (FField_FF_z_adj[k1*k_step,:,:] + FField_FF_z_adj[(k1+1)*k_step,:,:])
              
        dum3 = dum3 + \
          0.5*(zgrid_macro[(k1+1)*k_step]-zgrid_macro[k1*k_step]) * \
              (FField_FF_z_adj_abs[k1*k_step,:,:] + FField_FF_z_adj_abs[(k1+1)*k_step,:,:])


FField_FF_int = dum
FField_FF_int_adj = dum2
FField_FF_int_adj_abs = dum3

diff_full = (FField_FF_int_adj_abs - FField_FF_int)/np.max(FField_FF_int)
diff_disp = (FField_FF_int_adj - FField_FF_int)/np.max(FField_FF_int)

diff_full_a = (abs(FField_FF_int_adj_abs) - abs(FField_FF_int))/np.max(FField_FF_int)
diff_disp_a = (abs(FField_FF_int_adj) - abs(FField_FF_int))/np.max(FField_FF_int)


## test the function
def dispersion_function(omega):

    f1_values_test = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')

    nXUV_test  = 1.0 - rho0_init*units.r_electron_classical * \
              ((lambdaSI**2)*f1_values_test/(2.0*np.pi))              
            
    phase_velocities_XUV_test  = units.c_light / nXUV_test 

    return ((1./group_velocity_IR) - (1./phase_velocities_XUV_test))

def absorption_function(omega):

    f2_values_test = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')

    beta_factor_test = rho0_init*units.r_electron_classical * \
                  ((lambdaSI**2)*f2_values_test/(2.0*np.pi))
                  

    
    return beta_factor_test / units.c_light

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
                                          data = np.stack((FField_FF_int.real, FField_FF_int.imag),axis=-1)
                                          )
    # grp.create_dataset('Spectrum_on_screen_abs',
    #                    data = np.stack((FField_FF_int.real, FField_FF_int.imag),axis=-1)
    #                    )
    
    grp.create_dataset('Spectrum_on_screen_disp',
                       data = np.stack((FField_FF_int_adj.real, FField_FF_int_adj.imag),axis=-1)
                       ) 

    grp.create_dataset('Spectrum_on_screen_disp+abs',
                       data = np.stack((FField_FF_int_adj_abs.real, FField_FF_int_adj_abs.imag),axis=-1)
                       )   
    grp.create_dataset('Spectrum_on_screen_test',
                                          data = np.stack((FField_FF_int_test.real, FField_FF_int_test.imag),axis=-1)
                                          )    

Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_int.T)**2);
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
FF_spectrum_logscale = np.log10(abs(FField_FF_int_adj.T)**2);
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
FF_spectrum_logscale = np.log10(abs(FField_FF_int_adj_abs.T)**2);
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
FF_spectrum_logscale = np.log10(abs(FField_FF_int_adj_abs.T)**2);
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
FF_spectrum_logscale = np.log10(abs(diff_full.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Error, disp + abs, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(diff_disp.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Error, disp, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()


# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(diff_full_a.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Error, disp + abs, log, a')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(diff_disp_a.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Error, disp, log, a')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
if showplots: plt.show()
# plt.close(fig)
# sys.exit()


# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(FField_FF_z_sum[-1].T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), integrated, log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)
# # sys.exit()

    
# FField_FF_z_sum = [] # trapezoid
# for k1 in range(len(kz_steps)):
#     N_integral = len(range(0,Nz_max_sum,kz_steps[k1])) - 1
#     k_step = kz_steps[k1]
#     # FField_FF_z_sum.append(FField_FF_z[0,:,:])
    
#     for k2 in range(N_integral):
#         print((k2+1)*k_step)
#         if (k2 == 0):
#             dum = 0.5*(zgrid_macro[(k2+1)*k_step]-zgrid_macro[k2*k_step]) * \
#                   (FField_FF_z[k2*k_step,:,:] + FField_FF_z[(k2+1)*k_step,:,:])
#         else:
#             dum = dum + \
#               0.5*(zgrid_macro[(k2+1)*k_step]-zgrid_macro[k2*k_step]) * \
#                   (FField_FF_z[k2*k_step,:,:] + FField_FF_z[(k2+1)*k_step,:,:])
#     # dum = dum / len(range(0,Nz_max_sum,kz_steps[k1]))
#     FField_FF_z_sum.append(dum)


# Hankel_errors = []
# for k1 in range(len(kz_steps)-1):
#     Hankel_errors.append(
#                          (FField_FF_z_sum[k1]-FField_FF_z_sum[-1])/np.max(abs(FField_FF_z_sum[-1]))
#                          )
    
# Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()  
# Hankel_errors_logscale = np.log10(abs(Hankel_errors[0].T))
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,Hankel_errors_logscale, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# dz_string = ', '+"{:.1f}".format(1e6*(zgrid_macro[39]-zgrid_macro[31])) +' mum'
# plt.title('Error' + dz_string)
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()  
# Hankel_errors_logscale = np.log10(abs(Hankel_errors[1].T))
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,Hankel_errors_logscale, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# dz_string = ', '+"{:.1f}".format(1e6*(zgrid_macro[39]-zgrid_macro[35])) +' mum'
# plt.title('Error' + dz_string)
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()  
# Hankel_errors_logscale = np.log10(abs(Hankel_errors[2].T))
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,Hankel_errors_logscale, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# dz_string = ', '+"{:.1f}".format(1e6*(zgrid_macro[39]-zgrid_macro[37])) +' mum'
# plt.title('Error' + dz_string)
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(np.squeeze(FField_FF_z[0,:,:]).T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), entry, log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)


# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(np.squeeze(FField_FF_z[Nz_max_sum-1,:,:]).T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), 40pl, log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(FField_FF_z_sum[-1].T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), integrated, log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)
# # sys.exit()

# # fig = plt.figure()
# # plt.plot(rgrid_FF,abs(FField_FF_z_sum[0]))
# # plt.plot(rgrid_FF,abs(FField_FF_z_sum[1]))
# # plt.plot(rgrid_FF,abs(FField_FF_z_sum[2]))
# # plt.plot(rgrid_FF,abs(FField_FF_z_sum[3]))
# # plt.show()

# # fig = plt.figure()
# # dum = np.max(abs(FField_FF_z_sum[-1]))
# # plt.semilogy(rgrid_FF,abs(FField_FF_z_sum[0]-FField_FF_z_sum[-1])/dum)
# # plt.semilogy(rgrid_FF,abs(FField_FF_z_sum[1]-FField_FF_z_sum[-1])/dum)
# # plt.semilogy(rgrid_FF,abs(FField_FF_z_sum[2]-FField_FF_z_sum[-1])/dum)
# # plt.show()

# # fig = plt.figure()
# # plt.plot(rgrid_FF,FField_FF_z[0,0,:].real)
# # plt.plot(rgrid_FF,FField_FF_z[0,0,:].imag)
# # plt.plot(rgrid_FF,abs(FField_FF_z[0,0,:]))
# # plt.show()

# # fig = plt.figure()
# # plt.plot(rgrid_FF,FField_FF_z[-1,0,:].real)
# # plt.plot(rgrid_FF,FField_FF_z[-1,0,:].imag)
# # plt.plot(rgrid_FF,abs(FField_FF_z[-1,0,:]))
# # plt.show()

# # fig = plt.figure()
# # for k1 in range(35):
# #     plt.plot(rgrid_FF,FField_FF_z[k1,0,:].real)
# # plt.show()

# # fig = plt.figure()
# # for k1 in range(35):
# #     plt.plot(rgrid_FF,FField_FF_z[k1,0,:].imag)
# # plt.show()

# # fig = plt.figure()
# # for k1 in range(35):
# #     plt.plot(rgrid_FF,np.unwrap(np.angle(FField_FF_z[k1,0,:])))
# # plt.show()

# # fig = plt.figure()
# # for k1 in range(35):
# #     plt.plot(rgrid_FF,abs(FField_FF_z[k1,0,:]))
# # plt.show()
