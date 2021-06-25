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


# f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(q*res.omega0, 'omegaSI', 'eV'))
f1 = XUV_index.getf1(gas_type+'_' + XUV_table_type_diffraction, mn.ConvertPhoton(20e-9, 'lambdaSI', 'eV'))

f_abs = XUV_index.getf2(gas_type+'_' + XUV_table_type_absorption, mn.ConvertPhoton(20e-9, 'lambdaSI', 'eV'))

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
else:
    results_CUPRAD = os.path.join("D:\data", "Discharges", "TDSE", "t6")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSEH1")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE40planes1")
    results_TDSE = os.path.join("C:\data", "Discharges", "TDSE", "TDSE50planes1") # fine up to 46
    
    results_TDSE = os.path.join("G:\data_ELI", "Discharges", "TDSE", "TDSE_2mm") # fine up to 46
    
    # results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE10planes1")


file_CUPRAD = 'results_1.h5'
file_TDSE = 'results_merged_Fsource.h5' # 'hdf5_temp_0000000.h5'
file_TDSE = 'results_merged.h5' # 'hdf5_temp_0000000.h5'


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
   
   tgrid_CUPRAD = InputArchiveCUPRAD['outputs/tgrid'][:]
   
   # SourceTerm_TDSE = InputArchiveTDSE['SourceTerm'][0,0,:]
   # FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,:,0] + \
   #                    1j*InputArchiveTDSE['FSourceTerm'][:,:,:,1]
   tgrid = InputArchiveTDSE['tgrid'][:]
   rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
   zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]
   Plasma = InputArchiveTDSE['PopTot'][:]
    
   # PopTot_TDSE = InputArchiveTDSE['PopTot'][0,0,:]
   # PopInt_TDSE = InputArchiveTDSE['PopInt'][0,0,:]
   # expval_x_TDSE = InputArchiveTDSE['expval_x'][0,0,:]
   
   # GS_init = InputArchiveTDSE['ground_state'][:,0] + 1j*InputArchiveTDSE['ground_state'][:,1]

print('data loaded:')

Delta_t_TDSE = tgrid[-1] - tgrid[0]
Delta_t_CUPRAD = tgrid_CUPRAD[-1] - tgrid_CUPRAD[0]

a1 = Delta_t_CUPRAD/Delta_t_TDSE
b1 = tgrid_CUPRAD[-1] - a1*tgrid[-1]

tgrid_TDSE_resc = a1*tgrid + b1

index_t0 = mn.FindInterval(tgrid_TDSE_resc, 0.0e-15)




plasma_map = Plasma[:,:,index_t0]



# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
map1 = ax.pcolor(zgrid_macro, rgrid_macro, 100*(1.0-plasma_map), shading='auto', cmap='plasma')
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('plasma, t=0 fs')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
plt.show()
# plt.close(fig)


plasma_map = Plasma[:,:,-1]
# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
map1 = ax.pcolor(zgrid_macro, rgrid_macro, 100*(1.0-plasma_map), shading='auto', cmap='plasma')
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('plasma, t_end')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
plt.show()
# plt.close(fig)







# FField_FF = []
# for k1 in range(Nz_max_sum):
    
#     # FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
    
#     # if (k1 == 0):  # to enforce correct dimensions
#     #      dum = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_indices[0]:H_indices[1]:ko_step]).T
#     #      FSourceTerm_select = np.zeros((1,)+dum.shape, dtype=np.cdouble)
#     #      FSourceTerm_select = dum # np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
#     # else:
#     FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_indices[0]:H_indices[1]:ko_step]).T
         
#     # FSourceTerm_select[np.newaxis] # to enforce correct dimensions
#     # FSourceTerm_select.reshape((1,) + FSourceTerm_select.shape) # to enforce correct dimensions
    
   
    
#     FField_FF = Hfn2.HankelTransform(ogrid_select_SI,
#                                      rgrid_macro[0:Nr_max:kr_step],
#                                      FSourceTerm_select,
#                                      0.3-zgrid_macro[k1],
#                                      rgrid_FF)
    
#     if (k1 == 0): FField_FF_z = np.zeros( (Nz_max_sum,) + FField_FF.shape,dtype=np.cdouble)  
#     FField_FF_z[k1,:,:] = FField_FF                    



# # FField_FF_z_sum = []
# # for k1 in range(len(kz_steps)):
# #     FField_FF_z_sum.append(FField_FF_z[0,:,:])
# #     for k2 in range(kz_steps[k1],Nz_max_sum,kz_steps[k1]):
# #         FField_FF_z_sum[k1] = FField_FF_z_sum[k1] + FField_FF_z[k2,:,:]
# #     FField_FF_z_sum[k1] = FField_FF_z_sum[k1] / len(range(0,Nz_max_sum,kz_steps[k1]))
    
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
