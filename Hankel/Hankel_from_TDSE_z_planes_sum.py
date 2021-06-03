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




arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    results_CUPRAD = os.path.join("D:\data", "Discharges", "TDSE", "t6")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSEH1")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE10planes4")
    
    # results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE10planes1")


file_CUPRAD = 'results_1.h5'
file_TDSE = 'results_merged.h5' # 'hdf5_temp_0000000.h5'


file_CUPRAD = os.path.join(results_CUPRAD,file_CUPRAD)
file_TDSE = os.path.join(results_TDSE,file_TDSE)


print('processing:', file_CUPRAD, file_TDSE)             
with h5py.File(file_CUPRAD, 'r') as InputArchiveCUPRAD, h5py.File(file_TDSE, 'r') as InputArchiveTDSE:
    # load data
   omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchiveCUPRAD,'/inputs/laser_wavelength','N'),'lambdaSI','omegaau')
   
   
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


Nr_max = 235 #470; 235; 155-still fine
kr_step = 2 # descending order, tha last is "the most accurate"
ko_step = 2

rmax_FF = 8*1e-4
Nr_FF = 800

FF_orders_plot = 10

omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
ogridSI = omega_au2SI * ogrid

Hgrid = ogrid/omega0
Hvalue = 17 # [14, 36]


H_index = mn.FindInterval(Hgrid,Hvalue)

rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)
ogrid_select_SI = [ogridSI[H_index]]


FField_FF = []
for k1 in range(len(zgrid_macro)):
    
    # FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
    
    if (k1 == 0):  # to enforce correct dimensions
         dum = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
         FSourceTerm_select = np.zeros((1,)+dum.shape, dtype=np.cdouble)
         FSourceTerm_select[0,:] = dum # np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
    else:
         FSourceTerm_select[0,:] = np.squeeze(FSourceTerm[0:Nr_max:kr_step,k1,H_index]).T
         
    # FSourceTerm_select[np.newaxis] # to enforce correct dimensions
    # FSourceTerm_select.reshape((1,) + FSourceTerm_select.shape) # to enforce correct dimensions
    
   
    
    FField_FF = Hfn2.HankelTransform(ogrid_select_SI,
                                     rgrid_macro[0:Nr_max:kr_step],
                                     FSourceTerm_select,
                                     0.3,
                                     rgrid_FF)
    
    if (k1 == 0): FField_FF_z = np.zeros( (len(zgrid_macro),) + FField_FF.shape,dtype=np.cdouble)  
    FField_FF_z[k1,:,:] = FField_FF                    



fig = plt.figure()
plt.plot(rgrid_FF,FField_FF_z[0,0,:].real)
plt.plot(rgrid_FF,FField_FF_z[0,0,:].imag)
plt.plot(rgrid_FF,abs(FField_FF_z[0,0,:]))
plt.show()

fig = plt.figure()
plt.plot(rgrid_FF,FField_FF_z[-1,0,:].real)
plt.plot(rgrid_FF,FField_FF_z[-1,0,:].imag)
plt.plot(rgrid_FF,abs(FField_FF_z[-1,0,:]))
plt.show()

fig = plt.figure()
for k1 in range(6):
    plt.plot(rgrid_FF,FField_FF_z[k1,0,:].real)
plt.show()

fig = plt.figure()
for k1 in range(6):
    plt.plot(rgrid_FF,FField_FF_z[k1,0,:].imag)
plt.show()

fig = plt.figure()
for k1 in range(6):
    plt.plot(rgrid_FF,abs(FField_FF_z[k1,0,:]))
plt.show()




# Hankel_errors = []
# for k1 in range(len(kr_steps)-1):
#     Hankel_errors.append(
#                          (FField_FF[k1+1]-FField_FF[k1])/np.max(abs(FField_FF[k1]))
#                          )

# Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]

# # vmin = np.max(np.log(Gaborr))-6.
# fig = plt.figure()
# plt.pcolor(Hgrid_select,rgrid_FF,abs(FField_FF.T)**2, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# plt.title('Far-field spectrum')
# plt.show()
# # plt.close(fig)

# # sys.exit()


# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log(abs(FField_FF[-1].T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)


# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()  
# Hankel_errors_logscale = np.log(abs(Hankel_errors[0].T))
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,Hankel_errors_logscale, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Error')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()  
# Hankel_errors_logscale = np.log(abs(Hankel_errors[1].T))
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,Hankel_errors_logscale, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Error')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()  
# Hankel_errors_logscale = np.log(abs(Hankel_errors[2].T))
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,Hankel_errors_logscale, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Error')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # plt.close(fig)


# # sys.exit()