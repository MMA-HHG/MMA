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
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE40planes1")
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE50planes1") # fine up to 46
    
    # results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE10planes1")


file_CUPRAD = 'results_1.h5'
file_TDSE = 'results_merged_Fsource.h5' # 'hdf5_temp_0000000.h5'


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

Nz_max_sum = 41
kz_steps = [8,4,2,1] # descending order, tha last is "the most accurate"

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


FField_FF_z_sum = []
for k1 in range(len(kz_steps)):
    FField_FF_z_sum.append(FField_FF_z[0,0,:])
    for k2 in range(kz_steps[k1],Nz_max_sum,kz_steps[k1]):
        FField_FF_z_sum[k1] = FField_FF_z_sum[k1] + FField_FF_z[k2,0,:]
    FField_FF_z_sum[k1] = FField_FF_z_sum[k1] / len(range(0,Nz_max_sum,kz_steps[k1]))


fig = plt.figure()
plt.plot(rgrid_FF,abs(FField_FF_z_sum[0]))
plt.plot(rgrid_FF,abs(FField_FF_z_sum[1]))
plt.plot(rgrid_FF,abs(FField_FF_z_sum[2]))
plt.plot(rgrid_FF,abs(FField_FF_z_sum[3]))
plt.show()

fig = plt.figure()
dum = np.max(abs(FField_FF_z_sum[-1]))
plt.semilogy(rgrid_FF,abs(FField_FF_z_sum[0]-FField_FF_z_sum[-1])/dum)
plt.semilogy(rgrid_FF,abs(FField_FF_z_sum[1]-FField_FF_z_sum[-1])/dum)
plt.semilogy(rgrid_FF,abs(FField_FF_z_sum[2]-FField_FF_z_sum[-1])/dum)
plt.show()

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
for k1 in range(35):
    plt.plot(rgrid_FF,FField_FF_z[k1,0,:].real)
plt.show()

fig = plt.figure()
for k1 in range(35):
    plt.plot(rgrid_FF,FField_FF_z[k1,0,:].imag)
plt.show()

fig = plt.figure()
for k1 in range(35):
    plt.plot(rgrid_FF,np.unwrap(np.angle(FField_FF_z[k1,0,:])))
plt.show()

fig = plt.figure()
for k1 in range(35):
    plt.plot(rgrid_FF,abs(FField_FF_z[k1,0,:]))
plt.show()
