import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
sys.path.append('D:\git\python_modules')
import units
import mynumerics as mn

# import mynumerics as mn
import matplotlib.pyplot as plt


#
arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
    results_path = os.path.join("D:\data", "Discharges")
    results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "foci")
    results_path = os.path.join("D:\data", "Discharges", "f_scan")
    results_path = os.path.join("D:\data", "Discharges", "TDSE", "t6")


file = 'results_1.h5' # 'results_Ar_vac.h5', 'Ar_vac_long.h5' 'results_3.h5' 'results_1.h5'

file_TDSE = 'results_merged_t4.h5' # 'hdf5_temp_0000000.h5'

file_path = os.path.join(results_path,file)
print('processing:', file_path)             
with h5py.File(file_path, 'r') as InputArchiveCUPRAD, h5py.File(file_TDSE, 'r') as InputArchiveTDSE:
    # load data
   Efield_CUPRAD = InputArchiveCUPRAD['/outputs/output_field'][:,0,0]
   tgrid_CUPRAD = InputArchiveCUPRAD['/outputs/tgrid'][:]
   omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchiveCUPRAD,'/inputs/laser_wavelength','N'),'lambdaSI','omegaau')
   
   Efield_TDSE = InputArchiveTDSE['Efield'][0,0,:]
   FEfield_TDSE = InputArchiveTDSE['FEfield'][0,0,:,0] + \
                      1j*InputArchiveTDSE['FEfield'][0,0,:,1]
   tgrid_TDSE = InputArchiveTDSE['tgrid'][:]
   
   SourceTerm_TDSE = InputArchiveTDSE['SourceTerm'][0,0,:]
   FSourceTerm_TDSE = InputArchiveTDSE['FSourceTerm'][0,0,:,0] + \
                      1j*InputArchiveTDSE['FSourceTerm'][0,0,:,1]
   ogrid_TDSE = InputArchiveTDSE['omegagrid'][:]
   PopTot_TDSE = InputArchiveTDSE['PopTot'][0,0,:]
   PopInt_TDSE = InputArchiveTDSE['PopInt'][0,0,:]
   expval_x_TDSE = InputArchiveTDSE['expval_x'][0,0,:]
   
   GS_init = InputArchiveTDSE['ground_state'][:,0] + 1j*InputArchiveTDSE['ground_state'][:,1]
   xgrid_micro = InputArchiveTDSE['xgrid_micro'][:]


fig = plt.figure()
plt.plot(tgrid_CUPRAD,Efield_CUPRAD/units.EFIELDau)
# plt.plot(1e3*zgrid_Fluence,1e6*radius_inv_e2)
plt.title('CUPRAD')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_TDSE*units.TIMEau,Efield_TDSE)
plt.title('TDSE')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_CUPRAD-tgrid_CUPRAD[0],Efield_CUPRAD)
plt.plot(tgrid_TDSE*units.TIMEau,Efield_TDSE*units.EFIELDau)
plt.title('TDSE, CUPRAD')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_TDSE*units.TIMEau,SourceTerm_TDSE)
plt.title('SourceTerm')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_TDSE,
         mn.apply_filter(SourceTerm_TDSE, mn.filter_box, tgrid_TDSE, [6750,tgrid_TDSE[-1]])
         )
plt.title('SourceTerm, filtered')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(ogrid_TDSE/omega0,abs(FSourceTerm_TDSE))
plt.title('FSourceTerm')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.semilogy(ogrid_TDSE/omega0,abs(FSourceTerm_TDSE))
plt.title('FSourceTerm')
plt.show()
# plt.close(fig)

# filter spectrum
ogrid_FE, FE_filter, Nt = mn.fft_t_nonorm(
                            tgrid_TDSE,
                            mn.apply_filter(SourceTerm_TDSE, mn.filter_box, tgrid_TDSE, [6750,tgrid_TDSE[-1]])
                            )

fig = plt.figure()
plt.semilogy(ogrid_FE/omega0,abs(FE_filter))
plt.title('FS_filt')
plt.show()
# plt.close(fig)


fig = plt.figure()
plt.semilogy(ogrid_TDSE/omega0,abs(FEfield_TDSE))
plt.title('FEfield')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_TDSE,PopTot_TDSE)
plt.title('PopTot')
plt.show()

fig = plt.figure()
plt.plot(tgrid_TDSE,PopInt_TDSE)
plt.title('PopInt')
plt.show()

fig = plt.figure()
plt.plot(tgrid_TDSE,expval_x_TDSE)
plt.title('<x>')
plt.show()
# plt.close(fig)

# fig = plt.figure()
# plt.semilogy(tgrid_TDSE,PopTot_TDSE)
# plt.title('PopTot')
# plt.show()
# # plt.close(fig)

fig = plt.figure()
plt.plot(xgrid_micro,abs(GS_init))
plt.title('GS')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.semilogy(xgrid_micro,abs(GS_init))
plt.title('GS')
plt.show()
# plt.close(fig)





