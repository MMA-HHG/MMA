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

file_TDSE = 'results_merged.h5' # 'hdf5_temp_0000000.h5'

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
   ogrid_TDSE = InputArchiveTDSE['/omegagrid'][:]


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
plt.plot(ogrid_TDSE/omega0,abs(FSourceTerm_TDSE))
plt.title('FSourceTerm')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.semilogy(ogrid_TDSE/omega0,abs(FSourceTerm_TDSE))
plt.title('FSourceTerm')
plt.show()
# plt.close(fig)


fig = plt.figure()
plt.semilogy(ogrid_TDSE/omega0,abs(FEfield_TDSE))
plt.title('FEfield')
plt.show()
# plt.close(fig)




