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


file_path = os.path.join(results_path,file)
print('processing:', file_path)             
with h5py.File(file_path, 'r') as InputArchiveCUPRAD, h5py.File('hdf5_temp_0000000.h5', 'r') as InputArchiveTDSE:
    # load data
   Efield_CUPRAD = InputArchiveCUPRAD['/outputs/output_field'][:,0,0]
   tgrid_CUPRAD = InputArchiveCUPRAD['/outputs/tgrid'][:]
   
   Efield_TDSE = InputArchiveTDSE['Efield'][:,0]
   tgrid_TDSE = InputArchiveTDSE['tgrid'][:]
            


fig = plt.figure()
plt.plot(tgrid_CUPRAD,Efield_CUPRAD)
# plt.plot(1e3*zgrid_Fluence,1e6*radius_inv_e2)
plt.title('CUPRAD')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_TDSE*units.TIMEau,Efield_TDSE*units.EFIELDau)
plt.title('TDSE')
plt.show()
# plt.close(fig)

fig = plt.figure()
plt.plot(tgrid_CUPRAD-tgrid_CUPRAD[0],Efield_CUPRAD)
plt.plot(tgrid_TDSE*units.TIMEau,Efield_TDSE*units.EFIELDau)
plt.title('CUPRAD')
plt.show()
# plt.close(fig)






