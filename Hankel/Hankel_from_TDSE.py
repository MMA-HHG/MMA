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


Nr_max = 470;
kr_step = 10;
ko_step = 2;

omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
ogridSI = omega_au2SI * ogrid

Hgrid = ogrid/omega0
Hrange = [16, 20]
H_indices = [mn.FindInterval(Hgrid,Hvalue) for Hvalue in Hrange]
rgrid_FF = np.linspace(0.0, 1e-4, 100)
ogrid_select_SI = ogridSI[H_indices[0]:H_indices[1]:ko_step]
FSourceTerm_select = np.squeeze(FSourceTerm[0:Nr_max:kr_step,:,H_indices[0]:H_indices[1]:ko_step]).T

FField_FF = Hfn2.HankelTransform(ogrid_select_SI,
                                 rgrid_macro[0:Nr_max:kr_step],
                                 FSourceTerm_select,
                                 0.3,
                                 rgrid_FF)


Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]

# vmin = np.max(np.log(Gaborr))-6.
fig = plt.figure()
plt.pcolor(Hgrid_select,rgrid_FF,abs(FField_FF.T)**2, shading='auto')
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
plt.title('Far-field spectrum')
plt.show()
# plt.close(fig)

# sys.exit()


# vmin = np.max(np.log(Gaborr))-6.
fig = plt.figure()
FF_spectrum_logscale = np.log(abs(FField_FF.T)**2);
vmin = np.max(FF_spectrum_logscale)-6.
plt.pcolor(Hgrid_select,rgrid_FF,np.log(abs(FField_FF.T)**2), shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
plt.title('Far-field spectrum, log')
plt.show()
# plt.close(fig)

# sys.exit()