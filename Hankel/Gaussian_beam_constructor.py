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
import HHG

# import mynumerics as mn
import matplotlib.pyplot as plt
import plot_presets as pp  

import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index

from scipy import interpolate



def Gaussian_phase_map(z,r,w0,lambd,n=1.0,vacuum_frame=True):
    k = 2.0*np.pi*n/(lambd)
    zR = np.pi*(w0**2)*n/lambd
    inv_curv_radius = z/(zR**2+z**2)
    phi_G = np.arctan(z/zR)
    if vacuum_frame:
        k_corr = 2.0*np.pi*(n-1.0)/(lambd)
        phase = k_corr*z + 0.5*k*(r**2)*inv_curv_radius - phi_G
    else:
        phase = k*z + 0.5*k*(r**2)*inv_curv_radius - phi_G
    
    return phase


results_TDSE = os.path.join("D:\data", "TDSE_list", "Maker4")
file_TDSE = 'results_merged.h5'
file_TDSE = os.path.join(results_TDSE,file_TDSE)
print('processing:', file_TDSE)             
with h5py.File(file_TDSE, 'r') as InputArchiveTDSE:   
    # FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,0] + \
    #                    1j*InputArchiveTDSE['FSourceTerm'][:,:,1]   
    # dum = InputArchiveTDSE['grids_for_scans/varying_params'][:]; Np = len(dum)
    # varying_params = []
    # for k1 in range(Np): varying_params.append(dum[k1].decode())    
    # E0_grid = InputArchiveTDSE[ 'grids_for_scans/param_'+str(varying_params.index('E0')) ][:]   
    # ogrid = InputArchiveTDSE['omegagrid'][:]
    omega0 = InputArchiveTDSE['grids_for_scans/omega0'][()]
    # Hgrid = ogrid/omega0

w0 = 120e-6 # 120e-6 #25e-6
gas_type = 'Ar'
pressure = 100e-3

susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omega0,'omegaau','lambdaSI'))
n_IR = np.sqrt(1.+pressure * susc_IR)

zgr = np.linspace(-0.1,0.1,3)
rgr = np.linspace(0.,w0,2)

zm, rm = np.meshgrid(zgr,rgr)

phase_map = Gaussian_phase_map(zm,rm,w0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'))