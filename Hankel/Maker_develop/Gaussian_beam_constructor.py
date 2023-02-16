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



def Gaussian_phase_map(z,r,w0,lambd,n=1.0,vacuum_frame=True,
                       incl_curv = True, incl_Gouy = True,
                       incl_lin = True):
    k = 2.0*np.pi*n/(lambd)
    zR = np.pi*(w0**2)*n/lambd
    inv_curv_radius = z/(zR**2+z**2)
    phi_G = np.arctan(z/zR)
    if vacuum_frame: k_corr = 2.0*np.pi*(n-1.0)/(lambd)
    else: k_corr = k
        
    phase = 0.
    if incl_lin: phase += k_corr*z
    if incl_curv: phase += 0.5*k*(r**2)*inv_curv_radius 
    if incl_Gouy: phase += -phi_G
    
    return phase


def Gaussian_E0_map(z,r,w0,lambd,n=1.0,vacuum_frame=True,
                       incl_curv = True, incl_Gouy = True,
                       incl_lin = True):
    k = 2.0*np.pi*n/(lambd)
    zR = np.pi*(w0**2)*n/lambd
    inv_curv_radius = z/(zR**2+z**2)
    phi_G = np.arctan(z/zR)
    if vacuum_frame: k_corr = 2.0*np.pi*(n-1.0)/(lambd)
    else: k_corr = k
        
    phase = 0.
    if incl_lin: phase += k_corr*z
    if incl_curv: phase += 0.5*k*(r**2)*inv_curv_radius 
    if incl_Gouy: phase += -phi_G
    
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


zR = np.pi*(w0**2)/mn.ConvertPhoton(omega0, 'omegaau', 'lambdaSI')

susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omega0,'omegaau','lambdaSI'))
n_IR = np.sqrt(1.+pressure * susc_IR)

zgr = np.linspace(-zR,zR,2000)
rgr = np.linspace(0.,1.5*w0,1000)

zm, rm = np.meshgrid(zgr,rgr)





included_eff = {'incl_Gouy' : True,
                'incl_curv' : True,
                'incl_lin'  :True}


phase_map = Gaussian_phase_map(zm,rm,w0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'),
                               **included_eff)

phase_map_ref = Gaussian_phase_map(zm,rm,w0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'),
                               n = n_IR,
                               **included_eff)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]

image.sf[0].args = [phase_map[:,0]]
image.sf[1].args = [phase_map[:,-1]]
# image.sf[0].method = plt.pcolormesh

# image.sf[0].colorbar.show = True
pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [phase_map[:,0]]
# image.sf[1].args = [phase_map_ref[:,0]]
# pp.plot_preset(image)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.sf[0].args = [zgr,phase_map[0,:]]
image.sf[1].args = [zgr,phase_map[-1,:]]
pp.plot_preset(image)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.sf[0].args = [zgr,phase_map[0,:]]
image.sf[1].args = [zgr,phase_map_ref[0,:]]
pp.plot_preset(image)


image = pp.figure_driver()   
image.sf = [pp.plotter() for k2 in range(16)]
image.sf[0].args = [phase_map[:,:]]

image.sf[0].method = plt.pcolormesh

image.sf[0].colorbar.show = True
pp.plot_preset(image)