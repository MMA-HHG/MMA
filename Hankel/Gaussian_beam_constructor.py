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

from scipy import interpolate



def Gaussian_phase_map(r,z,w0,lambd,n=1.0,vacuum_frame=True):
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

