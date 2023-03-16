import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
# import sys
# import units
# import mynumerics as mn
# import Hfn
# import Hfn2
import HHG

# import mynumerics as mn
import matplotlib.pyplot as plt
import plot_presets as pp  

# import XUV_refractive_index as XUV_index
# import IR_refractive_index as IR_index

# from scipy import interpolate
# from scipy import integrate

# import XUV_signal_computation as XUV_sig



# fun = lambda x : np.ones(x.shape)

results_TDSE = os.path.join("C:\data", "TDSE_list", "Maker4")
file_TDSE = 'results_merged.h5'
file_TDSE = os.path.join(results_TDSE,file_TDSE)
print('processing:', file_TDSE)             
with h5py.File(file_TDSE, 'r') as InputArchiveTDSE: 
    dum = InputArchiveTDSE['grids_for_scans/varying_params'][:]; Np = len(dum)
    varying_params = []
    for k1 in range(Np): varying_params.append(dum[k1].decode()) 
    
    omega0 = InputArchiveTDSE['grids_for_scans/omega0'][()]
    Ip_TDSE_au = -InputArchiveTDSE['Energy_of_the_ground_state'][()]
    
    pl_step_I = 3
    pl_step_omega = 2
    
    image = pp.figure_driver()    
    image.sf = [pp.plotter() for k1 in range(16)]
    
    image.sf[0].args = [HHG.ComputeCutoff(InputArchiveTDSE[ 'grids_for_scans/param_'+str(varying_params.index('E0')) ][::pl_step_I]**2,omega0, Ip_TDSE_au)[1],
                        InputArchiveTDSE['omegagrid'][::pl_step_omega]/omega0,
                        np.abs(InputArchiveTDSE['FSourceTerm'][::pl_step_I,::pl_step_omega,0] + 1j*InputArchiveTDSE['FSourceTerm'][::pl_step_I,::pl_step_omega,1]).T ]
    image.sf[0].method = plt.pcolormesh
    
    pp.plot_preset(image)
       

    # FSourceTerm =    
    
    # E0_grid = InputArchiveTDSE[ 'grids_for_scans/param_'+str(varying_params.index('E0')) ][:]   
    # ogrid = InputArchiveTDSE['omegagrid'][:]
    # omega0 = InputArchiveTDSE['grids_for_scans/omega0'][()]; omega0SI = mn.ConvertPhoton(omega0, 'omegaau', 'omegaSI')
    # Hgrid = ogrid/omega0
    
   




# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [HHG.ComputeCutoff(E0_grid**2, omega0, Ip_TDSE_au)[1],
#                     Hgrid, np.abs(FSourceTerm).T ]
# image.sf[0].method = plt.pcolormesh

# pp.plot_preset(image)

