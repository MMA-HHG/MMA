import numpy as np
import os
import time
import copy
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn
import re
import glob

import warnings


import matplotlib.pyplot as plt
import plot_presets as pp  

# import XUV_refractive_index as XUV_index
# import IR_refractive_index as IR_index

from matplotlib.lines import Line2D

import HHG
  


arguments = sys.argv



# results_path = [os.path.join("D:\data", "Discharges", "TDSE","scan3"),
#                  os.path.join("D:\data", "Discharges", "TDSE","scan4")]

results_path = os.path.join("D:\data", "TDSE_list", "convergence1")





### load results
results_fnames = [os.path.join("ref", "results_merged.h5"),
                  os.path.join("dx2", "results_merged.h5"),
                  os.path.join("dt2", "results_merged.h5"),
                  os.path.join("2Nx", "results_merged.h5"),
                  os.path.join("2Nxdx2", "results_merged.h5"),
                  os.path.join("2Nxdt2", "results_merged.h5"),
                  os.path.join("4Nx", "results_merged.h5"),
                  os.path.join("4Nxdx2", "results_merged.h5"),
                  os.path.join("4Nxdt2", "results_merged.h5")]

# fname = os.path.join(results_path, results_fname)
tgrid=[]; ogrid=[]; SourceTerm=[]; FSourceTerm=[]; expval_x=[]; PopTot=[];
PopInt=[]; Efields=[]; xgrid_m=[]; ground_state=[]; omega0=[]; E0s=[];
Hgrid=[]
for k1 in range(len(results_fnames)):
    fname = os.path.join(results_path, results_fnames[k1]) 
    with h5py.File(fname,'r') as f:
        omega0.append(f['grids_for_scans/omega0'][()])
        tgrid.append(f['tgrid'][:])
        ogrid.append(f['omegagrid'][:]); Hgrid.append(ogrid[k1]/omega0[k1])
        SourceTerm.append(f['SourceTerm'][:,:])
        FSourceTerm.append(f['FSourceTerm'][:,:,0] + 1j*f['FSourceTerm'][:,:,1])
        expval_x.append(f['expval_x'][:])
        PopTot.append(f['PopTot'][:])
        PopInt.append(f['PopInt'][:])
        Efields.append(f['Efield'][:])
        xgrid_m.append(f['xgrid_micro'][:])
        ground_state.append(f['ground_state'][:,0]+1j*f['ground_state'][:,1])
        
        
        
        E0s.append(f['grids_for_scans/param_1'][:])

Ip = 0.5792




image = pp.figure_driver()    
image.sf = [pp.plotter() for k1 in range(16)]

image.sf[0].args = [Hgrid[3], abs(FSourceTerm[6][15,:])]
image.sf[0].method = plt.semilogy

image.sf[1].args = [Hgrid[4], abs(FSourceTerm[7][15,:])]
image.sf[1].method = plt.semilogy

# image.sf[2].args = [Hgrid[5], abs(FSourceTerm[8][15,:])]
# image.sf[2].method = plt.semilogy

# image.sf[3].args = [Hgrid[6], abs(FSourceTerm[6][15,:])]
# image.sf[3].method = plt.semilogy

# image.sf[4].args = [Hgrid[6], abs(FSourceTerm[6][1,:])]
# image.sf[4].method = plt.semilogy

# image.sf[5].args = [Hgrid[5], abs(FSourceTerm[5][0,:])]
# image.sf[5].method = plt.semilogy


# image.sf[3].args = [Hgrid[3], abs(FSourceTerm[3][15,:])]
# image.sf[3].method = plt.semilogy

# image.sf[1].args = [Hgrid, abs(FSourceTerm[1,:])]
# image.sf[1].method = plt.semilogy

# image.sf[2].args = [Hgrid, abs(FSourceTerm[2,:])]
# image.sf[2].method = plt.semilogy

# image.sf[3].args = [Hgrid, abs(FSourceTerm[12,:])]
# image.sf[3].method = plt.semilogy

# image.sf[4].args = [Hgrid, abs(FSourceTerm[13,:])]
# image.sf[4].method = plt.plot

# image.sf[5].args = [Hgrid, abs(FSourceTerm[14,:])]
# image.sf[5].method = plt.plot

# image.sf[6].args = [Hgrid, abs(FSourceTerm[15,:])]
# image.sf[6].method = plt.plot

pp.plot_preset(image)




# print(HHG.ComputeCutoff(E0s[0]**2, omega0, Ip))

# print(HHG.ComputeCutoff(E0s[3]**2, omega0, Ip))

# print(HHG.ComputeCutoff(E0s[15]**2, omega0, Ip))


# Hgrid = ogrid/omega0
# # source = ST

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(4)]

# image.sf[0].args = [xgrid_m, abs(ground_state)]
# image.sf[0].method = plt.plot

# pp.plot_preset(image)



# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(4)]

# image.sf[0].args = [SourceTerm[0,:]]
# image.sf[0].method = plt.plot

# image.sf[1].args = [SourceTerm[3,:]]
# image.sf[1].method = plt.plot

# image.sf[2].args = [SourceTerm[15,:]]
# image.sf[2].method = plt.plot

# pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [Hgrid, abs(FSourceTerm[0,:])]
# image.sf[0].method = plt.semilogy

# image.sf[1].args = [Hgrid, abs(FSourceTerm[1,:])]
# image.sf[1].method = plt.semilogy

# image.sf[2].args = [Hgrid, abs(FSourceTerm[2,:])]
# image.sf[2].method = plt.semilogy

# image.sf[3].args = [Hgrid, abs(FSourceTerm[12,:])]
# image.sf[3].method = plt.semilogy

# image.sf[4].args = [Hgrid, abs(FSourceTerm[13,:])]
# image.sf[4].method = plt.plot

# image.sf[5].args = [Hgrid, abs(FSourceTerm[14,:])]
# image.sf[5].method = plt.plot

# image.sf[6].args = [Hgrid, abs(FSourceTerm[15,:])]
# image.sf[6].method = plt.plot

# pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(4)]

# image.sf[0].args = [tgrid, PopTot[0,:]]
# image.sf[0].method = plt.plot

# image.sf[1].args = [tgrid, PopTot[3,:]]
# image.sf[1].method = plt.plot

# image.sf[2].args = [tgrid, PopTot[15,:]]
# image.sf[2].method = plt.plot

# pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [tgrid, expval_x[0,:]]
# image.sf[0].method = plt.plot

# image.sf[1].args = [tgrid, expval_x[3,:]]
# image.sf[1].method = plt.plot

# image.sf[2].args = [tgrid, expval_x[13,:]]
# image.sf[2].method = plt.plot

# image.sf[3].args = [tgrid, expval_x[14,:]]
# image.sf[3].method = plt.plot

# image.sf[4].args = [tgrid, expval_x[15,:]]
# image.sf[4].method = plt.plot

# pp.plot_preset(image)

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(7)]

# image.sf[0].args = [tgrid, PopInt[0,:]]
# image.sf[0].method = plt.plot

# image.sf[1].args = [tgrid, PopInt[3,:]]
# image.sf[1].method = plt.plot

# image.sf[2].args = [tgrid, PopInt[15,:]]
# image.sf[2].method = plt.plot

# image.sf[3].args = [tgrid, PopTot[0,:],'--']
# image.sf[3].method = plt.plot

# image.sf[4].args = [tgrid, PopTot[3,:],'--']
# image.sf[4].method = plt.plot

# image.sf[5].args = [tgrid, PopTot[15,:],'--']
# image.sf[5].method = plt.plot

# pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(4)]

# image.sf[0].args = [tgrid, Efields[0,:]]
# image.sf[0].method = plt.plot

# image.sf[1].args = [tgrid, Efields[3,:]]
# image.sf[1].method = plt.plot

# image.sf[2].args = [tgrid, Efields[15,:]]
# image.sf[2].method = plt.plot

# pp.plot_preset(image)








# # store results
# if os.path.exists(OutPath) and os.path.isdir(OutPath):
#   shutil.rmtree(OutPath)
#   print('deleted previous results')
# os.mkdir(OutPath)

# os.chdir(OutPath)


# ## plot H17 + ionisations
# k1 = 1 # index to access given harmonic plot


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k2 in range(11)]



# image.sf[0].args = [p_grid, XUV_energy_pp[0][:,0,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
# image.sf[1].args = [p_grid, XUV_energy_pp[0][:,1,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

# image.sf[2].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge/2'}
# # image.sf[3].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.max(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}



# # A0 = A_norm(np.mean(XUV_energy_pp[0][:,0,k1]),Hgrid_study[k1],omegaSI)

# # image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}



# # image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
# # image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    

# for k2 in range(3,8): image.sf[k2].axis_location = 'right'
# image.sf[3].args = [p_grid, ionisations[0]['half_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
# image.sf[4].args = [p_grid, ionisations[0]['half'], 'b--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
# image.sf[5].args = [p_grid, ionisations[1]['half_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
# image.sf[6].args = [p_grid, ionisations[1]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 

# image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal(Hgrid_study[k1])*np.ones(len(p_grid)), 'g--']
  
# # # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  




# # preions = ionisation_ratio(A0,XUV_energy_pp[0][:,1,k1],Hgrid_study[k1],omegaSI)
# # image.sf[7].args = [p_grid, 100*preions[0][:], 'c--'];
# # image.sf[8].args = [p_grid, 100*preions[1][:], 'c--'];
# # image.sf[7].method = None
# # image.sf[8].method = None

# ## custom legend
# # custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
# #                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
# #                 Line2D([0], [0], color="tab:blue", lw=3),
# #                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
# #                 Line2D([0], [0], color="tab:green", lw=3),                
# #                 Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

# custom_lines = [Line2D([1], [0], color="k"),
#                 Line2D([0], [0], color="tab:grey", linestyle=":"),
#                 Line2D([0], [0], color="b"),
#                 Line2D([0], [0], color="tab:grey", linestyle="--"),
#                 Line2D([0], [0], color="r"),                
#                 Line2D([0], [0], color="g", linestyle="--")]

# # ax.legend(custom_lines, [I0s_leg[0],
# #                          pressures_leg[0],
# #                          I0s_leg[1],
# #                          pressures_leg[1],
# #                          I0s_leg[2],
# #                          pressures_leg[2]],
# #           loc=1, ncol=3)

# image.legend_args = [custom_lines,['no preion.', r'$\eta_0$', '40 A', r'$\eta_{las.}$','50 A', '$\eta_{opt.}$']]
# image.legend_kwargs = {'loc': 1, 'ncol': 3}

# # image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
# image.xlabel = r'$p$ [mbar]'; image.ylabel = r'$I_{\mathrm{XUV}}$ [arb. u.]'; image.right_ylabel = 'ionisation [%]'

# image.title = r'$H_{'+str(Hgrid_study[k1]) + r'}$, $T_{\mathrm{discharge}}/2$'

# image.savefig_args = ['compare1.pdf']
# image.savefig_kwargs = {'dpi' : 600,'bbox_inches' : 'tight'}

# image.set_fontsizes = 'doublet+'

# pp.plot_preset(image)


