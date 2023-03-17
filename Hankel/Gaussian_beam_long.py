import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import functools
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
from scipy import integrate

import XUV_signal_computation as XUV_sig







# fun = lambda x : np.ones(x.shape)

results_TDSE = os.path.join("C:\data", "TDSE_list", "Maker4")
file_TDSE = 'results_merged.h5'
file_TDSE = os.path.join(results_TDSE,file_TDSE)
print('processing:', file_TDSE)             
with h5py.File(file_TDSE, 'r') as InputArchiveTDSE:   
    FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,0] + \
                        1j*InputArchiveTDSE['FSourceTerm'][:,:,1]   
    dum = InputArchiveTDSE['grids_for_scans/varying_params'][:]; Np = len(dum)
    varying_params = []
    for k1 in range(Np): varying_params.append(dum[k1].decode())    
    E0_grid = InputArchiveTDSE[ 'grids_for_scans/param_'+str(varying_params.index('E0')) ][:]   
    ogrid = InputArchiveTDSE['omegagrid'][:]
    omega0 = InputArchiveTDSE['grids_for_scans/omega0'][()]; omega0SI = mn.ConvertPhoton(omega0, 'omegaau', 'omegaSI')
    Hgrid = ogrid/omega0
    
    Ip_TDSE_au = -InputArchiveTDSE['Energy_of_the_ground_state'][()]

XUV_table_type_absorption = 'Henke' # {Henke, NIST}    
XUV_table_type_dispersion = 'NIST'

w0 = 120e-6 # 120e-6 #25e-6
gas_type = 'Ar'
pressure = 10e-3
E0 = 1.


zR = np.pi*(w0**2)/mn.ConvertPhoton(omega0, 'omegaau', 'lambdaSI')

susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omega0,'omegaau','lambdaSI'))

n_IR = np.sqrt(1.+pressure * susc_IR)

# zgr = np.linspace(0,0.1*zR,1000)
# zgr = np.linspace(0,zR,10000)
zgr = np.linspace(-zR,zR,20000)
rgr = np.linspace(0.,1.5*w0,1000)


Hlimit = [16, 18]




## Field interpolator
domega = ogrid[1]-ogrid[0]
dE0 = E0_grid[1]-E0_grid[0]
k_omega_sel = list(mn.FindInterval(Hgrid,Hlimit))
ogrid_sel = ogrid[k_omega_sel[0]:k_omega_sel[1]]
ogrid_sel_SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')*ogrid_sel
No_sel = len(ogrid_sel)

Hgrid_sel = ogrid_sel/omega0

FSourceTerm_sel = FSourceTerm[:,k_omega_sel[0]:k_omega_sel[1]]
FSourceTerm_interpE0 = interpolate.interp1d( E0_grid, FSourceTerm_sel ,axis=0)



zm, rm = np.meshgrid(zgr,rgr)




Nz = len(zgr)
ogrid = [17*omega0SI]

FSource = np.ones((1,Nz))

sig = Hfn2.Signal_cum_integrator(ogrid, zgr, FSource)

sig_abs = np.abs(sig)**2

image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.sf[0].args = [sig_abs[0,:]]
# image.sf[1].args = [zgr,phase_map_ref[0,:]]
pp.plot_preset(image)









I0_start = 12.5e17/units.INTENSITYau
I0_end = 37.5e17/units.INTENSITYau#E0_grid[-1]**2
# I0_grid = np.linspace(I0_start,E0_grid[-1]**2,N_I0)

I0_grid = E0_grid**2

H_intens = 25

I0_comp = HHG.ComputeInvCutoff(H_intens, omega0, Ip_TDSE_au)

# E0_sel = np.sqrt((0.75*(I0_grid[0]+I0_grid[-1])))
E0_sel = np.sqrt(I0_comp)


## signal including IR phase from the phase factor
def dispersion_function(omega):
    return XUV_index.dispersion_function(omega, pressure, gas_type+'_' + XUV_table_type_dispersion, n_IR=n_IR)
# dispersion_function = functools.partial(XUV_index.dispersion_function, pressure, gas_type+'_' + XUV_table_type_dispersion, n_IR=n_IR) #dispersion_function_disp_def
No = len(ogrid_sel)
dispersion_factor = np.empty(No)
for k1 in range(No):
    dispersion_factor[k1] = ogrid_sel_SI[k1]*dispersion_function(ogrid_sel_SI[k1])   
    # dispersion_factor[k1] = np.pi/2.1777380065358176e-07 # ogrid[k1]*dispersion_function(ogrid[k1])   

factor_e1 = pressure*np.exp(1j*np.outer(zgr,dispersion_factor))


Fsource_long_prof_disp1 =  FSourceTerm_interpE0( E0_sel * 
                                          mn.Gaussian_E0_map(zgr,np.asarray([0]),w0,1.0,
                                                          mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                                                          n=n_IR,
                                                          incl_z_profile = True,
                                                          incl_radial_wz_profile = True))

sig_long_prof_disp1 = Hfn2.Signal_cum_integrator(ogrid_sel, zgr, (factor_e1 * np.conj(Fsource_long_prof_disp1)).T)


## signal including IR phase from the Gaussian beam constructor
# dispersion_function = dispersion_function_def
# dispersion_function = functools.partial(XUV_index.dispersion_function, pressure, gas_type+'_' + XUV_table_type_dispersion)
def dispersion_function(omega):
    return XUV_index.dispersion_function(omega, pressure, gas_type+'_' + XUV_table_type_dispersion)

No = len(ogrid_sel)
dispersion_factor = np.empty(No)
for k1 in range(No):
    dispersion_factor[k1] = ogrid_sel_SI[k1]*dispersion_function(ogrid_sel_SI[k1])   
    # dispersion_factor[k1] = np.pi/2.1777380065358176e-07 # ogrid[k1]*dispersion_function(ogrid[k1])   

factor_e_vac = np.exp(1j*np.outer(zgr,dispersion_factor))


Fsource_long_prof_disp2 =  FSourceTerm_interpE0( E0_sel * 
                                          mn.Gaussian_E0_map(zgr,np.asarray([0]),w0,1.0,
                                                          mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                                                          n=n_IR,
                                                          incl_z_profile = True,
                                                          incl_radial_wz_profile = True))

factor_e_Gauss = np.exp(1j*np.outer((ogrid_sel_SI/omega0SI),
                                    mn.Gaussian_phase_map(
                                        zgr,0,w0,mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                                        n=n_IR,
                                        vacuum_frame=True, incl_curv = False, incl_Gouy = False, incl_lin = True)
                                    ).T
                        )

factor_e_Gauss_geom = np.exp(1j*np.outer((ogrid_sel_SI/omega0SI),
                                    mn.Gaussian_phase_map(
                                        zgr,0,w0,mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                                        n=n_IR,
                                        vacuum_frame=True, incl_curv = False, incl_Gouy = True, incl_lin = True)
                                    ).T
                        )



factor_e2 = pressure*factor_e_vac * factor_e_Gauss
factor_e3 = pressure*factor_e_vac * factor_e_Gauss_geom

sig_long_prof_disp2 = Hfn2.Signal_cum_integrator(ogrid_sel, zgr, (factor_e2 * Fsource_long_prof_disp2).T)
sig_long_prof_disp3 = Hfn2.Signal_cum_integrator(ogrid_sel, zgr, (factor_e3 * Fsource_long_prof_disp2).T)




## compare phases
k_Hsel = mn.FindInterval(Hgrid_sel, 17)
image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'full model phases'
image.sf[0].args = [zgr,np.angle(factor_e1.T[k_Hsel,:])]
image.sf[1].args = [zgr,np.angle(factor_e2.T[k_Hsel,:]), '--']
pp.plot_preset(image)







image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'from TDSE incl. all'
# image.sf[0].args = [zgr,np.abs(sig_long[k_Hsel,:])**2]
# image.sf[1].args = [zgr,np.abs(sig_long_prof[k_Hsel,:])**2]
image.sf[2].args = [zgr,np.abs(sig_long_prof_disp1[k_Hsel,:])**2, '-']
image.sf[3].args = [zgr,np.abs(sig_long_prof_disp2[k_Hsel,:])**2, '-.']
image.sf[4].args = [zgr,np.abs(sig_long_prof_disp3[k_Hsel,:])**2, ':']
pp.plot_preset(image)






































# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [HHG.ComputeCutoff(E0_grid**2, omega0, Ip_TDSE_au)[1],
#                     Hgrid, np.abs(FSourceTerm).T ]
# image.sf[0].method = plt.pcolormesh

# pp.plot_preset(image)


### old stuff 
# included_eff = {'incl_Gouy' : True,
#                 'incl_curv' : True,
#                 'incl_lin'  :True}


# phase_map = Gaussian_phase_map(zm,rm,w0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'),
#                                **included_eff)

# phase_map_ref = Gaussian_phase_map(zm,rm,w0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'),
#                                n = n_IR,
#                                **included_eff)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k2 in range(16)]

# image.sf[0].args = [phase_map[:,0]]
# image.sf[1].args = [phase_map[:,-1]]
# # image.sf[0].method = plt.pcolormesh

# # image.sf[0].colorbar.show = True
# pp.plot_preset(image)


# # image = pp.figure_driver()    
# # image.sf = [pp.plotter() for k2 in range(16)]
# # image.sf[0].args = [phase_map[:,0]]
# # image.sf[1].args = [phase_map_ref[:,0]]
# # pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [zgr,phase_map[0,:]]
# image.sf[1].args = [zgr,phase_map[-1,:]]
# pp.plot_preset(image)


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [zgr,phase_map[0,:]]
# image.sf[1].args = [zgr,phase_map_ref[0,:]]
# pp.plot_preset(image)


# image = pp.figure_driver()   
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [zgr,rgr,phase_map[:,:]]

# image.sf[0].method = plt.pcolormesh

# image.sf[0].colorbar.show = True
# pp.plot_preset(image)


# image = pp.figure_driver()   
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [Gaussian_E0_map(zm,rm,w0,E0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'))]

# image.sf[0].method = plt.pcolormesh

# image.sf[0].colorbar.show = True
# pp.plot_preset(image)

# image = pp.figure_driver()   
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [Gaussian_E0_map(zm,rm,w0,E0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'),incl_z_profile=False)]

# image.sf[0].method = plt.pcolormesh

# image.sf[0].colorbar.show = True
# pp.plot_preset(image)


# image = pp.figure_driver()   
# image.sf = [pp.plotter() for k2 in range(16)]
# image.sf[0].args = [Gaussian_E0_map(zm,rm,w0,E0,mn.ConvertPhoton(omega0,'omegaau','lambdaSI'),incl_radial_wz_profile=False)]

# image.sf[0].method = plt.pcolormesh

# image.sf[0].colorbar.show = True
# pp.plot_preset(image)
