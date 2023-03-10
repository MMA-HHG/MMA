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
from scipy import integrate

import XUV_signal_computation as XUV_sig



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


def Gaussian_E0_map(z,r,w0,E0,lambd,n=1.0,
                        incl_z_profile = True,
                        incl_radial_wz_profile = True):
    
    zR = np.pi*(w0**2)*n/lambd
    w_z = lambda z: w0*np.sqrt(1+(z/zR)**2)
        
    E0_rz = E0
    if incl_z_profile: E0_rz *= (w0/w_z(z))
    
    if incl_radial_wz_profile: E0_rz *= np.exp(-(r/w_z(z))**2)
    else: E0_rz *= np.exp(-(r/w0)**2)
    
    return E0_rz


def Signal_cum_integrator(ogrid, zgrid, FSourceTerm,
                         integrator = integrate.cumulative_trapezoid):
    
    No = len(ogrid); Nz = len(zgrid)
    signal = np.zeros((No,Nz), dtype=np.cdouble)
    integrand = FSourceTerm
    for k1 in range(No):
        signal[k1,1:] = integrator(integrand[k1,:],zgrid)
    return signal

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
print('susIR code:', susc_IR)
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
No_sel = len(ogrid_sel)

Hgrid_sel = ogrid_sel/omega0

FSourceTerm_sel = FSourceTerm[:,k_omega_sel[0]:k_omega_sel[1]]
FSourceTerm_interpE0 = interpolate.interp1d( E0_grid, FSourceTerm_sel ,axis=0)



zm, rm = np.meshgrid(zgr,rgr)




Nz = len(zgr)
ogrid = [17*omega0SI]

FSource = np.ones((1,Nz))

sig = Signal_cum_integrator(ogrid, zgr, FSource)

sig_abs = np.abs(sig)**2

image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.sf[0].args = [sig_abs[0,:]]
# image.sf[1].args = [zgr,phase_map_ref[0,:]]
pp.plot_preset(image)


N_atm = 2.7e19 * 1e6
def f1_funct(E):
    return XUV_index.getf1(gas_type+'_' + XUV_table_type_dispersion, E)

def dispersion_function_vac_def(omega):
    f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    nXUV     = 1.0 - pressure*N_atm*units.r_electron_classical * ((lambdaSI**2)*f1_value/(2.0*np.pi))           
    phase_velocity_XUV  = units.c_light / nXUV
    # print(phase_velocity_XUV)
    # print(((1./units.c_light) - (1./phase_velocity_XUV)))
    return ((1./units.c_light) - (1./phase_velocity_XUV))




def dispersion_function_disp_def(omega):
    f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    nXUV     = 1.0 - pressure*N_atm*units.r_electron_classical * ((lambdaSI**2)*f1_value/(2.0*np.pi))           
    phase_velocity_XUV  = units.c_light / nXUV
    
    phase_velocity_IR = units.c_light / n_IR
    # print(phase_velocity_XUV)
    # print(((1./units.c_light) - (1./phase_velocity_XUV)))
    # print('nXUV, nIR, code:', nXUV, n_IR)
    return ((1./phase_velocity_IR) - (1./phase_velocity_XUV))






dispersion_function = dispersion_function_vac_def
No = len(ogrid)
dispersion_factor = np.empty(No)
for k1 in range(No):
    dispersion_factor[k1] = ogrid[k1]*dispersion_function(ogrid[k1])   
    # dispersion_factor[k1] = np.pi/2.1777380065358176e-07 # ogrid[k1]*dispersion_function(ogrid[k1])   

factor_e = np.exp(1j*np.outer(zgr,dispersion_factor))
sig3 = pressure*Signal_cum_integrator(ogrid, zgr, (factor_e.T) * FSource)



image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'numerical vs. analytic'
image.sf[1].args = [zgr,np.abs(sig3[0,:])**2,'--']
pp.plot_preset(image)



Fsource_IR_mod = np.exp(1j*(ogrid/omega0SI)*Gaussian_phase_map(
                        zgr,0,w0,mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                        n=n_IR,
                        vacuum_frame=True, incl_curv = False, incl_Gouy = False, incl_lin = True))

Fsource_IR_mod2 = np.exp(1j*(ogrid/omega0SI)*Gaussian_phase_map(
                        zgr,0,w0,mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                        n=n_IR,
                        vacuum_frame=True, incl_curv = False, incl_Gouy = True, incl_lin = True))

Fsource_IR_mod3 = np.exp(1j*(ogrid/omega0SI)*Gaussian_phase_map(
                        zgr,0,w0,mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                        n=1.0,
                        vacuum_frame=True, incl_curv = False, incl_Gouy = True, incl_lin = True))


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'phase 1'
image.sf[1].args = [zgr,np.angle(Fsource_IR_mod)]
image.sf[1].args = [zgr,np.angle(Fsource_IR_mod2)]
pp.plot_preset(image)





sig_anal = XUV_sig.compute_S1_abs(pressure, 0.0, 0.0, zgr - zgr[0], 17,
                                  {'omegaSI': omega0SI, #ogrid[0],
                                   'XUV_table_type_absorption': XUV_table_type_absorption,
                                   'XUV_table_type_dispersion': XUV_table_type_dispersion,
                                   'gas_type': gas_type,
                                   'Aq': 1.},
                                   include_absorption=False)



sig_comp = pressure*Signal_cum_integrator(ogrid, zgr, (factor_e.T) * Fsource_IR_mod)
sig_comp2 = pressure*Signal_cum_integrator(ogrid, zgr, (factor_e.T) * Fsource_IR_mod2)
sig_comp3 = pressure*Signal_cum_integrator(ogrid, zgr, np.reshape(Fsource_IR_mod3,(1,len(Fsource_IR_mod3))))

image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'phase fast'
image.sf[1].args = [zgr,np.angle( ((factor_e.T) * Fsource_IR_mod))[0,:] ]
image.sf[2].args = [zgr,np.angle( ((factor_e.T) * Fsource_IR_mod2))[0,:] ]
pp.plot_preset(image)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'numerical vs. analytic'
image.sf[0].args = [zgr,np.abs(sig_anal[0])**2]
image.sf[1].args = [zgr,np.abs(sig_comp[0,:])**2,'--']
image.sf[2].args = [zgr,np.abs(sig_comp2[0,:])**2,':']
# image.sf[3].args = [zgr,np.abs(sig_comp3[0,:])**2,'-.']
pp.plot_preset(image)



I0_start = 12.5e17/units.INTENSITYau
I0_end = 37.5e17/units.INTENSITYau#E0_grid[-1]**2
# I0_grid = np.linspace(I0_start,E0_grid[-1]**2,N_I0)

I0_grid = E0_grid**2

H_intens = 21

I0_comp = HHG.ComputeInvCutoff(H_intens, omega0, Ip_TDSE_au)

# E0_sel = np.sqrt((0.75*(I0_grid[0]+I0_grid[-1])))
E0_sel = np.sqrt(I0_comp)


print('sel cutoff', HHG.ComputeCutoff(E0_sel**2, omega0, Ip_TDSE_au)[1])
Fsource_long =  FSourceTerm_interpE0( E0_sel * np.ones((Nz,)) )
Fsource_long_prof =  FSourceTerm_interpE0( E0_sel * 
                                          Gaussian_E0_map(zgr,np.asarray([0]),w0,1.0,
                                                          mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                                                          n=1.0,
                                                          incl_z_profile = True,
                                                          incl_radial_wz_profile = True))


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'profile'
image.sf[0].args = [zgr,Gaussian_E0_map(zgr,np.asarray([0]),w0,1.0,
                mn.ConvertPhoton(omega0SI, 'omegaSI', 'lambdaSI'),
                n=1.0,
                incl_z_profile = True,
                incl_radial_wz_profile = True)]
pp.plot_preset(image)

image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'sources'
image.sf[0].args = [HHG.ComputeCutoff(E0_grid**2, omega0, Ip_TDSE_au)[1], np.abs(FSourceTerm[:,mn.FindInterval(Hgrid, 21)]) ]
pp.plot_preset(image)




sig_long = pressure*Signal_cum_integrator(ogrid_sel, zgr, Fsource_long.T)
sig_long_prof = pressure*Signal_cum_integrator(ogrid_sel, zgr, Fsource_long_prof.T)

sig_long_prof_phase = pressure*Signal_cum_integrator(ogrid_sel, zgr, Fsource_long_prof.T/np.abs(Fsource_long_prof.T))
sig_long_prof_ampl = pressure*Signal_cum_integrator(ogrid_sel, zgr, np.abs(Fsource_long_prof.T))


k_Hsel = mn.FindInterval(Hgrid_sel, 17)

image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'integrands'
image.sf[0].args = [zgr,np.abs(Fsource_long.T[k_Hsel,:])]
image.sf[1].args = [zgr,np.abs(Fsource_long_prof.T[k_Hsel,:])]
image.sf[2].args = [zgr,np.abs(np.abs(Fsource_long_prof.T)[k_Hsel,:])]
pp.plot_preset(image)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'integrands phase'
image.sf[0].args = [zgr,np.angle(Fsource_long.T[k_Hsel,:])]
image.sf[1].args = [zgr,np.angle(Fsource_long_prof.T[k_Hsel,:])]
image.sf[2].args = [zgr,np.angle(np.abs(Fsource_long_prof.T)[k_Hsel,:])]
pp.plot_preset(image)





image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'from TDSE'
image.sf[0].args = [zgr,np.abs(sig_long[k_Hsel,:])**2]
image.sf[1].args = [zgr,np.abs(sig_long_prof[k_Hsel,:])**2]
image.sf[2].args = [zgr,np.abs(sig_long_prof_ampl[k_Hsel,:])**2]
pp.plot_preset(image)


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]
image.title = 'from TDSE phase'
image.sf[0].args = [zgr,np.abs(sig_long_prof_phase[k_Hsel,:])**2]
pp.plot_preset(image)


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
