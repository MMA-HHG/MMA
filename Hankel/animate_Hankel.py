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
from matplotlib import animation
import plot_presets as pp  

import XUV_refractive_index as XUV_index

from scipy import interpolate







file_Hankel = 'Hankel_M5.h5'


# load data
          
with h5py.File(file_Hankel, 'r') as Hankel_results:
    HHG_onscreen = Hankel_results['XUV/Spectrum_on_screen'][:,:,:,0] + 1j*Hankel_results['XUV/Spectrum_on_screen'][:,:,:,1]
    I0_grid = Hankel_results['XUV/I0_grid'][:]
    rgrid_FF = Hankel_results['XUV/rgrid_FF'][:]
    Hgrid = Hankel_results['XUV/Hgrid_sel'][:]
   
    # FSourceTerm = InputArchiveTDSE['FSourceTerm'][:,:,0] + \
    #                    1j*InputArchiveTDSE['FSourceTerm'][:,:,1]
   
    # dum = InputArchiveTDSE['grids_for_scans/varying_params'][:]; Np = len(dum)
    # varying_params = []
    # for k1 in range(Np): varying_params.append(dum[k1].decode())
    
    # E0_grid = InputArchiveTDSE[ 'grids_for_scans/param_'+str(varying_params.index('E0')) ][:]
             
    # # E0grid =  InputArchiveTDSE['grids_for_scans'][:]
   
    # ogrid = InputArchiveTDSE['omegagrid'][:]
    # omega0 = InputArchiveTDSE['grids_for_scans/omega0'][()]
    # Hgrid = ogrid/omega0
    # # rgrid_macro = InputArchiveTDSE['rgrid_coarse'][:]
    # # zgrid_macro = InputArchiveTDSE['zgrid_coarse'][:]

print('data loaded:')



# plot the first spectrum



# fig, ax = plt.subplots()


image = pp.figure_driver()    
image.sf = [pp.plotter() for k1 in range(16)]

image.sf[0].args = [Hgrid, rgrid_FF, np.abs(HHG_onscreen[-1,:,:].T) ]
image.sf[0].method = plt.pcolormesh

pp.plot_preset(image)


k1 = mn.FindInterval(Hgrid, 17.0)
image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(16)]

image.sf[0].args = [units.INTENSITYau*I0_grid, rgrid_FF, np.abs(HHG_onscreen[:,k1,:].T) ]
image.sf[0].method = plt.pcolormesh

image.sf[0].colorbar.show = True
pp.plot_preset(image)


fig, ax = plt.subplots()
# cax = ax.pcolormesh(Hgrid, rgrid_FF, np.abs(HHG_onscreen[0,:,:].T)) #, vmin=-1, vmax=1, cmap='Blues')
cax = ax.pcolormesh(Hgrid, rgrid_FF, np.abs(HHG_onscreen[0,:,:].T), vmax=0.6e-9)
# ax.set_title(str(I0_grid[-1]))
ax.set_title("{:.2E}".format(1e-4*units.INTENSITYau*I0_grid[0]) + ' W/cm2')
fig.colorbar(cax)

def animate(i):
    cax.set_array(np.abs(HHG_onscreen[i,:,:].T)) #G[:-1, :-1, i].flatten())
    # ax.set_title(str(I0_grid[i]))
    # ax.set_title("{:10.4f}".format(units.INTENSITYau*I0_grid[i]))
    ax.set_title("{:.2E}".format(1e-4*units.INTENSITYau*I0_grid[i]) + ' W/cm2')
    # fig.colorbar(cf, cax=cax)
    # fig.colorbar(cax)
    # fig.colorbar(fig, cax=cax)
   
# anim = animation.FuncAnimation(fig, animate, interval=100, frames=len(I0_grid))
anim = animation.FuncAnimation(fig, animate, interval=250, frames=len(I0_grid))

plt.show()

# ## Laser
# NI0 = 300
# I0_start = 5e17/units.INTENSITYau
# I0_end = 35e17/units.INTENSITYau#E0_grid[-1]**2
# I0_grid = np.linspace(I0_start,E0_grid[-1]**2,NI0)
# w0 = 120e-6 #25e-6

# Gaussian_E_r = lambda r : np.exp(-(r/w0)**2)

# Nr = 200
# rgrid = np.linspace(0, 1.2*w0, Nr)


# Hlimit = [10, 36]
# # Hlimit = [24, 26]
# Hlimit = [15, 26]
# Hlimit = [16, 18]




# ## construct sources

# domega = ogrid[1]-ogrid[0]
# dE0 = E0_grid[1]-E0_grid[0]

# k_omega_sel = list(mn.FindInterval(Hgrid,Hlimit))

# ogrid_sel = ogrid[k_omega_sel[0]:k_omega_sel[1]]
# No_sel = len(ogrid_sel)

# FSourceTerm_sel = FSourceTerm[:,k_omega_sel[0]:k_omega_sel[1]]
# FSourceTerm_interpE0 = interpolate.interp1d( E0_grid, FSourceTerm_sel ,axis=0)


# # the sources
# # FSource_interp = FSourceTerm_interpE0( 0.75*np.sqrt(I0_grid[-1]) * Gaussian_E_r(rgrid) )
# FSource_interp = FSourceTerm_interpE0( np.sqrt(((2e18/units.INTENSITYau))) * Gaussian_E_r(rgrid) )
# # FSource_interp = FSourceTerm_interpE0( np.sqrt(I0_grid[-1]) * Gaussian_E_r(rgrid) )

# Nr_FF = 100
# rmax_FF = 0.012

# rgrid_FF = np.linspace(0,rmax_FF,Nr_FF)
# ## Hankel
# distance = 3.0 # 1.0
# omega_convert = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')


# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [rgrid ,ogrid_sel/omega0, np.log(np.abs(FSource_interp.T)) ]
# image.sf[0].method = plt.pcolormesh

# pp.plot_preset(image)

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [rgrid ,ogrid_sel/omega0, np.abs(FSource_interp.T) ]
# image.sf[0].method = plt.pcolormesh

# pp.plot_preset(image)

# HHG_onscreen = []
# for k1 in range(len(I0_grid)):
#     FSource_interp = FSourceTerm_interpE0( np.sqrt(I0_grid[k1]) * Gaussian_E_r(rgrid) )
#     HHG_onscreen.append(
#         Hfn2.HankelTransform(omega_convert * ogrid_sel, rgrid, FSource_interp.T, distance, rgrid_FF)
#         )

# HHG_onscreen = np.asarray(HHG_onscreen)
# ## reference plots

# # for k1 in range(len(I0_grid)):
# #     image = pp.figure_driver()    
# #     image.sf = [pp.plotter() for k1 in range(16)]

# #     image.sf[0].args = [ogrid_sel/omega0, rgrid_FF, np.log(np.abs(HHG_onscreen[k1,:,:].T)) ]
# #     image.sf[0].method = plt.pcolormesh


# #     # image.sf[1].args = [Hgrid[4], abs(FSourceTerm[7][15,:])]
# #     # image.sf[1].method = plt.semilogy


# #     pp.plot_preset(image)


# k1 = mn.FindInterval(ogrid_sel/omega0, 17.0)

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [units.INTENSITYau*I0_grid, rgrid_FF, np.log(np.abs(HHG_onscreen[:,k1,:].T)) ]
# image.sf[0].method = plt.pcolormesh

# pp.plot_preset(image)

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(16)]

# image.sf[0].args = [units.INTENSITYau*I0_grid, rgrid_FF, np.abs(HHG_onscreen[:,k1,:].T) ]
# image.sf[0].method = plt.pcolormesh

# pp.plot_preset(image)

# # image = pp.figure_driver()    
# # image.sf = [pp.plotter() for k1 in range(16)]

# # image.sf[0].args = [ogrid_sel/omega0, rgrid_FF, np.abs(HHG_onscreen.T) ]
# # image.sf[0].method = plt.pcolormesh


# # # image.sf[1].args = [Hgrid[4], abs(FSourceTerm[7][15,:])]
# # # image.sf[1].method = plt.semilogy


# # pp.plot_preset(image)





# # image = pp.figure_driver()    
# # image.sf = [pp.plotter() for k1 in range(16)]

# # image.sf[0].args = [ogrid_sel/omega0, np.real(FSourceTerm_sel[2500,:])]
# # image.sf[0].method = plt.semilogy

# # image.sf[1].args = [ogrid_sel/omega0, np.real(FSourceTerm_sel[2501,:])]
# # image.sf[1].method = plt.semilogy

# # image.sf[2].args = [ogrid_sel/omega0, np.real(FSourceTerm_interpE0(E0_grid[2500]+dE0/2.0))]
# # image.sf[2].method = plt.semilogy

# # # image.sf[1].args = [Hgrid[4], abs(FSourceTerm[7][15,:])]
# # # image.sf[1].method = plt.semilogy


# # pp.plot_preset(image)



# # image = pp.figure_driver()    
# # image.sf = [pp.plotter() for k1 in range(16)]

# # image.sf[0].args = [Hgrid, abs(FSourceTerm[-1,:])]
# # image.sf[0].method = plt.semilogy

# # image.sf[1].args = [Hgrid, abs(FSourceTerm[-2,:])]
# # image.sf[1].method = plt.semilogy

# # image.sf[2].args = [Hgrid, abs(FSourceTerm[-100,:])]
# # image.sf[2].method = plt.semilogy

# # pp.plot_preset(image)

# #######################################
# #######################################
# #######################################

# # def Efield_r(r):
# #     E0_r = Gaussian_E_r(r) # corresponding peak intensity




# # omega_au2SI = mn.ConvertPhoton(1.0, 'omegaau', 'omegaSI')
# # ogridSI = omega_au2SI * ogrid
# # Hgrid = ogrid/omega0
# # H_indices = [mn.FindInterval(Hgrid,Hvalue) for Hvalue in Hrange]

# # try:
# #     os.remove(out_h5name)
# #     print("previous results deleted")
# # except:
# #     print("no files deleted")  

# # rgrid_FF = np.linspace(0.0, rmax_FF, Nr_FF)
# # ogrid_select_SI = ogridSI[H_indices[0]:H_indices[1]:ko_step]



# # # Here are the fuction to obtain the phase factors in SI units: exp(i*omega*function(omega))
# # def f1_funct(E):
# #     return XUV_index.getf1(gas_type+'_' + XUV_table_type_diffraction, E)
# # def f2_funct(E):
# #     return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)


# # def dispersion_function_def(omega):
# #     f1_value = f1_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))    
# #     lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
# #     nXUV     = 1.0 - rho0_init*units.r_electron_classical * \
# #                ((lambdaSI**2)*f1_value/(2.0*np.pi))           
# #     phase_velocity_XUV  = units.c_light / nXUV
# #     return ((1./group_velocity_IR) - (1./phase_velocity_XUV))

# # def absorption_function_def(omega):
# #     f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
# #     lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
# #     beta_factor = rho0_init*units.r_electron_classical * \
# #                   ((lambdaSI**2)*f2_value/(2.0*np.pi))
# #     return beta_factor / units.c_light

# # if ('dispersion' in apply_diffraction): dispersion_function = dispersion_function_def
# # else: dispersion_function = None

# # if ('absorption' in apply_diffraction): absorption_function = absorption_function_def
# # else: dispersion_function = None

# # ## create subintervals to analyse the intensities
# # Hgrid_I_study = mn.get_odd_interior_points(Hrange)
# # omega_I_study_intervals = []
# # omega0SI = mn.ConvertPhoton(omega0, 'omegaau', 'omegaSI')
# # for k1 in Hgrid_I_study:
# #     omega_I_study_intervals.append(
# #         [np.max([ogrid_select_SI[0],omega0SI*(k1-0.5)]), np.min([ogrid_select_SI[-1],omega0SI*(k1+0.5)])]
# #         )

# # # The main integration
# # FField_FF_integrated, source_maxima = Hfn2.HankelTransform_long(
# #                                                ogrid_select_SI,
# #                                                rgrid_macro[0:Nr_max:kr_step],
# #                                                zgrid_macro[:Nz_max_sum],
# #                                                FSourceTerm[0:Nr_max:kr_step,:Nz_max_sum,H_indices[0]:H_indices[1]:ko_step],
# #                                                distance_FF,
# #                                                rgrid_FF,
# #                                                dispersion_function = dispersion_function, # None, #dispersion_function,
# #                                                absorption_function = absorption_function,
# #                                                frequencies_to_trace_maxima = omega_I_study_intervals)


# # # Save the data
# # Hgrid_select = Hgrid[H_indices[0]:H_indices[1]:ko_step]
# # with h5py.File(out_h5name,'w') as OutFile:
# #     grp = OutFile.create_group('XUV')
# #     grp.create_dataset('Spectrum_on_screen',
# #                                           data = np.stack((FField_FF_integrated.real, FField_FF_integrated.imag),axis=-1)
# #                                           )
# #     grp.create_dataset('Maxima_of_planes',
# #                                           data = np.asarray(source_maxima)
# #                                           )
# #     grp.create_dataset('Maxima_Hgrid',
# #                                           data = Hgrid_I_study
# #                                           )
# #     grp.create_dataset('rgrid_FF',
# #                                           data = rgrid_FF
# #                                           )    
# #     grp.create_dataset('Hgrid_select',
# #                                           data = Hgrid_select
# #                                           )
        


# # # vmin = np.max(np.log(Gaborr))-6.
# # fig, ax = plt.subplots()   
# # FF_spectrum_logscale = np.log10(abs(FField_FF_integrated.T)**2);
# # vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# # map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# # fig.colorbar(map1)
# # plt.title('Far-field spectrum (30 cm), integrated, log')
# # plt.xlabel('H [-]')
# # plt.ylabel('r [m]')
# # if showplots: plt.show()
# # # plt.close(fig)
# # # sys.exit()
