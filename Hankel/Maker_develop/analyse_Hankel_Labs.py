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

import XUV_refractive_index as XUV_index

  
# filename = 'Hankel.h5'
# filename = 'Hankel_dx2.h5'
# filename = 'Hankel_dt2.h5'
# filename = 'Hankel_2Nx.h5'
# filename = '60pl/Hankel.h5'
# filename = '60pl/Hankel_dr2.h5'
filename = 'PoIs/Hankel_1250pl.h5'
# filename = 'PoIs/Hankel_500pl.h5'
# filename = 'PoIs/Hankel_1000pl.h5'

results_CUPRAD = os.path.join("D:\data", "Discharges", "TDSE", "t6")
file_CUPRAD = 'results_1.h5'
file_CUPRAD = os.path.join(results_CUPRAD,file_CUPRAD)

FF_orders_plot = 4
         
with h5py.File(filename, 'r') as InputArchive, h5py.File(file_CUPRAD, 'r') as InputArchiveCUPRAD:
    # load data
   Maxima = InputArchive['XUV/Maxima_of_planes'][:] #/np.pi
   Maxima_Hgrid = InputArchive['XUV/Maxima_Hgrid'][:]
   Phases_onax = InputArchive['XUV/Phase_on_axis'][:] #/np.pi
   Phases_first = InputArchive['XUV/Phase_first_plane'][:]
   FField_FF = InputArchive['XUV/Spectrum_on_screen'][:,:,0] + \
               1j*InputArchive['XUV/Spectrum_on_screen'][:,:,1]
   Hgrid = InputArchive['XUV/Hgrid_select'][:]
   rgrid_FF = InputArchive['XUV/rgrid_FF'][:]
   zgrid_integration = InputArchive['XUV/zgrid_integration'][:]
   
   omega0SI = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchiveCUPRAD,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
   # rho0_init = 1e6 * mn.readscalardataset(InputArchiveCUPRAD, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N') # SI
   rho0_init = 1e6 * 8.1e17 # SI

gas_type = 'Kr'
XUV_table_type_absorption = 'Henke' # {Henke, NIST}    
def f2_funct(E):
    return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)

def L_abs(omega):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 1.0 / (2.0 * rho0_init * units.r_electron_classical * lambdaSI * f2_value) 
    
Maxima = Maxima/np.max(Maxima)  

L_abs_Hgrid = np.zeros((len(Hgrid),))
for k1 in range(len(L_abs_Hgrid)):
    L_abs_Hgrid[k1] = L_abs(omega0SI*Hgrid[k1])

L_abs_Maxima_Hgrid = np.zeros((len(Maxima_Hgrid),))
for k1 in range(len(L_abs_Maxima_Hgrid)):
    L_abs_Maxima_Hgrid[k1] = L_abs(omega0SI*Maxima_Hgrid[k1])


fig, ax = plt.subplots()     
plt.plot(zgrid_integration,Maxima[0,:])
plt.plot(zgrid_integration,
         np.exp((zgrid_integration-zgrid_integration[-1])/(2.0*L_abs_Maxima_Hgrid[0]))
         )
plt.show()


fig, ax = plt.subplots()     
plt.semilogy(zgrid_integration,Maxima[0,:])
plt.semilogy(zgrid_integration,
         np.exp((zgrid_integration-zgrid_integration[-1])/(2.0*L_abs_Maxima_Hgrid[0]))
         )
plt.show()
    

fig, ax = plt.subplots()     
plt.plot(Maxima[0,:])
plt.plot(Maxima[1,:])
plt.plot(Maxima[2,:])
plt.plot(Maxima[3,:])
plt.show()

fig, ax = plt.subplots()     
plt.semilogy(Maxima[0,:])
plt.semilogy(Maxima[1,:])
plt.semilogy(Maxima[2,:])
plt.semilogy(Maxima[3,:])
plt.show()

fig, ax = plt.subplots()     
plt.plot(Phases_onax[0,:])
plt.plot(Phases_onax[1,:])
plt.plot(Phases_onax[2,:])
plt.plot(Phases_onax[3,:])
plt.show()

fig, ax = plt.subplots()     
plt.plot(Phases_first[0,:])
plt.plot(Phases_first[1,:])
plt.plot(Phases_first[2,:])
plt.plot(Phases_first[3,:])
plt.show()

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
plt.show()
# if showplots: plt.show()
# plt.close(fig)
# sys.exit()

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(FField_FF.T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
map1 = ax.pcolor(Hgrid,rgrid_FF,abs(FField_FF.T)**2, shading='auto')
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
plt.show()
# if showplots: plt.show()
# plt.close(fig)
# sys.exit()

   
fig, ax = plt.subplots()     
plt.plot(Hgrid,1e3*L_abs_Hgrid)
plt.show()   

# LabsH17 = 1e3*L_abs(omega0SI*17)
   
   
   
   
   
