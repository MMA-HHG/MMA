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

import warnings

# import mynumerics as mn
import matplotlib.pyplot as plt

  
# filename = 'Hankel.h5'
# filename = 'Hankel_dx2.h5'
# filename = 'Hankel_dt2.h5'
# filename = 'Hankel_2Nx.h5'
# filename = '60pl/Hankel.h5'
# filename = '60pl/Hankel_dr2.h5'
# filename = 'PoIs/Hankel_1250pl.h5'
# filename = 'PoIs/Hankel_500pl.h5'
# filename = 'PoIs/Hankel_1000pl.h5'
# filename = 'PoIs/Hankel_all_cummulative1.h5'
# filename = 'PoIs/Hankel_all_cummulative2.h5'

# filename = 'PoIs/Hankel_all_cummulative1.h5'


results_path = os.path.join("D:\data", "Discharges", "TDSE","scan1")

filename = 'Hankel_all_cummulative_20_8.h5'

FF_orders_plot = 4

filename_path = os.path.join(results_path,filename)
         
with h5py.File(filename_path, 'r') as InputArchive:
    # load data
   data_group = InputArchive['XUV']
   available_data = list(data_group.keys())
    
   Maxima = InputArchive['XUV/Maxima_of_planes'][:] #/np.pi
   Phases_onax = InputArchive['XUV/Phase_on_axis'][:] #/np.pi
   Phases_first = InputArchive['XUV/Phase_first_plane'][:]
   FField_FF = InputArchive['XUV/Spectrum_on_screen'][:,:,0] + \
               1j*InputArchive['XUV/Spectrum_on_screen'][:,:,1]
   Hgrid = InputArchive['XUV/Hgrid_select'][:]; NH = len(Hgrid)
   rgrid_FF = InputArchive['XUV/rgrid_FF'][:]
   zgrid_integration = InputArchive['XUV/zgrid_integration'][:]; Nz = len(zgrid_integration)
   
   Hgrid_study = InputArchive['XUV/Maxima_Hgrid'][:]
   
   cummulative_spectrum = ('Spectrum_on_screen_cummulative' in available_data)
   if cummulative_spectrum:
       FField_FF_cummulative = InputArchive['XUV/Spectrum_on_screen_cummulative'][:,:,:,0] + \
               1j*InputArchive['XUV/Spectrum_on_screen_cummulative'][:,:,:,1]
               
   L_abs_analyse = ('L_abs_Hgrid' in available_data)
   if L_abs_analyse:
       L_abs_Hgrid = InputArchive['XUV/L_abs_Hgrid'][:]
               
 
zgrid_integration_midpoints = 0.5*(zgrid_integration[1:]+zgrid_integration[:-1])
   
Maxima = Maxima/np.max(Maxima)  
   
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


# Compute integrated spectrum
dE_dH = np.zeros((NH,))
for k1 in range(NH):
    dE_dH[k1] = np.trapz(np.abs(FField_FF[k1,:])**2)


fig, ax = plt.subplots()     
plt.plot(Hgrid,dE_dH)
plt.xlabel('H [-]')
plt.ylabel('dE/dH [arb. u.]')
plt.show()



print(np.max(abs(FField_FF.T)**2))

Hs_to_trace_maxima = []
for k1 in range(len(Hgrid_study)):
    Hs_to_trace_maxima.append([Hgrid_study[k1]-0.5 , Hgrid_study[k1]+0.5])

H_indices = []
planes_maxima = []
XUV_beams_energy = []
for H_list in Hs_to_trace_maxima:
    try:
        H_indices.append(mn.FindInterval(Hgrid, H_list))
        planes_maxima.append([]); XUV_beams_energy.append([])
    except:
        warnings.warn("A frequency from frequencies_to_trace_maxima doesn't match ogrid.")

dE_dH_z = np.zeros((NH,))
    
if (len(H_indices)>0):
    for k1 in range(Nz-1):
        for k2 in range(NH): dE_dH_z[k2] = np.trapz(np.abs(FField_FF_cummulative[k1,k2,:])**2) 
        for k2 in range(len(H_indices)):
            planes_maxima[k2].append(np.max(abs(
                              FField_FF_cummulative[k1,H_indices[k2][0]:H_indices[k2][1],:]
                                    )))
            XUV_beams_energy[k2].append(
                              mn.integrate_subinterval(dE_dH_z,Hgrid,[Hgrid_study[k2]-0.5 , Hgrid_study[k2]+0.5])
                                    )
            
    for k1 in range(len(H_indices)):
        planes_maxima[k1] = np.asarray(planes_maxima[k1])
        XUV_beams_energy[k1] = np.asarray(XUV_beams_energy[k1])

planes_maxima = np.asarray(planes_maxima)
XUV_beams_energy = np.asarray(XUV_beams_energy)

# find maxima in the case we need them
fig, ax = plt.subplots()     
plt.plot(planes_maxima[0,:])
plt.plot(planes_maxima[1,:])
plt.plot(planes_maxima[2,:])
plt.plot(planes_maxima[3,:])
plt.show()

fig, ax = plt.subplots()     
plt.plot(zgrid_integration_midpoints, planes_maxima[0,:])
plt.plot(zgrid_integration_midpoints, planes_maxima[1,:])
plt.plot(zgrid_integration_midpoints, planes_maxima[2,:])
plt.plot(zgrid_integration_midpoints, planes_maxima[3,:])
plt.show()

fig, ax = plt.subplots()     
plt.plot(zgrid_integration_midpoints, XUV_beams_energy[0,:])
plt.plot(zgrid_integration_midpoints, XUV_beams_energy[1,:])
plt.plot(zgrid_integration_midpoints, XUV_beams_energy[2,:])
plt.plot(zgrid_integration_midpoints, XUV_beams_energy[3,:])
plt.show()
   
   
   
if L_abs_analyse:  
    fig, ax = plt.subplots()     
    plt.plot(Hgrid,1e3*L_abs_Hgrid)
    plt.xlabel('H [-]')
    plt.ylabel('Labs [mm]')
    plt.show()
   
   
