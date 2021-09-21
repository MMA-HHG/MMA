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

filename = 'PoIs/Hankel_all_cummulative_noabs2.h5'

FF_orders_plot = 4
         
with h5py.File(filename, 'r') as InputArchive:
    # load data
   data_group = InputArchive['XUV']
   available_data = list(data_group.keys())
    
   Maxima = InputArchive['XUV/Maxima_of_planes'][:] #/np.pi
   Phases_onax = InputArchive['XUV/Phase_on_axis'][:] #/np.pi
   Phases_first = InputArchive['XUV/Phase_first_plane'][:]
   FField_FF = InputArchive['XUV/Spectrum_on_screen'][:,:,0] + \
               1j*InputArchive['XUV/Spectrum_on_screen'][:,:,1]
   Hgrid = InputArchive['XUV/Hgrid_select'][:]
   rgrid_FF = InputArchive['XUV/rgrid_FF'][:]
   zgrid_integration = InputArchive['XUV/zgrid_integration'][:]; Nz = len(zgrid_integration)
   
   Hgrid_study = InputArchive['XUV/Maxima_Hgrid'][:]
   
   cummulative_spectrum = ('Spectrum_on_screen_cummulative' in available_data)
   if cummulative_spectrum:
       FField_FF_cummulative = InputArchive['XUV/Spectrum_on_screen_cummulative'][:,:,:,0] + \
               1j*InputArchive['XUV/Spectrum_on_screen_cummulative'][:,:,:,1]
   
 
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

if (len(H_indices)>0):
    for k1 in range(Nz-1):
        for k2 in range(len(H_indices)):
            planes_maxima[k2].append(np.max(abs(
                              FField_FF_cummulative[k1,H_indices[k2][0]:H_indices[k2][1],:]
                                    )))
            
    for k1 in range(len(H_indices)):
        planes_maxima[k1] = np.asarray(planes_maxima[k1])

planes_maxima = np.asarray(planes_maxima)

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
   
   
   
   
   
   
