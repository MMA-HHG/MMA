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
from contextlib import ExitStack

# import mynumerics as mn
import matplotlib.pyplot as plt

  
# filename = 'Hankel.h5'
# filename = 'Hankel_dx2.h5'
# filename = 'Hankel_dt2.h5'
filename = 'Hankel_2Nx.h5'

cwd = os.getcwd()
results_path = cwd

FF_orders_plot = 4

files = ['Hankel.h5','Hankel_dx2.h5','Hankel_dt2.h5','Hankel_2Nx.h5']
Nfiles = len(files)
Maxima = []; FField_FF = []; Hgrid = []; rgrid_FF = []

with ExitStack() as stack:    
    InArch = [stack.enter_context(h5py.File(os.path.join(results_path,fname), 'r')) for fname in files]
    for k1 in range(Nfiles):
        # load data 
        Maxima.append(InArch[k1]['XUV/Maxima_of_planes'][:])
        FField_FF.append(InArch[k1]['XUV/Spectrum_on_screen'][:,:,0] + \
                         1j*InArch[k1]['XUV/Spectrum_on_screen'][:,:,1])
        Hgrid.append(InArch[k1]['XUV/Hgrid_select'][:])
        rgrid_FF.append(InArch[k1]['XUV/rgrid_FF'][:])

sys.exit()
            
# with h5py.File(filename, 'r') as InputArchive:
#     # load data
#    Maxima = InputArchive['XUV/Maxima_of_planes'][:]
#    FField_FF = InputArchive['XUV/Spectrum_on_screen'][:,:,0] + \
#                1j*InputArchive['XUV/Spectrum_on_screen'][:,:,1]
#    Hgrid = InputArchive['XUV/Hgrid_select'][:]
#    rgrid_FF = InputArchive['XUV/rgrid_FF'][:]
   
   
   
fig, ax = plt.subplots()     
plt.plot(Maxima[0,:])
plt.plot(Maxima[1,:])
plt.plot(Maxima[2,:])
plt.plot(Maxima[3,:])
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

   
   
   
   
   
   
   
