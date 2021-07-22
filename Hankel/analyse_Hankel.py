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

  

FF_orders_plot = 4
         
with h5py.File('Hankel.h5', 'r') as InputArchive:
    # load data
   Maxima = InputArchive['XUV/Maxima_of_planes'][:]
   FField_FF = InputArchive['XUV/Spectrum_on_screen'][:,:,0] + \
               1j*InputArchive['XUV/Spectrum_on_screen'][:,:,1]
   
   
   
   
plt.plot(Maxima[0,:])
plt.show()

# vmin = np.max(np.log(Gaborr))-6.
fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF.T)**2);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
map1 = ax.pcolor(FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field spectrum (30 cm), integrated, log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
plt.show()
# if showplots: plt.show()
# plt.close(fig)
# sys.exit()

   
   
   
   
   
   
   
