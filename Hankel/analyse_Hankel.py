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

           
with h5py.File('Hankel.h5', 'r') as InputArchiveCUPRAD:
    # load data
   Maxima = InputArchiveCUPRAD['XUV/Maxima_of_planes'][:]
   
   
   
   
plt.plot(Maxima[0,:])
plt.show()
   
   
   
   
   
   
   
   
   # PopTot_TDSE = InputArchiveTDSE['PopTot'][0,0,:]
   # PopInt_TDSE = InputArchiveTDSE['PopInt'][0,0,:]
   # expval_x_TDSE = InputArchiveTDSE['expval_x'][0,0,:]
   
   # GS_init = InputArchiveTDSE['ground_state'][:,0] + 1j*InputArchiveTDSE['ground_state'][:,1]

# print('data loaded:')

# Delta_t_TDSE = tgrid[-1] - tgrid[0]
# Delta_t_CUPRAD = tgrid_CUPRAD[-1] - tgrid_CUPRAD[0]

# a1 = Delta_t_CUPRAD/Delta_t_TDSE
# b1 = tgrid_CUPRAD[-1] - a1*tgrid[-1]

# tgrid_TDSE_resc = a1*tgrid + b1

# index_t0 = mn.FindInterval(tgrid_TDSE_resc, 0.0e-15)




# plasma_map = Plasma[:,:,index_t0]



# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# map1 = ax.pcolor(1e3*zgrid_macro, 1e6*rgrid_macro, 100*(1.0-plasma_map), shading='auto', cmap='plasma')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('plasma, t=0 fs')
# plt.xlabel('z [mm]')
# plt.ylabel('r [mum]')
# plt.ylim([0, 130])
# plt.show()
# # plt.close(fig)


# plasma_map = Plasma[:,:,-1]
# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# map1 = ax.pcolor(1e3*zgrid_macro, 1e6*rgrid_macro, 100*(1.0-plasma_map), shading='auto', cmap='plasma')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('plasma, t_end')
# plt.xlabel('z [mm]')
# plt.ylabel('r [mum]')
# plt.ylim([0, 130])
# plt.show()
# # plt.close(fig)