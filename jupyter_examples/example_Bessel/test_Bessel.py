import numpy as np
import h5py
import mynumerics as mn
import os

import dataformat_CUPRAD as dfC

import matplotlib.pyplot as plt
import plot_presets as pp

h5file1 = os.path.join('D:\sharepoint', 'OneDrive - ELI Beamlines', 'data', 'Sunrise', 'demos', 'Bessel', 'results', 'results_Bessel_4.h5')

h5filename = 'results_Bessel.h5'

with h5py.File(h5file1,'r') as f:
    Efield_r = f['CUPRAD/pre-processed/startfield_r'][:]
    Efield_i = f['CUPRAD/pre-processed/startfield_i'][:]
    CUPRAD_res = dfC.get_data(f)
    
    
# fig, ax = plt.subplots()
# ax.plot(Efield_r[0,:])
# ax.plot(Efield_i[0,:], ':')
# plt.show()  


# Test complexification
Nr, Nt = np.shape(Efield_r)
Efield_cmplx_pyth = np.empty((Nr, Nt),dtype=np.cdouble)

for k1 in range(Nr): 
    Efield_cmplx_pyth[k1,:] = mn.complexify_fft(Efield_r[k1,:],convention='-')
    
    
# fig, ax = plt.subplots()
# ax.plot(Efield_cmplx_pyth[0,:].real)
# ax.plot(Efield_cmplx_pyth[0,:].imag, ':')
# plt.show()  

# fig, ax = plt.subplots()
# ax.plot(Efield_i[0,:])
# ax.plot(Efield_cmplx_pyth[0,:].imag, ':')
# plt.show()  

fig, ax = plt.subplots()
ax.pcolormesh(1e15*CUPRAD_res.tgrid, 1e6*CUPRAD_res.rgrid, CUPRAD_res.E_zrt[CUPRAD_res.Nz//2,:,:], shading='auto')
plt.show()  

fig, ax = plt.subplots()
ax.pcolormesh(1e15*CUPRAD_res.tgrid, 1e6*CUPRAD_res.rgrid, CUPRAD_res.E_zrt[-1,:,:]**2, shading='auto')
plt.show()  


fig, ax = plt.subplots()
ax.plot(1e15*CUPRAD_res.tgrid,Efield_r[0,:])
ax.plot(1e15*CUPRAD_res.tgrid,CUPRAD_res.E_zrt[0,0,:], ':')
ax.plot(1e15*CUPRAD_res.tgrid,CUPRAD_res.E_zrt[CUPRAD_res.Nz//2,0,:],)
plt.show()  


xxx = (Efield_i-Efield_cmplx_pyth.imag)



# fig, ax = plt.subplots()
# ax.plot(Efield_i[0,:])
# ax.plot(xxx[0,:])
# plt.show()  

# np.testing.assert_allclose(Efield_i,Efield_cmplx_pyth.imag)

