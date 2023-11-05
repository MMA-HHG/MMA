import numpy as np
import h5py
import mynumerics as mn


# modulation function
def Gaussian_jet(zgrid_, z0_, a_, peak_max = 1.):
    return peak_max * np.exp(-((zgrid_ - z0_)/a_)**2)


# prepare zgrid
zmin = 0.
zmax = 1e-2
dz = 100e-6
zgrid = np.ogrid[zmin:zmax+dz:dz]


# add to file
with h5py.File('results.h5', 'a') as f:
    mn.adddataset(f, 'density_mod/zgrid', zgrid, '[m]')
    mn.adddataset(f, 'density_mod/table', Gaussian_jet(zgrid, 0.5*zmax, 5e-3), '[-]')
    # f['density_mod/zgrid'] = zgrid
    # f['density_mod/table'] = Gaussian_jet(zgrid, 0.5*zmax, 5e-3)
    


print('done')
