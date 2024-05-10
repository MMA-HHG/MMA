import numpy as np
import h5py
# import sys
import mynumerics as mn

import MMA_administration as MMA

# import plot_presets as pp

# import HHG
# import units


# modulation function
def Gaussian_jet(zgrid_, z0_, a_, peak_max = 1.):
    return peak_max * np.exp(-((zgrid_ - z0_)/a_)**2)

# Gaussian train
def Gaussian_train(zgrid_, p_max_, a_, positions_):
    res = p_max_ * np.exp(-((zgrid_-positions_[0])/a_)**2)
    for position in positions_[1:]:  res += p_max_ * np.exp(-((zgrid_-position)/a_)**2)
    return res

def Gauss_train_stride(zgrid_, p_max_, a_, pos0_, stride_, N_jets_): # the defining parameters of the train are now the position of the first jet and the 'stride'
    return Gaussian_train(zgrid_, p_max_, a_, [pos0_ + k1*stride_ for k1 in range(N_jets_)])


# prepare zgrid
zmin = 0.
zmax = 2e-3 # float(sys.argv[1]) # 1e-2
dz = 100e-8
zgrid = np.ogrid[zmin:zmax+dz:dz]

a_Gauss = 500e-6



# image = pp.figure_driver()
    
# image.sf = [pp.plotter() for k1 in range(32)]

# image.sf[0].args = [zgrid, Gaussian_jet(zgrid, zmax/2., a_Gauss)]  


# pp.plot_preset(image)  

# print(np.trapz(Gaussian_jet(zgrid, zmax/2., a_Gauss),zgrid)/zmax)

# add to file
with h5py.File('results.h5', 'a') as f:
    mn.adddataset(f, MMA.paths['global_inputs']+'/density_mod/zgrid', zgrid, '[m]')
    mn.adddataset(f, MMA.paths['global_inputs']+'/density_mod/table', Gaussian_jet(zgrid, zmax/2., a_Gauss), '[-]')
    
    # f['density_mod/zgrid'] = zgrid
    # f['density_mod/table'] = Gaussian_jet(zgrid, 0.5*zmax, 5e-3)
    
    
    # outf.require_group(h5path)
    


print('done')
