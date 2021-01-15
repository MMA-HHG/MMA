import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
sys.path.append('D:\git\python_modules')
import units
import mynumerics as mn
# import mynumerics as mn
import matplotlib.pyplot as plt

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")

filename = "results_7.h5"

file_path = os.path.join(results_path,filename)

with h5py.File(file_path, 'r') as InputArchive:

    print(1e-2*InputArchive['/inputs/laser_wavelength'][()])
    omega0 = mn.ConvertPhoton(1e-2*InputArchive['/inputs/laser_wavelength'][()],'lambdaSI','omegaSI')
    print(omega0)
    tgrid = InputArchive['/outputs/tgrid'][:]
    rgrid = InputArchive['/outputs/rgrid'][:]
    zgrid = InputArchive['/outputs/zgrid'][:]
    electron_density_map = InputArchive['/outputs/output_plasma'][:]
    Efield = InputArchive['/outputs/output_field'][:]

    plasma_frequency_map = np.sqrt((units.elcharge ** 2) / (units.eps0 * units.elmass)) * np.sqrt(electron_density_map)


    # n_plasma_map = np.sqrt(1-(plasma_frequency_map/omega0)**2)
    #
    # k_plas_map = 25*(2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))*(n_plasma_map-1)
    #
    # plt.pcolor(tgrid,rgrid,k_plas_map[:,:,0].T)
    # plt.colorbar()
    # plt.savefig('kmap.png', dpi = 600)
    # plt.show()
    #
    # plt.pcolor(tgrid,rgrid,abs(np.pi/k_plas_map[:,:,15].T))
    # plt.colorbar()
    # plt.savefig('Lcohmap.png', dpi = 600)
    # plt.show()
    #
    # xxx = electron_density_map[0:3,0,0]
    # yyy = electron_density_map[0, 0:3, 0]

    plt.plot(tgrid, Efield[:,0,1],linewidth=0.2)
    plt.plot(tgrid, Efield[:, 0, 2],linewidth=0.2)
    plt.savefig('Efield.png', dpi = 600)
    plt.show()


    plt.pcolor(tgrid,rgrid,plasma_frequency_map[:,:,0].T)
    plt.colorbar()
    plt.savefig('map_zfix.png', dpi = 600)
    plt.show()


    plt.pcolor(zgrid,tgrid,np.squeeze(plasma_frequency_map[:,0,:]))
    plt.colorbar()
    plt.savefig('map_onax.png', dpi = 600)
    plt.show()

    plt.pcolor(zgrid,rgrid,np.squeeze(plasma_frequency_map[512,:,:]))
    plt.colorbar()
    plt.savefig('map_thalf.png', dpi = 600)
    plt.show()

print('test1')


## refractive index contribs


#### THE MAIN PROGRAM #####


# class NumericalParamsClass:
#   def __init__(self):
#     pass
#
# NumericalParams.z_medium = ParamFile['inputs/'+'jetpositions'][()]

# test matplot

plt.plot([1, 2, 3, 4])
plt.ylabel('some numbers')
plt.show()

