import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
sys.path.append('D:\git\python_modules')
import units
# import mynumerics as mn
import matplotlib.pyplot as plt

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")

filename = "results_7.h5"

file_path = os.path.join(results_path,filename)

with h5py.File(file_path, 'r') as InputArchive:

    print(InputArchive['/inputs/laser_wavelength'][()])
    tgrid = InputArchive['/outputs/tgrid'][:]
    rgrid = InputArchive['/outputs/rgrid'][:]
    zgrid = InputArchive['/outputs/zgrid'][:]
    electron_density_map = InputArchive['/outputs/output_plasma'][:]
    xxx = electron_density_map[0:3,0,0]
    yyy = electron_density_map[0, 0:3, 0]

    plasma_frequency_map = np.sqrt((units.elcharge**2)/(units.eps0*units.elmass)) * np.sqrt(electron_density_map)

    plt.pcolor(tgrid,rgrid,plasma_frequency_map[:,:,0].T)
    plt.colorbar()
    plt.savefig('map.png', dpi = 600)
    plt.show()




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

