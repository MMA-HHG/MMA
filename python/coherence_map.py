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

filename = "results_1.h5"

file_path = os.path.join(results_path,filename)

with h5py.File(file_path, 'r') as InputArchive:

    print(InputArchive['/inputs/laser_wavelength'][()])

    electron_density_map = 5
    plasma_frequency_map = np.sqrt((units.elcharge**2)/(units.eps0*units.elmass)) * np.sqrt(electron_density_map)


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

