import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py

import matplotlib.pyplot as plt



   




with h5py.File(os.path.join("D:\data", "JZ","density_mod","series10","test1_modT2","results.h5"),'r') as myfile:   # hdf5: ztr order
    
    
    field = myfile['/outputs/output_field'][-1,:,0]
    
    plt.figure()
    plt.plot(field)
    plt.show()
    
