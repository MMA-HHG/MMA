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
   
   
   
   
   
   
   
