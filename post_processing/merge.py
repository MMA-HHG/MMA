# https://stackoverflow.com/questions/17198319/how-to-configure-custom-pythonpath-with-vm-and-pycharm

# from scipy import special
# from scipy import integrate
# from scipy import interpolate
import numpy as np
# import struct
# import array
import os
# import time
# # import ray
# # import matlab.engine
# # import string
# import multiprocessing as mp
import math
# # import joblib
# # from mpi4py import MPI
# # import oct2py
import shutil
import h5py
import sys
import units
import mynumerics as mn
import glob

print(mn.IsPowerOf2(4))

files = glob.glob('hdf5_temp_*.h5')
print(files)



k1 = 0;
for fname in files:
    print("file")
    with h5py.File(fname,'r') as f:
        # if first: prepare
        # else: print

        print(list(f.keys()))


reflist = ["a", "car", "dog"]
print("car" in reflist)
print("a" in reflist)
print("b" in reflist)


# scan datasets and then fill it according to keys and lengths (extend each dataset by 2 (resolve dimensions...))