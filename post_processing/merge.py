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
outfname = "merged.h5"
available_outputs_list = ['Efield', 'Gabor']

print(files)

def prepare_ouput_file(f,outf,dset_list):
    joint_rz_shape = (mn.readscalardataset(f,'Nr_orig','N')[0], mn.readscalardataset(f,'Nz_orig','N')[0])
    for dsetname in dset_list:
        if (dsetname in available_outputs_list):
            dset = f[dsetname]
            newshape = dset.shape[0:-1] + joint_rz_shape
            newdset = outf.create_dataset(dsetname, newshape,'d')







with h5py.File(outfname,'w') as outf:
    firstrun = True;
    for fname in files:
        print("file")
        with h5py.File(fname,'r') as f:
            dset_list = list(f.keys())
            if firstrun:
                prepare_ouput_file(f, outf, dset_list)
                firstrun = False;
            else:
                pass
            nsim_loc = mn.readscalardataset(f,'number_of_local_simulations','N')[0]
            print(nsim_loc)
            data = f["Efield"][()]
            print(data[:,0:nsim_loc])
            # if first: prepare
            # else: print
            print(list(f.keys()))
            myshape = f["FEfield"].shape
            print(myshape[0:-1])


reflist = ["a", "car", "dog"]
print("car" in reflist)
print("a" in reflist)
print("b" in reflist)


# scan datasets and then fill it according to keys and lengths (extend each dataset by 2 (resolve dimensions...))