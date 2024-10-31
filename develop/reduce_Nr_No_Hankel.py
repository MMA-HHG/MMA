import sys
import h5py
import numpy as np
import MMA_administration as MMA

with open('msg.tmp','r') as msg_file:
    results_file = msg_file.readline()[:-1] # need to strip the last character due to Fortran msg.tmp

arguments = sys.argv    
with h5py.File(results_file, 'a') as h5f:
    h5f[MMA.paths['Hankel_inputs']+'/Nr_FF'][()] = arguments[1]    # 20
    h5f[MMA.paths['Hankel_inputs']+'/rmax_FF'][()] = arguments[2]  # 0.005
    h5f[MMA.paths['Hankel_inputs']+'/ko_step'][()] = arguments[3]  # 3

print("Done")
