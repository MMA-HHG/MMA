import sys
import h5py
import MMA_administration as MMA

with open('msg.tmp','r') as msg_file:
    results_file = msg_file.readline()[:-1] # need to strip the last character due to Fortran msg.tmp

 
with h5py.File('results_Hankel.h5', 'r') as f_src, h5py.File(results_file, 'a') as f_dest:
    f_src.copy(f_src[MMA.paths['Hankel_outputs']],f_dest[MMA.paths['Hankel']])

print("Done")
