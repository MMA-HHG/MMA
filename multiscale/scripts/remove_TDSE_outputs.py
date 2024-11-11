import h5py
import MMA_administration as MMA

with open('msg.tmp','r') as msg_file:
    results_file = msg_file.readline()[:-1] # need to strip the last character due to Fortran msg.tmp

 
with h5py.File(results_file, 'a') as f:
    del f[MMA.paths['CTDSE_outputs']]

print("Done")
