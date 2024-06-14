import h5py
import numpy as np
import MMA_administration as MMA

with open('msg.tmp','r') as msg_file:
    results_file = msg_file.readline()
    
with h5py.File(results_file, 'a') as h5f:
    Nz = len(h5f[MMA.paths['CUPRAD_outputs']+'/zgrid'][:])
    
    dset_id = h5f.create_dataset(MMA.paths['CTDSE_inputs']+'/Nz_max', data=Nz)
    dset_id.attrs['units'] = np.string_('[-]')

print("Added the maximal Nz length ("+str(Nz)+") for CTDSE.")
