####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2020)
#
# The purpose of this code is to merge outputs generated from 1DTDSE. It uses a general
# procedure for assigning indices based on keys, it thus can be easily generalised for more
# dimensions.



# finding mn from external location
# https://stackoverflow.com/questions/17198319/how-to-configure-custom-pythonpath-with-vm-and-pycharm

import numpy as np
import os
import math
import shutil
import h5py
import sys
import units
import mynumerics as mn
import glob


files = glob.glob('hdf5_temp_*.h5') # filter all the single-proc files
outfname = "results_merged.h5"
available_outputs_list = ['Efield', 'FEfield', 'SourceTerm', 'FSourceTerm', 'FEfieldM2', 'FSourceTermM2', 'PopTot', 'Gabor',
                          'PopInt', 'expval_x'] # Gabor is not implemented, it's here to test an extra argument
available_further_data = ['tgrid', 'omegagrid', 'Energy_of_the_ground_state', 'xgrid_micro', 'ground_state', 'zgrid_coarse',
                          'rgrid_coarse'] # these are in all the files and supposed to be same, e.g. grids
precision = 'd'

def prepare_ouput_file(f,outf,dset_list,nsim_tot):
    for dsetname in dset_list:

        if (dsetname in available_outputs_list):
            dset = f[dsetname]
            newshape =  (nsim_tot,) + dset.shape[0:-1]
            newdset = outf.create_dataset(dsetname, newshape,precision)

        if (dsetname in available_further_data):
            f.copy(dsetname,outf)


def print_ouput_file(f,outf,dset_list):
    Nr = mn.readscalardataset(f, 'Nr_orig', 'N')[0]
    Nz = mn.readscalardataset(f, 'Nz_orig', 'N')[0]
    nsim_loc = mn.readscalardataset(f, 'number_of_local_simulations', 'N')[0]
    ks_loc = f['keys'][()]
    for dsetname in dset_list:
        if (dsetname in available_outputs_list):
            dset = f[dsetname]
            ndims = len(dset.shape)
            newdset = outf[dsetname]
            data = dset[()]
            for k1 in range(nsim_loc):
                if (ks_loc[k1] >= 0): # only for used keys
                    # kr, kz = mn.n1n2mapping_inv(ks_loc[k1], Nr)
                    # print(kr,kz)
                    if (ndims == 2): # for reals
                        newdset[ks_loc[k1],:] = dset[:,k1]
                    elif (ndims == 3):  # for complex
                        newdset[ks_loc[k1],:,:] = dset[:,:,k1]
                    else:
                        print('warning, dataset with unsupported dimension: ' + dsetname + ', nothing done')# general procedure should be possible
    return nsim_loc

# total number of simulations

nsim_tot = 0
for fname in files:
    with h5py.File(fname,'r') as f:
        nsim_tot += mn.readscalardataset(f, 'number_of_local_simulations', 'N')[0]
            
# nsim_tot = 0
with h5py.File(outfname,'w') as outf:
    firstrun = True;
    for fname in files:
        with h5py.File(fname,'r') as f:
            dset_list = list(f.keys())
            if firstrun:
                prepare_ouput_file(f, outf, dset_list, nsim_tot)
                firstrun = False
            nsim_tot = nsim_tot + print_ouput_file(f, outf, dset_list) # here is the printing

        # shutil.move(fname, 'temp/'+fname) #moving files

print("total number of TDSE's merged:",nsim_tot)
print('Done')