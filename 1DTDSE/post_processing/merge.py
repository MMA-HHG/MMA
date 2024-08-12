####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2020)
#
# The purpose of this code is to merge outputs generated from 1DTDSE. It uses a general
# procedure for assigning indices based on keys, it thus can be easily generalised for more
# dimensions.



# finding mn from external location
# https://stackoverflow.com/questions/17198319/how-to-configure-custom-pythonpath-with-vm-and-pycharm

import sys
import os
import numpy as np
import h5py
import mynumerics as mn
import glob
# import argparse
import MMA_administration as MMA

# ap = argparse.ArgumentParser()
# ap.add_argument("-p", "--printdata", nargs='+', required=False, help="Select data to print. Available " \
#                 "options are 'Efield', 'FEfield', 'SourceTerm', 'FSourceTerm', 'FEfieldM2', " \
#                 "'FSourceTermM2', 'PopTot', " \
#                 "'PopInt', 'expval_x'")
# args = vars(ap.parse_args())

arguments = sys.argv

if not(os.path.isfile('msg.tmp')): raise FileNotFoundError("'msg.tmp' required to pass the name off the archive.")
with open('msg.tmp') as f: outfname = f.readline().rstrip()  # "results_merged.h5"

outfname = 'merge_test.h5'

h5path = MMA.paths['CTDSE_outputs']

available_outputs = {'Efield'          : '[a.u.]',
                     'FEfield'         : '[a.u.]',
                     'SourceTerm'      : '[a.u.]',
                     'FSourceTerm'     : '[a.u.]',
                     'FEfieldM2'       : '[a.u.]',
                     'FSourceTermM2'   : '[a.u.]',
                     'PopInt'          : '[-]', 
                     'expval_x'        : '[a.u.]'}

available_further_data = {'tgrid'                       : '[a.u.]',
                          'omegagrid'                   : '[a.u.]',
                          'Energy_of_the_ground_state'  : '[a.u.]',
                          'xgrid_micro'                 : '[a.u.]',
                          'ground_state'                : '[a.u.]',
                          'zgrid_coarse'                : '[a.u.]',
                          'rgrid_coarse'                : '[a.u.]',
                          'trg_a'                       : '[a.u.]'} 

# if args['printdata'] != None: 
#     args = args['printdata']
#     available_outputs_list = set(args).intersection(available_outputs_list)
#     if available_outputs_list == set():
#         available_outputs_list = ['Efield', 'FEfield', 'SourceTerm', 'FSourceTerm', 'FEfieldM2', 'FSourceTermM2', 'PopTot',
#                           'PopInt', 'expval_x'] 
#     print("Datasets to be printed: ")
#     print(available_outputs_list)


files = glob.glob('hdf5_temp_*.h5') # filter all the single-proc files
# outfname = "results_merged.h5"

precision = 'd'

def prepare_ouput_file(f,outf,dset_list):
    joint_zr_shape = (mn.readscalardataset(f,'Nz_orig','N')[0],mn.readscalardataset(f,'Nr_orig','N')[0])
    for dsetname in dset_list:

        if (dsetname in available_outputs.keys()):
            dset = f[dsetname]
            newshape =  joint_zr_shape + dset.shape[0:-1]
            newdset = outf.create_dataset(h5path+'/'+dsetname, newshape,precision)
            newdset.attrs['units'] = np.string_(available_outputs[dsetname]) # add units

        if (dsetname in available_further_data.keys()):
            f.copy(dsetname,outf[h5path])
            outf[h5path+'/'+dsetname].attrs['units'] = np.string_(available_further_data[dsetname]) # add units
            


def print_ouput_file(f,outf,dset_list):
    Nr = mn.readscalardataset(f, 'Nr_orig', 'N')[0]
    Nz = mn.readscalardataset(f, 'Nz_orig', 'N')[0]
    nsim_loc = mn.readscalardataset(f, 'number_of_local_simulations', 'N')[0]
    ks_loc = f['keys'][()]
    for dsetname in dset_list:
        if (dsetname in available_outputs.keys()):
            dset = f[dsetname]
            ndims = len(dset.shape)
            newdset = outf[h5path+'/'+dsetname]
            data = dset[()]
            for k1 in range(nsim_loc):
                if (ks_loc[k1] >= 0): # only for used keys
                    kr, kz = mn.n1n2mapping_inv(ks_loc[k1], Nr)
                    # print(kr,kz)
                    if (ndims == 2): # for reals
                        newdset[kz,kr,:] = dset[:,k1]
                    elif (ndims == 3):  # for complex
                        newdset[kz,kr,:,:] = dset[:,:,k1]
                    else:
                        print('warning, dataset with unsupported dimension: ' + dsetname + ', nothing done')# general procedure should be possible
    return nsim_loc

nsim_tot = 0
with h5py.File(outfname,'a') as outf:
    firstrun = True
    for fname in files:
        with h5py.File(fname,'r') as f:
            dset_list = list(f.keys())
            if firstrun:
                outf.require_group(h5path)
                prepare_ouput_file(f, outf, dset_list)
                firstrun = False
            nsim_tot = nsim_tot + print_ouput_file(f, outf, dset_list) # here is the printing

print("total number of TDSE's merged:",nsim_tot)


if not('-keep-files' in arguments):
    for f in files: os.remove(f)
    print("Intermediate files deleted.")
else:
     print("Intermediate files kept.")
     
print('Done')