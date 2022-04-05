import numpy as np
import os
import math
import shutil
import h5py
import sys


import glob


# files = glob.glob('*.h5') # filter all available files
# results_path = 'exit'

fname1 = 'results.h5'
fname2 = 'results_merged.h5'
outfile = 'test.h5'

precision = 'd'


with h5py.File(fname1,'r') as inpf, h5py.File(fname2,'r') as inpf2, h5py.File(outfile,'w') as outf:
    inpf.copy('inputs',outf) # copy inps
    
    
    inpgrp = inpf['outputs']       
    outgrp = outf.create_group('outputs')
    zgrid = inpgrp['zgrid']; Nz = len(zgrid)
    # infield = inpgrp['output_field'][Nz-1,:,:]
    outfield = inpgrp['output_field'][[0,Nz-1],:,:]
    outplasma = inpgrp['output_plasma'][[0,Nz-1],:,:]
    
    outgrp.create_dataset('field_ee', data = outfield)
    outgrp.create_dataset('plasma_ee', data = outplasma)
    outgrp.create_dataset('zgrid', data = zgrid[[0,Nz-1]])
    inpgrp.copy('rgrid',outgrp)
    inpgrp.copy('tgrid',outgrp)
    
    outgrp = outf.create_group('TDSE')
    zgrid = ['zgrid_coarse']; Nz = len(zgrid)
    for dataset in ['tgrid','omegagrid','rgrid_coarse', 'ground_state', 'xgrid_micro']:
        inpf2.copy(dataset,outgrp)
    
    # omega domain
    print([0,Nz-1])
    dum = inpf2['FSourceTerm'][:,[0,Nz-1],:,:]
    outgrp.create_dataset('FSourceTerm', data = dum)
    
    # time domain
    for dataset in ['PopInt','PopTot','expval_x']:
        dum = inpf2[dataset][:,[0,Nz-1],:]
        outgrp.create_dataset(dataset, data = dum)
    
        
        # inpgrp.copy
        
        
        
        # dset_list = list(inpgrp.keys())
        # dset_list.remove('Spectrum_on_screen_cummulative')
        # for dset in dset_list:
        #     inpgrp.copy(dset,outgrp)



print('Done')