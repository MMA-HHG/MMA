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


files = glob.glob('*.h5') # filter all available files
results_path = 'exit'


for fname in files:
    with h5py.File(fname,'r') as inpf, h5py.File(os.path.join(results_path,fname),'w') as outf:
        inpgrp = inpf['XUV']
        outgrp = outf.create_group('XUV')
        dset_list = list(inpgrp.keys())
        dset_list.remove('Spectrum_on_screen_cummulative')
        for dset in dset_list:
            inpgrp.cop(dset,outgrp)



print('Done')