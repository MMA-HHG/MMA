import numpy as np
import os
import math
import shutil
import h5py
import sys
import units
import mynumerics as mn
import glob


folders = glob.glob('*_end') + glob.glob('*_half') + glob.glob('*_no') # filter all available files

print(folders)

results_path = 'Hankels'


for folder in folders:
    # path = os.path.join(folder, 'TDSEs')
    Hankels = glob.glob(os.path.join(folder, 'TDSEs','Hankel_all*.h5'))
    print(Hankels)
    for Hankel in Hankels:
        print(os.path.basename(Hankel).replace('.h5','_'+folder+'.h5'))
        new_name = os.path.basename(Hankel).replace('.h5','_'+folder+'.h5')
        shutil.copyfile(Hankel, os.path.join(results_path, new_name))

# for fname in files:
#     with h5py.File(fname,'r') as inpf, h5py.File(os.path.join(results_path,fname),'w') as outf:
#         inpgrp = inpf['XUV']
#         outgrp = outf.create_group('XUV')
#         dset_list = list(inpgrp.keys())
#         dset_list.remove('Spectrum_on_screen_cummulative')
#         for dset in dset_list:
#             inpgrp.copy(dset,outgrp)



print('Done')