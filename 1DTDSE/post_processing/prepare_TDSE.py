import numpy as np
import h5py
import sys
import os
import shutil
import re
import glob
import subprocess







results_file = glob.glob('../results_*.h5')
print(results_file[0])

with h5py.File(results_file[0], 'r') as InputArchive:
    Nz = len(InputArchive['/outputs/zgrid'][:])
    # Nz = 1500
    # pass


inp_part_path = os.environ['TDSE_1D_HOME']+'/processing/TDSE_scan1/FreeFormInputsTDSE_Xpl_stab.inp'
content = 'Nz_max\t'+str(Nz)+'\tI\t-\n' # Nz_max  1250 I   -

with open(inp_part_path,'r') as inp_part, open('TDSE.inp.tmp','w') as inp_tmp:
    inp_tmp.write(inp_part.read()+'\n'+content)


run_args = ['python3', os.environ['UNIV_INPUT_PATH']+'/create_universal_HDF5.py',
            '-i', 'TDSE.inp.tmp',
            '-ihdf5', results_file[0],
            '-ohdf5', 'results.h5',
            '-g', 'TDSE_inputs']

# print(run_args)
subprocess.run(run_args)

print('done')